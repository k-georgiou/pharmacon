Periodic Boundary Conditions (PBC)
====================================

Most MD engines write trajectories in which solute molecules diffuse across
the edges of the simulation box, get wrapped back through the opposite face,
or end up split between periodic images. Distances, angles, and contacts
measured on such raw frames are physically meaningless unless the box
geometry is taken into account. Pharmacon handles this in two complementary
ways: **geometry-aware calculations** that consult the box at every frame,
and an **opt-in transformation pipeline** that rewrites coordinates before
the analysis sees them.

PBC-aware calculations
----------------------

The trajectory analyzers that measure geometric quantities — ``distances``,
``angles``, ``h-bonds``, ``pl-interactions``, ``pp-interactions`` — read the
per-frame unit-cell dimensions from the trajectory and use the
**minimum-image convention** when box information is available. Concretely,
this means a hydrogen-bond donor in the primary image and an acceptor in a
neighbouring image are still detected as a single contact at the correct
distance, without any user intervention.

If the trajectory has no box dimensions (e.g. a stripped ``.dcd`` or a
gas-phase run), Pharmacon falls back to plain Euclidean geometry and logs
the fact at ``DEBUG`` level.

The ``-at / --add-transformations`` flag
-----------------------------------------

PBC-aware distances handle the *measurement* side, but they cannot help when
a molecule is visually broken across box faces (e.g. a protein straddling
the box boundary), when downstream tools require whole molecules (RDKit
SMARTS matching, RMSD fitting, average-structure export), or when the
ligand drifts out of the protein's primary image. For these cases most
trajectory subcommands that consume coordinates (all except ``rmsd`` and
``rmsf``) accept:

.. code-block:: text

   -at, --add-transformations    Apply PBC unwrapping transformations to the
                                 MDAnalysis Universe before analysis.
                                 (default: False)

When the flag is set, Pharmacon attaches an **on-the-fly transformation
pipeline** to the MDAnalysis ``Universe`` in
``pharmacon.utils.mda.create_universe()``. The pipeline is built and applied
in a fixed order:

1. **unwrap** — ``MDAnalysis.transformations.unwrap(u.atoms)``. Makes every
   molecule whole by reversing the wrapping that the MD engine applied on
   write. Required before any per-molecule operation (RMSD fit, COM/COG
   distance, SMARTS perception).
2. **center_in_box** — ``center_in_box(ref, center="geometry")``.
   Re-centers the system on the geometric centre of the protein (falling
   back to all atoms if no ``protein`` selection exists). This keeps the
   receptor stationary in the primary image and prevents it from drifting
   toward a box face over long simulations.
3. **wrap** — ``MDAnalysis.transformations.wrap(u.atoms, compound="atoms")``.
   After unwrapping and centering, wraps everything back into the primary
   image so that downstream distance/contact searches stay numerically
   well-behaved.

The transformations are **lazy**: they are attached once to
``u.trajectory`` and re-evaluated on every frame iteration, so memory usage
stays flat regardless of trajectory length. If the pipeline cannot be
applied (e.g. missing bond topology required by ``unwrap``), Pharmacon logs
a warning and continues with the untransformed trajectory rather than
aborting the run — the flag is treated as a convenience, not a hard
requirement.

When to enable ``--add-transformations``
-----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Situation
     - Recommendation
   * - Trajectory already imaged / centered by the MD engine
     - Leave **off** (default)
   * - Protein visibly split across box edges in the input trajectory
     - Enable
   * - Ligand drifts out of the primary box during a long run
     - Enable
   * - PCA or average-structure on a periodic system
     - Enable
   * - Protein–ligand or protein–protein interactions with periodic ligands
     - Enable
   * - Box dimensions are missing from the trajectory
     - Has no effect — pipeline is silently skipped

The flag is *idempotent* on already-imaged trajectories: re-running with
``-at`` on a clean trajectory simply re-applies the same identity pass and
produces the same results.

Recommended: pre-image the trajectory yourself
-----------------------------------------------

The on-the-fly ``unwrap → center → wrap`` pipeline runs on **every** frame
iteration, every time you launch an analysis. On long trajectories or
workflows that touch the same trajectory repeatedly (RMSD, PCA, then
PL-interactions, then H-bonds, …) this overhead is paid over and over again
and can **heavily affect calculation speed** — sometimes dominating the
total runtime.

For best performance we **strongly recommend pre-processing the trajectory
once with your MD engine's native tooling** and then feeding the cleaned-up
trajectory to Pharmacon *without* ``-at``. Native tools are implemented in
C/Fortran, operate frame-by-frame in a single streaming pass, and write the
result to disk so the cost is paid exactly once.

Examples of native equivalents:

- **GROMACS** — ``gmx trjconv -pbc whole``, then
  ``gmx trjconv -center -pbc mol -ur compact`` (or ``-pbc nojump`` for
  unwrapped trajectories).
- **Amber / cpptraj** — ``autoimage`` followed by ``center`` and
  ``image familiar``.
- **CHARMM / NAMD (VMD)** — the ``pbctools`` plugin: ``pbc unwrap``,
  ``pbc wrap -center com -centersel "protein" -compound res``.
- **MDAnalysis (Python)** — write a short script that mirrors the same
  ``unwrap → center_in_box → wrap`` chain and dumps the result with an
  ``MDAnalysis.coordinates.XTC.XTCWriter``.

.. note::

   If you're not sure how to do this with your engine of choice, it is
   perfectly fine to just pass ``-at`` and let Pharmacon handle it for you
   — the result is numerically identical, you just have to wait a bit
   longer for the analysis to finish. The flag exists exactly for this
   case.

Subcommands that accept ``-at``
--------------------------------

The flag is exposed on every trajectory subcommand whose output depends on
coordinate geometry:

- ``pharmacon trajectory distances``
- ``pharmacon trajectory angles``
- ``pharmacon trajectory h-bonds``
- ``pharmacon trajectory pl-interactions``
- ``pharmacon trajectory pp-interactions``
- ``pharmacon trajectory pca``
- ``pharmacon trajectory average-st``

Example:

.. code-block:: bash

   pharmacon trajectory pl-interactions \
       -p out.dms -x out.xtc -o pli.pta \
       -prt "chainid R" -lig "resname LIG" \
       -w "resname T3P and around 5 resname LIG" \
       --add-transformations \
       --workers 10

The ``add_transformations`` setting is persisted in the resulting ``.pta``
file's root-level metadata, so a downstream consumer can verify which PBC
treatment a given artifact was produced under without re-reading the
original trajectory.