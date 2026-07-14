trajectory
==========

Analyze MD trajectories frame by frame.  Every subcommand reads a topology
and a trajectory file and writes results to a signed ``.pta`` HDF5 artifact.

**Supported topology formats:** ``.tpr``, ``.prmtop``, ``.parm7``, ``.psf``, ``.dms``

**Supported trajectory formats:** ``.xtc``, ``.trr``, ``.dcd``, ``.nc``

.. tip::

   **Pre-image your trajectory with the MD engine instead of relying on**
   ``-at``. Several subcommands accept ``-at / --add-transformations`` to
   unwrap → center → wrap a periodic system on the fly, but that pipeline is
   **lazy**: MDAnalysis re-evaluates it on *every* frame of *every* run, so the
   cost is paid again on each analysis and can dominate the total runtime.

   The strongly recommended workflow is to image the trajectory **once** with
   your engine's native, C/Fortran tooling and feed that cleaned trajectory to
   Pharmacon **without** ``-at`` — it is dramatically faster and numerically
   identical:

   - **GROMACS** — ``gmx trjconv -pbc whole`` then
     ``gmx trjconv -center -pbc mol -ur compact``
   - **Amber / cpptraj** — ``autoimage`` then ``center`` + ``image familiar``
   - **CHARMM / NAMD (VMD)** — the ``pbctools`` plugin (``pbc unwrap`` /
     ``pbc wrap``)

   Reach for ``-at`` only as a convenience when you cannot pre-process. See
   :doc:`/pbc` for the full rationale and engine-specific commands.

.. contents:: Subcommands
   :local:
   :depth: 1

----

rmsd
----

Calculate Root Mean Square Deviation of atomic positions relative to a
reference frame.  Supports multiple selections per run with a dedicated
fitting (alignment) group.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selections``
     - Yes
     - One or more MDAnalysis selection strings for RMSD calculation
   * - ``-n / --names``
     - Yes
     - Label for each selection (must match count of ``-sel``)
   * - ``-f / --fitting-group``
     - Yes
     - MDAnalysis selection used for structural alignment before RMSD
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``rmsd.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-r / --reference-frame``
     - No
     - Frame index used as RMSD reference (default: 0)
   * - ``-b / --begin``
     - No
     - First frame to process (default: 0)
   * - ``-e / --end``
     - No
     - Last frame to process (default: last)
   * - ``-s / --step``
     - No
     - Process every Nth frame (default: 1)
   * - ``-l / --log``
     - No
     - Log file (default: ``rmsd.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

RMSD of C-alpha atoms aligned on C-alpha:

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p topol.tpr \
       -x traj.xtc \
       -o rmsd.pta \
       -sel "protein and name CA" \
       -f  "protein and name CA" \
       -n  calpha

Multiple selections in a single run (backbone, C-alpha, and ligand):

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p topol.tpr \
       -x traj.xtc \
       -o rmsd.pta \
       -sel "protein and name CA" "backbone" "resname LIG and not name H*" \
       -f  "protein and name CA" \
       -n  calpha backbone ligand_heavy

Process frames 100–500 with stride 5, using frame 100 as reference:

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p topol.tpr \
       -x traj.xtc \
       -o rmsd_subset.pta \
       -sel "protein and name CA" \
       -f  "protein and name CA" \
       -n  calpha \
       -b 100 -e 500 -s 5 -r 100

Full run with debug logging:

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p protein.prmtop \
       -x md.nc \
       -sel "protein and name CA" "resname LIG and not name H*" \
       -n  bb_ca ligand_heavy \
       -f  "protein and name CA" \
       -r 0 -b 0 -e 10000 -s 10 \
       -o  rmsd.pta \
       -l  rmsd.log -fl DEBUG -tl WARNING

----

rmsf
----

Calculate per-atom Root Mean Square Fluctuation (RMSF) over a trajectory.
The trajectory is aligned in-memory on the fitting group, then ``rms.RMSF``
is computed for each selection.  By default alignment is two-pass (align to
the initial frame, build the average structure, then re-align to it); pass
``-r`` to align to a single reference frame instead.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selections``
     - Yes
     - One or more MDAnalysis selection strings for per-atom RMSF
   * - ``-n / --names``
     - Yes
     - Label for each selection (must match count of ``-sel``)
   * - ``-f / --fitting-group``
     - Yes
     - MDAnalysis selection used for in-memory alignment before RMSF
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``rmsf.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-r / --reference-frame``
     - No
     - Frame index used as the alignment reference; omit to align to the
       average structure (two-pass, recommended)
   * - ``-b / --begin``
     - No
     - First frame to process (default: 0)
   * - ``-e / --end``
     - No
     - Last frame to process (default: last)
   * - ``-s / --step``
     - No
     - Process every Nth frame (default: 1)
   * - ``-l / --log``
     - No
     - Log file (default: ``rmsf.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Per-atom RMSF of C-alpha atoms (two-pass alignment to the average structure):

.. code-block:: bash

   pharmacon trajectory rmsf \
       -p topol.tpr \
       -x traj.xtc \
       -o rmsf.pta \
       -sel "protein and name CA" \
       -f  "protein and name CA" \
       -n  calpha

Multiple selections aligned to a single reference frame:

.. code-block:: bash

   pharmacon trajectory rmsf \
       -p topol.tpr \
       -x traj.xtc \
       -o rmsf.pta \
       -sel "protein and name CA" "backbone" \
       -f  "protein and name CA" \
       -n  calpha bb \
       -r 0

----

distances
---------

Compute pairwise or group-to-group distance time series between user-defined
atom selections.  Supports minimum-distance, maximum-distance,
centre-of-mass, and centre-of-geometry modes with PBC awareness.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel1 / --selection1``
     - Yes
     - First atom group selection(s)
   * - ``-sel2 / --selection2``
     - Yes
     - Second atom group selection(s)
   * - ``-m / --method``
     - Yes
     - Distance method per pair: ``MIN``, ``MAX``, ``COM``, ``COG``
   * - ``-n / --names``
     - Yes
     - Label for each distance pair
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``distances.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-l / --log``
     - No
     - Log file (default: ``distances.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

All four lists (``-sel1``, ``-sel2``, ``-m``, ``-n``) must have the same
number of entries.

Distance methods:

- **MIN** — minimum distance between any atom pair across the two groups
- **MAX** — maximum distance between any atom pair across the two groups
- **COM** — distance between the centres of mass
- **COG** — distance between the centres of geometry

**Examples**

Centre-of-mass distance between ligand and C-alpha atoms:

.. code-block:: bash

   pharmacon trajectory distances \
       -p topol.tpr \
       -x traj.xtc \
       -o distances.pta \
       -sel1 "resname LIG" \
       -sel2 "protein and name CA" \
       -m   COM \
       -n   lig_to_ca

Multiple distance pairs in one run:

.. code-block:: bash

   pharmacon trajectory distances \
       -p topol.tpr \
       -x traj.xtc \
       -o distances.pta \
       -sel1 "resname LIG" "resname LIG" \
       -sel2 "resid 45"    "resid 102" \
       -m    MIN            MIN \
       -n    lig_res45      lig_res102

With PBC unwrapping enabled:

.. code-block:: bash

   pharmacon trajectory distances \
       -p topol.tpr \
       -x traj.xtc \
       -o distances_pbc.pta \
       -sel1 "protein" \
       -sel2 "resname LIG" \
       -m   COM \
       -n   prot_lig \
       --add-transformations

----

angles
------

Measure geometric angles across a trajectory.  Three measurement types are
supported:

- **Three-atom angle** — select three atoms; Pharmacon measures the angle at the central atom.
  Example selection: ``"index 10 or index 11 or index 12"``
- **Vector angle** — define two vectors using the ``->`` syntax; Pharmacon measures the angle between the two vectors.
  Example: ``"index 10 or index 100 -> index 300 or index 301"``
- **Dihedral / torsion angle** — select four atoms; Pharmacon measures the torsion around the central bond.
  Example: ``"index 10 or index 11 or index 12 or index 13"``

All calculations are PBC-aware when box dimensions are available.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selections``
     - Yes
     - Selection strings defining each angle (see syntax above)
   * - ``-n / --names``
     - Yes
     - Label for each angle measurement
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``angles.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``-l / --log``
     - No
     - Log file (default: ``angles.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Vector angle relative to the x-axis:

.. code-block:: bash

   pharmacon trajectory angles \
       -p topol.tpr \
       -x traj.xtc \
       -o angles.pta \
       -sel "index 4 or index 5 -> x-axis" \
       -n  vector1

Dihedral angle of a backbone torsion:

.. code-block:: bash

   pharmacon trajectory angles \
       -p topol.tpr \
       -x traj.xtc \
       -o angles.pta \
       -sel "index 10 or index 11 or index 12 or index 13" \
       -n  phi_backbone

Multiple angle measurements in one run:

.. code-block:: bash

   pharmacon trajectory angles \
       -p topol.tpr \
       -x traj.xtc \
       -o angles.pta \
       -sel "index 4 or index 5 -> x-axis" "index 10 or index 11 -> z-axis" \
       -n  vec1 vec2

----

h-bonds
-------

Detect and count hydrogen bonds over a trajectory using donor–acceptor
geometry with distance and angle criteria.  Produces a focused output of
hydrogen-bond records only.  Useful when the hydrogen-bonding network is the
sole interest; for full interaction profiling see :ref:`pl-interactions`.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selection``
     - Yes
     - Atom group to analyse for hydrogen bonds
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``hbonds.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``--workers``
     - No
     - Number of parallel worker processes (default: 1)
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``-l / --log``
     - No
     - Log file (default: ``hbonds.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Hydrogen bonds within the protein:

.. code-block:: bash

   pharmacon trajectory h-bonds \
       -p topol.tpr \
       -x traj.xtc \
       -o hbonds.pta \
       -sel "protein"

Use 4 parallel workers for faster processing:

.. code-block:: bash

   pharmacon trajectory h-bonds \
       -p topol.tpr \
       -x traj.xtc \
       -o hbonds.pta \
       -sel "protein" \
       --workers 4

Hydrogen bonds in chain A, every 10th frame:

.. code-block:: bash

   pharmacon trajectory h-bonds \
       -p topol.tpr \
       -x traj.xtc \
       -o hbonds_chainA.pta \
       -sel "chainid A" \
       -s 10

----

.. _pl-interactions:

pl-interactions
---------------

Detect and quantify nine types of protein–ligand non-covalent interactions on
a per-frame basis.  Uses RDKit SMARTS pattern matching for accurate atom-type
classification.

**Detected interaction types:**

1. **Hydrophobic contacts** — distance-based between non-polar atoms
2. **Hydrogen bonds** — donor–hydrogen–acceptor geometry
3. **Ionic / salt-bridge** — oppositely charged groups within cutoff
4. **Halogen bonds** — C–X···A geometry with angle validation
5. **Metal contacts** — metal ion coordination to electronegative atoms
6. **Water bridges (1st degree)** — single bridging water mediating an H-bond chain
7. **Water bridges (2nd degree)** — two bridging waters (not yet implemented)
8. **Pi–cation** — aromatic ring to charged group
9. **Pi–stacking** — parallel or T-shaped aromatic ring pairs

Each record includes full atom metadata (name, residue, chain, segment) for
both partners plus interaction-specific geometric details.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-prt / --protein``
     - Yes
     - MDAnalysis selection for the protein / receptor
   * - ``-lig / --ligand``
     - Yes
     - MDAnalysis selection for the ligand
   * - ``-w / --water``
     - No
     - MDAnalysis selection for water molecules (required for water-bridge detection)
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``pl_interactions.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``--workers``
     - No
     - Number of parallel worker processes (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``--disable-hydrophobic``
     - No
     - Skip hydrophobic contact detection
   * - ``--disable-hbonds``
     - No
     - Skip hydrogen-bond detection
   * - ``--disable-pi-stacking``
     - No
     - Skip pi-stacking detection
   * - ``--disable-pi-cation``
     - No
     - Skip pi-cation detection
   * - ``--disable-ionic``
     - No
     - Skip ionic / salt-bridge detection
   * - ``--disable-water-bridges``
     - No
     - Skip water-bridge detection
   * - ``--disable-halogen``
     - No
     - Skip halogen-bond detection
   * - ``--disable-metal``
     - No
     - Skip metal-contact detection
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-l / --log``
     - No
     - Log file (default: ``pl_interactions.log``)

.. caution::

   **Always scope the water selection dynamically around the ligand.**
   Prefer ``-w "resname WAT and around 5 resname LIG"`` over a global
   ``-w "resname WAT"``.  A global selection forces Pharmacon to test
   **every** water molecule in the box for bridging on **every** frame, which
   can **dramatically increase runtime**.  Restricting waters to the
   neighbourhood of the binding site yields identical results for the relevant
   bridges at a fraction of the cost.  (Adjust ``WAT``/``LIG`` and the
   ``around`` cutoff to match your system.)

**Examples**

Full interaction profile with water bridges:

.. code-block:: bash

   pharmacon trajectory pl-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pl_interactions.pta \
       -prt "protein" \
       -lig "resname LIG" \
       -w  "resname WAT and around 5 resname LIG"

Specific chain against the ligand, PBC-corrected:

.. code-block:: bash

   pharmacon trajectory pl-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pl_interactions.pta \
       -prt "chainid A" \
       -lig "resname LIG" \
       --add-transformations

Hydrophobic and hydrogen bonds only (all others disabled):

.. code-block:: bash

   pharmacon trajectory pl-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pl_hbonds_hydro.pta \
       -prt "protein" \
       -lig "resname LIG" \
       --disable-pi-stacking \
       --disable-pi-cation \
       --disable-ionic \
       --disable-water-bridges \
       --disable-halogen \
       --disable-metal

Parallel processing with 8 workers:

.. code-block:: bash

   pharmacon trajectory pl-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pl_interactions.pta \
       -prt "protein" \
       -lig "resname LIG" \
       --workers 8

----

pp-interactions
---------------

Detect protein–protein non-covalent interactions between two chains or groups
over a trajectory.  Supports the same nine interaction types as
:ref:`pl-interactions`.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-prt1 / --protein1``
     - Yes
     - MDAnalysis selection for the first protein / chain
   * - ``-prt2 / --protein2``
     - Yes
     - MDAnalysis selection for the second protein / chain
   * - ``-w / --water``
     - No
     - MDAnalysis selection for water molecules (required for water-bridge detection)
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``pp_interactions.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``--workers``
     - No
     - Number of parallel worker processes (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``--disable-hydrophobic``
     - No
     - Skip hydrophobic contact detection
   * - ``--disable-hbonds``
     - No
     - Skip hydrogen-bond detection
   * - ``--disable-pi-stacking``
     - No
     - Skip pi-stacking detection
   * - ``--disable-pi-cation``
     - No
     - Skip pi-cation detection
   * - ``--disable-ionic``
     - No
     - Skip ionic / salt-bridge detection
   * - ``--disable-water-bridges``
     - No
     - Skip water-bridge detection
   * - ``--disable-halogen``
     - No
     - Skip halogen-bond detection
   * - ``--disable-metal``
     - No
     - Skip metal-contact detection
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-l / --log``
     - No
     - Log file (default: ``pp_interactions.log``)

.. caution::

   **Always scope the water selection dynamically around the interface.**
   Prefer ``-w "resname WAT and around 5 chainid B"`` over a global
   ``-w "resname WAT"``.  A global selection forces Pharmacon to test
   **every** water molecule in the box for bridging on **every** frame, which
   can **dramatically increase runtime**.  Restricting waters to the
   neighbourhood of the protein–protein interface yields identical results for
   the relevant bridges at a fraction of the cost.  (Adjust the resname and the
   ``around`` cutoff to match your system.)

**Examples**

Interactions between chain A and chain B:

.. code-block:: bash

   pharmacon trajectory pp-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pp_interactions.pta \
       -prt1 "chainid A" \
       -prt2 "chainid B"

Interface between two specific residue ranges:

.. code-block:: bash

   pharmacon trajectory pp-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pp_interface.pta \
       -prt1 "resid 1:50" \
       -prt2 "resid 51:100"

Parallel run keeping only H-bond and ionic interactions:

.. code-block:: bash

   pharmacon trajectory pp-interactions \
       -p topol.tpr \
       -x traj.xtc \
       -o pp_hbond_ionic.pta \
       -prt1 "chainid A" \
       -prt2 "chainid B" \
       --disable-hydrophobic \
       --disable-pi-stacking \
       --disable-pi-cation \
       --disable-water-bridges \
       --disable-halogen \
       --disable-metal \
       --workers 4

----

pca
---

Principal Component Analysis on atomic positions over a trajectory.  Outputs
eigenvalues (explained variance ratios), eigenvectors (principal components),
and per-frame projections.  Results support scatter plots, time-series
projections, variance scree plots, free-energy surface (FES) heatmaps, and
probability density heatmaps via ``pharmacon plot pta``.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selection``
     - Yes
     - MDAnalysis selection of atoms to include in PCA
   * - ``-c / --components``
     - No
     - Number of principal components to retain (default: 3)
   * - ``-o / --output``
     - No
     - Output ``.pta`` file (default: ``pca.pta``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``--parallel``
     - No
     - Number of parallel worker processes (default: 1)
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``-l / --log``
     - No
     - Log file (default: ``pca.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

PCA on C-alpha atoms, 5 components:

.. code-block:: bash

   pharmacon trajectory pca \
       -p topol.tpr \
       -x traj.xtc \
       -o pca.pta \
       -sel "protein and name CA" \
       -c 5

PCA on backbone, 10 components, every 5th frame:

.. code-block:: bash

   pharmacon trajectory pca \
       -p topol.tpr \
       -x traj.xtc \
       -o pca_bb.pta \
       -sel "backbone" \
       -c 10 \
       -s 5

Then generate all PCA plots:

.. code-block:: bash

   pharmacon plot pta -i pca.pta -o pca_plots/ --overwrite

----

average-st
----------

Compute the per-atom mean coordinates over a specified frame range, with
optional alignment to a reference frame.  Writes the average structure as a
coordinate file for visualisation or further analysis.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file
   * - ``-x / --trajectory``
     - Yes
     - Input trajectory file
   * - ``-sel / --selection``
     - Yes
     - MDAnalysis selection of atoms to average
   * - ``-o / --output``
     - No
     - Output coordinate file (default: ``average_st.pdb``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-r / --reference-frame``
     - No
     - Reference frame for alignment before averaging (default: 0)
   * - ``-b / --begin``
     - No
     - First frame (default: 0)
   * - ``-e / --end``
     - No
     - Last frame (default: last)
   * - ``-s / --step``
     - No
     - Frame stride (default: 1)
   * - ``-at / --add-transformations``
     - No
     - Apply PBC wrap/unwrap transformations
   * - ``-l / --log``
     - No
     - Log file (default: ``average_st.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Average C-alpha structure over the full trajectory:

.. code-block:: bash

   pharmacon trajectory average-st \
       -p topol.tpr \
       -x traj.xtc \
       -o average.pdb \
       -sel "protein and name CA" \
       -r 0

Average over the last 500 frames:

.. code-block:: bash

   pharmacon trajectory average-st \
       -p topol.tpr \
       -x traj.xtc \
       -o average_last500.pdb \
       -sel "protein" \
       -b 500

Average every 10th frame using frame 50 as reference:

.. code-block:: bash

   pharmacon trajectory average-st \
       -p topol.tpr \
       -x traj.xtc \
       -o average_stride.pdb \
       -sel "protein and name CA" \
       -r 50 -s 10
