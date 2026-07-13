File Formats
============

Pharmacon persists every analysis result as a **signed, self-describing HDF5
artifact**.  There are two concrete formats — ``.pta`` for trajectory analysis
and ``.psa`` for structure analysis — both built on top of the shared
``PharmaconHDF5File`` base class.

Schema version: ``1.0``

Both formats are plain HDF5 files and can be opened with ``h5py``, ``h5dump``,
``HDFView``, or any HDF5-aware library.  Pharmacon-specific metadata lives in
HDF5 attributes on the root group and on each analysis group.

.. contents::
   :local:
   :depth: 1

----

Common root-level attributes
-----------------------------

Every Pharmacon artifact (``.pta`` or ``.psa``) carries the following metadata
block on the root group, written automatically at creation time.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Attribute
     - Type
     - Meaning
   * - ``pharmacon_version``
     - str
     - Pharmacon version that wrote the file (e.g. ``"1.0.0"``)
   * - ``schema_version``
     - str
     - HDF5 layout version (``"1.0"``)
   * - ``file_type``
     - str
     - ``TRAJECTORY_ANALYSIS`` or ``STRUCTURE_ANALYSIS``
   * - ``command``
     - str
     - Top-level command (e.g. ``"Trajectory Analysis"``)
   * - ``subcommand``
     - str
     - Specific subcommand (e.g. ``"rmsd"``, ``"pl-interactions"``)
   * - ``description``
     - str
     - Human-readable one-liner describing the run
   * - ``created_at``
     - str
     - ISO-8601 UTC timestamp at which the file was opened for writing
   * - ``signature``
     - str
     - ``PHARMACON::<COMMAND>::<SUBCOMMAND>::<chunked-SHA256-prefix>`` identity token
   * - ``fingerprint``
     - str
     - 128-bit BLAKE2b hex digest — tamper-check hint for the identity block
   * - ``blueprint``
     - str
     - Deterministic hash of the run inputs (topology, selections, flags).
       Used by ``pharmacon merge results`` to refuse incompatible files.
   * - ``artifact_token``
     - str
     - HMAC-style token bound to the blueprint
   * - ``artifact_token_version``
     - str
     - Token scheme version (currently ``"1"``)
   * - ``artifact_status``
     - str
     - ``"SUCCESS"`` if the run completed cleanly
   * - ``artifact_status_code``
     - str
     - Exit code (``"0"`` on success)
   * - ``completed``
     - str
     - ``"True"`` once the writer has closed the file with all groups populated
   * - ``is_merged``
     - str
     - ``"True"`` if the file was produced by ``pharmacon merge results``

In addition, each subcommand appends its own run-specific attributes (e.g.
``begin``, ``end``, ``step``, ``total_frames``, ``ligand``, ``protein``,
``water``, ``fitting_group``, …) so a file is fully reproducible from its
metadata alone.

Signature, fingerprint, and blueprint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These three tokens look similar but answer different questions:

- **signature** — a stable *identity* string built from
  ``PHARMACON::<format>::<version>::<command>::<subcommand>::<salt>``,
  prefix-chunked for readability.  Two files with the same signature were
  produced by the same command/subcommand flavour of Pharmacon.
- **fingerprint** — a 128-bit BLAKE2b over the same payload.  Acts as a cheap
  tamper-check hint for the identity block.
- **blueprint** — hashes the *actual inputs* (topology, trajectory,
  selections, flags).  Two files share a blueprint only if they were run
  against equivalent inputs.  ``pharmacon merge results`` uses this to refuse
  merging incompatible replicates.

None of these are cryptographic signatures — they are integrity hints, not
authentication.

----

PTA — Pharmacon Trajectory Analysis
-------------------------------------

Extension: ``.pta``

Stores results from trajectory-based analyses.  Each subcommand writes one or
more **analysis groups** at the root, containing per-frame datasets.

**General layout:**

.. code-block:: text

   results.pta
   ├── @ pharmacon_version, signature, fingerprint, blueprint,
   │     command, subcommand, begin, end, step, total_frames, …
   └── <analysis_group>/
       ├── @ completed=True
       ├── frame_0/
       │   └── <dataset>
       ├── frame_1/
       │   └── <dataset>
       └── …

**Group names per subcommand:**

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Subcommand
     - Root group
     - Per-frame dataset(s)
   * - ``rmsd``
     - ``rmsd``
     - ``frame_<N>/rmsd``
   * - ``distances``
     - ``distances``
     - ``frame_<N>/distances``
   * - ``angles``
     - ``angles``
     - ``frame_<N>/angles``
   * - ``h-bonds``
     - ``hydrogen_bonds``
     - ``frame_<N>/interactions``
   * - ``pl-interactions``
     - ``pl_interactions``
     - ``frame_<N>/interactions``
   * - ``pp-interactions``
     - ``pp_interactions``
     - ``frame_<N>/interactions``
   * - ``pca``
     - ``pca_<selection>`` (one group per selection; bare ``pca`` accepted for back-compat)
     - ``components``, ``projections``, ``variance_ratio``

.. _interaction-schema:

Interaction row schema
^^^^^^^^^^^^^^^^^^^^^^

The ``interactions`` dataset for ``pl-interactions``, ``pp-interactions``,
and ``h-bonds`` is a variable-length array of rows with a fixed 21-column
core schema:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``frame_number``
     - Frame index in the source trajectory
   * - ``interaction``
     - Interaction type tag (see table below)
   * - ``atom1_index``
     - Atom index (first partner)
   * - ``atom1_name``
     - Atom name
   * - ``atom1_id``
     - Atom ID
   * - ``atom1_type``
     - Atom type
   * - ``atom1_element``
     - Element symbol
   * - ``atom1_resname``
     - Residue name
   * - ``atom1_resid``
     - Residue ID
   * - ``atom1_chainid``
     - Chain ID
   * - ``atom1_segid``
     - Segment ID
   * - ``atom2_*``
     - Same 9 fields for the second partner
   * - ``details``
     - Dict of interaction-specific geometry fields (see below)

**Details fields per interaction type:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Interaction tag
     - ``details`` fields
   * - ``HYDROPHOBIC``
     - ``distance``, ``is_hydrogen``
   * - ``HYDROGEN-BOND``
     - ``distance``, ``angle_dha``, ``angle_hax``, ``orientation``
   * - ``IONIC``
     - ``distance``, ``orientation``
   * - ``HALOGEN-BOND``
     - ``distance``, ``angle_cxa``, ``angle_xay``, ``orientation``
   * - ``METAL-CONTACT``
     - ``distance``, ``role``, ``orientation``
   * - ``WATER-BRIDGE-1``
     - 16 fields: bridging water + two distances + two angle pairs + orientation
   * - ``WATER-BRIDGE-2``
     - 28 fields: two bridging waters + three distances + three angle pairs + orientation *(reserved — not yet implemented)*
   * - ``PI-CATION``
     - ``distance``, ``theta``, ``orientation``
   * - ``PI-STACKING``
     - ``distance``, ``angle``, ``stacking_type``, ``orientation``

Aggregated interaction views
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For interaction analyses, Pharmacon additionally writes **mode tables**
alongside the per-frame data:

- ``modes/mode<N>/table`` — per-file aggregated tables (residue × frame,
  residue × interaction type, residue × ligand atom, …) computed at write
  time so plotters do not need to re-walk every frame.
- ``modes_merged/<mode_name>/table`` — present only on files produced by
  ``pharmacon merge results``.  Contains the aggregate across all merged
  replicates.  Merged files drop per-frame datasets in favour of these tables.

----

PSA — Pharmacon Structure Analysis
-------------------------------------

Extension: ``.psa``

Stores results from static structure analyses.  Inherits the entire root
attribute schema from ``.pta`` and adds structure-specific groups.

**Sequence layout:**

.. code-block:: text

   sequence.psa
   ├── @ pharmacon_version, signature, …, command="Structure Analysis", subcommand="sequence"
   └── sequence/
       ├── @ completed=True
       └── A/                    ← one subgroup per chain
           ├── aa1_seq           ← single-letter sequence
           ├── aa3_list          ← three-letter residue list
           └── resid_seq         ← residue ID sequence

**Properties layout:**

.. code-block:: text

   properties.psa
   ├── @ pharmacon_version, signature, …
   └── properties/
       ├── @ completed=True
       ├── properties_table      ← scalar descriptors per molecule
       └── fingerprints/
           ├── morgan
           ├── topological_torsion
           ├── atom_pair
           └── maccs_keys

----

Inspecting artifacts
---------------------

Using the Pharmacon CLI (rich formatted panels):

.. code-block:: bash

   pharmacon dump pta -i run.pta
   pharmacon dump psa -i structure.psa

Using h5py directly:

.. code-block:: python

   import h5py

   with h5py.File("run.pta", "r") as f:
       print(dict(f.attrs))       # root metadata
       print(list(f.keys()))      # analysis groups

----

Exporting artifacts
--------------------

.. code-block:: bash

   pharmacon export pta -i rmsd.pta            -f csv   -o ./results/
   pharmacon export pta -i pl_interactions.pta -f tsv   -o ./results/
   pharmacon export psa -i sequence.psa        -f fasta -o ./results/
   pharmacon export psa -i properties.psa      -f csv   -o ./results/
