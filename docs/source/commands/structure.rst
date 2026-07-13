structure
=========

Analyze a single topology or structure file without requiring a trajectory.
Results are written to a ``.psa`` (Pharmacon Structure Analysis) HDF5 artifact.

**Supported formats:** ``.mol2``, ``.smi`` (``properties``); ``.pdb``, ``.gro``,
``.crd`` (``sequence``).  ``.sdf`` is accepted by the parser but the reader is
not yet implemented — use ``.mol2`` or ``.smi``.

.. contents:: Subcommands
   :local:
   :depth: 1

----

properties
----------

Compute structural and chemical descriptors for small molecules using RDKit.

**Computed properties:**

- Molecular weight, LogP, TPSA (topological polar surface area)
- Rotatable bond count, ring count, aromaticity
- Stereocenter count, net charge
- Element composition, molecular volume, fragment counts

**Computed fingerprints:**

- Morgan / ECFP
- Topological torsion
- Atom pair
- MACCS keys

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-i / --input``
     - Yes
     - Input small-molecule file (``.mol2`` or ``.smi``; ``.sdf`` not yet supported)
   * - ``-o / --output``
     - No
     - Output ``.psa`` file (default: ``properties.psa``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-l / --log``
     - No
     - Log file (default: ``properties.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Compute properties for a single MOL2 ligand:

.. code-block:: bash

   pharmacon structure properties \
       -i ligand.mol2 \
       -o ligand_props.psa

Process a SMILES file with multiple molecules:

.. code-block:: bash

   pharmacon structure properties \
       -i ligands.smi \
       -o ligand_props.psa

Then export to CSV:

.. code-block:: bash

   pharmacon export psa -i ligand_props.psa -f csv -o ./results/

Overwrite an existing output:

.. code-block:: bash

   pharmacon structure properties \
       -i ligands.smi \
       -o ligand_props.psa \
       --overwrite

----

sequence
--------

Extract the amino-acid sequence from a protein topology, organized by chain.
Outputs single-letter and three-letter representations with residue ID
mappings.  Supports FASTA export via ``pharmacon export psa``.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-p / --topology``
     - Yes
     - Input topology file (``.pdb``, ``.gro``, ``.crd``)
   * - ``-o / --output``
     - No
     - Output ``.psa`` file (default: ``sequence.psa``)
   * - ``--overwrite``
     - No
     - Overwrite existing output file
   * - ``-l / --log``
     - No
     - Log file (default: ``sequence.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Extract sequence from a PDB file:

.. code-block:: bash

   pharmacon structure sequence \
       -p protein.pdb \
       -o sequence.psa

Extract from a GROMACS structure:

.. code-block:: bash

   pharmacon structure sequence \
       -p protein.gro \
       -o sequence.psa

Export the sequence to FASTA for downstream bioinformatics:

.. code-block:: bash

   pharmacon structure sequence -p protein.pdb -o sequence.psa
   pharmacon export psa -i sequence.psa -f fasta -o ./results/
