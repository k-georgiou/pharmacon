Quick Start
===========

CLI structure
-------------

Pharmacon follows a **git-style** command structure:

.. code-block:: text

   pharmacon <command> <subcommand> [options]

Getting help at every level
-----------------------------

.. code-block:: bash

   pharmacon --help                            # list all commands
   pharmacon trajectory --help                 # list trajectory subcommands
   pharmacon trajectory rmsd --help            # show all RMSD flags
   pharmacon plot pta --help                   # show all plot flags

Supported MD engines
---------------------

Pharmacon reads trajectories produced by **Amber**, **GROMACS**, **CHARMM**,
**NAMD**, and **OpenMM** through MDAnalysis.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File type
     - Supported formats
   * - Topology
     - ``.tpr``, ``.prmtop``, ``.parm7``, ``.psf``, ``.dms``
   * - Trajectory
     - ``.xtc``, ``.trr``, ``.dcd``, ``.nc``
   * - Structure (static)
     - ``.pdb``, ``.gro``, ``.crd``, ``.mol2``, ``.smi`` (``.sdf`` planned, not yet supported)

Common flags (all trajectory subcommands)
------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Flag
     - Description
   * - ``-p / --topology``
     - Topology file
   * - ``-x / --trajectory``
     - Trajectory file
   * - ``-o / --output``
     - Output file or directory
   * - ``--overwrite``
     - Overwrite existing output
   * - ``-b / --begin``
     - First frame (default: 0)
   * - ``-e / --end``
     - Last frame (default: last)
   * - ``-s / --step``
     - Process every Nth frame (default: 1)
   * - ``-l / --log``
     - Log file path
   * - ``-fl / --file-logging-level``
     - File log verbosity (default: DEBUG)
   * - ``-tl / --terminal-logging-level``
     - Terminal log verbosity (default: INFO)

Example 1 — RMSD analysis and plot
------------------------------------

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p protein.prmtop \
       -x md.nc \
       -sel "protein and name CA" "resname LIG and not name H*" \
       -n  bb_ca ligand_heavy \
       -f  "protein and name CA" \
       -r 0 -b 0 -e 10000 -s 10 \
       -o  rmsd.pta

   pharmacon plot pta -i rmsd.pta -o plots/ --overwrite

Example 2 — Protein–ligand interactions across replicates
----------------------------------------------------------

.. code-block:: bash

   # Analyse each replicate
   for rep in rep1 rep2 rep3; do
     pharmacon trajectory pl-interactions \
         -p $rep/topol.tpr \
         -x $rep/traj.xtc \
         -prt "protein" \
         -lig "resname LIG" \
         -w  "resname WAT" \
         -o  $rep/pli.pta
   done

   # Merge replicates into one artifact
   pharmacon merge results \
       -i rep1/pli.pta rep2/pli.pta rep3/pli.pta \
       -o pli_merged.pta

   # Plot the merged artifact
   pharmacon plot pta \
       -i pli_merged.pta \
       -o plots/merged/ \
       --overwrite \
       -c my_theme.ini

Example 3 — PCA
----------------

.. code-block:: bash

   pharmacon trajectory pca \
       -p protein.prmtop \
       -x md.nc \
       -sel "protein and name CA" \
       -c 5 \
       -o pca.pta

   pharmacon plot pta -i pca.pta -o pca_plots/ --overwrite
   # Produces: pca_timeseries, pca_scatter, pca_variance_ratio,
   #           pca_probability, pca_heatmap_fes

Example 4 — Structure properties
----------------------------------

.. code-block:: bash

   pharmacon structure properties -i ligands.smi  -o ligand_props.psa
   pharmacon export psa -i ligand_props.psa -f csv -o ./results/

Example 5 — Inspect and export results
----------------------------------------

.. code-block:: bash

   pharmacon dump   pta -i rmsd.pta           # pretty-print metadata
   pharmacon export pta -i rmsd.pta -f csv -o ./results/

Debug mode
----------

Set ``PHARMACON_DEBUG=1`` to get a full traceback on unexpected errors:

.. code-block:: bash

   PHARMACON_DEBUG=1 pharmacon trajectory rmsd ...
