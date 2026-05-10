Command Reference
=================

Pharmacon exposes six top-level commands, each grouping related subcommands.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Description
   * - :doc:`trajectory`
     - Analyze MD trajectories: RMSD, distances, angles, hydrogen bonds, protein–ligand and protein–protein interactions, PCA, average structure
   * - :doc:`structure`
     - Analyze static structures: molecular properties and amino-acid sequence extraction
   * - :doc:`dump`
     - Print metadata and group tree of ``.pta`` / ``.psa`` files to stdout
   * - :doc:`export`
     - Export ``.pta`` / ``.psa`` data to CSV, TSV, or FASTA
   * - :doc:`plot`
     - Render publication-ready plots from ``.pta`` files using optional INI configuration
   * - :doc:`merge`
     - Merge multiple ``.pta`` replicates from the same topology into one aggregated artifact

General syntax
--------------

.. code-block:: text

   pharmacon <command> <subcommand> [options]

.. code-block:: bash

   pharmacon --help                          # top-level help
   pharmacon trajectory --help               # trajectory subcommands
   pharmacon trajectory rmsd --help          # RMSD flags

Adding a new subcommand
-----------------------

Pharmacon uses auto-discovery.  Drop a ``.py`` file into the appropriate
``src/pharmacon/command_line/<command>/`` directory.  The module must export:

.. code-block:: python

   SUBCOMMAND_NAME: str          # e.g. "my-analysis"
   SUMMARY: str                  # one-line description

   def build_parser(subparsers, parents): ...
   def validate(args): ...
   def run(args): ...

No registration step is required — the registry discovers it at import time.
See :doc:`/architecture` for the full development pattern.
