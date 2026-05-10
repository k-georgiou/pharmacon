Architecture
============

Pharmacon is a pure-Python command-line package centred on a single
executable ``pharmacon``.  It follows a **git-style** command/subcommand
structure and is organised around the following core modules.

Module overview
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Module
     - Responsibility
   * - **CLI Dispatcher** (``pharmacon.command_line``)
     - Routes ``pharmacon <command> <subcommand>`` to the correct analysis,
       plotting, or export module.  Commands and subcommands are discovered
       automatically at import time — no registration required.
   * - **Analyzer** (``pharmacon.analyzer``)
     - Computational core.  Stateless routines that operate on MDAnalysis
       selections and return raw numerical results.
       See :doc:`analyzer` for details.
   * - **File I/O** (``pharmacon.fileio``)
     - Reads and writes signed, self-describing HDF5 artifacts (``.pta`` /
       ``.psa``).  See :doc:`file_formats` for the full schema.
   * - **Plotter** (``pharmacon.plotter``)
     - Generates publication-ready figures (heatmaps, stacked columns, pie
       charts, scatter plots, time series) directly from analysis artifacts.
   * - **Plot Settings** (``pharmacon.constants.plots``)
     - Fully customisable settings system driven by INI configuration files.
       See :doc:`plot_configuration`.
   * - **Logger** (``pharmacon.logger``)
     - Structured terminal and file logging with multiple verbosity levels
       and Rich-coloured output.  See :doc:`logging`.
   * - **Workspace** (``pharmacon.utils.workspace``)
     - Manages working and temporary directories with automatic cleanup on
       exit.  Configurable via environment variables.
   * - **Fingerprint** (``pharmacon.utils.fingerprint``)
     - Generates deterministic identity tokens (signature, fingerprint,
       blueprint) and tamper-detection checksums for every output file.
   * - **Merge** (``pharmacon.command_line.merge``)
     - Verifies blueprint compatibility across replicates and combines
       ``.pta`` files into one aggregated artifact.

Adding a new subcommand
-----------------------

Follow these four steps:

1. Create ``src/pharmacon/command_line/<group>/<name>.py`` exporting
   ``SUBCOMMAND_NAME``, ``SUMMARY``, ``build_parser``, ``validate``, ``run``.
2. Put the computational logic in ``pharmacon.analyzer.<name>`` so it remains
   unit-testable independently of the CLI.
3. If the output is a new HDF5 artifact, extend ``pharmacon.fileio``.
4. If it needs plots, add a settings dataclass to
   ``pharmacon.constants.plots`` and a rendering function to
   ``pharmacon.plotter``.

The registry auto-discovers the new file on next import — no other changes
are needed.

Dependencies
------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Library
     - Role
   * - MDAnalysis ≥ 2.7
     - Topology/trajectory parsing and atom selection
   * - RDKit ≥ 2024.03.1
     - SMARTS-based atom-type detection, molecular fingerprints and properties
   * - NumPy ≥ 1.26
     - Numerical arrays and linear algebra
   * - SciPy ≥ 1.11
     - Scientific algorithms (PCA, statistics)
   * - Matplotlib ≥ 3.8
     - Figure rendering
   * - h5py ≥ 3.10
     - HDF5 file I/O
   * - NetworkX ≥ 3.2
     - Graph-based analysis (interaction networks)
   * - Rich ≥ 13.0
     - Terminal output formatting
   * - configobj ≥ 5.0
     - INI plot-settings parsing
   * - pyarrow ≥ 14.0
     - Parquet fingerprint export
   * - periodictable ≥ 1.7
     - Element data
   * - psutil ≥ 5.9
     - Process and memory monitoring
