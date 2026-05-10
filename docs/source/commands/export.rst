export
======

Convert Pharmacon HDF5 artifacts into human-readable flat files for use in
spreadsheets, R, pandas, or bioinformatics pipelines.

.. contents:: Subcommands
   :local:
   :depth: 1

----

export pta
----------

Export a Pharmacon Trajectory Analysis (``.pta``) file to CSV or TSV.  All
datasets inside the file are written as separate files into the specified
output directory.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-i / --input``
     - Yes
     - Input ``.pta`` file
   * - ``-f / --format``
     - Yes
     - Output format: ``csv`` or ``tsv``
   * - ``-o / --output``
     - Yes
     - Output directory where exported files will be written
   * - ``--overwrite``
     - No
     - Overwrite existing output files
   * - ``-mw / --max-warnings``
     - No
     - Maximum number of coercion warnings to tolerate
   * - ``-l / --log``
     - No
     - Log file (default: ``export_pta.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Export RMSD results to CSV:

.. code-block:: bash

   pharmacon export pta \
       -i rmsd.pta \
       -f csv \
       -o ./results/

Export interaction data to TSV:

.. code-block:: bash

   pharmacon export pta \
       -i pl_interactions.pta \
       -f tsv \
       -o ./exports/

Export a merged artifact to CSV and overwrite existing output:

.. code-block:: bash

   pharmacon export pta \
       -i pli_merged.pta \
       -f csv \
       -o ./results/ \
       --overwrite

Quick inspect then export workflow:

.. code-block:: bash

   pharmacon dump   pta -i rmsd.pta
   pharmacon export pta -i rmsd.pta -f csv -o ./results/

----

export psa
----------

Export a Pharmacon Structure Analysis (``.psa``) file to CSV, TSV, or FASTA.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-i / --input``
     - Yes
     - Input ``.psa`` file
   * - ``-f / --format``
     - Yes
     - Output format: ``csv``, ``tsv``, or ``fasta``
   * - ``-o / --output``
     - Yes
     - Output directory where exported files will be written
   * - ``--overwrite``
     - No
     - Overwrite existing output files
   * - ``-l / --log``
     - No
     - Log file (default: ``export_psa.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Export molecular properties to CSV:

.. code-block:: bash

   pharmacon export psa \
       -i properties.psa \
       -f csv \
       -o ./results/

Export amino-acid sequence to FASTA (for downstream bioinformatics):

.. code-block:: bash

   pharmacon export psa \
       -i sequence.psa \
       -f fasta \
       -o ./results/

Export to TSV and overwrite existing output:

.. code-block:: bash

   pharmacon export psa \
       -i properties.psa \
       -f tsv \
       -o ./results/ \
       --overwrite
