dump
====

Print metadata and the internal group/dataset tree of Pharmacon HDF5 artifacts
to the terminal as rich formatted panels.  Useful for verifying a file is
intact, checking which analyses it contains, and inspecting run parameters
without opening an HDF5 viewer or writing any Python code.

.. contents:: Subcommands
   :local:
   :depth: 1

----

dump pta
--------

Inspect a Pharmacon Trajectory Analysis (``.pta``) file.

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
   * - ``--disable-file-metadata``
     - No
     - Suppress the file-level metadata panel
   * - ``--disable-group-meta``
     - No
     - Suppress group-level metadata
   * - ``--disable-dataset-meta``
     - No
     - Suppress dataset-level metadata
   * - ``--verbose``
     - No
     - Print extended info including dataset shapes and dtypes

**Examples**

Basic inspection showing metadata, group tree, and datasets:

.. code-block:: bash

   pharmacon dump pta -i rmsd.pta

Print everything in verbose mode (shapes, dtypes):

.. code-block:: bash

   pharmacon dump pta -i rmsd.pta --verbose

Show only file-level metadata (hide groups and datasets):

.. code-block:: bash

   pharmacon dump pta -i rmsd.pta \
       --disable-group-meta \
       --disable-dataset-meta

Inspect a merged artifact:

.. code-block:: bash

   pharmacon dump pta -i pli_merged.pta --verbose

----

dump psa
--------

Inspect a Pharmacon Structure Analysis (``.psa``) file.

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
   * - ``--disable-file-metadata``
     - No
     - Suppress the file-level metadata panel
   * - ``--disable-group-meta``
     - No
     - Suppress group-level metadata
   * - ``--disable-dataset-meta``
     - No
     - Suppress dataset-level metadata
   * - ``--verbose``
     - No
     - Print extended info including dataset shapes and dtypes

**Examples**

Basic inspection of a properties file:

.. code-block:: bash

   pharmacon dump psa -i properties.psa

Verbose inspection of a sequence file:

.. code-block:: bash

   pharmacon dump psa -i sequence.psa --verbose

Show only file-level metadata:

.. code-block:: bash

   pharmacon dump psa -i properties.psa \
       --disable-group-meta \
       --disable-dataset-meta
