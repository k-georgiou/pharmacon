Logging & Workspaces
=====================

Logging
-------

Every Pharmacon subcommand except ``dump`` and ``export`` (which have no
logging flags) accepts three logging flags:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Flag
     - Description
   * - ``-l / --log FILE``
     - Path to the log file (default varies per subcommand, e.g. ``rmsd.log``)
   * - ``-fl / --file-logging-level LEVEL``
     - Verbosity written to the log file (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level LEVEL``
     - Verbosity printed to the terminal (default: ``INFO``)

**Available log levels** (from most to least verbose):

.. code-block:: text

   TRACE  →  DEBUG  →  INFO  →  WARNING  →  ERROR  →  CRITICAL

``TRACE`` is a Pharmacon-specific level below ``DEBUG`` used for very detailed
per-frame output (e.g. logging each processed frame).

**Example** — suppress terminal output except errors, keep a full DEBUG log:

.. code-block:: bash

   pharmacon trajectory rmsd \
       -p topol.tpr -x traj.xtc \
       -sel "protein and name CA" -f "protein and name CA" -n ca \
       -o rmsd.pta \
       -l  rmsd.log \
       -fl DEBUG \
       -tl ERROR

Debug mode
----------

Set the ``PHARMACON_DEBUG`` environment variable to get a full Python
traceback on unexpected errors instead of the default one-liner:

.. code-block:: bash

   PHARMACON_DEBUG=1 pharmacon trajectory rmsd ...

----

Workspaces
----------

``PharmaconWorkspace`` manages the working directory and an optional temporary
directory for each run.  Two environment variables control its behaviour.

PHARMACON_WORKDIR
^^^^^^^^^^^^^^^^^

Sets the **working directory** used by Pharmacon for the run.

**Resolution order:**

1. Value of ``PHARMACON_WORKDIR`` environment variable
2. ``.`` — the shell's current working directory (default)

The directory is **created automatically** (``mkdir -p``) if it does not
already exist.

.. code-block:: bash

   # Run with a dedicated project directory
   export PHARMACON_WORKDIR=/data/my_project
   pharmacon trajectory rmsd -p topol.tpr -x traj.xtc ...

   # Or inline for a single command
   PHARMACON_WORKDIR=/data/my_project pharmacon trajectory rmsd ...

PHARMACON_TEMPDIR
^^^^^^^^^^^^^^^^^

Sets the **base directory** inside which Pharmacon creates a per-run scratch
directory when one is needed.

**Resolution order:**

1. Value of ``PHARMACON_TEMPDIR`` environment variable
2. ``None`` — the OS default temporary location (e.g. ``/tmp`` on Linux) is
   used by ``tempfile.mkdtemp``

**Important:** ``PHARMACON_TEMPDIR`` is the *parent* for the scratch space, not
the scratch space itself.  The actual directory created is a uniquely named
subdirectory:

.. code-block:: text

   $PHARMACON_TEMPDIR/pharmacon-<random>/

The base directory is **created automatically** if it does not exist.

.. note::

   Not every subcommand creates a temporary directory.  The scratch space is
   only allocated when the specific analysis requires it
   (``is_tmp_dir_needed=True`` internally).  Most standard subcommands
   (``rmsd``, ``distances``, ``pca``, etc.) do not use one.

**Automatic cleanup:**

The temporary directory is deleted automatically via ``atexit`` when the
Pharmacon process exits — including on unhandled exceptions and
``KeyboardInterrupt``.  Only a hard kill (``SIGKILL``) bypasses cleanup.

.. code-block:: bash

   # Place all scratch space under /scratch instead of /tmp
   export PHARMACON_TEMPDIR=/scratch/pharmacon_tmp
   pharmacon trajectory pl-interactions -p topol.tpr -x traj.xtc ...
   # Creates: /scratch/pharmacon_tmp/pharmacon-a3f9c2b1/

Environment variable summary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Variable
     - Default
     - Effect
   * - ``PHARMACON_WORKDIR``
     - ``.`` (current directory)
     - Working directory for the run; created automatically if absent
   * - ``PHARMACON_TEMPDIR``
     - OS default (``/tmp``)
     - Parent directory for per-run scratch space; a ``pharmacon-<random>/``
       subdirectory is created inside it; auto-cleaned on exit
   * - ``PHARMACON_DEBUG``
     - unset
     - Print full Python tracebacks on unexpected errors (same effect as the ``--debug`` flag)
