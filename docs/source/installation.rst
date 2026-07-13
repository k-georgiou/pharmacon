Installation
============

Pharmacon requires **Python ≥ 3.12**.

Option 1 — PyPI (recommended)
-------------------------------

The simplest way to install Pharmacon. Downloads a pre-built wheel from the
Python Package Index together with all required dependencies.  Creating a
virtual environment first is strongly recommended.

.. code-block:: bash

   python3.12 -m venv venv
   source venv/bin/activate      # Windows: venv\Scripts\activate

   pip install pharmacon

Option 2 — Conda
-----------------

Recommended if you manage binary dependencies (especially RDKit and
MDAnalysis) through Conda.  An ``environment.yml`` is provided in the
repository.

.. code-block:: bash

   conda env create -f environment.yml
   conda activate pharmacon

   pip install pharmacon

Option 3 — Source / development
---------------------------------

Clone the repository and install in editable mode.  Use this if you want to
contribute or stay on the latest development version.

.. code-block:: bash

   git clone https://github.com/k-georgiou/pharmacon.git
   cd pharmacon
   python -m venv venv && source venv/bin/activate
   pip install -e ".[dev]"

The ``dev`` extra pulls in ``pytest``.  To run the full test suite:

.. code-block:: bash

   python -m pytest

Optional: MPI-aware logging
---------------------------

Pharmacon's parallel subcommands use Python ``multiprocessing`` and require no
extra dependencies.  The optional ``mpi`` extra only installs ``mpi4py``,
which Pharmacon uses — when present — to tag per-rank log files if you launch
it under an MPI runner.  Install it only in that case:

.. code-block:: bash

   pip install "pharmacon[mpi]"

Optional: development tools only (no editable install)
-------------------------------------------------------

If you already have Pharmacon installed and only need the test dependencies:

.. code-block:: bash

   pip install "pharmacon[dev]"

Verifying the installation
--------------------------

.. code-block:: bash

   pharmacon --help

You should see the Pharmacon banner followed by the top-level help listing
all available commands.
