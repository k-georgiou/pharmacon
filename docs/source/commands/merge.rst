merge
=====

Combine multiple ``.pta`` files produced from the same topology into a single
aggregated artifact.  The standard use case is pooling independent replicate
runs before plotting or exporting.

Before merging, Pharmacon requires that all inputs share the same **command**
and **subcommand** — this is a hard requirement and mismatched files are
refused.  It also compares each file's **blueprint** (a deterministic hash of
the run inputs — topology, selections, flags); a differing blueprint (e.g. a
different topology) is reported as a warning and, under the default
zero-warning budget (``-mw 0``), aborts the merge.

Merged files replace per-frame datasets with aggregate summary tables stored
under ``modes_merged/``.  The ``is_merged`` root attribute is set to ``True``.

.. contents:: Subcommands
   :local:
   :depth: 1

----

merge results
-------------

Merge two or more ``.pta`` files into one aggregated artifact.

**Arguments**

.. list-table::
   :header-rows: 1
   :widths: 35 10 55

   * - Flag
     - Required
     - Description
   * - ``-i / --input``
     - Yes
     - Two or more input ``.pta`` files (space-separated)
   * - ``-o / --output``
     - Yes
     - Output merged ``.pta`` file
   * - ``--overwrite``
     - No
     - Overwrite an existing output file
   * - ``-mw / --max-warnings``
     - No
     - Maximum number of coercion warnings to tolerate
   * - ``-l / --log``
     - No
     - Log file (default: ``merge_results.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Merge two RMSD files from consecutive trajectory segments:

.. code-block:: bash

   pharmacon merge results \
       -i seg1/rmsd.pta seg2/rmsd.pta \
       -o merged_rmsd.pta

Merge three replicate interaction files:

.. code-block:: bash

   pharmacon merge results \
       -i rep1/pli.pta rep2/pli.pta rep3/pli.pta \
       -o pli_merged.pta

Merge and overwrite an existing output:

.. code-block:: bash

   pharmacon merge results \
       -i rep1/pli.pta rep2/pli.pta rep3/pli.pta \
       -o pli_merged.pta \
       --overwrite

Typical replicate workflow (analyse → merge → plot):

.. code-block:: bash

   # Step 1 — run each replicate
   for rep in rep1 rep2 rep3; do
     pharmacon trajectory pl-interactions \
         -p $rep/topol.tpr -x $rep/traj.xtc \
         -prt "protein" -lig "resname LIG" \
         -w  "resname WAT" \
         -o  $rep/pli.pta
   done

   # Step 2 — merge
   pharmacon merge results \
       -i rep1/pli.pta rep2/pli.pta rep3/pli.pta \
       -o pli_merged.pta

   # Step 3 — plot
   pharmacon plot pta \
       -i pli_merged.pta \
       -o plots/merged/ \
       --overwrite
