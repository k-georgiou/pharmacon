plot
====

Render publication-ready plots directly from Pharmacon artifacts.  The plotter
auto-discovers which analysis groups are present in the file and dispatches to
the appropriate rendering functions.  Plot appearance is fully customisable via
INI configuration files — see :doc:`/plot_configuration` for the full
reference.

.. contents:: Subcommands
   :local:
   :depth: 1

----

plot pta
--------

Generate all compatible plots for every group found in a ``.pta`` file.
The type of figure produced depends on the analysis stored in the file:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Analysis in file
     - Plots produced
   * - ``rmsd``, ``distances``, ``angles``
     - Time-series line plot (``PTA-UNIFIED``)
   * - ``rmsf``
     - Per-atom RMSF profile
   * - ``pl_interactions``
     - Stacked column (×2), heatmap (×2), pie charts, ligand monitor heatmap
   * - ``pp_interactions``
     - Timeline heatmap, residue heatmap, stacked column
   * - ``hbonds``
     - Residue×residue heatmap, occupancy ranking, H-bond network graph, timeline, count-per-frame
   * - ``pca``
     - Time series, scatter, variance-ratio scree, FES heatmap, probability heatmap

Plots that require per-frame data are automatically skipped on merged files
(``is_merged=True``).

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
   * - ``-o / --output``
     - Yes
     - Output directory where plot files will be written
   * - ``-c / --config``
     - No
     - Pharmacon INI file with plot-settings overrides (see :doc:`/plot_configuration`)
   * - ``--overwrite``
     - No
     - Overwrite existing output directory contents
   * - ``-mw / --maxwarnings``
     - No
     - Maximum coercion warnings per plot before the plot is skipped (default: 0)
   * - ``-l / --log``
     - No
     - Log file (default: ``plot_pta.log``)
   * - ``-fl / --file-logging-level``
     - No
     - File log verbosity (default: ``DEBUG``)
   * - ``-tl / --terminal-logging-level``
     - No
     - Terminal log verbosity (default: ``INFO``)

**Examples**

Plot RMSD results:

.. code-block:: bash

   pharmacon plot pta \
       -i rmsd.pta \
       -o ./plots/

Plot protein–ligand interactions:

.. code-block:: bash

   pharmacon plot pta \
       -i pl_interactions.pta \
       -o ./plots/

Plot a merged artifact with a custom INI theme:

.. code-block:: bash

   pharmacon plot pta \
       -i pli_merged.pta \
       -o ./plots/ \
       -c my_theme.ini \
       --overwrite

Use the master INI that covers every plot type:

.. code-block:: bash

   pharmacon plot pta \
       -i run.pta \
       -o ./plots/ \
       -c examples/plot_ini/all_plots.ini

Generate all PCA plots:

.. code-block:: bash

   pharmacon trajectory pca -p topol.tpr -x traj.xtc -sel "protein and name CA" -c 5 -o pca.pta
   pharmacon plot pta -i pca.pta -o pca_plots/ --overwrite
   # Produces: pca_timeseries, pca_scatter, pca_variance_ratio,
   #            pca_probability, pca_heatmap_fes
