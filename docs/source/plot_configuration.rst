Plot Configuration (INI Files)
==============================

Plot appearance in Pharmacon is fully customisable through
`configobj <https://configobj.readthedocs.io/>`_-style INI files.
Pass any INI file to ``pharmacon plot pta`` via the ``-c / --config`` flag.

.. code-block:: bash

   pharmacon plot pta -i run.pta -o ./plots/ -c my_settings.ini

How it works
------------

Each INI section header must match the **alias** of a plot-settings dataclass.
Sections not recognised by the invoked subcommand are silently ignored, so a
single master INI file can safely be reused across different ``plot`` commands.

Unknown field values are **coerced** with a warning.  The ``-mw / --max-warnings``
flag controls how many coercion warnings a plot may accumulate before it is
skipped entirely (default: 0).

Example INI
-----------

.. code-block:: ini

   # ── Shared figure theme ─────────────────────────────────────────────────
   [PTA-UNIFIED]
   fig_dpi          = 300
   fig_format       = svg
   font_family      = DejaVu Sans
   enable_grid      = true
   grid_style       = dashed
   grid_alpha       = 0.4

   # ── Per-plot overrides ───────────────────────────────────────────────────
   [PLI-STACKED-COLUMN-1]
   bar_width         = 0.9
   bar_alpha         = 0.85
   color_hydrophobic = "#b39ddb"

   [PLI-HEATMAP-2]
   cmap              = viridis
   threshold         = 0.05

   [PCA-SCATTER]
   cmap              = plasma
   scatter_size      = 14
   disable_colorbar  = false

Available plot sections and aliases
-------------------------------------

Ready-to-copy template INI files for every plot type are in
``examples/plot_ini/``.

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - INI template file
     - Plot type
     - Section alias
   * - ``all_plots.ini``
     - Master file (all plots)
     - *(all aliases)*
   * - ``pta_unified.ini``
     - Universal PTA time series (RMSD, distances, angles, …)
     - ``PTA-UNIFIED``
   * - ``pli_stacked_column_1.ini``
     - PLI stacked column (per interaction type)
     - ``PLI-STACKED-COLUMN-1``
   * - ``pli_stacked_column_2.ini``
     - PLI stacked column (backbone vs side-chain)
     - ``PLI-STACKED-COLUMN-2``
   * - ``pli_heatmap_1.ini``
     - PLI heatmap (residue × frame)
     - ``PLI-HEATMAP-1``
   * - ``pli_heatmap_2.ini``
     - PLI heatmap (interaction × frame)
     - ``PLI-HEATMAP-2``
   * - ``pli_pie_charts_1.ini``
     - PLI pie charts (collage)
     - ``PLI-PIE-CHARTS-1``
   * - ``pli_ligand_monitor.ini``
     - Residue × ligand-atom contact heatmap
     - ``PLI-LIGAND-MONITOR``
   * - ``ppi_timeline_pairs.ini``
     - PPI timeline heatmap for top residue pairs
     - ``PPI-TIMELINE-PAIRS``
   * - ``ppi_heatmap.ini``
     - PPI residue × residue heatmap
     - ``PPI-HEATMAP``
   * - ``ppi_stacked_column.ini``
     - PPI per-pair stacked column
     - ``PPI-STACKED-COLUMN``
   * - ``pca_timeseries.ini``
     - PCA projections vs time
     - ``PCA-TIMESERIES``
   * - ``pca_scatter.ini``
     - PCA scatter of two components
     - ``PCA-SCATTER``
   * - ``pca_variance_ratio.ini``
     - PCA explained-variance scree plot
     - ``PCA-VARIANCE-RATIO``
   * - ``pca_fes_heatmap.ini``
     - PCA Free-Energy Surface heatmap
     - ``PCA-FES-HEATMAP``
   * - ``pca_probability_heatmap.ini``
     - PCA probability density heatmap
     - ``PCA-PROBABILITY-HEATMAP``

Using template files
---------------------

Point any ``pharmacon plot`` command at a template directly:

.. code-block:: bash

   # Use the PTA unified settings for an RMSD plot
   pharmacon plot pta -i rmsd.pta -c examples/plot_ini/pta_unified.ini -o ./plots/

   # Use the master file for PLI plots
   pharmacon plot pta -i pl_interactions.pta -c examples/plot_ini/all_plots.ini -o ./plots/

   # The master file is safe to use for any plot type —
   # unrecognised sections are silently ignored
   pharmacon plot pta -i pca.pta -c examples/plot_ini/all_plots.ini -o ./pca_plots/
