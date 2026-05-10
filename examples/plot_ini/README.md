# Pharmacon — Examples

This directory contains ready-to-use example configuration files for the
Pharmacon CLI.

## `plot_ini/`

Plot configuration files written in INI format (parsed by `configobj`).
Every supported plot type has two representations:

| File                          | Purpose                                                                   |
|-------------------------------|---------------------------------------------------------------------------|
| `all_plots.ini`               | Single master file containing one section per plot type (minimal fields). |
| `pli_stacked_column_1.ini`    | PLI stacked column plot (per-residue × interaction type)                  |
| `pli_stacked_column_2.ini`    | PLI stacked column plot (backbone vs side-chain)                          |
| `pli_heatmap_1.ini`           | PLI heatmap (residue × frame)                                             |
| `pli_heatmap_2.ini`           | PLI heatmap (interaction × frame, auto-sized)                             |
| `pli_pie_charts_1.ini`        | Per-residue pie charts + collage                                          |
| `pli_ligand_monitor.ini`      | Residue × ligand-atom contact heatmap                                     |
| `ppi_timeline_pairs.ini`      | PPI timeline heatmap for the top residue-pairs                            |
| `ppi_heatmap.ini`             | PPI residue × residue contact heatmap                                     |
| `ppi_stacked_column.ini`      | PPI per-pair stacked column                                               |
| `pta_unified.ini`             | Universal PTA time series (RMSD/Hbonds/…)                                 |
| `pca_timeseries.ini`          | PCA projections vs time                                                   |
| `pca_scatter.ini`             | PCA scatter of two components                                             |
| `pca_variance_ratio.ini`      | PCA explained-variance scree plot                                         |
| `pca_fes_heatmap.ini`         | PCA Free-Energy Surface heatmap                                           |
| `pca_probability_heatmap.ini` | PCA probability density heatmap                                           |

Individual files list **every variable** Pharmacon exposes for each plot,
together with their default values and inline hints. Copy any one of them,
edit, and pass it to `pharmacon plot ...` via `-c`.

## Usage

```bash
# Use a single-plot config
pharmacon plot pta -i run.pta -c examples/plot_ini/pta_unified.ini

# Or a master config containing every section
pharmacon plot pta -i run.pta -c examples/plot_ini/all_plots.ini
pharmacon plot pca -i run.pta -c examples/plot_ini/all_plots.ini
```

## Rules / Tips

- **Section headers are case-insensitive** and accept any of the aliases listed
  at the top of each file (e.g. `[PLI-HEATMAP-1]`, `[PLI_HEATMAP_1]`, and
  `[PL-INTERACTIONS-HEATMAP-1]` are equivalent).
- **Unknown keys are silently ignored**, so a master file with every section
  can safely be reused across different `plot` subcommands.
- **Empty values mean "use default / auto"**. For instance, leaving
  `vmax =` blank lets the plot compute its own maximum.
- **List-valued fields** (e.g. `line_colors`, `pcs`, `components`) accept a
  comma-separated list.
- **Invalid values are coerced back to the default** and a warning is emitted
  on stderr.
