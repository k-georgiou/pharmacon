<div align="center">

# Pharmacon — Plot INI Reference

</div>

This directory holds ready-to-edit plot configuration files for the
`pharmacon plot ...` subcommands. Every plot Pharmacon ships is controlled
by an INI file — copy one of the examples below, tweak, and pass it via
`-c`:

```bash
pharmacon plot pta -i run.pta -o ./plots -c examples/plot_ini/rmsf.ini
pharmacon plot pta -i merged.pta -o ./plots -c examples/plot_ini/all_plots.ini
```

> [!TIP]
> If you don't want to manage many small files, **`all_plots.ini`** is a
> single master file with one section per plot type — pass it to any `plot`
> subcommand and only the relevant sections are read; the rest are silently
> ignored.

---

## 1. Quick reference — what's in this directory

| File                          | Plot family     | Primary alias             | What it produces                                                 |
|-------------------------------|-----------------|---------------------------|------------------------------------------------------------------|
| `all_plots.ini`               | (master)        | (every alias)             | One section per plot type — pass to any `plot` subcommand        |
| `pta_unified.ini`             | PTA time-series | `PTA-UNIFIED`             | RMSD / angles / distances time-series                            |
| `rmsf.ini`                    | RMSF            | `PTA-RMSF`                | Per-atom RMSF profile, one curve per selection                   |
| `rmsf_ligand.ini`             | RMSF            | `PTA-RMSF`                | Per-atom RMSF tuned for a ligand (atom_index + atom_name labels) |
| `pca_timeseries.ini`          | PCA             | `PCA-TIMESERIES`          | PC projections vs time                                           |
| `pca_scatter.ini`             | PCA             | `PCA-SCATTER`             | 2-component PC scatter                                           |
| `pca_variance_ratio.ini`      | PCA             | `PCA-VARIANCE-RATIO`      | Scree (explained-variance) plot                                  |
| `pca_fes_heatmap.ini`         | PCA             | `PCA-FES-HEATMAP`         | Free-energy surface −kT·ln P                                     |
| `pca_probability_heatmap.ini` | PCA             | `PCA-PROBABILITY-HEATMAP` | 2D probability density map                                       |
| `pli_stacked_column_1.ini`    | PLI             | `PLI-STACKED-COLUMN-1`    | Per-residue × interaction-type stacked bar                       |
| `pli_stacked_column_2.ini`    | PLI             | `PLI-STACKED-COLUMN-2`    | Per-residue backbone vs side-chain stacked bar                   |
| `pli_heatmap_1.ini`           | PLI             | `PLI-HEATMAP-1`           | Residue × frame heatmap                                          |
| `pli_heatmap_2.ini`           | PLI             | `PLI-HEATMAP-2`           | Interaction × frame auto-sized heatmap                           |
| `pli_pie_charts_1.ini`        | PLI             | `PLI-PIE-CHARTS-1`        | Per-residue pie charts + optional collage                        |
| `pli_ligand_monitor.ini`      | PLI             | `PLI-LIGAND-MONITOR`      | Residue × ligand-atom contact heatmap                            |
| `ppi_timeline_pairs.ini`      | PPI             | `PPI-TIMELINE-PAIRS`      | Top residue-pair timeline heatmap                                |
| `ppi_heatmap.ini`             | PPI             | `PPI-HEATMAP`             | Residue × residue contact frequency heatmap                      |
| `ppi_stacked_column.ini`      | PPI             | `PPI-STACKED-COLUMN`      | Per-pair stacked column                                          |

---

## 2. How INI files work

### Sections and aliases

Each plot type has one or more **aliases**. A section header in the INI
file uses any of them — case-insensitive, hyphens or underscores both
accepted:

```ini
[PTA-UNIFIED]    ; same plot
[PTA_UNIFIED]    ; same plot
[pta-unified]    ; same plot
```

Every plot family's aliases are listed in the per-plot reference below.

### Unknown keys are silently ignored

The parser only consumes keys it recognises for a given plot class.
This is *deliberate* — it's what lets one master file (`all_plots.ini`)
hold sections for every plot type without errors when only one
subcommand reads it. A small downside: a typo in a key name is silent.
If a plot looks suspiciously default-styled, double-check that you spelt
the keys correctly.

### Unknown sections produce a warning

If a `[SECTION]` doesn't match any registered alias, you'll see a
console warning at the end of the run listing the unrecognised sections.

### Values: scalars, lists, and "auto"

- Numeric, string, and boolean values are parsed as you'd expect
  (`True`/`False`, `true`/`false`, `1`/`0`).
- **Lists** are comma-separated:
  ```ini
  line_colors = #1f77b4, #ff7f0e, #2ca02c
  pcs         = 1, 2, 3
  ```
- **Empty** values mean "use default" (e.g. `vmax =` lets the plot auto-pick).
- **`auto`** is accepted for axis-limit fields (`x_min`, `x_max`,
  `y_min`, `y_max` on the RMSF plot).
- **Hex colours must be quoted.** `#` starts an inline comment in INI
  syntax, so an unquoted hex value gets stripped to an empty string and
  the plot silently falls back to its default colour. Always write:
  ```ini
  color_hydrophobic = "#b39ddb"
  line_colors       = "#1f77b4", "#ff7f0e", "#2ca02c"
  colors_by_label   = calpha:"#d62728", backbone:"#1f77b4"
  ```
  Named colours (`red`, `dodgerblue`) and rgb tuples don't need quoting.

### Validation & coercion

Pharmacon's settings classes validate every field. Invalid values
(out-of-range numbers, mistyped enums, malformed colours, bad
templates) **do not raise** — they emit a console warning and fall back
to the documented default. This lets a long INI file with one bad value
still render. Each subcommand also has a `--maxwarnings N` flag that
will *skip* a plot whose settings produced more than N coercion warnings.

### Editable text in SVG / PDF

Both plotter modules force matplotlib's rcParams to `svg.fonttype="none"`
and `pdf.fonttype=42` at import time, so **every SVG/PDF Pharmacon writes
has real `<text>` elements / TrueType fonts** — open them in Illustrator
or Inkscape and you can retype labels, recolour them, restyle them
without redrawing the figure. No INI switch needed.

### Universal output settings — fig_format and dpi

| Key                                 | Accepted values                                 | Notes                                                                                      |
|-------------------------------------|-------------------------------------------------|--------------------------------------------------------------------------------------------|
| `fig_format`                        | `png`, `jpg`/`jpeg`, `svg`, `pdf`, `tif`/`tiff` | Use `svg` or `pdf` for publication / Illustrator                                           |
| `fig_dpi`                           | 50 – 2000                                       | 300 for documents, 600 for posters, 150 for slides                                         |
| `fig_size_width`, `fig_size_height` | inches (typically 6–14)                         | matplotlib's natural unit                                                                  |
| `fig_transparent`                   | bool                                            | `True` strips the white background — useful when pasting onto coloured slides              |
| `tight_layout`                      | bool                                            | Auto-fit ticks/labels. Off by default for some stacked plots that hand-position the legend |
| `bg_color`                          | colour                                          | Default `white`                                                                            |

---

## 3. Variables shared by *most* plots

These appear across nearly every settings class. Per-plot tables below
only list the *additional* / *plot-specific* variables — assume the
shared block is also available.

### 3.1 Figure
| Key                                  | Default | Notes                                            |
|--------------------------------------|---------|--------------------------------------------------|
| `fig_size_width` / `fig_size_height` | varies  | inches                                           |
| `fig_dpi`                            | 300–800 | clamped to 50–2000                               |
| `fig_basename`                       | varies  | output stem (extension comes from `fig_format`)  |
| `fig_format`                         | `png`   | png / jpg / svg / pdf / tif                      |
| `fig_transparent`                    | `False` | strip white background                           |
| `tight_layout`                       | varies  | auto-fit margins                                 |
| `bg_color`                           | `white` | any colour string                                |

### 3.2 Title and axis labels
| Key                                  | Default | Notes                               |
|--------------------------------------|---------|-------------------------------------|
| `fig_title`                          | varies  | empty string = auto / generic title |
| `x_label`, `y_label`                 | varies  | leave empty for auto                |
| `disable_title`                      | `False` | omit the title entirely             |
| `disable_x_label`, `disable_y_label` | `False` | omit axis labels                    |
| `disable_x_axis`, `disable_y_axis`   | `False` | omit the whole axis spine + ticks   |
| `disable_ticks`                      | `False` | omit ticks but keep the axis        |

### 3.3 Fonts
| Key                                             | Default       | Notes                                                                      |
|-------------------------------------------------|---------------|----------------------------------------------------------------------------|
| `font_family`                                   | `dejavu sans` | Any matplotlib-available font name; falls back to `dejavu sans` if unknown |
| `font_size_title`                               | 12–18         |                                                                            |
| `font_size_label`                               | 10–16         | axis labels                                                                |
| `font_size_ticks`                               | 8–10          | tick numbers                                                               |
| `font_size_legend`                              | 8–10          |                                                                            |
| `font_size_x`, `font_size_y`                    | 8–10          | per-axis tick fonts (heatmaps)                                             |
| `font_size_cbar`                                | 8–12          | colourbar tick font                                                        |
| `font_weight_title`                             | `bold`        | one of: `normal`, `bold`, `light`, `medium`, `heavy`, `semibold`, `book` … |
| `font_weight_label`, `font_weight_legend`, etc. | `normal`      |                                                                            |

### 3.4 Grid
| Key              | Default             | Notes                                                                |
|------------------|---------------------|----------------------------------------------------------------------|
| `enable_grid`    | varies              |                                                                      |
| `grid_style`     | `dashed`            | one of: `solid`, `dashed`, `dashdot`, `dotted`, `-`, `--`, `-.`, `:` |
| `grid_alpha`     | 0.3–0.5             | 0=invisible, 1=opaque                                                |
| `grid_color`     | `lightgray`/`black` | heatmaps only                                                        |
| `grid_linewidth` | 0.2                 | heatmaps only                                                        |

### 3.5 Legend
| Key                    | Default | Notes                                                  |
|------------------------|---------|--------------------------------------------------------|
| `disable_legend`       | `False` |                                                        |
| `legend_loc`           | `best`  | matplotlib loc strings: `upper right`, `lower left`, … |
| `legend_frame`         | `True`  | draw a frame around the legend                         |
| `legend_alpha`         | 1.0     | frame transparency                                     |
| `legend_n_col`         | varies  | columns; useful for stacked-column legends             |
| `legend_bbox_y`        | varies  | vertical offset for off-axis legends                   |
| `legend_margin_bottom` | varies  | adds bottom margin so an off-axis legend isn't clipped |

### 3.6 Residue-label rewriting (PLI / PPI / PPI-heatmap)

Several plots show residue identifiers. You can rewrite them via:

| Key                                                                     | Default                | Notes                                                                                                |
|-------------------------------------------------------------------------|------------------------|------------------------------------------------------------------------------------------------------|
| `representation` *or* `x_axis_representation` / `y_axis_representation` | `chainid:resnameresid` | Pattern combining `resname`, `resid`, `chainid`, `segid` separated by `:`, `-`, `_`, or no separator |
| `aa3_to_aa1`                                                            | varies                 | `ALA` → `A`, `LYS` → `K`, …                                                                          |
| `renumber`                                                              | `False`                | Renumber residues starting at `renumber_int` (or 0)                                                  |
| `renumber_int`                                                          | `0`                    | start integer for renumbering                                                                        |
| `alter_chains`                                                          | `False`                | rename chains; format defined in `alter_chains_str` (`"A>X, B>Y"`)                                   |
| `alter_segments`                                                        | `False`                | rename segments via `alter_segments_str`                                                             |

Representation pattern examples:
- `chainid:resnameresid` → `A:ARG52`
- `resnameresid` → `ARG52`
- `resid-resname` → `52-ARG`
- `chainid:resid` → `A:52`

---

## 4. PTA Universal time-series — `[PTA-UNIFIED]`

**Used for:** RMSD, angles, distances. Anything stored under
`/<group>/frame_<N>` with a scalar value per frame.

**Aliases:** `PTA-UNIFIED`, `PTA_UNIFIED`

**Example file:** `pta_unified.ini`

### Specific variables
| Key              | Default       | Notes                                                 |
|------------------|---------------|-------------------------------------------------------|
| `x_axis`         | `time_ns`     | `time_ns`, `time_ps`, or `frame`                      |
| `line_width`     | 1.5           |                                                       |
| `line_alpha`     | 1.0           |                                                       |
| `line_style`     | `solid`       | solid / dashed / dashdot / dotted                     |
| `line_colors`    | tab10 palette | comma-separated list; cycled per series               |
| `cycle_colors`   | `True`        | False = use only the first colour                     |
| `show_std_band`  | `True`        | mean ± std band; only renders on **merged** PTA files |
| `std_band_alpha` | 0.25          |                                                       |
| `plot_every_n`   | 1             | downsample frames                                     |
| `plot_multiple`  | `False`       | `True` = one figure per series, `False` = overlay     |

### Best practices

> [!TIP]
> **Overlay first, split if needed.** `plot_multiple=False` (default)
> gives you one plot per analysis (e.g. one RMSD figure with all
> selections overlaid). Switch to `True` only when curves visually
> collide.
>
> **Pick `time_ns` for ns-scale runs, `time_ps` for ps-scale.** Wrong
> unit on the x-axis is the most common readability bug.
>
> **`show_std_band=True` is harmless on non-merged files** — there's no
> std to plot, so the band silently doesn't render. Leave it on.

---

## 5. RMSF profile — `[PTA-RMSF]`

**Used for:** per-atom RMSF data produced by `pharmacon trajectory
rmsf`. One curve per selection. Different axis logic than the universal
time-series plotter because RMSF has no frame axis.

**Aliases:** `PTA-RMSF`, `PTA_RMSF`

**Example files:** `rmsf.ini` (general), `rmsf_ligand.ini` (per-atom
ligand-style labels)

### Specific variables
| Key                                | Default        | Notes                                                                                              |
|------------------------------------|----------------|----------------------------------------------------------------------------------------------------|
| `x_axis`                           | `resid`        | `resid` / `atom_index` / `position` / `atom_name` (see below)                                      |
| `xtick_format`                     | `""`           | Python format string — fields: `{atom_index}`, `{resid}`, `{resname}`, `{atom_name}`, `{position}` |
| `xtick_rotation`                   | `auto`         | `auto` (90° if >20 labels, else 0°), or any angle                                                  |
| `xtick_max_labels`                 | 200            | thin xticks via stride when above this                                                             |
| `x_min`, `x_max`, `y_min`, `y_max` | `auto`         | `auto` or numeric — clamp axis limits                                                              |
| `line_colors`, `cycle_colors`      | tab10 / `True` | per-series palette                                                                                 |
| `colors_by_label`                  | `""`           | `label1:#hex, label2:#hex` — pin a colour to a named selection                                     |
| `show_std_band`                    | `True`         | renders only on merged RMSF (per-atom `mean`/`std`)                                                |
| `std_band_alpha`                   | 0.25           |                                                                                                    |
| `shading`                          | `""`           | semicolon-separated regions `start,end,color[,alpha[,label]]`                                      |
| `shading_alpha`                    | 0.25           | default alpha when a region omits one                                                              |
| `shading_show_legend`              | `False`        | include labeled shading regions in the legend                                                      |
| `plot_multiple`                    | `False`        | `False` = overlay all selections; `True` = one figure per selection                                |

### `x_axis` modes

| Mode         | Placement               | Default tick labels                | Use when                                                  |
|--------------|-------------------------|------------------------------------|-----------------------------------------------------------|
| `resid`      | residue number          | matplotlib auto                    | One-atom-per-residue selections (`name CA`)               |
| `atom_index` | topology atom index     | matplotlib auto                    | Multi-atom-per-residue selections (`backbone`); monotonic |
| `position`   | 0..N-1 dataset position | matplotlib auto                    | "Just plot them in order"                                 |
| `atom_name`  | 0..N-1 dataset position | atom names (e.g. `C1`, `C2`, `N4`) | Ligands or single-residue zooms                           |

### `xtick_format` examples

| Format                         | Produces                 | Use case                                                   |
|--------------------------------|--------------------------|------------------------------------------------------------|
| `{atom_name}`                  | `C1`, `C2`, `N4`         | Ligand atom-by-atom view (default when `x_axis=atom_name`) |
| `{atom_index} {atom_name}`     | `4523 C1`, `4524 C2`     | Ligand with topology indices for cross-referencing         |
| `{resname}{resid}-{atom_name}` | `LIG901-C1`, `LIG901-C2` | Fully-qualified atom labels                                |
| `{resid}`                      | `1`, `2`, `3`            | Force resid labels even on per-position placement          |

### `shading` syntax

```ini
shading = 0,10,#888888; 11,100,#ffa500
shading = 50,75,#d62728,0.4,active site
```

Each region: `start, end, color[, alpha[, label]]`. The `start`/`end`
values are in the **same units as `x_axis`** — residues for `resid`,
atoms for `atom_index`, positions for `position`/`atom_name`.

### Best practices

> [!TIP]
> **Default to one figure per analysis (`plot_multiple=False`)**, even
> with several selections. Comparison across selections is the main
> reason you run RMSF on more than one.
>
> **For ligand RMSF, copy `rmsf_ligand.ini` as a starting point.**
> It pre-configures `x_axis = atom_name` and `xtick_format =
> "{atom_index} {atom_name}"` which is the convention biologists expect.
>
> **Pin specific colours to specific selections** with
> `colors_by_label = calpha:#d62728, backbone:#1f77b4` so figures stay
> consistent across runs/papers even if selection order changes.
>
> **Use `shading` for domain annotation.** Pair it with
> `shading_show_legend = True` and label each region so reviewers can
> read off boundaries directly from the plot.
>
> **Set `y_min = 0`** if you're comparing RMSF magnitude across runs —
> matplotlib's auto-y-min varies by data and makes comparison harder.

---

## 6. PCA — five plot types

PCA results are stored at `/pca/...`. **PCA is not mergeable**: the
`merge results` subcommand rejects PCA files, and `plot pta` skips PCA
on merged files with a console message.

### 6.1 PC time series — `[PCA-TIMESERIES]`

**Aliases:** `PCA-TIMESERIES`, `PCA_TIMESERIES`
**Example:** `pca_timeseries.ini`

| Key                        | Default   | Notes                                  |
|----------------------------|-----------|----------------------------------------|
| `pcs`                      | `1, 2, 3` | List of components to plot (1-indexed) |
| `line_width`, `line_alpha` | 1.5 / 0.9 |                                        |

### 6.2 PC scatter — `[PCA-SCATTER]`

**Aliases:** `PCA-SCATTER`, `PCA_SCATTER`
**Example:** `pca_scatter.ini`

| Key                | Default   | Notes                                |
|--------------------|-----------|--------------------------------------|
| `pc_x`, `pc_y`     | 1 / 2     | Component to put on each axis        |
| `scatter_size`     | 20        | matplotlib `s=`                      |
| `scatter_alpha`    | 0.8       |                                      |
| `cmap`             | `viridis` | Colour by frame; any matplotlib cmap |
| `disable_colorbar` | `False`   |                                      |

### 6.3 Explained-variance / scree — `[PCA-VARIANCE-RATIO]`

**Aliases:** `PCA-VARIANCE-RATIO`, `PCA_VARIANCE_RATIO`, `PCA-SCREE`, `PCA_SCREE`
**Example:** `pca_variance_ratio.ini`

| Key                     | Default | Notes                            |
|-------------------------|---------|----------------------------------|
| `bar_alpha`             | 0.8     |                                  |
| `cumulative_line`       | `True`  | overlay cumulative variance line |
| `cumulative_line_width` | 2.0     |                                  |
| `cumulative_alpha`      | 0.9     |                                  |

### 6.4 Free-energy surface heatmap — `[PCA-FES-HEATMAP]`

Computes a 2D free-energy surface from a histogram of PC projections.

**Aliases:** `PCA-FES-HEATMAP`, `PCA_FES_HEATMAP`, `PCA FES HEATMAP`, `PCA FREE ENERGY HEATMAP`
**Example:** `pca_fes_heatmap.ini`

| Key                            | Default         | Notes                                    |
|--------------------------------|-----------------|------------------------------------------|
| `components`                   | `(1, 2)`        | Tuple or list of pairs                   |
| `plot_multiple`, `plot_in_one` | `True` / `True` | Layout of multiple component pairs       |
| `bins`                         | 60              | Histogram resolution (10–500)            |
| `smooth_sigma`                 | 1.0             | Gaussian smoothing on the density (0–10) |
| `temperature`                  | 310.0           | Kelvin — used in −kT·ln P                |
| `cmap`                         | `plasma`        | sequential cmap                          |
| `n_levels`                     | 12              | Contour levels                           |

### 6.5 Probability heatmap — `[PCA-PROBABILITY-HEATMAP]`

Same shape as the FES plot but plots P directly (no temperature scaling).

**Aliases:** `PCA-PROBABILITY-HEATMAP`, `PCA_PROBABILITY_HEATMAP`, `PCA PROBABILITY HEATMAP`
**Example:** `pca_probability_heatmap.ini`

Same variables as FES minus `temperature`; default `cmap = viridis`.

### PCA best practices

> [!TIP]
> **Look at the variance ratio first.** If PC1+PC2 explain <50% of the
> motion, the 2D scatter / FES will be misleading.
>
> **`smooth_sigma = 1.0` is a sensible starting point.** Larger values
> smear basins together; smaller values leave granular spikes. Tune on
> the heatmap, not the scatter.
>
> **PCA scatter is the right plot for showing trajectories**; FES /
> probability heatmaps are for showing *populations*. Pick the right one
> for the question you're answering.
>
> **Use `cmap = plasma` or `viridis`** — both perceptually-uniform.
> Avoid `jet` and `rainbow` for publication.

---

## 7. Protein–Ligand Interactions (PLI) — six plot types

All PLI plots share these per-interaction-type colours (override any
to match your house style):

```
color_hydrophobic      = #b39ddb
color_hydrogen_bonds   = #f1c40f
color_pi_cation        = #e67e22
color_pi_stacking      = #2ecc71
color_water_bridge_1   = #3498db
color_ionic            = #ff3333
color_halogen          = #f49ac2
color_metal_contact    = #95a5a6
```

### 7.1 Stacked column type 1 — `[PLI-STACKED-COLUMN-1]`

Per-residue × interaction-type stacked bars.

**Aliases:** `PLI-STACKED-COLUMN-1`, `PLI_STACKED_COLUMN_1`,
`PL-INTERACTIONS-STACKED-COLUMN-1`, `PL-INTERACTIONS_STACKED_COLUMN_1`,
plus `MERGED-…` variants of all four.
**Example:** `pli_stacked_column_1.ini`

| Key                                               | Default                | Notes                                    |
|---------------------------------------------------|------------------------|------------------------------------------|
| `threshold`                                       | 0.02                   | min frequency to include a residue (0–1) |
| `bar_width`                                       | 0.8                    | bar width                                |
| `bar_alpha`                                       | 0.8                    |                                          |
| `bar_edge_color`, `bar_edge_width`                | black / 0.5            |                                          |
| `y_limit_min`, `y_limit_max`                      | auto                   |                                          |
| `error_bars`                                      | `True`                 | only renders on merged files             |
| `error_bars_capsize` / `_alpha` / `_color` / etc. | sensible defaults      |                                          |
| `x_axis_representation`                           | `chainid:resnameresid` | residue label format                     |

### 7.2 Stacked column type 2 — `[PLI-STACKED-COLUMN-2]`

Per-residue backbone vs side-chain stacked bars.

**Aliases:** `PLI-STACKED-COLUMN-2`, `PLI_STACKED_COLUMN_2`,
`PL-INTERACTIONS-STACKED-COLUMN-2`, `PL-INTERACTIONS_STACKED_COLUMN_2`
**Example:** `pli_stacked_column_2.ini`

| Key                | Default   | Notes                       |
|--------------------|-----------|-----------------------------|
| `color_backbone`   | `#b28dff` |                             |
| `color_side_chain` | `#aff8db` |                             |
| `aa3_to_aa1`       | `True`    | one-letter codes by default |
| `threshold`        | 0.02      |                             |

### 7.3 Heatmap type 1 — `[PLI-HEATMAP-1]`

Residue × frame frequency heatmap. Y = residue, X = frame.

**Aliases:** `PLI-HEATMAP-1`, `PLI_HEATMAP_1`, `PL-INTERACTIONS-HEATMAP-1`, `PL_INTERACTIONS_HEATMAP_1`
**Example:** `pli_heatmap_1.ini`

| Key                       | Default                | Notes                              |
|---------------------------|------------------------|------------------------------------|
| `threshold`               | 0.02                   | drop residues below this frequency |
| `cmap`                    | `viridis`              | any matplotlib cmap                |
| `vmin`, `vmax`            | 0.0 / auto             | clamp the colourbar                |
| `interpolation`           | `nearest`              | matplotlib imshow interp           |
| `cbar_orientation`        | `vertical`             | `vertical` or `horizontal`         |
| `cbar_shrink`, `cbar_pad` | 1.0 / 0.04             | colourbar geometry                 |
| `y_axis_representation`   | `chainid:resnameresid` |                                    |

### 7.4 Heatmap type 2 — `[PLI-HEATMAP-2]`

Interaction-type × frame auto-sized heatmap. Width auto-scales with
frame count so the bars stay legible.

**Aliases:** `PLI-HEATMAP-2`, `PLI_HEATMAP_2`, `PL-INTERACTIONS-HEATMAP-2`, `PL_INTERACTIONS_HEATMAP_2`
**Example:** `pli_heatmap_2.ini`

| Key                                                  | Default         | Notes                                        |
|------------------------------------------------------|-----------------|----------------------------------------------|
| `threshold`                                          | 0.02            |                                              |
| `normalize`                                          | `none`          | `none`, `by_frame`, `max1`                   |
| `xtick_max`                                          | 50              | maximum number of x-ticks to show            |
| `drop_empty_rows`                                    | `True`          | omit rows below threshold entirely           |
| `per_frame_in`, `per_inter_in`                       | 0.10 / 0.85     | inches per frame / interaction (auto-sizing) |
| `min_width`, `max_width`, `min_height`, `max_height` | sensible bounds | auto-size envelope                           |

### 7.5 Pie charts — `[PLI-PIE-CHARTS-1]`

One pie chart per residue plus an optional collage figure.

**Aliases:** `PLI-PIE-CHARTS-1`, `PLI_PIE_CHARTS_1`, `PL-INTERACTIONS-PIE-CHARTS-1`, `PL_INTERACTIONS_PIE_CHARTS_1`
**Example:** `pli_pie_charts_1.ini`

| Key                                | Default | Notes                             |
|------------------------------------|---------|-----------------------------------|
| `top_n`                            | 20      | number of top residues to render  |
| `collage`                          | `True`  | also produce one combined figure  |
| `collage_cols`                     | 5       | columns in the collage            |
| `collage_pad`                      | 0.5     | inter-pie spacing                 |
| `make_pdf`                         | `False` | also dump a multi-page PDF        |
| `make_overall`                     | `True`  | also write an overall summary pie |
| `font_size_title`, `font_size_pct` | 28 / 16 |                                   |

### 7.6 Ligand atom monitor — `[PLI-LIGAND-MONITOR]`

Heatmap of residue (rows) × ligand atom (columns) contact frequency.

**Aliases:** `PLI-LIGAND-MONITOR`, `PLI_LIGAND_MONITOR`, `PL-INTERACTIONS-LIGAND-MONITOR`, `PL_INTERACTIONS_LIGAND_MONITOR`
**Example:** `pli_ligand_monitor.ini`

| Key                                  | Default         | Notes |
|--------------------------------------|-----------------|-------|
| `threshold`                          | 0.02            |       |
| `cmap`                               | `viridis`       |       |
| `vmin`, `vmax`                       | 0.0 / 1.0       |       |
| `drop_empty_rows`, `drop_empty_cols` | `True` / `True` |       |

### PLI best practices

> [!TIP]
> **Set `threshold` based on simulation length.** 0.02 (2% of frames)
> is a reasonable default for >500-ns runs; raise to 0.05–0.10 if you
> have very short trajectories and want to filter noise.
>
> **Use `aa3_to_aa1 = True` for stacked columns**, off for heatmaps.
> One-letter codes are compact for narrow bars; three-letter codes are
> clearer when there's plenty of vertical space.
>
> **Keep the per-interaction colour palette consistent across your
> paper.** It's worth setting `color_*` overrides in `all_plots.ini`
> once and referencing the same file for every figure.
>
> **For ligand-atom monitors, set both `vmin = 0.0` and `vmax = 1.0`**
> when you want to compare runs — the auto-vmax otherwise rescales per
> figure.

---

## 8. Protein–Protein Interactions (PPI) — three plot types

### 8.1 Timeline pairs — `[PPI-TIMELINE-PAIRS]`

Top-N residue-pair × frame timeline heatmap.

**Aliases:** `PPI-TIMELINE-PAIRS`, `PPI_TIMELINE_PAIRS`, `PP-TIMELINE`, `PP_TIMELINE`
**Example:** `ppi_timeline_pairs.ini`

| Key               | Default   | Notes                                     |
|-------------------|-----------|-------------------------------------------|
| `top_n`           | 50        | top N pairs by frequency                  |
| `threshold`       | 0.02      | minimum total frequency to include a pair |
| `drop_empty_rows` | `True`    |                                           |
| `cmap`            | `gnuplot` |                                           |
| `vmin`, `vmax`    | 0.0 / 1.0 |                                           |
| `interpolation`   | `nearest` |                                           |

### 8.2 Residue × residue heatmap — `[PPI-HEATMAP]`

Symmetric residue × residue contact-frequency map.

**Aliases:** `PPI-HEATMAP`, `PPI_HEATMAP`, `PP-HEATMAP`, `PP_HEATMAP`
**Example:** `ppi_heatmap.ini`

| Key            | Default    | Notes                                                     |
|----------------|------------|-----------------------------------------------------------|
| `threshold`    | 0.02       |                                                           |
| `min_total`    | 0.0        | drop rows/columns whose row-sum is below this             |
| `top_n`        | 0          | 0 = no limit; >0 keeps only the N most-contacted residues |
| `symmetric`    | `True`     | mirror upper/lower triangles                              |
| `cmap`         | `viridis`  |                                                           |
| `vmin`, `vmax` | 0.0 / auto |                                                           |

### 8.3 Stacked column — `[PPI-STACKED-COLUMN]`

Per residue-pair stacked column (one bar per pair, stacked by
interaction type).

**Aliases:** `PPI-STACKED-COLUMN`, `PPI_STACKED_COLUMN`, `PP-STACKED`, `PP_STACKED`
**Example:** `ppi_stacked_column.ini`

| Key                      | Default             | Notes                          |
|--------------------------|---------------------|--------------------------------|
| `threshold`              | 0.05                |                                |
| `bar_width`, `bar_alpha` | 0.9 / 1.0           |                                |
| `error_bars`             | `True`              | only renders on merged files   |
| `legend_n_col`           | 5                   | columns in the off-axis legend |
| All `color_*` overrides  | same palette as PLI | hydrophobic, hydrogen-bonds, … |

### PPI best practices

> [!TIP]
> **Default to `[PPI-HEATMAP]` when you don't know which protein–protein
> contacts matter.** It surfaces hotspots quickly. Then drill into
> specific pairs with `[PPI-TIMELINE-PAIRS]`.
>
> **`top_n` is a soft filter, `threshold` is a hard one.** Set
> `threshold = 0.02` to drop noise, then `top_n = 50` (timeline) or
> `top_n = 30` (heatmap) to crop the leading set.
>
> **For asymmetric systems** (e.g. complex monomer A on top, monomer B
> on bottom), set `symmetric = False` on the heatmap.

---

## 9. Common gotchas & troubleshooting

### "I changed a key but nothing happened"

- **Typo in the key name** — unknown keys are silently ignored. Compare
  against the per-plot tables above.
- **Wrong section alias** — verify you used a valid alias for the
  intended plot (case-insensitive, but the alias has to match a
  registered one). A typo here is *not* silent: you'll see an
  "unknown section" warning at the end of the run.
- **Wrong plot subcommand** — `plot pta` only reads sections for
  PTA/PCA/RMSF/PLI/PPI; for example `plot pca`-only sections won't render
  through `plot pta`. (In practice `plot pta` handles every section, so
  this is rare.)

### Coercion warnings

Every settings class falls back to defaults when a value is invalid
(bad colour, malformed enum, out-of-range number). Warnings go to
stderr — read them, they're how you spot bad keys without crashing the
whole run. Use `--maxwarnings N` to skip plots whose warnings exceed a
threshold.

### Editable text in Illustrator / Inkscape

Always-on for SVG and PDF — no INI switch. Implementation: module-level
`plt.rcParams["svg.fonttype"]="none"` and `plt.rcParams["pdf.fonttype"]=42`
at the top of `pharmacon/plotter/{universal,interactions}.py`. The
renderer must have the font installed; the default `dejavu sans` ships
with matplotlib and is broadly available.

### Merged-file behaviour

| Plot                      | On merged files                                   |
|---------------------------|---------------------------------------------------|
| PTA-UNIFIED               | Draws mean line + std band (`show_std_band=True`) |
| PTA-RMSF                  | Draws mean line + std band per atom               |
| PCA-*                     | **Skipped** (PCA is not mergeable)                |
| PLI / PPI stacked columns | Draw error bars (`error_bars=True`)               |
| PLI / PPI heatmaps        | No special change                                 |

### Hex colours — quote them

`#` is the INI comment character, so an unquoted hex code like
`#1f77b4` becomes an empty value (everything after `#` is stripped).
The bundled example files quote all hex values; do the same in your own
configs:

```ini
color_hydrophobic   = "#b39ddb"           # single value — quote
line_colors         = "#1f77b4", "#ff7f0e", "#2ca02c"  # list — quote each
colors_by_label     = calpha:"#d62728", backbone:"#1f77b4"
shading             = 0,10,"#888888"; 11,100,"#ffa500"  # quote the colour part
```

Symptom of forgetting: a warning like `line_colors empty, restoring
defaults`, and a plot that uses only the first default colour.

### Booleans vs integers

The INI parser tries numeric coercion **before** boolean literal
matching, so values like `y_min = 0` and `bins = 1` stay as integers
instead of becoming `False` / `True`. Word literals (`true`, `false`,
`yes`, `no`, `on`, `off`) still parse as booleans. Bool fields also
accept native int `0` / `1`, so `flag = 0` works either way.

### Where settings classes live

If a key isn't in the tables above and you need to dig deeper:

| Family                      | Source                                       |
|-----------------------------|----------------------------------------------|
| Universal time-series       | `src/pharmacon/constants/plots/universal.py` |
| RMSF                        | `src/pharmacon/constants/plots/rmsf.py`      |
| PCA (5)                     | `src/pharmacon/constants/plots/pca.py`       |
| PLI (6)                     | `src/pharmacon/constants/plots/pli.py`       |
| PPI (3)                     | `src/pharmacon/constants/plots/ppi.py`       |
| Validation helpers, aliases | `src/pharmacon/constants/plots/_base.py`     |

---

## 10. Editing workflow tips

1. **Start from an example.** Don't write an INI from scratch — copy
   the closest example and edit. Defaults are sensible.
2. **Use the master file (`all_plots.ini`) for batch runs.** One file,
   one source of truth, works with every `plot` subcommand.
3. **Render to SVG, finish in Illustrator/Inkscape.** Editable text and
   vector lines mean you can polish typography and add overlays without
   re-rendering.
4. **Pair `fig_format = svg` with `fig_dpi = 300`** — DPI affects
   raster fallbacks (when SVG embeds raster elements like heatmap
   images) but vector content is resolution-independent.
5. **Override colours in one place.** Set `color_*` or `colors_by_label`
   in `all_plots.ini` so every figure picks up the same palette.
6. **Re-render after a `merge`.** The `is_merged` flag in the file
   metadata changes plot behaviour (std bands, error bars), so the same
   INI produces a more informative figure on a merged input.

---

## Usage recap

```bash
# Single-plot config
pharmacon plot pta -i run.pta -c examples/plot_ini/rmsf.ini

# Master config (one file, every plot family)
pharmacon plot pta -i run.pta -c examples/plot_ini/all_plots.ini
pharmacon plot pta -i merged.pta -c examples/plot_ini/all_plots.ini -o ./plots/

# Render only files you care about by deleting/commenting out sections
# from your INI; unmatched sections produce a single end-of-run warning.
```
