<div align="center">

# Pharmacon — Plot INI Reference

</div>

This directory holds ready-to-edit plot configuration files for the
`pharmacon plot ...` subcommands. Every plot Pharmacon ships is controlled
by an INI file — copy one of the examples below, tweak, and pass it via `-c`:

```bash
pharmacon plot pta -i run.pta    -o ./plots -c examples/plot_ini/rmsf.ini
pharmacon plot pta -i merged.pta -o ./plots -c examples/plot_ini/all_plots.ini
```

> [!TIP]
> If you don't want to manage many small files, **`all_plots.ini`** is a
> single master file with one section per plot type — pass it to any `plot`
> subcommand and only the relevant sections are read; the rest are silently
> ignored.

> [!IMPORTANT]
> **The single most common mistake is quoting.** The INI parser splits any
> value containing a comma into a *list*, and treats an unquoted `#` as the
> start of a comment. So hex colours, and any multi-value field, must be
> quoted. See **[§2.4 Quoting rules](#24-quoting-rules--read-this)** — it will
> save you a confusing afternoon.

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
| `pl_interactions_all.ini`     | PLI (master)    | (every PLI alias)         | All six PLI plot types in one minimalistic file                  |
| `ppi_timeline_pairs.ini`      | PPI             | `PPI-TIMELINE-PAIRS`      | Top residue-pair timeline heatmap                                |
| `ppi_heatmap.ini`             | PPI             | `PPI-HEATMAP`             | Residue × residue contact frequency heatmap                      |
| `ppi_stacked_column.ini`      | PPI             | `PPI-STACKED-COLUMN`      | Per-pair stacked column                                          |
| `hbonds_heatmap.ini`          | H-bonds         | `HBONDS-HEATMAP`          | Residue × residue H-bond frequency heatmap                       |
| `hbonds_timeline_pairs.ini`   | H-bonds         | `HBONDS-TIMELINE-PAIRS`   | Per-pair H-bond occupancy timeline                               |
| `hbonds_count_per_frame.ini`  | H-bonds         | `HBONDS-COUNT-PER-FRAME`  | Number of H-bonds vs frame/time                                  |
| `hbonds_occupancy.ini`        | H-bonds         | `HBONDS-OCCUPANCY`        | Ranked per-pair occupancy bar                                    |
| `hbonds_network.ini`          | H-bonds         | `HBONDS-NETWORK`          | Residue H-bond network graph                                     |

> [!NOTE]
> Every plot is produced by `pharmacon plot pta`. It auto-discovers the
> analysis groups in your `.pta` file and renders whatever it finds — there is
> no separate `plot rmsf` / `plot pca` subcommand. The INI just styles what
> gets drawn.

---

## 2. How INI files work

### 2.1 Sections and aliases

Each plot type has one or more **aliases**. A section header uses any of them —
case-insensitive, hyphens or underscores both accepted:

```ini
[PTA-UNIFIED]    ; same plot
[PTA_UNIFIED]    ; same plot
[pta-unified]    ; same plot
```

Every plot family's aliases are listed in the per-plot reference below.

### 2.2 Unknown keys are silently ignored

The parser only consumes keys it recognises for a given plot class. This is
*deliberate* — it lets one master file (`all_plots.ini`) hold sections for
every plot type without errors. The downside: **a typo in a key name is
silent.** If a plot looks stubbornly default-styled, re-check your key spelling
against the per-plot tables.

### 2.3 Unknown sections produce a warning

If a `[SECTION]` doesn't match any registered alias, you'll see an
"unknown section" warning at the end of the run listing them. (Unlike key
typos, section typos are *not* silent.)

### 2.4 Quoting rules — read this

> [!WARNING]
> The INI loader (`configobj`) does two things that trip people up:
> 1. **A comma splits a value into a list.** `a, b, c` becomes `["a","b","c"]`.
> 2. **An unquoted `#` starts a comment.** Everything after it is discarded.
>
> Hex colours contain `#`; several fields take comma-separated payloads.
> Get the quoting wrong and the value is silently mangled.

The rules, by field shape:

| Field shape | Example fields | How to quote |
|-------------|----------------|--------------|
| **Single colour** | `bg_color`, `color_hydrophobic`, `node_color`, `bar_color` | Quote the value: `bg_color = "#ffffff"` |
| **List of colours** | `line_colors` | Quote **each** item: `line_colors = "#1f77b4", "#ff7f0e"` |
| **Multi-value string** (commas/`#` *inside* one value) | `shading`, `colors_by_label`, `alter_chains_str`, `alter_segments_str` | Quote the **whole** value: `shading = "0,10,#888888; 11,100,#ffa500"` |

```ini
# Single colours — quote
color_hydrophobic = "#b39ddb"
bg_color          = "#ffffff"

# A list — quote EACH element (each is its own list item)
line_colors       = "#1f77b4", "#ff7f0e", "#2ca02c"

# Multi-value strings — quote the WHOLE thing (commas/# live inside one value)
colors_by_label   = "calpha:#d62728, backbone:#1f77b4"
shading           = "0,10,#888888; 11,100,#ffa500"
alter_chains_str  = "{'A': (1, 300)}"
```

> [!CAUTION]
> Do **not** rely on per-item quoting inside a multi-value string
> (e.g. `shading = 0,10,"#888888"`). The comma still splits it into a list
> *before* the quotes are interpreted, and the value is mangled. Quote the
> **entire** value instead. Named colours (`red`, `dodgerblue`) avoid the `#`
> problem but still get comma-split, so multi-value strings always need the
> whole-value quote.

### 2.5 Values: scalars, lists, and "auto"

- Numeric, string, and boolean values parse as expected (`True`/`False`,
  `true`/`false`, `1`/`0`, `yes`/`no`, `on`/`off`).
- **Empty value = "use the default", silently.** `vmax =` (blank) lets the
  plot auto-pick, with no warning. This is the documented way to say "auto"
  for any numeric/limit/size field.
- **`auto`** is also accepted explicitly on the RMSF axis-limit fields
  (`x_min`, `x_max`, `y_min`, `y_max`).

### 2.6 Validation & coercion

Every settings field is validated. An **invalid** value (out-of-range number,
mistyped enum, malformed colour, unknown colormap) **does not crash** — it
emits a warning and falls back to the documented default, so one bad value in a
long file doesn't abort the whole render.

> [!IMPORTANT]
> The distinction matters for `--maxwarnings`:
> - **Empty / unset** field → silent default (does **not** count as a warning).
> - **Invalid** value (e.g. `cmap = nope`, `fig_format = xyz`) → one coercion
>   warning + default.
>
> Each subcommand has `--maxwarnings N` (default `0`). A plot whose settings
> produce **more than N** coercion warnings is *skipped*. With the default `0`,
> a single invalid value will skip that plot — read the warnings, they tell you
> exactly which key to fix. Empty "auto" fields are safe under `--maxwarnings 0`.

### 2.7 Editable text in SVG / PDF

The plotter modules set matplotlib's `svg.fonttype="none"` and
`pdf.fonttype=42`, so **every SVG/PDF Pharmacon writes keeps real `<text>` /
TrueType text** — open in Illustrator or Inkscape and retype/recolour labels
without redrawing. No INI switch needed (the renderer just needs the font
installed; the default `dejavu sans` ships with matplotlib).

### 2.8 Universal output settings

| Key                                 | Accepted values                                 | Notes                                                            |
|-------------------------------------|-------------------------------------------------|------------------------------------------------------------------|
| `fig_format`                        | `png`, `jpg`/`jpeg`, `svg`, `pdf`, `tif`/`tiff` | Use `svg`/`pdf` for publication; invalid → coerced to `png`      |
| `fig_dpi`                           | 50 – 2000                                       | 300 for documents, 600 for posters, 150 for slides              |
| `fig_size_width`, `fig_size_height` | inches                                          | matplotlib's natural unit                                        |
| `fig_transparent`                   | bool                                            | `True` strips the background. **JPG has no alpha**, so it's ignored there |
| `tight_layout`                      | bool                                            | Auto-fit ticks/labels; off by default for some hand-laid-out plots |
| `bg_color`                          | colour                                          | Default `white`                                                  |

---

## 3. Variables shared by *most* plots

These appear across nearly every settings class. Per-plot tables below list
only the *additional* / *plot-specific* variables — assume this shared block is
also available.

### 3.1 Figure
| Key                                  | Default | Notes                                            |
|--------------------------------------|---------|--------------------------------------------------|
| `fig_size_width` / `fig_size_height` | varies  | inches                                           |
| `fig_dpi`                            | 300–800 | clamped 50–2000                                  |
| `fig_basename`                       | varies  | output stem (extension from `fig_format`)        |
| `fig_format`                         | `png`   | png / jpg / svg / pdf / tif                      |
| `fig_transparent`                    | `False` | strip background                                 |
| `tight_layout`                       | varies  | auto-fit margins                                 |
| `bg_color`                           | `white` | any colour string                                |

### 3.2 Title and axis labels
| Key                                  | Default | Notes                              |
|--------------------------------------|---------|------------------------------------|
| `fig_title`                          | varies  | per-plot sensible default          |
| `x_label`, `y_label`                 | varies  | per-plot sensible default          |
| `disable_title`                      | `False` | omit the title                     |
| `disable_x_label`, `disable_y_label` | `False` | omit axis labels                   |
| `disable_x_axis`, `disable_y_axis`   | `False` | omit the whole axis spine + ticks  |
| `disable_ticks`                      | `False` | omit ticks but keep the axis       |

### 3.3 Fonts
| Key                                             | Default       | Notes                                                                      |
|-------------------------------------------------|---------------|----------------------------------------------------------------------------|
| `font_family`                                   | `dejavu sans` | Any installed font; unknown → falls back to `dejavu sans` (with a warning) |
| `font_size_title`                               | 12–18         |                                                                            |
| `font_size_label`                               | 10–16         | axis labels                                                                |
| `font_size_ticks`                               | 8–10          | tick numbers                                                               |
| `font_size_legend`                              | 8–10          |                                                                            |
| `font_size_x`, `font_size_y`                    | 8–10          | per-axis tick fonts (heatmaps)                                             |
| `font_size_cbar`                                | 8–12          | colourbar tick font                                                        |
| `font_weight_title`                             | `bold`        | `normal`, `bold`, `light`, `medium`, `heavy`, `semibold`, `book`, …        |
| `font_weight_label`, `font_weight_legend`, etc. | `normal`      |                                                                            |

### 3.4 Grid
| Key              | Default             | Notes                                                          |
|------------------|---------------------|----------------------------------------------------------------|
| `enable_grid`    | varies              |                                                                |
| `grid_style`     | `dashed`            | `solid`, `dashed`, `dashdot`, `dotted`, `-`, `--`, `-.`, `:`   |
| `grid_alpha`     | 0.3–0.5             | 0=invisible, 1=opaque                                          |
| `grid_color`     | `lightgray`/`black` | heatmaps only                                                  |
| `grid_linewidth` | 0.2                 | heatmaps only                                                  |

### 3.5 Legend
| Key                    | Default | Notes                                                  |
|------------------------|---------|--------------------------------------------------------|
| `disable_legend`       | `False` |                                                        |
| `legend_loc`           | `best`  | matplotlib loc strings: `upper right`, `lower left`, … |
| `legend_frame`         | `True`  | frame around the legend                                |
| `legend_alpha`         | 1.0     | frame transparency                                     |
| `legend_n_col`         | varies  | columns; useful for stacked-column legends             |
| `legend_bbox_y`        | varies  | vertical offset for off-axis legends                   |
| `legend_margin_bottom` | varies  | bottom margin so an off-axis legend isn't clipped      |

### 3.6 Residue-label rewriting (PLI / PPI / H-bonds)

Several plots show residue identifiers. Rewrite them via:

| Key                                                                     | Default                | Notes                                                                          |
|-------------------------------------------------------------------------|------------------------|--------------------------------------------------------------------------------|
| `representation` *or* `x_axis_representation` / `y_axis_representation`  | `chainid:resnameresid` | Combine `resname`, `resid`, `chainid`, `segid` with `:`, `-`, `_`, or nothing  |
| `aa3_to_aa1`                                                            | varies                 | `ALA` → `A`, `LYS` → `K`, …                                                     |
| `renumber`                                                              | `False`                | Reassign residue numbers sequentially from `renumber_int` (first-seen order)   |
| `renumber_int`                                                          | `0`                    | start integer for renumbering                                                  |
| `alter_chains` / `alter_segments`                                       | `False`                | rename chains/segments by resid range (see below)                              |
| `alter_chains_str` / `alter_segments_str`                               | `""`                   | **dict** of `new_id → (lo, hi)` resid ranges — quote the whole value           |

Representation examples:
`chainid:resnameresid` → `A:ARG52`; `resnameresid` → `ARG52`;
`resid-resname` → `52-ARG`; `chainid:resid` → `A:52`.

> [!CAUTION]
> **`alter_chains_str` is a dict of resid ranges, not a rename map**, and it
> must be quoted as a whole (it contains commas and `{}`):
> ```ini
> alter_chains     = True
> alter_chains_str = "{'A': (1, 300), 'B': (301, 600)}"   # resid 1–300 → chain A, 301–600 → chain B
> ```
> Writing `A,B,C` or `A>X` will not work. Same shape for `alter_segments_str`.

---

## 4. PTA Universal time-series — `[PTA-UNIFIED]`

**Used for:** RMSD, angles, distances — any per-frame scalar stored under
`/<group>/frame_<N>`.
**Aliases:** `PTA-UNIFIED`, `PTA_UNIFIED`
**Example:** `pta_unified.ini`

### Specific variables
| Key              | Default       | Notes                                                                          |
|------------------|---------------|--------------------------------------------------------------------------------|
| `x_axis`         | `time_ns`     | `time_ns`, `time_ps`, `time_us`, or `frame_index` (`frame` accepted as alias)  |
| `line_width`     | 1.5           |                                                                                |
| `line_alpha`     | 1.0           |                                                                                |
| `line_style`     | `solid`       | solid / dashed / dashdot / dotted                                              |
| `line_colors`    | tab10 palette | comma-separated list, cycled per series — **quote each hex**                   |
| `cycle_colors`   | `True`        | False = use only the first colour                                              |
| `show_std_band`  | `True`        | mean ± std band; renders only on **merged** PTA files                          |
| `std_band_alpha` | 0.25          |                                                                                |
| `plot_every_n`   | 1             | downsample: keep every Nth frame                                               |
| `plot_multiple`  | `False`       | `True` = one figure per series, `False` = overlay                              |

> [!NOTE]
> Y-axis units are auto-appended when the data carries a `units` attribute and
> you keep the default `y_label` (e.g. angles → `Value (degrees)`). Set
> `y_label` explicitly to override.

> [!TIP]
> - **Overlay first, split if needed.** Switch `plot_multiple=True` only when
>   curves visually collide.
> - **Pick the right time unit.** Wrong unit on the x-axis is the most common
>   readability bug; use `time_ns` for ns-scale runs.
> - **`show_std_band=True` is harmless on non-merged files** — there's no std
>   to draw, so it silently does nothing. Leave it on.

---

## 5. RMSF profile — `[PTA-RMSF]`

**Used for:** per-atom RMSF from `pharmacon trajectory rmsf`. One curve per
selection; no frame axis.
**Aliases:** `PTA-RMSF`, `PTA_RMSF`
**Examples:** `rmsf.ini` (general), `rmsf_ligand.ini` (ligand atom-by-atom)

### Specific variables
| Key                                | Default        | Notes                                                                                              |
|------------------------------------|----------------|----------------------------------------------------------------------------------------------------|
| `x_axis`                           | `resid`        | `resid` / `atom_index` / `position` / `atom_name`                                                  |
| `xtick_format`                     | `""`           | Python format string — fields: `{atom_index}`, `{resid}`, `{resname}`, `{atom_name}`, `{position}` |
| `xtick_rotation`                   | `auto`         | `auto` (90° if many labels) or any angle                                                           |
| `xtick_max_labels`                 | 200            | thin xticks via stride above this                                                                  |
| `x_min`, `x_max`, `y_min`, `y_max` | `auto`         | `auto` or numeric                                                                                  |
| `line_colors`, `cycle_colors`      | tab10 / `True` | per-series palette — **quote each hex**                                                            |
| `colors_by_label`                  | `""`           | pin colours to selections — **quote the whole value** (see below)                                  |
| `show_std_band`, `std_band_alpha`  | `True` / 0.25  | renders only on merged RMSF                                                                        |
| `shading`                          | `""`           | shaded regions — **quote the whole value** (see below)                                             |
| `shading_alpha`                    | 0.25           | default alpha when a region omits one                                                              |
| `shading_show_legend`              | `False`        | include labelled regions in the legend                                                            |
| `plot_multiple`                    | `False`        | `False` = overlay selections; `True` = one figure each                                            |

### `x_axis` modes
| Mode         | Placement               | Default tick labels                | Use when                                       |
|--------------|-------------------------|------------------------------------|------------------------------------------------|
| `resid`      | residue number          | matplotlib auto                    | One-atom-per-residue selections (`name CA`)    |
| `atom_index` | topology atom index     | matplotlib auto                    | Multi-atom-per-residue selections (`backbone`) |
| `position`   | 0..N-1 dataset position | matplotlib auto                    | "Just plot them in order"                      |
| `atom_name`  | 0..N-1 dataset position | atom names (`C1`, `C2`, `N4`)      | Ligands or single-residue zooms                |

### `shading` & `colors_by_label` — quote the whole value
```ini
# region: start,end,color[,alpha[,label]] ; semicolon-separated rows
shading         = "0,10,#888888; 50,75,#d62728,0.4,active site"
colors_by_label = "calpha:#d62728, backbone:#1f77b4"
```
`start`/`end` are in the **same units as `x_axis`** (residues for `resid`,
atoms for `atom_index`, positions otherwise).

> [!TIP]
> - **Default to overlay (`plot_multiple=False`)** — comparing selections is
>   usually the point.
> - **For ligands, start from `rmsf_ligand.ini`** (`x_axis=atom_name`,
>   atom-index+name tick labels).
> - **Pin colours with `colors_by_label`** so figures stay consistent across
>   runs even if selection order changes.
> - **Set `y_min = 0`** when comparing RMSF magnitude across runs.

---

## 6. PCA — five plot types

PCA results are stored one group per selection: `/pca_<selection>/…`
(e.g. `pca_protein`). `plot pta` renders every `pca_*` group it finds.

> [!WARNING]
> **PCA is not mergeable.** `merge results` rejects PCA files, and `plot pta`
> skips PCA on merged files with a console message.

### 6.1 PC time series — `[PCA-TIMESERIES]`
**Aliases:** `PCA-TIMESERIES`, `PCA_TIMESERIES` · **Example:** `pca_timeseries.ini`

| Key                        | Default   | Notes                                  |
|----------------------------|-----------|----------------------------------------|
| `pcs`                      | `1, 2, 3` | components to plot (1-indexed)         |
| `line_width`, `line_alpha` | 1.5 / 0.9 |                                        |

### 6.2 PC scatter — `[PCA-SCATTER]`
**Aliases:** `PCA-SCATTER`, `PCA_SCATTER` · **Example:** `pca_scatter.ini`

| Key                | Default   | Notes                                |
|--------------------|-----------|--------------------------------------|
| `pc_x`, `pc_y`     | 1 / 2     | component per axis                   |
| `scatter_size`     | 20        | matplotlib `s=`                      |
| `scatter_alpha`    | 0.8       |                                      |
| `cmap`             | `viridis` | colour by frame; invalid → coerced   |
| `disable_colorbar` | `False`   |                                      |

### 6.3 Explained-variance / scree — `[PCA-VARIANCE-RATIO]`
**Aliases:** `PCA-VARIANCE-RATIO`, `PCA_VARIANCE_RATIO`, `PCA-SCREE`, `PCA_SCREE` · **Example:** `pca_variance_ratio.ini`

| Key                     | Default | Notes                            |
|-------------------------|---------|----------------------------------|
| `bar_alpha`             | 0.8     |                                  |
| `cumulative_line`       | `True`  | overlay cumulative variance line |
| `cumulative_line_width` | 2.0     |                                  |
| `cumulative_alpha`      | 0.9     |                                  |

### 6.4 Free-energy surface heatmap — `[PCA-FES-HEATMAP]`
**Aliases:** `PCA-FES-HEATMAP`, `PCA_FES_HEATMAP` · **Example:** `pca_fes_heatmap.ini`

| Key                            | Default         | Notes                                    |
|--------------------------------|-----------------|------------------------------------------|
| `components`                   | `(1, 2)`        | component pair                           |
| `plot_multiple`, `plot_in_one` | `True` / `True` | layout of multiple pairs                 |
| `bins`                         | 60              | histogram resolution (10–500)            |
| `smooth_sigma`                 | 1.0             | Gaussian smoothing (0–10)                |
| `temperature`                  | 310.0           | Kelvin — used in −kT·ln P                |
| `cmap`                         | `plasma`        | sequential cmap                          |
| `n_levels`                     | 12              | contour levels                           |

### 6.5 Probability heatmap — `[PCA-PROBABILITY-HEATMAP]`
**Aliases:** `PCA-PROBABILITY-HEATMAP`, `PCA_PROBABILITY_HEATMAP` · **Example:** `pca_probability_heatmap.ini`
Same as FES minus `temperature`; default `cmap = viridis`.

> [!TIP]
> - **Read the variance ratio first.** If PC1+PC2 explain <50%, the 2D
>   scatter/FES is misleading.
> - **Scatter = trajectory, FES/probability = populations.** Pick the plot for
>   the question.
> - **Use `viridis`/`plasma`** (perceptually uniform); avoid `jet`/`rainbow`.

---

## 7. Protein–Ligand Interactions (PLI) — six plot types

All PLI plots share the per-interaction-type colour palette (override to match
your house style — **quote each hex**):

```ini
color_hydrophobic    = "#b39ddb"
color_hydrogen_bonds = "#f1c40f"
color_pi_cation      = "#e67e22"
color_pi_stacking    = "#2ecc71"
color_water_bridge_1 = "#3498db"
color_ionic          = "#ff3333"
color_halogen        = "#f49ac2"
color_metal_contact  = "#95a5a6"
```

> [!TIP]
> Prefer one file for the whole family? **`pl_interactions_all.ini`** carries a
> minimalistic section for all six PLI plots — copy it and tweak.

> [!NOTE]
> PLI plots are built from **interaction modes** (`mode1`/`mode2`/`mode3`):
> mode1 counts every occurrence, mode2 de-duplicates per frame (occupancy),
> mode3 is a hybrid. Stacked-column-1 renders once **per mode** (merged and
> non-merged); the other PLI plots render once and **only on non-merged files**
> — a merged `.pta` produces just stacked-column-1.

### 7.1 Stacked column type 1 — `[PLI-STACKED-COLUMN-1]`
Per-residue × interaction-type stacked bars.
**Aliases:** `PLI-STACKED-COLUMN-1`, `PLI_STACKED_COLUMN_1`, `PL-INTERACTIONS-STACKED-COLUMN-1`, … plus `MERGED-…` variants.
**Example:** `pli_stacked_column_1.ini`

| Key                                | Default                | Notes                                    |
|------------------------------------|------------------------|------------------------------------------|
| `threshold`                        | 0.02                   | min frequency to include a residue (0–1) |
| `bar_width`, `bar_alpha`           | 0.8 / 0.8              |                                          |
| `bar_edge_color`, `bar_edge_width` | black / 0.5            |                                          |
| `y_limit_min`, `y_limit_max`       | auto                   |                                          |
| `error_bars` (+ `_capsize`, …)     | `True`                 | only renders on merged files             |
| `x_axis_representation`            | `chainid:resnameresid` | residue label format                     |

### 7.2 Stacked column type 2 — `[PLI-STACKED-COLUMN-2]`
Per-residue backbone vs side-chain stacked bars.
**Aliases:** `PLI-STACKED-COLUMN-2`, `PLI_STACKED_COLUMN_2`, … · **Example:** `pli_stacked_column_2.ini`

| Key                | Default   | Notes                       |
|--------------------|-----------|-----------------------------|
| `color_backbone`   | `#b28dff` | quote it                    |
| `color_side_chain` | `#aff8db` | quote it                    |
| `aa3_to_aa1`       | `True`    | one-letter codes by default |
| `threshold`        | 0.02      |                             |

### 7.3 Heatmap type 1 — `[PLI-HEATMAP-1]`
Residue × frame frequency heatmap.
**Aliases:** `PLI-HEATMAP-1`, `PLI_HEATMAP_1`, … · **Example:** `pli_heatmap_1.ini`

| Key                       | Default                | Notes                              |
|---------------------------|------------------------|------------------------------------|
| `threshold`               | 0.02                   | drop residues below this frequency |
| `cmap`                    | `viridis`              | invalid → coerced                  |
| `vmin`, `vmax`            | 0.0 / auto             | clamp the colourbar                |
| `interpolation`           | `nearest`              | imshow interpolation; invalid → coerced |
| `cbar_orientation`        | `vertical`             | `vertical` or `horizontal`         |
| `cbar_shrink`, `cbar_pad` | 1.0 / 0.04             |                                    |
| `y_axis_representation`   | `chainid:resnameresid` |                                    |

### 7.4 Heatmap type 2 — `[PLI-HEATMAP-2]`
Interaction-type × frame auto-sized heatmap.
**Aliases:** `PLI-HEATMAP-2`, `PLI_HEATMAP_2`, … · **Example:** `pli_heatmap_2.ini`

| Key                                                  | Default         | Notes                                        |
|------------------------------------------------------|-----------------|----------------------------------------------|
| `threshold`                                          | 0.02            |                                              |
| `normalize`                                          | `none`          | `none`, `by_frame`, `max1` (invalid → `none` + warning) |
| `xtick_max`                                          | 50              | max x-ticks to show                          |
| `drop_empty_rows`                                    | `True`          | omit empty rows                              |
| `per_frame_in`, `per_inter_in`                       | 0.10 / 0.85     | inches per frame / interaction (auto-sizing) |
| `min_width`, `max_width`, `min_height`, `max_height` | bounds          | auto-size envelope                           |

### 7.5 Pie charts — `[PLI-PIE-CHARTS-1]`
One pie per residue + optional collage.
**Aliases:** `PLI-PIE-CHARTS-1`, `PLI_PIE_CHARTS_1`, … · **Example:** `pli_pie_charts_1.ini`

| Key                                | Default | Notes                                                       |
|------------------------------------|---------|-------------------------------------------------------------|
| `fig_size`                         | 5.0     | **single square size** (pie charts use this, not `fig_size_width/height`) |
| `top_n`                            | 20      | top residues to render (0 = all)                            |
| `collage`, `collage_cols`, `collage_pad` | `True` / 5 / 0.5 | combined figure layout                       |
| `make_pdf`                         | `False` | also dump a multi-page PDF                                  |
| `make_overall`                     | `True`  | also write an overall summary pie                          |
| `disable_labels`, `disable_autopct`| `False` | hide slice labels / percentages                            |
| `font_size_title`, `font_size_pct` | 28 / 16 |                                                            |

> [!CAUTION]
> Pie charts need backbone/side-chain (`SC`/`BB`) data, which requires the
> **ligand atoms to have names** in the topology. Ligands parametrised from
> SMILES/Maestro often arrive nameless — if so the pies are skipped
> ("No SC/BB data available"). Re-prepare the ligand with atom names.

### 7.6 Ligand atom monitor — `[PLI-LIGAND-MONITOR]`
Residue (rows) × ligand atom (cols) contact-frequency heatmap.
**Aliases:** `PLI-LIGAND-MONITOR`, `PLI_LIGAND_MONITOR`, … · **Example:** `pli_ligand_monitor.ini`

| Key                                  | Default         | Notes |
|--------------------------------------|-----------------|-------|
| `threshold`                          | 0.02            |       |
| `cmap`                               | `viridis`       |       |
| `vmin`, `vmax`                       | 0.0 / 1.0       | set both to compare across runs |
| `drop_empty_rows`, `drop_empty_cols` | `True` / `True` |       |

> [!TIP]
> - **Scale `threshold` with run length** — 0.02 suits long runs; raise to
>   0.05–0.10 to filter noise on short trajectories.
> - **`aa3_to_aa1=True` for stacked columns, `False` for heatmaps.**
> - **Set `vmin=0, vmax=1`** on the ligand monitor to compare runs (auto-vmax
>   rescales per figure otherwise).

---

## 8. Protein–Protein Interactions (PPI) — three plot types

### 8.1 Timeline pairs — `[PPI-TIMELINE-PAIRS]`
Top-N residue-pair × frame timeline heatmap.
**Aliases:** `PPI-TIMELINE-PAIRS`, `PPI_TIMELINE_PAIRS`, `PP-TIMELINE`, `PP_TIMELINE` · **Example:** `ppi_timeline_pairs.ini`

| Key               | Default   | Notes                                     |
|-------------------|-----------|-------------------------------------------|
| `top_n`           | 50        | top N pairs by frequency                  |
| `threshold`       | 0.02      | min total frequency to include a pair     |
| `drop_empty_rows` | `True`    |                                           |
| `xtick_max`       | 50        | **thin frame ticks to ~this many** (avoids overlapping x-labels on long runs) |
| `cmap`            | `gnuplot` |                                           |
| `vmin`, `vmax`    | 0.0 / 1.0 |                                           |

### 8.2 Residue × residue heatmap — `[PPI-HEATMAP]`
**Aliases:** `PPI-HEATMAP`, `PPI_HEATMAP`, `PP-HEATMAP`, `PP_HEATMAP` · **Example:** `ppi_heatmap.ini`

| Key            | Default    | Notes                                                     |
|----------------|------------|-----------------------------------------------------------|
| `threshold`    | 0.02       |                                                           |
| `min_total`    | 0.0        | drop residues whose row-sum is below this                 |
| `top_n`        | 0          | 0 = no limit; >0 keeps the N most-contacted residues      |
| `symmetric`    | `True`     | mirror upper/lower triangles                              |
| `cmap`         | `viridis`  |                                                           |
| `vmin`, `vmax` | 0.0 / auto |                                                           |

### 8.3 Stacked column — `[PPI-STACKED-COLUMN]`
Per residue-pair stacked column.
**Aliases:** `PPI-STACKED-COLUMN`, `PPI_STACKED_COLUMN`, `PP-STACKED`, `PP_STACKED` · **Example:** `ppi_stacked_column.ini`

| Key                      | Default             | Notes                          |
|--------------------------|---------------------|--------------------------------|
| `threshold`              | 0.05                |                                |
| `bar_width`, `bar_alpha` | 0.9 / 1.0           |                                |
| `error_bars`             | `True`              | only renders on merged files   |
| `legend_n_col`           | 5                   | columns in the off-axis legend |
| All `color_*` overrides  | same palette as PLI |                                |

> [!TIP]
> - **Start with `[PPI-HEATMAP]`** to find hotspots, then drill into pairs with
>   `[PPI-TIMELINE-PAIRS]`.
> - **`threshold` is the hard filter, `top_n` the soft crop.**
> - **`symmetric = False`** for genuinely directional/asymmetric systems.

---

## 9. Hydrogen Bonds (H-bonds) — five plot types

Produced by `pharmacon trajectory h-bonds`. H-bonds are residue↔residue
contacts of a single interaction type, so the residue-pair plots reuse the PPI
machinery with H-bond titles. The analysis writes only **mode1** (count-all)
and **mode2** (once-per-frame ⇒ true occupancy) — there is no mode3 for
H-bonds. Per-frame datasets also store `time_ps`, enabling time-based x-axes.

### 9.1 Residue × residue heatmap — `[HBONDS-HEATMAP]`
**Aliases:** `HBONDS-HEATMAP`, `HBONDS_HEATMAP`, `HB-HEATMAP`, `HB_HEATMAP` · **Example:** `hbonds_heatmap.ini`
Same knobs as `[PPI-HEATMAP]` (`threshold`, `min_total`, `top_n`, `symmetric`,
`cmap`, `vmin`/`vmax`, representation). Default title "Hydrogen Bond Contact
Frequency".

### 9.2 Pair occupancy timeline — `[HBONDS-TIMELINE-PAIRS]`
The classic H-bond persistence plot: pairs (rows) × frame (cols).
**Aliases:** `HBONDS-TIMELINE-PAIRS`, `HB-TIMELINE`, … · **Example:** `hbonds_timeline_pairs.ini`
Same knobs as `[PPI-TIMELINE-PAIRS]`, including **`xtick_max`** (thin frame
ticks). Non-merged only.

### 9.3 H-bonds per frame — `[HBONDS-COUNT-PER-FRAME]`
Number of H-bonds vs frame/time (line).
**Aliases:** `HBONDS-COUNT-PER-FRAME`, `HB-COUNT`, … · **Example:** `hbonds_count_per_frame.ini`

| Key            | Default       | Notes                                                                   |
|----------------|---------------|-------------------------------------------------------------------------|
| `x_axis`       | `frame_index` | `frame_index`, `time_ns`, `time_ps`, `time_us` — auto-labels the x-axis  |
| `line_width`, `line_style`, `line_colors` | universal | quote each hex                                       |
| `enable_grid`, `grid_style`, `grid_alpha` | grid styling |                                                       |

### 9.4 Occupancy ranking — `[HBONDS-OCCUPANCY]`
Ranked horizontal bar of per-pair occupancy (built from mode2).
**Aliases:** `HBONDS-OCCUPANCY`, `HB-OCCUPANCY`, … · **Example:** `hbonds_occupancy.ini`

| Key                                          | Default   | Notes                            |
|----------------------------------------------|-----------|----------------------------------|
| `top_n`                                      | 25        | most-occupied pairs (0 = all)    |
| `threshold`                                  | 0.0       | drop pairs below this occupancy  |
| `bar_color`, `bar_alpha`, `bar_edge_color`   | styling   | quote `bar_color` hex            |
| `representation`, `aa3_to_aa1`               | labels    |                                  |

### 9.5 Network graph — `[HBONDS-NETWORK]`
Residue H-bond network: nodes = residues, edges = H-bonds weighted by
occupancy. **Built to stay readable on large proteins.**
**Aliases:** `HBONDS-NETWORK`, `HB-NETWORK`, … · **Example:** `hbonds_network.ini`

| Key                              | Default   | Notes                                                                   |
|----------------------------------|-----------|-------------------------------------------------------------------------|
| `mode`                           | `mode2`   | source mode (mode2 ⇒ occupancy in 0–1)                                  |
| `threshold`                      | 0.3       | keep only bonds at/above this occupancy — **the main de-clutter lever** |
| `top_n`                          | 200       | keep the N strongest edges (0 = all)                                    |
| `min_seq_sep`                    | 0         | drop \|i−j\| ≤ this within a chain (e.g. 2 hides sequential backbone)   |
| `largest_component`              | `False`   | keep only the largest connected network                                 |
| `label_top_n`                    | 15        | label only the N highest-degree hubs (0 = none, −1 = all)               |
| `layout`                         | `spring`  | `spring`, `kamada_kawai`, `circular`                                    |
| `seed`                           | 42        | reproducible spring layout                                              |
| `node_size_min`, `node_size_max` | 50 / 600 | node size scales with degree                                            |
| `node_color`, `node_edge_color`  | colours  | quote hex                                                               |
| `edge_cmap`                      | `viridis` | edges coloured by occupancy                                             |
| `edge_width_min`, `edge_width_max`, `edge_alpha` | styling | width scales with occupancy                             |

> [!WARNING]
> On a large protein the **raw** network is a hairball. Tame it with
> `threshold` (try 0.3–0.5 for "stable" bonds), `top_n`, `largest_component=True`,
> and `min_seq_sep=2` to hide expected α-helix backbone bonds. `label_top_n`
> keeps labels legible by tagging only the hubs.

> [!TIP]
> The **heatmap** is the matrix view of the same data (most scalable, no
> hairball); the **network** adds topology (hubs/clusters). The **occupancy**
> bar and **count-per-frame** line are the two summaries reviewers expect.

---

## 10. Common gotchas & troubleshooting

### "I changed a key but nothing happened"
- **Typo in the key name** — unknown keys are silently ignored. Check the
  per-plot table.
- **Wrong section alias** — must match a registered alias (an "unknown
  section" warning prints at the end if not).
- **You quoted wrong** — see §2.4. A mangled `line_colors`/`shading`/
  `colors_by_label` is the usual culprit.

### "My plot got skipped"
- **`--maxwarnings` (default 0).** One *invalid* value skips the plot. Read the
  warnings — they name the offending key. Empty "auto" fields don't count.
- **Empty data after filtering.** A high `threshold`/`top_n` can filter
  everything out; the plot is skipped with a clear message.

### Hex colours & multi-value strings — quote them
`#` is the INI comment character and commas split into lists (§2.4). Symptom of
forgetting: a `line_colors empty, restoring defaults` warning, default colours,
or mangled shading/labels.

### Booleans vs integers
Numeric coercion runs **before** boolean matching, so `y_min = 0` and `bins = 1`
stay integers (not `False`/`True`). Word literals (`true`, `false`, `yes`,
`no`, `on`, `off`) still parse as booleans; bool fields also accept `0`/`1`.

### Merged-file behaviour
| Plot                        | On merged files                                   |
|-----------------------------|---------------------------------------------------|
| PTA-UNIFIED                 | Mean line + std band (`show_std_band=True`)       |
| PTA-RMSF                    | Mean line + std band per atom                     |
| PCA-*                       | **Skipped** (PCA is not mergeable)                |
| PLI stacked-column-1        | Renders per mode, with error bars                 |
| PLI heatmaps / pie / monitor / stacked-column-2 | **Skipped** — non-merged only |
| PPI stacked column          | Draws error bars (`error_bars=True`)              |
| PPI / H-bond heatmaps       | Render (no special change)                        |
| H-bond occupancy / network  | Render                                            |
| H-bond timeline / count     | Per-frame plots — non-merged only                 |

### Where settings classes live
| Family                      | Source                                       |
|-----------------------------|----------------------------------------------|
| Universal time-series       | `src/pharmacon/constants/plots/universal.py` |
| RMSF                        | `src/pharmacon/constants/plots/rmsf.py`      |
| PCA (5)                     | `src/pharmacon/constants/plots/pca.py`       |
| PLI (6)                     | `src/pharmacon/constants/plots/pli.py`       |
| PPI (3)                     | `src/pharmacon/constants/plots/ppi.py`       |
| H-bonds (5)                 | `src/pharmacon/constants/plots/hbonds.py`    |
| Validation helpers, aliases | `src/pharmacon/constants/plots/_base.py`     |

---

## 11. Editing workflow tips

1. **Start from an example** — copy the closest file and edit; defaults are
   sensible.
2. **Use `all_plots.ini` for batch runs** — one source of truth for every plot
   family.
3. **Render to SVG, finish in Illustrator/Inkscape** — text stays editable.
4. **Pair `fig_format = svg` with `fig_dpi = 300`** — DPI only affects raster
   fallbacks (e.g. embedded heatmap images); vector content is
   resolution-independent.
5. **Override colours in one place** — set `color_*` / `colors_by_label` in
   `all_plots.ini` so every figure shares the palette.
6. **Re-render after a `merge`** — the `is_merged` metadata flips std-bands /
   error-bars on, so the same INI yields a richer figure.

---

## Usage recap

```bash
# Single-plot config
pharmacon plot pta -i run.pta -o ./plots -c examples/plot_ini/hbonds_network.ini

# Master config (one file, every plot family)
pharmacon plot pta -i run.pta    -o ./plots -c examples/plot_ini/all_plots.ini
pharmacon plot pta -i merged.pta -o ./plots -c examples/plot_ini/all_plots.ini

# Be strict: skip any plot whose settings emit a coercion warning
pharmacon plot pta -i run.pta -o ./plots -c my.ini --maxwarnings 0
```
