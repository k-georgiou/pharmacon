"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Module :mod:`pharmacon.plotter.universal`.
"""
import csv
import json
import math
import warnings
from itertools import combinations
from pathlib import Path
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, minimum_filter


from pharmacon.logger import get_logger, PharmaconLogger
from pharmacon.constants.plots import (PlotUniversalSettings, RMSFPlotSettings,
                                       PCAPlotTimeSeriesSettings, PCAPlotScatterSettings,
                                       PCAPlotVarianceRatioSettings, PCAPlotProbabilityHeatmapSettings,
                                       PCAPlotFESHeatmapSettings)




__all__ = [
    "logger",
    "plot_pta_timeseries_from_file",
    "plot_pta_rmsf_from_file",
    "plot_pca_timeseries_from_file",
    "plot_pca_scatter_from_file",
    "plot_pca_variance_ratio_from_file",
    "plot_pca_probability_from_file",
    "plot_pca_fes_from_file",
    "export_pca_summary",
]


warnings.filterwarnings("ignore")


# Make text in SVG/PDF outputs remain editable in Illustrator/Inkscape
# (real <text> elements / TrueType) rather than outlined to paths. Affects
# every savefig in this module; raster formats ignore these rcParams.
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42


logger: PharmaconLogger = get_logger(__name__)



def _get_x_value(attrs, frame: int, x_axis: str) -> float:
    """
    Determines the x-axis value based on the provided x_axis type and attributes.

    The function processes and converts the x-axis value depending on whether it is
    a frame index or time-based (picoseconds, nanoseconds, microseconds). It first
    attempts to use the direct attribute from the provided attrs dictionary. If that
    is not available, it falls back to "time_ps" and calculates the required time value
    by converting it to the appropriate unit. If no valid x_axis type is provided, a
    ValueError is raised.

    Supported x_axis:
      - "frame_index"
      - "time_ps"
      - "time_ns"
      - "time_us"

    Logic:
      1) If x_axis == frame_index -> frame
      2) If attrs has the requested unit -> use it
      3) Else fall back to time_ps and convert
      4) Else fall back to frame

    :param attrs: Dictionary of attributes containing possible x-axis values or fallback
                  time_ps attribute
    :type attrs: dict
    :param frame: The frame index to use as the x-axis value if x_axis is "frame_index"
    :type frame: int
    :param x_axis: The type of x-axis data to use; must be one of "time_ps", "time_ns",
                   "time_us", or "frame_index"
    :type x_axis: str
    :return: The computed x-axis value based on the provided settings
    :rtype: float
    :raises ValueError: If the x_axis type is invalid or conversion is not possible
    """
    x_axis = (x_axis or "time_ps").strip().lower()

    if x_axis == "frame_index":
        return float(frame)

    if x_axis in ("time_ps", "time_ns", "time_us"):
        # try direct attribute first (if writer stored it)
        direct = attrs.get(x_axis, None)
        if direct is not None:
            try:
                return float(direct)
            except Exception:
                pass

        # fall back to time_ps then convert
        t_ps = attrs.get("time_ps", None)
        if t_ps is None:
            return float(frame)

        t_ps = float(t_ps)
        if x_axis == "time_ps":
            return t_ps
        if x_axis == "time_ns":
            return t_ps / 1_000.0
        if x_axis == "time_us":
            return t_ps / 1_000_000.0

    raise ValueError("settings.x_axis must be one of: time_ps, time_ns, time_us, frame_index")


def plot_pta_timeseries_from_file(pta_file,
                                  *,
                                  group_name: str,
                                  settings: PlotUniversalSettings,
                                  out_dir: Path,
                                  is_merged: bool = False
                                  ) -> None:
    """
    Generate and save time series plots from PTA dataset files, with customizable
    plotting settings for specific groups, such as "distances", "angles", or
    other data categories available in the PTA file.

    This function selectively processes data based on filtering criteria in
    the `settings` parameter and groups data into series based on specific
    keys ("label", "kind", "method", etc.) depending on the group type. The plots
    are saved as image files into the specified output directory.

    Unified timeseries plotter for rmsd / angles / distances PTA groups.

    STORAGE (non-merged):
      /<group>/frame_<i>/<group>  -> rows:
        rmsd:      {label, value}
        angles:    {label, kind, value}
        distances: {label, method, distance}

    STORAGE (merged):
      /<group>/frame_<i>/<group>  -> rows:
        rmsd:      {label, mean, std, n}
        angles:    {label, kind, mean, std, n}
        distances: {label, method, mean, std, n}

    :param pta_file: PTA file containing the dataset with group and frame hierarchy.
    :param group_name: Specific group name within the PTA file to process, such as "distances", "angles", etc.
    :param settings: An instance of PlotUniversalSettings containing user-defined plotting parameters.
    :param out_dir: Output directory where generated plot files will be saved.
    :param is_merged: Boolean flag indicating whether the data represents merged/averaged statistics.

    :return: None
    """
    logger.debug(
        "plot_pta_timeseries_from_file: group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )
    settings.validate()

    if group_name not in pta_file.file:
        raise ValueError(f"Group '{group_name}' not found in PTA file.")

    # new: choose x axis (default time_ps)
    x_axis = getattr(settings, "x_axis", "time_ps")
    logger.debug("x_axis='%s' plot_every_n=%d plot_multiple=%s",
                 x_axis, settings.plot_every_n, settings.plot_multiple)

    data: dict[str, dict[str, list]] = {}  # label -> {"x":[], "y":[], "std":[] or None}

    frames = list(pta_file._iter_frames(group_name))
    if not frames:
        raise ValueError(f"No frames found under group '{group_name}'")
    logger.debug("Discovered %d frame(s) under '%s'", len(frames), group_name)

    frames_used = 0
    rows_total = 0
    for frame in frames:
        if frame % settings.plot_every_n != 0:
            continue

        dset_path = f"{group_name}/frame_{frame}/{group_name}"
        if dset_path not in pta_file.file:
            logger.trace("Missing dataset '{}' (frame {}), skipping", dset_path, frame)
            continue

        frames_used += 1
        dset = pta_file.file[dset_path]
        attrs = dset.attrs

        # new: compute x based on requested axis
        x_val = _get_x_value(attrs, frame, x_axis)

        for row in dset:
            rec = json.loads(row)
            rows_total += 1

            # Build a stable label key. Merged records store the averaged
            # value under "mean"; non-merged records use the original
            # per-subcommand field names ("value" or "distance").
            if group_name == "angles":
                base = str(rec.get("label", "")).strip()
                kind = str(rec.get("kind", "")).strip()
                series_key = f"{base}:{kind}" if kind else base

                y = float(rec.get("mean", 0.0)) if is_merged else float(rec.get("value", 0.0))
                std = float(rec.get("std", 0.0)) if is_merged else None

            elif group_name == "distances":
                base = str(rec.get("label", "")).strip()
                method = str(rec.get("method", "")).strip()
                series_key = f"{base}:{method}" if method else base

                y = float(rec.get("mean", 0.0)) if is_merged else float(rec.get("distance", 0.0))
                std = float(rec.get("std", 0.0)) if is_merged else None

            else:
                # rmsd (and anything else with {label, value})
                series_key = str(rec.get("label", "")).strip()
                y = float(rec.get("mean", 0.0)) if is_merged else float(rec.get("value", 0.0))
                std = float(rec.get("std", 0.0)) if is_merged else None

            if not series_key:
                continue

            if series_key not in data:
                data[series_key] = {"x": [], "y": [], "std": [] if is_merged else None}

            data[series_key]["x"].append(x_val)
            data[series_key]["y"].append(y)

            if is_merged:
                data[series_key]["std"].append(std)

    if not data:
        raise ValueError("No PTA data collected (empty after filtering).")

    logger.debug(
        "Collected %d row(s) across %d frame(s) → %d series: %s",
        rows_total, frames_used, len(data), sorted(data.keys()),
    )

    def _plot(labels: list[str], suffix: str = "") -> None:
        fig, ax = plt.subplots(
            figsize=(settings.fig_size_width, settings.fig_size_height),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
        )
        plt.rcParams["font.family"] = settings.font_family

        n_colors = max(1, len(settings.line_colors))

        for i, label in enumerate(labels):
            content = data[label]
            x = np.asarray(content["x"], dtype=float)
            y = np.asarray(content["y"], dtype=float)

            order = np.argsort(x)
            x = x[order]
            y = y[order]

            color = settings.line_colors[(i % n_colors) if settings.cycle_colors else 0]

            ax.plot(
                x, y,
                color=color,
                linewidth=settings.line_width,
                alpha=settings.line_alpha,
                linestyle=settings.line_style,
                label=label,
            )

            if is_merged and settings.show_std_band:
                std = np.asarray(content["std"], dtype=float)[order]
                ax.fill_between(
                    x, y - std, y + std,
                    color=color,
                    alpha=settings.std_band_alpha,
                    linewidth=0,
                )

        if not settings.disable_title:
            ax.set_title(
                settings.fig_title,
                fontsize=settings.font_size_title,
                fontweight=settings.font_weight_title,
            )

        ax.set_xlabel(settings.x_label, fontsize=settings.font_size_label, fontweight=settings.font_weight_label)
        ax.set_ylabel(settings.y_label, fontsize=settings.font_size_label, fontweight=settings.font_weight_label)
        ax.tick_params(labelsize=settings.font_size_ticks)

        if settings.enable_grid:
            ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

        if not settings.disable_legend:
            ax.legend(
                fontsize=settings.font_size_legend,
                frameon=settings.legend_frame,
                framealpha=settings.legend_alpha,
            )

        if settings.tight_layout:
            fig.tight_layout()

        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{group_name}_timeseries{suffix}.{settings.fig_format}"
        fig.savefig(out_path, dpi=settings.fig_dpi, transparent=settings.fig_transparent)
        plt.close(fig)

        logger.info(f"Saved PTA plot → {out_path}")

    keys = sorted(data.keys())

    if settings.plot_multiple:
        for k in keys:
            safe = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in k)[:120]
            _plot([k], suffix=f"_{safe}")
    else:
        _plot(keys, suffix="")


def plot_pta_rmsf_from_file(pta_file,
                            *,
                            group_name: str = "rmsf",
                            settings: RMSFPlotSettings,
                            out_dir: Path,
                            is_merged: bool = False
                            ) -> None:
    """
    Render per-atom RMSF profiles, one curve per selection.

    Unlike :func:`plot_pta_timeseries_from_file`, this plotter operates on
    the per-selection layout produced by the ``trajectory rmsf`` subcommand:
    there is no frame axis. Each selection's records are read from
    ``/<group>/selection_<label>/atoms``.

    STORAGE (non-merged):
      /<group>/selection_<label>/atoms -> rows:
        {atom_index, resid, resname, atom_name, rmsf}

    STORAGE (merged):
      /<group>/selection_<label>/atoms -> rows:
        {atom_index, resid, resname, atom_name, mean, std, n}

    :param pta_file: Open PTA file handle.
    :param group_name: HDF5 group containing RMSF data (default ``"rmsf"``).
    :param settings: :class:`RMSFPlotSettings` with figure/font/line options.
    :param out_dir: Directory to write the plot into.
    :param is_merged: ``True`` if the PTA file was produced by
        ``merge results`` (records carry ``mean``/``std`` instead of
        ``rmsf``).
    :return: None
    """
    logger.debug(
        "plot_pta_rmsf_from_file: group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )
    settings.validate()

    if group_name not in pta_file.file:
        raise ValueError(f"Group '{group_name}' not found in PTA file.")

    x_axis = getattr(settings, "x_axis", "resid")
    logger.debug("x_axis='%s' plot_multiple=%s is_merged=%s",
                 x_axis, settings.plot_multiple, is_merged)

    # Per-label container:
    #   "x":         list[float]     – x-axis values
    #   "y":         list[float]     – RMSF values (or merged mean)
    #   "std":       list[float]|None
    #   "atom_info": list[dict]      – per-row {atom_index, resid, resname,
    #                                  atom_name, position} for xtick label formatting
    data: dict[str, dict[str, list]] = {}

    labels = list(pta_file._iter_selections(group_name))
    if not labels:
        raise ValueError(f"No RMSF selections found under group '{group_name}'.")
    logger.debug("Discovered %d selection(s) under '%s'", len(labels), group_name)

    rows_total = 0
    for label in labels:
        dset_path = f"{group_name}/selection_{label}/atoms"
        if dset_path not in pta_file.file:
            logger.trace("Missing dataset '{}' for selection {}",
                         dset_path, label)
            continue
        dset = pta_file.file[dset_path]

        xs: list = []
        ys: list = []
        stds: list = []
        atom_info: list[dict] = []

        for idx, row in enumerate(dset):
            rec = json.loads(row)
            rows_total += 1

            # x value depends on requested axis
            if x_axis == "resid":
                x_val = float(rec.get("resid", 0))
            elif x_axis == "atom_index":
                x_val = float(rec.get("atom_index", 0))
            else:  # "position" or "atom_name" (both place atoms at 0..N-1)
                x_val = float(idx)

            # y value: merged records carry "mean" (+ "std"); non-merged carry "rmsf"
            if is_merged:
                y_val = float(rec.get("mean", 0.0))
                s_val = float(rec.get("std", 0.0))
                stds.append(s_val)
            else:
                y_val = float(rec.get("rmsf", 0.0))

            xs.append(x_val)
            ys.append(y_val)
            atom_info.append({
                "atom_index": rec.get("atom_index", idx),
                "resid":      rec.get("resid", 0),
                "resname":    rec.get("resname", ""),
                "atom_name":  rec.get("atom_name", ""),
                "position":   idx,
            })

        if not xs:
            continue

        data[label] = {
            "x": xs, "y": ys,
            "std": stds if is_merged else None,
            "atom_info": atom_info,
        }

    if not data:
        raise ValueError("No RMSF data collected (empty after filtering).")

    logger.debug("Collected %d row(s) across %d selection(s)",
                 rows_total, len(data))

    # Resolve the effective xtick label template:
    #   - explicit `xtick_format` always wins
    #   - else x_axis=atom_name → "{atom_name}" as a sensible default
    #   - else "" → matplotlib's automatic numeric ticks
    effective_format = settings.xtick_format
    if not effective_format and x_axis == "atom_name":
        effective_format = "{atom_name}"

    def _plot(plot_labels: list[str], suffix: str = "") -> None:
        fig, ax = plt.subplots(
            figsize=(settings.fig_size_width, settings.fig_size_height),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
        )
        plt.rcParams["font.family"] = settings.font_family

        n_colors = max(1, len(settings.line_colors))

        # Shaded background regions. Drawn first so they sit behind everything
        # else; zorder is set explicitly because tight_layout / fill_between
        # can otherwise reshuffle layering.
        shading_regions = getattr(settings, "shading_regions", []) or []
        for region in shading_regions:
            ax.axvspan(
                region["start"], region["end"],
                color=region["color"],
                alpha=region["alpha"],
                label=region["label"] if settings.shading_show_legend else None,
                zorder=0,
                linewidth=0,
            )

        # For custom xtick labels we need the per-atom info, sorted by x;
        # collect it from the first plotted label (or union over labels if
        # overlay mode). In practice atom_name labels are used with a single
        # selection (one ligand, one fragment) — we use the first label as the
        # source-of-truth and assume other overlaid selections share positions.
        tick_positions: list[float] = []
        tick_labels: list[str] = []

        colors_by_label = getattr(settings, "colors_by_label_map", {}) or {}

        for i, label in enumerate(plot_labels):
            content = data[label]
            x = np.asarray(content["x"], dtype=float)
            y = np.asarray(content["y"], dtype=float)
            info = content["atom_info"]

            order = np.argsort(x, kind="stable")
            x_sorted = x[order]
            y_sorted = y[order]

            # Per-label override wins; otherwise fall back to the cycling palette.
            if label in colors_by_label:
                color = colors_by_label[label]
            else:
                color = settings.line_colors[(i % n_colors) if settings.cycle_colors else 0]

            ax.plot(
                x_sorted, y_sorted,
                color=color,
                linewidth=settings.line_width,
                alpha=settings.line_alpha,
                linestyle=settings.line_style,
                label=label,
            )

            if is_merged and settings.show_std_band and content["std"] is not None:
                std = np.asarray(content["std"], dtype=float)[order]
                ax.fill_between(
                    x_sorted, y_sorted - std, y_sorted + std,
                    color=color,
                    alpha=settings.std_band_alpha,
                    linewidth=0,
                )

            # Build tick metadata from the first selection only.
            if i == 0 and effective_format:
                info_sorted = [info[int(j)] for j in order]
                tick_positions = [float(x_sorted[k]) for k in range(len(x_sorted))]
                tick_labels = []
                for rec_info in info_sorted:
                    try:
                        tick_labels.append(effective_format.format(**rec_info))
                    except Exception:
                        # Fallback if format misbehaves on a specific row.
                        tick_labels.append(str(rec_info.get("atom_name", "")))

        if not settings.disable_title:
            ax.set_title(
                settings.fig_title,
                fontsize=settings.font_size_title,
                fontweight=settings.font_weight_title,
            )

        ax.set_xlabel(settings.x_label,
                      fontsize=settings.font_size_label,
                      fontweight=settings.font_weight_label)
        ax.set_ylabel(settings.y_label,
                      fontsize=settings.font_size_label,
                      fontweight=settings.font_weight_label)
        ax.tick_params(labelsize=settings.font_size_ticks)

        # Apply custom xtick labels if a template is in effect.
        if effective_format and tick_positions:
            n_ticks = len(tick_positions)
            max_labels = max(1, int(settings.xtick_max_labels))
            if n_ticks > max_labels:
                stride = max(1, n_ticks // max_labels)
                tick_positions = tick_positions[::stride]
                tick_labels = tick_labels[::stride]
                logger.debug(
                    "Thinning xtick labels: %d → %d (stride=%d)",
                    n_ticks, len(tick_positions), stride,
                )

            if settings.xtick_rotation == "auto":
                rotation = 90.0 if len(tick_labels) > 20 else 0.0
            else:
                rotation = float(settings.xtick_rotation)

            ax.set_xticks(tick_positions)
            ax.set_xticklabels(
                tick_labels,
                rotation=rotation,
                ha=("right" if 0 < rotation < 180 else "center"),
            )

        if settings.enable_grid:
            ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

        # Axis limits ("auto" leaves matplotlib's autoscaling untouched).
        cur_xmin, cur_xmax = ax.get_xlim()
        cur_ymin, cur_ymax = ax.get_ylim()
        new_xmin = settings.x_min if settings.x_min != "auto" else cur_xmin
        new_xmax = settings.x_max if settings.x_max != "auto" else cur_xmax
        new_ymin = settings.y_min if settings.y_min != "auto" else cur_ymin
        new_ymax = settings.y_max if settings.y_max != "auto" else cur_ymax
        ax.set_xlim(new_xmin, new_xmax)
        ax.set_ylim(new_ymin, new_ymax)

        if not settings.disable_legend:
            ax.legend(
                fontsize=settings.font_size_legend,
                frameon=settings.legend_frame,
                framealpha=settings.legend_alpha,
            )

        if settings.tight_layout:
            fig.tight_layout()

        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{group_name}_profile{suffix}.{settings.fig_format}"
        fig.savefig(out_path, dpi=settings.fig_dpi, transparent=settings.fig_transparent)
        plt.close(fig)

        logger.info(f"Saved RMSF plot → {out_path}")

    keys = sorted(data.keys())

    if settings.plot_multiple:
        for k in keys:
            safe = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in k)[:120]
            _plot([k], suffix=f"_{safe}")
    else:
        _plot(keys, suffix="")


def plot_pca_timeseries_from_file(pta_file,
                                  *,
                                  group_name: str = "pca",
                                  settings: PCAPlotTimeSeriesSettings,
                                  out_dir: Path
                                  ) -> None:
    """
    Generate and save a PCA time-series plot from a specified PTA file group.

    This function reads PCA time-series data from a specific group in a PTA file
    and produces a customizable plot. The plot's visual style, dimensions, and
    behavior can be adjusted via the provided settings. The generated plot is
    then saved as an image in the specified output directory.

    :param pta_file: PTA file providing access to the data. The file must support
        operations to check group existence and iterate through frames.
    :type pta_file: PTAFileHandler
    :param group_name: Name of the group in the PTA file containing PCA data. Defaults to "pca".
    :type group_name: str
    :param settings: Configuration object for customizing the appearance and layout
        of the plot, including line styles, dimensions, and labels.
    :type settings: PCAPlotTimeSeriesSettings
    :param out_dir: Output directory where the generated plot image will be saved.
        The directory will be created if it does not already exist.
    :type out_dir: Path
    :return: This function does not return a value. It only saves the plot to a file.
    :rtype: None
    """

    logger.debug(
        "plot_pca_timeseries_from_file: group='%s' pcs=%s out_dir='%s'",
        group_name, settings.pcs, out_dir,
    )
    settings._validate_fields()

    if not pta_file.group_exists(group_name):
        raise ValueError(f"Group '{group_name}' not found in PTA file")

    data = {}

    frame_count = 0
    for frame in pta_file._iter_frames(group_name):
        frame_count += 1

        dset_path = f"{group_name}/frame_{frame}/pca"
        dset = pta_file.file[dset_path]
        time_ps = float(dset.attrs.get("time_ps", frame))

        for row in dset:
            rec = json.loads(row)

            pc = int(rec["pc"])
            if pc not in settings.pcs:
                continue

            if pc not in data:
                data[pc] = {"x": [], "y": []}

            data[pc]["x"].append(time_ps)
            data[pc]["y"].append(float(rec["value"]))

    logger.debug(
        "Collected %d frame(s) → %d PC series %s",
        frame_count, len(data), sorted(data.keys()),
    )

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    for pc, content in sorted(data.items()):
        ax.plot(
            content["x"],
            content["y"],
            linewidth=settings.line_width,
            alpha=settings.line_alpha,
            label=f"PC{pc}",
        )

    if not settings.disable_title:
        ax.set_title(settings.fig_title, fontsize=settings.font_size_title)

    ax.set_xlabel(settings.x_label, fontsize=settings.font_size_label)
    ax.set_ylabel(settings.y_label, fontsize=settings.font_size_label)

    if settings.enable_grid:
        ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

    if not settings.disable_legend:
        ax.legend(fontsize=settings.font_size_legend)

    if settings.tight_layout:
        fig.tight_layout()

    out_path = Path(out_dir) / f"{settings.fig_basename}.{settings.fig_format}"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_path, dpi=settings.fig_dpi)
    plt.close(fig)
    logger.info(f"Saved PCA timeseries → {out_path}")


def plot_pca_scatter_from_file(pta_file,
                               *,
                               group_name: str = "pca",
                               settings: PCAPlotScatterSettings,
                               out_dir: Path
                               ) -> None:
    """
    Generates a PCA scatter plot from the provided file, based on specified settings.
    Reads PCA data associated with particular principal components and plots them with
    visual properties outline by the settings. Saves the resulting visualization to the
    specified output directory.

    :param pta_file: The input file containing PCA data and frame-related information.
    :param group_name: Optional group name in the file to fetch frames and PCA data
                       (default: "pca").
    :param settings: Configuration object containing scatter plot attributes and
                     rendering specifications.
    :param out_dir: Directory path where the generated scatter plot image file will be saved.
    :return: None
    """

    logger.debug(
        "plot_pca_scatter_from_file: group='%s' pc_x=%d pc_y=%d out_dir='%s'",
        group_name, settings.pc_x, settings.pc_y, out_dir,
    )
    settings._validate_fields()

    x_vals = []
    y_vals = []
    times = []

    for frame in pta_file._iter_frames(group_name):

        dset = pta_file.file[f"{group_name}/frame_{frame}/pca"]
        time_ps = float(dset.attrs.get("time_ps", frame))

        pc_values = {}

        for row in dset:
            rec = json.loads(row)
            pc_values[int(rec["pc"])] = float(rec["value"])

        if settings.pc_x in pc_values and settings.pc_y in pc_values:
            x_vals.append(pc_values[settings.pc_x])
            y_vals.append(pc_values[settings.pc_y])
            times.append(time_ps)

    logger.debug("Collected %d point(s) for PC%d vs PC%d",
                 len(x_vals), settings.pc_x, settings.pc_y)

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    sc = ax.scatter(
        x_vals,
        y_vals,
        c=times,
        cmap=settings.cmap,
        s=settings.scatter_size,
        alpha=settings.scatter_alpha,
    )

    if not settings.disable_colorbar:
        plt.colorbar(sc, ax=ax, label="Time (ps)")

    ax.set_xlabel(f"PC{settings.pc_x}")
    ax.set_ylabel(f"PC{settings.pc_y}")

    if not settings.disable_title:
        ax.set_title(settings.fig_title)

    if settings.enable_grid:
        ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

    if settings.tight_layout:
        fig.tight_layout()

    out_path = Path(out_dir) / f"{settings.fig_basename}.{settings.fig_format}"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_path, dpi=settings.fig_dpi)
    plt.close(fig)
    logger.info(f"Saved PCA scatter → {out_path}")


def plot_pca_variance_ratio_from_file(pta_file,
                                      *,
                                      group_name: str = "pca",
                                      settings: PCAPlotVarianceRatioSettings,
                                      out_dir: Path
                                      ) -> None:
    """
    Plots the PCA variance ratio from a provided file and saves the plot to an output
    directory.

    This function reads PCA variance ratios and principal components from a specified
    input file, constructs a bar plot, overlays a cumulative line (if enabled), and
    saves the resulting plot as an image file. The visual appearance of the plot can
    be customized using the provided settings.

    :param pta_file: The input file containing PCA data. Must support an iterator
                     for frames and provide access to datasets by group and frame.
    :type pta_file: Any
    :param group_name: The group name in the input file where PCA data is stored.
                       Defaults to "pca".
    :type group_name: str, optional
    :param settings: The configuration object containing the layout and appearance
                     settings for the PCA variance-ratio plot. Includes attributes
                     like figure size, DPI, axis labels, and grid options.
    :type settings: PCAPlotVarianceRatioSettings
    :param out_dir: The directory where the generated plot will be saved. The method
                    ensures the directory is created if it doesn't exist.
    :type out_dir: Path
    :return: None
    """

    logger.debug(
        "plot_pca_variance_ratio_from_file: group='%s' out_dir='%s'",
        group_name, out_dir,
    )
    settings._validate_fields()

    first_frame = next(pta_file._iter_frames(group_name))
    dset = pta_file.file[f"{group_name}/frame_{first_frame}/pca"]
    logger.debug("Reading variance ratios from first frame dataset '%s'",
                 f"{group_name}/frame_{first_frame}/pca")

    pcs = []
    ratios = []

    for row in dset:
        rec = json.loads(row)
        pcs.append(int(rec["pc"]))
        ratios.append(float(rec["variance_ratio"]))

    logger.debug("Loaded %d component(s) with variance ratios", len(pcs))

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    ax.bar(pcs, ratios, alpha=settings.bar_alpha)

    if settings.cumulative_line:
        cumulative = np.cumsum(ratios)
        ax.plot(
            pcs,
            cumulative,
            linewidth=settings.cumulative_line_width,
            alpha=settings.cumulative_alpha,
        )

    ax.set_xlabel(settings.x_label)
    ax.set_ylabel(settings.y_label)

    if not settings.disable_title:
        ax.set_title(settings.fig_title)

    if settings.enable_grid:
        ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

    if settings.tight_layout:
        fig.tight_layout()

    out_path = Path(out_dir) / f"{settings.fig_basename}.{settings.fig_format}"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_path, dpi=settings.fig_dpi)
    plt.close(fig)
    logger.info(f"Saved PCA variance-ratio plot → {out_path}")


def _get_all_pca_components_from_file(pta_file, *, group_name: str = "pca") -> List[int]:
    """
    Extracts all PCA components from a given file under a specific group.

    This function iterates over frames in the specified group within a file-like
    object, extracts PCA components from a dataset in a JSON format, and ensures
    that at least two unique PCA components are available for further use. If
    no frames, datasets, or insufficient PCA components are found, exceptions
    are raised.

    :param pta_file: The file-like object to extract PCA components from.
    :param group_name: The group under which to search for PCA-related data, default is "pca".
    :type group_name: str
    :return: A sorted list of unique PCA components.
    :rtype: List[int]
    :raises ValueError: If frames cannot be found under the specified group.
    :raises ValueError: If the PCA dataset cannot be located in the file.
    :raises ValueError: If fewer than two PCA components are found.
    """
    frames = list(pta_file._iter_frames(group_name))
    if not frames:
        raise ValueError(f"No frames found under group '{group_name}'")

    first_frame = frames[0]
    dset_path = f"{group_name}/frame_{first_frame}/pca"

    if dset_path not in pta_file.file:
        raise ValueError(f"PCA dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]

    pcs: List[int] = []
    for row in dset:
        rec = json.loads(row)
        pcs.append(int(rec["pc"]))

    pcs = sorted(set(pcs))

    if len(pcs) < 2:
        raise ValueError("At least two PCA components are required for heatmap plotting")

    return pcs


def _normalize_pca_component_pairs(components,
                                   *,
                                   allow_multiple: bool,
                                   all_components: List[int],
                                   ) -> List[Tuple[int, int]]:
    """
    Normalize PCA component pairs.

    This function processes and validates a given set of PCA components to produce
    a list of PCA component pairs. It supports a variety of input formats, including
    a single tuple, a list of tuples, or a list of integers, and ensures the output
    is a consistent list of valid component pairs.

    Normalize PCA component specifications into explicit (pc_x, pc_y) pairs.

    Accepted forms:
        - None                      -> all available component pairs if allow_multiple=True,
                                      otherwise the first two components
        - (1, 2)
        - [1, 2, 3]                 -> pairwise combinations if allow_multiple=True,
                                      otherwise only (1, 2)
        - [(1, 2), (1, 3), (2, 3)]

    :param components: A specification of PCA components to process, which can be
        None, a tuple of two integers, a list of tuples of two integers, or a list
        of integers depending on the desired input.
    :param allow_multiple: A boolean flag indicating whether multiple PCA component
        pairs are allowed in the result.
    :param all_components: A list of all PCA components that can be used to
        generate pairs in the case where `components` is None.
    :return: A list of normalized and validated PCA component pairs as tuples of
        two integers.
    :rtype: List[Tuple[int, int]]
    :raises ValueError: If the input components cannot be parsed or do not meet the
        required criteria.
    """

    if components is None:
        if allow_multiple:
            return list(combinations(all_components, 2))
        return [(all_components[0], all_components[1])]

    if (
        isinstance(components, tuple)
        and len(components) == 2
        and all(isinstance(x, int) for x in components)
    ):
        return [(int(components[0]), int(components[1]))]

    comps_list = list(components)

    if not comps_list:
        raise ValueError("No PCA components specified")

    if all(isinstance(x, tuple) and len(x) == 2 for x in comps_list):
        pairs = [(int(a), int(b)) for a, b in comps_list]
        if not allow_multiple and len(pairs) > 1:
            return [pairs[0]]
        return pairs

    if all(isinstance(x, int) for x in comps_list):
        ints = [int(x) for x in comps_list]
        if len(ints) < 2:
            raise ValueError("At least two PCA components are required")
        if allow_multiple:
            return list(combinations(ints, 2))
        return [(ints[0], ints[1])]

    raise ValueError(
        "Invalid components specification. Use None, (1, 2), [1, 2, 3], or [(1, 2), (1, 3)]."
    )


def _collect_pca_xy_from_file(pta_file,
                              *,
                              group_name: str,
                              pc_x: int,
                              pc_y: int
                              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extracts PCA x and y values for specified principal components (PCs) from a given
    group within a PTA file. The method iterates through all frames in the specified
    group, reading the relevant PCA data for the provided PC indices, and aggregates
    them into separate arrays for PCx and PCy. Any frames that do not contain PCA data
    for the specified PCs are ignored.

    :param pta_file: The PTA file object to extract PCA values from.
    :param group_name: The name of the group within the PTA file where PCA data is located.
    :param pc_x: The index of the principal component for the x-axis values.
    :param pc_y: The index of the principal component for the y-axis values.
    :return: A tuple containing numpy arrays of the x-axis and y-axis PCA values.
    :rtype: Tuple[np.ndarray, np.ndarray]
    :raises ValueError: When no PCA values are found for the specified PCs.
    """
    x_vals: List[float] = []
    y_vals: List[float] = []

    for frame in pta_file._iter_frames(group_name):
        dset_path = f"{group_name}/frame_{frame}/pca"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]

        pc_values: Dict[int, float] = {}
        for row in dset:
            rec = json.loads(row)
            pc_values[int(rec["pc"])] = float(rec["value"])

        if pc_x in pc_values and pc_y in pc_values:
            x_vals.append(pc_values[pc_x])
            y_vals.append(pc_values[pc_y])

    if not x_vals or not y_vals:
        raise ValueError(f"No PCA values found for PC{pc_x} vs PC{pc_y}")

    return np.asarray(x_vals, dtype=float), np.asarray(y_vals, dtype=float)


def _build_pca_probability_surface(x: np.ndarray,
                                   y: np.ndarray,
                                   *,
                                   bins: int,
                                   smooth_sigma: float
                                   ) -> Dict[str, np.ndarray]:
    """
    Build a probability distribution surface using PCA (Principal Component Analysis)
    coordinates. This function calculates a 2D histogram from the input data and
    optionally applies Gaussian smoothing to the histogram to generate a
    density-normalized probability distribution.

    :param x: PCA component values for the first dimension.
    :type x: np.ndarray
    :param y: PCA component values for the second dimension.
    :type y: np.ndarray
    :param bins: Number of bins to use for the 2D histogram.
    :type bins: int
    :param smooth_sigma: Sigma for Gaussian smoothing. If set to 0, no smoothing is
        applied.
    :type smooth_sigma: float
    :return: A dictionary containing the following keys:
        - **H** (*np.ndarray*): The normalized 2D histogram (probability surface).
        - **X** (*np.ndarray*): X-coordinates of the histogram grid.
        - **Y** (*np.ndarray*): Y-coordinates of the histogram grid.
    :rtype: Dict[str, np.ndarray]
    """
    from scipy.ndimage import gaussian_filter

    H, xedges, yedges = np.histogram2d(x, y, bins=bins, density=True)

    if smooth_sigma > 0:
        H = gaussian_filter(H, sigma=smooth_sigma)

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        norm = np.nansum(H) * dx * dy
        if norm > 0:
            H = H / norm

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcenters, ycenters, indexing="ij")

    return {
        "H": H,
        "X": X,
        "Y": Y,
    }


def _build_pca_fes_surface(x: np.ndarray,
                           y: np.ndarray,
                           *,
                           bins: int,
                           smooth_sigma: float,
                           temperature: float
                           ) -> Dict[str, np.ndarray]:
    """
    Constructs a free energy surface (FES) based on principal component analysis (PCA).
    This function utilizes a 2D histogram of input data, applies optional Gaussian smoothing,
    and calculates free energy values from normalized probability densities.

    The input data `x` and `y` are binned according to the specified number of bins.
    If a positive `smooth_sigma` is provided, Gaussian smoothing is applied to the
    histogram. The free energy is then calculated using the Boltzmann relationship and
    the specified temperature in Kelvin.

    F = -RT ln(P), shifted so min(F) = 0

    :param x: 1D array-like, Numerical data for the x-axis in PCA space.
    :param y: 1D array-like, Numerical data for the y-axis in PCA space.
    :param bins: int, Number of bins for the 2D histogram.
    :param smooth_sigma: float, Standard deviation for Gaussian smoothing in histogram
        computations. A value less than or equal to 0 disables smoothing.
    :param temperature: float, Absolute temperature in Kelvin used for free energy
        calculations.
    :return: A dictionary containing the following:
        - "F": 2D numpy array of free energy values [kcal/mol].
        - "X": 2D numpy array of x-axis grid coordinates corresponding to the free
          energy surface.
        - "Y": 2D numpy array of y-axis grid coordinates corresponding to the free
          energy surface.
    """
    from scipy.ndimage import gaussian_filter

    H, xedges, yedges = np.histogram2d(x, y, bins=bins, density=True)

    if smooth_sigma > 0:
        H = gaussian_filter(H, sigma=smooth_sigma)

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        norm = np.nansum(H) * dx * dy
        if norm > 0:
            H = H / norm

    H_safe = H.copy()
    H_safe[H_safe <= 0] = np.nan

    R = 0.0019872041  # kcal/mol/K
    F = -R * temperature * np.log(H_safe)
    F = F - np.nanmin(F)

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcenters, ycenters, indexing="ij")

    return {
        "F": F,
        "X": X,
        "Y": Y,
    }


def plot_pca_probability_from_file(pta_file,
                                   *,
                                   group_name: str = "pca",
                                   settings: PCAPlotProbabilityHeatmapSettings,
                                   out_dir: Path
                                   ) -> None:
    """
    Plots PCA probability heatmaps from a specified file. This function reads
    principal component data from a given PTA file, creates density surfaces
    for specified or all principal component pairs, and generates individual
    or combined heatmap plots based on the specified settings.

    :param pta_file: The file containing Principal Component Analysis (PCA) data.
    :type pta_file: PTAFile
    :param group_name: The group name under which PCA data is structured in the
                       file. Defaults to "pca".
    :type group_name: str
    :param settings: An instance of PCAPlotProbabilityHeatmapSettings providing
                     configuration for the plotting process, such as plot
                     aesthetics and component pairs.
    :type settings: PCAPlotProbabilityHeatmapSettings
    :param out_dir: The output directory to save the generated plots.
                    If the directory does not exist, it will be created.
    :type out_dir: Path
    :return: None
    """
    logger.debug(
        "plot_pca_probability_from_file: group='%s' components=%s "
        "plot_multiple=%s plot_in_one=%s bins=%d smooth_sigma=%s out_dir='%s'",
        group_name, settings.components, settings.plot_multiple,
        settings.plot_in_one, settings.bins, settings.smooth_sigma, out_dir,
    )
    settings._validate_fields()

    if not pta_file.group_exists(group_name):
        raise ValueError(f"Group '{group_name}' not found in PTA file")

    plt.rcParams["font.family"] = settings.font_family

    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_components = _get_all_pca_components_from_file(
        pta_file,
        group_name=group_name,
    )
    logger.debug("Discovered %d PCA component(s): %s",
                 len(all_components), all_components)

    component_pairs = _normalize_pca_component_pairs(
        settings.components,
        allow_multiple=settings.plot_multiple,
        all_components=all_components,
    )
    logger.debug("Resolved %d component pair(s) for density surfaces: %s",
                 len(component_pairs), component_pairs)

    surfaces: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}

    for pc_x, pc_y in component_pairs:
        x, y = _collect_pca_xy_from_file(
            pta_file,
            group_name=group_name,
            pc_x=pc_x,
            pc_y=pc_y,
        )
        logger.trace("Pair PC{} vs PC{}: {} point(s)", pc_x, pc_y, len(x))
        surfaces[(pc_x, pc_y)] = _build_pca_probability_surface(
            x,
            y,
            bins=settings.bins,
            smooth_sigma=settings.smooth_sigma,
        )

    def _levels(H: np.ndarray) -> np.ndarray:
        hmax = np.nanmax(H)
        if not np.isfinite(hmax) or hmax <= 0:
            hmax = 1.0
        return np.linspace(0.0, hmax, settings.n_levels)

    def _draw(ax, surf: Dict[str, np.ndarray], pc_x: int, pc_y: int) -> None:
        levels = _levels(surf["H"])

        cf = ax.contourf(
            surf["X"],
            surf["Y"],
            surf["H"],
            levels=levels,
            cmap=settings.cmap,
        )
        ax.contour(
            surf["X"],
            surf["Y"],
            surf["H"],
            levels=levels,
            linewidths=0.8,
        )

        ax.set_xlabel(settings.x_label or f"Projection on eigenvector {pc_x} (PC{pc_x})")
        ax.set_ylabel(settings.y_label or f"Projection on eigenvector {pc_y} (PC{pc_y})")

        if not settings.disable_title:
            title = settings.fig_title or f"PCA Probability Density: PC{pc_x} vs PC{pc_y}"
            ax.set_title(
                title,
                fontsize=settings.font_size_title,
                fontweight=settings.font_weight_title,
            )

        ax.tick_params(labelsize=settings.font_size_ticks)

        if settings.enable_grid:
            ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

        plt.colorbar(cf, ax=ax, label="Probability density")

    # Always save individual plots for all resolved pairs
    for pc_x, pc_y in component_pairs:
        fig, ax = plt.subplots(
            figsize=(settings.fig_size_width, settings.fig_size_height),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
        )

        _draw(ax, surfaces[(pc_x, pc_y)], pc_x, pc_y)

        if settings.tight_layout:
            fig.tight_layout()

        out_path = out_dir / f"{settings.fig_basename}_pc{pc_x}_vs_pc{pc_y}.{settings.fig_format}"
        fig.savefig(
            out_path,
            dpi=settings.fig_dpi,
            transparent=settings.fig_transparent,
        )
        plt.close(fig)

        logger.info(f"Saved PCA probability heatmap → {out_path}")

    # Optionally also save one combined figure
    if settings.plot_in_one:
        fig, axes = plt.subplots(
            nrows=len(component_pairs),
            ncols=1,
            figsize=(settings.fig_size_width, settings.fig_size_height * len(component_pairs)),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
            squeeze=False,
        )

        for i, (pc_x, pc_y) in enumerate(component_pairs):
            _draw(axes[i][0], surfaces[(pc_x, pc_y)], pc_x, pc_y)

        if settings.tight_layout:
            fig.tight_layout()

        out_path = out_dir / f"{settings.fig_basename}.{settings.fig_format}"
        fig.savefig(
            out_path,
            dpi=settings.fig_dpi,
            transparent=settings.fig_transparent,
        )
        plt.close(fig)

        logger.info(f"Saved combined PCA probability heatmap → {out_path}")


def plot_pca_fes_from_file(pta_file,
                           *,
                           group_name: str = "pca",
                           settings: PCAPlotFESHeatmapSettings,
                           out_dir: Path
                           ) -> None:
    """
    Plots and saves PCA-based Free Energy Surface (FES) heatmaps from data within a specified PTA file.
    The function supports generating individual component pair plots as well as a combined figure
    containing all component pairs resolved. All generated plots are saved to the specified output directory.

    :param pta_file: An object representing the PTA file containing PCA data.
    :param group_name: Name of the group within the PTA file holding PCA-related information. Default is "pca".
    :param settings: A configuration object defining plot settings (e.g., heatmap appearance, labels).
    :param out_dir: Directory path where the resulting plots will be saved.
    :return: None
    """
    logger.debug(
        "plot_pca_fes_from_file: group='%s' components=%s "
        "plot_multiple=%s plot_in_one=%s bins=%d smooth_sigma=%s T=%sK out_dir='%s'",
        group_name, settings.components, settings.plot_multiple,
        settings.plot_in_one, settings.bins, settings.smooth_sigma,
        settings.temperature, out_dir,
    )
    settings._validate_fields()

    if not pta_file.group_exists(group_name):
        raise ValueError(f"Group '{group_name}' not found in PTA file")

    plt.rcParams["font.family"] = settings.font_family

    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_components = _get_all_pca_components_from_file(
        pta_file,
        group_name=group_name,
    )
    logger.debug("Discovered %d PCA component(s): %s",
                 len(all_components), all_components)

    component_pairs = _normalize_pca_component_pairs(
        settings.components,
        allow_multiple=settings.plot_multiple,
        all_components=all_components,
    )
    logger.debug("Resolved %d component pair(s) for FES surfaces: %s",
                 len(component_pairs), component_pairs)

    surfaces: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}

    for pc_x, pc_y in component_pairs:
        x, y = _collect_pca_xy_from_file(
            pta_file,
            group_name=group_name,
            pc_x=pc_x,
            pc_y=pc_y,
        )
        logger.trace("Pair PC{} vs PC{}: {} point(s)", pc_x, pc_y, len(x))
        surfaces[(pc_x, pc_y)] = _build_pca_fes_surface(
            x,
            y,
            bins=settings.bins,
            smooth_sigma=settings.smooth_sigma,
            temperature=settings.temperature,
        )

    def _levels(F: np.ndarray) -> np.ndarray:
        finite = F[np.isfinite(F)]
        if finite.size == 0:
            return np.linspace(0.0, 1.0, settings.n_levels)

        fmax = np.nanpercentile(finite, 95.0)
        if not np.isfinite(fmax) or fmax <= 0:
            fmax = np.nanmax(finite)
        if not np.isfinite(fmax) or fmax <= 0:
            fmax = 1.0

        return np.linspace(0.0, fmax, settings.n_levels)

    def _draw(ax, surf: Dict[str, np.ndarray], pc_x: int, pc_y: int) -> None:
        levels = _levels(surf["F"])

        cf = ax.contourf(
            surf["X"],
            surf["Y"],
            surf["F"],
            levels=levels,
            cmap=settings.cmap,
        )
        ax.contour(
            surf["X"],
            surf["Y"],
            surf["F"],
            levels=levels,
            linewidths=0.8,
        )

        ax.set_xlabel(settings.x_label or f"Projection on eigenvector {pc_x} (PC{pc_x})")
        ax.set_ylabel(settings.y_label or f"Projection on eigenvector {pc_y} (PC{pc_y})")

        if not settings.disable_title:
            title = settings.fig_title or f"PCA Free Energy Surface: PC{pc_x} vs PC{pc_y}"
            ax.set_title(
                title,
                fontsize=settings.font_size_title,
                fontweight=settings.font_weight_title,
            )

        ax.tick_params(labelsize=settings.font_size_ticks)

        if settings.enable_grid:
            ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

        plt.colorbar(cf, ax=ax, label="Free energy (kcal/mol)")

    # Always save individual plots for all resolved pairs
    for pc_x, pc_y in component_pairs:
        fig, ax = plt.subplots(
            figsize=(settings.fig_size_width, settings.fig_size_height),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
        )

        _draw(ax, surfaces[(pc_x, pc_y)], pc_x, pc_y)

        if settings.tight_layout:
            fig.tight_layout()

        out_path = out_dir / f"{settings.fig_basename}_pc{pc_x}_vs_pc{pc_y}.{settings.fig_format}"
        fig.savefig(
            out_path,
            dpi=settings.fig_dpi,
            transparent=settings.fig_transparent,
        )
        plt.close(fig)

        logger.info(f"Saved PCA FES heatmap → {out_path}")

    # Optionally also save one combined figure
    if settings.plot_in_one:
        fig, axes = plt.subplots(
            nrows=len(component_pairs),
            ncols=1,
            figsize=(settings.fig_size_width, settings.fig_size_height * len(component_pairs)),
            dpi=settings.fig_dpi,
            facecolor=settings.bg_color,
            squeeze=False,
        )

        for i, (pc_x, pc_y) in enumerate(component_pairs):
            _draw(axes[i][0], surfaces[(pc_x, pc_y)], pc_x, pc_y)

        if settings.tight_layout:
            fig.tight_layout()

        out_path = out_dir / f"{settings.fig_basename}.{settings.fig_format}"
        fig.savefig(
            out_path,
            dpi=settings.fig_dpi,
            transparent=settings.fig_transparent,
        )
        plt.close(fig)

        logger.info(f"Saved combined PCA FES heatmap → {out_path}")


def export_pca_summary(pta_file,
                       *,
                       group_name: str = "pca",
                       out_dir: Path,
                       landscape_components: tuple[int, int] = (1, 2),
                       bins: int = 120,
                       smooth_sigma: float = 1.0,
                       temperature: float = 310.0,
                       top_n_global: int = 10,
                       top_n_per_basin: int = 5,
                       max_basins: int = 10,
                       min_basin_population_fraction: float = 0.01,
                       min_basin_population_count: int = 10,
                       export_json: bool = True,
                       export_frame_scores_csv: bool = True,
                       export_basins_csv: bool = True,
                       export_representatives_csv: bool = True,
                       export_extremes_csv: bool = True,
                       export_frame_lists_txt: bool = True,
                       json_filename: str = "pca_summary.json",
                       frame_scores_csv_filename: str = "pca_frame_scores.csv",
                       basins_csv_filename: str = "pca_basins.csv",
                       representatives_csv_filename: str = "pca_representatives.csv",
                       extremes_csv_filename: str = "pca_extremes.csv",
                       representative_frames_txt_filename: str = "pca_representative_frames.txt",
                       extreme_frames_txt_filename: str = "pca_extreme_frames.txt") -> None:
    """
    Exports a PCA summary with detailed statistical analysis and visualization of principal components. The
    summary includes information about the PCA landscape, basins, frame scores, representatives, and
    extremes. Additionally, visualizations can be generated for easier interpretation.

    :param pta_file: The input PTA file containing PCA data.
    :param group_name: Name of the group in the PTA file where PCA data resides. Defaults to "pca".
    :param out_dir: Output directory where all the resulting files will be saved.
    :param landscape_components: Tuple representing the principal components for the landscape (e.g., (1, 2)).
    :param bins: Number of bins for visualizing the PCA density.
    :param smooth_sigma: Sigma value used for density smoothing. Must be non-negative.
    :param temperature: The temperature value (in Kelvin) used in free energy calculations.
    :param top_n_global: Number of top frames globally to include in the analysis.
    :param top_n_per_basin: Number of top frames per basin to include in further calculations.
    :param max_basins: Maximum number of basins to consider during the PCA analysis.
    :param min_basin_population_fraction: Minimum fraction of total population required for a basin to be considered.
    :param min_basin_population_count: Minimum absolute count of data points required for a basin to be considered.
    :param export_json: Whether to export the summary data as a JSON file.
    :param export_frame_scores_csv: Whether to export frame scores into a CSV file.
    :param export_basins_csv: Whether to export basin-related statistics into a CSV file.
    :param export_representatives_csv: Whether to export representative samples per basin into a CSV file.
    :param export_extremes_csv: Whether to export extremely low and high values into a CSV file.
    :param export_frame_lists_txt: Whether to export relevant frame indices into text files.
    :param json_filename: Name of the JSON file to save PCA summary data.
    :param frame_scores_csv_filename: Name of the CSV file to save frame scores.
    :param basins_csv_filename: Name of the CSV file to save basin statistics.
    :param representatives_csv_filename: Name of the CSV file to save representative frames.
    :param extremes_csv_filename: Name of the CSV file to save extreme frames.
    :param representative_frames_txt_filename: Name of the text file to save representative frames list.
    :param extreme_frames_txt_filename: Name of the text file to save extreme frames list.
    :return: None
    """
    if not pta_file.group_exists(group_name):
        raise ValueError(f"Group '{group_name}' not found in PTA file")

    if bins < 10:
        raise ValueError("bins must be >= 10")
    if smooth_sigma < 0:
        raise ValueError("smooth_sigma must be >= 0")
    if temperature <= 0:
        raise ValueError("temperature must be > 0")
    if top_n_global < 1:
        raise ValueError("top_n_global must be >= 1")
    if top_n_per_basin < 1:
        raise ValueError("top_n_per_basin must be >= 1")
    if max_basins < 1:
        raise ValueError("max_basins must be >= 1")

    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    pc_x, pc_y = int(landscape_components[0]), int(landscape_components[1])

    # Helpers
    def _safe_float(x, default: float = 0.0) -> float:
        """
        Convert a value to a floating-point number safely.

        This function attempts to convert the input value to a float. If the
        conversion fails due to an exception, a provided default value is returned.

        :param x: The value to be converted into a float.
        :type x: Any
        :param default: The default floating-point value to return in case the
            conversion fails.
        :return: The floating-point representation of the input value, or the
            default value if the conversion fails.
        :rtype: float
        """
        try:
            return float(x)
        except Exception:
            return float(default)

    def _nearest_center(x0: float, y0: float, centers: list[tuple[float, float, int]]) -> tuple[int | None, float]:
        """
        Finds the nearest center to a given point (x0, y0) from a list of centers and calculates
        the distance to that center.

        :param x0: X coordinate of the point
        :param y0: Y coordinate of the point
        :param centers: A list of tuples where each tuple contains the center's X coordinate,
            Y coordinate, and associated basin ID
        :return: A tuple containing the ID of the nearest center (or None if no centers are
            provided) and the distance to the nearest center
        """
        if not centers:
            return None, float("nan")
        best_id = None
        best_d = float("inf")
        for cx, cy, basin_id in centers:
            d = math.hypot(x0 - cx, y0 - cy)
            if d < best_d:
                best_d = d
                best_id = basin_id
        return best_id, best_d

    def _first_pc_index_for_threshold(cum_percentages: np.ndarray, threshold: float) -> int:
        """
        Calculates the index for the first principal component that meets or exceeds a given
        threshold in cumulative percentages. This function determines the smallest index in
        the cumulative percentage array such that the value at that index is greater than or
        equal to the given threshold, adjusted to a 1-based index.

        :param cum_percentages: A 1D numpy array representing cumulative percentages of
            explained variance. Should be sorted in ascending order.
        :param threshold: A float representing the threshold value to determine the index.
        :return: The 1-based index (int) that corresponds to the first principal component
            meeting or exceeding the threshold.
        """
        idx = np.searchsorted(cum_percentages, threshold, side="left")
        idx = min(idx, len(cum_percentages) - 1)
        return int(idx + 1)

    def _top_rows(rows: list[dict], key: str, n: int, reverse: bool = True, finite_only: bool = False) -> list[dict]:
        """
        Filters and ranks the top `n` rows from a list of dictionaries based on the value of a
        specified key. Optionally allows filtering for finite values only and can sort in
        ascending or descending order.

        :param rows:
            A list of dictionaries where each dictionary represents a row of data.
        :param key:
            The key in each dictionary whose value will be used for ranking.
        :param n:
            The maximum number of rows to return.
        :param reverse:
            If True, sorts the rows in descending order. If False, sorts in ascending order.
            Defaults to True.
        :param finite_only:
            If True, rows with non-finite values (e.g., NaN, inf) in the specified key will
            be excluded. Defaults to False.
        :return:
            A list of dictionaries containing the top `n` rows, each with an additional key
            "rank" indicating the rank of the row.
        """
        data = rows
        if finite_only:
            data = [r for r in rows if np.isfinite(r.get(key, float("nan")))]
        ranked = sorted(data, key=lambda r: r.get(key, float("-inf") if reverse else float("inf")), reverse=reverse)
        out = []
        for i, rec in enumerate(ranked[:n], start=1):
            row = dict(rec)
            row["rank"] = i
            out.append(row)
        return out

    def _frame_record_subset(rec: dict, basin_id_key: str = "basin_id") -> dict:
        """
        Processes a given record dictionary and returns a subset of the record containing
        selected keys transformed to appropriate types. The function handles cases where
        certain keys may contain non-finite values or be missing, assigning `None` to such cases.

        :param rec: The input dictionary representing the record with numerical and other
            metadata.
        :type rec: dict
        :param basin_id_key: The key in the record dictionary used to extract the `basin_id`
            value. Defaults to "basin_id".
        :type basin_id_key: str
        :return: A subset of the record with transformed fields for selected keys. Non-finite
            or missing values are replaced with `None`.
        :rtype: dict
        """
        return {
            "frame_index": int(rec["frame_index"]),
            "time_ps": float(rec["time_ps"]),
            "time_ns": float(rec["time_ns"]),
            "time_us": float(rec["time_us"]),
            "pc_x": float(rec["pc_x"]),
            "pc_y": float(rec["pc_y"]),
            "free_energy": None if not np.isfinite(rec["free_energy"]) else float(rec["free_energy"]),
            "probability_density": None if not np.isfinite(rec["probability_density"]) else float(rec["probability_density"]),
            "basin_id": None if rec.get(basin_id_key) is None else int(rec[basin_id_key]),
        }

    # Read PCA table
    frames = list(pta_file._iter_frames(group_name))
    if not frames:
        raise ValueError(f"No frames found under group '{group_name}'")

    frame_rows: list[dict] = []
    variances_by_pc: dict[int, float] = {}
    ratios_by_pc: dict[int, float] = {}
    projection_values_by_pc: dict[int, list[float]] = {}

    for frame in frames:
        dset_path = f"{group_name}/frame_{frame}/pca"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        attrs = dset.attrs

        time_ps = _safe_float(attrs.get("time_ps", frame))
        time_ns = _safe_float(attrs.get("time_ns", time_ps / 1_000.0))
        time_us = _safe_float(attrs.get("time_us", time_ps / 1_000_000.0))

        row = {
            "frame_index": int(frame),
            "time_ps": time_ps,
            "time_ns": time_ns,
            "time_us": time_us,
        }

        for raw in dset:
            rec = json.loads(raw)
            pc = int(rec["pc"])
            value = float(rec["value"])
            variance = float(rec["variance"])
            variance_ratio = float(rec["variance_ratio"])

            row[pc] = value
            variances_by_pc[pc] = variance
            ratios_by_pc[pc] = variance_ratio
            projection_values_by_pc.setdefault(pc, []).append(value)

        frame_rows.append(row)

    if not frame_rows:
        raise ValueError("No PCA rows collected")

    all_pcs = sorted(variances_by_pc.keys())
    if pc_x not in all_pcs or pc_y not in all_pcs:
        raise ValueError(
            f"Landscape components {landscape_components} not present in PCA data. "
            f"Available PCs: {all_pcs}"
        )

    # Global metadata
    n_frames = len(frame_rows)
    n_components = len(all_pcs)

    time_ps_values = np.asarray([float(r["time_ps"]) for r in frame_rows], dtype=float)
    time_ps_sorted = np.sort(time_ps_values)
    if len(time_ps_sorted) >= 2:
        diffs = np.diff(time_ps_sorted)
        positive_diffs = diffs[diffs > 0]
        time_step_ps_estimate = float(np.median(positive_diffs)) if positive_diffs.size else 0.0
    else:
        time_step_ps_estimate = 0.0

    pcs_arr = np.asarray(all_pcs, dtype=int)
    variances_arr = np.asarray([variances_by_pc[pc] for pc in all_pcs], dtype=float)
    ratios_arr = np.asarray([ratios_by_pc[pc] for pc in all_pcs], dtype=float)
    percentages_arr = ratios_arr * 100.0
    cumulative_percentages_arr = np.cumsum(percentages_arr)

    # Effective dimensionality (entropy-based)
    safe_ratios = ratios_arr[ratios_arr > 0]
    if safe_ratios.size:
        entropy = -float(np.sum(safe_ratios * np.log(safe_ratios)))
        effective_dimensionality = float(np.exp(entropy))
    else:
        effective_dimensionality = 0.0

    # Per-component stats
    pca_components_summary: list[dict] = []
    for i, pc in enumerate(all_pcs):
        vals = np.asarray(projection_values_by_pc.get(pc, []), dtype=float)
        if vals.size == 0:
            continue
        q1, q3 = np.percentile(vals, [25, 75])
        pca_components_summary.append(
            {
                "pc": int(pc),
                "variance": float(variances_arr[i]),
                "variance_ratio": float(ratios_arr[i]),
                "percentage": float(percentages_arr[i]),
                "cumulative_percentage": float(cumulative_percentages_arr[i]),
                "projection_min": float(np.min(vals)),
                "projection_max": float(np.max(vals)),
                "projection_mean": float(np.mean(vals)),
                "projection_std": float(np.std(vals)),
                "projection_median": float(np.median(vals)),
                "projection_iqr": float(q3 - q1),
            }
        )

    coverage_thresholds = {
        "50_percent": _first_pc_index_for_threshold(cumulative_percentages_arr, 50.0),
        "75_percent": _first_pc_index_for_threshold(cumulative_percentages_arr, 75.0),
        "90_percent": _first_pc_index_for_threshold(cumulative_percentages_arr, 90.0),
        "95_percent": _first_pc_index_for_threshold(cumulative_percentages_arr, 95.0),
    }

    # Build PCx / PCy landscape
    landscape_rows: list[dict] = []
    x_vals = []
    y_vals = []

    for r in frame_rows:
        if pc_x in r and pc_y in r:
            row = {
                "frame_index": int(r["frame_index"]),
                "time_ps": float(r["time_ps"]),
                "time_ns": float(r["time_ns"]),
                "time_us": float(r["time_us"]),
                "pc_x": float(r[pc_x]),
                "pc_y": float(r[pc_y]),
            }
            landscape_rows.append(row)
            x_vals.append(float(r[pc_x]))
            y_vals.append(float(r[pc_y]))

    if not landscape_rows:
        raise ValueError(f"No frames contain both PC{pc_x} and PC{pc_y}")

    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)

    H, xedges, yedges = np.histogram2d(x, y, bins=bins, density=True)

    if smooth_sigma > 0:
        H = gaussian_filter(H, sigma=smooth_sigma)
        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        norm = np.nansum(H) * dx * dy
        if norm > 0:
            H = H / norm

    H_safe = H.copy()
    H_safe[H_safe <= 0] = np.nan

    R = 0.0019872041  # kcal/mol/K
    F = -R * temperature * np.log(H_safe)
    F = F - np.nanmin(F)

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])

    # Annotate frames with density / free energy by bin lookup
    for row in landscape_rows:
        xi = int(np.searchsorted(xedges, row["pc_x"], side="right") - 1)
        yi = int(np.searchsorted(yedges, row["pc_y"], side="right") - 1)

        xi = max(0, min(xi, H.shape[0] - 1))
        yi = max(0, min(yi, H.shape[1] - 1))

        density = float(H[xi, yi])
        free_energy = float(F[xi, yi]) if np.isfinite(F[xi, yi]) else float("nan")

        row["probability_density"] = density
        row["free_energy"] = free_energy
        row["bin_x"] = xi
        row["bin_y"] = yi

    # Detect basins on the FES grid
    basin_candidates: list[dict] = []
    F_work = F.copy()
    finite_mask = np.isfinite(F_work)

    if np.any(finite_mask):
        local_min_mask = (F_work == minimum_filter(F_work, size=3, mode="nearest")) & finite_mask

        candidate_indices = np.argwhere(local_min_mask)
        for bx, by in candidate_indices:
            basin_candidates.append(
                {
                    "bin_x": int(bx),
                    "bin_y": int(by),
                    "center_pc_x": float(xcenters[bx]),
                    "center_pc_y": float(ycenters[by]),
                    "minimum_free_energy": float(F_work[bx, by]),
                }
            )

        basin_candidates = sorted(
            basin_candidates,
            key=lambda r: r["minimum_free_energy"]
        )[: max_basins * 5]  # extra headroom before population filtering

    # Assign frames to nearest candidate basin
    candidate_centers = [
        (b["center_pc_x"], b["center_pc_y"], i + 1)
        for i, b in enumerate(basin_candidates)
    ]

    for row in landscape_rows:
        basin_id, nearest_dist = _nearest_center(
            row["pc_x"],
            row["pc_y"],
            candidate_centers,
        )
        row["nearest_basin_id_raw"] = basin_id
        row["distance_to_nearest_basin_center_raw"] = nearest_dist

    # Population counts per raw basin
    raw_counts: dict[int, int] = {}
    for row in landscape_rows:
        bid = row.get("nearest_basin_id_raw")
        if bid is not None:
            raw_counts[bid] = raw_counts.get(bid, 0) + 1

    retained_raw_basin_ids: list[int] = []
    for raw_id, count in sorted(raw_counts.items(), key=lambda kv: kv[1], reverse=True):
        frac = count / float(len(landscape_rows))
        if count >= min_basin_population_count and frac >= min_basin_population_fraction:
            retained_raw_basin_ids.append(int(raw_id))

    retained_raw_basin_ids = retained_raw_basin_ids[:max_basins]

    # Remap retained basins to contiguous ids
    raw_to_final_basin_id = {raw_id: i + 1 for i, raw_id in enumerate(retained_raw_basin_ids)}

    final_basin_candidates: list[dict] = []
    for raw_id in retained_raw_basin_ids:
        b = basin_candidates[raw_id - 1]
        final_basin_candidates.append(
            {
                "basin_id": raw_to_final_basin_id[raw_id],
                "raw_basin_id": raw_id,
                "center_pc_x": b["center_pc_x"],
                "center_pc_y": b["center_pc_y"],
                "minimum_free_energy": b["minimum_free_energy"],
            }
        )

    final_centers = [
        (b["center_pc_x"], b["center_pc_y"], b["basin_id"])
        for b in final_basin_candidates
    ]

    # Final per-frame basin assignment
    for row in landscape_rows:
        if final_centers:
            basin_id, nearest_dist = _nearest_center(
                row["pc_x"],
                row["pc_y"],
                final_centers,
            )
        else:
            basin_id, nearest_dist = (None, float("nan"))

        row["basin_id"] = basin_id
        row["distance_to_nearest_basin_center"] = nearest_dist
        row["distance_to_basin_center"] = nearest_dist if basin_id is not None else float("nan")

    # Global frame rankings
    finite_fe_rows = [r for r in landscape_rows if np.isfinite(r["free_energy"])]
    finite_pd_rows = [r for r in landscape_rows if np.isfinite(r["probability_density"])]

    top_lowest_energy_frames = _top_rows(finite_fe_rows, "free_energy", top_n_global, reverse=False, finite_only=True)
    top_highest_density_frames = _top_rows(finite_pd_rows, "probability_density", top_n_global, reverse=True, finite_only=True)
    top_highest_free_energy_frames = _top_rows(finite_fe_rows, "free_energy", top_n_global, reverse=True, finite_only=True)
    top_lowest_density_frames = _top_rows(finite_pd_rows, "probability_density", top_n_global, reverse=False, finite_only=True)
    top_largest_abs_pc_x_frames = _top_rows(
        [{**r, "abs_pc_x": abs(r["pc_x"])} for r in landscape_rows],
        "abs_pc_x",
        top_n_global,
        reverse=True,
        finite_only=True,
    )
    top_largest_abs_pc_y_frames = _top_rows(
        [{**r, "abs_pc_y": abs(r["pc_y"])} for r in landscape_rows],
        "abs_pc_y",
        top_n_global,
        reverse=True,
        finite_only=True,
    )
    top_farthest_from_nearest_basin_frames = _top_rows(
        [r for r in landscape_rows if np.isfinite(r["distance_to_nearest_basin_center"])],
        "distance_to_nearest_basin_center",
        top_n_global,
        reverse=True,
        finite_only=True,
    )

    global_minimum_frame = top_lowest_energy_frames[0] if top_lowest_energy_frames else None
    highest_density_frame = top_highest_density_frames[0] if top_highest_density_frames else None

    # Basin summaries
    basin_summaries: list[dict] = []
    basin_representatives: list[int] = []
    basin_lowest_energy_frames: list[int] = []
    basin_highest_density_frames: list[int] = []
    basin_farthest_frames: list[int] = []

    for basin in final_basin_candidates:
        basin_id = int(basin["basin_id"])
        cx = float(basin["center_pc_x"])
        cy = float(basin["center_pc_y"])

        members = [r for r in landscape_rows if r.get("basin_id") == basin_id]
        if not members:
            continue

        for m in members:
            m["distance_to_basin_center"] = math.hypot(m["pc_x"] - cx, m["pc_y"] - cy)

        population_count = len(members)
        population_fraction = population_count / float(len(landscape_rows))

        representative = min(members, key=lambda r: r["distance_to_basin_center"])
        members_finite_fe = [r for r in members if np.isfinite(r["free_energy"])]
        members_finite_pd = [r for r in members if np.isfinite(r["probability_density"])]

        lowest_energy = min(members_finite_fe, key=lambda r: r["free_energy"]) if members_finite_fe else representative
        highest_density = max(members_finite_pd, key=lambda r: r["probability_density"]) if members_finite_pd else representative
        farthest = max(members, key=lambda r: r["distance_to_basin_center"])

        distances = np.asarray([m["distance_to_basin_center"] for m in members], dtype=float)
        radius_estimate = float(np.percentile(distances, 90)) if distances.size else 0.0

        top_frames_basin = [
            int(r["frame_index"])
            for r in sorted(members_finite_fe, key=lambda r: r["free_energy"])[:top_n_per_basin]
        ]

        basin_summary = {
            "basin_id": basin_id,
            "center_pc_x": cx,
            "center_pc_y": cy,
            "minimum_free_energy": float(np.nanmin([m["free_energy"] for m in members_finite_fe])) if members_finite_fe else float("nan"),
            "mean_free_energy": float(np.nanmean([m["free_energy"] for m in members_finite_fe])) if members_finite_fe else float("nan"),
            "population_count": int(population_count),
            "population_fraction": float(population_fraction),
            "representative_frame": int(representative["frame_index"]),
            "representative_time_ps": float(representative["time_ps"]),
            "representative_time_ns": float(representative["time_ns"]),
            "representative_time_us": float(representative["time_us"]),
            "representative_pc_x": float(representative["pc_x"]),
            "representative_pc_y": float(representative["pc_y"]),
            "lowest_energy_frame": int(lowest_energy["frame_index"]),
            "highest_density_frame": int(highest_density["frame_index"]),
            "farthest_from_center_frame": int(farthest["frame_index"]),
            "farthest_from_center_distance": float(farthest["distance_to_basin_center"]),
            "radius_estimate": radius_estimate,
            "top_frames": top_frames_basin,
        }
        basin_summaries.append(basin_summary)

        basin_representatives.append(int(representative["frame_index"]))
        basin_lowest_energy_frames.append(int(lowest_energy["frame_index"]))
        basin_highest_density_frames.append(int(highest_density["frame_index"]))
        basin_farthest_frames.append(int(farthest["frame_index"]))

    # Flag per-frame annotations
    global_min_frame_idx = int(global_minimum_frame["frame_index"]) if global_minimum_frame else None
    highest_density_frame_idx = int(highest_density_frame["frame_index"]) if highest_density_frame else None

    extreme_high_energy_ids = {int(r["frame_index"]) for r in top_highest_free_energy_frames}
    extreme_low_density_ids = {int(r["frame_index"]) for r in top_lowest_density_frames}
    extreme_far_from_wells_ids = {int(r["frame_index"]) for r in top_farthest_from_nearest_basin_frames}
    basin_rep_ids = set(basin_representatives)
    basin_low_ids = set(basin_lowest_energy_frames)
    basin_dense_ids = set(basin_highest_density_frames)
    basin_far_ids = set(basin_farthest_frames)

    for row in landscape_rows:
        fid = int(row["frame_index"])
        row["is_global_minimum"] = fid == global_min_frame_idx
        row["is_highest_density"] = fid == highest_density_frame_idx
        row["is_basin_representative"] = fid in basin_rep_ids
        row["is_basin_lowest_energy"] = fid in basin_low_ids
        row["is_basin_highest_density"] = fid in basin_dense_ids
        row["is_basin_farthest"] = fid in basin_far_ids
        row["is_extreme_high_energy"] = fid in extreme_high_energy_ids
        row["is_extreme_low_density"] = fid in extreme_low_density_ids
        row["is_extreme_far_from_wells"] = fid in extreme_far_from_wells_ids

    # Build final JSON payload
    summary = {
        "metadata": {
            "group_name": str(group_name),
            "n_frames": int(n_frames),
            "n_components": int(n_components),
            "landscape_components": [int(pc_x), int(pc_y)],
            "bins": int(bins),
            "smooth_sigma": float(smooth_sigma),
            "temperature": float(temperature),
            "time_start_ps": float(np.min(time_ps_values)),
            "time_end_ps": float(np.max(time_ps_values)),
            "time_step_ps_estimate": float(time_step_ps_estimate),
        },
        "pca": {
            "effective_dimensionality": float(effective_dimensionality),
            "dominant_component": int(all_pcs[int(np.argmax(ratios_arr))]) if ratios_arr.size else None,
            "coverage_thresholds": coverage_thresholds,
            "components": pca_components_summary,
        },
        "landscape": {
            "components_used": [int(pc_x), int(pc_y)],
            "n_detected_basins_raw": int(len(basin_candidates)),
            "n_retained_basins": int(len(basin_summaries)),
            "global_minimum_frame": _frame_record_subset(global_minimum_frame) if global_minimum_frame else None,
            "highest_density_frame": _frame_record_subset(highest_density_frame) if highest_density_frame else None,
            "top_lowest_energy_frames": [_frame_record_subset(r) for r in top_lowest_energy_frames],
            "top_highest_density_frames": [_frame_record_subset(r) for r in top_highest_density_frames],
        },
        "basins": basin_summaries,
        "representatives": {
            "global_minimum_frame": global_min_frame_idx,
            "highest_density_frame": highest_density_frame_idx,
            "basin_representatives": basin_representatives,
            "basin_lowest_energy_frames": basin_lowest_energy_frames,
            "basin_highest_density_frames": basin_highest_density_frames,
            "basin_farthest_frames": basin_farthest_frames,
        },
        "extremes": {
            "highest_free_energy_frames": [
                {
                    "rank": int(r["rank"]),
                    **_frame_record_subset(r),
                }
                for r in top_highest_free_energy_frames
            ],
            "lowest_density_frames": [
                {
                    "rank": int(r["rank"]),
                    **_frame_record_subset(r),
                }
                for r in top_lowest_density_frames
            ],
            "farthest_from_nearest_basin_center_frames": [
                {
                    "rank": int(r["rank"]),
                    "frame_index": int(r["frame_index"]),
                    "time_ps": float(r["time_ps"]),
                    "time_ns": float(r["time_ns"]),
                    "time_us": float(r["time_us"]),
                    "pc_x": float(r["pc_x"]),
                    "pc_y": float(r["pc_y"]),
                    "distance_to_nearest_basin_center": float(r["distance_to_nearest_basin_center"]),
                    "nearest_basin_id": None if r["basin_id"] is None else int(r["basin_id"]),
                }
                for r in top_farthest_from_nearest_basin_frames
            ],
            "largest_absolute_pc_x_frames": [
                {
                    "rank": int(r["rank"]),
                    "frame_index": int(r["frame_index"]),
                    "time_ps": float(r["time_ps"]),
                    "time_ns": float(r["time_ns"]),
                    "time_us": float(r["time_us"]),
                    "pc_x": float(r["pc_x"]),
                }
                for r in top_largest_abs_pc_x_frames
            ],
            "largest_absolute_pc_y_frames": [
                {
                    "rank": int(r["rank"]),
                    "frame_index": int(r["frame_index"]),
                    "time_ps": float(r["time_ps"]),
                    "time_ns": float(r["time_ns"]),
                    "time_us": float(r["time_us"]),
                    "pc_y": float(r["pc_y"]),
                }
                for r in top_largest_abs_pc_y_frames
            ],
        },
        "files": {
            "frame_scores_csv": frame_scores_csv_filename if export_frame_scores_csv else None,
            "basins_csv": basins_csv_filename if export_basins_csv else None,
            "representatives_csv": representatives_csv_filename if export_representatives_csv else None,
            "extremes_csv": extremes_csv_filename if export_extremes_csv else None,
            "representative_frames_txt": representative_frames_txt_filename if export_frame_lists_txt else None,
            "extreme_frames_txt": extreme_frames_txt_filename if export_frame_lists_txt else None,
        },
    }

    # Write JSON
    if export_json:
        json_path = out_dir / json_filename
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=4)
        logger.info(f"Saved PCA summary JSON → {json_path}")

    # Write frame scores CSV
    if export_frame_scores_csv:
        csv_path = out_dir / frame_scores_csv_filename
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "frame_index",
                "time_ps",
                "time_ns",
                "time_us",
                "pc_x",
                "pc_y",
                "probability_density",
                "free_energy",
                "basin_id",
                "distance_to_basin_center",
                "distance_to_nearest_basin_center",
                "is_global_minimum",
                "is_highest_density",
                "is_basin_representative",
                "is_basin_lowest_energy",
                "is_basin_highest_density",
                "is_basin_farthest",
                "is_extreme_high_energy",
                "is_extreme_low_density",
                "is_extreme_far_from_wells",
            ])
            for r in sorted(landscape_rows, key=lambda x: x["frame_index"]):
                writer.writerow([
                    int(r["frame_index"]),
                    f"{float(r['time_ps']):.6f}",
                    f"{float(r['time_ns']):.6f}",
                    f"{float(r['time_us']):.6f}",
                    f"{float(r['pc_x']):.6f}",
                    f"{float(r['pc_y']):.6f}",
                    "" if not np.isfinite(r["probability_density"]) else f"{float(r['probability_density']):.12g}",
                    "" if not np.isfinite(r["free_energy"]) else f"{float(r['free_energy']):.12g}",
                    "" if r["basin_id"] is None else int(r["basin_id"]),
                    "" if not np.isfinite(r["distance_to_basin_center"]) else f"{float(r['distance_to_basin_center']):.12g}",
                    "" if not np.isfinite(r["distance_to_nearest_basin_center"]) else f"{float(r['distance_to_nearest_basin_center']):.12g}",
                    int(bool(r["is_global_minimum"])),
                    int(bool(r["is_highest_density"])),
                    int(bool(r["is_basin_representative"])),
                    int(bool(r["is_basin_lowest_energy"])),
                    int(bool(r["is_basin_highest_density"])),
                    int(bool(r["is_basin_farthest"])),
                    int(bool(r["is_extreme_high_energy"])),
                    int(bool(r["is_extreme_low_density"])),
                    int(bool(r["is_extreme_far_from_wells"])),
                ])
        logger.info(f"Saved PCA frame scores CSV → {csv_path}")

    # Write basins CSV
    if export_basins_csv:
        csv_path = out_dir / basins_csv_filename
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "basin_id",
                "center_pc_x",
                "center_pc_y",
                "minimum_free_energy",
                "mean_free_energy",
                "population_count",
                "population_fraction",
                "representative_frame",
                "lowest_energy_frame",
                "highest_density_frame",
                "farthest_from_center_frame",
                "farthest_from_center_distance",
                "radius_estimate",
            ])
            for b in basin_summaries:
                writer.writerow([
                    int(b["basin_id"]),
                    f"{float(b['center_pc_x']):.6f}",
                    f"{float(b['center_pc_y']):.6f}",
                    "" if not np.isfinite(b["minimum_free_energy"]) else f"{float(b['minimum_free_energy']):.12g}",
                    "" if not np.isfinite(b["mean_free_energy"]) else f"{float(b['mean_free_energy']):.12g}",
                    int(b["population_count"]),
                    f"{float(b['population_fraction']):.6f}",
                    int(b["representative_frame"]),
                    int(b["lowest_energy_frame"]),
                    int(b["highest_density_frame"]),
                    int(b["farthest_from_center_frame"]),
                    f"{float(b['farthest_from_center_distance']):.12g}",
                    f"{float(b['radius_estimate']):.12g}",
                ])
        logger.info(f"Saved PCA basins CSV → {csv_path}")

    # Write representatives CSV
    if export_representatives_csv:
        csv_path = out_dir / representatives_csv_filename
        rows_out = []

        if global_minimum_frame:
            rows_out.append(("global_minimum", None, global_minimum_frame))
        if highest_density_frame:
            rows_out.append(("highest_density", None, highest_density_frame))

        for b in basin_summaries:
            basin_id = int(b["basin_id"])
            basin_members = [r for r in landscape_rows if r.get("basin_id") == basin_id]
            if not basin_members:
                continue
            rep = next(r for r in basin_members if int(r["frame_index"]) == int(b["representative_frame"]))
            low = next(r for r in basin_members if int(r["frame_index"]) == int(b["lowest_energy_frame"]))
            dense = next(r for r in basin_members if int(r["frame_index"]) == int(b["highest_density_frame"]))
            far = next(r for r in basin_members if int(r["frame_index"]) == int(b["farthest_from_center_frame"]))

            rows_out.extend([
                ("basin_representative", basin_id, rep),
                ("basin_lowest_energy", basin_id, low),
                ("basin_highest_density", basin_id, dense),
                ("basin_farthest_from_center", basin_id, far),
            ])

        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "category",
                "basin_id",
                "frame_index",
                "time_ps",
                "time_ns",
                "time_us",
                "pc_x",
                "pc_y",
                "free_energy",
                "probability_density",
                "distance_to_basin_center",
            ])
            for category, basin_id, r in rows_out:
                writer.writerow([
                    category,
                    "" if basin_id is None else int(basin_id),
                    int(r["frame_index"]),
                    f"{float(r['time_ps']):.6f}",
                    f"{float(r['time_ns']):.6f}",
                    f"{float(r['time_us']):.6f}",
                    f"{float(r['pc_x']):.6f}",
                    f"{float(r['pc_y']):.6f}",
                    "" if not np.isfinite(r["free_energy"]) else f"{float(r['free_energy']):.12g}",
                    "" if not np.isfinite(r["probability_density"]) else f"{float(r['probability_density']):.12g}",
                    "" if not np.isfinite(r["distance_to_basin_center"]) else f"{float(r['distance_to_basin_center']):.12g}",
                ])
        logger.info(f"Saved PCA representatives CSV → {csv_path}")

    # Write extremes CSV
    if export_extremes_csv:
        csv_path = out_dir / extremes_csv_filename
        rows_out = []

        for r in top_highest_free_energy_frames:
            rows_out.append(("highest_free_energy", r))
        for r in top_lowest_density_frames:
            rows_out.append(("lowest_density", r))
        for r in top_farthest_from_nearest_basin_frames:
            rows_out.append(("farthest_from_nearest_basin", r))
        for r in top_largest_abs_pc_x_frames:
            rows_out.append(("largest_absolute_pc_x", r))
        for r in top_largest_abs_pc_y_frames:
            rows_out.append(("largest_absolute_pc_y", r))

        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "category",
                "rank",
                "frame_index",
                "time_ps",
                "time_ns",
                "time_us",
                "pc_x",
                "pc_y",
                "free_energy",
                "probability_density",
                "basin_id",
                "distance_to_nearest_basin_center",
            ])
            for category, r in rows_out:
                writer.writerow([
                    category,
                    int(r["rank"]),
                    int(r["frame_index"]),
                    f"{float(r['time_ps']):.6f}",
                    f"{float(r['time_ns']):.6f}",
                    f"{float(r['time_us']):.6f}",
                    f"{float(r.get('pc_x', float('nan'))):.6f}" if "pc_x" in r else "",
                    f"{float(r.get('pc_y', float('nan'))):.6f}" if "pc_y" in r else "",
                    "" if not np.isfinite(r.get("free_energy", float("nan"))) else f"{float(r['free_energy']):.12g}",
                    "" if not np.isfinite(r.get("probability_density", float("nan"))) else f"{float(r['probability_density']):.12g}",
                    "" if r.get("basin_id") is None else int(r["basin_id"]),
                    "" if not np.isfinite(r.get("distance_to_nearest_basin_center", float("nan"))) else f"{float(r['distance_to_nearest_basin_center']):.12g}",
                ])
        logger.info(f"Saved PCA extremes CSV → {csv_path}")

    # Write frame lists TXT
    if export_frame_lists_txt:
        rep_txt_path = out_dir / representative_frames_txt_filename
        ext_txt_path = out_dir / extreme_frames_txt_filename

        rep_ids = []
        if global_minimum_frame:
            rep_ids.append(int(global_minimum_frame["frame_index"]))
        if highest_density_frame:
            rep_ids.append(int(highest_density_frame["frame_index"]))
        rep_ids.extend(basin_representatives)
        rep_ids.extend(basin_lowest_energy_frames)
        rep_ids.extend(basin_highest_density_frames)
        rep_ids.extend(basin_farthest_frames)
        rep_ids = sorted(set(rep_ids))

        ext_ids = []
        ext_ids.extend([int(r["frame_index"]) for r in top_highest_free_energy_frames])
        ext_ids.extend([int(r["frame_index"]) for r in top_lowest_density_frames])
        ext_ids.extend([int(r["frame_index"]) for r in top_farthest_from_nearest_basin_frames])
        ext_ids.extend([int(r["frame_index"]) for r in top_largest_abs_pc_x_frames])
        ext_ids.extend([int(r["frame_index"]) for r in top_largest_abs_pc_y_frames])
        ext_ids = sorted(set(ext_ids))

        with open(rep_txt_path, "w", encoding="utf-8") as fh:
            for fid in rep_ids:
                fh.write(f"{fid}\n")

        with open(ext_txt_path, "w", encoding="utf-8") as fh:
            for fid in ext_ids:
                fh.write(f"{fid}\n")

        logger.info(f"Saved PCA representative frame list → {rep_txt_path}")
        logger.info(f"Saved PCA extreme frame list → {ext_txt_path}")
