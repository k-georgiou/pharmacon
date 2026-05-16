"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Module :mod:`pharmacon.plotter.interactions`.
"""
import re
import math
import json

import warnings


import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import ast
from typing import Dict, Tuple


from typing_extensions import Final, Set, List
from matplotlib.colors import is_color_like

import networkx as nx
from pharmacon.logger import get_logger, PharmaconLogger

from pharmacon.constants import (PLIStackedSettings1,  PLIStackedSettings2, PLIHeatmapSettings1, PLIHeatmapSettings2,
                                 PLIPieChartsSettings1, PLILigandMonitorSettings,

                                 PPITimelinePairsSettings, PPIHeatmapSettings, PPIStackedColumnSettings,

                                 AA3_to_AA1)


# Canonical bottom->top stack order for PLI stacked-column plots.
# Labels are the normalized lowercase forms produced by `norm_interaction`
# inside the build_pli_* functions.
_PLI_STACK_ORDER: Final[Tuple[str, ...]] = (
    "hydrophobic",
    "hydrogen_bonds",
    "pi_cation",
    "pi_stacking",
    "water_bridge",
    "ionic",
    "halogen",
    "metal_contact",
)
_PLI_STACK_RANK: Final[Dict[str, int]] = {k: i for i, k in enumerate(_PLI_STACK_ORDER)}







__all__ = [
    "parse_range_dict",
    "logger",
    "plot_protein_ligand_interactions_stacked_column_1_from_file",
    "plot_protein_ligand_interactions_stacked_column_2_from_file",
    "plot_protein_ligand_interactions_heatmap_1_from_file",
    "plot_protein_ligand_interactions_heatmap_2_from_file",
    "plot_protein_ligand_interactions_pie_charts_from_file",
    "plot_protein_ligand_interactions_ligand_monitor_from_file",
    "parse_frame_interaction_record",
    "build_pli_merged_stacked_data",
    "build_pli_normal_data",
    "plot_protein_protein_timeline_pairs_from_file",
    "plot_protein_protein_heatmap_freq_from_file",
    "plot_protein_protein_interactions_stacked_column_from_file",
]


def parse_range_dict(obj) -> Dict[str, Tuple[int, int]]:
    """
    Parses an input object representing a range mapping into a dictionary where keys
    are strings and values are tuples of two integers.

    Detailed behavior:
    - If the input is already a dictionary, it is used directly.
    - If the input is a string, it is evaluated as a dictionary using `ast.literal_eval`.
    - The resulting dictionary's keys must be strings, and its values must be tuples of
      length 2 containing non-negative integers, where the first integer is less than
      or equal to the second.

    Raises detailed exceptions if the input doesn't conform to these expectations.

    :param obj: The input object to parse. Can be a string representation of a
                dictionary or a dictionary itself.
    :type obj: Union[str, dict]
    :return: A dictionary mapping string keys to tuples of two integers.
    :rtype: Dict[str, Tuple[int, int]]
    :raises ValueError: If the input is not a valid dictionary or its structure does
                        not align with the described format.
    """

    # If already dict (because literal_eval parsed it), use directly
    if isinstance(obj, dict):
        raw = obj
    elif isinstance(obj, str):
        obj = obj.strip()
        if not obj:
            return {}
        try:
            raw = ast.literal_eval(obj)
        except Exception as e:
            raise ValueError(f"Invalid range mapping syntax: {obj}") from e
    else:
        raise ValueError("Alter mapping must be dict or string.")

    if not isinstance(raw, dict):
        raise ValueError("Alter mapping must be a dict.")

    result: Dict[str, Tuple[int, int]] = {}

    for key, value in raw.items():

        if not isinstance(key, str):
            raise ValueError(f"Invalid key type: {key} (must be string)")

        if not isinstance(value, (tuple, list)) or len(value) != 2:
            raise ValueError(
                f"Invalid value for '{key}': {value} (must be tuple of length 2)"
            )

        lo = int(value[0])
        hi = int(value[1])

        if lo < 0 or hi < 0:
            raise ValueError(
                f"Range values must be non-negative for '{key}': {value}"
            )

        if lo > hi:
            raise ValueError(
                f"Lower bound > upper bound for '{key}': {value}"
            )

        result[key] = (lo, hi)

    return result


warnings.filterwarnings("ignore")


logger: PharmaconLogger = get_logger(__name__)


def plot_protein_ligand_interactions_stacked_column_1_from_file(pta_file, *,
                                                                group_name: str,
                                                                mode_name: str,
                                                                settings: PLIStackedSettings1,
                                                                out_dir: Path = Path("."),
                                                                attach_to_name: str = "mode",
                                                                is_merged: bool = False) -> None:
    """
    Generates stacked column bar plots representing protein-ligand interactions
    from a specified dataset file. The function supports custom styling, axis
    labels, legends, and error bars as defined by the provided settings.

    The function handles both normal and merged datasets, processes input data
    for required sequences of residues and interactions, applies chain and
    segment renaming if specified, and generates the x-axis labels accordingly.
    Customizations such as figure size, colors, fonts, and gridlines are applied
    based on the preferences provided in the settings object. The final plot
    is saved to the specified output directory.

    :param pta_file: The protein-ligand interaction file used as the data source.
    :type pta_file: CustomFileType
    :param group_name: Name of the receptor group used for dataset extraction.
    :type group_name: str
    :param mode_name: Identifier for the ligand mode being analyzed.
    :type mode_name: str
    :param settings: The configuration object that defines figure appearance
        and data processing options.
    :type settings: PLIStackedSettings1
    :param out_dir: Path to the directory where the plot image will be saved.
    :type out_dir: pathlib.Path
    :param attach_to_name: Designates a prefix to append to the final saved image
        file name. Defaults to "mode".
    :type attach_to_name: str
    :param is_merged: Flag indicating whether the dataset is merged.
    :type is_merged: bool
    :return: None
    """
    logger.debug(
        "plot_protein_ligand_interactions_stacked_column_1_from_file: "
        "group='%s' mode='%s' is_merged=%s out_dir='%s'",
        group_name, mode_name, is_merged, out_dir,
    )

    # DATASET PATH
    if is_merged:
        dset_path = f"{group_name}/modes_merged/{mode_name}/table"
        fig_basename = f"merged_{settings.fig_basename}"
    else:
        dset_path = f"{group_name}/modes/{mode_name}/table"
        fig_basename = settings.fig_basename

    logger.debug("Dataset path: '%s'", dset_path)
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError(f"Dataset empty: {dset_path}")
    logger.debug("Dataset rows: %d", int(getattr(dset, "size", 0)))

    # OUTPUT PATH
    ext = f".{settings.fig_format}"
    out_path = Path(out_dir) / f"{attach_to_name}-{fig_basename}"
    out_path = out_path.with_suffix(ext)

    # BUILD DATA
    if is_merged:
        ordered_residues, interactions, values, errors = \
            build_pli_merged_stacked_data(
                pta_file,
                group_name=group_name,
                mode_name=mode_name,
                threshold=settings.threshold,
                aa3_to_aa1=settings.aa3_to_aa1,
                renumber=settings.renumber,
                renumber_int=settings.renumber_int,
                debug=False,
            )
        error_bars = settings.error_bars
    else:
        ordered_residues, interactions, values = \
            build_pli_normal_data(
                pta_file,
                group_name=group_name,
                mode_name=mode_name,
                threshold=settings.threshold,
                aa3_to_aa1=settings.aa3_to_aa1,
                renumber=settings.renumber,
                renumber_int=settings.renumber_int,
                debug=False,
            )
        errors = None
        error_bars = False

    if values.size == 0:
        logger.warning(f"No data available after processing. Skipping {mode_name}.")
        return

    # ALTER CHAINS / SEGMENTS
    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    segment_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    x_labels = []

    for resname, resid, chainid, segid in ordered_residues:

        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= int(resid) <= hi:
                    chainid = new_chain
                    break

        if settings.alter_segments:
            for new_seg, (lo, hi) in segment_map.items():
                if lo <= int(resid) <= hi:
                    segid = new_seg
                    break

        label = settings.x_axis_representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))

        x_labels.append(label)


    # FIGURE SETUP
    plt.rcParams["font.family"] = settings.font_family
    plt.rcParams["axes.labelsize"] = settings.font_size_label
    plt.rcParams["axes.labelweight"] = settings.font_weight_label

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    x = np.arange(len(x_labels))
    cumulative = np.zeros(len(x_labels), dtype=float)

    interaction_colors = {
        "hydrophobic": settings.color_hydrophobic,
        "hydrogen_bonds": settings.color_hydrogen_bonds,
        "pi_cation": settings.color_pi_cation,
        "pi_stacking": settings.color_pi_stacking,
        "water_bridge": settings.color_water_bridge_1,
        "ionic": settings.color_ionic,
        "halogen": settings.color_halogen,
        "metal_contact": settings.color_metal_contact,
    }

    colors = [interaction_colors.get(k, "gray") for k in interactions]

    # STACKED BARS
    for j, (ikey, col) in enumerate(zip(interactions, colors)):

        ax.bar(
            x,
            values[:, j],
            bottom=cumulative,
            width=settings.bar_width,
            color=col,
            alpha=settings.bar_alpha,
            edgecolor=settings.bar_edge_color,
            linewidth=settings.bar_edge_width,
            label=ikey.replace("_", " ").title(),
        )

        if error_bars and errors is not None:
            ax.errorbar(
                x,
                cumulative + values[:, j],
                yerr=errors[:, j],
                fmt="none",
                ecolor=settings.error_bars_color,
                elinewidth=settings.error_bars_line_width,
                alpha=settings.error_bars_alpha,
                capsize=settings.error_bars_capsize,
                linestyle=settings.error_bars_line_style,
            )

        cumulative += values[:, j]

    # AXES
    if not settings.disable_x_axis:
        ax.set_xticks(x)
        ax.set_xticklabels(
            x_labels,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_ticks,
            fontweight=settings.font_weight_ticks,
        )
    else:
        ax.get_xaxis().set_visible(False)

    if not settings.disable_y_axis:
        ax.tick_params(axis="y", labelsize=settings.font_size_ticks)

        if not settings.disable_ticks:
            for tick in ax.get_yticklabels():
                tick.set_rotation(settings.y_tick_rotation)
                tick.set_fontweight(settings.font_weight_ticks)
    else:
        ax.get_yaxis().set_visible(False)

    if not settings.disable_x_label:
        ax.set_xlabel(
            settings.x_label,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )

    if not settings.disable_y_label:
        ax.set_ylabel(
            settings.y_label,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )

    if not settings.disable_title:
        ax.set_title(
            settings.fig_title,
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    if settings.disable_ticks:
        ax.tick_params(left=False, bottom=False)

    # Y LIMITS
    if settings.y_limit_min is not None or settings.y_limit_max is not None:
        ax.set_ylim(settings.y_limit_min, settings.y_limit_max)

    # GRID
    if settings.enable_grid:
        ax.grid(
            linestyle=settings.grid_style,
            color=settings.grid_color,
            alpha=settings.grid_alpha,
        )

    # LEGEND
    if not settings.disable_legend:
        leg = ax.legend(
            loc=settings.legend_loc,
            ncol=settings.legend_n_col,
            fontsize=settings.font_size_legend,
            frameon=settings.legend_frame,
            bbox_to_anchor=(0.5, settings.legend_bbox_y),
        )

        if leg:
            leg.get_frame().set_alpha(settings.legend_alpha)
            for text in leg.get_texts():
                text.set_fontweight(settings.font_weight_legend)
            leg.set_in_layout(False)

    # LAYOUT
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.92,
        bottom=0.25
    )

    if settings.tight_layout:
        fig.tight_layout()

    # SAVE

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        format=settings.fig_format,
        transparent=settings.fig_transparent,
    )

    plt.close(fig)

    logger.info(f"Saved stacked column type 1 to {out_path}")


def plot_protein_ligand_interactions_stacked_column_2_from_file(pta_file, *,
                                                                group_name: str,
                                                                settings: PLIStackedSettings2,
                                                                out_dir: Path,
                                                                attach_to_name: str = "",
                                                                is_merged: bool = False) -> None:
    """
    Plots a stacked column chart representing protein-ligand interactions for a given PTA file.
    The chart displays the frequency of interactions categorized into backbone and sidechain
    regions, computed across all frames extracted from the specified group.

    :param pta_file: The input file containing the protein-ligand interaction data as PTA format.
    :type pta_file: PTAFile
    :param group_name: The dataset group from which interaction frames are extracted.
    :type group_name: str
    :param settings: Configuration settings for plotting the stacked column chart.
    :type settings: PLIStackedSettings2
    :param out_dir: Output directory where the resulting figure will be saved.
    :type out_dir: Path
    :param attach_to_name: An optional string to attach to the output file name.
    :type attach_to_name: str, optional
    :param is_merged: Flag indicating whether the PTA file is merged.
                      Only disjoint PTA files are supported in this method.
    :type is_merged: bool, optional
    :raises RuntimeError: If provided with merged PTA files, or fails to find interaction data
                          or frames for the specified group, or if all residues are filtered
                          by the threshold.
    :return: None. Outputs the generated figure to the specified directory.
    """
    logger.debug(
        "plot_protein_ligand_interactions_stacked_column_2_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError(
            "Stacked Column Type 2 is not supported for merged PTA files."
        )

    # Frame discovery
    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    frame_indices = list(pta_file._iter_frames(group_name))
    if not frame_indices:
        raise RuntimeError("No frames found.")
    logger.debug("Discovered %d interaction frame(s)", len(frame_indices))

    n_frames = float(len(frame_indices))
    residue_data = {}

    # Helper
    def _extract_protein_residue_from_row(rec):
        """
        Extracts protein residue information from a given record (row).

        This function processes the input record based on the label provided in the record
        and extracts specific properties such as atom type, residue name, residue ID,
        chain ID, and segment ID. The logic for extracting these properties varies
        depending on the label ("PI-STACKING", "PI-CATION", or other labels) in the
        record.

        :param rec: List containing a record to extract protein residue information from.
            The structure and indices of this list depend on the specific label
            provided in the record.
        :type rec: list
        :return: A tuple containing the atom type, residue name, residue ID,
            chain ID, and segment ID.
        :rtype: tuple
        """

        label = rec[0]

        if label == "PI-STACKING":
            atom_type = rec[5].split(",")[0].upper()
            resname = rec[6].split(",")[0]
            resid = int(rec[7].split(",")[0])
            chainid = rec[8].split(",")[0]
            segid = rec[9].split(",")[0]

        elif label == "PI-CATION":
            atom_type = str(rec[6]).upper()
            resname = str(rec[7]).split(",")[0]
            resid = int(str(rec[8]).split(",")[0])
            chainid = str(rec[9]).split(",")[0]
            segid = str(rec[10]).split(",")[0]

        else:
            atom_type = str(rec[5]).upper()
            resname = str(rec[6])
            resid = int(rec[7])
            chainid = str(rec[8])
            segid = str(rec[9])

        return atom_type, resname, resid, chainid, segid

    # Stream frames
    for frame in frame_indices:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        seen_in_frame = set()

        for row in dset:
            rec = json.loads(row)

            atom_type, resname, resid, chainid, segid = \
                _extract_protein_residue_from_row(rec)

            bbsc = "Backbone" if atom_type == "BB" else "Sidechain"

            residue_key = (resname, resid, chainid, segid)

            unique_frame_key = (residue_key, bbsc)
            if unique_frame_key in seen_in_frame:
                continue

            seen_in_frame.add(unique_frame_key)

            if residue_key not in residue_data:
                residue_data[residue_key] = {
                    "Backbone": set(),
                    "Sidechain": set(),
                }

            residue_data[residue_key][bbsc].add(frame)

    if not residue_data:
        raise RuntimeError("No interaction data collected.")

    # Frequency computation
    processed = []

    for (resname, resid, chainid, segid), r in residue_data.items():

        bb_freq = len(r["Backbone"]) / n_frames
        sc_freq = len(r["Sidechain"]) / n_frames

        if (bb_freq + sc_freq) < settings.threshold:
            continue

        processed.append(
            (resname, resid, chainid, segid, bb_freq, sc_freq)
        )

    if not processed:
        raise RuntimeError("All residues filtered by threshold.")

    # Sorting
    processed.sort(key=lambda x: x[1])  # by resid

    # Label construction
    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    segment_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    labels = []
    backbone_vals = []
    sidechain_vals = []

    for resname, resid, chainid, segid, bb_freq, sc_freq in processed:

        # AA3 → AA1
        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        # Renumber
        if settings.renumber and settings.renumber_int is not None:
            resid = resid + settings.renumber_int

        # Alter chains
        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= int(resid) <= hi:
                    chainid = new_chain
                    break

        # Alter segments
        if settings.alter_segments:
            for new_seg, (lo, hi) in segment_map.items():
                if lo <= int(resid) <= hi:
                    segid = new_seg
                    break

        # Build label using representation
        label = settings.x_axis_representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))

        labels.append(label)
        backbone_vals.append(bb_freq)
        sidechain_vals.append(sc_freq)

    backbone_vals = np.array(backbone_vals, dtype=float)
    sidechain_vals = np.array(sidechain_vals, dtype=float)
    x = np.arange(len(labels))

    # Output path
    ext = "." + settings.fig_format
    out_name = f"{attach_to_name}{settings.fig_basename}"
    out_path = (Path(out_dir) / out_name).with_suffix(ext)

    # Figure
    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    # Bars
    ax.bar(
        x,
        backbone_vals,
        label="Backbone",
        color=settings.color_backbone,
        edgecolor=settings.bar_edge_color,
        width=settings.bar_width,
        alpha=settings.bar_alpha,
        linewidth=settings.bar_edge_width,
    )

    ax.bar(
        x,
        sidechain_vals,
        bottom=backbone_vals,
        label="Sidechain",
        color=settings.color_side_chain,
        edgecolor=settings.bar_edge_color,
        width=settings.bar_width,
        alpha=settings.bar_alpha,
        linewidth=settings.bar_edge_width,
    )

    # Title
    if not settings.disable_title:
        ax.set_title(
            settings.fig_title,
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    # Axis labels (NOW using font_size_x/y)
    if not settings.disable_x_axis and not settings.disable_x_label:
        ax.set_xlabel(
            settings.x_label,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )

    if not settings.disable_y_axis and not settings.disable_y_label:
        ax.set_ylabel(
            settings.y_label,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )

    # Ticks
    if settings.disable_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.set_xticks(x)
        ax.set_xticklabels(
            labels,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_ticks,
            fontweight=settings.font_weight_ticks,
        )

        ax.tick_params(axis="y", labelsize=settings.font_size_ticks)

        for tick in ax.get_yticklabels():
            tick.set_rotation(settings.y_tick_rotation)
            tick.set_fontweight(settings.font_weight_ticks)

    # Y limits
    if settings.y_limit_min is not None or settings.y_limit_max is not None:
        ax.set_ylim(settings.y_limit_min, settings.y_limit_max)

    # Grid
    if settings.enable_grid:
        ax.grid(
            True,
            axis="y",
            linestyle=settings.grid_style,
            color=settings.grid_color,
            alpha=settings.grid_alpha,
        )

    # Legend
    if not settings.disable_legend:

        handles, labels_ = ax.get_legend_handles_labels()

        leg = ax.legend(
            handles,
            labels_,
            loc=settings.legend_loc,
            bbox_to_anchor=(0.5, settings.legend_bbox_y)
            if settings.legend_loc in ("lower center", "upper center")
            else None,
            ncol=min(settings.legend_n_col, len(labels_)),
            frameon=settings.legend_frame,
            framealpha=settings.legend_alpha,
            fontsize=settings.font_size_legend,
        )

        for txt in leg.get_texts():
            txt.set_fontweight(settings.font_weight_legend)

    if settings.tight_layout:
        plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
    )

    plt.close(fig)

    logger.info(f"Saved stacked column type 2 to {out_path}")

    return None


def plot_protein_ligand_interactions_heatmap_1_from_file(pta_file, *,
                                                         group_name: str,
                                                         settings: PLIHeatmapSettings1,
                                                         out_dir: Path,
                                                         attach_to_name: str = "",
                                                         is_merged: bool = False) -> None:
    """
    Plots a protein-ligand interaction heatmap (type 1) from the given PTA file. The heatmap visualizes
    interactions between residues and ligands over frames, filtered and transformed based on user-defined
    settings. Results are saved to the specified output directory.

    :param pta_file: PTA file object containing interaction data
    :param group_name: Group name within the PTA file to extract data from
    :param settings: Settings object defining visualization and threshold parameters
    :param out_dir: Directory to save the resulting heatmap
    :param attach_to_name: Optional string to attach to the output file name
    :param is_merged: Boolean indicating if the PTA file represents merged data
    :return: None
    """
    logger.debug(
        "plot_protein_ligand_interactions_heatmap_1_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError(
            "Heatmap Type 1 is not supported for merged PTA files."
        )

    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    # Frame discovery
    frame_indices = list(pta_file._iter_frames(group_name))
    if not frame_indices:
        raise RuntimeError("No frames found.")
    logger.debug("Discovered %d interaction frame(s)", len(frame_indices))

    residue_frame_counts = {}
    n_frames = len(frame_indices)

    # Extract data (FRAME BASED)
    for frame in frame_indices:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        seen_in_frame = set()

        for row in dset:
            rec = json.loads(row)

            label = rec[0]

            if label == "PI-STACKING":
                resname = rec[6].split(",")[0]
                resid = int(rec[7].split(",")[0])
                chainid = rec[8].split(",")[0]
                segid = rec[9].split(",")[0]
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue
            elif label == "PI-CATION":
                resname = str(rec[7]).split(",")[0]
                resid = int(str(rec[8]).split(",")[0])
                chainid = str(rec[9]).split(",")[0]
                segid = str(rec[10]).split(",")[0]
                atom_type = rec[6]
                # Ligand atom
                atom2_name = str(rec[12]).strip()
                if not atom2_name:
                    continue
            else:
                resname = str(rec[6])
                resid = int(rec[7])
                chainid = str(rec[8])
                segid = str(rec[9])
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue

            residue_key = (resname, resid, chainid, segid)

            if residue_key not in residue_frame_counts:
                residue_frame_counts[residue_key] = set()

            if residue_key in seen_in_frame:
                continue

            seen_in_frame.add(residue_key)
            residue_frame_counts[residue_key].add(frame)

    if not residue_frame_counts:
        raise RuntimeError("No interaction data collected.")

    # Apply threshold
    threshold = settings.threshold
    filtered = []

    for key, frames in residue_frame_counts.items():
        freq = len(frames) / float(n_frames)
        if freq >= threshold:
            filtered.append((key, frames))

    if not filtered:
        raise RuntimeError("All residues filtered by threshold.")

    # Sort residues by frequency
    filtered.sort(key=lambda x: len(x[1]), reverse=True)

    # Build matrix
    frame_list_sorted = sorted(frame_indices)
    frame_index_map = {f: i for i, f in enumerate(frame_list_sorted)}

    heatmap_matrix = np.zeros(
        (len(filtered), len(frame_list_sorted)),
        dtype=float,
    )

    labels = []

    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    segment_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    for row_idx, ((resname, resid, chainid, segid), frames) in enumerate(filtered):

        # AA3 → AA1
        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        # Renumber
        if settings.renumber and settings.renumber_int is not None:
            resid = resid + settings.renumber_int

        # Alter chains
        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= int(resid) <= hi:
                    chainid = new_chain
                    break

        # Alter segments
        if settings.alter_segments:
            for new_seg, (lo, hi) in segment_map.items():
                if lo <= int(resid) <= hi:
                    segid = new_seg
                    break

        # Build label
        label = settings.y_axis_representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))

        labels.append(label)

        for frame in frames:
            col_idx = frame_index_map[frame]
            heatmap_matrix[row_idx, col_idx] = 1.0

    # Output path
    ext = "." + settings.fig_format
    out_name = f"{attach_to_name}{settings.fig_basename}"
    out_path = (Path(out_dir) / out_name).with_suffix(ext)

    # Figure
    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    im = ax.imshow(
        heatmap_matrix,
        aspect="auto",
        interpolation=settings.interpolation,
        cmap=settings.cmap,
        vmin=settings.vmin,
        vmax=settings.vmax,
    )

    # Title
    if not settings.disable_title:
        ax.set_title(
            settings.fig_title,
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    # X ticks (adaptive reduction)
    n_frames_total = len(frame_list_sorted)

    if settings.disable_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

    else:
        # Always keep full matrix — reduce only tick labels
        max_ticks = 50  # safe upper bound for readability

        if n_frames_total <= max_ticks:
            tick_indices = np.arange(n_frames_total)
        else:
            tick_indices = np.linspace(
                0,
                n_frames_total - 1,
                max_ticks,
                dtype=int
            )

        ax.set_xticks(tick_indices)
        ax.set_xticklabels(
            [frame_list_sorted[i] for i in tick_indices],
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )

        # Y ticks unchanged
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels(
            labels,
            rotation=settings.y_tick_rotation,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )

    # Axis labels
    if not settings.disable_x_axis and not settings.disable_x_label:
        ax.set_xlabel(
            settings.x_label,
            fontsize=settings.font_size_label,
            fontweight=settings.font_weight_label,
        )

    if not settings.disable_y_axis and not settings.disable_y_label:
        ax.set_ylabel(
            settings.y_label,
            fontsize=settings.font_size_label,
            fontweight=settings.font_weight_label,
        )

    # Colorbar
    if not settings.disable_legend:
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation=settings.cbar_orientation,
            shrink=settings.cbar_shrink,
            pad=settings.cbar_pad,
        )

        cbar.ax.tick_params(labelsize=settings.font_size_cbar)

    # Layout
    if settings.tight_layout:
        plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
    )

    plt.close(fig)

    logger.info(f"Saved heatmap type 1 to {out_path}")

    return None


def plot_protein_ligand_interactions_heatmap_2_from_file(pta_file, *,
                                                         group_name: str,
                                                         settings: PLIHeatmapSettings2,
                                                         out_dir: Path,
                                                         attach_to_name: str = "",
                                                         is_merged: bool = False) -> None:
    """
    Generates and saves a heatmap visualization for protein-ligand interactions from a PTA file.

    This function processes a PTA file to visualize interactions of residues over a series of
    frames for a given group. It supports several normalization modes and provides customizable
    settings to control the appearance of the heatmap figure, such as figure size, color map,
    labels, and gridline options.

    :title param pta_file: An object representing the parsed PTA file containing interaction data.
    :type pta_file: PTAFile
    :param group_name: The name of the group in the PTA file to visualize.
    :param settings: Configuration settings specifying appearance and normalization options
                     for the heatmap. Must be an instance of PLIHeatmapSettings2.
    :param out_dir: The directory path where the heatmap figure should be saved.
    :param attach_to_name: An optional string to prepend to the figure filename.
    :param is_merged: A boolean flag indicating whether the PTA file contains merged frames.
                      This mode is not supported for this heatmap type.
    :return: This method does not return any value; it saves the generated heatmap to the
             specified output directory.

    :raises RuntimeError: If the heatmap type is used in merged PTA files, the specified
                          group is missing in the PTA file, there are no frames or interaction
                          data, or the data is filtered out based on the provided settings.
    """
    logger.debug(
        "plot_protein_ligand_interactions_heatmap_2_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError("Heatmap Type 2 is not supported for merged PTA files.")

    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    # Frame discovery
    frame_indices = list(pta_file._iter_frames(group_name))
    if not frame_indices:
        raise RuntimeError("No frames found.")
    logger.debug("Discovered %d interaction frame(s)", len(frame_indices))

    frame_list_sorted = sorted(frame_indices)
    frame_index_map = {f: i for i, f in enumerate(frame_list_sorted)}

    # Build interaction × frame structure
    interaction_frame_residues = {}

    for frame in frame_list_sorted:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]

        for row in dset:
            rec = json.loads(row)
            interaction = rec[0]

            if interaction == "PI-STACKING":
                resname = rec[6].split(",")[0]
                resid = int(rec[7].split(",")[0])
                chainid = rec[8].split(",")[0]
                segid = rec[9].split(",")[0]
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue
            elif interaction == "PI-CATION":
                resname = str(rec[7]).split(",")[0]
                resid = int(str(rec[8]).split(",")[0])
                chainid = str(rec[9]).split(",")[0]
                segid = str(rec[10]).split(",")[0]
                atom_type = rec[6]
                # Ligand atom
                atom2_name = str(rec[12]).strip()
                if not atom2_name:
                    continue
            else:
                resname = str(rec[6])
                resid = int(rec[7])
                chainid = str(rec[8])
                segid = str(rec[9])
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue

            residue_label = f"{resname}:{resid}"

            interaction_frame_residues.setdefault(interaction, {})
            interaction_frame_residues[interaction].setdefault(frame, set())
            interaction_frame_residues[interaction][frame].add(residue_label)

    if not interaction_frame_residues:
        raise RuntimeError("No interaction data collected.")

    # Convert to matrix
    interactions_sorted = sorted(interaction_frame_residues.keys())

    heatmap_matrix = np.zeros(
        (len(interactions_sorted), len(frame_list_sorted)),
        dtype=float,
    )

    for i, inter in enumerate(interactions_sorted):
        for frame, residues in interaction_frame_residues[inter].items():
            j = frame_index_map[frame]
            heatmap_matrix[i, j] = float(len(residues))

    # Drop empty interactions
    if settings.drop_empty_rows:
        mask = heatmap_matrix.sum(axis=1) > 0
        heatmap_matrix = heatmap_matrix[mask]
        interactions_sorted = [
            inter for inter, keep in zip(interactions_sorted, mask) if keep
        ]

    if heatmap_matrix.size == 0:
        raise RuntimeError("No data left after filtering.")

    # Normalization
    if settings.normalize == "by_frame":
        col_max = heatmap_matrix.max(axis=0)
        col_max[col_max == 0] = 1.0
        heatmap_matrix = heatmap_matrix / col_max
        vmin, vmax = 0.0, 1.0
        cbar_label = "Normalized per frame (0–1)"

    elif settings.normalize == "max1":
        gmax = heatmap_matrix.max() or 1.0
        heatmap_matrix = heatmap_matrix / gmax
        vmin, vmax = 0.0, 1.0
        cbar_label = "Normalized (global max = 1)"

    else:
        vmin, vmax = settings.vmin, settings.vmax
        cbar_label = "Unique residues per frame"

    # Fixed figure size (NO SCALING)
    n_inter, n_frames = heatmap_matrix.shape

    fig_size_w = settings.fig_size_width or 10
    fig_size_h = settings.fig_size_height or 6

    # Output path

    ext = "." + settings.fig_format
    out_name = f"{attach_to_name}{settings.fig_basename}"
    out_path = (Path(out_dir) / out_name).with_suffix(ext)

    # Plot
    fig, ax = plt.subplots(
        figsize=(fig_size_w, fig_size_h),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    im = ax.imshow(
        heatmap_matrix,
        aspect="auto",
        interpolation=settings.interpolation,
        cmap=settings.cmap,
        vmin=vmin,
        vmax=vmax,
    )

    # Title
    if not settings.disable_title:
        ax.set_title(
            settings.fig_title,
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
            pad=20,
        )

    # Adaptive X ticks (NO overcrowding)
    if settings.disable_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        max_ticks = settings.xtick_max

        if n_frames <= max_ticks:
            tick_indices = np.arange(n_frames)
        else:
            tick_indices = np.linspace(
                0,
                n_frames - 1,
                max_ticks,
                dtype=int,
            )

        ax.set_xticks(tick_indices)
        ax.set_xticklabels(
            [frame_list_sorted[i] for i in tick_indices],
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )

        ax.set_yticks(np.arange(n_inter))
        ax.set_yticklabels(
            interactions_sorted,
            fontsize=settings.font_size_y,
            rotation=settings.y_tick_rotation,
            fontweight=settings.font_weight_y,
        )

    # Axis labels
    if not settings.disable_x_axis and not settings.disable_x_label:
        ax.set_xlabel(
            settings.x_label,
            fontsize=settings.font_size_label,
            fontweight=settings.font_weight_label,
        )

    if not settings.disable_y_axis and not settings.disable_y_label:
        ax.set_ylabel(
            settings.y_label,
            fontsize=settings.font_size_label,
            fontweight=settings.font_weight_label,
        )

    # Colorbar
    if not settings.disable_legend:
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation=settings.cbar_orientation,
            shrink=settings.cbar_shrink,
            pad=settings.cbar_pad,
        )
        cbar.ax.tick_params(labelsize=settings.font_size_cbar)
        cbar.set_label(
            cbar_label,
            fontsize=settings.font_size_cbar,
            fontweight=settings.font_weight_label,
        )

    # Gridlines
    if settings.enable_grid:
        ax.set_xticks(np.arange(-0.5, n_frames, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_inter, 1), minor=True)
        ax.grid(
            which="minor",
            linestyle="-",
            linewidth=settings.grid_linewidth,
            color=settings.grid_color,
        )

    # Layout & Save
    if settings.tight_layout:
        fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
    )

    plt.close(fig)

    logger.info(f"Saved heatmap type 2 to {out_path}")

    return None


def plot_protein_ligand_interactions_pie_charts_from_file(pta_file, *,
                                                          group_name: str,
                                                          settings: PLIPieChartsSettings1,
                                                          out_dir: Path,
                                                          attach_to_name: str = "",
                                                          is_merged: bool = False) -> None:
    """
    Plots pie charts, based on protein-ligand interaction (PLI) data, for visualization of
    side-chain (SC) and backbone (BB) interactions. Processes a Protein-Ligand
    Interaction (PTA) file for a specified group and generates pie charts and summaries
    of the interactions using provided plot settings. It creates individual pie charts
    for residues or atom groups, optionally a collage of these charts, and an overall
    summary pie chart for all interactions.

    :param pta_file: A PTA file object containing protein-ligand interaction data.
    :type pta_file: PTAFile
    :param group_name: The name of the group within the PTA file to analyze.
    :type group_name: str
    :param settings: Configuration for generating pie charts, including plot formatting
        preferences such as colors, layout, and labels.
    :type settings: PLIPieChartsSettings1
    :param out_dir: Directory where the resulting plots and charts are stored.
    :type out_dir: Path
    :param attach_to_name: A string to append to the generated filenames for the pie charts
        for context, e.g., molecule name or other details. Default is an empty string.
    :type attach_to_name: str, optional
    :param is_merged: Boolean flag determining if the analysis is on merged PTA files,
        which is not supported for pie charts. Default is False.
    :type is_merged: bool, optional
    :return: This function does not return a value. It saves generated pie charts to the
        specified output directory.
    :rtype: None
    """
    logger.debug(
        "plot_protein_ligand_interactions_pie_charts_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError("Pie charts are not supported for merged PTA files.")

    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    # Basic settings
    dpi = settings.fig_dpi
    figsize = settings.fig_size
    bg_color = settings.bg_color
    transparent = settings.fig_transparent
    output_format = settings.fig_format

    tight_layout = settings.tight_layout
    top_n = settings.top_n

    collage = settings.collage
    collage_cols = settings.collage_cols
    collage_pad = settings.collage_pad

    make_pdf = settings.make_pdf
    make_overall = settings.make_overall

    disable_title = settings.disable_title
    disable_labels = settings.disable_labels
    disable_autopct = settings.disable_autopct
    disable_overall_title = settings.disable_overall_title

    fontsize_title = settings.font_size_title
    fontsize_pct = settings.font_size_pct
    fontweight_title = settings.font_weight_title

    color_sidechain = settings.color_side_chain
    color_backbone = settings.color_backbone

    # Frame discovery
    frame_indices = list(pta_file._iter_frames(group_name))
    if not frame_indices:
        raise RuntimeError("No frames found.")
    logger.debug("Discovered %d interaction frame(s)", len(frame_indices))

    residue_data = {}

    for frame in frame_indices:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        seen_in_frame = set()

        for row in dset:
            rec = json.loads(row)
            interaction = rec[0]

            if interaction == "PI-STACKING":
                resname = rec[6].split(",")[0]
                resid = int(rec[7].split(",")[0])
                chainid = rec[8].split(",")[0]
                segid = rec[9].split(",")[0]
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue
            elif interaction == "PI-CATION":
                resname = str(rec[7]).split(",")[0]
                resid = int(str(rec[8]).split(",")[0])
                chainid = str(rec[9]).split(",")[0]
                segid = str(rec[10]).split(",")[0]
                atom_type = rec[6]
                # Ligand atom
                atom2_name = str(rec[12]).strip()
                if not atom2_name:
                    continue
            else:
                resname = str(rec[6])
                resid = int(rec[7])
                chainid = str(rec[8])
                segid = str(rec[9])
                atom_type = rec[5]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue

            atom_type = str(atom_type).upper()
            if atom_type not in ("SC", "BB"):
                continue

            label = settings.x_axis_representation
            label = label.replace("resname", str(resname))
            label = label.replace("resid", str(resid))
            label = label.replace("chainid", str(chainid))
            label = label.replace("segid", str(segid))

            residue_key = (label, atom_type)
            if residue_key in seen_in_frame:
                continue
            seen_in_frame.add(residue_key)

            residue_data.setdefault(label, {"SC": set(), "BB": set()})
            residue_data[label][atom_type].add(frame)

    if not residue_data:
        raise RuntimeError("No SC/BB data available for pie charts.")

    # Build summary
    summary = []
    for residue, data in residue_data.items():
        n_SC = len(data["SC"])
        n_BB = len(data["BB"])
        total = n_SC + n_BB
        if total == 0:
            continue

        pct_SC = (n_SC / total) * 100.0
        pct_BB = (n_BB / total) * 100.0

        summary.append((residue, n_SC, n_BB, total, pct_SC, pct_BB))

    summary.sort(key=lambda x: x[3], reverse=True)

    if top_n > 0:
        summary = summary[:top_n]

    if not summary:
        raise RuntimeError("No data to plot pie charts.")

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1 Individual Pies
    for i, (residue, n_SC, n_BB, total, pct_SC, pct_BB) in enumerate(summary, start=1):

        fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi, facecolor=bg_color)

        labels = None if disable_labels else ["Sidechain", "Backbone"]
        autopct = None if disable_autopct else (lambda p: f"{p:.1f}%")

        ax.pie(
            [pct_SC, pct_BB],
            labels=labels,
            colors=[color_sidechain, color_backbone],
            autopct=autopct,
            startangle=90,
            counterclock=False,
            textprops={"fontsize": fontsize_pct},
        )
        ax.axis("equal")

        if not disable_title:
            ax.set_title(
                f"{residue} (n={total})",
                fontsize=fontsize_title,
                fontweight=fontweight_title,
                pad=14,
            )

        fpath = out_dir / f"{i:03d}_{residue}.{output_format}"
        if tight_layout:
            plt.tight_layout()

        fig.savefig(fpath, dpi=dpi, transparent=transparent, bbox_inches="tight")
        plt.close(fig)

    # 2️ Collage
    if collage and len(summary) > 1:

        n = len(summary)
        cols = max(1, collage_cols)
        rows = math.ceil(n / cols)

        fig, axes = plt.subplots(
            rows,
            cols,
            figsize=(cols * figsize, rows * figsize),
            dpi=dpi,
            facecolor=bg_color,
        )

        axes = np.array(axes).reshape(rows, cols)
        it = iter(summary)

        for r in range(rows):
            for c in range(cols):
                ax = axes[r, c]
                try:
                    residue, n_SC, n_BB, total, pct_SC, pct_BB = next(it)
                except StopIteration:
                    ax.axis("off")
                    continue

                ax.pie(
                    [pct_SC, pct_BB],
                    labels=None if disable_labels else ["SC", "BB"],
                    colors=[color_sidechain, color_backbone],
                    autopct=None if disable_autopct else (lambda p: f"{p:.0f}%"),
                    startangle=90,
                    counterclock=False,
                )
                ax.axis("equal")

                if not disable_title:
                    ax.set_title(
                        f"{residue} (n={total})",
                        fontsize=max(10, fontsize_title - 4),
                        pad=8,
                    )

        plt.tight_layout(pad=collage_pad)

        fig.savefig(
            out_dir / f"pli_pie_collage.{output_format}",
            dpi=dpi,
            transparent=transparent,
            bbox_inches="tight",
        )
        plt.close(fig)

    # 3️ Overall Pie
    if make_overall:

        total_SC = sum(x[1] for x in summary)
        total_BB = sum(x[2] for x in summary)
        total_sum = total_SC + total_BB

        sc_pct = (total_SC / total_sum) * 100 if total_sum else 0
        bb_pct = (total_BB / total_sum) * 100 if total_sum else 0

        fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi, facecolor=bg_color)

        autopct = None if disable_autopct else (lambda p: f"{p:.1f}%")

        ax.pie(
            [sc_pct, bb_pct],
            labels=None if disable_labels else ["Sidechain", "Backbone"],
            colors=[color_sidechain, color_backbone],
            autopct=autopct,
            startangle=90,
            counterclock=False,
        )
        ax.axis("equal")

        if not disable_overall_title:
            ax.set_title(
                f"Overall Composition (n={int(total_sum)})",
                fontsize=fontsize_title,
                fontweight=fontweight_title,
                pad=14,
            )

        fig.savefig(
            out_dir / f"pli_pie_overall.{output_format}",
            dpi=dpi,
            transparent=transparent,
            bbox_inches="tight",
        )
        plt.close(fig)

    # 4️ PDF Export
    if make_pdf:
        from matplotlib.backends.backend_pdf import PdfPages

        pdf_path = out_dir / "pli_pie_charts.pdf"
        with PdfPages(pdf_path) as pdf:

            for residue, n_SC, n_BB, total, pct_SC, pct_BB in summary:

                fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi, facecolor=bg_color)

                ax.pie(
                    [pct_SC, pct_BB],
                    labels=None if disable_labels else ["Sidechain", "Backbone"],
                    colors=[color_sidechain, color_backbone],
                    autopct=None if disable_autopct else (lambda p: f"{p:.1f}%"),
                    startangle=90,
                    counterclock=False,
                )
                ax.axis("equal")

                if not disable_title:
                    ax.set_title(
                        f"{residue} (n={total})",
                        fontsize=fontsize_title,
                        fontweight=fontweight_title,
                        pad=14,
                    )

                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    logger.info(f"Saved unified SC/BB pie charts to {out_dir}")



def plot_protein_ligand_interactions_ligand_monitor_from_file(pta_file, *,
                                                              group_name: str,
                                                              settings: PLILigandMonitorSettings,
                                                              out_dir: Path,
                                                              attach_to_name: str = "",
                                                              is_merged: bool = False) -> None:
    logger.debug(
        "plot_protein_ligand_interactions_ligand_monitor_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError("Ligand monitor is not supported for merged PTA files.")

    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    # Frame discovery
    frame_indices = list(pta_file._iter_frames(group_name))
    if not frame_indices:
        raise RuntimeError("No frames found.")

    total_frames = float(len(frame_indices))

    # Data structure:
    # residue_key -> ligand_atom -> set(frames)
    residue_atom_frames = {}

    for frame in frame_indices:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        seen_in_frame = set()

        for row in dset:
            rec = json.loads(row)

            # Extract protein residue
            label = rec[0]

            if label == "PI-STACKING":
                resname = rec[6].split(",")[0]
                resid = int(rec[7].split(",")[0])
                chainid = rec[8].split(",")[0]
                segid = rec[9].split(",")[0]
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue
            elif label == "PI-CATION":
                resname = str(rec[7]).split(",")[0]
                resid = int(str(rec[8]).split(",")[0])
                chainid = str(rec[9]).split(",")[0]
                segid = str(rec[10]).split(",")[0]
                # Ligand atom
                atom2_name = str(rec[12]).strip()
                if not atom2_name:
                    continue
            else:
                resname = str(rec[6])
                resid = int(rec[7])
                chainid = str(rec[8])
                segid = str(rec[9])
                # Ligand atom
                atom2_name = str(rec[11]).strip()
                if not atom2_name:
                    continue

            residue_key = (resname, resid, chainid, segid)
            unique_key = (residue_key, atom2_name)

            if unique_key in seen_in_frame:
                continue
            seen_in_frame.add(unique_key)

            if residue_key not in residue_atom_frames:
                residue_atom_frames[residue_key] = {}

            if atom2_name not in residue_atom_frames[residue_key]:
                residue_atom_frames[residue_key][atom2_name] = set()

            residue_atom_frames[residue_key][atom2_name].add(frame)

    if not residue_atom_frames:
        logger.warning("No interaction data collected.")
        return

    # Apply residue transformations

    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    segment_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    residue_labels = []
    ligand_atoms = set()

    processed = []

    for (resname, resid, chainid, segid), atom_dict in residue_atom_frames.items():

        # AA3 → AA1
        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        # Renumber
        if settings.renumber and settings.renumber_int is not None:
            resid = resid + settings.renumber_int

        # Alter chains
        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= int(resid) <= hi:
                    chainid = new_chain
                    break

        # Alter segments
        if settings.alter_segments:
            for new_seg, (lo, hi) in segment_map.items():
                if lo <= int(resid) <= hi:
                    segid = new_seg
                    break

        # Build residue label
        label = settings.y_axis_representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))

        residue_labels.append(label)
        processed.append((label, atom_dict))

        ligand_atoms.update(atom_dict.keys())

    ligand_atoms = sorted(ligand_atoms)

    # Build matrix (frequency)
    n_res = len(processed)
    n_lig = len(ligand_atoms)

    matrix = np.zeros((n_res, n_lig), dtype=float)

    for i, (res_label, atom_dict) in enumerate(processed):
        for j, atom in enumerate(ligand_atoms):
            frames = atom_dict.get(atom, set())
            freq = len(frames) / total_frames
            if freq >= settings.threshold:
                matrix[i, j] = freq

    # Drop empty rows/cols
    if settings.drop_empty_rows:
        row_mask = matrix.sum(axis=1) > 0
        matrix = matrix[row_mask]
        residue_labels = [r for r, keep in zip(residue_labels, row_mask) if keep]

    if settings.drop_empty_cols:
        col_mask = matrix.sum(axis=0) > 0
        matrix = matrix[:, col_mask]
        ligand_atoms = [c for c, keep in zip(ligand_atoms, col_mask) if keep]

    if matrix.size == 0:
        logger.warning("Matrix empty after filtering; skipping plot.")
        return

    # Figure
    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    im = ax.imshow(
        matrix,
        aspect="auto",
        interpolation=settings.interpolation,
        cmap=settings.cmap,
        vmin=settings.vmin,
        vmax=settings.vmax,
    )

    # X axis
    if not settings.disable_x_axis:
        ax.set_xticks(np.arange(len(ligand_atoms)))
        ax.set_xticklabels(
            ligand_atoms,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )
        if not settings.disable_x_label:
            ax.set_xlabel(
                settings.x_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )
    else:
        ax.set_xticks([])

    # Y axis
    if not settings.disable_y_axis:
        ax.set_yticks(np.arange(len(residue_labels)))
        ax.set_yticklabels(
            residue_labels,
            rotation=settings.y_tick_rotation,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )
        if not settings.disable_y_label:
            ax.set_ylabel(
                settings.y_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )
    else:
        ax.set_yticks([])

    # Title
    if not settings.disable_title:
        ax.set_title(
            f"{settings.fig_title} (≥ {settings.threshold:.2f})",
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    # Colorbar
    if not settings.disable_colorbar:
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation=settings.cbar_orientation,
            shrink=settings.cbar_shrink,
            pad=settings.cbar_pad,
        )
        cbar.ax.tick_params(labelsize=settings.font_size_cbar)

    # Gridlines
    if settings.enable_grid:
        ax.set_xticks(np.arange(-0.5, len(ligand_atoms), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(residue_labels), 1), minor=True)
        ax.grid(which="minor",
                color=settings.grid_color,
                linewidth=settings.grid_linewidth)

    if settings.tight_layout:
        fig.tight_layout()

    # Save
    ext = "." + settings.fig_format
    out_name = f"{attach_to_name}{settings.fig_basename}"
    out_path = (Path(out_dir) / out_name).with_suffix(ext)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
    )

    plt.close(fig)

    logger.info(f"Saved ligand monitor heatmap to {out_path}")


def parse_frame_interaction_record(rec: List, *, debug: bool = False,) -> Dict[str, object]:
    """
    Parse a single frame-level interaction record from PTA file.

    Returns:
        {
            "interaction": str,
            "protein": {
                "resname": str,
                "resid": int,
                "chainid": str,
                "segid": str,
            },
            "ligand_atoms": List[str]
        }
    """

    interaction = str(rec[0]).strip().upper()

    # -------------------------------------------------
    # Default containers
    # -------------------------------------------------
    resname = None
    resid = None
    chainid = None
    segid = None
    ligand_atoms: List[str] = []

    # -------------------------------------------------
    # PI-STACKING
    # -------------------------------------------------
    if interaction == "PI-STACKING":

        # Protein residue info
        resname = rec[6].split(",")[0]
        resid = int(rec[7].split(",")[0])
        chainid = rec[8].split(",")[0]
        segid = rec[9].split(",")[0]

        # Ligand atom names are at index 12
        if len(rec) > 12 and rec[12]:
            ligand_atoms = [
                a.strip()
                for a in str(rec[12]).split(",")
                if a.strip()
            ]

    # -------------------------------------------------
    # PI-CATION
    # -------------------------------------------------
    elif interaction == "PI-CATION":

        resname = str(rec[7])
        resid = int(rec[8])
        chainid = str(rec[9])
        segid = str(rec[10])

        # Ligand atom names at index 12
        if len(rec) > 12 and rec[12]:
            ligand_atoms = [
                a.strip()
                for a in str(rec[12]).split(",")
                if a.strip()
            ]

    # -------------------------------------------------
    # NORMAL INTERACTIONS
    # -------------------------------------------------
    else:

        resname = str(rec[6])
        resid = int(rec[7])
        chainid = str(rec[8])
        segid = str(rec[9])

        # Single ligand atom at index 11
        if len(rec) > 11 and rec[11]:
            ligand_atoms = [str(rec[11]).strip()]

    # -------------------------------------------------
    # Debug output
    # -------------------------------------------------
    if debug:
        print("-------------------------------------------------")
        print("Interaction:", interaction)
        print("Protein:", resname, resid, chainid, segid)
        print("Ligand atoms:", ligand_atoms)
        print("-------------------------------------------------")

    return {
        "interaction": interaction,
        "protein": {
            "resname": str(resname),
            "resid": int(resid),
            "chainid": str(chainid),
            "segid": str(segid),
        },
        "ligand_atoms": ligand_atoms,
    }



def build_pli_merged_stacked_data(pta_file, *, group_name: str, mode_name: str, threshold: float = 0.0,
                                  aa3_to_aa1: bool = True, renumber: bool = False, renumber_int: int | None = None,
                                  debug: bool = True
                                  ) -> Tuple[List[Tuple[str, int, str, str]], List[str], np.ndarray, np.ndarray]:
    """
    Builds a merged and stacked data representation used in protein-ligand interaction analysis.

    This function processes a dataset of protein-ligand interactions, normalizes interaction
    labels, filters records based on the provided threshold, optionally applies residue renumbering,
    and generates matrices representing interaction values and errors. The resulting data includes an
    ordered list of residues, a list of interaction types, and corresponding values and errors.

    :param pta_file: File object representing the dataset containing protein-ligand interaction data.
    :param group_name: Name of the group in the dataset structure where interaction data is stored.
    :param mode_name: Specific mode name within the group to locate the interaction data.
    :param threshold: Minimum mean frequency value to include an interaction in the results.
    :param aa3_to_aa1: Flag to determine whether to convert residues from 3-letter to 1-letter codes.
    :param renumber: Flag indicating whether residue renumbering should be applied.
    :param renumber_int: Starting integer for residue renumbering, required if renumber is True.
    :param debug: Boolean flag for enabling or disabling debug logging.
    :return: A tuple containing the processed data:
        - List of ordered residues represented as tuples (resname, resid, chainid, segid).
        - List of normalized interaction types.
        - Numpy array of interaction values.
        - Numpy array of interaction errors.
    """

    # normalize labels so plotting can color them
    def norm_interaction(s: str) -> str:
        """
        Normalizes interaction type names by transforming them into a standard
        format. This function takes an input string, removes any leading or
        trailing whitespace, converts it to uppercase, and maps it against
        a predefined dictionary of interaction types. If the input string
        matches a key in the dictionary, the corresponding normalized string
        is returned. Otherwise, it transforms the input to lowercase, replaces
        hyphens with underscores, and returns it.

        :param s: The interaction type as a string. Must be non-empty and may
                  contain leading/trailing spaces or hyphens.
        :type s: str
        :return: A normalized representation of the interaction type. Returns
                 the corresponding value from the mapping dictionary or a
                 transformed version of the input in lowercase with underscores
                 replacing hyphens.
        :rtype: str
        """
        t = str(s).strip().upper()
        return {
            "HYDROPHOBIC": "hydrophobic",
            "HYDROGEN-BOND": "hydrogen_bonds",
            "PI-CATION": "pi_cation",
            "PI-STACKING": "pi_stacking",
            "WATER-BRIDGE": "water_bridge",
            "WATER-BRIDGE-1": "water_bridge",
            "IONIC": "ionic",
            "HALOGEN-BOND": "halogen",
            "METAL-CONTACT": "metal_contact",
        }.get(t, t.lower().replace("-", "_"))

    logger.debug(
        "build_pli_merged_stacked_data: group='%s' mode='%s' threshold=%s "
        "aa3_to_aa1=%s renumber=%s",
        group_name, mode_name, threshold, aa3_to_aa1, renumber,
    )

    dset_path = f"{group_name}/modes_merged/{mode_name}/table"
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError("Merged dataset is empty.")

    logger.debug("Merged dataset path='%s' rows=%d", dset_path, int(dset.size))
    if debug:
        logger.trace("---- DATASET INFO ----")
        logger.trace("Dataset path: {}", dset_path)
        logger.trace("Rows: {}", dset.size)

    resid_map = {}
    if renumber:
        if renumber_int is None:
            raise RuntimeError("renumber=True but renumber_int=None")

        current = renumber_int
        for row in dset:
            rec = json.loads(row)
            r_protein, _, _ = rec["key"]
            resid = r_protein[1]
            if resid not in resid_map:
                resid_map[resid] = current
                current += 1

        if debug:
            logger.trace("---- RENUMBER MAP ----")
            logger.trace("{}", resid_map)

    data = {}
    residue_order = {}

    for i, row in enumerate(dset):
        rec = json.loads(row)

        r_protein, _, interaction_raw = rec["key"]
        interaction = norm_interaction(interaction_raw)

        mean = float(rec["mean_frequency"])
        std = float(rec["std_frequency"])

        if mean < threshold:
            continue

        resname, resid, chainid, segid = r_protein

        if renumber:
            resid = resid_map[resid]

        if aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        residue_key = (str(resname), int(resid), str(chainid).strip(), str(segid).strip())

        if debug and i < 5:
            logger.trace("Row {}", i)
            logger.trace("Raw key: {}", rec["key"])
            logger.trace("Protein selected: {}", residue_key)
            logger.trace("Interaction normalized: {} -> {}", interaction_raw, interaction)

        if residue_key not in data:
            data[residue_key] = {}
            residue_order[residue_key] = int(resid)

        if interaction in data[residue_key]:
            prev_m, prev_s = data[residue_key][interaction]
            data[residue_key][interaction] = (prev_m + mean, prev_s + std)
        else:
            data[residue_key][interaction] = (mean, std)

    if not data:
        raise RuntimeError("No data collected.")

    ordered_residues = [k for k, _ in sorted(residue_order.items(), key=lambda x: x[1])]

    interactions = []
    for v in data.values():
        for ikey in v.keys():
            if ikey not in interactions:
                interactions.append(ikey)

    # Sort into canonical bottom->top order; unknown labels go to the end (stable).
    interactions.sort(key=lambda k: (_PLI_STACK_RANK.get(k, len(_PLI_STACK_ORDER)), k))

    values = np.zeros((len(ordered_residues), len(interactions)))
    errors = np.zeros_like(values)

    for i, rkey in enumerate(ordered_residues):
        for j, ikey in enumerate(interactions):
            if ikey in data[rkey]:
                values[i, j], errors[i, j] = data[rkey][ikey]

    if debug:
        logger.trace("---- INTERACTIONS (normalized) ----")
        logger.trace("{}", interactions)

    logger.debug(
        "build_pli_merged_stacked_data → %d residue(s) × %d interaction(s), "
        "matrix shape=%s",
        len(ordered_residues), len(interactions), values.shape,
    )
    return ordered_residues, interactions, values, errors


def build_pli_normal_data(pta_file, *, group_name: str, mode_name: str, threshold: float = 0.0,
                          aa3_to_aa1: bool = True, renumber: bool = False, renumber_int: int | None = None,
                          debug: bool = False) -> Tuple[List[Tuple[str, int, str, str]], List[str], np.ndarray]:
    logger.debug(
        "build_pli_normal_data: group='%s' mode='%s' threshold=%s "
        "aa3_to_aa1=%s renumber=%s",
        group_name, mode_name, threshold, aa3_to_aa1, renumber,
    )

    def norm_interaction(s: str) -> str:
        t = str(s).strip().upper()
        return {
            "HYDROPHOBIC": "hydrophobic",
            "HYDROGEN-BOND": "hydrogen_bonds",
            "PI-CATION": "pi_cation",
            "PI-STACKING": "pi_stacking",
            "WATER-BRIDGE": "water_bridge",
            "WATER-BRIDGE-1": "water_bridge",
            "IONIC": "ionic",
            "HALOGEN-BOND": "halogen",
            "METAL-CONTACT": "metal_contact",
        }.get(t, t.lower().replace("-", "_"))

    dset_path = f"{group_name}/modes/{mode_name}/table"
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError("Dataset is empty.")

    resid_map = {}
    if renumber:
        if renumber_int is None:
            raise RuntimeError("renumber=True but renumber_int=None")

        current = renumber_int
        for row in dset:
            rec = json.loads(row)
            r_protein, _, _ = rec["key"]
            resid = r_protein[1]
            if resid not in resid_map:
                resid_map[resid] = current
                current += 1

    data = {}
    residue_order = {}

    for i, row in enumerate(dset):
        rec = json.loads(row)

        r_protein, _, interaction_raw = rec["key"]
        interaction = norm_interaction(interaction_raw)

        freq = float(rec["frequency"])

        if freq < threshold:
            continue

        resname, resid, chainid, segid = r_protein

        if renumber:
            resid = resid_map[resid]

        if aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        residue_key = (
            str(resname),
            int(resid),
            str(chainid).strip(),
            str(segid).strip(),
        )

        if residue_key not in data:
            data[residue_key] = {}
            residue_order[residue_key] = int(resid)

        if interaction in data[residue_key]:
            data[residue_key][interaction] += freq
        else:
            data[residue_key][interaction] = freq

    if not data:
        raise RuntimeError("No data collected.")

    ordered_residues = [
        k for k, _ in sorted(residue_order.items(), key=lambda x: x[1])
    ]

    interactions = []
    for v in data.values():
        for ikey in v.keys():
            if ikey not in interactions:
                interactions.append(ikey)

    # Sort into canonical bottom->top order; unknown labels go to the end (stable).
    interactions.sort(key=lambda k: (_PLI_STACK_RANK.get(k, len(_PLI_STACK_ORDER)), k))

    values = np.zeros((len(ordered_residues), len(interactions)))

    for i, rkey in enumerate(ordered_residues):
        for j, ikey in enumerate(interactions):
            if ikey in data[rkey]:
                values[i, j] = data[rkey][ikey]

    logger.debug(
        "build_pli_normal_data → %d residue(s) × %d interaction(s), "
        "matrix shape=%s",
        len(ordered_residues), len(interactions), values.shape,
    )
    return ordered_residues, interactions, values


# PPI
def _collect_ppi_pairs(pta_file, group_name: str, settings):

    logger.debug("_collect_ppi_pairs: group='%s'", group_name)

    frame_indices = sorted(list(pta_file._iter_frames(group_name)))
    total_frames = float(len(frame_indices))
    logger.debug("_collect_ppi_pairs: %d frame(s) to scan", len(frame_indices))

    pair_frames = {}

    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    seg_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    def build_label(resname, resid, chainid, segid):

        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        if settings.renumber and settings.renumber_int is not None:
            resid = resid + settings.renumber_int

        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= resid <= hi:
                    chainid = new_chain
                    break

        if settings.alter_segments:
            for new_seg, (lo, hi) in seg_map.items():
                if lo <= resid <= hi:
                    segid = new_seg
                    break

        label = settings.representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))
        return label

    for frame in frame_indices:

        dset_path = f"{group_name}/frame_{frame}/interactions"
        if dset_path not in pta_file.file:
            continue

        dset = pta_file.file[dset_path]
        seen_this_frame = set()

        for row in dset:
            rec = json.loads(row)
            label = rec[0]

            # -------- π-CATION ----------
            if label == "PI-CATION":
                role = rec[1]
                if role == "ring":
                    r1 = (rec[7].split(",")[0],
                          int(rec[8].split(",")[0]),
                          rec[9].split(",")[0],
                          rec[10].split(",")[0])

                    r2 = (rec[16],
                          int(rec[17]),
                          rec[18],
                          rec[19])
                else:
                    r1 = (rec[7],
                          int(rec[8]),
                          rec[9],
                          rec[10])

                    r2 = (rec[16].split(",")[0],
                          int(rec[17].split(",")[0]),
                          rec[18].split(",")[0],
                          rec[19].split(",")[0])

            # -------- π-STACKING ----------
            elif label == "PI-STACKING":
                r1 = (rec[6].split(",")[0],
                      int(rec[7].split(",")[0]),
                      rec[8].split(",")[0],
                      rec[9].split(",")[0])

                r2 = (rec[15].split(",")[0],
                      int(rec[16].split(",")[0]),
                      rec[17].split(",")[0],
                      rec[18].split(",")[0])

            # -------- STANDARD ----------
            else:
                r1 = (rec[6], int(rec[7]), rec[8], rec[9])
                r2 = (rec[15], int(rec[16]), rec[17], rec[18])

            res1 = build_label(*r1)
            res2 = build_label(*r2)

            pair = "_".join(sorted([res1, res2]))

            if pair not in pair_frames:
                pair_frames[pair] = set()

            if pair not in seen_this_frame:
                pair_frames[pair].add(frame)
                seen_this_frame.add(pair)

    return pair_frames, frame_indices, total_frames


def plot_protein_protein_timeline_pairs_from_file(pta_file, *,
                                                  group_name: str,
                                                  settings: PPITimelinePairsSettings,
                                                  out_dir: Path,
                                                  attach_to_name: str = "",
                                                  is_merged: bool = False) -> None:
    logger.debug(
        "plot_protein_protein_timeline_pairs_from_file: "
        "group='%s' is_merged=%s out_dir='%s'",
        group_name, is_merged, out_dir,
    )

    if is_merged:
        raise RuntimeError("Timeline pairs not supported for merged PTA files.")

    pair_frames, frame_indices, total_frames = _collect_ppi_pairs(
        pta_file, group_name, settings
    )
    logger.debug("Collected %d unique pair(s) across %d frame(s)",
                 len(pair_frames), len(frame_indices))

    if not pair_frames:
        logger.warning("No PPI data collected.")
        return

    pairs_sorted = sorted(
        pair_frames.keys(),
        key=lambda p: len(pair_frames[p]),
        reverse=True,
    )

    if settings.top_n > 0:
        pairs_sorted = pairs_sorted[:settings.top_n]

    n_pairs = len(pairs_sorted)
    n_frames = len(frame_indices)

    matrix = np.zeros((n_pairs, n_frames))

    for i, pair in enumerate(pairs_sorted):
        for j, frame in enumerate(frame_indices):
            if frame in pair_frames[pair]:
                matrix[i, j] = 1.0

    if settings.threshold > 0.0:
        min_frames = settings.threshold * total_frames
        mask = matrix.sum(axis=1) >= min_frames
        matrix = matrix[mask]
        pairs_sorted = [p for p, keep in zip(pairs_sorted, mask) if keep]

    if settings.drop_empty_rows:
        mask = matrix.sum(axis=1) > 0
        matrix = matrix[mask]
        pairs_sorted = [p for p, keep in zip(pairs_sorted, mask) if keep]

    if matrix.size == 0:
        logger.warning("Matrix empty after filtering.")
        return

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    im = ax.imshow(
        matrix,
        aspect="auto",
        interpolation=settings.interpolation,
        cmap=settings.cmap,
        vmin=settings.vmin,
        vmax=settings.vmax,
    )

    if not settings.disable_colorbar:
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation=settings.cbar_orientation,
            shrink=settings.cbar_shrink,
            pad=settings.cbar_pad,
        )
        cbar.ax.tick_params(labelsize=settings.font_size_cbar)

    if not settings.disable_x_axis:
        ax.set_xticks(np.arange(n_frames))
        ax.set_xticklabels(
            frame_indices,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )
        if not settings.disable_x_label:
            ax.set_xlabel(
                settings.x_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )

    if not settings.disable_y_axis:
        ax.set_yticks(np.arange(len(pairs_sorted)))
        ax.set_yticklabels(
            pairs_sorted,
            rotation=settings.y_tick_rotation,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )
        if not settings.disable_y_label:
            ax.set_ylabel(
                settings.y_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )

    if settings.disable_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

    if not settings.disable_title:
        ax.set_title(
            settings.fig_title,
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    if settings.enable_grid:
        ax.set_xticks(np.arange(-0.5, n_frames, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(pairs_sorted), 1), minor=True)
        ax.grid(which="minor",
                color=settings.grid_color,
                linewidth=settings.grid_linewidth)

    if settings.tight_layout:
        fig.tight_layout()

    out_name = f"{attach_to_name}{settings.fig_basename}" if attach_to_name else settings.fig_basename
    out_path = (Path(out_dir) / out_name).with_suffix("." + settings.fig_format)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_path, dpi=settings.fig_dpi,
                transparent=settings.fig_transparent)
    plt.close(fig)

    logger.info(f"Saved PPI timeline pairs to {out_path}")


def plot_protein_protein_heatmap_freq_from_file(pta_file,
                                                *,
                                                group_name: str,
                                                mode_name: str,
                                                settings: PPIHeatmapSettings,
                                                out_dir: Path,
                                                attach_to_name: str = "mode",
                                                is_merged: bool = False) -> None:
    """
    Generates residue–residue interaction heatmap from PTA file.
    Fully settings-driven. Supports merged and normal modes.
    """
    logger.debug(
        "plot_protein_protein_heatmap_freq_from_file: "
        "group='%s' mode='%s' is_merged=%s out_dir='%s'",
        group_name, mode_name, is_merged, out_dir,
    )

    # -------------------------------------------------
    # DATASET PATH (same logic as PLI stacked)
    # -------------------------------------------------
    if is_merged:
        dset_path = f"{group_name}/modes_merged/{mode_name}/table"
        fig_basename = f"merged_{settings.fig_basename}"
    else:
        dset_path = f"{group_name}/modes/{mode_name}/table"
        fig_basename = settings.fig_basename

    logger.debug("Dataset path: '%s'", dset_path)
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError(f"Dataset empty: {dset_path}")
    logger.debug("Dataset rows: %d", int(getattr(dset, "size", 0)))

    # -------------------------------------------------
    # OUTPUT PATH
    # -------------------------------------------------
    ext = f".{settings.fig_format}"
    out_path = Path(out_dir) / f"{attach_to_name}-{fig_basename}"
    out_path = out_path.with_suffix(ext)

    # -------------------------------------------------
    # BUILD DATA
    # -------------------------------------------------
    if is_merged:
        ordered_residues, interactions, values, _ = \
            build_pli_merged_stacked_data(
                pta_file,
                group_name=group_name,
                mode_name=mode_name,
                threshold=settings.threshold,
                aa3_to_aa1=settings.aa3_to_aa1,
                renumber=settings.renumber,
                renumber_int=settings.renumber_int,
                debug=False,
            )
    else:
        ordered_residues, interactions, values = \
            build_pli_normal_data(
                pta_file,
                group_name=group_name,
                mode_name=mode_name,
                threshold=settings.threshold,
                aa3_to_aa1=settings.aa3_to_aa1,
                renumber=settings.renumber,
                renumber_int=settings.renumber_int,
                debug=False,
            )

    if values.size == 0:
        logger.warning(f"No data available after processing. Skipping {mode_name}.")
        return

    # -------------------------------------------------
    # Convert stacked → symmetric contact matrix
    # -------------------------------------------------
    # ordered_residues = list of (resname, resid, chainid, segid)
    residue_labels = []
    for r in ordered_residues:
        label = settings.representation
        label = label.replace("resname", str(r[0]))
        label = label.replace("resid", str(r[1]))
        label = label.replace("chainid", str(r[2]))
        label = label.replace("segid", str(r[3]))
        residue_labels.append(label)

    n = len(residue_labels)
    matrix = np.zeros((n, n))

    # values shape: (n_residues, n_interactions)
    # collapse interactions → total frequency
    totals = values.sum(axis=1)

    for i in range(n):
        matrix[i, i] = totals[i]

    if settings.symmetric:
        matrix = matrix + matrix.T - np.diag(matrix.diagonal())

    # -------------------------------------------------
    # min_total filtering
    # -------------------------------------------------
    if settings.min_total > 0:
        mask = matrix.sum(axis=1) >= settings.min_total
        matrix = matrix[mask][:, mask]
        residue_labels = [r for r, keep in zip(residue_labels, mask) if keep]

    # -------------------------------------------------
    # top_n filtering
    # -------------------------------------------------
    if settings.top_n > 0 and len(residue_labels) > settings.top_n:
        totals = matrix.sum(axis=1)
        idx = np.argsort(totals)[-settings.top_n:]
        matrix = matrix[np.ix_(idx, idx)]
        residue_labels = [residue_labels[i] for i in idx]

    if matrix.size == 0:
        logger.warning(f"No data left after filtering. Skipping {mode_name}.")
        return

    # -------------------------------------------------
    # PLOT
    # -------------------------------------------------
    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    vmax = settings.vmax if settings.vmax is not None else matrix.max()

    im = ax.imshow(
        matrix,
        cmap=settings.cmap,
        vmin=settings.vmin,
        vmax=vmax,
        aspect="auto",
        interpolation=settings.interpolation,
    )

    # ---------------- Title ----------------
    if not settings.disable_title:
        ax.set_title(
            f"{settings.fig_title} ({mode_name})",
            fontsize=settings.font_size_title,
            fontweight=settings.font_weight_title,
        )

    # ---------------- X Axis ----------------
    if not settings.disable_x_axis:
        ax.set_xticks(np.arange(len(residue_labels)))
        ax.set_xticklabels(
            residue_labels,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
            fontweight=settings.font_weight_x,
        )
        if not settings.disable_x_label:
            ax.set_xlabel(
                settings.x_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )
    else:
        ax.set_xticks([])

    # ---------------- Y Axis ----------------
    if not settings.disable_y_axis:
        ax.set_yticks(np.arange(len(residue_labels)))
        ax.set_yticklabels(
            residue_labels,
            rotation=settings.y_tick_rotation,
            fontsize=settings.font_size_y,
            fontweight=settings.font_weight_y,
        )
        if not settings.disable_y_label:
            ax.set_ylabel(
                settings.y_label,
                fontsize=settings.font_size_label,
                fontweight=settings.font_weight_label,
            )
    else:
        ax.set_yticks([])

    # ---------------- Gridlines ----------------
    if settings.enable_grid:
        ax.set_xticks(np.arange(-0.5, len(residue_labels), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(residue_labels), 1), minor=True)
        ax.grid(
            which="minor",
            color=settings.grid_color,
            linewidth=settings.grid_linewidth,
        )

    # ---------------- Colorbar ----------------
    if not settings.disable_colorbar:
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation=settings.cbar_orientation,
            shrink=settings.cbar_shrink,
            pad=settings.cbar_pad,
        )
        cbar.ax.tick_params(
            labelsize=settings.font_size_cbar,
        )

    # ---------------- Layout ----------------
    if settings.tight_layout:
        fig.tight_layout()

    # ---------------- Save ----------------
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
        bbox_inches="tight",
    )

    plt.close(fig)

    logger.info(f"Saved PPI heatmap → {out_path}")



def plot_protein_protein_interactions_stacked_column_from_file(
    pta_file,
    *,
    group_name: str,
    mode_name: str,
    settings: PPIStackedColumnSettings,
    out_dir: Path,
    attach_to_name: str = "mode",
    is_merged: bool = False,
) -> None:
    logger.debug(
        "plot_protein_protein_interactions_stacked_column_from_file: "
        "group='%s' mode='%s' is_merged=%s out_dir='%s'",
        group_name, mode_name, is_merged, out_dir,
    )

    # -------------------------------------------------
    # DATASET PATH
    # -------------------------------------------------
    if is_merged:
        dset_path = f"{group_name}/modes_merged/{mode_name}/table"
        fig_basename = f"merged_{settings.fig_basename}"
    else:
        dset_path = f"{group_name}/modes/{mode_name}/table"
        fig_basename = settings.fig_basename

    logger.debug("Dataset path: '%s'", dset_path)
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError("Dataset empty.")

    # -------------------------------------------------
    # INTERACTION NORMALIZATION (LIKE PLI)
    # -------------------------------------------------
    def norm_interaction(s: str) -> str:
        t = str(s).strip().upper()
        return {
            "HYDROPHOBIC": "hydrophobic",
            "HYDROGEN-BOND": "hydrogen_bond",
            "PI-CATION": "pi_cation",
            "PI-STACKING": "pi_stacking",
            "WATER-BRIDGE": "water_bridge_1",
            "WATER-BRIDGE-1": "water_bridge_1",
            "IONIC": "ionic",
            "HALOGEN-BOND": "halogen",
            "METAL-CONTACT": "metal_contact",
        }.get(t, t.lower().replace("-", "_"))

    # -------------------------------------------------
    # REPRESENTATION LOGIC
    # -------------------------------------------------
    chain_map = parse_range_dict(settings.alter_chains_str) \
        if settings.alter_chains and settings.alter_chains_str else {}

    seg_map = parse_range_dict(settings.alter_segments_str) \
        if settings.alter_segments and settings.alter_segments_str else {}

    def build_label(resname, resid, chainid, segid):

        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)

        if settings.renumber and settings.renumber_int is not None:
            resid = resid + settings.renumber_int

        if settings.alter_chains:
            for new_chain, (lo, hi) in chain_map.items():
                if lo <= resid <= hi:
                    chainid = new_chain
                    break

        if settings.alter_segments:
            for new_seg, (lo, hi) in seg_map.items():
                if lo <= resid <= hi:
                    segid = new_seg
                    break

        label = settings.representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))

        return label

    # -------------------------------------------------
    # COLLECT DATA
    # -------------------------------------------------
    pair_data = {}
    pair_errors = {}

    for row in dset:
        rec = json.loads(row)
        r1, r2, interaction_raw = rec["key"]

        if is_merged:
            value = float(rec["mean_frequency"])
            std = float(rec["std_frequency"])
        else:
            value = float(rec["frequency"])
            std = 0.0

        if value < settings.threshold:
            continue

        res1 = build_label(*r1)
        res2 = build_label(*r2)

        pair = f"{res1}–{res2}"
        ikey = norm_interaction(interaction_raw)

        pair_data.setdefault(pair, {})
        pair_data[pair][ikey] = pair_data[pair].get(ikey, 0.0) + value

        if is_merged:
            pair_errors.setdefault(pair, {})
            pair_errors[pair][ikey] = pair_errors[pair].get(ikey, 0.0) + std

    if not pair_data:
        logger.warning("No data after filtering.")
        return

    # -------------------------------------------------
    # BUILD MATRIX
    # -------------------------------------------------
    pairs = sorted(pair_data.keys())

    interactions = sorted({
        ikey
        for v in pair_data.values()
        for ikey in v.keys()
    })

    matrix = np.zeros((len(pairs), len(interactions)))
    errors = np.zeros_like(matrix)

    for i, pair in enumerate(pairs):
        for j, ikey in enumerate(interactions):
            matrix[i, j] = pair_data[pair].get(ikey, 0.0)
            if is_merged:
                errors[i, j] = pair_errors[pair].get(ikey, 0.0)

    # Drop empty rows
    mask = matrix.sum(axis=1) > 0
    matrix = matrix[mask]
    errors = errors[mask]
    pairs = [p for p, keep in zip(pairs, mask) if keep]

    if matrix.size == 0:
        logger.warning("No non-zero pairs remain.")
        return

    # -------------------------------------------------
    # STANDARD FIGSIZE (STRICT SETTINGS)
    # -------------------------------------------------
    fig_height = max(settings.fig_size_height, 10.0)  # hard minimum protection

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, fig_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )

    plt.rcParams["font.family"] = settings.font_family

    x = np.arange(len(pairs))
    bottom = np.zeros(len(pairs))

    # -------------------------------------------------
    # COLOR MAP FROM SETTINGS
    # -------------------------------------------------
    interaction_colors = {
        "hydrophobic": settings.color_hydrophobic,
        "hydrogen_bond": settings.color_hydrogen_bonds,
        "pi_cation": settings.color_pi_cation,
        "pi_stacking": settings.color_pi_stacking,
        "water_bridge_1": settings.color_water_bridge_1,
        "ionic": settings.color_ionic,
        "halogen": settings.color_halogen,
        "metal_contact": settings.color_metal_contact,
    }

    # -------------------------------------------------
    # STACKED BARS
    # -------------------------------------------------
    for j, ikey in enumerate(interactions):

        ax.bar(
            x,
            matrix[:, j],
            bottom=bottom,
            width=settings.bar_width,
            color=interaction_colors.get(ikey, "gray"),
            edgecolor=settings.bar_edge_color,
            alpha=settings.bar_alpha,
            linewidth=settings.bar_edge_width,
            label=ikey.replace("_", " ").title(),
        )

        bottom += matrix[:, j]

    # -------------------------------------------------
    # ERROR BARS
    # -------------------------------------------------
    if is_merged and settings.error_bars:
        total_errors = errors.sum(axis=1)

        ax.errorbar(
            x,
            bottom,
            yerr=total_errors,
            fmt="none",
            ecolor=settings.error_bars_color,
            elinewidth=settings.error_bars_line_width,
            capsize=settings.error_bars_capsize,
            alpha=settings.error_bars_alpha,
            linestyle=settings.error_bars_line_style,
        )

    # -------------------------------------------------
    # LABELS & TICKS
    # -------------------------------------------------
    if not settings.disable_title:
        ax.set_title(settings.fig_title,
                     fontsize=settings.font_size_title,
                     fontweight=settings.font_weight_title)

    if not settings.disable_x_axis:
        ax.set_xlabel(settings.x_label,
                      fontsize=settings.font_size_label,
                      fontweight=settings.font_weight_label)

    if not settings.disable_y_axis:
        ax.set_ylabel(settings.y_label,
                      fontsize=settings.font_size_label,
                      fontweight=settings.font_weight_label)

    if settings.disable_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.set_xticks(x)
        ax.set_xticklabels(
            pairs,
            rotation=settings.x_tick_rotation,
            fontsize=settings.font_size_x,
        )
        ax.tick_params(axis="y", labelsize=settings.font_size_y)

    if settings.enable_grid:
        ax.grid(True,
                linestyle=settings.grid_style,
                alpha=settings.grid_alpha,
                axis="y")

    if not settings.disable_legend:
        ax.legend(
            loc=settings.legend_loc,
            bbox_to_anchor=(0.5, settings.legend_bbox_y),
            ncol=min(settings.legend_n_col, len(interactions)),
            fontsize=settings.font_size_legend,
            frameon=settings.legend_frame,
            framealpha=settings.legend_alpha,
        )
        fig.subplots_adjust(bottom=settings.legend_margin_bottom)

    if settings.tight_layout:
        fig.tight_layout()

    # -------------------------------------------------
    # SAVE
    # -------------------------------------------------
    ext = "." + settings.fig_format
    out_path = Path(out_dir) / f"{attach_to_name}-{fig_basename}"
    out_path = out_path.with_suffix(ext)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_path,
        dpi=settings.fig_dpi,
        transparent=settings.fig_transparent,
        bbox_inches="tight",
    )

    plt.close(fig)

    logger.info(f"Saved PPI stacked column plot → {out_path}")
