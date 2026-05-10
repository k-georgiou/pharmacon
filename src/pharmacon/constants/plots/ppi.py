"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Protein–Protein Interaction plot settings.

- ``ProteinProteinInteractionsTimelinePairs`` — residue-pair × frame timeline
  heatmap.
- ``ProteinProteinInteractionsHeatmap`` — residue × residue contact
  frequency heatmap.
- ``ProteinProteinInteractionsStackedColumnSettings`` — per residue-pair
  stacked column.
"""
import re

from ._base import (
    PlotSettingsBase,
    VALID_EXTENSIONS,
    VALID_FONT_WEIGHTS,
    VALID_LINE_STYLES,
    VALID_LEGEND_LOCS,
    AVAILABLE_FONTS,
    Namespace,
    ClassVar,
    Tuple,
    List,
    dataclass,
    field,
    fields,
    asdict,
    logger,
)
from matplotlib.colors import is_color_like


__all__ = [
    "ProteinProteinInteractionsTimelinePairs",
    "ProteinProteinInteractionsHeatmap",
    "ProteinProteinInteractionsStackedColumnSettings",
]


@dataclass
class ProteinProteinInteractionsTimelinePairs(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PPI-TIMELINE-PAIRS",
        "PPI_TIMELINE_PAIRS",
        "PP-TIMELINE",
        "PP_TIMELINE",
    )

    # Figure basics
    fig_size_width: float = 14
    fig_size_height: float = 8
    fig_dpi: int = 300
    fig_basename: str = "ppi_timeline_pairs"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Labels & Title
    fig_title: str = "Top Interaction Pairs – Timeline"
    x_label: str = "Frame"
    y_label: str = "Residue Pair"

    # Fonts
    font_size_x: int = 10
    font_weight_x: str = "normal"

    font_size_y: int = 10
    font_weight_y: str = "normal"

    font_size_label: int = 14
    font_weight_label: str = "bold"

    font_size_title: int = 18
    font_weight_title: str = "bold"

    font_size_cbar: int = 10
    font_weight_cbar: str = "normal"

    font_family: str = "dejavu sans"

    # Tick rotation
    x_tick_rotation: float = 0.0
    y_tick_rotation: float = 0.0

    # Data controls
    top_n: int = 50
    threshold: float = 0.02
    drop_empty_rows: bool = True

    # Representation system
    representation: str = "chainid:resnameresid"

    aa3_to_aa1: bool = False
    renumber: bool = False
    renumber_int: int | None = None

    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    # Colormap (Binary by default)
    cmap: str = "gnuplot"
    vmin: float = 0.0
    vmax: float = 1.0
    interpolation: str = "nearest"

    # Gridlines
    enable_grid: bool = False
    grid_color: str = "black"
    grid_linewidth: float = 0.2

    # Colorbar
    disable_colorbar: bool = False
    cbar_orientation: str = "vertical"
    cbar_shrink: float = 1.0
    cbar_pad: float = 0.04

    # Disable switches
    disable_title: bool = False
    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False
    disable_ticks: bool = False


@dataclass
class ProteinProteinInteractionsHeatmap(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PPI-HEATMAP",
        "PPI_HEATMAP",
        "PP-HEATMAP",
        "PP_HEATMAP",
    )

    # -------------------------------------------------
    # Figure basics
    # -------------------------------------------------
    fig_size_width: float = 16
    fig_size_height: float = 14
    fig_dpi: int = 300
    fig_basename: str = "ppi_heatmap"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # -------------------------------------------------
    # Title & Labels
    # -------------------------------------------------
    fig_title: str = "Protein–Protein Contact Frequency"
    x_label: str = "Residue"
    y_label: str = "Residue"

    font_size_title: int = 20
    font_weight_title: str = "bold"

    font_size_label: int = 16
    font_weight_label: str = "bold"

    font_size_x: int = 8
    font_weight_x: str = "normal"

    font_size_y: int = 8
    font_weight_y: str = "normal"

    font_size_cbar: int = 12
    font_weight_cbar: str = "normal"

    font_family: str = "dejavu sans"

    # -------------------------------------------------
    # Tick rotation
    # -------------------------------------------------
    x_tick_rotation: float = 90.0
    y_tick_rotation: float = 0.0

    # -------------------------------------------------
    # Data controls
    # -------------------------------------------------
    threshold: float = 0.02
    min_total: float = 0.0
    top_n: int = 0
    symmetric: bool = True

    # -------------------------------------------------
    # Representation
    # -------------------------------------------------
    representation: str = "chainid:resnameresid"
    aa3_to_aa1: bool = False
    renumber: bool = False
    renumber_int: int | None = None

    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    # -------------------------------------------------
    # Colormap
    # -------------------------------------------------
    cmap: str = "viridis"
    vmin: float = 0.0
    vmax: float | None = None  # None = auto-compute from data
    interpolation: str = "nearest"

    # -------------------------------------------------
    # Grid
    # -------------------------------------------------
    enable_grid: bool = False
    grid_color: str = "black"
    grid_linewidth: float = 0.2

    # -------------------------------------------------
    # Colorbar
    # -------------------------------------------------
    disable_colorbar: bool = False
    cbar_orientation: str = "vertical"
    cbar_shrink: float = 1.0
    cbar_pad: float = 0.04

    # -------------------------------------------------
    # Disable switches
    # -------------------------------------------------
    disable_title: bool = False
    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False
    disable_ticks: bool = False


@dataclass
class ProteinProteinInteractionsStackedColumnSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PPI-STACKED-COLUMN",
        "PPI_STACKED_COLUMN",
        "PP-STACKED",
        "PP_STACKED",
    )

    # Figure
    fig_size_width: float = 30
    fig_size_height: float = 4
    fig_dpi: int = 300
    fig_basename: str = "ppi_stacked_column"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = False
    bg_color: str = "white"

    # Labels
    fig_title: str = "Protein–Protein Interactions"
    x_label: str = "Residue–Residue Pair"
    y_label: str = "Interaction Frequency"

    # Fonts
    font_size_x: int = 10
    font_size_y: int = 12
    font_size_label: int = 10
    font_size_title: int = 24
    font_size_legend: int = 8
    font_weight_title: str = "bold"
    font_weight_label: str = "bold"
    font_family: str = "dejavu sans"

    # Tick rotation
    x_tick_rotation: float = 90.0
    y_tick_rotation: float = 0.0

    # Threshold
    threshold: float = 0.05

    # Axis limits
    y_limit_min: float | None = None
    y_limit_max: float | None = None

    # Bars
    bar_width: float = 0.9
    bar_edge_color: str = "black"
    bar_edge_width: float | None = None
    bar_alpha: float = 1.0

    # Grid
    enable_grid: bool = True
    grid_style: str = "--"
    grid_alpha: float = 0.5

    # Legend
    disable_legend: bool = True
    legend_loc: str = "upper center"
    legend_bbox_y: float = -0.20
    legend_n_col: int = 5
    legend_frame: bool = False
    legend_alpha: float = 1.0
    legend_margin_bottom: float = 0.28

    # Representation
    representation: str = "chainid:resnameresid"
    aa3_to_aa1: bool = True
    renumber: bool = False
    renumber_int: int | None = None
    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    # Disable switches
    disable_title: bool = False
    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_ticks: bool = False

    # for merged:
    error_bars: bool = True
    error_bars_capsize: float = 3
    error_bars_color: str = "black"
    error_bars_alpha: float = 0.5
    error_bars_line_width: float = 0.5
    error_bars_line_style: str = "solid"

    color_hydrophobic: str = "#b39ddb"
    color_hydrogen_bonds: str = "#f1c40f"
    color_pi_cation: str = "#e67e22"
    color_pi_stacking: str = "#2ecc71"
    color_water_bridge_1: str = "#3498db"
    color_ionic: str = "#ff3333"
    color_halogen: str = "#f49ac2"
    color_metal_contact: str = "#95a5a6"


