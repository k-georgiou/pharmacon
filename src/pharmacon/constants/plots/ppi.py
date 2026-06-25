"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

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


# Token grammar shared by every PPI ``representation`` string, e.g.
# "chainid:resnameresid". Mirrors the PLI x_axis_representation validator.
_REP_PATTERN = re.compile(
    r"^(resname|resid|chainid|segid)"
    r"(([:\-_]?)(resname|resid|chainid|segid))*$"
)


def _validate_representation(settings: "PlotSettingsBase") -> str:
    """Validate ``settings.representation`` against the token grammar,
    coercing an invalid value to the default with a warning."""
    rep = str(settings.representation).strip().lower()
    if not _REP_PATTERN.fullmatch(rep):
        settings._warn(
            f"Invalid representation '{settings.representation}', "
            "using 'chainid:resnameresid'"
        )
        return "chainid:resnameresid"
    return rep


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
    xtick_max: int = 50            # thin frame ticks to ~this many when exceeded

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

    def _validate_fields(self) -> None:
        # Figure basics
        self.fig_size_width = self._safe_float(self.fig_size_width, 14.0, 1.0, 200.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 8.0, 1.0, 200.0)
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # Labels
        self.fig_title = str(self.fig_title).strip()
        self.x_label = str(self.x_label).strip()
        self.y_label = str(self.y_label).strip()

        # Font sizes
        self.font_size_x = self._safe_int(self.font_size_x, 10, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 10, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 14, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 18, 1)
        self.font_size_cbar = self._safe_int(self.font_size_cbar, 10, 1)

        # Font weights
        self.font_weight_x = self._safe_font_weight(self.font_weight_x, "normal")
        self.font_weight_y = self._safe_font_weight(self.font_weight_y, "normal")
        self.font_weight_label = self._safe_font_weight(self.font_weight_label, "bold")
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_weight_cbar = self._safe_font_weight(self.font_weight_cbar, "normal")

        # Font family
        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")

        # Rotations
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 0.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # Data controls
        self.top_n = self._safe_int(self.top_n, 50, 0)
        self.threshold = self._safe_float(self.threshold, 0.02, 0.0, 1.0)
        self.drop_empty_rows = self._safe_bool(self.drop_empty_rows, True)
        self.xtick_max = self._safe_int(self.xtick_max, 50, 1)

        # Representation / labeling
        self.representation = _validate_representation(self)
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, False)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)
        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        # Colormap
        self.cmap = self._safe_cmap(self.cmap, "gnuplot")
        self.vmin = self._safe_float(self.vmin, 0.0)
        self.vmax = self._safe_float(self.vmax, 1.0)
        self.interpolation = self._safe_interpolation(self.interpolation, "nearest")

        # Grid
        self.enable_grid = self._safe_bool(self.enable_grid, False)
        self.grid_color = self._safe_color(self.grid_color, "black")
        self.grid_linewidth = self._safe_float(self.grid_linewidth, 0.2, 0.0, 10.0)

        # Colorbar
        self.disable_colorbar = self._safe_bool(self.disable_colorbar, False)
        self.cbar_orientation = str(self.cbar_orientation).strip().lower() or "vertical"
        if self.cbar_orientation not in {"vertical", "horizontal"}:
            self._warn("Invalid cbar_orientation, using 'vertical'")
            self.cbar_orientation = "vertical"
        self.cbar_shrink = self._safe_float(self.cbar_shrink, 1.0, 0.1, 2.0)
        self.cbar_pad = self._safe_float(self.cbar_pad, 0.04, 0.0, 1.0)

        # Disable switches
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


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

    def _validate_fields(self) -> None:
        # Figure basics
        self.fig_size_width = self._safe_float(self.fig_size_width, 16.0, 1.0, 200.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 14.0, 1.0, 200.0)
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # Labels
        self.fig_title = str(self.fig_title).strip()
        self.x_label = str(self.x_label).strip()
        self.y_label = str(self.y_label).strip()

        # Font sizes
        self.font_size_title = self._safe_int(self.font_size_title, 20, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 16, 1)
        self.font_size_x = self._safe_int(self.font_size_x, 8, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 8, 1)
        self.font_size_cbar = self._safe_int(self.font_size_cbar, 12, 1)

        # Font weights
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_weight_label = self._safe_font_weight(self.font_weight_label, "bold")
        self.font_weight_x = self._safe_font_weight(self.font_weight_x, "normal")
        self.font_weight_y = self._safe_font_weight(self.font_weight_y, "normal")
        self.font_weight_cbar = self._safe_font_weight(self.font_weight_cbar, "normal")

        # Font family
        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")

        # Rotations
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # Data controls
        self.threshold = self._safe_float(self.threshold, 0.02, 0.0, 1.0)
        self.min_total = self._safe_float(self.min_total, 0.0, 0.0)
        self.top_n = self._safe_int(self.top_n, 0, 0)
        self.symmetric = self._safe_bool(self.symmetric, True)

        # Representation / labeling
        self.representation = _validate_representation(self)
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, False)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)
        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        # Colormap (vmax empty/unset → None = auto from data)
        self.cmap = self._safe_cmap(self.cmap, "viridis")
        self.vmin = self._safe_float(self.vmin, 0.0)
        self.vmax = self._safe_float(self.vmax, None)
        self.interpolation = self._safe_interpolation(self.interpolation, "nearest")

        # Grid
        self.enable_grid = self._safe_bool(self.enable_grid, False)
        self.grid_color = self._safe_color(self.grid_color, "black")
        self.grid_linewidth = self._safe_float(self.grid_linewidth, 0.2, 0.0, 10.0)

        # Colorbar
        self.disable_colorbar = self._safe_bool(self.disable_colorbar, False)
        self.cbar_orientation = str(self.cbar_orientation).strip().lower() or "vertical"
        if self.cbar_orientation not in {"vertical", "horizontal"}:
            self._warn("Invalid cbar_orientation, using 'vertical'")
            self.cbar_orientation = "vertical"
        self.cbar_shrink = self._safe_float(self.cbar_shrink, 1.0, 0.1, 2.0)
        self.cbar_pad = self._safe_float(self.cbar_pad, 0.04, 0.0, 1.0)

        # Disable switches
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


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

    def _validate_fields(self) -> None:
        # Figure basics
        self.fig_size_width = self._safe_float(self.fig_size_width, 30.0, 1.0, 200.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 4.0, 1.0, 200.0)
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, False)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # Labels
        self.fig_title = str(self.fig_title).strip()
        self.x_label = str(self.x_label).strip()
        self.y_label = str(self.y_label).strip()

        # Font sizes
        self.font_size_x = self._safe_int(self.font_size_x, 10, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 12, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 10, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 24, 1)
        self.font_size_legend = self._safe_int(self.font_size_legend, 8, 1)

        # Font weights
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_weight_label = self._safe_font_weight(self.font_weight_label, "bold")

        # Font family
        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")

        # Rotations
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # Threshold
        self.threshold = self._safe_float(self.threshold, 0.05, 0.0, 1.0)

        # Axis limits (empty/unset → None = auto)
        self.y_limit_min = self._safe_float(self.y_limit_min, None)
        self.y_limit_max = self._safe_float(self.y_limit_max, None)

        # Bars (bar_edge_width empty/unset → None = matplotlib default)
        self.bar_width = self._safe_float(self.bar_width, 0.9, 0.0)
        self.bar_edge_color = self._safe_color(self.bar_edge_color, "black")
        self.bar_edge_width = self._safe_float(self.bar_edge_width, None)
        self.bar_alpha = self._safe_float(self.bar_alpha, 1.0, 0.0, 1.0)

        # Grid
        self.enable_grid = self._safe_bool(self.enable_grid, True)
        self.grid_style = str(self.grid_style).strip().lower()
        if self.grid_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid grid_style '{self.grid_style}', using '--'")
            self.grid_style = "--"
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.5, 0.0, 1.0)

        # Legend
        self.disable_legend = self._safe_bool(self.disable_legend, True)
        self.legend_loc = str(self.legend_loc).strip().lower()
        if self.legend_loc not in VALID_LEGEND_LOCS:
            self._warn(f"Invalid legend_loc '{self.legend_loc}', using 'upper center'")
            self.legend_loc = "upper center"
        self.legend_bbox_y = self._safe_float(self.legend_bbox_y, -0.20, -2.0, 2.0)
        self.legend_n_col = self._safe_int(self.legend_n_col, 5, 1, 20)
        self.legend_frame = self._safe_bool(self.legend_frame, False)
        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0.0, 1.0)
        self.legend_margin_bottom = self._safe_float(self.legend_margin_bottom, 0.28, 0.0, 1.0)

        # Representation / labeling
        self.representation = _validate_representation(self)
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, True)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)
        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        # Disable switches
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

        # Error bars (merged only)
        self.error_bars = self._safe_bool(self.error_bars, True)
        self.error_bars_capsize = self._safe_float(self.error_bars_capsize, 3.0, 0.0, 20.0)
        self.error_bars_color = self._safe_color(self.error_bars_color, "black")
        self.error_bars_alpha = self._safe_float(self.error_bars_alpha, 0.5, 0.0, 1.0)
        self.error_bars_line_width = self._safe_float(self.error_bars_line_width, 0.5, 0.0, 5.0)
        self.error_bars_line_style = str(self.error_bars_line_style).strip().lower()
        if self.error_bars_line_style not in VALID_LINE_STYLES:
            self._warn(
                f"Invalid error_bars_line_style '{self.error_bars_line_style}', using 'solid'"
            )
            self.error_bars_line_style = "solid"

        # Interaction colors
        self.color_hydrophobic = self._safe_color(self.color_hydrophobic, "#b39ddb")
        self.color_hydrogen_bonds = self._safe_color(self.color_hydrogen_bonds, "#f1c40f")
        self.color_pi_cation = self._safe_color(self.color_pi_cation, "#e67e22")
        self.color_pi_stacking = self._safe_color(self.color_pi_stacking, "#2ecc71")
        self.color_water_bridge_1 = self._safe_color(self.color_water_bridge_1, "#3498db")
        self.color_ionic = self._safe_color(self.color_ionic, "#ff3333")
        self.color_halogen = self._safe_color(self.color_halogen, "#f49ac2")
        self.color_metal_contact = self._safe_color(self.color_metal_contact, "#95a5a6")

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


