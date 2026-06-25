"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Universal PTA time-series plot settings used for RMSD, angles, distances
(and any other per-frame scalar dataset stored under ``/<group>/frame_<N>``).

For per-atom data with no frame axis, see
:class:`pharmacon.constants.plots.rmsf.RMSFPlotSettings`.
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


__all__ = ["PlotUniversalSettings"]


@dataclass
class PlotUniversalSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PTA-UNIFIED",
        "PTA_UNIFIED",
    )

    # ---------------- FIGURE ----------------
    fig_size_width: float = 8
    fig_size_height: float = 5
    fig_dpi: int = 300
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    fig_title: str = "PTA Timeseries"
    x_label: str = "Time (ns)"
    y_label: str = "Value"
    x_axis: str = "time_ns"

    # ---------------- FONTS ----------------
    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_weight_title: str = "bold"
    font_size_label: int = 10
    font_weight_label: str = "normal"
    font_size_ticks: int = 8
    font_size_legend: int = 8
    font_weight_legend: str = "normal"

    # ---------------- LINE STYLE ----------------
    line_width: float = 1.5
    line_alpha: float = 1.0
    line_style: str = "solid"

    line_colors: List[str] = field(
        default_factory=lambda: [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
    )

    cycle_colors: bool = True

    # ---------------- STD BAND ----------------
    show_std_band: bool = True
    std_band_alpha: float = 0.25

    # ---------------- GRID ----------------
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.4

    # ---------------- LEGEND ----------------
    disable_legend: bool = False
    legend_frame: bool = True
    legend_alpha: float = 1.0

    # ---------------- BEHAVIOR ----------------
    plot_every_n: int = 1
    plot_multiple: bool = False
    disable_title: bool = False

    # VALIDATION
    def _validate_fields(self) -> None:

        # ---- figure ----
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 60)
        self.fig_size_height = self._safe_float(self.fig_size_height, 5, 2, 60)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # ---- labels ----
        self.fig_title = str(self.fig_title).strip()
        self.x_label = str(self.x_label).strip()
        self.y_label = str(self.y_label).strip()

        # ---- x axis ('frame' is accepted as an alias for 'frame_index') ----
        x_axis = str(self.x_axis).strip().lower()
        if x_axis == "frame":
            x_axis = "frame_index"
        if x_axis not in {"time_ps", "time_ns", "time_us", "frame_index"}:
            self._warn(f"Invalid x_axis '{self.x_axis}', using 'time_ns'")
            x_axis = "time_ns"
        self.x_axis = x_axis

        # ---- fonts ----
        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")
        self.font_size_title = self._safe_int(self.font_size_title, 14, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 10, 1)
        self.font_size_ticks = self._safe_int(self.font_size_ticks, 8, 1)
        self.font_size_legend = self._safe_int(self.font_size_legend, 8, 1)
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_weight_label = self._safe_font_weight(self.font_weight_label, "normal")
        self.font_weight_legend = self._safe_font_weight(self.font_weight_legend, "normal")

        # ---- lines ----
        self.line_width = self._safe_float(self.line_width, 1.5, 0.1, 10)
        self.line_alpha = self._safe_float(self.line_alpha, 1.0, 0.0, 1.0)
        self.std_band_alpha = self._safe_float(self.std_band_alpha, 0.25, 0.0, 1.0)

        self.line_style = str(self.line_style).strip().lower()
        if self.line_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid line_style '{self.line_style}', using 'solid'")
            self.line_style = "solid"

        # ---- grid ----
        self.enable_grid = self._safe_bool(self.enable_grid, True)
        self.grid_style = str(self.grid_style).strip().lower()
        if self.grid_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid grid_style '{self.grid_style}', using 'dashed'")
            self.grid_style = "dashed"
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.4, 0.0, 1.0)

        # ---- legend ----
        self.disable_legend = self._safe_bool(self.disable_legend, False)
        self.legend_frame = self._safe_bool(self.legend_frame, True)
        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0.0, 1.0)

        # ---- behaviour ----
        self.show_std_band = self._safe_bool(self.show_std_band, True)
        self.plot_multiple = self._safe_bool(self.plot_multiple, False)
        self.cycle_colors = self._safe_bool(self.cycle_colors, True)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.plot_every_n = self._safe_int(self.plot_every_n, 1, 1, 100000)

        # ---- line colors ----
        # A single color from the INI arrives as a plain string; split it so we
        # don't iterate the string character-by-character.
        if isinstance(self.line_colors, str):
            self.line_colors = [c.strip() for c in self.line_colors.split(",") if c.strip()]

        validated = []
        for c in self.line_colors:
            validated.append(self._safe_color(c, "#1f77b4"))

        if not validated:
            self._warn("line_colors empty, restoring defaults")
            validated = ["#1f77b4"]

        self.line_colors = validated

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


