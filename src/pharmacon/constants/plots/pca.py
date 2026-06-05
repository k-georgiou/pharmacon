"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Principal Component Analysis plot settings.

- ``PCAPlotTimeSeriesSettings``        — PC projections vs time.
- ``PCAPlotScatterSettings``           — 2-component scatter plot.
- ``PCAPlotVarianceRatioSettings``     — explained-variance scree plot.
- ``PCAPlotFESHeatmapSettings``        — free-energy surface (−kT·ln P) map.
- ``PCAPlotProbabilityHeatmapSettings``— probability density map.
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
    "PCAPlotTimeSeriesSettings",
    "PCAPlotScatterSettings",
    "PCAPlotVarianceRatioSettings",
    "PCAPlotFESHeatmapSettings",
    "PCAPlotProbabilityHeatmapSettings",
]


@dataclass
class PCAPlotTimeSeriesSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PCA-TIMESERIES",
        "PCA_TIMESERIES",

    )

    # ---- figure ----
    fig_size_width: float = 10
    fig_size_height: float = 6
    fig_dpi: int = 300
    fig_basename: str = "pca_timeseries"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # ---- labels ----
    fig_title: str = "PCA Time Series"
    x_label: str = "Time (ps)"
    y_label: str = "Projection Value"

    # ---- fonts ----
    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_ticks: int = 10
    font_size_legend: int = 10
    font_weight_title: str = "bold"

    # ---- line styling ----
    line_width: float = 1.5
    line_alpha: float = 0.9

    # ---- PC selection ----
    pcs: List[int] = field(default_factory=lambda: [1, 2, 3])

    # ---- legend ----
    disable_legend: bool = False
    legend_loc: str = "best"
    legend_frame: bool = True
    legend_alpha: float = 1.0

    # ---- grid ----
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.3

    disable_title: bool = False

    # -----------------------------------------------------

    def _validate_fields(self) -> None:

        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 10, 2, 40)
        self.fig_size_height = self._safe_float(self.fig_size_height, 6, 2, 40)

        self.fig_format = str(self.fig_format).lower()
        if f".{self.fig_format}" not in VALID_EXTENSIONS:
            self._warn("Invalid fig_format — using png")
            self.fig_format = "png"

        self.bg_color = self._safe_color(self.bg_color, "white")

        self.font_family = (
            self.font_family.lower()
            if self.font_family.lower() in AVAILABLE_FONTS
            else "dejavu sans"
        )

        self.font_size_title = self._safe_int(self.font_size_title, 14, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 12, 1)
        self.font_size_ticks = self._safe_int(self.font_size_ticks, 10, 1)
        self.font_size_legend = self._safe_int(self.font_size_legend, 10, 1)

        self.line_width = self._safe_float(self.line_width, 1.5, 0.1, 10)
        self.line_alpha = self._safe_float(self.line_alpha, 0.9, 0, 1)

        self.legend_loc = (
            self.legend_loc.lower()
            if self.legend_loc.lower() in VALID_LEGEND_LOCS
            else "best"
        )

        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0, 1)
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.3, 0, 1)

        if self.grid_style not in VALID_LINE_STYLES:
            self.grid_style = "dashed"

        self.pcs = [int(pc) for pc in self.pcs if int(pc) > 0]

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


@dataclass
class PCAPlotScatterSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PCA-SCATTER",
        "PCA_SCATTER",
    )

    # ---- figure ----
    fig_size_width: float = 8
    fig_size_height: float = 8
    fig_dpi: int = 300
    fig_basename: str = "pca_scatter"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # ---- labels ----
    fig_title: str = "PCA Scatter"
    x_label: str = "PC1"
    y_label: str = "PC2"

    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_ticks: int = 10

    font_weight_title: str = "bold"

    # ---- PCA components ----
    pc_x: int = 1
    pc_y: int = 2

    # ---- scatter styling ----
    scatter_size: float = 20
    scatter_alpha: float = 0.8
    cmap: str = "viridis"

    disable_colorbar: bool = False

    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.3

    disable_title: bool = False

    # -----------------------------------------------------

    def _validate_fields(self) -> None:

        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 40)
        self.fig_size_height = self._safe_float(self.fig_size_height, 8, 2, 40)

        self.bg_color = self._safe_color(self.bg_color, "white")

        self.scatter_size = self._safe_float(self.scatter_size, 20, 1, 1000)
        self.scatter_alpha = self._safe_float(self.scatter_alpha, 0.8, 0, 1)

        self.pc_x = max(1, int(self.pc_x))
        self.pc_y = max(1, int(self.pc_y))

        if self.grid_style not in VALID_LINE_STYLES:
            self.grid_style = "dashed"

        self.grid_alpha = self._safe_float(self.grid_alpha, 0.3, 0, 1)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


@dataclass
class PCAPlotVarianceRatioSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PCA-VARIANCE-RATIO",
        "PCA_VARIANCE_RATIO",
        "PCA-SCREE",
        "PCA_SCREE",
    )

    # ---- figure ----
    fig_size_width: float = 8
    fig_size_height: float = 6
    fig_dpi: int = 300
    fig_basename: str = "pca_variance_ratio"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # ---- labels ----
    fig_title: str = "PCA Explained Variance"
    x_label: str = "Principal Component"
    y_label: str = "Explained Variance Ratio"

    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_ticks: int = 10
    font_weight_title: str = "bold"

    # ---- style ----
    bar_alpha: float = 0.8
    cumulative_line: bool = True
    cumulative_line_width: float = 2.0
    cumulative_alpha: float = 0.9

    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.3

    disable_title: bool = False

    # -----------------------------------------------------

    def _validate_fields(self) -> None:

        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 40)
        self.fig_size_height = self._safe_float(self.fig_size_height, 6, 2, 40)

        self.bg_color = self._safe_color(self.bg_color, "white")

        self.bar_alpha = self._safe_float(self.bar_alpha, 0.8, 0, 1)
        self.cumulative_line_width = self._safe_float(self.cumulative_line_width, 2.0, 0.1, 10)
        self.cumulative_alpha = self._safe_float(self.cumulative_alpha, 0.9, 0, 1)

        if self.grid_style not in VALID_LINE_STYLES:
            self.grid_style = "dashed"

        self.grid_alpha = self._safe_float(self.grid_alpha, 0.3, 0, 1)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))



@dataclass
class PCAPlotFESHeatmapSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PCA-FES-HEATMAP",
        "PCA_FES_HEATMAP",
        "PCA FES HEATMAP",
        "PCA FREE ENERGY HEATMAP",
    )

    # ---- components ----
    components: Tuple[int, int] | list[int] | list[tuple[int, int]] | None = None
    plot_multiple: bool = True
    plot_in_one: bool = True

    # ---- calculation ----
    bins: int = 60
    smooth_sigma: float = 1.0
    temperature: float = 310.0

    # ---- contour / heatmap ----
    cmap: str = "plasma"
    n_levels: int = 12

    # ---- figure ----
    fig_size_width: float = 8
    fig_size_height: float = 6
    fig_dpi: int = 600
    fig_basename: str = "pca_heatmap_fes"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # ---- labels ----
    fig_title: str = ""
    x_label: str = ""
    y_label: str = ""

    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_ticks: int = 10
    font_weight_title: str = "bold"

    # ---- style ----
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.3

    disable_title: bool = False

    # -----------------------------------------------------

    def _validate_fields(self) -> None:

        # ---- figure ----
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 40)
        self.fig_size_height = self._safe_float(self.fig_size_height, 6, 2, 40)

        self.bg_color = self._safe_color(self.bg_color, "white")

        # ---- bins / smoothing ----
        self.bins = self._safe_int(self.bins, 60, 10, 500)
        self.smooth_sigma = self._safe_float(self.smooth_sigma, 1.0, 0.0, 10.0)
        self.temperature = self._safe_float(self.temperature, 300.0, 1.0, 1000.0)

        # ---- contour ----
        self.n_levels = self._safe_int(self.n_levels, 12, 3, 100)

        # ---- grid ----
        if self.grid_style not in VALID_LINE_STYLES:
            self.grid_style = "dashed"

        self.grid_alpha = self._safe_float(self.grid_alpha, 0.3, 0, 1)

        # ---- components sanity ----
        if isinstance(self.components, tuple):
            if len(self.components) != 2:
                self.components = (1, 2)

        # ---- safety ----
        if not isinstance(self.plot_multiple, bool):
            self.plot_multiple = False

        if not isinstance(self.plot_in_one, bool):
            self.plot_in_one = False

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))



@dataclass
class PCAPlotProbabilityHeatmapSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PCA-PROBABILITY-HEATMAP",
        "PCA_PROBABILITY_HEATMAP",
        "PCA PROBABILITY HEATMAP",
    )

    # ---- components ----
    components: Tuple[int, int] | list[int] | list[tuple[int, int]] | None = None
    plot_multiple: bool = True
    plot_in_one: bool = True

    # ---- calculation ----
    bins: int = 60
    smooth_sigma: float = 1.0

    # ---- contour / heatmap ----
    cmap: str = "viridis"
    n_levels: int = 12

    # ---- figure ----
    fig_size_width: float = 8
    fig_size_height: float = 6
    fig_dpi: int = 600
    fig_basename: str = "pca_heatmap_probability"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # ---- labels ----
    fig_title: str = ""
    x_label: str = ""
    y_label: str = ""

    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_ticks: int = 10
    font_weight_title: str = "bold"

    # ---- style ----
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.3

    disable_title: bool = False

    # -----------------------------------------------------

    def _validate_fields(self) -> None:

        # ---- figure ----
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 40)
        self.fig_size_height = self._safe_float(self.fig_size_height, 6, 2, 40)

        self.bg_color = self._safe_color(self.bg_color, "white")

        # ---- bins / smoothing ----
        self.bins = self._safe_int(self.bins, 60, 10, 500)
        self.smooth_sigma = self._safe_float(self.smooth_sigma, 1.0, 0.0, 10.0)

        # ---- contour ----
        self.n_levels = self._safe_int(self.n_levels, 12, 3, 100)

        # ---- grid ----
        if self.grid_style not in VALID_LINE_STYLES:
            self.grid_style = "dashed"

        self.grid_alpha = self._safe_float(self.grid_alpha, 0.3, 0, 1)

        # ---- components sanity ----
        if isinstance(self.components, tuple):
            if len(self.components) != 2:
                self.components = (1, 2)

        # ---- safety ----
        if not isinstance(self.plot_multiple, bool):
            self.plot_multiple = False

        if not isinstance(self.plot_in_one, bool):
            self.plot_in_one = False

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))