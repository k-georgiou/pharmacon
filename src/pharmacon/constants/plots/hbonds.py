"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Hydrogen-bond plot settings.

H-bonds are residue↔residue contacts of a single interaction type, so the
residue-pair plots reuse the protein–protein machinery (subclassed here only
to give H-bond-specific aliases, titles and basenames). Two H-bond-specific
summaries are added on top:

- ``HBondsCountPerFrame``  — number of H-bonds vs frame (time series)
- ``HBondsOccupancy``      — ranked per-pair occupancy (horizontal bar)
"""
from ._base import (
    PlotSettingsBase,
    VALID_EXTENSIONS,
    VALID_LINE_STYLES,
    Namespace,
    ClassVar,
    Tuple,
    dataclass,
    asdict,
)

from .ppi import (
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsTimelinePairs,
)
from .universal import PlotUniversalSettings


__all__ = [
    "HBondsHeatmap",
    "HBondsTimelinePairs",
    "HBondsCountPerFrame",
    "HBondsOccupancy",
    "HBondsNetwork",
]


# ---------------------------------------------------------------------------
# Residue×residue H-bond frequency heatmap (reuses the PPI heatmap plotter).
# ---------------------------------------------------------------------------
@dataclass
class HBondsHeatmap(ProteinProteinInteractionsHeatmap):

    alias: ClassVar[Tuple[str, ...]] = (
        "HBONDS-HEATMAP",
        "HBONDS_HEATMAP",
        "HB-HEATMAP",
        "HB_HEATMAP",
    )

    fig_basename: str = "hbonds_heatmap"
    fig_title: str = "Hydrogen Bond Contact Frequency"
    x_label: str = "Residue"
    y_label: str = "Residue"


# ---------------------------------------------------------------------------
# Per-pair occupancy timeline (reuses the PPI timeline-pairs plotter).
# ---------------------------------------------------------------------------
@dataclass
class HBondsTimelinePairs(ProteinProteinInteractionsTimelinePairs):

    alias: ClassVar[Tuple[str, ...]] = (
        "HBONDS-TIMELINE-PAIRS",
        "HBONDS_TIMELINE_PAIRS",
        "HB-TIMELINE",
        "HB_TIMELINE",
    )

    fig_basename: str = "hbonds_timeline_pairs"
    fig_title: str = "Hydrogen Bond Pairs – Timeline"
    y_label: str = "Residue Pair"


# ---------------------------------------------------------------------------
# Number of H-bonds per frame (time series; reuses universal line settings).
# ---------------------------------------------------------------------------
@dataclass
class HBondsCountPerFrame(PlotUniversalSettings):

    alias: ClassVar[Tuple[str, ...]] = (
        "HBONDS-COUNT-PER-FRAME",
        "HBONDS_COUNT_PER_FRAME",
        "HB-COUNT",
        "HB_COUNT",
    )

    fig_title: str = "Hydrogen Bonds per Frame"
    # H-bond frame datasets carry no time metadata, so default to the frame
    # index (the universal x_axis validator accepts 'frame_index').
    x_axis: str = "frame_index"
    x_label: str = "Frame"
    y_label: str = "Number of H-bonds"


# ---------------------------------------------------------------------------
# Ranked per-pair occupancy (horizontal bar).
# ---------------------------------------------------------------------------
@dataclass
class HBondsOccupancy(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "HBONDS-OCCUPANCY",
        "HBONDS_OCCUPANCY",
        "HB-OCCUPANCY",
        "HB_OCCUPANCY",
    )

    # Figure
    fig_size_width: float = 8
    fig_size_height: float = 10
    fig_dpi: int = 300
    fig_basename: str = "hbonds_occupancy"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Title / labels
    fig_title: str = "Hydrogen Bond Occupancy"
    x_label: str = "Occupancy (fraction of frames)"
    y_label: str = "Residue Pair"

    # Fonts
    font_family: str = "dejavu sans"
    font_size_title: int = 14
    font_weight_title: str = "bold"
    font_size_label: int = 12
    font_weight_label: str = "bold"
    font_size_ticks: int = 8

    # Data selection
    top_n: int = 25            # keep the N most-occupied pairs (0 = all)
    threshold: float = 0.0     # drop pairs below this occupancy

    # Bars
    bar_color: str = "#3498db"
    bar_alpha: float = 0.9
    bar_edge_color: str = "black"
    bar_edge_width: float = 0.5

    # Grid
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.4

    # Labeling
    representation: str = "chainid:resnameresid"
    aa3_to_aa1: bool = True

    # Disable switches
    disable_title: bool = False
    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_ticks: bool = False

    def _validate_fields(self) -> None:
        self.fig_size_width = self._safe_float(self.fig_size_width, 8.0, 1.0, 200.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 10.0, 1.0, 200.0)
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

        self.fig_title = str(self.fig_title)
        self.x_label = str(self.x_label)
        self.y_label = str(self.y_label)

        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")
        self.font_size_title = self._safe_int(self.font_size_title, 14, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 12, 1)
        self.font_size_ticks = self._safe_int(self.font_size_ticks, 8, 1)
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_weight_label = self._safe_font_weight(self.font_weight_label, "bold")

        self.top_n = self._safe_int(self.top_n, 25, 0)
        self.threshold = self._safe_float(self.threshold, 0.0, 0.0, 1.0)

        self.bar_color = self._safe_color(self.bar_color, "#3498db")
        self.bar_alpha = self._safe_float(self.bar_alpha, 0.9, 0.0, 1.0)
        self.bar_edge_color = self._safe_color(self.bar_edge_color, "black")
        self.bar_edge_width = self._safe_float(self.bar_edge_width, 0.5, 0.0, 5.0)

        self.enable_grid = self._safe_bool(self.enable_grid, True)
        self.grid_style = str(self.grid_style).strip().lower()
        if self.grid_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid grid_style '{self.grid_style}', using 'dashed'")
            self.grid_style = "dashed"
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.4, 0.0, 1.0)

        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, True)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


# ---------------------------------------------------------------------------
# H-bond network graph (nodes = residues, edges = H-bonds weighted by occupancy).
# ---------------------------------------------------------------------------
@dataclass
class HBondsNetwork(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "HBONDS-NETWORK",
        "HBONDS_NETWORK",
        "HB-NETWORK",
        "HB_NETWORK",
    )

    # Figure
    fig_size_width: float = 12
    fig_size_height: float = 12
    fig_dpi: int = 300
    fig_basename: str = "hbonds_network"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Title
    fig_title: str = "Hydrogen Bond Network"
    disable_title: bool = False

    # Fonts
    font_family: str = "dejavu sans"
    font_size_title: int = 16
    font_weight_title: str = "bold"
    font_size_labels: int = 8       # node (residue) labels
    font_size_cbar: int = 8

    # Source mode (mode2 = once-per-frame ⇒ occupancy in [0, 1])
    mode: str = "mode2"

    # Data selection (defaults tuned for large proteins)
    threshold: float = 0.3          # keep only bonds at/above this occupancy
    top_n: int = 200                # keep the N strongest edges (0 = all)
    min_seq_sep: int = 0            # drop |i-j| <= this within a chain (0 = off)
    largest_component: bool = False  # keep only the largest connected component

    # Labeling
    representation: str = "chainid:resnameresid"
    aa3_to_aa1: bool = True
    label_top_n: int = 15           # label the N highest-degree hubs (0 = none, -1 = all)

    # Layout
    layout: str = "spring"          # spring | kamada_kawai | circular
    seed: int = 42

    # Nodes (size scales with degree)
    node_size_min: float = 50
    node_size_max: float = 600
    node_color: str = "#cfe8ff"
    node_edge_color: str = "#1f4e79"

    # Edges (colored by occupancy, width scales with occupancy)
    edge_cmap: str = "viridis"
    edge_width_min: float = 0.5
    edge_width_max: float = 4.0
    edge_alpha: float = 0.7

    # Colorbar
    disable_colorbar: bool = False

    def _validate_fields(self) -> None:
        self.fig_size_width = self._safe_float(self.fig_size_width, 12.0, 1.0, 200.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 12.0, 1.0, 200.0)
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

        self.fig_title = str(self.fig_title)
        self.disable_title = self._safe_bool(self.disable_title, False)

        self.font_family = self._safe_font_family(self.font_family, "dejavu sans")
        self.font_size_title = self._safe_int(self.font_size_title, 16, 1)
        self.font_weight_title = self._safe_font_weight(self.font_weight_title, "bold")
        self.font_size_labels = self._safe_int(self.font_size_labels, 8, 1)
        self.font_size_cbar = self._safe_int(self.font_size_cbar, 8, 1)

        self.mode = str(self.mode).strip() or "mode2"

        self.threshold = self._safe_float(self.threshold, 0.3, 0.0, 1.0)
        self.top_n = self._safe_int(self.top_n, 200, 0)
        self.min_seq_sep = self._safe_int(self.min_seq_sep, 0, 0)
        self.largest_component = self._safe_bool(self.largest_component, False)

        rep = str(self.representation).strip().lower()
        self.representation = rep or "chainid:resnameresid"
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, True)
        self.label_top_n = self._safe_int(self.label_top_n, 15, -1)

        self.layout = str(self.layout).strip().lower()
        if self.layout not in {"spring", "kamada_kawai", "circular"}:
            self._warn(f"Invalid layout '{self.layout}', using 'spring'")
            self.layout = "spring"
        self.seed = self._safe_int(self.seed, 42)

        self.node_size_min = self._safe_float(self.node_size_min, 50.0, 0.0)
        self.node_size_max = self._safe_float(self.node_size_max, 600.0, 0.0)
        self.node_color = self._safe_color(self.node_color, "#cfe8ff")
        self.node_edge_color = self._safe_color(self.node_edge_color, "#1f4e79")

        self.edge_cmap = self._safe_cmap(self.edge_cmap, "viridis")
        self.edge_width_min = self._safe_float(self.edge_width_min, 0.5, 0.0, 20.0)
        self.edge_width_max = self._safe_float(self.edge_width_max, 4.0, 0.0, 20.0)
        self.edge_alpha = self._safe_float(self.edge_alpha, 0.7, 0.0, 1.0)

        self.disable_colorbar = self._safe_bool(self.disable_colorbar, False)

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))
