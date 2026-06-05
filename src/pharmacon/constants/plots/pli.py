"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Protein–Ligand Interaction plot settings.

Contains one dataclass per PLI plot variant:

- ``ProteinLigandInteractionsStackedColumn1PlotSettings`` — stacked bar
  per-residue × interaction type.
- ``ProteinLigandInteractionsStackedColumn2PlotSettings`` — stacked bar
  per-residue split by backbone / side-chain.
- ``ProteinLigandInteractionsHeatmap1PlotSettings`` — residue × frame
  aggregated-intensity heatmap.
- ``ProteinLigandInteractionsHeatmap2PlotSettings`` — interaction × frame
  auto-sized heatmap.
- ``ProteinLigandInteractionsPieCharts1PlotSettings`` — per-residue pie
  charts with optional collage.
- ``ProteinLigandInteractionsLigandMonitorSettings`` — residue × ligand
  atom contact heatmap.
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
    "ProteinLigandInteractionsStackedColumn1PlotSettings",
    "ProteinLigandInteractionsStackedColumn2PlotSettings",
    "ProteinLigandInteractionsHeatmap1PlotSettings",
    "ProteinLigandInteractionsHeatmap2PlotSettings",
    "ProteinLigandInteractionsPieCharts1PlotSettings",
    "ProteinLigandInteractionsLigandMonitorSettings",
]


@dataclass
class ProteinLigandInteractionsStackedColumn1PlotSettings(PlotSettingsBase):
    """
    Settings related to the visualization of protein-ligand interaction data in
    stacked column chart format.

    The `ProteinLigandInteractionsStackedColumn1PlotSettings` class provides numerous
    configuration options to customize the appearance and properties of plots for protein-ligand
    interactions rendered as stacked column charts. These settings include options for figure
    size, fonts, colors, and other stylistic customizations necessary for clear and detailed
    representation of interaction data.

    Its attributes control visual elements including figure dimensions, axis labels, tick
    mark appearance, title formatting, legend options, grid styles, bar styles, and color
    themes for various interaction types. Additional options allow toggling labels, legend,
    grid visibility, and other rendering features.

    :ivar alias: A tuple of valid aliases for this type of plot setting.
    :type alias: ClassVar[Tuple[str, ...]]
    :ivar fig_size_width: The width of the figure in inches.
    :ivar fig_size_height: The height of the figure in inches.
    :ivar fig_dpi: The dots per inch (DPI) resolution for the figure.
    :ivar fig_basename: Basename for saving the figure.
    :ivar fig_format: File format for the figure (e.g. png, jpg).
    :ivar fig_transparent: Flag to save the figure with a transparent background.
    :ivar tight_layout: Whether to use the tight_layout option for better spacing.
    :ivar bg_color: Background color of the figure.

    :ivar fig_title: Title displayed on the figure.
    :ivar x_label: Label for the x-axis.
    :ivar y_label: Label for the y-axis.

    :ivar font_size_x: Font size for x-axis label.
    :ivar font_weight_x: Font weight for x-axis label (e.g., normal, bold).
    :ivar font_size_y: Font size for y-axis label.
    :ivar font_weight_y: Font weight for y-axis label (e.g., normal, bold).
    :ivar font_size_title: Font size for the figure title.
    :ivar font_weight_title: Font weight for the figure title.
    :ivar font_size_legend: Font size for the legend.
    :ivar font_weight_legend: Font weight for the legend.
    :ivar font_size_label: Font size for labels other than axes.
    :ivar font_weight_label: Font weight for labels other than axes.
    :ivar font_size_ticks: Font size for tick marks on axes.
    :ivar font_weight_ticks: Font weight for tick marks on axes.

    :ivar font_family: Name of the font family used for text in the figure.

    :ivar x_tick_rotation: Rotation angle for x-axis tick labels, in degrees.
    :ivar y_tick_rotation: Rotation angle for y-axis tick labels, in degrees.

    :ivar threshold: Minimum frequency threshold for displaying data.

    :ivar disable_x_axis: Flag to disable rendering the x-axis.
    :ivar disable_y_axis: Flag to disable rendering the y-axis.
    :ivar disable_x_label: Flag to disable rendering the x-axis label.
    :ivar disable_y_label: Flag to disable rendering the y-axis label.
    :ivar disable_title: Flag to disable rendering the title.
    :ivar disable_legend: Flag to disable rendering the legend.
    :ivar disable_ticks: Flag to disable rendering tick marks on axes.

    :ivar enable_grid: Flag to enable a grid on the plot.
    :ivar grid_style: Style of the grid lines (e.g., dashed, solid).
    :ivar grid_color: Color of the grid lines.
    :ivar grid_alpha: Transparency level for the grid lines.

    :ivar bar_width: Width of the bars in the stacked column chart.
    :ivar bar_edge_color: Edge color of the bars.
    :ivar bar_alpha: Transparency level for the bars.
    :ivar bar_edge_width: Edge width of the bars.

    :ivar y_limit_min: Minimum value for the y-axis scale.
    :ivar y_limit_max: Maximum value for the y-axis scale.

    :ivar legend_loc: Location of the legend on the plot.
    :ivar legend_bbox_y: Vertical position offset for the legend's bounding box.
    :ivar legend_n_col: Number of columns in the legend.
    :ivar legend_frame: Flag to display a frame around the legend.
    :ivar legend_alpha: Transparency level for the legend background.
    :ivar legend_margin_bottom: Additional margin at the bottom of the legend box.

    :ivar color_hydrophobic: Color for bars representing hydrophobic interactions.
    :ivar color_hydrogen_bonds: Color for bars representing hydrogen bonds.
    :ivar color_pi_cation: Color for bars representing pi-cation interactions.
    :ivar color_pi_stacking: Color for bars representing pi-stacking interactions.
    :ivar color_water_bridge_1: Color for bars representing water bridge interactions.
    :ivar color_ionic: Color for bars representing ionic interactions.
    :ivar color_halogen: Color for bars representing halogen interactions.
    :ivar color_metal_contact: Color for bars representing metal contact interactions.

    :ivar aa3_to_aa1: Flag to convert three-letter amino acid codes to one-letter codes.
    :ivar renumber: Flag to renumber residues instead of using original numbering.
    :ivar renumber_int: Starting integer for renumbering residues.
    :ivar alter_chains: Flag to modify chain identifiers in the plot.
    :ivar alter_segments: Flag to modify segment identifiers in the plot.
    :ivar alter_chains_str: String specification for altering chain identifiers.
    :ivar alter_segments_str: String specification for altering segment identifiers.

    :ivar x_axis_representation: Format for representing residue identifiers on the x-axis.

    :ivar error_bars: Flag to enable error bars in the plot.
    :ivar error_bars_capsize: Length of caps at the ends of error bars.
    :ivar error_bars_color: Color of error bars.
    :ivar error_bars_alpha: Transparency level for error bars.
    :ivar error_bars_line_width: Line width for error bars.
    :ivar error_bars_line_style: Line style for error bars (e.g., solid, dashed).
    """
    alias: ClassVar[Tuple[str, ...]] = ("PLI-STACKED-COLUMN-1",
                                        "PLI_STACKED_COLUMN_1",
                                        "PL-INTERACTIONS-STACKED-COLUMN-1",
                                        "PL-INTERACTIONS_STACKED_COLUMN_1",
                                        "MERGED-PLI-STACKED-COLUMN-1",
                                        "MERGED-PLI_STACKED_COLUMN_1",
                                        "MERGED-PL-INTERACTIONS-STACKED-COLUMN-1",
                                        "MERGED-PL-INTERACTIONS_STACKED_COLUMN_1")
    fig_size_width: float = 6
    fig_size_height: float = 4
    fig_dpi: int = 800
    fig_basename: str = "pli_stacked_column_1"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = False
    bg_color: str = "white"

    fig_title: str = "Protein-Ligand Interactions Stacked Column 1"
    x_label: str = "Residues"
    y_label: str = "Frequency"

    font_size_x: int = 8
    font_weight_x: str = "normal"

    font_size_y: int = 8
    font_weight_y: str = "normal"

    font_size_title: int = 14
    font_weight_title: str = "bold"

    font_size_legend: int = 8
    font_weight_legend: str = "normal"

    font_size_label: int = 8
    font_weight_label: str = "normal"

    font_size_ticks: int = 8
    font_weight_ticks: str = "normal"

    font_family: str = "dejavu sans"

    x_tick_rotation: int | float = 90
    y_tick_rotation: int | float = 0

    threshold: float = 0.02

    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False

    disable_title: bool = False
    disable_legend: bool = False
    disable_ticks: bool = False

    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_color: str = "lightgray"
    grid_alpha: float = 0.5

    bar_width: float = 0.8
    bar_edge_color: str = "black"
    bar_alpha: float = 0.8
    bar_edge_width: float = 0.5

    y_limit_min: float | None = None
    y_limit_max: float | None = None

    legend_loc: str = "lower center"
    legend_bbox_y: float = -0.37
    legend_n_col: int = 4
    legend_frame: bool = True
    legend_alpha: float = 1.0
    legend_margin_bottom: float = 0.02

    color_hydrophobic: str = "#b39ddb"
    color_hydrogen_bonds: str = "#f1c40f"
    color_pi_cation: str = "#e67e22"
    color_pi_stacking: str = "#2ecc71"
    color_water_bridge_1: str = "#3498db"
    color_ionic: str = "#ff3333"
    color_halogen: str = "#f49ac2"
    color_metal_contact: str = "#95a5a6"

    aa3_to_aa1: bool = True
    renumber: bool = False
    renumber_int: int | None = None
    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    x_axis_representation: str = "chainid:resnameresid"

    # for merged:
    error_bars: bool = True
    error_bars_capsize: float = 3
    error_bars_color: str = "black"
    error_bars_alpha: float = 0.5
    error_bars_line_width: float = 0.5
    error_bars_line_style: str = "solid"

    def _validate_fields(self) -> None:
        """
        Validate and sanitize various configuration fields to ensure integrity and consistency.

        This method enforces limits, sets defaults, and will replace invalid entries with
        appropriate fallback values for graphical and display parameters such as figure
        dimensions, labels, font styles, grid display, bar styles, and legend settings.
        Additionally, it processes color definitions, error bar display attributes, and
        other related visualization options.

        :raises Warning: If invalid values are found during field validations, warnings
                         are raised explaining the fallback changes made.

        :return: None
        """

        # FIGURE BASICS
        self.fig_dpi = self._safe_int(self.fig_dpi, 800, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 6.0, 1.0, 40.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 4.0, 1.0, 40.0)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            ext = ".png"
        self.fig_format = ext.lstrip(".")

        self.fig_basename = str(self.fig_basename).strip() or "pli_stacked_column_1"

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, False)

        self.bg_color = self._safe_color(self.bg_color, "white")


        # LABELS
        self.fig_title = str(self.fig_title).strip() or \
                         "Protein-Ligand Interactions Stacked Column 1"

        self.x_label = str(self.x_label).strip() or "Residues"
        self.y_label = str(self.y_label).strip() or "Frequency"

        # FONT SIZES
        self.font_size_x = self._safe_int(self.font_size_x, 8, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 8, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 14, 1)
        self.font_size_legend = self._safe_int(self.font_size_legend, 8, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 8, 1)
        self.font_size_ticks = self._safe_int(self.font_size_ticks, 8, 1)

        # FONT WEIGHTS
        def _fw(value: str, default: str) -> str:
            """
            Converts and validates a given font weight value, ensuring it complies with the
            set of valid font weights. If the provided value is not valid, it issues a warning
            and falls back to a default value.

            :param value: The input font weight as a string.
            :param default: The default font weight to use in case the input value is invalid.
            :return: A validated and sanitized font weight string.
            """
            v = str(value).strip().lower()
            if v not in VALID_FONT_WEIGHTS:
                self._warn(f"Invalid font weight '{value}', using '{default}'")
                return default
            return v

        self.font_weight_x = _fw(self.font_weight_x, "normal")
        self.font_weight_y = _fw(self.font_weight_y, "normal")
        self.font_weight_title = _fw(self.font_weight_title, "bold")
        self.font_weight_legend = _fw(self.font_weight_legend, "normal")
        self.font_weight_label = _fw(self.font_weight_label, "normal")
        self.font_weight_ticks = _fw(self.font_weight_ticks, "normal")

        # FONT FAMILY
        font = str(self.font_family).strip().lower()
        if font not in AVAILABLE_FONTS:
            self._warn(f"Font '{font}' not found, using 'dejavu sans'")
            self.font_family = "dejavu sans"
        else:
            self.font_family = font

        # ROTATIONS
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # THRESHOLD
        self.threshold = self._safe_float(self.threshold, 0.02, 0.0, 1.0)

        # BOOLEAN FLAGS
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_legend = self._safe_bool(self.disable_legend, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)
        self.enable_grid = self._safe_bool(self.enable_grid, True)

        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, True)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)
        self.error_bars = self._safe_bool(self.error_bars, True)

        # GRID
        self.grid_style = str(self.grid_style).strip().lower()
        if self.grid_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid grid_style '{self.grid_style}', using 'dashed'")
            self.grid_style = "dashed"

        self.grid_color = self._safe_color(self.grid_color, "lightgray")
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.5, 0.0, 1.0)

        # BARS
        self.bar_width = self._safe_float(self.bar_width, 0.8, 0.1, 2.0)
        self.bar_alpha = self._safe_float(self.bar_alpha, 0.8, 0.0, 1.0)
        self.bar_edge_width = self._safe_float(self.bar_edge_width, 0.5, 0.0, 5.0)
        self.bar_edge_color = self._safe_color(self.bar_edge_color, "black")

        # Y LIMITS
        if self.y_limit_min is not None:
            self.y_limit_min = self._safe_float(self.y_limit_min, None)

        if self.y_limit_max is not None:
            self.y_limit_max = self._safe_float(self.y_limit_max, None)

        if self.y_limit_min is not None and self.y_limit_max is not None:
            if self.y_limit_min > self.y_limit_max:
                self._warn("y_limit_min > y_limit_max, resetting limits")
                self.y_limit_min = None
                self.y_limit_max = None

        # LEGEND
        self.legend_loc = str(self.legend_loc).strip().lower()
        if self.legend_loc not in VALID_LEGEND_LOCS:
            self._warn(f"Invalid legend_loc '{self.legend_loc}', using 'lower center'")
            self.legend_loc = "lower center"

        self.legend_n_col = self._safe_int(self.legend_n_col, 4, 1, 20)
        self.legend_frame = self._safe_bool(self.legend_frame, True)
        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0.0, 1.0)
        self.legend_bbox_y = self._safe_float(self.legend_bbox_y, -0.37, -2.0, 2.0)
        self.legend_margin_bottom = self._safe_float(self.legend_margin_bottom, 0.02, 0.0, 1.0)

        # COLORS
        self.color_hydrophobic = self._safe_color(self.color_hydrophobic, "#b39ddb")
        self.color_hydrogen_bonds = self._safe_color(self.color_hydrogen_bonds, "#f1c40f")
        self.color_pi_cation = self._safe_color(self.color_pi_cation, "#e67e22")
        self.color_pi_stacking = self._safe_color(self.color_pi_stacking, "#2ecc71")
        self.color_water_bridge_1 = self._safe_color(self.color_water_bridge_1, "#3498db")
        self.color_ionic = self._safe_color(self.color_ionic, "#ff3333")
        self.color_halogen = self._safe_color(self.color_halogen, "#f49ac2")
        self.color_metal_contact = self._safe_color(self.color_metal_contact, "#95a5a6")

        self.error_bars_color = self._safe_color(self.error_bars_color, "black")

        # ERROR BARS
        self.error_bars_capsize = self._safe_float(self.error_bars_capsize, 3.0, 0.0, 20.0)
        self.error_bars_alpha = self._safe_float(self.error_bars_alpha, 0.5, 0.0, 1.0)
        self.error_bars_line_width = self._safe_float(self.error_bars_line_width, 0.5, 0.0, 5.0)

        self.error_bars_line_style = str(self.error_bars_line_style).strip().lower()
        if self.error_bars_line_style not in VALID_LINE_STYLES:
            self._warn(
                f"Invalid error_bars_line_style '{self.error_bars_line_style}', using 'solid'"
            )
            self.error_bars_line_style = "solid"

        # RENUMBER
        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        # X AXIS REPRESENTATION
        pattern = re.compile(
            r"^(resname|resid|chainid|segid)"
            r"(([:\-_]?)(resname|resid|chainid|segid))*$"
        )

        rep = str(self.x_axis_representation).strip().lower()

        if not pattern.fullmatch(rep):
            self._warn(
                f"Invalid x_axis_representation '{self.x_axis_representation}', "
                "using 'chainid:resnameresid'"
            )
            self.x_axis_representation = "chainid:resnameresid"
        else:
            self.x_axis_representation = rep

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


@dataclass
class ProteinLigandInteractionsStackedColumn2PlotSettings(PlotSettingsBase):
    alias: ClassVar[Tuple[str, ...]] = ("PLI-STACKED-COLUMN-2",
                                        "PLI_STACKED_COLUMN_2",
                                        "PL-INTERACTIONS-STACKED-COLUMN-2",
                                        "PL-INTERACTIONS_STACKED_COLUMN_2")

    fig_size_width: float = 6
    fig_size_height: float = 4
    fig_dpi: int = 800
    fig_basename: str = "pli_stacked_column_2"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    fig_title: str = "Protein-Ligand Interactions Stacked Column 2"
    x_label: str = "Residues"
    y_label: str = "Frequency"

    font_size_x: int = 8
    font_weight_x: str = "normal"

    font_size_y: int = 8
    font_weight_y: str = "normal"

    font_size_title: int = 14
    font_weight_title: str = "bold"

    font_size_legend: int = 8
    font_weight_legend: str = "normal"

    font_size_label: int = 8
    font_weight_label: str = "normal"

    font_size_ticks: int = 8
    font_weight_ticks: str = "normal"

    font_family: str = "dejavu sans"

    x_tick_rotation: int | float = 90
    y_tick_rotation: int | float = 0

    threshold: float = 0.02

    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False

    disable_title: bool = False
    disable_legend: bool = False
    disable_ticks: bool = False

    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_color: str = "lightgray"
    grid_alpha: float = 0.5

    bar_width: float = 0.8
    bar_edge_color: str = "black"
    bar_alpha: float = 0.8
    bar_edge_width: float = 0.5

    y_limit_min: float | None = None
    y_limit_max: float | None = None

    legend_loc: str = "upper right"
    legend_bbox_y: float = 0
    legend_n_col: int = 2
    legend_frame: bool = True
    legend_alpha: float = 1.0
    legend_margin_bottom: float = 0.02

    color_backbone: str = "#b28dff"
    color_side_chain: str = "#aff8db"

    aa3_to_aa1: bool = True
    renumber: bool = False
    renumber_int: int | None = None
    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    x_axis_representation: str = "chainid:resnameresid"

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))

    def _validate_fields(self) -> None:
        """
        Validates and normalizes the configuration properties of a plotting or visualization
        object to ensure their correctness, completeness, and compliance with predefined
        constraints. Performs tasks such as type adjustments, boundary checks, default
        assignments, and warning notifications when values are invalid or missing.

        :param self: Reference to the object on which the field validation and normalization
                     operations are performed.

        :raises Warning: If specific fields contain invalid or unsupported values, warnings
                         may be raised to notify the user of the default value adjustments.

        :return: None
        """

        # FIGURE BASICS
        self.fig_dpi = self._safe_int(self.fig_dpi, 800, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 6.0, 1.0, 40.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 4.0, 1.0, 40.0)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            ext = ".png"
        self.fig_format = ext.lstrip(".")

        self.fig_basename = str(self.fig_basename).strip() or "pli_stacked_column_2"

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)

        self.bg_color = self._safe_color(self.bg_color, "white")

        # LABELS
        self.fig_title = str(self.fig_title).strip() or \
                         "Protein-Ligand Interactions Stacked Column 2"

        self.x_label = str(self.x_label).strip() or "Residues"
        self.y_label = str(self.y_label).strip() or "Frequency"

        # FONT SIZES
        self.font_size_x = self._safe_int(self.font_size_x, 8, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 8, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 14, 1)
        self.font_size_legend = self._safe_int(self.font_size_legend, 8, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 8, 1)
        self.font_size_ticks = self._safe_int(self.font_size_ticks, 8, 1)

        # FONT WEIGHTS
        def _fw(value: str, default: str) -> str:
            """
            Validates and adjusts the font weight value based on predefined acceptable font weights.

            The method checks the given font weight against a list of valid font weights. If the
            value is invalid, it logs a warning and uses a default font weight instead. The function
            ensures that the returned value is properly formatted by stripping whitespace and
            converting it to lowercase.

            :param value: The font weight to validate.
            :param default: The default font weight to use if the given value is not valid.
            :return: The valid or default font weight.
            """
            v = str(value).strip().lower()
            if v not in VALID_FONT_WEIGHTS:
                self._warn(f"Invalid font weight '{value}', using '{default}'")
                return default
            return v

        self.font_weight_x = _fw(self.font_weight_x, "normal")
        self.font_weight_y = _fw(self.font_weight_y, "normal")
        self.font_weight_title = _fw(self.font_weight_title, "bold")
        self.font_weight_legend = _fw(self.font_weight_legend, "normal")
        self.font_weight_label = _fw(self.font_weight_label, "normal")
        self.font_weight_ticks = _fw(self.font_weight_ticks, "normal")

        # FONT FAMILY
        font = str(self.font_family).strip().lower()
        if font not in AVAILABLE_FONTS:
            self._warn(f"Font '{font}' not found, using 'dejavu sans'")
            self.font_family = "dejavu sans"
        else:
            self.font_family = font

        # ROTATIONS
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # THRESHOLD
        self.threshold = self._safe_float(self.threshold, 0.02, 0.0, 1.0)

        # BOOLEAN FLAGS
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_legend = self._safe_bool(self.disable_legend, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)
        self.enable_grid = self._safe_bool(self.enable_grid, True)

        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, True)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)

        # GRID
        self.grid_style = str(self.grid_style).strip().lower()
        if self.grid_style not in VALID_LINE_STYLES:
            self._warn(f"Invalid grid_style '{self.grid_style}', using 'dashed'")
            self.grid_style = "dashed"

        self.grid_color = self._safe_color(self.grid_color, "lightgray")
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.5, 0.0, 1.0)

        # BARS
        self.bar_width = self._safe_float(self.bar_width, 0.8, 0.1, 2.0)
        self.bar_alpha = self._safe_float(self.bar_alpha, 0.8, 0.0, 1.0)
        self.bar_edge_width = self._safe_float(self.bar_edge_width, 0.5, 0.0, 5.0)
        self.bar_edge_color = self._safe_color(self.bar_edge_color, "black")

        # Y LIMITS
        if self.y_limit_min is not None:
            self.y_limit_min = self._safe_float(self.y_limit_min, None)

        if self.y_limit_max is not None:
            self.y_limit_max = self._safe_float(self.y_limit_max, None)

        if self.y_limit_min is not None and self.y_limit_max is not None:
            if self.y_limit_min > self.y_limit_max:
                self._warn("y_limit_min > y_limit_max, resetting limits")
                self.y_limit_min = None
                self.y_limit_max = None

        # LEGEND
        self.legend_loc = str(self.legend_loc).strip().lower()
        if self.legend_loc not in VALID_LEGEND_LOCS:
            self._warn(f"Invalid legend_loc '{self.legend_loc}', using 'upper right'")
            self.legend_loc = "upper right"

        self.legend_n_col = self._safe_int(self.legend_n_col, 2, 1, 20)
        self.legend_frame = self._safe_bool(self.legend_frame, True)
        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0.0, 1.0)
        self.legend_bbox_y = self._safe_float(self.legend_bbox_y, 0.0, -2.0, 2.0)
        self.legend_margin_bottom = self._safe_float(self.legend_margin_bottom, 0.02, 0.0, 1.0)

        # COLORS
        self.color_backbone = self._safe_color(self.color_backbone, "#b28dff")
        self.color_side_chain = self._safe_color(self.color_side_chain, "#aff8db")

        if self.color_backbone == self.color_side_chain:
            self._warn("Backbone and Sidechain colors identical. Resetting to defaults.")
            self.color_backbone = "#b28dff"
            self.color_side_chain = "#aff8db"

        # RENUMBER
        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        # X AXIS REPRESENTATION
        pattern = re.compile(
            r"^(resname|resid|chainid|segid)"
            r"(([:\-_]?)(resname|resid|chainid|segid))*$"
        )

        rep = str(self.x_axis_representation).strip().lower()

        if not pattern.fullmatch(rep):
            self._warn(
                f"Invalid x_axis_representation '{self.x_axis_representation}', "
                "using 'chainid:resnameresid'"
            )
            self.x_axis_representation = "chainid:resnameresid"
        else:
            self.x_axis_representation = rep


@dataclass
class ProteinLigandInteractionsHeatmap1PlotSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PLI-HEATMAP-1",
        "PLI_HEATMAP_1",
        "PL-INTERACTIONS-HEATMAP-1",
        "PL_INTERACTIONS_HEATMAP_1",
    )

    # Figure
    fig_size_width: float = 10.0
    fig_size_height: float = 6.0
    fig_dpi: int = 600
    fig_basename: str = "pli_heatmap_1"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Labels & Title
    fig_title: str = "Protein-Ligand Interaction Heatmap"
    x_label: str = "Frame"
    y_label: str = "Residue"

    font_size_x: int = 8
    font_weight_x: str = "normal"

    font_size_y: int = 8
    font_weight_y: str = "normal"

    font_size_label: int = 10
    font_weight_label: str = "normal"

    font_size_title: int = 12
    font_weight_title: str = "bold"

    font_size_cbar: int = 8
    font_weight_cbar: str = "normal"

    font_family: str = "dejavu sans"

    # Ticks
    x_tick_rotation: float = 90.0
    y_tick_rotation: float = 0.0

    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False
    disable_title: bool = False
    disable_legend: bool = False
    disable_ticks: bool = False

    # Heatmap behavior
    threshold: float = 0.02
    cmap: str = "viridis"
    vmin: float | None = 0.0
    vmax: float | None = None
    interpolation: str = "nearest"

    cbar_orientation: str = "vertical"
    cbar_shrink: float = 1.0
    cbar_pad: float = 0.04

    # Residue transformation
    aa3_to_aa1: bool = False
    renumber: bool = False
    renumber_int: int | None = None

    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    y_axis_representation: str = "chainid:resnameresid"

    # Validation
    def _validate_fields(self) -> None:
        """
        Validates and sanitizes fields and attributes to ensure they fall within
        expected ranges or values. Wrong or invalid values are corrected or replaced
        with defaults as necessary. This function ensures consistency and proper setup
        of properties such as figure configurations, text settings, axes options,
        heatmap configurations, and residue transforms.

        :return: None
        """

        # FIGURE
        self.fig_dpi = self._safe_int(self.fig_dpi, 600, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 10.0, 1.0, 100.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 6.0, 1.0, 100.0)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using 'png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_basename = str(self.fig_basename).strip() or "pli_heatmap_1"
        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # TEXT
        self.fig_title = str(self.fig_title).strip() or \
                         "Protein-Ligand Interaction Heatmap"

        self.x_label = str(self.x_label).strip() or "Frame"
        self.y_label = str(self.y_label).strip() or "Residue"

        self.font_size_x = self._safe_int(self.font_size_x, 8, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 8, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 10, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 12, 1)
        self.font_size_cbar = self._safe_int(self.font_size_cbar, 8, 1)

        def _fw(value: str, default: str) -> str:
            v = str(value).strip().lower()
            if v not in VALID_FONT_WEIGHTS:
                self._warn(f"Invalid font weight '{value}', using '{default}'")
                return default
            return v

        self.font_weight_x = _fw(self.font_weight_x, "normal")
        self.font_weight_y = _fw(self.font_weight_y, "normal")
        self.font_weight_label = _fw(self.font_weight_label, "normal")
        self.font_weight_title = _fw(self.font_weight_title, "bold")
        self.font_weight_cbar = _fw(self.font_weight_cbar, "normal")

        font = str(self.font_family).strip().lower()
        if font not in AVAILABLE_FONTS:
            self._warn(f"Font '{font}' not found, using 'dejavu sans'")
            self.font_family = "dejavu sans"
        else:
            self.font_family = font

        # TICKS
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_legend = self._safe_bool(self.disable_legend, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

        # HEATMAP
        self.threshold = self._safe_float(self.threshold, 0.0, 0.0, 1.0)

        self.cmap = str(self.cmap).strip() or "viridis"

        if self.vmin is not None:
            self.vmin = self._safe_float(self.vmin, 0.0)

        if self.vmax is not None:
            self.vmax = self._safe_float(self.vmax, None)

        self.interpolation = str(self.interpolation).strip().lower() or "nearest"

        self.cbar_orientation = str(self.cbar_orientation).strip().lower()
        if self.cbar_orientation not in {"vertical", "horizontal"}:
            self._warn("Invalid cbar_orientation, using 'vertical'")
            self.cbar_orientation = "vertical"

        self.cbar_shrink = self._safe_float(self.cbar_shrink, 1.0, 0.1, 5.0)
        self.cbar_pad = self._safe_float(self.cbar_pad, 0.04, 0.0, 1.0)

        # RESIDUE TRANSFORMS
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, False)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)

        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        pattern = re.compile(
            r"^(resname|resid|chainid|segid)"
            r"(([:\-_]?)(resname|resid|chainid|segid))*$"
        )

        rep = str(self.y_axis_representation).strip().lower()

        if not pattern.fullmatch(rep):
            self._warn(
                f"Invalid y_axis_representation '{self.y_axis_representation}', "
                "using 'chainid:resnameresid'"
            )
            self.y_axis_representation = "chainid:resnameresid"
        else:
            self.y_axis_representation = rep

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


@dataclass
class ProteinLigandInteractionsHeatmap2PlotSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PLI-HEATMAP-2",
        "PLI_HEATMAP_2",
        "PL-INTERACTIONS-HEATMAP-2",
        "PL_INTERACTIONS_HEATMAP_2",
    )

    # Figure basics
    fig_size_width: float | None = None
    fig_size_height: float | None = None
    fig_dpi: int = 800
    fig_basename: str = "pli_heatmap_2"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Labels & title
    fig_title: str = "Protein–Ligand Interaction Frequency Heatmap"
    x_label: str = "Frame"
    y_label: str = "Interaction type"

    font_size_x: int = 10
    font_weight_x: str = "normal"

    font_size_y: int = 10
    font_weight_y: str = "normal"

    font_size_label: int = 10
    font_weight_label: str = "normal"

    font_size_title: int = 12
    font_weight_title: str = "bold"

    font_size_cbar: int = 10

    font_family: str = "dejavu sans"

    x_tick_rotation: float = 90.0
    y_tick_rotation: float = 0.0

    # Filtering / normalization
    drop_empty_rows: bool = True
    threshold: float = 0.02
    normalize: str = "none"   # "none", "by_frame", "max1"

    # Tick control
    xtick_max: int = 50
    disable_x_axis: bool = False
    disable_y_axis: bool = False
    disable_x_label: bool = False
    disable_y_label: bool = False
    disable_title: bool = False
    disable_legend: bool = False
    disable_ticks: bool = False

    # Axis limits
    x_limit_min: float | None = None
    x_limit_max: float | None = None

    # Auto scaling
    per_frame_in: float = 0.10
    per_inter_in: float = 0.85

    min_width: float = 24.0
    max_width: float = 42.0
    min_height: float = 18.0
    max_height: float = 100.0

    # Gridlines
    enable_grid: bool = False
    grid_color: str = "black"
    grid_linewidth: float = 0.2

    # Colormap
    cmap: str = "viridis"
    vmin: float | None = None
    vmax: float | None = None
    interpolation: str = "nearest"

    # Colorbar
    cbar_orientation: str = "vertical"
    cbar_shrink: float = 1.0
    cbar_pad: float = 0.04

    # Validation
    def _validate_fields(self) -> None:

        # Figure basics
        self.fig_dpi = self._safe_int(self.fig_dpi, 800, 50, 2000)

        if self.fig_size_width is not None:
            self.fig_size_width = self._safe_float(self.fig_size_width, None, 1.0, 200.0)

        if self.fig_size_height is not None:
            self.fig_size_height = self._safe_float(self.fig_size_height, None, 1.0, 200.0)

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



@dataclass
class ProteinLigandInteractionsPieCharts1PlotSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PLI-PIE-CHARTS-1",
        "PLI_PIE_CHARTS_1",
        "PL-INTERACTIONS-PIE-CHARTS-1",
        "PL_INTERACTIONS_PIE_CHARTS_1",
    )

    # Figure basics
    fig_size: float = 5.0
    fig_dpi: int = 150
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Residue filtering
    top_n: int = 20

    # Collage
    collage: bool = True
    collage_cols: int = 5
    collage_pad: float = 0.5

    # Optional outputs
    make_pdf: bool = False
    make_overall: bool = True

    # Disable flags
    disable_title: bool = False
    disable_labels: bool = False
    disable_autopct: bool = False
    disable_overall_title: bool = False

    # Typography
    font_size_title: int = 28
    font_size_pct: int = 16
    font_weight_title: str = "bold"

    # Colors
    color_backbone: str = "#b28dff"
    color_side_chain: str = "#aff8db"

    # Residue label representation (SC2 style)
    x_axis_representation: str = "chainid:resnameresid"

    # Validation
    def _validate_fields(self) -> None:

        # Figure
        self.fig_size = self._safe_float(self.fig_size, 5.0, 1.0, 20.0)
        self.fig_dpi = self._safe_int(self.fig_dpi, 150, 50, 2000)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using 'png'")
            self.fig_format = "png"
        else:
            self.fig_format = ext.lstrip(".")

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)
        self.bg_color = self._safe_color(self.bg_color, "white")

        # Top-N
        self.top_n = self._safe_int(self.top_n, 20, 0)

        # Collage
        self.collage = self._safe_bool(self.collage, True)
        self.collage_cols = self._safe_int(self.collage_cols, 5, 1, 20)
        self.collage_pad = self._safe_float(self.collage_pad, 0.5, 0.0, 5.0)

        # Optional outputs
        self.make_pdf = self._safe_bool(self.make_pdf, False)
        self.make_overall = self._safe_bool(self.make_overall, True)

        # Disable flags
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_labels = self._safe_bool(self.disable_labels, False)
        self.disable_autopct = self._safe_bool(self.disable_autopct, False)
        self.disable_overall_title = self._safe_bool(self.disable_overall_title, False)

        # Typography
        self.font_size_title = self._safe_int(self.font_size_title, 28, 1)
        self.font_size_pct = self._safe_int(self.font_size_pct, 16, 1)

        fw = str(self.font_weight_title).strip().lower()
        if fw not in VALID_FONT_WEIGHTS:
            self._warn(f"Invalid font_weight_title '{self.font_weight_title}', using 'bold'")
            self.font_weight_title = "bold"
        else:
            self.font_weight_title = fw

        # Colors
        self.color_side_chain = self._safe_color(self.color_side_chain, "#aff8db")
        self.color_backbone = self._safe_color(self.color_backbone, "#b28dff")

        if self.color_side_chain == self.color_backbone:
            self._warn("Sidechain and Backbone colors identical. Resetting to defaults.")
            self.color_side_chain = "#aff8db"
            self.color_backbone = "#b28dff"

        # X axis representation pattern validation
        pattern = re.compile(
            r"^(resname|resid|chainid|segid)"
            r"(([:\-_]?)(resname|resid|chainid|segid))*$"
        )

        rep = str(self.x_axis_representation).strip().lower()

        if not pattern.fullmatch(rep):
            self._warn(
                f"Invalid x_axis_representation '{self.x_axis_representation}', "
                "using 'chainid:resnameresid'"
            )
            self.x_axis_representation = "chainid:resnameresid"
        else:
            self.x_axis_representation = rep

    # Namespace export
    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))



@dataclass
class ProteinLigandInteractionsLigandMonitorSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PLI-LIGAND-MONITOR",
        "PLI_LIGAND_MONITOR",
        "PL-INTERACTIONS-LIGAND-MONITOR",
        "PL_INTERACTIONS_LIGAND_MONITOR",
    )

    # Figure basics
    fig_size_width: float = 10
    fig_size_height: float = 8
    fig_dpi: int = 300
    fig_basename: str = "pli_ligand_monitor"
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    # Labels & Title
    fig_title: str = "Residue ↔ Ligand Atom Contact Frequency"
    x_label: str = "Ligand atom (Group 2)"
    y_label: str = "Protein residue (Group 1)"

    # Fonts
    font_size_x: int = 10
    font_weight_x: str = "normal"

    font_size_y: int = 10
    font_weight_y: str = "normal"

    font_size_label: int = 12
    font_weight_label: str = "bold"

    font_size_title: int = 16
    font_weight_title: str = "bold"

    font_size_cbar: int = 10
    font_weight_cbar: str = "normal"

    font_family: str = "dejavu sans"

    # Tick rotation
    x_tick_rotation: float = 90.0
    y_tick_rotation: float = 0.0

    # Filtering
    threshold: float = 0.02
    drop_empty_rows: bool = True
    drop_empty_cols: bool = True

    # Representation system
    y_axis_representation: str = "chainid:resnameresid"

    aa3_to_aa1: bool = False
    renumber: bool = False
    renumber_int: int | None = None
    alter_chains: bool = False
    alter_segments: bool = False
    alter_chains_str: str | None = None
    alter_segments_str: str | None = None

    # Colormap
    cmap: str = "viridis"
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

    # Validation
    def _validate_fields(self) -> None:

        # Figure basics
        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 10.0, 2.0, 60.0)
        self.fig_size_height = self._safe_float(self.fig_size_height, 8.0, 2.0, 60.0)

        ext = f".{str(self.fig_format).strip().lower()}"
        if ext not in VALID_EXTENSIONS:
            self._warn(f"Invalid fig_format '{self.fig_format}', using '.png'")
            ext = ".png"
        self.fig_format = ext.lstrip(".")

        self.fig_basename = str(self.fig_basename).strip() or "pli_ligand_monitor"

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)

        self.bg_color = self._safe_color(self.bg_color, "white")

        # Fonts
        self.font_size_x = self._safe_int(self.font_size_x, 10, 1)
        self.font_size_y = self._safe_int(self.font_size_y, 10, 1)
        self.font_size_label = self._safe_int(self.font_size_label, 12, 1)
        self.font_size_title = self._safe_int(self.font_size_title, 16, 1)
        self.font_size_cbar = self._safe_int(self.font_size_cbar, 10, 1)

        def _fw(value: str, default: str) -> str:
            v = str(value).strip().lower()
            if v not in VALID_FONT_WEIGHTS:
                self._warn(f"Invalid font weight '{value}', using '{default}'")
                return default
            return v

        self.font_weight_x = _fw(self.font_weight_x, "normal")
        self.font_weight_y = _fw(self.font_weight_y, "normal")
        self.font_weight_label = _fw(self.font_weight_label, "bold")
        self.font_weight_title = _fw(self.font_weight_title, "bold")
        self.font_weight_cbar = _fw(self.font_weight_cbar, "normal")

        # Font family
        font = str(self.font_family).strip().lower()
        if font not in AVAILABLE_FONTS:
            self._warn(f"Font '{font}' not found, using 'dejavu sans'")
            self.font_family = "dejavu sans"
        else:
            self.font_family = font

        # Rotations
        self.x_tick_rotation = self._safe_float(self.x_tick_rotation, 90.0, -360, 360)
        self.y_tick_rotation = self._safe_float(self.y_tick_rotation, 0.0, -360, 360)

        # Filtering
        self.threshold = self._safe_float(self.threshold, 0.0, 0.0, 1.0)
        self.drop_empty_rows = self._safe_bool(self.drop_empty_rows, True)
        self.drop_empty_cols = self._safe_bool(self.drop_empty_cols, True)

        # Representation
        self.aa3_to_aa1 = self._safe_bool(self.aa3_to_aa1, False)
        self.renumber = self._safe_bool(self.renumber, False)
        self.alter_chains = self._safe_bool(self.alter_chains, False)
        self.alter_segments = self._safe_bool(self.alter_segments, False)

        if self.renumber:
            self.renumber_int = self._safe_int(self.renumber_int, 0)

        pattern = re.compile(
            r"^(resname|resid|chainid|segid)"
            r"(([:\-_]?)(resname|resid|chainid|segid))*$"
        )

        rep = str(self.y_axis_representation).strip().lower()
        if not pattern.fullmatch(rep):
            self._warn(
                f"Invalid y_axis_representation '{self.y_axis_representation}', "
                "using 'chainid:resnameresid'"
            )
            self.y_axis_representation = "chainid:resnameresid"
        else:
            self.y_axis_representation = rep

        # Colormap
        self.cmap = str(self.cmap).strip()
        self.vmin = self._safe_float(self.vmin, 0.0)
        self.vmax = self._safe_float(self.vmax, 1.0)
        self.interpolation = str(self.interpolation).strip().lower()

        # Grid
        self.enable_grid = self._safe_bool(self.enable_grid, False)
        self.grid_color = self._safe_color(self.grid_color, "black")
        self.grid_linewidth = self._safe_float(self.grid_linewidth, 0.2, 0.0, 10.0)

        # Colorbar
        self.disable_colorbar = self._safe_bool(self.disable_colorbar, False)
        self.cbar_orientation = str(self.cbar_orientation).strip().lower()
        self.cbar_shrink = self._safe_float(self.cbar_shrink, 1.0, 0.1, 2.0)
        self.cbar_pad = self._safe_float(self.cbar_pad, 0.04, 0.0, 1.0)

        # Disable switches
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.disable_x_axis = self._safe_bool(self.disable_x_axis, False)
        self.disable_y_axis = self._safe_bool(self.disable_y_axis, False)
        self.disable_x_label = self._safe_bool(self.disable_x_label, False)
        self.disable_y_label = self._safe_bool(self.disable_y_label, False)
        self.disable_ticks = self._safe_bool(self.disable_ticks, False)

    # Namespace export
    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))


