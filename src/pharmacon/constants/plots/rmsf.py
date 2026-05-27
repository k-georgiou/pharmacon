"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

RMSF plot settings — per-atom Root Mean Square Fluctuation profiles.

Unlike the universal time-series settings (which assume a frame axis), this
class is shaped for the per-selection / per-atom layout produced by the
``trajectory rmsf`` subcommand and merged in by ``merge results``.
"""
from ._base import (
    PlotSettingsBase,
    Namespace,
    ClassVar,
    Tuple,
    List,
    dataclass,
    field,
    asdict,
)


__all__ = ["RMSFPlotSettings"]


_VALID_X_AXES = frozenset({"resid", "atom_index", "position", "atom_name"})
_XTICK_FORMAT_FIELDS = frozenset({"atom_index", "resid", "resname", "atom_name", "position"})


@dataclass
class RMSFPlotSettings(PlotSettingsBase):

    alias: ClassVar[Tuple[str, ...]] = (
        "PTA-RMSF",
        "PTA_RMSF",
    )

    # ---------------- FIGURE ----------------
    fig_size_width: float = 8
    fig_size_height: float = 5
    fig_dpi: int = 300
    fig_format: str = "png"
    fig_transparent: bool = False
    tight_layout: bool = True
    bg_color: str = "white"

    fig_title: str = "RMSF Profile"
    x_label: str = "Residue"
    y_label: str = "RMSF (Å)"
    x_axis: str = "resid"  # resid | atom_index | position | atom_name

    # ---------------- AXIS LIMITS ----------------
    # "auto" leaves matplotlib's autoscaling; numeric overrides it.
    x_min: str = "auto"
    x_max: str = "auto"
    y_min: str = "auto"
    y_max: str = "auto"

    # ---------------- X-TICK LABELS ----------------
    # Custom per-atom tick labels. Empty string = auto:
    #   x_axis=atom_name → "{atom_name}" labels
    #   other modes      → matplotlib's automatic numeric ticks.
    # Accepts a Python format string with any of:
    #   {atom_index} {resid} {resname} {atom_name} {position}
    # Examples:
    #   xtick_format = {atom_name}
    #   xtick_format = {atom_index} {atom_name}
    #   xtick_format = {resname}{resid}-{atom_name}
    xtick_format: str = ""
    xtick_rotation: str = "auto"   # "auto" or a numeric angle (e.g. "45")
    xtick_max_labels: int = 200    # thin ticks down to ~this many when exceeded

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

    # Per-selection color overrides. Format: "label1:#hex, label2:#hex, ...".
    # Selections not listed fall back to the cycling `line_colors`.
    colors_by_label: str = ""

    # ---------------- STD BAND ----------------
    show_std_band: bool = True
    std_band_alpha: float = 0.25

    # ---------------- SHADED REGIONS ----------------
    # Semicolon-separated list of "start,end,color[,alpha[,label]]" entries.
    # The start/end values are interpreted in the SAME units as `x_axis`
    # (resid → residue numbers; atom_index → topology indices; position /
    # atom_name → 0..N-1 dataset positions).
    # Examples:
    #   shading = 0,10,#888888; 11,100,#ffa500
    #   shading = 50,75,#d62728,0.4,active site
    shading: str = ""
    shading_alpha: float = 0.25         # default alpha when a region omits its own
    shading_show_legend: bool = False   # include labeled regions in the legend

    # ---------------- GRID ----------------
    enable_grid: bool = True
    grid_style: str = "dashed"
    grid_alpha: float = 0.4

    # ---------------- LEGEND ----------------
    disable_legend: bool = False
    legend_frame: bool = True
    legend_alpha: float = 1.0

    # ---------------- BEHAVIOR ----------------
    plot_multiple: bool = False
    disable_title: bool = False

    # VALIDATION
    def _validate_fields(self) -> None:

        self.fig_dpi = self._safe_int(self.fig_dpi, 300, 50, 2000)
        self.fig_size_width = self._safe_float(self.fig_size_width, 8, 2, 60)
        self.fig_size_height = self._safe_float(self.fig_size_height, 5, 2, 60)

        self.fig_transparent = self._safe_bool(self.fig_transparent, False)
        self.tight_layout = self._safe_bool(self.tight_layout, True)

        self.line_width = self._safe_float(self.line_width, 1.5, 0.1, 10)
        self.line_alpha = self._safe_float(self.line_alpha, 1.0, 0.0, 1.0)
        self.std_band_alpha = self._safe_float(self.std_band_alpha, 0.25, 0.0, 1.0)
        self.shading_alpha = self._safe_float(self.shading_alpha, 0.25, 0.0, 1.0)
        self.grid_alpha = self._safe_float(self.grid_alpha, 0.4, 0.0, 1.0)
        self.legend_alpha = self._safe_float(self.legend_alpha, 1.0, 0.0, 1.0)

        self.enable_grid = self._safe_bool(self.enable_grid, True)
        self.show_std_band = self._safe_bool(self.show_std_band, True)
        self.disable_legend = self._safe_bool(self.disable_legend, False)
        self.legend_frame = self._safe_bool(self.legend_frame, True)
        self.plot_multiple = self._safe_bool(self.plot_multiple, False)
        self.cycle_colors = self._safe_bool(self.cycle_colors, True)
        self.disable_title = self._safe_bool(self.disable_title, False)
        self.shading_show_legend = self._safe_bool(self.shading_show_legend, False)

        # Parse shading string into a list of dicts on the instance.
        self.shading_regions = self._parse_shading(
            str(self.shading), self.shading_alpha,
        )

        # x_axis: enum-like
        x_axis_val = str(self.x_axis).strip().lower()
        if x_axis_val not in _VALID_X_AXES:
            self._warn(
                f"Invalid x_axis '{self.x_axis}', "
                f"using default 'resid' (valid: {sorted(_VALID_X_AXES)})"
            )
            x_axis_val = "resid"
        self.x_axis = x_axis_val

        # xtick_format: empty string = auto. Otherwise reject templates that
        # reference unknown placeholder fields (catch typos early).
        fmt = str(self.xtick_format)
        if fmt.strip():
            import string
            try:
                referenced = {
                    name for _, name, _, _ in string.Formatter().parse(fmt)
                    if name
                }
            except ValueError:
                self._warn(f"Malformed xtick_format '{fmt}', disabling custom labels.")
                fmt = ""
                referenced = set()
            unknown = referenced - _XTICK_FORMAT_FIELDS
            if unknown:
                self._warn(
                    f"xtick_format references unknown field(s) {sorted(unknown)}; "
                    f"valid fields: {sorted(_XTICK_FORMAT_FIELDS)}. Disabling custom labels."
                )
                fmt = ""
        self.xtick_format = fmt

        # xtick_rotation: "auto" or a numeric string in [-360, 360].
        rot = str(self.xtick_rotation).strip().lower()
        if rot == "auto":
            self.xtick_rotation = "auto"
        else:
            try:
                angle = float(rot)
                if -360 <= angle <= 360:
                    self.xtick_rotation = rot
                else:
                    self._warn(
                        f"xtick_rotation '{self.xtick_rotation}' out of [-360, 360], "
                        f"using 'auto'."
                    )
                    self.xtick_rotation = "auto"
            except ValueError:
                self._warn(
                    f"Invalid xtick_rotation '{self.xtick_rotation}', using 'auto'."
                )
                self.xtick_rotation = "auto"

        self.xtick_max_labels = self._safe_int(self.xtick_max_labels, 200, 5, 100000)

        # Accept comma-separated strings (INI loader delivers them this way
        # for list-typed fields). Already-listy inputs pass through.
        if isinstance(self.line_colors, str):
            self.line_colors = [c.strip() for c in self.line_colors.split(",") if c.strip()]

        # Validate colors
        validated = []
        for c in self.line_colors:
            validated.append(self._safe_color(c, "#1f77b4"))

        if not validated:
            self._warn("line_colors empty, restoring defaults")
            validated = ["#1f77b4"]

        self.line_colors = validated

        # Per-label color overrides: "label:hex, label:hex"
        self.colors_by_label_map = self._parse_colors_by_label(str(self.colors_by_label))

        # Axis limits: "auto" or numeric
        self.x_min = self._parse_axis_limit(self.x_min, "x_min")
        self.x_max = self._parse_axis_limit(self.x_max, "x_max")
        self.y_min = self._parse_axis_limit(self.y_min, "y_min")
        self.y_max = self._parse_axis_limit(self.y_max, "y_max")

    def _parse_axis_limit(self, value, field_name: str):
        """
        Accept "auto" (case-insensitive) or a numeric string/value. Returns
        either the literal string "auto" or a float. Native int / float / bool
        inputs are coerced numerically so an INI value of ``0`` (even if the
        loader hands it through as False / True) still resolves correctly.
        """
        if value is None:
            return "auto"
        if isinstance(value, bool):
            return float(value)
        if isinstance(value, (int, float)):
            return float(value)
        s = str(value).strip().lower()
        if s in {"", "auto", "none"}:
            return "auto"
        try:
            return float(s)
        except ValueError:
            self._warn(f"Invalid {field_name} '{value}', using 'auto'.")
            return "auto"

    def _parse_colors_by_label(self, raw: str) -> dict:
        """
        Parse a "label:hex, label:hex" string into a {label: color} dict.
        Skips malformed entries with a warning.
        """
        out: dict = {}
        s = (raw or "").strip()
        if not s:
            return out
        for entry in s.split(","):
            entry = entry.strip()
            if not entry:
                continue
            if ":" not in entry:
                self._warn(
                    f"colors_by_label entry {entry!r} missing ':' separator; skipping."
                )
                continue
            label, color = entry.split(":", 1)
            label = label.strip()
            color = self._safe_color(color.strip(), "#1f77b4")
            if not label:
                self._warn("colors_by_label entry has empty label; skipping.")
                continue
            out[label] = color
        return out

    def _parse_shading(self, raw: str, default_alpha: float) -> list:
        """
        Parse the ``shading`` field into a list of region dicts.

        Each row of `raw` is "start,end,color[,alpha[,label]]" separated by
        semicolons. Bad rows are warned about and skipped — they do not
        block parsing of the rest.
        """
        regions: list = []
        s = (raw or "").strip()
        if not s:
            return regions

        rows = [r.strip() for r in s.split(";") if r.strip()]
        if len(rows) > 100:
            self._warn(
                f"shading has {len(rows)} regions; only the first 100 are used."
            )
            rows = rows[:100]

        for i, row in enumerate(rows, start=1):
            parts = [p.strip() for p in row.split(",")]
            if len(parts) < 3:
                self._warn(
                    f"shading entry #{i} needs 'start,end,color' "
                    f"(got: {row!r}); skipping."
                )
                continue
            try:
                start = float(parts[0])
                end = float(parts[1])
            except ValueError:
                self._warn(
                    f"shading entry #{i}: start/end must be numeric "
                    f"(got: {row!r}); skipping."
                )
                continue
            if end < start:
                start, end = end, start
            color = self._safe_color(parts[2], "#cccccc")

            alpha = default_alpha
            if len(parts) >= 4 and parts[3]:
                try:
                    a = float(parts[3])
                    if 0.0 <= a <= 1.0:
                        alpha = a
                    else:
                        self._warn(
                            f"shading entry #{i}: alpha {a} out of [0,1], "
                            f"using default {default_alpha}."
                        )
                except ValueError:
                    self._warn(
                        f"shading entry #{i}: alpha {parts[3]!r} not numeric, "
                        f"using default {default_alpha}."
                    )

            label = parts[4] if len(parts) >= 5 and parts[4] else None

            regions.append({
                "start": start, "end": end,
                "color": color, "alpha": alpha,
                "label": label,
            })

        return regions

    def namespace(self) -> Namespace:
        return Namespace(**asdict(self))
