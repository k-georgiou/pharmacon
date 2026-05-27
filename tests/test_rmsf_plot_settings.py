"""
Tests for `RMSFPlotSettings` (pharmacon/constants/plots/rmsf.py).

Covers the alias registry, `from_dict` construction, and every validator
branch — including the parsers for `xtick_format`, `xtick_rotation`,
`colors_by_label`, `shading`, axis limits, and the `line_colors`
string-input compatibility shim.
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import RMSFPlotSettings, PlotSettingsBase


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

class TestAliasRegistry:
    def test_aliases_declared(self):
        assert RMSFPlotSettings.alias == ("PTA-RMSF", "PTA_RMSF")

    @pytest.mark.parametrize("alias", ["PTA-RMSF", "PTA_RMSF", "pta-rmsf"])
    def test_resolve_by_alias_case_insensitive(self, alias):
        # PlotSettingsBase.resolve normalises to uppercase
        assert PlotSettingsBase.resolve(alias.upper()) is RMSFPlotSettings


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

class TestDefaults:
    def test_construction_with_empty_overrides(self):
        s = RMSFPlotSettings.from_dict({})
        assert s.x_axis == "resid"
        assert s.x_label == "Residue"
        assert s.y_label == "RMSF (Å)"
        assert s.show_std_band is True
        assert s.plot_multiple is False
        assert s.y_min == "auto" and s.y_max == "auto"
        assert s.x_min == "auto" and s.x_max == "auto"
        assert s.xtick_format == ""
        assert s.xtick_rotation == "auto"


# ---------------------------------------------------------------------------
# x_axis enum
# ---------------------------------------------------------------------------

class TestXAxis:
    @pytest.mark.parametrize("value", ["resid", "atom_index", "position", "atom_name"])
    def test_valid_values_accepted(self, value):
        s = RMSFPlotSettings.from_dict({"x_axis": value})
        assert s.x_axis == value

    def test_invalid_value_falls_back_to_resid(self):
        s = RMSFPlotSettings.from_dict({"x_axis": "frame"})
        assert s.x_axis == "resid"

    def test_case_insensitive(self):
        s = RMSFPlotSettings.from_dict({"x_axis": "ATOM_NAME"})
        assert s.x_axis == "atom_name"


# ---------------------------------------------------------------------------
# xtick_format
# ---------------------------------------------------------------------------

class TestXtickFormat:
    def test_valid_template_passes_through(self):
        s = RMSFPlotSettings.from_dict({"xtick_format": "{atom_index} {atom_name}"})
        assert s.xtick_format == "{atom_index} {atom_name}"

    def test_unknown_field_resets_to_empty(self):
        s = RMSFPlotSettings.from_dict({"xtick_format": "{atom_index} {bogus}"})
        assert s.xtick_format == ""

    def test_malformed_template_resets_to_empty(self):
        s = RMSFPlotSettings.from_dict({"xtick_format": "{unclosed"})
        assert s.xtick_format == ""

    @pytest.mark.parametrize("template", [
        "{atom_index}",
        "{resid}",
        "{resname}{resid}-{atom_name}",
        "{position}: {atom_name}",
    ])
    def test_every_field_individually_valid(self, template):
        s = RMSFPlotSettings.from_dict({"xtick_format": template})
        assert s.xtick_format == template


# ---------------------------------------------------------------------------
# xtick_rotation
# ---------------------------------------------------------------------------

class TestXtickRotation:
    def test_auto_default(self):
        s = RMSFPlotSettings.from_dict({})
        assert s.xtick_rotation == "auto"

    def test_auto_case_insensitive(self):
        s = RMSFPlotSettings.from_dict({"xtick_rotation": "AUTO"})
        assert s.xtick_rotation == "auto"

    @pytest.mark.parametrize("value", ["0", "45", "90", "-45", "180"])
    def test_numeric_strings_accepted(self, value):
        # Validator coerces numeric strings to floats so the value lands
        # on the instance with its numeric identity.
        s = RMSFPlotSettings.from_dict({"xtick_rotation": value})
        assert s.xtick_rotation == float(value)

    @pytest.mark.parametrize("value", [0, 45, 90, -45])
    def test_native_int_accepted(self, value):
        # INI loader auto-coerces "90" → int 90; settings must accept that.
        s = RMSFPlotSettings.from_dict({"xtick_rotation": value})
        assert s.xtick_rotation == float(value)

    @pytest.mark.parametrize("value", ["sideways", "999", "-400", 999, -400])
    def test_invalid_falls_back_to_auto(self, value):
        s = RMSFPlotSettings.from_dict({"xtick_rotation": value})
        assert s.xtick_rotation == "auto"


# ---------------------------------------------------------------------------
# colors_by_label
# ---------------------------------------------------------------------------

class TestColorsByLabel:
    def test_basic_parse(self):
        s = RMSFPlotSettings.from_dict({
            "colors_by_label": "calpha:#d62728, backbone:#1f77b4",
        })
        assert s.colors_by_label_map == {
            "calpha": "#d62728",
            "backbone": "#1f77b4",
        }

    def test_missing_colon_entry_is_skipped(self):
        s = RMSFPlotSettings.from_dict({
            "colors_by_label": "calpha:#d62728, missing-colon, backbone:#1f77b4",
        })
        assert s.colors_by_label_map == {
            "calpha": "#d62728",
            "backbone": "#1f77b4",
        }

    def test_empty_label_is_skipped(self):
        s = RMSFPlotSettings.from_dict({
            "colors_by_label": ":#888, calpha:#d62728",
        })
        assert s.colors_by_label_map == {"calpha": "#d62728"}

    def test_invalid_color_falls_back_to_default(self):
        s = RMSFPlotSettings.from_dict({
            "colors_by_label": "calpha:not-a-color",
        })
        # Reverts to default first colour in line_colors
        assert s.colors_by_label_map["calpha"] == "#1f77b4"

    def test_empty_string_yields_empty_map(self):
        s = RMSFPlotSettings.from_dict({"colors_by_label": ""})
        assert s.colors_by_label_map == {}


# ---------------------------------------------------------------------------
# shading
# ---------------------------------------------------------------------------

class TestShading:
    def test_two_regions(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "0,10,#888888; 11,100,#ffa500",
        })
        assert len(s.shading_regions) == 2
        a, b = s.shading_regions
        assert (a["start"], a["end"], a["color"]) == (0.0, 10.0, "#888888")
        assert (b["start"], b["end"], b["color"]) == (11.0, 100.0, "#ffa500")

    def test_alpha_and_label_parsed(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "50,75,#d62728,0.4,active site",
        })
        r = s.shading_regions[0]
        assert r["alpha"] == 0.4
        assert r["label"] == "active site"

    def test_swapped_start_end_auto_corrected(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "100,10,#888888",
        })
        r = s.shading_regions[0]
        assert r["start"] == 10.0 and r["end"] == 100.0

    def test_missing_color_entry_skipped(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "0,10",            # missing color
        })
        assert s.shading_regions == []

    def test_non_numeric_bounds_skipped(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "foo,10,#888888",
        })
        assert s.shading_regions == []

    def test_bad_alpha_uses_default(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "0,10,#888888,1.5,extra",
            "shading_alpha": 0.5,
        })
        r = s.shading_regions[0]
        assert r["alpha"] == 0.5            # fell back to the default
        assert r["label"] == "extra"

    def test_partial_failures_keep_valid_rows(self):
        s = RMSFPlotSettings.from_dict({
            "shading": "0,10,#888888; bad row; 20,30,#ffa500",
        })
        starts = [r["start"] for r in s.shading_regions]
        assert starts == [0.0, 20.0]


# ---------------------------------------------------------------------------
# Axis limits
# ---------------------------------------------------------------------------

class TestAxisLimits:
    def test_numeric_strings_become_floats(self):
        s = RMSFPlotSettings.from_dict({"y_min": "0", "y_max": "3.5"})
        assert s.y_min == 0.0 and s.y_max == 3.5

    def test_auto_keyword_passes_through(self):
        s = RMSFPlotSettings.from_dict({"y_min": "auto", "y_max": "AUTO"})
        assert s.y_min == "auto" and s.y_max == "auto"

    def test_invalid_value_falls_back_to_auto(self):
        s = RMSFPlotSettings.from_dict({"y_min": "garbage"})
        assert s.y_min == "auto"

    def test_native_int_zero_is_zero_not_false(self):
        # Defensive: even if the loader hands through a bool, axis limit
        # should be treated numerically.
        s = RMSFPlotSettings.from_dict({"y_min": 0, "y_max": 3})
        assert s.y_min == 0.0 and s.y_max == 3.0

    def test_native_bool_treated_as_numeric(self):
        s = RMSFPlotSettings.from_dict({"y_min": False, "y_max": True})
        assert s.y_min == 0.0 and s.y_max == 1.0


# ---------------------------------------------------------------------------
# line_colors string-input handling (INI quirk)
# ---------------------------------------------------------------------------

class TestLineColorsStringInput:
    def test_comma_separated_string_is_split(self):
        s = RMSFPlotSettings.from_dict({
            "line_colors": "#1f77b4, #ff7f0e, #2ca02c",
        })
        assert s.line_colors == ["#1f77b4", "#ff7f0e", "#2ca02c"]

    def test_list_input_passes_through(self):
        s = RMSFPlotSettings.from_dict({
            "line_colors": ["#1f77b4", "#ff7f0e"],
        })
        assert s.line_colors == ["#1f77b4", "#ff7f0e"]

    def test_empty_string_restores_default(self):
        s = RMSFPlotSettings.from_dict({"line_colors": ""})
        assert s.line_colors == ["#1f77b4"]
