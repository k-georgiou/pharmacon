"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for the Protein–Protein Interaction (PPI) settings
validation layer.

Previously none of the three PPI settings classes overrode
``_validate_fields``: every INI value reached the plotter raw, so invalid
enums/colors/numbers crashed with a hard matplotlib ValueError at render
time, and the shipped example INIs crashed outright on their empty "auto"
fields (`vmax =`, `bar_edge_width =`). These tests pin down that values are
now validated/coerced (with a warning) and that empty/unset fields restore
the default silently.
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import (
    PlotSettingsBase,
    ProteinProteinInteractionsTimelinePairs,
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsStackedColumnSettings,
)

TIMELINE = ProteinProteinInteractionsTimelinePairs
HEATMAP = ProteinProteinInteractionsHeatmap
STACKED = ProteinProteinInteractionsStackedColumnSettings

ALL_PPI = [TIMELINE, HEATMAP, STACKED]
HEATMAP_LIKE = [TIMELINE, HEATMAP]  # have cmap / interpolation / cbar


@pytest.mark.parametrize("cls", ALL_PPI)
def test_all_ppi_classes_validate(cls):
    """Every PPI class must now own a real validator (not the base no-op)."""
    assert "_validate_fields" in cls.__dict__


# --------------------------------------------------------------------------
# Invalid values are coerced + warn (no crash), not passed through raw.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cls", HEATMAP_LIKE)
class TestColormapValidation:
    def test_invalid_cmap_coerced(self, cls):
        s = cls.from_dict({"cmap": "not_a_cmap"})
        assert s.cmap in ("viridis", "gnuplot")  # each class's default
        assert s._current_warnings >= 1

    def test_valid_cmap_preserved(self, cls):
        s = cls.from_dict({"cmap": "magma"})
        assert s.cmap == "magma"

    def test_invalid_interpolation_coerced(self, cls):
        s = cls.from_dict({"interpolation": "bogus"})
        assert s.interpolation == "nearest"
        assert s._current_warnings >= 1

    def test_invalid_cbar_orientation_coerced(self, cls):
        s = cls.from_dict({"cbar_orientation": "sideways"})
        assert s.cbar_orientation == "vertical"
        assert s._current_warnings >= 1


@pytest.mark.parametrize("cls", ALL_PPI)
class TestSharedValidation:
    def test_invalid_bg_color_coerced(self, cls):
        s = cls.from_dict({"bg_color": "notacolor"})
        assert s.bg_color == "white"
        assert s._current_warnings >= 1

    def test_invalid_fig_format_coerced(self, cls):
        s = cls.from_dict({"fig_format": "xyz"})
        assert s.fig_format == "png"
        assert s._current_warnings >= 1

    def test_invalid_font_weight_coerced(self, cls):
        s = cls.from_dict({"font_weight_title": "superbold"})
        assert s.font_weight_title in ("bold", "normal")
        assert s._current_warnings >= 1

    def test_invalid_font_family_coerced(self, cls):
        s = cls.from_dict({"font_family": "NotARealFont"})
        assert s.font_family == "dejavu sans"
        assert s._current_warnings >= 1

    def test_invalid_representation_coerced(self, cls):
        s = cls.from_dict({"representation": "banana:wrongtoken"})
        assert s.representation == "chainid:resnameresid"
        assert s._current_warnings >= 1

    def test_valid_representation_preserved(self, cls):
        s = cls.from_dict({"representation": "resname-resid"})
        assert s.representation == "resname-resid"
        assert s._current_warnings == 0


class TestStackedColumnEnums:
    def test_invalid_grid_style_coerced(self):
        s = STACKED.from_dict({"grid_style": "wiggly"})
        assert s.grid_style == "--"
        assert s._current_warnings >= 1

    def test_invalid_legend_loc_coerced(self):
        s = STACKED.from_dict({"legend_loc": "nowhere"})
        assert s.legend_loc == "upper center"
        assert s._current_warnings >= 1

    def test_invalid_interaction_color_coerced(self):
        s = STACKED.from_dict({"color_hydrophobic": "notacolor"})
        assert s.color_hydrophobic == "#b39ddb"
        assert s._current_warnings >= 1

    def test_invalid_error_bars_line_style_coerced(self):
        s = STACKED.from_dict({"error_bars_line_style": "zigzag"})
        assert s.error_bars_line_style == "solid"
        assert s._current_warnings >= 1


# --------------------------------------------------------------------------
# Empty "auto" fields (the ones the example INIs leave blank) restore the
# default SILENTLY — this is what made the shipped examples crash before.
# --------------------------------------------------------------------------

class TestEmptyAutoFieldsSilent:
    def test_heatmap_empty_vmax_silent(self):
        # ppi_heatmap.ini: `vmax =`  → crashed before; now None (auto), no warn.
        s = HEATMAP.from_dict({"vmax": ""})
        assert s.vmax is None
        assert s._current_warnings == 0

    def test_stacked_empty_bar_edge_width_silent(self):
        # ppi_stacked_column.ini: `bar_edge_width =` → crashed before.
        s = STACKED.from_dict({"bar_edge_width": ""})
        assert s.bar_edge_width is None
        assert s._current_warnings == 0

    def test_stacked_empty_ylimits_silent(self):
        s = STACKED.from_dict({"y_limit_min": "", "y_limit_max": ""})
        assert s.y_limit_min is None and s.y_limit_max is None
        assert s._current_warnings == 0

    def test_garbage_numeric_still_warns(self):
        s = HEATMAP.from_dict({"vmax": "abc"})
        assert s.vmax is None
        assert s._current_warnings >= 1


# --------------------------------------------------------------------------
# Valid values from the shipped examples propagate unchanged & silently.
# --------------------------------------------------------------------------

class TestValidExampleValues:
    def test_timeline_defaults_silent(self):
        s = TIMELINE.from_dict({"cmap": "gnuplot", "vmin": 0.0, "vmax": 1.0,
                                "top_n": 50, "aa3_to_aa1": False})
        assert s._current_warnings == 0
        assert s.cmap == "gnuplot" and s.vmax == 1.0 and s.top_n == 50

    def test_heatmap_symmetric_and_topn(self):
        s = HEATMAP.from_dict({"symmetric": False, "top_n": 10, "min_total": 0.5})
        assert s.symmetric is False and s.top_n == 10 and s.min_total == 0.5
        assert s._current_warnings == 0

    def test_stacked_aa3_and_colors(self):
        s = STACKED.from_dict({"aa3_to_aa1": True, "color_hydrophobic": "#b39ddb"})
        assert s.aa3_to_aa1 is True and s.color_hydrophobic == "#b39ddb"
        assert s._current_warnings == 0
