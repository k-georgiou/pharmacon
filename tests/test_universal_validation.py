"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for the completed `PlotUniversalSettings` validator
(RMSD / angles / distances time-series).

The validator previously left many fields uncoerced — `fig_format`,
`bg_color`, `line_style`, `grid_style`, `x_axis`, the font sizes/weights,
`grid_alpha`, `legend_alpha` — so an invalid OR empty value crashed deep in
matplotlib at render time. They are now validated/coerced (with a warning),
empty/unset restores the default silently, and the documented `x_axis =
frame` is accepted as an alias for `frame_index`.
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import PlotUniversalSettings as U


class TestInvalidValuesCoerced:
    def test_fig_format(self):
        s = U.from_dict({"fig_format": "xyz"})
        assert s.fig_format == "png"
        assert s._current_warnings >= 1

    def test_bg_color(self):
        s = U.from_dict({"bg_color": "notacolor"})
        assert s.bg_color == "white"
        assert s._current_warnings >= 1

    def test_line_style(self):
        s = U.from_dict({"line_style": "wiggly"})
        assert s.line_style == "solid"
        assert s._current_warnings >= 1

    def test_grid_style(self):
        s = U.from_dict({"grid_style": "wiggly"})
        assert s.grid_style == "dashed"
        assert s._current_warnings >= 1

    def test_font_weight(self):
        s = U.from_dict({"font_weight_title": "superbold"})
        assert s.font_weight_title in ("bold", "normal")
        assert s._current_warnings >= 1

    def test_font_family(self):
        s = U.from_dict({"font_family": "NotARealFont"})
        assert s.font_family == "dejavu sans"
        assert s._current_warnings >= 1


class TestEmptyFieldsSilent:
    def test_empty_font_size(self):
        s = U.from_dict({"font_size_title": "", "font_size_label": ""})
        assert s.font_size_title == 14 and s.font_size_label == 10
        assert s._current_warnings == 0

    def test_empty_alphas(self):
        s = U.from_dict({"grid_alpha": "", "legend_alpha": ""})
        assert s.grid_alpha == 0.4 and s.legend_alpha == 1.0
        assert s._current_warnings == 0

    def test_garbage_alpha_still_warns(self):
        s = U.from_dict({"grid_alpha": "abc"})
        assert s.grid_alpha == 0.4
        assert s._current_warnings >= 1


class TestXAxisValidation:
    @pytest.mark.parametrize("value", ["time_ps", "time_ns", "time_us", "frame_index"])
    def test_valid_values_preserved(self, value):
        s = U.from_dict({"x_axis": value})
        assert s.x_axis == value
        assert s._current_warnings == 0

    def test_frame_alias_maps_to_frame_index(self):
        # The example INI documents `frame`; it must work, not crash.
        s = U.from_dict({"x_axis": "frame"})
        assert s.x_axis == "frame_index"
        assert s._current_warnings == 0

    def test_case_insensitive(self):
        s = U.from_dict({"x_axis": "TIME_NS"})
        assert s.x_axis == "time_ns"

    def test_invalid_coerced_to_time_ns(self):
        s = U.from_dict({"x_axis": "parsecs"})
        assert s.x_axis == "time_ns"
        assert s._current_warnings >= 1


class TestValidValuesPreserved:
    def test_defaults_silent(self):
        s = U.from_dict({})
        assert s._current_warnings == 0
        assert s.fig_format == "png" and s.x_axis == "time_ns"
        assert s.line_style == "solid" and s.grid_style == "dashed"

    def test_valid_overrides_silent(self):
        s = U.from_dict({
            "fig_format": "svg", "bg_color": "black",
            "line_style": "dashed", "grid_style": "dotted",
            "font_weight_title": "bold", "grid_alpha": 0.5, "legend_alpha": 0.8,
        })
        assert s._current_warnings == 0
        assert s.fig_format == "svg" and s.bg_color == "black"
        assert s.line_style == "dashed" and s.grid_style == "dotted"
        assert s.grid_alpha == 0.5 and s.legend_alpha == 0.8
