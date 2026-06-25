"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for "unset" override handling in the plot-settings
coercion layer (`PlotSettingsBase._safe_*`).

An empty/whitespace string or ``None`` is the documented "leave empty for
auto" / "use the default" sentinel — identical to omitting the key. It must
restore the default WITHOUT emitting a coercion warning, otherwise the
shipped example INIs (which use `vmax =`, `y_limit_min =`, ...) get their
plots skipped under the default `--maxwarnings 0`. A genuinely invalid
non-empty value (e.g. "abc") must still warn.
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import (
    PlotSettingsBase,
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
)


class _Probe(PlotSettingsBase):
    """Minimal harness to exercise the _safe_* helpers directly."""

    def __init__(self):
        self._current_warnings = 0


@pytest.fixture
def probe():
    return _Probe()


# --------------------------------------------------------------------------
# Unit-level: the _safe_* helpers treat "" / whitespace / None as unset.
# --------------------------------------------------------------------------

class TestIsUnset:
    @pytest.mark.parametrize("value", ["", "   ", "\t", None])
    def test_unset_values(self, value):
        assert PlotSettingsBase._is_unset(value) is True

    @pytest.mark.parametrize("value", ["abc", "0", 0, 0.0, False, "none"])
    def test_set_values(self, value):
        # 0 / False / "0" are real inputs, NOT unset.
        assert PlotSettingsBase._is_unset(value) is False


class TestSafeFloatUnset:
    @pytest.mark.parametrize("value", ["", "  ", None])
    def test_empty_returns_default_silently(self, probe, value):
        assert probe._safe_float(value, None) is None
        assert probe._current_warnings == 0

    def test_empty_with_numeric_default_silent(self, probe):
        assert probe._safe_float("", 1.0) == 1.0
        assert probe._current_warnings == 0

    def test_invalid_still_warns(self, probe):
        assert probe._safe_float("abc", None) is None
        assert probe._current_warnings == 1

    def test_valid_value_preserved(self, probe):
        assert probe._safe_float("0.5", None) == 0.5
        assert probe._current_warnings == 0


class TestSafeIntUnset:
    @pytest.mark.parametrize("value", ["", "  ", None])
    def test_empty_returns_default_silently(self, probe, value):
        assert probe._safe_int(value, 7) == 7
        assert probe._current_warnings == 0

    def test_invalid_still_warns(self, probe):
        assert probe._safe_int("abc", 7) == 7
        assert probe._current_warnings == 1


class TestSafeBoolColorUnset:
    def test_bool_empty_silent(self, probe):
        assert probe._safe_bool("", True) is True
        assert probe._current_warnings == 0

    def test_bool_false_is_not_unset(self, probe):
        assert probe._safe_bool("false", True) is False
        assert probe._current_warnings == 0

    def test_color_empty_silent(self, probe):
        assert probe._safe_color("", "white") == "white"
        assert probe._current_warnings == 0

    def test_color_invalid_warns(self, probe):
        assert probe._safe_color("notacolor", "white") == "white"
        assert probe._current_warnings == 1


# --------------------------------------------------------------------------
# Integration: empty overrides on real settings classes raise NO warnings,
# garbage does. Mirrors the shipped-example INI fields.
# --------------------------------------------------------------------------

class TestEmptyOverridesNoWarnings:
    def test_stacked_column_1_empty_ylimits_silent(self):
        s = ProteinLigandInteractionsStackedColumn1PlotSettings.from_dict({
            "y_limit_min": "",
            "y_limit_max": "",
        })
        assert s._current_warnings == 0

    def test_heatmap_1_empty_vmax_silent(self):
        s = ProteinLigandInteractionsHeatmap1PlotSettings.from_dict({"vmax": ""})
        assert s._current_warnings == 0
        assert s.vmax is None

    def test_heatmap_2_empty_numeric_fields_silent(self):
        # These are exactly the fields the example pli_heatmap_2.ini leaves blank.
        s = ProteinLigandInteractionsHeatmap2PlotSettings.from_dict({
            "fig_size_width": "",
            "fig_size_height": "",
            "vmin": "",
            "vmax": "",
            "x_limit_min": "",
            "x_limit_max": "",
        })
        assert s._current_warnings == 0
        assert s.vmin is None and s.vmax is None

    def test_heatmap_2_garbage_vmin_warns(self):
        s = ProteinLigandInteractionsHeatmap2PlotSettings.from_dict({"vmin": "abc"})
        assert s._current_warnings >= 1
        assert s.vmin is None

    def test_heatmap_2_normalize_invalid_warns(self):
        s = ProteinLigandInteractionsHeatmap2PlotSettings.from_dict({"normalize": "bogus"})
        assert s._current_warnings >= 1
        assert s.normalize == "none"

    @pytest.mark.parametrize("mode", ["none", "by_frame", "max1"])
    def test_heatmap_2_normalize_valid_preserved(self, mode):
        s = ProteinLigandInteractionsHeatmap2PlotSettings.from_dict({"normalize": mode})
        assert s.normalize == mode
        assert s._current_warnings == 0
