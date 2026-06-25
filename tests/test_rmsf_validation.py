"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for the completed `RMSFPlotSettings` validator.

The validator previously left `fig_format`, `bg_color`, `line_style`,
`grid_style`, `font_family`, the font sizes and the font weights uncoerced,
so an invalid OR empty value crashed deep in matplotlib at render time. They
are now validated/coerced (with a warning) and empty/unset restores the
default silently.

Rendering is exercised at BOTH the residue level (``x_axis=resid``) and the
atom level (``x_axis=atom_index`` / ``atom_name`` / ``position``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_rmsf_pta

from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.plotter.universal import plot_pta_rmsf_from_file
from pharmacon.constants.plots import RMSFPlotSettings as R


# --------------------------------------------------------------------------
# Settings-level: invalid values coerce + warn; empty values are silent.
# --------------------------------------------------------------------------

class TestInvalidValuesCoerced:
    def test_fig_format(self):
        s = R.from_dict({"fig_format": "xyz"})
        assert s.fig_format == "png"
        assert s._current_warnings >= 1

    def test_bg_color(self):
        s = R.from_dict({"bg_color": "notacolor"})
        assert s.bg_color == "white"
        assert s._current_warnings >= 1

    def test_line_style(self):
        s = R.from_dict({"line_style": "wiggly"})
        assert s.line_style == "solid"
        assert s._current_warnings >= 1

    def test_grid_style(self):
        s = R.from_dict({"grid_style": "wiggly"})
        assert s.grid_style == "dashed"
        assert s._current_warnings >= 1

    def test_font_weight(self):
        s = R.from_dict({"font_weight_title": "superbold"})
        assert s.font_weight_title in ("bold", "normal")
        assert s._current_warnings >= 1

    def test_font_family(self):
        s = R.from_dict({"font_family": "NotARealFont"})
        assert s.font_family == "dejavu sans"
        assert s._current_warnings >= 1


class TestEmptyFieldsSilent:
    def test_empty_font_sizes(self):
        s = R.from_dict({"font_size_title": "", "font_size_label": "",
                         "font_size_ticks": "", "font_size_legend": ""})
        assert (s.font_size_title, s.font_size_label,
                s.font_size_ticks, s.font_size_legend) == (14, 10, 8, 8)
        assert s._current_warnings == 0

    def test_empty_alphas_silent(self):
        # grid_alpha / legend_alpha go through _safe_float → empty is silent.
        s = R.from_dict({"grid_alpha": "", "legend_alpha": ""})
        assert s.grid_alpha == 0.4 and s.legend_alpha == 1.0
        assert s._current_warnings == 0


class TestValidValuesPreserved:
    def test_defaults_silent(self):
        s = R.from_dict({})
        assert s._current_warnings == 0
        assert s.fig_format == "png" and s.line_style == "solid"
        assert s.grid_style == "dashed" and s.bg_color == "white"

    def test_valid_overrides_silent(self):
        s = R.from_dict({"fig_format": "svg", "bg_color": "black",
                         "line_style": "dashed", "grid_style": "dotted",
                         "font_weight_title": "bold", "font_size_title": 20})
        assert s._current_warnings == 0
        assert s.fig_format == "svg" and s.bg_color == "black"
        assert s.line_style == "dashed" and s.grid_style == "dotted"
        assert s.font_size_title == 20


# --------------------------------------------------------------------------
# Render-level: previously-crashing settings now render at residue & atom level.
# --------------------------------------------------------------------------

RESIDUE_LEVEL = "resid"
ATOM_LEVELS = ["atom_index", "atom_name", "position"]


def _render(tmp_path, *, x_axis, **overrides):
    pta = build_rmsf_pta(tmp_path / "rmsf.pta")
    out = tmp_path / f"out_{x_axis}"
    out.mkdir(exist_ok=True)
    s = R.from_dict({"fig_dpi": 80, "fig_size_width": 4, "fig_size_height": 3,
                     "x_axis": x_axis, **overrides})
    with PharmaconPTAFile(pta, mode="r") as f:
        plot_pta_rmsf_from_file(f, group_name="rmsf", settings=s,
                                out_dir=out, is_merged=False)
    files = [p for p in out.iterdir() if p.is_file()]
    assert files and all(p.stat().st_size > 0 for p in files)


@pytest.mark.slow
@pytest.mark.parametrize("x_axis", [RESIDUE_LEVEL] + ATOM_LEVELS)
class TestRenderBothLevels:
    def test_baseline_renders(self, tmp_path, x_axis):
        _render(tmp_path, x_axis=x_axis)

    def test_invalid_settings_now_render(self, tmp_path, x_axis):
        # All of these crashed before the validator was completed.
        _render(tmp_path, x_axis=x_axis,
                fig_format="xyz", bg_color="notacolor",
                line_style="wiggly", grid_style="wiggly",
                font_weight_title="superbold")

    def test_empty_font_sizes_now_render(self, tmp_path, x_axis):
        _render(tmp_path, x_axis=x_axis,
                font_size_title="", font_size_label="")
