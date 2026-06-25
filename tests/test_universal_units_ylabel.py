"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for the units-aware y-label enhancement in the universal time-series
plotter.

When the dataset carries a ``units`` attribute (e.g. angles → "degrees") and
the user kept the default y-label ("Value"), the y-axis label becomes
"Value (<units>)". An explicit y-label is always respected, and groups
without a ``units`` attribute keep the plain default.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_universal_pta

from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.plotter.universal import plot_pta_timeseries_from_file
from pharmacon.constants.plots import PlotUniversalSettings


def _ylabel_of(svg_path: Path) -> str:
    import re
    import html
    txt = svg_path.read_text()
    texts = [html.unescape(re.sub("<[^>]+>", "", t))
             for t in re.findall(r"<text[^>]*>(.*?)</text>", txt, re.S)]
    return " ".join(texts)


def _render(tmp_path, group, series, *, units=None, fmt="svg", **overrides):
    pta = build_universal_pta(tmp_path / f"{group}.pta", group=group,
                              series=series, units=units)
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    settings = PlotUniversalSettings.from_dict({"fig_format": fmt,
                                                "fig_dpi": 80, **overrides})
    with PharmaconPTAFile(pta, mode="r") as f:
        plot_pta_timeseries_from_file(f, group_name=group, settings=settings,
                                      out_dir=out, is_merged=False)
    svg = next(p for p in out.iterdir() if p.suffix == ".svg")
    return _ylabel_of(svg)


ANGLE_SERIES = {"dihedral": {"dihedral": [10.0, 20.0, 30.0, 25.0, 15.0]}}
RMSD_SERIES = {"Calpha": [1.0, 1.2, 1.1, 1.3, 1.4]}


def test_units_appended_when_default_label(tmp_path):
    txt = _render(tmp_path, "angles", ANGLE_SERIES, units="degrees")
    assert "Value (degrees)" in txt


def test_no_units_attr_keeps_plain_default(tmp_path):
    # rmsd mock writes no 'units' attr → label stays the plain default.
    txt = _render(tmp_path, "rmsd", RMSD_SERIES, units=None)
    assert "Value" in txt
    assert "Value (" not in txt


def test_explicit_label_is_respected(tmp_path):
    txt = _render(tmp_path, "angles", ANGLE_SERIES, units="degrees",
                  y_label="Torsion (deg)")
    assert "Torsion (deg)" in txt
    assert "Value (degrees)" not in txt


def test_units_present_but_empty_string_no_suffix(tmp_path):
    txt = _render(tmp_path, "angles", ANGLE_SERIES, units="")
    assert "Value" in txt
    assert "Value (" not in txt
