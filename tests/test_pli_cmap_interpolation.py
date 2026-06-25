"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for `cmap` / `interpolation` coercion in the PLI heatmap-style
settings classes (heatmap_1, heatmap_2, ligand_monitor).

Previously an invalid `cmap` or `interpolation` passed settings validation
unchanged and only blew up with a hard `ValueError` deep inside matplotlib at
render time. They are now validated and coerced to safe defaults with a warning,
consistent with `fig_format` / `legend_loc` / `grid_style`.
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import (
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
)

HEATMAP_CLASSES = [
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
]


@pytest.mark.parametrize("cls", HEATMAP_CLASSES)
class TestCmapCoercion:
    def test_valid_cmap_preserved(self, cls):
        s = cls.from_dict({"cmap": "plasma"})
        assert s.cmap == "plasma"

    def test_invalid_cmap_falls_back_to_viridis(self, cls):
        s = cls.from_dict({"cmap": "not_a_real_cmap"})
        assert s.cmap == "viridis"
        assert s._current_warnings >= 1


@pytest.mark.parametrize("cls", HEATMAP_CLASSES)
class TestInterpolationCoercion:
    def test_valid_interpolation_preserved(self, cls):
        s = cls.from_dict({"interpolation": "bilinear"})
        assert s.interpolation == "bilinear"

    def test_invalid_interpolation_falls_back_to_nearest(self, cls):
        s = cls.from_dict({"interpolation": "bogus_interp"})
        assert s.interpolation == "nearest"
        assert s._current_warnings >= 1

    def test_interpolation_case_insensitive(self, cls):
        s = cls.from_dict({"interpolation": "BICUBIC"})
        assert s.interpolation == "bicubic"
