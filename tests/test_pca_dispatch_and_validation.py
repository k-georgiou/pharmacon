"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for PCA plotting:

1. Dispatch — the analyzer writes one group per selection named
   ``pca_<selection>`` (e.g. ``pca_protein``), but the plot dispatcher used to
   match the exact string ``"pca"`` only, so real PCA files rendered NOTHING.
   It now matches ``pca`` and any ``pca_*`` group.

2. Settings validation — `fig_format` (scatter / variance_ratio / probability /
   fes) and `cmap` (scatter / probability / fes) were unvalidated and crashed
   matplotlib at render time. They now coerce with a warning.
"""

from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_pca_pta

from pharmacon.command_line.plot.pta import run, _is_pca_group
from pharmacon.constants.plots import (
    PCAPlotScatterSettings, PCAPlotVarianceRatioSettings,
    PCAPlotProbabilityHeatmapSettings, PCAPlotFESHeatmapSettings,
)


# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------

class TestIsPcaGroup:
    @pytest.mark.parametrize("name", ["pca", "pca_protein", "pca_backbone",
                                      "pca_ligand", "pca_CA"])
    def test_matches(self, name):
        assert _is_pca_group(name) is True

    @pytest.mark.parametrize("name", ["pcamol", "rmsf", "pp_interactions",
                                      "pdca", "PCA"])
    def test_non_matches(self, name):
        # exact "pca" or "pca_" prefix only (case-sensitive group keys).
        assert _is_pca_group(name) is False


def _make_args(input_pta: Path, out_dir: Path) -> Namespace:
    return Namespace(
        input=str(input_pta),
        output=out_dir,
        is_merged=False,
        maxwarnings=9999,
        config=None,
        config_overrides={},
        overwrite=True,
        log=str(out_dir / "plot.log"),
        terminal_logging_level="INFO",
        file_logging_level="TRACE",
    )


@pytest.mark.slow
def test_real_pca_group_name_renders_all_five(tmp_path):
    # The production group naming convention — the case the suite was missing.
    pta = build_pca_pta(tmp_path / "pca.pta", n_frames=200, n_components=3,
                        group="pca_protein")
    out = tmp_path / "out"
    out.mkdir()
    run(_make_args(pta, out))

    pngs = list(out.glob("*.png"))
    names = " ".join(p.name for p in pngs)
    assert len(pngs) >= 5, f"expected >=5 PCA plots, got {len(pngs)}: {names}"
    # the five PCA plot kinds should each have produced output
    for kind in ("timeseries", "scatter", "variance", "probability", "fes"):
        assert kind in names, f"missing PCA plot '{kind}' in {names}"


# --------------------------------------------------------------------------
# Settings validation
# --------------------------------------------------------------------------

class TestFigFormatValidation:
    @pytest.mark.parametrize("cls", [
        PCAPlotScatterSettings, PCAPlotVarianceRatioSettings,
        PCAPlotProbabilityHeatmapSettings, PCAPlotFESHeatmapSettings,
    ])
    def test_invalid_fig_format_coerced(self, cls):
        s = cls.from_dict({"fig_format": "xyz"})
        assert s.fig_format == "png"
        assert s._current_warnings >= 1

    @pytest.mark.parametrize("cls", [
        PCAPlotScatterSettings, PCAPlotVarianceRatioSettings,
        PCAPlotProbabilityHeatmapSettings, PCAPlotFESHeatmapSettings,
    ])
    def test_valid_fig_format_preserved(self, cls):
        s = cls.from_dict({"fig_format": "svg"})
        assert s.fig_format == "svg"
        assert s._current_warnings == 0


class TestCmapValidation:
    @pytest.mark.parametrize("cls,default", [
        (PCAPlotScatterSettings, "viridis"),
        (PCAPlotProbabilityHeatmapSettings, "viridis"),
        (PCAPlotFESHeatmapSettings, "plasma"),
    ])
    def test_invalid_cmap_coerced(self, cls, default):
        s = cls.from_dict({"cmap": "not_a_cmap"})
        assert s.cmap == default
        assert s._current_warnings >= 1

    @pytest.mark.parametrize("cls", [
        PCAPlotScatterSettings, PCAPlotProbabilityHeatmapSettings,
        PCAPlotFESHeatmapSettings,
    ])
    def test_valid_cmap_preserved(self, cls):
        s = cls.from_dict({"cmap": "magma"})
        assert s.cmap == "magma"
        assert s._current_warnings == 0
