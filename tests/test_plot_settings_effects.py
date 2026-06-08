"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Layer 6 — settings-effect tests (every plot family).

The smoke layer (test_plotter_smoke.py) only checks that *a* file is
produced. These tests go one step further: they render with a specific
setting and then inspect the written artifact to prove the setting
actually changed the output.

A single CASES registry drives the whole suite so every plot type gets the
same battery:

  - fig_format        -> the file type that hits disk (PNG magic / SVG / PDF)
  - fig_dpi           -> higher dpi yields strictly more pixels (PNG)
  - fig_transparent   -> the saved PNG has a transparent background (or not)
  - fig_title /
    disable_title     -> the title string appears / disappears in the SVG
                         (skipped for pie charts, which title per-residue)
  - colour overrides  -> a pinned hex shows up in the SVG (where a single
                         colour applies; cmap-based plots are exempt)

Plus an RMSF-specific exact-pixel check (its savefig path has no
bbox_inches='tight', so pixels == inches * dpi).

These render via matplotlib's Agg backend (no display) on mock PTAs and
are marked @pytest.mark.slow alongside the other render tests.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib  # noqa: E402
matplotlib.use("Agg", force=True)
import matplotlib.image as mpimg  # noqa: E402

import pytest  # noqa: E402

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import (  # noqa: E402
    build_universal_pta, build_rmsf_pta, build_pca_pta,
    build_pli_frames_pta, build_ppi_frames_pta,
)

from pharmacon.fileio.pta import PharmaconPTAFile  # noqa: E402
import pharmacon.plotter.universal as U  # noqa: E402
import pharmacon.plotter.interactions as I  # noqa: E402
from pharmacon.constants.plots import (  # noqa: E402
    PlotUniversalSettings,
    RMSFPlotSettings,
    PCAPlotTimeSeriesSettings,
    PCAPlotScatterSettings,
    PCAPlotVarianceRatioSettings,
    PCAPlotProbabilityHeatmapSettings,
    PCAPlotFESHeatmapSettings,
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinLigandInteractionsStackedColumn2PlotSettings,
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsPieCharts1PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsStackedColumnSettings,
    ProteinProteinInteractionsTimelinePairs,
)


pytestmark = pytest.mark.slow

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
MARKER_HEX = "#123456"          # distinctive; unlikely to be any default
MARKER_TITLE = "ZZ_TITLE_MARKER_ZZ"


# --------------------------------------------------------------------------
# Case registry — one entry per plot type
# --------------------------------------------------------------------------

@dataclass
class Case:
    cid: str
    build: Callable[[Path], Path]          # (dir) -> pta path
    settings_cls: type
    render: Callable                       # (open_file, out_dir, settings) -> None
    title: bool = True                     # fig_title reliably lands in the SVG
    color: Optional[Tuple[str, object]] = None  # (settings_key, hex/list) or None


def _uni(p):  return build_universal_pta(p / "u.pta", group="rmsd",
                                         series={"bb": [1.0, 1.2, 1.1, 1.3, 1.4, 1.2, 1.0, 1.1]})
def _rmsf(p): return build_rmsf_pta(p / "r.pta")
def _pca(p):  return build_pca_pta(p / "pca.pta", n_frames=200, n_components=3)
def _pli(p):  return build_pli_frames_pta(p / "pli.pta", n_frames=100)
def _ppi(p):  return build_ppi_frames_pta(p / "ppi.pta", n_frames=100)


CASES = [
    Case("universal_timeseries", _uni, PlotUniversalSettings,
         lambda f, o, s: U.plot_pta_timeseries_from_file(
             f, group_name="rmsd", settings=s, out_dir=o, is_merged=False),
         color=("line_colors", [MARKER_HEX])),

    Case("rmsf", _rmsf, RMSFPlotSettings,
         lambda f, o, s: U.plot_pta_rmsf_from_file(
             f, group_name="rmsf", settings=s, out_dir=o, is_merged=False),
         color=("colors_by_label", f"calpha:{MARKER_HEX}, backbone:#abcdef")),

    Case("pca_timeseries", _pca, PCAPlotTimeSeriesSettings,
         lambda f, o, s: U.plot_pca_timeseries_from_file(
             f, group_name="pca", settings=s, out_dir=o)),

    Case("pca_scatter", _pca, PCAPlotScatterSettings,
         lambda f, o, s: U.plot_pca_scatter_from_file(
             f, group_name="pca", settings=s, out_dir=o)),

    Case("pca_variance_ratio", _pca, PCAPlotVarianceRatioSettings,
         lambda f, o, s: U.plot_pca_variance_ratio_from_file(
             f, group_name="pca", settings=s, out_dir=o)),

    Case("pca_probability", _pca, PCAPlotProbabilityHeatmapSettings,
         lambda f, o, s: U.plot_pca_probability_from_file(
             f, group_name="pca", settings=s, out_dir=o)),

    Case("pca_fes", _pca, PCAPlotFESHeatmapSettings,
         lambda f, o, s: U.plot_pca_fes_from_file(
             f, group_name="pca", settings=s, out_dir=o)),

    Case("pli_stacked_column_1", _pli, ProteinLigandInteractionsStackedColumn1PlotSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_stacked_column_1_from_file(
             f, group_name="pl_interactions", mode_name="mode1", settings=s,
             out_dir=o, attach_to_name="mode1", is_merged=False),
         color=("color_hydrophobic", MARKER_HEX)),

    Case("pli_stacked_column_2", _pli, ProteinLigandInteractionsStackedColumn2PlotSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_stacked_column_2_from_file(
             f, group_name="pl_interactions", settings=s, out_dir=o, is_merged=False),
         color=("color_backbone", MARKER_HEX)),

    Case("pli_heatmap_1", _pli, ProteinLigandInteractionsHeatmap1PlotSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_heatmap_1_from_file(
             f, group_name="pl_interactions", settings=s, out_dir=o, is_merged=False)),

    Case("pli_heatmap_2", _pli, ProteinLigandInteractionsHeatmap2PlotSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_heatmap_2_from_file(
             f, group_name="pl_interactions", settings=s, out_dir=o, is_merged=False)),

    Case("pli_pie_charts", _pli, ProteinLigandInteractionsPieCharts1PlotSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_pie_charts_from_file(
             f, group_name="pl_interactions", settings=s, out_dir=o, is_merged=False),
         title=False,  # each pie is titled by its residue, not fig_title
         color=("color_backbone", MARKER_HEX)),  # pies colour by SC/BB, not interaction

    Case("pli_ligand_monitor", _pli, ProteinLigandInteractionsLigandMonitorSettings,
         lambda f, o, s: I.plot_protein_ligand_interactions_ligand_monitor_from_file(
             f, group_name="pl_interactions", settings=s, out_dir=o, is_merged=False)),

    Case("ppi_heatmap", _ppi, ProteinProteinInteractionsHeatmap,
         lambda f, o, s: I.plot_protein_protein_heatmap_freq_from_file(
             f, group_name="pp_interactions", mode_name="mode1", settings=s,
             out_dir=o, attach_to_name="mode1", is_merged=False)),

    Case("ppi_stacked_column", _ppi, ProteinProteinInteractionsStackedColumnSettings,
         lambda f, o, s: I.plot_protein_protein_interactions_stacked_column_from_file(
             f, group_name="pp_interactions", mode_name="mode1", settings=s,
             out_dir=o, attach_to_name="mode1", is_merged=False),
         color=("color_hydrophobic", MARKER_HEX)),

    Case("ppi_timeline_pairs", _ppi, ProteinProteinInteractionsTimelinePairs,
         lambda f, o, s: I.plot_protein_protein_timeline_pairs_from_file(
             f, group_name="pp_interactions", settings=s, out_dir=o, is_merged=False)),
]

CASES_BY_ID = {c.cid: c for c in CASES}
ALL_IDS = [c.cid for c in CASES]
TITLE_IDS = [c.cid for c in CASES if c.title]
COLOR_IDS = [c.cid for c in CASES if c.color is not None]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

_BASE = {"fig_dpi": 90, "fig_size_width": 5, "fig_size_height": 4, "fig_format": "png"}


def _render(case: Case, tmp_path: Path, out_name: str, **overrides) -> list[Path]:
    """Build the case's mock PTA, render with overrides, return output files."""
    src = tmp_path / "src"
    src.mkdir(exist_ok=True)
    pta = case.build(src)

    out_dir = tmp_path / out_name
    out_dir.mkdir(parents=True, exist_ok=True)

    settings = case.settings_cls.from_dict({**_BASE, **overrides})
    with PharmaconPTAFile(pta, mode="r") as f:
        case.render(f, out_dir, settings)

    files = sorted((p for p in out_dir.rglob("*") if p.is_file()),
                   key=lambda p: p.name)
    assert files, f"{case.cid}: no output produced"
    return files


def _first(files: list[Path], suffix: str) -> Path:
    matches = [p for p in files if p.suffix == suffix]
    assert matches, f"no {suffix} file among {[p.name for p in files]}"
    return matches[0]


# --------------------------------------------------------------------------
# fig_format -> file type on disk
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cid", ALL_IDS)
def test_fig_format_png(tmp_path, cid):
    case = CASES_BY_ID[cid]
    files = _render(case, tmp_path, "png", fig_format="png")
    png = _first(files, ".png")
    assert png.read_bytes()[:8] == PNG_MAGIC


@pytest.mark.parametrize("cid", ALL_IDS)
def test_fig_format_svg(tmp_path, cid):
    case = CASES_BY_ID[cid]
    files = _render(case, tmp_path, "svg", fig_format="svg")
    svg = _first(files, ".svg")
    assert b"<svg" in svg.read_bytes()


# --------------------------------------------------------------------------
# fig_dpi -> more pixels
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cid", ALL_IDS)
def test_higher_dpi_yields_more_pixels(tmp_path, cid):
    case = CASES_BY_ID[cid]
    low = mpimg.imread(str(_first(
        _render(case, tmp_path, "lo", fig_format="png", fig_dpi=60), ".png")))
    high = mpimg.imread(str(_first(
        _render(case, tmp_path, "hi", fig_format="png", fig_dpi=120), ".png")))
    # Strictly more pixels in at least one dimension; never fewer in either.
    assert high.shape[0] >= low.shape[0] and high.shape[1] >= low.shape[1]
    assert high.shape[0] * high.shape[1] > low.shape[0] * low.shape[1]


# --------------------------------------------------------------------------
# fig_transparent -> background transparency
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cid", ALL_IDS)
def test_fig_transparent_controls_background(tmp_path, cid):
    case = CASES_BY_ID[cid]
    transp = mpimg.imread(str(_first(
        _render(case, tmp_path, "t", fig_format="png", fig_transparent=True), ".png")))
    opaque = mpimg.imread(str(_first(
        _render(case, tmp_path, "o", fig_format="png", fig_transparent=False), ".png")))

    # Transparent render: RGBA with at least one fully see-through pixel.
    assert transp.shape[2] == 4, f"{cid}: transparent PNG should be RGBA"
    assert float(transp[..., 3].min()) == pytest.approx(0.0, abs=1e-3)

    # Opaque render: if it carries an alpha channel, nothing is transparent.
    if opaque.shape[2] == 4:
        assert float(opaque[..., 3].min()) == pytest.approx(1.0, abs=1e-3)


# --------------------------------------------------------------------------
# fig_title / disable_title -> title text in the SVG
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cid", TITLE_IDS)
def test_fig_title_appears_in_svg(tmp_path, cid):
    case = CASES_BY_ID[cid]
    files = _render(case, tmp_path, "out", fig_format="svg",
                    fig_title=MARKER_TITLE, disable_title=False)
    assert any(MARKER_TITLE in p.read_text() for p in files if p.suffix == ".svg")


@pytest.mark.parametrize("cid", TITLE_IDS)
def test_disable_title_removes_title_from_svg(tmp_path, cid):
    case = CASES_BY_ID[cid]
    files = _render(case, tmp_path, "out", fig_format="svg",
                    fig_title=MARKER_TITLE, disable_title=True)
    assert not any(MARKER_TITLE in p.read_text() for p in files if p.suffix == ".svg")


# --------------------------------------------------------------------------
# colour overrides -> pinned hex shows up in the SVG
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cid", COLOR_IDS)
def test_colour_override_appears_in_svg(tmp_path, cid):
    case = CASES_BY_ID[cid]
    key, value = case.color
    files = _render(case, tmp_path, "out", fig_format="svg", **{key: value})
    assert any(MARKER_HEX in p.read_text() for p in files if p.suffix == ".svg"), (
        f"{cid}: setting {key}={value!r} did not reach the SVG"
    )


# --------------------------------------------------------------------------
# RMSF-specific exact pixel check (no bbox_inches='tight' on this savefig)
# --------------------------------------------------------------------------

def test_rmsf_fig_size_times_dpi_sets_exact_png_pixels(tmp_path):
    case = CASES_BY_ID["rmsf"]
    png = _first(_render(
        case, tmp_path, "exact",
        fig_format="png", fig_size_width=5, fig_size_height=3, fig_dpi=100), ".png")
    arr = mpimg.imread(str(png))   # (height, width, channels)
    assert abs(arr.shape[1] - 500) <= 2, f"width {arr.shape[1]} != 5in*100dpi"
    assert abs(arr.shape[0] - 300) <= 2, f"height {arr.shape[0]} != 3in*100dpi"
