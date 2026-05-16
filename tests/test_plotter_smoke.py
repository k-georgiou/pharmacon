"""Layer 5 — end-to-end render smoke tests.

Each test renders one plot per family using matplotlib's Agg backend
(no display) on a mock PTA, then asserts the output file exists and is
non-empty. Marked @pytest.mark.slow so default `pytest` runs skip them
unless `pytest -m slow` is used.

Set MPLBACKEND=Agg in case the backend was already chosen by a previous
matplotlib import.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib  # noqa: E402
matplotlib.use("Agg", force=True)

import pytest  # noqa: E402

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import (  # noqa: E402
    build_pli_pta, build_ppi_pta, build_pca_pta, build_universal_pta,
)

from pharmacon.fileio.pta import PharmaconPTAFile  # noqa: E402
from pharmacon.constants.plots import (  # noqa: E402
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinProteinInteractionsHeatmap,
    PCAPlotScatterSettings,
    PlotUniversalSettings,
)


pytestmark = pytest.mark.slow


def _assert_nonempty_file_in(out_dir: Path) -> None:
    """Helper — at least one file was produced in out_dir."""
    files = [p for p in out_dir.iterdir() if p.is_file()]
    assert files, f"no output produced in {out_dir}"
    for p in files:
        assert p.stat().st_size > 0, f"empty output file: {p}"


# ---------------------------------------------------------------------------
# PLI stacked column 1 (modes table; non-merged)
# ---------------------------------------------------------------------------

def test_render_pli_stacked_column_1(tmp_path):
    from pharmacon.plotter.interactions import (
        plot_protein_ligand_interactions_stacked_column_1_from_file,
    )
    pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
    out_dir = tmp_path / "out_pli"
    out_dir.mkdir()

    settings = ProteinLigandInteractionsStackedColumn1PlotSettings.from_dict({
        "fig_dpi": 80,
        "fig_size_width": 4,
        "fig_size_height": 3,
        "fig_format": "png",
    })

    with PharmaconPTAFile(pta, mode="r") as f:
        plot_protein_ligand_interactions_stacked_column_1_from_file(
            f, group_name="pl_interactions", mode_name="mode1",
            settings=settings, out_dir=out_dir,
            attach_to_name="mode1", is_merged=False,
        )
    _assert_nonempty_file_in(out_dir)


# ---------------------------------------------------------------------------
# PPI heatmap (modes table; non-merged)
# ---------------------------------------------------------------------------

def test_render_ppi_heatmap(tmp_path):
    from pharmacon.plotter.interactions import (
        plot_protein_protein_heatmap_freq_from_file,
    )
    pta = build_ppi_pta(tmp_path / "ppi.pta", n_frames=80)
    out_dir = tmp_path / "out_ppi"
    out_dir.mkdir()

    settings = ProteinProteinInteractionsHeatmap.from_dict({
        "fig_dpi": 80,
        "fig_size_width": 4,
        "fig_size_height": 3,
        "fig_format": "png",
    })

    with PharmaconPTAFile(pta, mode="r") as f:
        plot_protein_protein_heatmap_freq_from_file(
            f, group_name="pp_interactions", mode_name="mode1",
            settings=settings, out_dir=out_dir,
            attach_to_name="mode1", is_merged=False,
        )
    _assert_nonempty_file_in(out_dir)


# ---------------------------------------------------------------------------
# PCA scatter
# ---------------------------------------------------------------------------

def test_render_pca_scatter(tmp_path):
    from pharmacon.plotter.universal import plot_pca_scatter_from_file
    pta = build_pca_pta(tmp_path / "pca.pta", n_frames=40, n_components=3)
    out_dir = tmp_path / "out_pca"
    out_dir.mkdir()

    settings = PCAPlotScatterSettings.from_dict({
        "fig_dpi": 80,
        "fig_size_width": 4,
        "fig_size_height": 3,
        "fig_format": "png",
        "pc_x": 1,
        "pc_y": 2,
    })

    with PharmaconPTAFile(pta, mode="r") as f:
        plot_pca_scatter_from_file(
            f, group_name="pca", settings=settings, out_dir=out_dir,
        )
    _assert_nonempty_file_in(out_dir)


# ---------------------------------------------------------------------------
# Universal timeseries (RMSD)
# ---------------------------------------------------------------------------

def test_render_universal_rmsd_timeseries(tmp_path):
    from pharmacon.plotter.universal import plot_pta_timeseries_from_file
    series = {"backbone": [1.0, 1.2, 1.1, 1.3, 1.4, 1.2, 1.1, 1.0]}
    pta = build_universal_pta(tmp_path / "rmsd.pta", group="rmsd", series=series)
    out_dir = tmp_path / "out_rmsd"
    out_dir.mkdir()

    settings = PlotUniversalSettings.from_dict({
        "fig_dpi": 80,
        "fig_size_width": 4,
        "fig_size_height": 3,
        "fig_format": "png",
    })

    with PharmaconPTAFile(pta, mode="r") as f:
        plot_pta_timeseries_from_file(
            f, group_name="rmsd", settings=settings, out_dir=out_dir,
            is_merged=False,
        )
    _assert_nonempty_file_in(out_dir)
