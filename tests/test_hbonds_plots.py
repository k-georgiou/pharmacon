"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for hydrogen-bond plotting.

H-bond `.pta` files were previously discovered by `plot pta` but never
dispatched (no plotter, no INIs). They now render five plots: a residue×residue
frequency heatmap (per mode), a pair-occupancy timeline, an H-bonds-per-frame
time series, and a ranked occupancy bar.
"""

from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_hbonds_frames_pta, DEFAULT_PPI_ROWS

from pharmacon.command_line.plot.pta import run
from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.plotter.hbonds import (
    plot_hbonds_count_per_frame_from_file,
    plot_hbonds_occupancy_from_file,
    plot_hbonds_network_from_file,
)
from pharmacon.constants.plots import (
    PlotSettingsBase, HBondsHeatmap, HBondsTimelinePairs,
    HBondsCountPerFrame, HBondsOccupancy, HBondsNetwork,
)


# H-bond rows are residue↔residue HYDROGEN-BOND pairs — same record shape as PPI.
HBOND_ROWS = tuple(
    ("HYDROGEN-BOND", rn1, rid1, ch1, rn2, rid2, ch2, freq)
    for (_label, rn1, rid1, ch1, rn2, rid2, ch2, freq) in DEFAULT_PPI_ROWS
)


def _build_hbonds(path: Path):
    # Writes the hbonds group with modes (from HBOND_ROWS) + per-frame records.
    return build_hbonds_frames_pta(path, rows=HBOND_ROWS, n_frames=60)


# --------------------------------------------------------------------------
# Settings: aliases register and titles are H-bond specific.
# --------------------------------------------------------------------------

class TestSettings:
    @pytest.mark.parametrize("alias,cls", [
        ("HBONDS-HEATMAP", HBondsHeatmap),
        ("HBONDS-TIMELINE-PAIRS", HBondsTimelinePairs),
        ("HBONDS-COUNT-PER-FRAME", HBondsCountPerFrame),
        ("HBONDS-OCCUPANCY", HBondsOccupancy),
        ("HBONDS-NETWORK", HBondsNetwork),
    ])
    def test_alias_resolves(self, alias, cls):
        assert PlotSettingsBase.resolve(alias) is cls

    def test_hbond_titles(self):
        assert HBondsHeatmap.from_dict({}).fig_title == "Hydrogen Bond Contact Frequency"
        assert HBondsTimelinePairs.from_dict({}).fig_title == "Hydrogen Bond Pairs – Timeline"
        assert HBondsCountPerFrame.from_dict({}).fig_title == "Hydrogen Bonds per Frame"
        assert HBondsOccupancy.from_dict({}).fig_title == "Hydrogen Bond Occupancy"
        assert HBondsNetwork.from_dict({}).fig_title == "Hydrogen Bond Network"

    def test_network_validates(self):
        s = HBondsNetwork.from_dict({"fig_format": "xyz", "layout": "bogus",
                                     "edge_cmap": "nope", "node_color": "notacolor"})
        assert s.fig_format == "png" and s.layout == "spring"
        assert s.edge_cmap == "viridis" and s.node_color == "#cfe8ff"
        assert s._current_warnings >= 1

    def test_count_defaults_to_frame_axis(self):
        assert HBondsCountPerFrame.from_dict({}).x_axis == "frame_index"

    def test_occupancy_validates(self):
        s = HBondsOccupancy.from_dict({"fig_format": "xyz", "bar_color": "nope",
                                       "grid_style": "wiggly"})
        assert s.fig_format == "png" and s.bar_color == "#3498db"
        assert s.grid_style == "dashed" and s._current_warnings >= 1


# --------------------------------------------------------------------------
# Plotters: the two H-bond-specific functions produce non-empty files.
# --------------------------------------------------------------------------

class TestPlotters:
    def test_count_per_frame_renders(self, tmp_path):
        pta = _build_hbonds(tmp_path / "hb.pta")
        out = tmp_path / "o1"; out.mkdir()
        s = HBondsCountPerFrame.from_dict({"fig_dpi": 80, "fig_size_width": 4,
                                           "fig_size_height": 3})
        with PharmaconPTAFile(pta, mode="r") as f:
            plot_hbonds_count_per_frame_from_file(f, group_name="hbonds",
                                                  settings=s, out_dir=out)
        files = list(out.glob("*.png"))
        assert files and files[0].stat().st_size > 0

    def test_occupancy_renders(self, tmp_path):
        pta = _build_hbonds(tmp_path / "hb.pta")
        out = tmp_path / "o2"; out.mkdir()
        s = HBondsOccupancy.from_dict({"fig_dpi": 80, "fig_size_width": 4,
                                       "fig_size_height": 4, "top_n": 10})
        with PharmaconPTAFile(pta, mode="r") as f:
            plot_hbonds_occupancy_from_file(f, group_name="hbonds",
                                            mode_name="mode2", settings=s, out_dir=out)
        files = list(out.glob("*.png"))
        assert files and files[0].stat().st_size > 0

    def test_network_renders(self, tmp_path):
        pta = _build_hbonds(tmp_path / "hb.pta")
        out = tmp_path / "o3"; out.mkdir()
        # low threshold so the mock's pairs survive into the graph
        s = HBondsNetwork.from_dict({"fig_dpi": 80, "fig_size_width": 5,
                                     "fig_size_height": 5, "threshold": 0.0})
        with PharmaconPTAFile(pta, mode="r") as f:
            plot_hbonds_network_from_file(f, group_name="hbonds",
                                          settings=s, out_dir=out)
        files = list(out.glob("*.png"))
        assert files and files[0].stat().st_size > 0

    def test_network_threshold_can_empty(self, tmp_path):
        pta = _build_hbonds(tmp_path / "hb.pta")
        out = tmp_path / "o4"; out.mkdir()
        # threshold above any mock occupancy → nothing to draw → RuntimeError.
        s = HBondsNetwork.from_dict({"fig_dpi": 80, "threshold": 1.0})
        with PharmaconPTAFile(pta, mode="r") as f:
            with pytest.raises(RuntimeError):
                plot_hbonds_network_from_file(f, group_name="hbonds",
                                              settings=s, out_dir=out)


# --------------------------------------------------------------------------
# End-to-end: `plot pta` dispatches an hbonds group to all five plots.
# --------------------------------------------------------------------------

def _make_args(input_pta: Path, out_dir: Path) -> Namespace:
    return Namespace(
        input=str(input_pta), output=out_dir, is_merged=False, maxwarnings=9999,
        config=None, config_overrides={}, overwrite=True,
        log=str(out_dir / "plot.log"),
        terminal_logging_level="INFO", file_logging_level="TRACE",
    )


class TestFrameTimeStorage:
    def test_time_ps_stored_when_provided(self, tmp_path):
        import h5py
        rec = ("HYDROGEN-BOND", 0, "N", 1, "N", "BB", "ALA", 1, "A", "S",
               10, "O", 2, "O", "BB", "LEU", 2, "A", "S", 2.0, True)
        path = tmp_path / "t.pta"
        with PharmaconPTAFile(path, overwrite=True, mode="a",
                              command="Trajectory Analysis", subcommand="h-bonds") as pta:
            pta.add_file_metadata({"is_merged": "False"})
            pta.create_group("hbonds")
            pta.write_frame_interactions(frame_index=0, group_name="hbonds",
                                         interactions=[rec], overwrite=True,
                                         time_ps=12.5)
            pta.write_frame_interactions(frame_index=1, group_name="hbonds",
                                         interactions=[rec], overwrite=True)
        f = h5py.File(path, "r")
        assert f["hbonds/frame_0/interactions"].attrs["time_ps"] == "12.5"
        assert "time_ps" not in f["hbonds/frame_1/interactions"].attrs

    def test_count_plot_time_axis_autolabels(self, tmp_path):
        import re
        import html
        rec = ("HYDROGEN-BOND", 0, "N", 1, "N", "BB", "ALA", 1, "A", "S",
               10, "O", 2, "O", "BB", "LEU", 2, "A", "S", 2.0, True)
        path = tmp_path / "t.pta"
        with PharmaconPTAFile(path, overwrite=True, mode="a",
                              command="Trajectory Analysis", subcommand="h-bonds") as pta:
            pta.add_file_metadata({"is_merged": "False"})
            pta.create_group("hbonds")
            for fr in range(5):
                pta.write_frame_interactions(frame_index=fr, group_name="hbonds",
                                             interactions=[rec], overwrite=True,
                                             time_ps=fr * 5.0)
        out = tmp_path / "o"; out.mkdir()
        s = HBondsCountPerFrame.from_dict({"x_axis": "time_ns", "fig_format": "svg",
                                           "fig_dpi": 80, "fig_size_width": 4,
                                           "fig_size_height": 3})
        with PharmaconPTAFile(path, mode="r") as f:
            plot_hbonds_count_per_frame_from_file(f, group_name="hbonds",
                                                  settings=s, out_dir=out)
        svg = next(out.glob("*.svg"))
        text = html.unescape(re.sub("<[^>]+>", "", " ".join(
            re.findall(r"<text[^>]*>(.*?)</text>", svg.read_text(), re.S))))
        assert "Time (ns)" in text


@pytest.mark.slow
def test_plot_pta_renders_all_hbond_plots(tmp_path):
    pta = _build_hbonds(tmp_path / "hb.pta")
    out = tmp_path / "out"; out.mkdir()
    run(_make_args(pta, out))

    names = " ".join(p.name for p in out.glob("*.png"))
    # heatmap (per mode), occupancy, timeline, count-per-frame, network
    for token in ("hbonds_heatmap", "hbonds_occupancy",
                  "hbonds_timeline_pairs", "hbonds_count_per_frame",
                  "hbonds_network"):
        assert token in names, f"missing '{token}' in {names}"
