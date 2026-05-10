"""
Tests for pharmacon.analyzer.frame_logger.log_every_frame

Note on MemoryReader: MDAnalysis.coordinates.memory.MemoryReader overrides
_apply_transformations as a no-op, so transformations stored in
_transformations are never invoked during normal iteration over a MemoryReader.
The frame logger works by hooking into this pipeline.

These tests work around this limitation by:
  1. Entering the context manager so _tap is registered in _transformations.
  2. Manually calling the _tap closure for each simulated frame.
  3. Verifying the logger is called with the expected arguments.

This correctly tests the frame-counting / slice / every logic in isolation.
"""
import pytest
import numpy as np
import MDAnalysis as Mda
from unittest.mock import Mock
from MDAnalysis.coordinates.memory import MemoryReader

from pharmacon.analyzer.frame_logger import log_every_frame


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_universe(n_frames: int) -> Mda.Universe:
    """Minimal single-atom Universe backed by MemoryReader."""
    u = Mda.Universe.empty(
        n_atoms=1, n_residues=1, n_segments=1,
        atom_resindex=[0], residue_segindex=[0], trajectory=True,
    )
    u.add_TopologyAttr("names", ["A0"])
    u.add_TopologyAttr("resnames", ["MOL"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    positions = np.zeros((n_frames, 1, 3), dtype=np.float32)
    u.load_new(positions, format=MemoryReader)
    return u


def _tap_fn(u: Mda.Universe) -> callable:
    """Return the _tap closure registered by the active log_every_frame context."""
    return u.trajectory._transformations[-1]


def _make_ts(frame: int, time: float | None = None):
    """Create a minimal mock timestep."""
    ts = Mock()
    ts.frame = frame
    ts.time = time
    return ts


def _simulate(tap, n_frames: int):
    """Drive tap through n_frames mock timesteps (mimics what __next__ would do)."""
    for f in range(n_frames):
        tap(_make_ts(f, float(f)))


# ---------------------------------------------------------------------------
# Context manager basics (do not depend on transformation pipeline)
# ---------------------------------------------------------------------------

class TestFrameLoggerContextManager:
    def test_no_exception_on_entry_exit(self):
        u = _make_universe(3)
        with log_every_frame(u, Mock()):
            pass

    def test_tap_is_registered_in_transformations(self):
        u = _make_universe(3)
        before = len(getattr(u.trajectory, "_transformations", []))
        with log_every_frame(u, Mock()):
            after = len(getattr(u.trajectory, "_transformations", []))
            assert after == before + 1

    def test_tap_is_callable(self):
        u = _make_universe(3)
        with log_every_frame(u, Mock()):
            tap = _tap_fn(u)
            assert callable(tap)

    def test_transformations_restored_after_context(self):
        u = _make_universe(3)
        original = list(getattr(u.trajectory, "_transformations", []))
        with log_every_frame(u, Mock()):
            pass
        assert list(getattr(u.trajectory, "_transformations", [])) == original

    def test_transformations_restored_on_exception(self):
        u = _make_universe(3)
        original = list(getattr(u.trajectory, "_transformations", []))
        try:
            with log_every_frame(u, Mock()):
                raise RuntimeError("boom")
        except RuntimeError:
            pass
        assert list(getattr(u.trajectory, "_transformations", [])) == original

    def test_no_debug_calls_without_frames(self):
        u = _make_universe(3)
        mock_log = Mock()
        with log_every_frame(u, mock_log):
            pass
        mock_log.debug.assert_not_called()


# ---------------------------------------------------------------------------
# Logging logic — every=1 (via manual _tap invocation)
# ---------------------------------------------------------------------------

class TestFrameLoggerEvery1:
    def test_4_frames_produces_4_debug_calls(self):
        # Frames 0-3 manually fed through _tap:
        # frame 0: first, no log
        # frame 1 accessed → done_f=0, count=0 (0%1==0) → log; count=1
        # frame 2 accessed → done_f=1, count=1 (1%1==0) → log; count=2
        # frame 3 accessed → done_f=2, count=2 (2%1==0) → log; count=3
        # finally: last_done_f=2 → DONE log
        # total = 3 normal + 1 DONE = 4
        u = _make_universe(4)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=1):
            _simulate(_tap_fn(u), 4)
        assert mock_log.debug.call_count == 4

    def test_done_message_in_last_call(self):
        u = _make_universe(3)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=1):
            _simulate(_tap_fn(u), 3)
        last = mock_log.debug.call_args_list[-1]
        assert "DONE" in last[0][0]

    def test_custom_prefix_appears_in_all_calls(self):
        # logger.debug("%s #%d", prefix, frame) → c[0][1] is the prefix
        u = _make_universe(3)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=1, prefix="myframe"):
            _simulate(_tap_fn(u), 3)
        for c in mock_log.debug.call_args_list:
            assert c[0][1] == "myframe"

    def test_single_frame_produces_no_log(self):
        # Only 1 frame: _tap called once (first frame), no "done" frame seen
        u = _make_universe(1)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=1):
            _simulate(_tap_fn(u), 1)
        # last_done_f is None → finally block doesn't log anything
        mock_log.debug.assert_not_called()


# ---------------------------------------------------------------------------
# every parameter
# ---------------------------------------------------------------------------

class TestFrameLoggerEvery:
    def test_every_2_with_5_frames(self):
        # Frames 0-4 (4 "done" frames: 0,1,2,3):
        # count=0 → frame 0 logged; count=1 → skip; count=2 → frame 2 logged; count=3 → skip
        # finally: last_done_f=3 → DONE
        # total = 2 normal + 1 DONE = 3
        u = _make_universe(5)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=2):
            _simulate(_tap_fn(u), 5)
        assert mock_log.debug.call_count == 3

    def test_every_10_with_3_frames(self):
        # 3 frames → 2 "done" frames (0 and 1):
        # count=0 → frame 0 logged (0%10==0); count=1 → skip
        # finally: last_done_f=1 → DONE
        # total = 1 normal + 1 DONE = 2
        u = _make_universe(3)
        mock_log = Mock()
        with log_every_frame(u, mock_log, every=10):
            _simulate(_tap_fn(u), 3)
        assert mock_log.debug.call_count == 2

    def test_every_1_vs_every_2_call_counts_differ(self):
        u1 = _make_universe(6)
        u2 = _make_universe(6)
        log1, log2 = Mock(), Mock()
        with log_every_frame(u1, log1, every=1):
            _simulate(_tap_fn(u1), 6)
        with log_every_frame(u2, log2, every=2):
            _simulate(_tap_fn(u2), 6)
        assert log1.debug.call_count > log2.debug.call_count


# ---------------------------------------------------------------------------
# start / stop slicing
# ---------------------------------------------------------------------------

class TestFrameLoggerSlice:
    def test_start_limits_logging_to_later_frames(self):
        # 6 frames, start=3: only frames 3,4 are "done" in-slice (frame 5 triggers frame 4 done)
        u = _make_universe(6)
        mock_log = Mock()
        with log_every_frame(u, mock_log, start=3, every=1):
            _simulate(_tap_fn(u), 6)
        normal_calls = [c for c in mock_log.debug.call_args_list if "DONE" not in c[0][0]]
        for c in normal_calls:
            frame_num = c[0][2]  # logger.debug("%s #%d", prefix, frame_num)
            assert frame_num >= 3, f"Expected frame >= 3, got {frame_num}"

    def test_stop_limits_logging_to_earlier_frames(self):
        # 6 frames, stop=3 (exclusive [0,3)): frames 0,1,2 in slice
        u = _make_universe(6)
        mock_log = Mock()
        with log_every_frame(u, mock_log, stop=3, every=1):
            _simulate(_tap_fn(u), 6)
        normal_calls = [c for c in mock_log.debug.call_args_list if "DONE" not in c[0][0]]
        for c in normal_calls:
            frame_num = c[0][2]  # logger.debug("%s #%d", prefix, frame_num)
            assert frame_num < 3, f"Expected frame < 3, got {frame_num}"

    def test_start_equals_stop_no_frames_in_slice(self):
        # start=stop → empty slice → no normal logs, no DONE (last_done_f stays None)
        u = _make_universe(5)
        mock_log = Mock()
        with log_every_frame(u, mock_log, start=2, stop=2, every=1):
            _simulate(_tap_fn(u), 5)
        mock_log.debug.assert_not_called()
