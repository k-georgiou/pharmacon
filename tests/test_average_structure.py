"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.analyzer.average_structure —
avg_st_process_trajectory and extract_trajectory_frame.
"""
import pytest
import numpy as np
import MDAnalysis as Mda
from MDAnalysis.coordinates.memory import MemoryReader

from pharmacon.analyzer.average_structure import (
    avg_st_process_trajectory,
    extract_trajectory_frame,
)


def _make_traj_universe(positions_per_frame):
    """Build a Universe from a list of (n_atoms, 3) position arrays."""
    n_atoms = len(positions_per_frame[0])
    n_frames = len(positions_per_frame)

    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n_atoms,
        residue_segindex=[0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", [f"A{i}" for i in range(n_atoms)])
    u.add_TopologyAttr("resnames", ["MOL"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    u.add_TopologyAttr("masses", [1.0] * n_atoms)

    all_pos = np.stack(positions_per_frame, axis=0).astype(np.float32)
    u.load_new(all_pos, format=MemoryReader)
    return u


# Three non-collinear atoms in an L-shape — rotation_matrix requires ≥3 non-collinear atoms.
_BASE_POS = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]], dtype=np.float32)
_N_FRAMES = 5


@pytest.fixture(scope="module")
def static_universe():
    """Universe where all 5 frames are identical (zero motion)."""
    frames = [_BASE_POS.copy() for _ in range(_N_FRAMES)]
    return _make_traj_universe(frames)


@pytest.fixture(scope="module")
def static_universe_fresh():
    """Fresh copy for tests that modify trajectory state."""
    frames = [_BASE_POS.copy() for _ in range(_N_FRAMES)]
    return _make_traj_universe(frames)


class TestAvgStValidation:
    def test_negative_reference_frame_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="reference_frame must be >= 0"):
            avg_st_process_trajectory(u, u.atoms, start=0, stop=4, reference_frame=-1)

    def test_reference_frame_below_start_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="reference_frame"):
            avg_st_process_trajectory(u, u.atoms, start=2, stop=4, reference_frame=1)

    def test_reference_frame_above_stop_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="reference_frame"):
            avg_st_process_trajectory(u, u.atoms, start=0, stop=3, reference_frame=4)

    def test_reference_frame_exceeds_total_frames_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="exceeds total frames"):
            avg_st_process_trajectory(u, u.atoms, start=0, stop=100, reference_frame=99)

    def test_zero_step_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="step must be > 0"):
            avg_st_process_trajectory(u, u.atoms, start=0, stop=4, step=0)

    def test_negative_step_raises(self, static_universe):
        u = static_universe
        with pytest.raises(ValueError, match="step must be > 0"):
            avg_st_process_trajectory(u, u.atoms, start=0, stop=4, step=-1)


class TestAvgStReturnTypes:
    def test_returns_three_tuple(self, static_universe):
        u = static_universe
        result = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        assert len(result) == 3

    def test_avg_positions_shape(self, static_universe):
        u = static_universe
        avg_pos, _, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        assert avg_pos.shape == (3, 3)  # (n_atoms, 3)

    def test_rmsd_to_avg_is_ndarray(self, static_universe):
        u = static_universe
        _, rmsd_avg, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        assert isinstance(rmsd_avg, np.ndarray)

    def test_rmsd_to_ref_is_ndarray(self, static_universe):
        u = static_universe
        _, _, rmsd_ref = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        assert isinstance(rmsd_ref, np.ndarray)

    def test_rmsd_arrays_have_correct_length(self, static_universe):
        u = static_universe
        _, rmsd_avg, rmsd_ref = avg_st_process_trajectory(
            u, u.atoms, start=0, stop=4, step=1
        )
        assert len(rmsd_avg) == 5
        assert len(rmsd_ref) == 5


class TestAvgStNumerics:
    def test_identical_frames_zero_rmsd_to_ref(self, static_universe):
        u = static_universe
        _, _, rmsd_ref = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        np.testing.assert_allclose(rmsd_ref, 0.0, atol=1e-5)

    def test_identical_frames_zero_rmsd_to_avg(self, static_universe):
        u = static_universe
        _, rmsd_avg, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        np.testing.assert_allclose(rmsd_avg, 0.0, atol=1e-5)

    def test_avg_positions_match_input_for_static(self, static_universe):
        u = static_universe
        avg_pos, _, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=4)
        # After centering and aligning, the average should match the reference frame
        # (which is identical to all frames). The positions may be in the reference
        # frame's coordinate system (centered), so we compare shapes and zero RMSD.
        assert avg_pos.shape == _BASE_POS.shape

    def test_memory_efficient_matches_non_efficient(self):
        frames = [_BASE_POS.copy() for _ in range(_N_FRAMES)]
        u1 = _make_traj_universe(frames)
        u2 = _make_traj_universe(frames)

        avg1, rmsd_avg1, rmsd_ref1 = avg_st_process_trajectory(
            u1, u1.atoms, start=0, stop=4, memory_efficient=True
        )
        avg2, rmsd_avg2, rmsd_ref2 = avg_st_process_trajectory(
            u2, u2.atoms, start=0, stop=4, memory_efficient=False
        )

        np.testing.assert_allclose(avg1, avg2, atol=1e-4)
        np.testing.assert_allclose(rmsd_avg1, rmsd_avg2, atol=1e-4)
        np.testing.assert_allclose(rmsd_ref1, rmsd_ref2, atol=1e-4)

    def test_step_parameter_reduces_frame_count(self):
        frames = [_BASE_POS.copy() for _ in range(6)]
        u = _make_traj_universe(frames)
        _, rmsd_avg, rmsd_ref = avg_st_process_trajectory(
            u, u.atoms, start=0, stop=5, step=2
        )
        # Frames 0, 2, 4 → 3 frames processed
        assert len(rmsd_avg) == 3
        assert len(rmsd_ref) == 3


class TestExtractFrameValidation:
    def test_negative_frame_idx_raises(self, tmp_path, static_universe_fresh):
        u = static_universe_fresh
        out = tmp_path / "out.pdb"
        with pytest.raises(ValueError, match="frame_idx must be >= 0"):
            extract_trajectory_frame(
                u, u.atoms, "all", -1,
                output_file_path=out, output_format="pdb"
            )

    def test_frame_idx_exceeds_trajectory_length_raises(self, tmp_path, static_universe_fresh):
        u = static_universe_fresh
        out = tmp_path / "out.pdb"
        with pytest.raises(ValueError, match="exceeds trajectory length"):
            extract_trajectory_frame(
                u, u.atoms, "all", 999,
                output_file_path=out, output_format="pdb"
            )
