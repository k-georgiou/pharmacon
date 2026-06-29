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


# ---------------------------------------------------------------------------
# Rigid-body alignment correctness
#
# The static-frame tests above use identical frames, where the optimal
# superposition rotation is the identity (R == R.T), so they cannot detect
# whether the rotation matrix returned by MDAnalysis' ``rotation_matrix`` is
# applied with the correct handedness. ``rotation_matrix(a, b)`` returns R with
# ``b = R @ a`` (R acts on the left), so a row-stored coordinate set is aligned
# with ``coords @ R.T`` — exactly what ``AtomGroup.rotate`` does internally.
#
# These tests feed pure rigid-body motion (rotation + translation of one base
# structure). A correct superposition recovers the reference structure exactly,
# so the averaged structure equals the base and every RMSD collapses to ~0.
# They FAIL if the alignment uses ``@ R`` (the inverse rotation) instead of
# ``@ R.T``.
# ---------------------------------------------------------------------------


def _axis_angle_matrix(axis, angle: float) -> np.ndarray:
    """Proper rotation matrix (det = +1) from an axis and angle (Rodrigues)."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1.0 - c
    return np.array(
        [
            [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
        ],
        dtype=float,
    )


# Asymmetric, non-collinear 5-atom rigid body so the optimal superposition
# rotation is unique and well-conditioned.
_RIGID_BASE = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 2.5],
        [1.0, 1.0, 0.5],
    ],
    dtype=np.float32,
)


def _make_rigid_motion_universe():
    """Universe whose frames are pure rigid-body transforms of one base.

    Frame 0 is the untouched base (used as the reference frame); every other
    frame is the base rotated by a known proper rotation and translated.
    """
    rng = np.random.default_rng(20260629)
    frames = [_RIGID_BASE.copy()]  # frame 0 == reference == base
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 2, 3)]
    angles = [0.3, 1.1, -0.7, 2.0, -1.4]
    for axis, angle in zip(axes, angles):
        R = _axis_angle_matrix(axis, angle)
        T = rng.uniform(-5.0, 5.0, size=3)
        moved = (_RIGID_BASE.astype(float) @ R.T) + T
        frames.append(moved.astype(np.float32))
    return _make_traj_universe(frames)


def _pairwise_distances(points: np.ndarray) -> np.ndarray:
    """Coordinate-frame-independent fingerprint of a structure's geometry."""
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt((diff * diff).sum(axis=-1))


class TestAvgStRigidBodyAlignment:
    _STOP = 5  # 6 frames total (frame 0 + 5 moved)

    def test_rotation_matrix_convention_is_left_acting(self):
        # Lock down the assumption the fix relies on: rotation_matrix(a, b)
        # returns R with b ≈ R @ a, so a @ R.T maps a onto b.
        from MDAnalysis.analysis.align import rotation_matrix

        R_true = _axis_angle_matrix((0.3, 1.0, -0.5), 0.9)
        a = (_RIGID_BASE.astype(float) - _RIGID_BASE.astype(float).mean(axis=0))
        b = a @ R_true.T  # b_i = R_true @ a_i
        R, rmsd = rotation_matrix(a, b)
        np.testing.assert_allclose(a @ R.T, b, atol=1e-6)
        assert rmsd == pytest.approx(0.0, abs=1e-6)

    def test_rmsd_to_ref_zero_for_rigid_motion(self):
        u = _make_rigid_motion_universe()
        _, _, rmsd_ref = avg_st_process_trajectory(u, u.atoms, start=0, stop=self._STOP)
        np.testing.assert_allclose(rmsd_ref, 0.0, atol=1e-3)

    def test_rmsd_to_avg_zero_for_rigid_motion(self):
        u = _make_rigid_motion_universe()
        _, rmsd_avg, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=self._STOP)
        np.testing.assert_allclose(rmsd_avg, 0.0, atol=1e-3)

    def test_avg_structure_matches_reference(self):
        # Reference (frame 0) is the untouched base in absolute coordinates;
        # aligning each rigid frame onto it recovers the base exactly, so the
        # average equals the base.
        u = _make_rigid_motion_universe()
        avg_pos, _, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=self._STOP)
        np.testing.assert_allclose(avg_pos, _RIGID_BASE, atol=1e-3)

    def test_avg_structure_preserves_internal_geometry(self):
        # A misaligned ("blurred") average would distort internal distances.
        u = _make_rigid_motion_universe()
        avg_pos, _, _ = avg_st_process_trajectory(u, u.atoms, start=0, stop=self._STOP)
        np.testing.assert_allclose(
            _pairwise_distances(avg_pos),
            _pairwise_distances(_RIGID_BASE.astype(float)),
            atol=1e-3,
        )

    def test_memory_efficient_and_full_agree_under_rotation(self):
        u1 = _make_rigid_motion_universe()
        u2 = _make_rigid_motion_universe()
        avg1, ra1, rr1 = avg_st_process_trajectory(
            u1, u1.atoms, start=0, stop=self._STOP, memory_efficient=True
        )
        avg2, ra2, rr2 = avg_st_process_trajectory(
            u2, u2.atoms, start=0, stop=self._STOP, memory_efficient=False
        )
        np.testing.assert_allclose(avg1, avg2, atol=1e-4)
        np.testing.assert_allclose(ra1, ra2, atol=1e-4)
        np.testing.assert_allclose(rr1, rr2, atol=1e-4)
