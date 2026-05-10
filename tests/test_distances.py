"""
Tests for pharmacon.analyzer.distances.calculate_frame_distances
"""
import pytest
import numpy as np
import MDAnalysis as Mda

from pharmacon.analyzer.distances import calculate_frame_distances


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_universe(positions, masses=None):
    """Build a minimal Universe with atoms placed at exact coordinates."""
    n = len(positions)
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n,
        residue_segindex=[0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", [f"A{i}" for i in range(n)])
    u.add_TopologyAttr("resnames", ["LIG"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    u.add_TopologyAttr("masses", [1.0] * n if masses is None else list(masses))
    u.atoms.positions = np.array(positions, dtype=np.float32)
    return u


@pytest.fixture(scope="module")
def two_atom_universe():
    """Two atoms separated by exactly 5.0 Å on the x-axis."""
    return _make_universe([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])


@pytest.fixture(scope="module")
def four_atom_universe():
    """Four atoms on the x-axis: 0, 1, 3, 5 Å."""
    return _make_universe([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
    ])


@pytest.fixture(scope="module")
def weighted_universe():
    """Two atoms with unequal masses to distinguish COM from COG."""
    # atom0 at (0,0,0) mass=3, atom1 at (4,0,0) mass=1
    # COG = (2,0,0) → dist = 2.0  (midpoint)
    # COM = (1,0,0) and (4,0,0) → dist between centers = 3.0
    # Actually: COM of (0,0,0)[m=3] + (4,0,0)[m=1] = (1,0,0)
    # But we need two SEPARATE atom groups to compute distance between them.
    # ag1 = atom0 at (0,0,0) m=3, ag2 = atom1 at (4,0,0) m=1
    # COM of ag1 = (0,0,0), COM of ag2 = (4,0,0), dist = 4
    # Both COM and COG are the same for single atoms.
    # Use 3-atom group: ag1 = atoms at (0,0,0)[m=3] and (2,0,0)[m=1]
    # COM of ag1 = (0*3 + 2*1)/4 = (0.5, 0, 0)
    # COG of ag1 = (0+2)/2 = (1, 0, 0)
    # ag2 = atom at (6, 0, 0)
    # COM dist ≈ 5.5, COG dist = 5.0
    return _make_universe(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [6.0, 0.0, 0.0]],
        masses=[3.0, 1.0, 1.0],
    )


# ---------------------------------------------------------------------------
# Return format
# ---------------------------------------------------------------------------

class TestReturnFormat:
    def test_returns_list_of_tuples(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["com"], ["pair1"], box=None
        )
        assert isinstance(result, list)
        assert len(result) == 1
        label, method, dist = result[0]
        assert label == "pair1"
        assert method == "com"
        assert isinstance(dist, float)

    def test_label_preserved_in_output(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["com"], ["my_label"], box=None
        )
        assert result[0][0] == "my_label"

    def test_method_normalized_to_lowercase(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["COM"], ["d"], box=None
        )
        assert result[0][1] == "com"

    def test_multiple_pairs(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]], u.atoms[[0]]],
            [u.atoms[[1]], u.atoms[[1]]],
            ["com", "cog"],
            ["p1", "p2"],
            box=None,
        )
        assert len(result) == 2
        assert result[0][0] == "p1"
        assert result[1][0] == "p2"

    def test_empty_input_returns_empty_list(self):
        result = calculate_frame_distances([], [], [], [], box=None)
        assert result == []


# ---------------------------------------------------------------------------
# Distance methods
# ---------------------------------------------------------------------------

class TestDistanceMethods:
    def test_com_single_atoms(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["com"], ["d"], box=None
        )
        np.testing.assert_allclose(result[0][2], 5.0, atol=1e-4)

    def test_cog_single_atoms(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["cog"], ["d"], box=None
        )
        np.testing.assert_allclose(result[0][2], 5.0, atol=1e-4)

    def test_com_equals_cog_for_equal_masses(self, two_atom_universe):
        u = two_atom_universe
        r_com = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["com"], ["d"], box=None
        )[0][2]
        r_cog = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["cog"], ["d"], box=None
        )[0][2]
        np.testing.assert_allclose(r_com, r_cog, atol=1e-6)

    def test_com_differs_from_cog_unequal_masses(self, weighted_universe):
        u = weighted_universe
        ag1 = u.atoms[[0, 1]]   # (0,0,0)[m=3] + (2,0,0)[m=1]
        ag2 = u.atoms[[2]]      # (6,0,0)
        r_com = calculate_frame_distances([ag1], [ag2], ["com"], ["d"], box=None)[0][2]
        r_cog = calculate_frame_distances([ag1], [ag2], ["cog"], ["d"], box=None)[0][2]
        # COM of ag1 = (0.5,0,0) → dist to (6,0,0) = 5.5
        # COG of ag1 = (1.0,0,0) → dist to (6,0,0) = 5.0
        assert abs(r_com - r_cog) > 0.1

    def test_min_distance_multi_atom_groups(self, four_atom_universe):
        u = four_atom_universe
        # ag1: (0,0,0) and (1,0,0); ag2: (3,0,0) and (5,0,0)
        # min dist = |(1,0,0)-(3,0,0)| = 2.0
        result = calculate_frame_distances(
            [u.atoms[[0, 1]]], [u.atoms[[2, 3]]], ["min"], ["d"], box=None
        )
        np.testing.assert_allclose(result[0][2], 2.0, atol=1e-4)

    def test_max_distance_multi_atom_groups(self, four_atom_universe):
        u = four_atom_universe
        # max dist = |(0,0,0)-(5,0,0)| = 5.0
        result = calculate_frame_distances(
            [u.atoms[[0, 1]]], [u.atoms[[2, 3]]], ["max"], ["d"], box=None
        )
        np.testing.assert_allclose(result[0][2], 5.0, atol=1e-4)

    def test_method_case_insensitive(self, two_atom_universe):
        u = two_atom_universe
        r_lower = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["com"], ["d"], box=None
        )[0][2]
        r_upper = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["COM"], ["d"], box=None
        )[0][2]
        np.testing.assert_allclose(r_lower, r_upper, atol=1e-6)

    def test_method_with_leading_trailing_whitespace(self, two_atom_universe):
        u = two_atom_universe
        result = calculate_frame_distances(
            [u.atoms[[0]]], [u.atoms[[1]]], ["  com  "], ["d"], box=None
        )
        np.testing.assert_allclose(result[0][2], 5.0, atol=1e-4)


# ---------------------------------------------------------------------------
# Validation / error paths
# ---------------------------------------------------------------------------

class TestDistancesValidation:
    def test_unequal_group1_length_raises(self, two_atom_universe):
        u = two_atom_universe
        with pytest.raises(ValueError, match="equal length"):
            calculate_frame_distances(
                [u.atoms[[0]], u.atoms[[0]]],
                [u.atoms[[1]]],
                ["com", "com"],
                ["a", "b"],
                box=None,
            )

    def test_unequal_methods_length_raises(self, two_atom_universe):
        u = two_atom_universe
        with pytest.raises(ValueError, match="equal length"):
            calculate_frame_distances(
                [u.atoms[[0]]],
                [u.atoms[[1]]],
                ["com", "cog"],  # two methods, one pair
                ["a"],
                box=None,
            )

    def test_empty_group1_raises_runtime_error(self, two_atom_universe):
        u = two_atom_universe
        empty = u.select_atoms("resname ZZZZZ")
        assert empty.n_atoms == 0
        with pytest.raises(RuntimeError, match="EMPTY GROUP1"):
            calculate_frame_distances(
                [empty], [u.atoms[[1]]], ["com"], ["d"], box=None
            )

    def test_empty_group2_raises_runtime_error(self, two_atom_universe):
        u = two_atom_universe
        empty = u.select_atoms("resname ZZZZZ")
        with pytest.raises(RuntimeError, match="EMPTY GROUP2"):
            calculate_frame_distances(
                [u.atoms[[0]]], [empty], ["com"], ["d"], box=None
            )

    def test_unknown_method_raises_runtime_error(self, two_atom_universe):
        u = two_atom_universe
        with pytest.raises(RuntimeError, match="Unknown distance method"):
            calculate_frame_distances(
                [u.atoms[[0]]], [u.atoms[[1]]], ["euclidean"], ["d"], box=None
            )
