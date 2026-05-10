"""
Test suite for pharmacon.analyzer.angles

Covers all public and private functions with strict assertions on return values,
edge cases, error paths, and the recently-fixed calculate_frame_angles loop.
"""

import pytest
import numpy as np
import MDAnalysis as Mda
from MDAnalysis.core.groups import AtomGroup

from pharmacon.analyzer.angles import (
    _displacement,
    _ZERO_TOL,
    calculate_angle_between_vectors,
    calculate_angle_3_atoms,
    calculate_dihedral_angle,
    _parse_vector_spec,
    create_mda_angle_selections,
    _ordered_names_from_selection_string,
    _reorder_ag_by_names,
    calculate_frame_angles,
    GeometrySpecError,
    AXIS_VECTORS,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_universe(positions, names=None, resnames=None):
    """Build a Universe with atoms at exact positions."""
    n = len(positions)
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n,
        residue_segindex=[0],
        trajectory=True,
    )
    if names is None:
        names = [f"A{i}" for i in range(n)]
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", resnames or ["LIG"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    u.atoms.positions = np.array(positions, dtype=np.float32)
    return u


# ---------------------------------------------------------------------------
# Fixtures — universes with known geometry
# ---------------------------------------------------------------------------


@pytest.fixture
def right_angle_universe():
    """Three atoms forming a 90-degree angle at the origin.
    A0(1,0,0)  A1(0,0,0)  A2(0,1,0)  -> angle at A1 = 90 deg
    """
    return _make_universe([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])


@pytest.fixture
def linear_universe():
    """Three atoms in a line -> 180-degree angle."""
    return _make_universe([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ])


@pytest.fixture
def four_atom_planar_universe():
    """Four atoms for dihedral = 0 deg (all in the xy-plane, cis configuration)."""
    return _make_universe([
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], names=["C1", "C2", "C3", "C4"])


@pytest.fixture
def four_atom_trans_universe():
    """Four atoms for dihedral = 180 deg (trans configuration)."""
    return _make_universe([
        [-1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], names=["C1", "C2", "C3", "C4"])


@pytest.fixture
def named_3atom_universe():
    """3 atoms with distinct names for selection/reorder tests."""
    return _make_universe(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        names=["CA", "CB", "CG"],
    )


@pytest.fixture
def named_4atom_universe():
    """4 atoms with distinct names for dihedral selection/reorder tests."""
    return _make_universe(
        [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        names=["CA", "CB", "CG", "CD"],
    )


@pytest.fixture
def two_atom_universe():
    """Two atoms along x-axis for vector spec tests."""
    return _make_universe(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        names=["N1", "N2"],
    )


@pytest.fixture
def coincident_universe():
    """Two atoms at the same position -> zero-length vector."""
    return _make_universe([
        [1.0, 2.0, 3.0],
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
    ])


# ---------------------------------------------------------------------------
# _displacement
# ---------------------------------------------------------------------------


class TestDisplacement:

    def test_basic_displacement(self, two_atom_universe):
        a0, a1 = two_atom_universe.atoms
        v = _displacement(a0, a1, None)
        np.testing.assert_allclose(v, [1.0, 0.0, 0.0], atol=1e-6)

    def test_displacement_direction(self, two_atom_universe):
        a0, a1 = two_atom_universe.atoms
        v_fwd = _displacement(a0, a1, None)
        v_rev = _displacement(a1, a0, None)
        np.testing.assert_allclose(v_fwd, -v_rev, atol=1e-6)

    def test_same_atom_zero_vector(self, two_atom_universe):
        a0 = two_atom_universe.atoms[0]
        v = _displacement(a0, a0, None)
        np.testing.assert_allclose(v, [0.0, 0.0, 0.0], atol=1e-12)

    def test_displacement_3d(self):
        u = _make_universe([[1.0, 2.0, 3.0], [4.0, 6.0, 8.0]])
        a0, a1 = u.atoms
        v = _displacement(a0, a1, None)
        np.testing.assert_allclose(v, [3.0, 4.0, 5.0], atol=1e-5)

    def test_displacement_with_box(self):
        u = _make_universe([[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]])
        a0, a1 = u.atoms
        box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], dtype=np.float32)
        v = _displacement(a0, a1, box)
        assert abs(v[0]) < 5.0

    def test_returns_float_array(self, two_atom_universe):
        a0, a1 = two_atom_universe.atoms
        v = _displacement(a0, a1, None)
        assert v.dtype == np.float64 or v.dtype == np.float32


# ---------------------------------------------------------------------------
# calculate_angle_between_vectors
# ---------------------------------------------------------------------------


class TestCalculateAngleBetweenVectors:

    def test_parallel_vectors_zero_degrees(self):
        u = _make_universe([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0], [3.0, 0.0, 0.0],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_angle_between_vectors(a1, a2, a3, a4)
        assert abs(angle - 0.0) < 1e-6

    def test_antiparallel_vectors_180_degrees(self):
        u = _make_universe([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0], [2.0, 0.0, 0.0],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_angle_between_vectors(a1, a2, a3, a4)
        assert abs(angle - 180.0) < 1e-6

    def test_perpendicular_vectors_90_degrees(self):
        u = _make_universe([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_angle_between_vectors(a1, a2, a3, a4)
        assert abs(angle - 90.0) < 1e-5

    def test_45_degree_angle(self):
        u = _make_universe([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], [1.0, 1.0, 0.0],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_angle_between_vectors(a1, a2, a3, a4)
        assert abs(angle - 45.0) < 1e-4

    def test_returns_float(self):
        u = _make_universe([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        ])
        a1, a2, a3, a4 = u.atoms
        result = calculate_angle_between_vectors(a1, a2, a3, a4)
        assert isinstance(result, float)

    def test_zero_vector_raises(self, coincident_universe):
        a0, a1, a2 = coincident_universe.atoms
        with pytest.raises(ValueError, match="Zero-length"):
            calculate_angle_between_vectors(a0, a1, a0, a2)

    def test_range_0_to_180(self):
        for angle_deg in [0, 30, 60, 90, 120, 150, 180]:
            rad = np.radians(angle_deg)
            u = _make_universe([
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0], [np.cos(rad), np.sin(rad), 0.0],
            ])
            a1, a2, a3, a4 = u.atoms
            result = calculate_angle_between_vectors(a1, a2, a3, a4)
            assert abs(result - angle_deg) < 0.1


# ---------------------------------------------------------------------------
# calculate_angle_3_atoms
# ---------------------------------------------------------------------------


class TestCalculateAngle3Atoms:

    def test_right_angle(self, right_angle_universe):
        a1, a2, a3 = right_angle_universe.atoms
        angle = calculate_angle_3_atoms(a1, a2, a3)
        assert abs(angle - 90.0) < 1e-4

    def test_linear_angle(self, linear_universe):
        a1, a2, a3 = linear_universe.atoms
        angle = calculate_angle_3_atoms(a1, a2, a3)
        assert abs(angle - 180.0) < 1e-4

    def test_60_degree_equilateral(self):
        u = _make_universe([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, np.sqrt(3) / 2, 0.0],
        ])
        a1, a2, a3 = u.atoms
        angle = calculate_angle_3_atoms(a1, a2, a3)
        assert abs(angle - 60.0) < 0.1

    def test_returns_float(self, right_angle_universe):
        a1, a2, a3 = right_angle_universe.atoms
        result = calculate_angle_3_atoms(a1, a2, a3)
        assert isinstance(result, float)

    def test_order_matters(self, right_angle_universe):
        a1, a2, a3 = right_angle_universe.atoms
        angle_123 = calculate_angle_3_atoms(a1, a2, a3)
        angle_321 = calculate_angle_3_atoms(a3, a2, a1)
        assert abs(angle_123 - angle_321) < 1e-6

    def test_vertex_is_middle_atom(self):
        u = _make_universe([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        a1, a2, a3 = u.atoms
        angle = calculate_angle_3_atoms(a1, a2, a3)
        assert abs(angle - 90.0) < 1e-4

    def test_with_box(self):
        u = _make_universe([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0], dtype=np.float32)
        a1, a2, a3 = u.atoms
        angle = calculate_angle_3_atoms(a1, a2, a3, box=box)
        assert abs(angle - 90.0) < 1e-3


# ---------------------------------------------------------------------------
# calculate_dihedral_angle
# ---------------------------------------------------------------------------


class TestCalculateDihedralAngle:

    def test_cis_dihedral_zero(self, four_atom_planar_universe):
        a1, a2, a3, a4 = four_atom_planar_universe.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4)
        assert abs(angle) < 1.0 or abs(abs(angle) - 360.0) < 1.0

    def test_trans_dihedral_180(self, four_atom_trans_universe):
        a1, a2, a3, a4 = four_atom_trans_universe.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4)
        assert abs(abs(angle) - 180.0) < 1.0

    def test_90_degree_dihedral(self):
        u = _make_universe([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4)
        assert abs(abs(angle) - 90.0) < 1.0

    def test_returns_float(self, four_atom_planar_universe):
        a1, a2, a3, a4 = four_atom_planar_universe.atoms
        result = calculate_dihedral_angle(a1, a2, a3, a4)
        assert isinstance(result, float)

    def test_dihedral_range(self, four_atom_planar_universe):
        a1, a2, a3, a4 = four_atom_planar_universe.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4)
        assert -180.0 <= angle <= 180.0

    def test_with_box(self, four_atom_planar_universe):
        box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0], dtype=np.float32)
        a1, a2, a3, a4 = four_atom_planar_universe.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4, box=box)
        assert isinstance(angle, float)


# ---------------------------------------------------------------------------
# _parse_vector_spec
# ---------------------------------------------------------------------------


class TestParseVectorSpec:

    def test_x_axis(self, two_atom_universe):
        kind, obj = _parse_vector_spec(two_atom_universe, "x-axis")
        assert kind == "axis"
        np.testing.assert_allclose(obj, [1.0, 0.0, 0.0])

    def test_y_axis(self, two_atom_universe):
        kind, obj = _parse_vector_spec(two_atom_universe, "y-axis")
        assert kind == "axis"
        np.testing.assert_allclose(obj, [0.0, 1.0, 0.0])

    def test_z_axis(self, two_atom_universe):
        kind, obj = _parse_vector_spec(two_atom_universe, "z-axis")
        assert kind == "axis"
        np.testing.assert_allclose(obj, [0.0, 0.0, 1.0])

    def test_axis_case_insensitive(self, two_atom_universe):
        kind, _ = _parse_vector_spec(two_atom_universe, "X-AXIS")
        assert kind == "axis"

    def test_axis_with_whitespace(self, two_atom_universe):
        kind, _ = _parse_vector_spec(two_atom_universe, "  x-axis  ")
        assert kind == "axis"

    def test_atom_selection_two_atoms(self, two_atom_universe):
        kind, obj = _parse_vector_spec(two_atom_universe, "name N1 or name N2")
        assert kind == "atoms"
        assert isinstance(obj, AtomGroup)
        assert obj.n_atoms == 2

    def test_atom_selection_wrong_count_raises(self, right_angle_universe):
        with pytest.raises(GeometrySpecError, match="exactly 2 atoms"):
            _parse_vector_spec(right_angle_universe, "all")

    def test_atom_selection_one_atom_raises(self, two_atom_universe):
        with pytest.raises(GeometrySpecError, match="exactly 2 atoms"):
            _parse_vector_spec(two_atom_universe, "name N1")

    def test_atom_selection_zero_atoms_raises(self, two_atom_universe):
        with pytest.raises(GeometrySpecError, match="exactly 2 atoms"):
            _parse_vector_spec(two_atom_universe, "name NONEXISTENT")


# ---------------------------------------------------------------------------
# create_mda_angle_selections
# ---------------------------------------------------------------------------


class TestCreateMdaAngleSelections:

    def test_three_atoms_returns_angle(self, named_3atom_universe):
        result = create_mda_angle_selections(named_3atom_universe, "name CA or name CB or name CG")
        assert result["type"] == "angle"
        assert result["atoms"].n_atoms == 3
        assert "spec_str" in result

    def test_four_atoms_returns_dihedral(self, named_4atom_universe):
        result = create_mda_angle_selections(named_4atom_universe, "name CA or name CB or name CG or name CD")
        assert result["type"] == "dihedral"
        assert result["atoms"].n_atoms == 4
        assert "spec_str" in result

    def test_vector_angle_axis_to_axis(self, two_atom_universe):
        result = create_mda_angle_selections(two_atom_universe, "x-axis -> y-axis")
        assert result["type"] == "vector-angle"
        assert result["left"][0] == "axis"
        assert result["right"][0] == "axis"

    def test_vector_angle_atoms_to_axis(self, two_atom_universe):
        result = create_mda_angle_selections(two_atom_universe, "name N1 or name N2 -> x-axis")
        assert result["type"] == "vector-angle"
        assert result["left"][0] == "atoms"
        assert result["right"][0] == "axis"

    def test_multiple_arrows_raises(self, two_atom_universe):
        with pytest.raises(GeometrySpecError, match="Only one"):
            create_mda_angle_selections(two_atom_universe, "x-axis -> y-axis -> z-axis")

    def test_two_atoms_raises(self, two_atom_universe):
        with pytest.raises(GeometrySpecError, match="3 atoms.*4 atoms"):
            create_mda_angle_selections(two_atom_universe, "name N1 or name N2")

    def test_five_atoms_raises(self):
        u = _make_universe(
            [[i, 0, 0] for i in range(5)],
            names=["A", "B", "C", "D", "E"],
        )
        with pytest.raises(GeometrySpecError, match="3 atoms.*4 atoms"):
            create_mda_angle_selections(u, "all")

    def test_zero_atoms_raises(self, two_atom_universe):
        with pytest.raises(GeometrySpecError, match="3 atoms.*4 atoms"):
            create_mda_angle_selections(two_atom_universe, "name NONEXISTENT")


# ---------------------------------------------------------------------------
# _ordered_names_from_selection_string
# ---------------------------------------------------------------------------


class TestOrderedNamesFromSelectionString:

    def test_basic_extraction(self):
        s = "(resname LIG and name C14) or (resname LIG and name C09)"
        result = _ordered_names_from_selection_string(s)
        assert result == ["C14", "C09"]

    def test_preserves_order(self):
        s = "name ZZ or name AA or name MM"
        result = _ordered_names_from_selection_string(s)
        assert result == ["ZZ", "AA", "MM"]

    def test_deduplicates(self):
        s = "name CA or name CB or name CA"
        result = _ordered_names_from_selection_string(s)
        assert result == ["CA", "CB"]

    def test_case_uppercased(self):
        s = "name ca or name cb"
        result = _ordered_names_from_selection_string(s)
        assert result == ["CA", "CB"]

    def test_four_names_for_dihedral(self):
        s = "(resname UNA and name C14) or (resname UNA and name C09) or (resname UNA and name N08) or (resname UNA and name N07)"
        result = _ordered_names_from_selection_string(s)
        assert result == ["C14", "C09", "N08", "N07"]

    def test_no_name_tokens_returns_empty(self):
        result = _ordered_names_from_selection_string("resname LIG")
        assert result == []

    def test_empty_string(self):
        result = _ordered_names_from_selection_string("")
        assert result == []

    def test_single_name(self):
        result = _ordered_names_from_selection_string("name CA")
        assert result == ["CA"]

    def test_alphanumeric_names(self):
        result = _ordered_names_from_selection_string("name C1 or name H2_A")
        assert result == ["C1", "H2_A"]


# ---------------------------------------------------------------------------
# _reorder_ag_by_names
# ---------------------------------------------------------------------------


class TestReorderAgByNames:

    def test_reorder_reverses(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        reordered = _reorder_ag_by_names(ag, ["CG", "CB", "CA"])
        names = [a.name for a in reordered]
        assert names == ["CG", "CB", "CA"]

    def test_reorder_same_order(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        reordered = _reorder_ag_by_names(ag, ["CA", "CB", "CG"])
        names = [a.name for a in reordered]
        assert names == ["CA", "CB", "CG"]

    def test_reorder_4_atoms(self, named_4atom_universe):
        ag = named_4atom_universe.atoms
        reordered = _reorder_ag_by_names(ag, ["CD", "CG", "CB", "CA"])
        names = [a.name for a in reordered]
        assert names == ["CD", "CG", "CB", "CA"]

    def test_count_mismatch_raises(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        with pytest.raises(GeometrySpecError, match="Cannot reorder"):
            _reorder_ag_by_names(ag, ["CA", "CB"])

    def test_missing_name_raises(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        with pytest.raises(GeometrySpecError, match="missing"):
            _reorder_ag_by_names(ag, ["CA", "CB", "XX"])

    def test_duplicate_names_raises(self):
        u = _make_universe(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
            names=["CA", "CA", "CB"],
        )
        with pytest.raises(GeometrySpecError, match="duplicate"):
            _reorder_ag_by_names(u.atoms, ["CA", "CA", "CB"])

    def test_returns_atomgroup(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        result = _reorder_ag_by_names(ag, ["CG", "CB", "CA"])
        assert isinstance(result, AtomGroup)

    def test_preserves_positions(self, named_3atom_universe):
        ag = named_3atom_universe.atoms
        original_pos = {a.name: a.position.copy() for a in ag}
        reordered = _reorder_ag_by_names(ag, ["CG", "CB", "CA"])
        for atom in reordered:
            np.testing.assert_allclose(atom.position, original_pos[atom.name], atol=1e-6)


# ---------------------------------------------------------------------------
# calculate_frame_angles
# ---------------------------------------------------------------------------


class TestCalculateFrameAngles:

    def test_single_angle_spec(self, right_angle_universe):
        spec = {"type": "angle", "atoms": right_angle_universe.atoms, "spec_str": ""}
        results = calculate_frame_angles([spec], ["test"], box=None)
        assert len(results) == 1
        label, kind, value = results[0]
        assert label == "test"
        assert kind == "angle"
        assert abs(value - 90.0) < 1e-3

    def test_single_dihedral_spec(self, four_atom_trans_universe):
        spec = {"type": "dihedral", "atoms": four_atom_trans_universe.atoms, "spec_str": ""}
        results = calculate_frame_angles([spec], ["dih"], box=None)
        assert len(results) == 1
        label, kind, value = results[0]
        assert label == "dih"
        assert kind == "dihedral"
        assert abs(abs(value) - 180.0) < 1.0

    def test_vector_angle_axis_to_axis(self):
        spec = {
            "type": "vector-angle",
            "left": ("axis", AXIS_VECTORS["x-axis"]),
            "right": ("axis", AXIS_VECTORS["y-axis"]),
        }
        results = calculate_frame_angles([spec], ["xy"], box=None)
        assert len(results) == 1
        label, kind, value = results[0]
        assert label == "xy"
        assert kind == "vector-angle"
        assert abs(value - 90.0) < 1e-6

    def test_vector_angle_same_axis_zero(self):
        spec = {
            "type": "vector-angle",
            "left": ("axis", AXIS_VECTORS["x-axis"]),
            "right": ("axis", AXIS_VECTORS["x-axis"]),
        }
        results = calculate_frame_angles([spec], ["xx"], box=None)
        assert abs(results[0][2] - 0.0) < 1e-6

    def test_vector_angle_atoms(self, two_atom_universe):
        ag = two_atom_universe.atoms
        spec = {
            "type": "vector-angle",
            "left": ("atoms", ag),
            "right": ("axis", AXIS_VECTORS["x-axis"]),
        }
        results = calculate_frame_angles([spec], ["atoms_vs_x"], box=None)
        assert abs(results[0][2] - 0.0) < 1e-4

    def test_multiple_specs_all_processed(self, right_angle_universe, four_atom_trans_universe):
        spec_angle = {"type": "angle", "atoms": right_angle_universe.atoms, "spec_str": ""}
        spec_dih = {"type": "dihedral", "atoms": four_atom_trans_universe.atoms, "spec_str": ""}
        spec_vec = {
            "type": "vector-angle",
            "left": ("axis", AXIS_VECTORS["x-axis"]),
            "right": ("axis", AXIS_VECTORS["y-axis"]),
        }
        results = calculate_frame_angles(
            [spec_angle, spec_dih, spec_vec],
            ["ang", "dih", "vec"],
            box=None,
        )
        assert len(results) == 3
        assert results[0][0] == "ang"
        assert results[1][0] == "dih"
        assert results[2][0] == "vec"

    def test_multiple_specs_correct_values(self, right_angle_universe, four_atom_trans_universe):
        spec_angle = {"type": "angle", "atoms": right_angle_universe.atoms, "spec_str": ""}
        spec_dih = {"type": "dihedral", "atoms": four_atom_trans_universe.atoms, "spec_str": ""}
        results = calculate_frame_angles(
            [spec_angle, spec_dih],
            ["ang", "dih"],
            box=None,
        )
        assert abs(results[0][2] - 90.0) < 1e-3
        assert abs(abs(results[1][2]) - 180.0) < 1.0

    def test_empty_specs_returns_empty(self):
        results = calculate_frame_angles([], [], box=None)
        assert results == []

    def test_mismatched_specs_labels_raises(self, right_angle_universe):
        spec = {"type": "angle", "atoms": right_angle_universe.atoms, "spec_str": ""}
        with pytest.raises(ValueError, match="same length"):
            calculate_frame_angles([spec], ["a", "b"], box=None)

    def test_unknown_type_raises(self):
        spec = {"type": "unknown-type"}
        with pytest.raises(RuntimeError, match="Unknown geometry type"):
            calculate_frame_angles([spec], ["x"], box=None)

    def test_vector_angle_zero_vector_raises(self):
        spec = {
            "type": "vector-angle",
            "left": ("axis", np.array([0.0, 0.0, 0.0])),
            "right": ("axis", AXIS_VECTORS["x-axis"]),
        }
        with pytest.raises(ValueError, match="Zero-length"):
            calculate_frame_angles([spec], ["zero"], box=None)

    def test_angle_with_reorder(self, named_3atom_universe):
        spec_str = "(name CG) or (name CB) or (name CA)"
        ag = named_3atom_universe.atoms
        spec = {"type": "angle", "atoms": ag, "spec_str": spec_str}
        results = calculate_frame_angles([spec], ["reordered"], box=None)
        assert len(results) == 1
        assert results[0][1] == "angle"

    def test_dihedral_with_reorder(self, named_4atom_universe):
        spec_str = "(name CD) or (name CG) or (name CB) or (name CA)"
        ag = named_4atom_universe.atoms
        spec = {"type": "dihedral", "atoms": ag, "spec_str": spec_str}
        results = calculate_frame_angles([spec], ["reordered_dih"], box=None)
        assert len(results) == 1
        assert results[0][1] == "dihedral"

    def test_with_box(self, right_angle_universe):
        box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0], dtype=np.float32)
        spec = {"type": "angle", "atoms": right_angle_universe.atoms, "spec_str": ""}
        results = calculate_frame_angles([spec], ["boxed"], box=box)
        assert abs(results[0][2] - 90.0) < 1e-2


# ---------------------------------------------------------------------------
# Edge cases & integration
# ---------------------------------------------------------------------------


class TestEdgeCases:

    def test_angle_result_always_0_to_180(self):
        for deg in [10, 45, 90, 135, 170]:
            rad = np.radians(deg)
            u = _make_universe([
                [np.cos(rad), np.sin(rad), 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ])
            a1, a2, a3 = u.atoms
            angle = calculate_angle_3_atoms(a1, a2, a3)
            assert 0.0 <= angle <= 180.0

    def test_dihedral_result_minus180_to_180(self):
        u = _make_universe([
            [1.0, 0.0, 0.5],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, -0.5],
        ])
        a1, a2, a3, a4 = u.atoms
        angle = calculate_dihedral_angle(a1, a2, a3, a4)
        assert -180.0 <= angle <= 180.0

    def test_axis_vectors_are_unit_length(self):
        for name, vec in AXIS_VECTORS.items():
            assert abs(np.linalg.norm(vec) - 1.0) < 1e-12, f"{name} is not unit length"

    def test_full_pipeline_angle(self, named_3atom_universe):
        sel_str = "name CA or name CB or name CG"
        spec = create_mda_angle_selections(named_3atom_universe, sel_str)
        assert spec["type"] == "angle"
        results = calculate_frame_angles([spec], ["pipeline"], box=None)
        assert len(results) == 1
        assert results[0][1] == "angle"
        assert 0.0 <= results[0][2] <= 180.0

    def test_full_pipeline_dihedral(self, named_4atom_universe):
        sel_str = "name CA or name CB or name CG or name CD"
        spec = create_mda_angle_selections(named_4atom_universe, sel_str)
        assert spec["type"] == "dihedral"
        results = calculate_frame_angles([spec], ["pipeline_dih"], box=None)
        assert len(results) == 1
        assert results[0][1] == "dihedral"

    def test_full_pipeline_vector_angle(self, two_atom_universe):
        sel_str = "name N1 or name N2 -> x-axis"
        spec = create_mda_angle_selections(two_atom_universe, sel_str)
        assert spec["type"] == "vector-angle"
        results = calculate_frame_angles([spec], ["vec_pipe"], box=None)
        assert len(results) == 1
        assert abs(results[0][2] - 0.0) < 1e-4

    def test_geometry_spec_error_is_value_error(self):
        assert issubclass(GeometrySpecError, ValueError)
