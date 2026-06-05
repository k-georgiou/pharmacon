"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Aggressive test suite for pharmacon.analyzer.interactions

Tests cover every public and private function with emphasis on:
- Exact tuple schema validation (index-by-index)
- Return type correctness
- Edge cases (empty groups, zero contacts, self-interactions)
- Geometric criteria (distances, angles)
- Deduplication logic
- extract_mode_key correctness
- Custom MDAnalysis universes with controlled geometry
"""

import numpy as np
import pytest
import MDAnalysis as Mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysisTests.datafiles import TPR as TPR_FILE, XTC as XTC_FILE

from pharmacon.analyzer.interactions import (
    _get_chain,
    _get_segid,
    _classify_type,
    _is_hydrogen,
    _heavy_neighbors,
    _ring_normal_cross,
    _ring_meta_str,
    calculate_hydrophobic_contacts,
    calculate_hydrogen_bond_contacts,
    calculate_first_degree_water_bridge_contacts,
    calculate_ionic_contacts,
    calculate_halogen_bond_contacts,
    calculate_pi_stacking_contacts,
    calculate_pi_cation_contacts,
    calculate_metal_contacts,
    deduplicate_interactions,
    extract_mode_key,
    interactions_process_frame,
    hbonds_process_frame,
)


@pytest.fixture(scope="module")
def tpr_universe():
    """Universe from the real TPR topology (no trajectory)."""
    return Mda.Universe(str(TPR_FILE))


@pytest.fixture(scope="module")
def tpr_xtc_universe():
    """Universe from the real TPR + XTC (with trajectory)."""
    return Mda.Universe(str(TPR_FILE), str(XTC_FILE))


@pytest.fixture(scope="module")
def two_residue_universe():
    """
    Two-residue universe: ALA1 (N, CA, C, O, CB, HN) and ALA2 (N, CA, C, O, CB, HN).
    Residues placed 3.0 Å apart along X so contacts are predictable.
    Bonds are added for H-bond angle calculations.
    """
    n_atoms = 12
    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=2,
        n_segments=1,
        atom_resindex=[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "CA", "C", "O", "CB", "HN",
                                  "N", "CA", "C", "O", "CB", "HN"])
    u.add_TopologyAttr("elements", ["N", "C", "C", "O", "C", "H",
                                     "N", "C", "C", "O", "C", "H"])
    u.add_TopologyAttr("resnames", ["ALA", "ALA"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["PROA"])
    u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

    # Positions: residue 1 centred at x=0, residue 2 at x=3.0
    positions = np.array([
        [0.0, 0.0, 0.0],    # N(1)
        [1.5, 0.0, 0.0],    # CA(1)
        [2.0, 1.0, 0.0],    # C(1)
        [2.0, 2.0, 0.0],    # O(1)
        [1.5, -1.0, 0.0],   # CB(1)
        [0.0, 1.0, 0.0],    # HN(1) — bonded to N(1)

        [5.0, 0.0, 0.0],    # N(2)
        [6.5, 0.0, 0.0],    # CA(2)
        [7.0, 1.0, 0.0],    # C(2)
        [7.0, 2.0, 0.0],    # O(2)
        [6.5, -1.0, 0.0],   # CB(2)
        [5.0, 1.0, 0.0],    # HN(2) — bonded to N(2)
    ], dtype=np.float32)
    u.atoms.positions = positions

    # Bonds:  N-HN, N-CA, CA-C, C-O, CA-CB  for each residue
    bonds = [
        (0, 5), (0, 1), (1, 2), (2, 3), (1, 4),   # Res 1
        (6, 11), (6, 7), (7, 8), (8, 9), (7, 10),  # Res 2
    ]
    u.add_TopologyAttr("bonds", bonds)

    return u


@pytest.fixture(scope="module")
def hbond_universe():
    """
    Two residues placed so that a hydrogen bond is geometrically valid:
      ALA1-N(donor)-HN ··· O(acceptor)-ALA2

    HN is 2.0 Å from O2, D-H-A angle ~180°, and H-A-X angle ~120°.
    """
    n_atoms = 12
    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=2,
        n_segments=1,
        atom_resindex=[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "CA", "C", "O", "CB", "HN",
                                  "N", "CA", "C", "O", "CB", "HN"])
    u.add_TopologyAttr("elements", ["N", "C", "C", "O", "C", "H",
                                     "N", "C", "C", "O", "C", "H"])
    u.add_TopologyAttr("resnames", ["ALA", "ALA"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["PROA"])
    u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

    # Geometry: donor N(1) at origin, HN pointing toward O(2)
    # N1(0,0,0) — HN1(1,0,0) ··· O2(3,0,0) — C2(4,0.5,0)
    positions = np.array([
        [0.0, 0.0, 0.0],    # N(1) donor
        [0.0, 1.5, 0.0],    # CA(1)
        [0.0, 2.5, 0.0],    # C(1)
        [0.0, 3.5, 0.0],    # O(1)
        [0.0, -1.0, 0.0],   # CB(1)
        [1.0, 0.0, 0.0],    # HN(1) -> pointing at O(2)

        [8.0, 0.0, 0.0],    # N(2)
        [4.5, 0.0, 0.0],    # CA(2)
        [4.0, 0.5, 0.0],    # C(2) -> bonded to O(2)
        [3.0, 0.0, 0.0],    # O(2) acceptor, 2.0 Å from HN(1)
        [5.0, -1.0, 0.0],   # CB(2)
        [8.0, 1.0, 0.0],    # HN(2)
    ], dtype=np.float32)
    u.atoms.positions = positions

    bonds = [
        (0, 5), (0, 1), (1, 2), (2, 3), (1, 4),
        (6, 11), (6, 7), (7, 8), (8, 9), (7, 10),
    ]
    u.add_TopologyAttr("bonds", bonds)

    return u


@pytest.fixture(scope="module")
def water_bridge_universe():
    """
    Three-body system:  Donor(ALA1)-H ··· O(water) ··· H–Acceptor(ALA2)
    Water oxygen has 2 bonded H atoms.
    """
    # ALA1: N, CA, C, O, CB, HN  (indices 0-5)
    # ALA2: N, CA, C, O, CB, HN  (indices 6-11)
    # WAT:  OW, HW1, HW2         (indices 12-14)
    n_atoms = 15
    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=3,
        n_segments=1,
        atom_resindex=[0]*6 + [1]*6 + [2]*3,
        residue_segindex=[0, 0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "CA", "C", "O", "CB", "HN",
                                  "N", "CA", "C", "O", "CB", "HN",
                                  "OW", "HW1", "HW2"])
    u.add_TopologyAttr("elements", ["N", "C", "C", "O", "C", "H",
                                     "N", "C", "C", "O", "C", "H",
                                     "O", "H", "H"])
    u.add_TopologyAttr("resnames", ["ALA", "ALA", "HOH"])
    u.add_TopologyAttr("resids", [1, 2, 3])
    u.add_TopologyAttr("segids", ["PROA"])
    u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

    # Water at centre, ALA1 on left, ALA2 on right.
    # HW1 points left toward O(1), HW2 points right toward O(2).
    positions = np.array([
        # ALA1 (left)
        [-5.0, 0.0, 0.0],   # N(1)
        [-4.0, 1.0, 0.0],   # CA(1)
        [-3.0, 1.5, 0.0],   # C(1)
        [-2.5, 0.0, 0.0],   # O(1)  — acceptor, close to HW1
        [-4.0, -1.0, 0.0],  # CB(1)
        [-5.0, 1.0, 0.0],   # HN(1)

        # ALA2 (right)
        [5.0, 0.0, 0.0],    # N(2)
        [4.0, 1.0, 0.0],    # CA(2)
        [3.0, 1.5, 0.0],    # C(2)
        [2.5, 0.0, 0.0],    # O(2)  — acceptor, close to HW2
        [4.0, -1.0, 0.0],   # CB(2)
        [5.0, 1.0, 0.0],    # HN(2)

        # Water
        [0.0, 0.0, 0.0],    # OW
        [-2.0, 0.0, 0.0],   # HW1 — 2.0 Å from OW, pointing at O(1)
        [2.0, 0.0, 0.0],    # HW2 — 2.0 Å from OW, pointing at O(2)
    ], dtype=np.float32)
    u.atoms.positions = positions

    bonds = [
        (0, 5), (0, 1), (1, 2), (2, 3), (1, 4),
        (6, 11), (6, 7), (7, 8), (8, 9), (7, 10),
        (12, 13), (12, 14),  # water O-H bonds
    ]
    u.add_TopologyAttr("bonds", bonds)

    return u


@pytest.fixture(scope="module")
def ring_universe():
    """
    Two 6-membered aromatic rings (benzene-like) in parallel stacking geometry.
    Ring1 in XY plane at z=0, Ring2 in XY plane at z=3.5 (within face-to-face cutoff).
    """
    n_ring = 6
    n_atoms = 2 * n_ring
    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=2,
        n_segments=1,
        atom_resindex=[0]*n_ring + [1]*n_ring,
        residue_segindex=[0, 0],
        trajectory=True,
    )
    names = [f"C{k}" for k in range(1, n_ring + 1)] * 2
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("elements", ["C"] * n_atoms)
    u.add_TopologyAttr("resnames", ["PHE", "TYR"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["PROA"])
    u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

    # Build two regular hexagons
    angles = np.linspace(0, 2 * np.pi, n_ring, endpoint=False)
    ring1 = np.column_stack([np.cos(angles) * 1.4, np.sin(angles) * 1.4,
                             np.zeros(n_ring)])
    ring2 = ring1.copy()
    ring2[:, 2] = 3.5  # parallel, 3.5 Å apart

    u.atoms.positions = np.vstack([ring1, ring2]).astype(np.float32)

    return u


@pytest.fixture(scope="module")
def metal_universe():
    """
    A metal ion (ZN) near an acceptor oxygen.
    ZN at origin, O of ALA at 2.5 Å away.
    """
    n_atoms = 8
    u = Mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=2,
        n_segments=1,
        atom_resindex=[0, 1, 1, 1, 1, 1, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["ZN", "N", "CA", "C", "O", "CB", "HN", "H"])
    u.add_TopologyAttr("elements", ["Zn", "N", "C", "C", "O", "C", "H", "H"])
    u.add_TopologyAttr("resnames", ["ZN", "ALA"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["ION"])
    u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

    positions = np.array([
        [0.0, 0.0, 0.0],    # ZN
        [5.0, 0.0, 0.0],    # N
        [4.0, 1.0, 0.0],    # CA
        [3.5, 1.5, 0.0],    # C
        [2.5, 0.0, 0.0],    # O — 2.5 Å from ZN
        [4.0, -1.0, 0.0],   # CB
        [5.0, 1.0, 0.0],    # HN
        [5.5, 0.5, 0.0],    # H
    ], dtype=np.float32)
    u.atoms.positions = positions

    bonds = [(1, 6), (1, 2), (2, 3), (3, 4), (2, 5), (1, 7)]
    u.add_TopologyAttr("bonds", bonds)

    return u


class TestGetChain:
    def test_returns_chain_id(self, tpr_universe):
        protein = tpr_universe.select_atoms("protein")
        result = _get_chain(protein, 0)
        assert isinstance(result, str)

    def test_returns_empty_string_on_missing_attr(self):
        """When chainIDs attr is absent, should return empty string."""
        u = Mda.Universe.empty(n_atoms=1, trajectory=True)
        u.add_TopologyAttr("names", ["X"])
        result = _get_chain(u.atoms, 0)
        assert result == ""


class TestGetSegid:
    def test_returns_segid(self, tpr_universe):
        protein = tpr_universe.select_atoms("protein")
        result = _get_segid(protein, 0)
        assert isinstance(result, str)

    def test_returns_empty_string_on_missing_attr(self):
        u = Mda.Universe.empty(n_atoms=1, trajectory=True)
        u.add_TopologyAttr("names", ["X"])
        result = _get_segid(u.atoms, 0)
        assert result == ""


class TestClassifyType:
    def test_water_residue(self):
        assert _classify_type("HOH", "O") == "W"
        assert _classify_type("WAT", "O") == "W"
        assert _classify_type("TIP3", "OH2") == "W"

    def test_backbone_atom(self):
        assert _classify_type("ALA", "N") == "BB"
        assert _classify_type("ALA", "CA") == "BB"
        assert _classify_type("ALA", "C") == "BB"
        assert _classify_type("ALA", "O") == "BB"

    def test_sidechain_atom(self):
        assert _classify_type("ALA", "CB") == "SC"
        assert _classify_type("PHE", "CG") == "SC"

    def test_unknown_residue(self):
        assert _classify_type("LIG", "C1") == "UN"
        assert _classify_type("UNK", "X") == "UN"

    def test_ring_mode_forces_sidechain(self):
        # In ring_mode, even backbone atoms of a protein residue → SC
        assert _classify_type("ALA", "CA", ring_mode=True) == "SC"

    def test_case_insensitivity(self):
        # resname and atom_name are uppercased internally
        assert _classify_type("ala", "ca") == "BB"


class TestIsHydrogen:
    def test_hydrogen_true(self, two_residue_universe):
        hn = two_residue_universe.atoms[5]  # HN of res 1
        assert _is_hydrogen(hn) is True

    def test_heavy_atom_false(self, two_residue_universe):
        n = two_residue_universe.atoms[0]  # N of res 1
        assert _is_hydrogen(n) is False

    def test_missing_element_attr(self):
        """When element attribute is empty, should return False."""

        class FakeAtom:
            element = ""

        assert _is_hydrogen(FakeAtom()) is False


class TestHeavyNeighbors:
    def test_returns_only_heavy(self, two_residue_universe):
        n_atom = two_residue_universe.atoms[0]  # N has bonds to HN and CA
        heavies = _heavy_neighbors(n_atom)
        for nb in heavies:
            assert not _is_hydrogen(nb)

    def test_excludes_hydrogen(self, two_residue_universe):
        n_atom = two_residue_universe.atoms[0]
        all_bonded = list(n_atom.bonded_atoms)
        heavies = _heavy_neighbors(n_atom)
        assert len(heavies) < len(all_bonded)  # HN excluded


class TestRingNormalCross:
    def test_planar_ring_has_z_normal(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        normal = _ring_normal_cross(ring1)
        assert normal.shape == (3,)
        # Hexagon in XY plane → normal along ±Z
        assert abs(abs(normal[2]) - 1.0) < 0.01

    def test_unit_length(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        normal = _ring_normal_cross(ring1)
        assert abs(np.linalg.norm(normal) - 1.0) < 1e-6

    def test_degenerate_returns_default(self):
        """A single atom or collinear atoms → fallback [0,0,1]."""
        u = Mda.Universe.empty(n_atoms=3, trajectory=True)
        u.add_TopologyAttr("names", ["A", "B", "C"])
        u.atoms.positions = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                                     dtype=np.float32)
        normal = _ring_normal_cross(u.atoms)
        assert np.allclose(normal, [0, 0, 1], atol=1e-6)


class TestRingMetaStr:
    def test_returns_9_strings(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        meta = _ring_meta_str(ring1)
        assert len(meta) == 9
        for field in meta:
            assert isinstance(field, str)

    def test_comma_separated_atom_count(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        meta = _ring_meta_str(ring1)
        # 6 atoms → 6 comma-separated tokens in fields that carry data.
        # Some fields (chain, segid) may be empty strings if attribute missing.
        for idx, field in enumerate(meta):
            if field == "":
                continue  # empty chain/segid is valid
            assert len(field.split(",")) == 6, f"Field {idx} has wrong token count: {field!r}"

    def test_indices_field(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        meta = _ring_meta_str(ring1)
        indices = [int(x) for x in meta[0].split(",")]
        expected = list(ring1.indices)
        assert indices == expected

    def test_elements_field(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        meta = _ring_meta_str(ring1)
        elements = meta[3].split(",")
        assert all(e == "C" for e in elements)

    def test_type_uses_ring_mode(self, ring_universe):
        ring1 = ring_universe.atoms[:6]
        meta = _ring_meta_str(ring1)
        types = meta[4].split(",")
        # PHE atoms with ring_mode=True → "SC"
        assert all(t == "SC" for t in types)


class TestHydrophobicContacts:
    def test_returns_tuple_of_tuples(self, two_residue_universe):
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        assert isinstance(result, tuple)
        if len(result) > 0:
            assert isinstance(result[0], tuple)

    def test_empty_group_returns_empty(self, two_residue_universe):
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("name NONEXISTENT")
        assert calculate_hydrophobic_contacts(g1, g2, box=None) == ()

    def test_label_is_hydrophobic(self, two_residue_universe):
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        for rec in result:
            assert rec[0] == "HYDROPHOBIC"

    def test_tuple_length_is_21(self, two_residue_universe):
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        for rec in result:
            assert len(rec) == 21, f"Expected 21 fields, got {len(rec)}: {rec}"

    def test_atom1_fields_types(self, two_residue_universe):
        """Validate index 1-9: atom1 metadata types."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        assert len(result) > 0
        rec = result[0]
        assert isinstance(rec[1], int)          # atom1_index
        assert isinstance(rec[2], str)          # atom1_name
        assert isinstance(rec[3], int)          # atom1_id
        assert isinstance(rec[4], str)          # atom1_element
        assert isinstance(rec[5], str)          # atom1_type
        assert isinstance(rec[6], (str, np.str_))  # atom1_resname
        assert isinstance(rec[7], int)          # atom1_resid
        assert isinstance(rec[8], str)          # atom1_chain
        assert isinstance(rec[9], str)          # atom1_segid

    def test_atom2_fields_types(self, two_residue_universe):
        """Validate index 10-18: atom2 metadata types."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        assert len(result) > 0
        rec = result[0]
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[15], (str, np.str_))
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)

    def test_distance_field(self, two_residue_universe):
        """Index 19 is distance (float), index 20 is bool."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        assert len(result) > 0
        rec = result[0]
        assert isinstance(rec[19], float)
        assert rec[19] > 0.0

    def test_is_hydrogen_involved_field(self, two_residue_universe):
        """Index 20 is is_hydrogen_involved (bool)."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        for rec in result:
            assert isinstance(rec[20], bool)
            # C-C pairs: no hydrogen involved
            assert rec[20] is False

    def test_hydrogen_involved_when_h_present(self, two_residue_universe):
        """Include H atoms → is_hydrogen_involved should be True for H pairs."""
        g1 = two_residue_universe.select_atoms("resid 1 and (element C or element H)")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None,
                                                cutoff_heavy=20.0, cutoff_hydrogen=20.0)
        h_recs = [r for r in result if r[4] == "H" or r[13] == "H"]
        for rec in h_recs:
            assert rec[20] is True

    def test_hh_excluded(self, two_residue_universe):
        """H-H pairs must be excluded."""
        g1 = two_residue_universe.select_atoms("resid 1 and element H")
        g2 = two_residue_universe.select_atoms("resid 2 and element H")
        result = calculate_hydrophobic_contacts(g1, g2, box=None,
                                                cutoff_heavy=20.0, cutoff_hydrogen=20.0)
        assert result == ()

    def test_cutoff_filtering(self, two_residue_universe):
        """Tight cutoff → fewer or zero contacts."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        tight = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=0.1)
        loose = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        assert len(tight) <= len(loose)

    def test_atom_index_matches_universe(self, two_residue_universe):
        """atom1_index (rec[1]) should be a valid universe atom index."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        for rec in result:
            atom = two_residue_universe.atoms[rec[1]]
            assert atom.name == rec[2]
            assert atom.element == rec[4]

    def test_resname_resid_match(self, two_residue_universe):
        """rec[6]/rec[7] should match the atom's actual resname/resid."""
        g1 = two_residue_universe.select_atoms("resid 1 and element C")
        g2 = two_residue_universe.select_atoms("resid 2 and element C")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=20.0)
        for rec in result:
            a1 = two_residue_universe.atoms[rec[1]]
            assert str(a1.resname) == str(rec[6])
            assert int(a1.resid) == rec[7]
            a2 = two_residue_universe.atoms[rec[10]]
            assert str(a2.resname) == str(rec[15])
            assert int(a2.resid) == rec[16]


class TestHydrogenBondContacts:
    def test_returns_tuple_of_tuples(self, hbond_universe):
        g1_acc = hbond_universe.select_atoms("resid 2 and name O")
        g2_don = hbond_universe.select_atoms("resid 1 and name N")
        g1_don = hbond_universe.select_atoms("resid 1 and name N")
        g2_acc = hbond_universe.select_atoms("resid 2 and name O")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, g2_acc, g1_don, g2_don, box=None
        )
        assert isinstance(result, tuple)

    def test_empty_groups_return_empty(self, hbond_universe):
        empty = hbond_universe.select_atoms("name NONEXISTENT")
        result = calculate_hydrogen_bond_contacts(empty, empty, empty, empty, box=None)
        assert result == ()

    def test_finds_hbond_with_loose_criteria(self, hbond_universe):
        """With generous cutoffs, should find the designed H-bond."""
        # O(2) is acceptor at index 9, N(1) is donor at index 0
        g1_acc = hbond_universe.select_atoms("index 9")   # O(2)
        g2_don = hbond_universe.select_atoms("index 0")   # N(1) has bonded HN
        result = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        assert len(result) > 0

    def test_label_is_hydrogen_bond(self, hbond_universe):
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        for rec in result:
            assert rec[0] == "HYDROGEN-BOND"

    def test_tuple_length_is_23(self, hbond_universe):
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        for rec in result:
            assert len(rec) == 23, f"Expected 23 fields, got {len(rec)}"

    def test_schema_field_types(self, hbond_universe):
        """Strict type check on every index."""
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        assert len(result) > 0
        rec = result[0]
        # rec[0] = label
        assert rec[0] == "HYDROGEN-BOND"
        # atom1 block
        assert isinstance(rec[1], int)       # index
        assert isinstance(rec[2], str)       # name
        assert isinstance(rec[3], int)       # id
        assert isinstance(rec[4], str)       # element
        assert isinstance(rec[5], str)       # type
        assert isinstance(rec[7], int)       # resid
        assert isinstance(rec[8], str)       # chain
        assert isinstance(rec[9], str)       # segid
        # atom2 block
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)
        # geometry
        assert isinstance(rec[19], float)    # distance
        assert isinstance(rec[20], float)    # DHA angle
        assert isinstance(rec[21], float)    # HAX angle
        assert isinstance(rec[22], str)      # direction label

    def test_direction_label_values(self, hbond_universe):
        """Direction label must be one of the two defined strings."""
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        g1_don = hbond_universe.select_atoms("index 0")
        g2_acc = hbond_universe.select_atoms("index 9")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, g2_acc, g1_don, g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        valid_labels = {"(G1)-A···H–D-(G2)", "(G1)-D–H···A-(G2)"}
        for rec in result:
            assert rec[22] in valid_labels, f"Unexpected label: {rec[22]}"

    def test_distance_within_cutoff(self, hbond_universe):
        cutoff = 5.0
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        result = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=cutoff
        )
        for rec in result:
            assert rec[19] <= cutoff

    def test_strict_angle_filters_contacts(self, hbond_universe):
        """Very strict angle → fewer/no contacts."""
        g1_acc = hbond_universe.select_atoms("index 9")
        g2_don = hbond_universe.select_atoms("index 0")
        strict = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=179.0, min_angleHAX=179.0, cutoff=5.0
        )
        loose = calculate_hydrogen_bond_contacts(
            g1_acc, hbond_universe.select_atoms("name NONE99"),
            hbond_universe.select_atoms("name NONE99"), g2_don,
            box=None, min_angleDHA=90.0, min_angleHAX=30.0, cutoff=5.0
        )
        assert len(strict) <= len(loose)


class TestWaterBridge1:
    def test_returns_tuple(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None
        )
        assert isinstance(result, tuple)

    def test_empty_water_returns_empty(self, two_residue_universe):
        """No water residues → empty result."""
        u = two_residue_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None
        )
        assert result == ()

    def test_label_is_water_bridge_1(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        for rec in result:
            assert rec[0] == "WATER-BRIDGE-1"

    def test_tuple_length_is_35(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        for rec in result:
            assert len(rec) == 35, f"Expected 35, got {len(rec)}"

    def test_schema_field_types(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        if not result:
            pytest.skip("No water bridges found")
        rec = result[0]
        assert rec[0] == "WATER-BRIDGE-1"
        # atom1 block (indices 1-9)
        assert isinstance(rec[1], int)
        assert isinstance(rec[2], str)
        assert isinstance(rec[3], int)
        assert isinstance(rec[4], str)
        assert isinstance(rec[5], str)
        assert isinstance(rec[7], int)
        assert isinstance(rec[8], str)
        assert isinstance(rec[9], str)
        # atom2 block (indices 10-18)
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)
        # water block (indices 19-27)
        assert isinstance(rec[19], int)
        assert isinstance(rec[20], str)
        assert isinstance(rec[21], int)
        assert isinstance(rec[22], str)   # element
        assert isinstance(rec[23], str)   # type → "W"
        assert isinstance(rec[25], int)   # resid
        assert isinstance(rec[26], str)   # chain
        assert isinstance(rec[27], str)   # segid
        # distances (28-29)
        assert isinstance(rec[28], float)
        assert isinstance(rec[29], float)
        # angles (30-33)
        assert isinstance(rec[30], float)
        assert isinstance(rec[31], float)
        assert isinstance(rec[32], float)
        assert isinstance(rec[33], float)
        # label (34)
        assert isinstance(rec[34], str)

    def test_water_type_is_W(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        for rec in result:
            assert rec[23] == "W", f"Water type should be 'W', got '{rec[23]}'"

    def test_direction_label_values(self, water_bridge_universe):
        u = water_bridge_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        valid = {
            "(G1)-A···H–O–H···A-(G2)",
            "(G1)-A···H–O···H-D-(G2)",
            "(G1)-D-H···O–H···A-(G2)",
            "(G1)-D-H···O···H-D-(G2)",
        }
        for rec in result:
            assert rec[34] in valid, f"Unexpected label: {rec[34]}"


class TestWaterBridge2:
    @pytest.fixture(scope="class")
    def two_water_universe(self):
        """
        Two waters bridging ALA1 and ALA2:
        ALA1-O ··· HW1a-OW1 ··· HW2b-OW2-HW2a ··· O-ALA2
        """
        # ALA1: N, CA, C, O, CB, HN  (0-5)
        # ALA2: N, CA, C, O, CB, HN  (6-11)
        # WAT1: OW1, HW1a, HW1b      (12-14)
        # WAT2: OW2, HW2a, HW2b      (15-17)
        n_atoms = 18
        u = Mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=4,
            n_segments=1,
            atom_resindex=[0]*6 + [1]*6 + [2]*3 + [3]*3,
            residue_segindex=[0, 0, 0, 0],
            trajectory=True,
        )
        u.add_TopologyAttr("names", [
            "N", "CA", "C", "O", "CB", "HN",
            "N", "CA", "C", "O", "CB", "HN",
            "OW", "HW1", "HW2",
            "OW", "HW1", "HW2",
        ])
        u.add_TopologyAttr("elements", [
            "N", "C", "C", "O", "C", "H",
            "N", "C", "C", "O", "C", "H",
            "O", "H", "H",
            "O", "H", "H",
        ])
        u.add_TopologyAttr("resnames", ["ALA", "ALA", "HOH", "HOH"])
        u.add_TopologyAttr("resids", [1, 2, 3, 4])
        u.add_TopologyAttr("segids", ["PROA"])
        u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

        # Linear layout:
        # O(ALA1) at x=-4, OW1 at x=-1, OW2 at x=1, O(ALA2) at x=4
        # HW1a points left (toward O-ALA1), HW1b points right (toward OW2)
        # HW2a points right (toward O-ALA2), HW2b points left (toward OW1)
        positions = np.array([
            [-7.0, 0.0, 0.0],   # N(1)
            [-6.0, 1.0, 0.0],   # CA(1)
            [-5.0, 1.5, 0.0],   # C(1)
            [-4.0, 0.0, 0.0],   # O(1) — acceptor
            [-6.0, -1.0, 0.0],  # CB(1)
            [-7.0, 1.0, 0.0],   # HN(1)

            [7.0, 0.0, 0.0],    # N(2)
            [6.0, 1.0, 0.0],    # CA(2)
            [5.0, 1.5, 0.0],    # C(2)
            [4.0, 0.0, 0.0],    # O(2) — acceptor
            [6.0, -1.0, 0.0],   # CB(2)
            [7.0, 1.0, 0.0],    # HN(2)

            [-1.0, 0.0, 0.0],   # OW1
            [-2.5, 0.0, 0.0],   # HW1a → pointing at O(1), 1.5 Å from OW1
            [0.0, 0.0, 0.0],    # HW1b → pointing right toward OW2, 1.0 Å

            [1.0, 0.0, 0.0],    # OW2
            [2.5, 0.0, 0.0],    # HW2a → pointing at O(2), 1.5 Å from OW2
            [0.0, 0.0, 0.0],    # HW2b → pointing left toward OW1, 1.0 Å
        ], dtype=np.float32)
        u.atoms.positions = positions

        bonds = [
            (0, 5), (0, 1), (1, 2), (2, 3), (1, 4),
            (6, 11), (6, 7), (7, 8), (8, 9), (7, 10),
            (12, 13), (12, 14),  # WAT1
            (15, 16), (15, 17),  # WAT2
        ]
        u.add_TopologyAttr("bonds", bonds)
        return u

    def test_second_degree_false_no_bridge2(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=False,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        assert len(bridge2) == 0

    def test_second_degree_true_finds_bridge2(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        assert len(bridge2) > 0, "Expected second-degree water bridges"

    def test_bridge2_label(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        for rec in bridge2:
            assert rec[0] == "WATER-BRIDGE-2"

    def test_bridge2_tuple_length_is_46(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        for rec in bridge2:
            assert len(rec) == 46, f"Expected 46, got {len(rec)}"

    def test_bridge2_schema_field_types(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        if not bridge2:
            pytest.skip("No second-degree bridges found")
        rec = bridge2[0]

        # atom1 block (1-9)
        assert isinstance(rec[1], int)
        assert isinstance(rec[2], str)
        assert isinstance(rec[3], int)
        assert isinstance(rec[4], str)
        assert isinstance(rec[5], str)
        assert isinstance(rec[7], int)
        assert isinstance(rec[8], str)
        assert isinstance(rec[9], str)

        # atom2 block (10-18)
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)

        # water1 block (19-27)
        assert isinstance(rec[19], int)
        assert isinstance(rec[22], str)  # element
        assert isinstance(rec[23], str)  # type
        assert isinstance(rec[25], int)  # resid

        # water2 block (28-36)
        assert isinstance(rec[28], int)
        assert isinstance(rec[31], str)  # element
        assert isinstance(rec[32], str)  # type
        assert isinstance(rec[34], int)  # resid

        # 3 distances (37-39)
        assert isinstance(rec[37], float)
        assert isinstance(rec[38], float)
        assert isinstance(rec[39], float)

        # angles (40-44)
        for i in range(40, 45):
            assert isinstance(rec[i], float)

        # label (45)
        assert isinstance(rec[45], str)

    def test_bridge2_no_h_reuse(self, two_water_universe):
        """H-atom accounting: same H must not be used for two donation links."""
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        # Verify that the two waters in any bridge are different residues
        for rec in bridge2:
            w1_resid = rec[25]
            w2_resid = rec[34]
            assert w1_resid != w2_resid, "Two water molecules must be different"

    def test_bridge2_direction_labels(self, two_water_universe):
        u = two_water_universe
        g1_acc = u.select_atoms("resid 1 and name O")
        g2_acc = u.select_atoms("resid 2 and name O")
        g1_don = u.select_atoms("resid 1 and name N")
        g2_don = u.select_atoms("resid 2 and name N")
        result = calculate_first_degree_water_bridge_contacts(
            u, g1_acc, g2_acc, g1_don, g2_don, box=None,
            second_degree=True,
            min_angleDHA=10.0, min_angleHAX=10.0, cutoff=5.0, cutoff_loose=20.0
        )
        bridge2 = [r for r in result if r[0] == "WATER-BRIDGE-2"]
        valid = {
            "(G1)-A···W1···W2···A-(G2)",
            "(G1)-A···W1···W2···D-(G2)",
            "(G1)-D···W1···W2···A-(G2)",
            "(G1)-D···W1···W2···D-(G2)",
        }
        for rec in bridge2:
            assert rec[45] in valid, f"Unexpected label: {rec[45]}"


class TestIonicContacts:
    @pytest.fixture(scope="class")
    def ionic_universe(self):
        """Na+ and Cl- 4.0 Å apart."""
        u = Mda.Universe.empty(
            n_atoms=2, n_residues=2, n_segments=1,
            atom_resindex=[0, 1],
            residue_segindex=[0, 0],
            trajectory=True,
        )
        u.add_TopologyAttr("names", ["NA", "CL"])
        u.add_TopologyAttr("elements", ["Na", "Cl"])
        u.add_TopologyAttr("resnames", ["NA", "CL"])
        u.add_TopologyAttr("resids", [1, 2])
        u.add_TopologyAttr("segids", ["ION"])
        u.add_TopologyAttr("ids", [1, 2])
        u.atoms.positions = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
                                     dtype=np.float32)
        return u

    def test_returns_tuple(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None)
        assert isinstance(result, tuple)

    def test_finds_contact_within_cutoff(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        assert len(result) > 0

    def test_no_contact_outside_cutoff(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=3.0)
        assert len(result) == 0

    def test_label_is_ionic(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        for rec in result:
            assert rec[0] == "IONIC"

    def test_tuple_length_is_21(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        for rec in result:
            assert len(rec) == 21, f"Expected 21, got {len(rec)}"

    def test_schema_field_types(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        assert len(result) > 0
        rec = result[0]
        assert rec[0] == "IONIC"
        assert isinstance(rec[1], int)
        assert isinstance(rec[2], str)
        assert isinstance(rec[3], int)
        assert isinstance(rec[4], str)
        assert isinstance(rec[5], str)
        assert isinstance(rec[7], int)
        assert isinstance(rec[8], str)
        assert isinstance(rec[9], str)
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)
        assert isinstance(rec[19], float)    # distance
        assert isinstance(rec[20], str)      # label

    def test_direction_label_values(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        valid = {"(G1+)···(G2-)", "(G1-)···(G2+)"}
        for rec in result:
            assert rec[20] in valid

    def test_distance_value(self, ionic_universe):
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None, cutoff=5.0)
        for rec in result:
            assert abs(rec[19] - 4.0) < 0.01

    def test_empty_groups(self, ionic_universe):
        empty = ionic_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(empty, empty, empty, empty, box=None)
        assert result == ()

    def test_both_directions(self, ionic_universe):
        """Both G1+→G2- and G1-→G2+ should be tested."""
        pos = ionic_universe.select_atoms("name NA")
        neg = ionic_universe.select_atoms("name CL")
        # pos as G1+, neg as G2-: should produce "(G1+)···(G2-)"
        # neg as G1-, pos as G2+: should produce "(G1-)···(G2+)"
        result = calculate_ionic_contacts(pos, neg, pos, neg, box=None, cutoff=5.0)
        labels = {rec[20] for rec in result}
        assert "(G1+)···(G2-)" in labels
        assert "(G1-)···(G2+)" in labels


class TestHalogenBondContacts:
    @pytest.fixture(scope="class")
    def halogen_universe(self):
        """
        Halogen donor (Cl bonded to C) near an acceptor (O bonded to C).
        Cl-C at left, O-C at right, 3.0 Å apart.
        Geometry set for ∠C-Cl···O ≈ 180° and ∠Cl···O-C ≈ 120°.
        """
        n_atoms = 4
        u = Mda.Universe.empty(
            n_atoms=n_atoms, n_residues=2, n_segments=1,
            atom_resindex=[0, 0, 1, 1],
            residue_segindex=[0, 0],
            trajectory=True,
        )
        u.add_TopologyAttr("names", ["CL", "C", "O", "C2"])
        u.add_TopologyAttr("elements", ["Cl", "C", "O", "C"])
        u.add_TopologyAttr("resnames", ["LIG", "ALA"])
        u.add_TopologyAttr("resids", [1, 2])
        u.add_TopologyAttr("segids", ["LIG"])
        u.add_TopologyAttr("ids", [1, 2, 3, 4])

        # Cl(0,0,0) - C(-1.7,0,0)  ...  O(3.0,0,0) - C2(4.0, 0.8, 0)
        positions = np.array([
            [0.0, 0.0, 0.0],     # Cl
            [-1.7, 0.0, 0.0],    # C (bonded to Cl)
            [3.0, 0.0, 0.0],     # O (acceptor)
            [4.0, 0.8, 0.0],     # C2 (bonded to O)
        ], dtype=np.float32)
        u.atoms.positions = positions
        u.add_TopologyAttr("bonds", [(0, 1), (2, 3)])
        return u

    def test_returns_tuple(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None
        )
        assert isinstance(result, tuple)

    def test_finds_halogen_bond(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None,
            donor_min_angle_as_donor=100.0, acceptor_min_angle_as_donor=10.0,
            cutoff=5.0
        )
        assert len(result) > 0

    def test_label_is_halogen_bond(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None,
            donor_min_angle_as_donor=100.0, acceptor_min_angle_as_donor=10.0,
            cutoff=5.0
        )
        for rec in result:
            assert rec[0] == "HALOGEN-BOND"

    def test_tuple_length_is_23(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None,
            donor_min_angle_as_donor=100.0, acceptor_min_angle_as_donor=10.0,
            cutoff=5.0
        )
        for rec in result:
            assert len(rec) == 23, f"Expected 23, got {len(rec)}"

    def test_schema_field_types(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None,
            donor_min_angle_as_donor=100.0, acceptor_min_angle_as_donor=10.0,
            cutoff=5.0
        )
        if not result:
            pytest.skip("No halogen bonds found")
        rec = result[0]
        assert rec[0] == "HALOGEN-BOND"
        # atom1 block
        assert isinstance(rec[1], int)
        assert isinstance(rec[2], str)
        assert isinstance(rec[3], int)
        assert isinstance(rec[4], str)
        assert isinstance(rec[5], str)
        assert isinstance(rec[7], int)
        assert isinstance(rec[8], str)
        assert isinstance(rec[9], str)
        # atom2 block
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)
        # geometry
        assert isinstance(rec[19], float)  # distance
        assert isinstance(rec[20], float)  # angle 1
        assert isinstance(rec[21], float)  # angle 2
        assert isinstance(rec[22], str)    # label

    def test_direction_label_values(self, halogen_universe):
        halogens = halogen_universe.select_atoms("element Cl")
        acceptors = halogen_universe.select_atoms("element O")
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, halogens, acceptors, empty, empty, empty, box=None,
            donor_min_angle_as_donor=100.0, acceptor_min_angle_as_donor=10.0,
            cutoff=5.0
        )
        valid = {
            "(G1)-A···X–Donor-(G2)",
            "(G1)-X(acceptor)···H–D-(G2)",
            "(G1)-Donor X–A···(G2)",
            "(G1)-D–H···X(acceptor)-(G2)",
        }
        for rec in result:
            assert rec[22] in valid

    def test_empty_groups(self, halogen_universe):
        empty = halogen_universe.select_atoms("name NONE99")
        result = calculate_halogen_bond_contacts(
            empty, empty, empty, empty, empty, empty, box=None
        )
        assert result == ()


class TestPiStackingContacts:
    def test_returns_tuple(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        assert isinstance(result, tuple)

    def test_parallel_stacking_detected(self, ring_universe):
        """Two parallel rings 3.5 Å apart → face-to-face."""
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        assert len(result) > 0
        # Check it's classified as parallel
        types = [r[21] for r in result]
        assert "parallel" in types

    def test_empty_groups(self):
        result = calculate_pi_stacking_contacts([], [], box=None)
        assert result == ()

    def test_label_is_pi_stacking(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        for rec in result:
            assert rec[0] == "PI-STACKING"

    def test_tuple_length_is_23(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        for rec in result:
            assert len(rec) == 23, f"Expected 23, got {len(rec)}"

    def test_schema_field_types(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        assert len(result) > 0
        rec = result[0]
        assert rec[0] == "PI-STACKING"
        # ring1 meta (9 comma-separated strings, indices 1-9)
        for i in range(1, 10):
            assert isinstance(rec[i], str)
        # ring2 meta (indices 10-18)
        for i in range(10, 19):
            assert isinstance(rec[i], str)
        # geometry
        assert isinstance(rec[19], float)   # distance
        assert isinstance(rec[20], float)   # theta angle
        assert isinstance(rec[21], str)     # stacking_type
        assert isinstance(rec[22], str)     # direction

    def test_stacking_type_values(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None,
                                                dmax_f2f=10.0, dmax_f2e=10.0)
        valid = {"parallel", "T-shaped"}
        for rec in result:
            assert rec[21] in valid

    def test_direction_label(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        for rec in result:
            assert rec[22] == "(G1)-ring···ring-(G2)"

    def test_ring_indices_match_atoms(self, ring_universe):
        ring1 = [ring_universe.atoms[:6]]
        ring2 = [ring_universe.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        assert len(result) > 0
        rec = result[0]
        # rec[1] = comma-separated ring1 indices
        indices_r1 = [int(x) for x in rec[1].split(",")]
        expected_r1 = list(ring_universe.atoms[:6].indices)
        assert indices_r1 == expected_r1

    def test_too_far_rings_no_contact(self, ring_universe):
        """Move ring2 far away."""
        u = ring_universe.copy()
        u.atoms[6:12].positions += [100.0, 0.0, 0.0]
        ring1 = [u.atoms[:6]]
        ring2 = [u.atoms[6:12]]
        result = calculate_pi_stacking_contacts(ring1, ring2, box=None)
        assert result == ()


class TestPiCationContacts:
    @pytest.fixture(scope="class")
    def pi_cation_universe(self):
        """Ring (6 C atoms in XY plane) + cation (Na+) 4.0 Å above ring centre."""
        n_ring = 6
        n_atoms = n_ring + 1
        u = Mda.Universe.empty(
            n_atoms=n_atoms, n_residues=2, n_segments=1,
            atom_resindex=[0]*n_ring + [1],
            residue_segindex=[0, 0],
            trajectory=True,
        )
        names = [f"C{k}" for k in range(1, n_ring + 1)] + ["NA"]
        u.add_TopologyAttr("names", names)
        u.add_TopologyAttr("elements", ["C"]*n_ring + ["Na"])
        u.add_TopologyAttr("resnames", ["PHE", "NA"])
        u.add_TopologyAttr("resids", [1, 2])
        u.add_TopologyAttr("segids", ["PROA"])
        u.add_TopologyAttr("ids", list(range(1, n_atoms + 1)))

        angles = np.linspace(0, 2*np.pi, n_ring, endpoint=False)
        ring = np.column_stack([np.cos(angles)*1.4, np.sin(angles)*1.4,
                                np.zeros(n_ring)])
        cation = np.array([[0.0, 0.0, 4.0]])
        u.atoms.positions = np.vstack([ring, cation]).astype(np.float32)
        return u

    def test_returns_tuple(self, pi_cation_universe):
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        assert isinstance(result, tuple)

    def test_finds_pi_cation(self, pi_cation_universe):
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        assert len(result) > 0

    def test_label_is_pi_cation(self, pi_cation_universe):
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        for rec in result:
            assert rec[0] == "PI-CATION"

    def test_tuple_length_is_23(self, pi_cation_universe):
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        for rec in result:
            assert len(rec) == 23, f"Expected 23, got {len(rec)}"

    def test_schema_field_types(self, pi_cation_universe):
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        assert len(result) > 0
        rec = result[0]
        assert rec[0] == "PI-CATION"
        assert isinstance(rec[1], str)     # g1_role
        assert rec[1] in ("ring", "cation")
        # g1_payload (9 fields, indices 2-10)
        # g2_payload (9 fields, indices 11-19)
        assert isinstance(rec[20], float)  # distance
        assert isinstance(rec[21], float)  # theta
        assert isinstance(rec[22], str)    # label = "cation-above-face"

    def test_g1_role_ring(self, pi_cation_universe):
        """When G1 supplies rings, g1_role should be 'ring'."""
        rings = [pi_cation_universe.atoms[:6]]
        cations = pi_cation_universe.select_atoms("name NA")
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts(rings, [], empty, cations, box=None)
        for rec in result:
            assert rec[1] == "ring"

    def test_empty_groups(self, pi_cation_universe):
        empty = pi_cation_universe.select_atoms("name NONE99")
        result = calculate_pi_cation_contacts([], [], empty, empty, box=None)
        assert result == ()


class TestMetalContacts:
    def test_returns_tuple(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        don = metal_universe.select_atoms("name N")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, don,
                                          box=None)
        assert isinstance(result, tuple)

    def test_finds_metal_contact(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        don = metal_universe.select_atoms("name N")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, don,
                                          box=None, cutoff=3.0)
        assert len(result) > 0

    def test_label_is_metal_contact(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, empty,
                                          box=None, cutoff=3.0)
        for rec in result:
            assert rec[0] == "METAL-CONTACT"

    def test_tuple_length_is_22(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, empty,
                                          box=None, cutoff=3.0)
        for rec in result:
            assert len(rec) == 22, f"Expected 22, got {len(rec)}"

    def test_schema_field_types(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, empty,
                                          box=None, cutoff=3.0)
        assert len(result) > 0
        rec = result[0]
        assert rec[0] == "METAL-CONTACT"
        # metal block (1-9)
        assert isinstance(rec[1], int)
        assert isinstance(rec[2], str)
        assert isinstance(rec[3], int)
        assert isinstance(rec[4], str)
        assert isinstance(rec[5], str)
        assert isinstance(rec[7], int)
        assert isinstance(rec[8], str)
        assert isinstance(rec[9], str)
        # partner block (10-18)
        assert isinstance(rec[10], int)
        assert isinstance(rec[11], str)
        assert isinstance(rec[12], int)
        assert isinstance(rec[13], str)
        assert isinstance(rec[14], str)
        assert isinstance(rec[16], int)
        assert isinstance(rec[17], str)
        assert isinstance(rec[18], str)
        # distance, role, label
        assert isinstance(rec[19], float)
        assert isinstance(rec[20], str)
        assert isinstance(rec[21], str)

    def test_role_values(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        don = metal_universe.select_atoms("name N")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, don,
                                          box=None, cutoff=10.0)
        roles = {rec[20] for rec in result}
        assert roles <= {"acceptor", "donor"}

    def test_direction_label_values(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        don = metal_universe.select_atoms("name N")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, don,
                                          box=None, cutoff=10.0)
        valid = {
            "(G1)-M···A-(G2)", "(G1)-M···D-(G2)",
            "(G1)-A···M-(G2)", "(G1)-D···M-(G2)",
        }
        for rec in result:
            assert rec[21] in valid

    def test_distance_value(self, metal_universe):
        metals = metal_universe.select_atoms("element Zn")
        acc = metal_universe.select_atoms("name O")
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(metals, empty, empty, empty, acc, empty,
                                          box=None, cutoff=3.0)
        for rec in result:
            assert abs(rec[19] - 2.5) < 0.1

    def test_empty_groups(self, metal_universe):
        empty = metal_universe.select_atoms("name NONE99")
        result = calculate_metal_contacts(empty, empty, empty, empty, empty, empty,
                                          box=None)
        assert result == ()


class TestDeduplicateInteractions:
    def _make_atom_rec(self, label, idx1, resname1, resid1, idx2, resname2, resid2):
        """Build a minimal atom-atom record (21 fields)."""
        return (
            label,
            idx1, "N", 1, "N", "BB", resname1, resid1, "A", "PROA",
            idx2, "O", 2, "O", "BB", resname2, resid2, "A", "PROA",
            3.5, "some_label"
        )

    def test_empty_input(self):
        assert deduplicate_interactions(()) == ()

    def test_single_record(self):
        rec = self._make_atom_rec("IONIC", 0, "ALA", 1, 10, "GLU", 2)
        result = deduplicate_interactions((rec,))
        assert len(result) == 1

    def test_removes_same_residue(self):
        """Same resname+resid+chain+segid → dropped."""
        rec = self._make_atom_rec("IONIC", 0, "ALA", 1, 10, "ALA", 1)
        result = deduplicate_interactions((rec,))
        assert len(result) == 0

    def test_removes_reverse_duplicate(self):
        """A-B and B-A → keep only the first.
        Both records must have identical per-atom blocks (just swapped)
        so the sorted key matches."""
        a1_block = (0, "N", 1, "N", "BB", "ALA", 1, "A", "PROA")
        a2_block = (10, "O", 2, "O", "BB", "GLU", 2, "A", "PROA")
        rec1 = ("IONIC", *a1_block, *a2_block, 3.5, "(G1+)···(G2-)")
        rec2 = ("IONIC", *a2_block, *a1_block, 3.5, "(G1-)···(G2+)")
        result = deduplicate_interactions((rec1, rec2))
        assert len(result) == 1

    def test_keeps_different_residues(self):
        rec1 = self._make_atom_rec("IONIC", 0, "ALA", 1, 10, "GLU", 2)
        rec2 = self._make_atom_rec("IONIC", 0, "ALA", 1, 20, "ASP", 3)
        result = deduplicate_interactions((rec1, rec2))
        assert len(result) == 2

    def test_pi_stacking_dedup(self):
        """PI-STACKING records with same ring identity but reversed → dedup."""
        r1_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                   "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        r2_meta = ("6,7,8", "C1,C2,C3", "7,8,9", "C,C,C", "SC,SC,SC",
                   "TYR,TYR,TYR", "2,2,2", "A,A,A", "PROA,PROA,PROA")
        rec1 = ("PI-STACKING", *r1_meta, *r2_meta, 3.5, 10.0, "parallel",
                "(G1)-ring···ring-(G2)")
        rec2 = ("PI-STACKING", *r2_meta, *r1_meta, 3.5, 10.0, "parallel",
                "(G1)-ring···ring-(G2)")
        result = deduplicate_interactions((rec1, rec2))
        assert len(result) == 1

    def test_pi_stacking_same_residue_dropped(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        rec = ("PI-STACKING", *r_meta, *r_meta, 0.0, 0.0, "parallel",
               "(G1)-ring···ring-(G2)")
        result = deduplicate_interactions((rec,))
        assert len(result) == 0

    def test_pi_cation_dedup(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        c_meta = (10, "NZ", 11, "N", "SC", "LYS", 2, "A", "PROA")
        rec1 = ("PI-CATION", "ring", *r_meta, *c_meta, 4.0, 15.0,
                "cation-above-face")
        rec2 = ("PI-CATION", "cation", *c_meta, *r_meta, 4.0, 15.0,
                "cation-above-face")
        result = deduplicate_interactions((rec1, rec2))
        assert len(result) == 1

    def test_preserves_order(self):
        """First seen record is kept."""
        rec1 = self._make_atom_rec("IONIC", 0, "ALA", 1, 10, "GLU", 2)
        rec2 = self._make_atom_rec("IONIC", 10, "GLU", 2, 0, "ALA", 1)
        result = deduplicate_interactions((rec1, rec2))
        assert result[0] is rec1


class TestExtractModeKey:
    def _make_atom_rec(self, label, resname1, resid1, chain1, segid1,
                       resname2, resid2, chain2, segid2):
        return (
            label,
            0, "N", 1, "N", "BB", resname1, resid1, chain1, segid1,
            10, "O", 2, "O", "BB", resname2, resid2, chain2, segid2,
            3.5, "some_label"
        )

    def test_hydrophobic_key(self):
        rec = self._make_atom_rec("HYDROPHOBIC", "ALA", 1, "A", "PROA",
                                  "LEU", 2, "A", "PROA")
        key = extract_mode_key(rec)
        assert key is not None
        assert key == (("ALA", 1, "A", "PROA"), ("LEU", 2, "A", "PROA"), "HYDROPHOBIC")

    def test_hydrogen_bond_key(self):
        rec = self._make_atom_rec("HYDROGEN-BOND", "ALA", 1, "A", "PROA",
                                  "GLU", 5, "B", "PROB")
        key = extract_mode_key(rec)
        assert key == (("ALA", 1, "A", "PROA"), ("GLU", 5, "B", "PROB"), "HYDROGEN-BOND")

    def test_ionic_key(self):
        rec = self._make_atom_rec("IONIC", "LYS", 3, "", "ION",
                                  "ASP", 10, "", "ION")
        key = extract_mode_key(rec)
        assert key[2] == "IONIC"

    def test_halogen_bond_key(self):
        rec = self._make_atom_rec("HALOGEN-BOND", "LIG", 1, "", "",
                                  "ALA", 2, "A", "PROA")
        key = extract_mode_key(rec)
        assert key[2] == "HALOGEN-BOND"

    def test_metal_contact_key(self):
        rec = self._make_atom_rec("METAL-CONTACT", "ZN", 1, "", "ION",
                                  "ALA", 2, "A", "PROA")
        key = extract_mode_key(rec)
        assert key[2] == "METAL-CONTACT"

    def test_water_bridge_key(self):
        # Water bridge record — uses label starting with "WATER-BRIDGE"
        rec = (
            "WATER-BRIDGE-1",
            0, "N", 1, "N", "BB", "ALA", 1, "A", "PROA",
            10, "O", 2, "O", "BB", "GLU", 5, "A", "PROA",
            20, "OW", 3, "O", "W", "HOH", 100, "", "",
            2.5, 2.7, 150.0, 155.0, 100.0, 95.0,
            "(G1)-A···H–O–H···A-(G2)"
        )
        key = extract_mode_key(rec)
        assert key is not None
        assert key[2] == "WATER-BRIDGE-1"
        assert key[0] == ("ALA", 1, "A", "PROA")
        assert key[1] == ("GLU", 5, "A", "PROA")

    def test_water_bridge_2_key(self):
        rec = (
            "WATER-BRIDGE-2",
            0, "N", 1, "N", "BB", "ALA", 1, "A", "PROA",
            10, "O", 2, "O", "BB", "GLU", 5, "A", "PROA",
            20, "OW", 3, "O", "W", "HOH", 100, "", "",
            30, "OW", 4, "O", "W", "HOH", 101, "", "",
            2.5, 2.0, 2.7, 150.0, 160.0, 155.0, 100.0, 95.0,
            "(G1)-A···W1···W2···A-(G2)"
        )
        key = extract_mode_key(rec)
        assert key is not None
        assert key[2] == "WATER-BRIDGE-2"

    def test_pi_stacking_key(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        r2_meta = ("6,7,8", "C1,C2,C3", "7,8,9", "C,C,C", "SC,SC,SC",
                   "TYR,TYR,TYR", "2,2,2", "A,A,A", "PROA,PROA,PROA")
        rec = ("PI-STACKING", *r_meta, *r2_meta, 3.5, 10.0, "parallel",
               "(G1)-ring···ring-(G2)")
        key = extract_mode_key(rec)
        assert key is not None
        assert key[2] == "PI-STACKING"
        assert key[0] == ("PHE", 1, "A", "PROA")
        assert key[1] == ("TYR", 2, "A", "PROA")

    def test_pi_stacking_self_returns_none(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        rec = ("PI-STACKING", *r_meta, *r_meta, 0.0, 0.0, "parallel",
               "(G1)-ring···ring-(G2)")
        key = extract_mode_key(rec)
        assert key is None

    def test_pi_cation_key_ring_first(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        c_meta = (10, "NZ", 11, "N", "SC", "LYS", 2, "A", "PROA")
        rec = ("PI-CATION", "ring", *r_meta, *c_meta, 4.0, 15.0,
               "cation-above-face")
        key = extract_mode_key(rec)
        assert key is not None
        assert key[2] == "PI-CATION"
        assert key[0] == ("PHE", 1, "A", "PROA")
        assert key[1] == ("LYS", 2, "A", "PROA")

    def test_pi_cation_key_cation_first(self):
        r_meta = ("0,1,2", "C1,C2,C3", "1,2,3", "C,C,C", "SC,SC,SC",
                  "PHE,PHE,PHE", "1,1,1", "A,A,A", "PROA,PROA,PROA")
        c_meta = (10, "NZ", 11, "N", "SC", "LYS", 2, "A", "PROA")
        rec = ("PI-CATION", "cation", *c_meta, *r_meta, 4.0, 15.0,
               "cation-above-face")
        key = extract_mode_key(rec)
        assert key is not None
        assert key[2] == "PI-CATION"
        # cation payload is first → key[0] is cation residue
        assert key[0] == ("LYS", 2, "A", "PROA")
        assert key[1] == ("PHE", 1, "A", "PROA")

    def test_unknown_label_returns_none(self):
        rec = ("UNKNOWN-TYPE", 0, 1, 2, 3, 4, 5, 6, 7, 8,
               10, 11, 12, 13, 14, 15, 16, 17, 18, 3.5, "x")
        key = extract_mode_key(rec)
        assert key is None

    def test_debug_mode_uses_logger_trace(self, caplog):
        """debug param is accepted (backwards compat); trace goes through logger."""
        import logging
        rec = (
            "HYDROPHOBIC",
            0, "N", 1, "N", "BB", "ALA", 1, "A", "PROA",
            10, "O", 2, "O", "BB", "LEU", 2, "A", "PROA",
            3.5, False
        )
        with caplog.at_level(5):  # TRACE = 5
            key = extract_mode_key(rec, debug=True)
        assert key is not None
        # logger.trace output captured at level 5
        assert any("extract_mode_key" in msg for msg in caplog.messages)


class TestHbondsProcessFrame:
    def test_returns_list(self, hbond_universe):
        result = hbonds_process_frame(
            box=None,
            hydrogen_bond_acceptor_atoms_group1=hbond_universe.select_atoms("index 9"),
            hydrogen_bond_acceptor_atoms_group2=hbond_universe.select_atoms("name NONE99"),
            hydrogen_bond_donor_atoms_group1=hbond_universe.select_atoms("name NONE99"),
            hydrogen_bond_donor_atoms_group2=hbond_universe.select_atoms("index 0"),
        )
        assert isinstance(result, list)
        assert len(result) == 1
        assert isinstance(result[0], tuple)


class TestRealTrajectory:
    def test_hydrophobic_on_real_system(self, tpr_universe):
        """Hydrophobic contacts between first two protein residues."""
        g1 = tpr_universe.select_atoms("resid 1 and element C")
        g2 = tpr_universe.select_atoms("resid 2 and element C")
        if len(g1) == 0 or len(g2) == 0:
            pytest.skip("No C atoms in first two residues")
        result = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=5.0)
        assert isinstance(result, tuple)
        for rec in result:
            assert len(rec) == 21
            assert rec[0] == "HYDROPHOBIC"

    def test_hbond_on_real_system(self, tpr_universe):
        """H-bond search on small selection to verify schema."""
        acc = tpr_universe.select_atoms("protein and element O and resid 1:5")
        don = tpr_universe.select_atoms("protein and element N and resid 1:5")
        result = calculate_hydrogen_bond_contacts(acc, acc, don, don, box=None)
        assert isinstance(result, tuple)
        for rec in result:
            assert len(rec) == 23
            assert rec[0] == "HYDROGEN-BOND"

    def test_ionic_on_real_system(self, tpr_universe):
        """Ionic contacts — may be empty if no ions near protein."""
        try:
            pos = tpr_universe.select_atoms("name NA or name K")
            neg = tpr_universe.select_atoms("name CL")
        except Exception:
            pytest.skip("No ions in system")
        empty = tpr_universe.select_atoms("name NONE99")
        result = calculate_ionic_contacts(pos, empty, empty, neg, box=None)
        assert isinstance(result, tuple)
        for rec in result:
            assert len(rec) == 21

    def test_deduplicate_on_real_contacts(self, tpr_universe):
        """Run dedup on real contacts — should not crash, result ≤ input."""
        g1 = tpr_universe.select_atoms("resid 1 and element C")
        g2 = tpr_universe.select_atoms("resid 2 and element C")
        contacts = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=10.0)
        deduped = deduplicate_interactions(contacts)
        assert len(deduped) <= len(contacts)
        assert isinstance(deduped, tuple)

    def test_extract_mode_key_on_real_contacts(self, tpr_universe):
        """extract_mode_key should produce valid keys for real records."""
        g1 = tpr_universe.select_atoms("resid 1 and element C")
        g2 = tpr_universe.select_atoms("resid 2 and element C")
        contacts = calculate_hydrophobic_contacts(g1, g2, box=None, cutoff_heavy=10.0)
        for rec in contacts:
            key = extract_mode_key(rec)
            assert key is not None
            assert len(key) == 3
            r1, r2, label = key
            assert len(r1) == 4  # (resname, resid, chain, segid)
            assert len(r2) == 4
            assert label == "HYDROPHOBIC"
