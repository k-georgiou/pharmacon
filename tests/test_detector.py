"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Test suite for pharmacon.analyzer.detector

Covers all public and private functions with strict assertions on return types,
sizes, edge cases, error paths, and mapping direction normalization.
"""

import pytest
import numpy as np
import MDAnalysis as Mda
from MDAnalysis.core.groups import AtomGroup
from rdkit import Chem
from unittest.mock import patch

from pharmacon.analyzer.detector import (
    _detect_aromatic_atoms,
    _detect_metal_atoms,
    _detect_atoms_by_smarts,
    _log_count,
    detect_hydrophobic_atoms,
    detect_hydrogen_bond_acceptor_atoms,
    detect_hydrogen_bond_donor_atoms,
    detect_positive_charge_atoms,
    detect_negative_charge_atoms,
    detect_halogen_atoms,
    detect_aromatic_atoms,
    detect_metal_atoms,
)
from pharmacon.constants.smarts import (
    AROMATIC_PATTERNS,
    METAL_ELEMENTS,
    HYDROPHOBIC_PATTERNS,
    HYDROPHOBIC_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_PATTERNS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_PATTERNS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    POSITIVELY_CHARGED_PATTERNS,
    NEGATIVELY_CHARGED_PATTERNS,
    HALOGEN_PATTERNS,
    HALOGEN_ELEMENTS,
)


def _universe_from_mol(mol: Chem.Mol) -> Mda.Universe:
    """Build a minimal MDA Universe that mirrors an RDKit Mol atom-for-atom."""
    n = mol.GetNumAtoms()
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n,
        residue_segindex=[0],
    )
    names = [atom.GetSymbol() + str(atom.GetIdx()) for atom in mol.GetAtoms()]
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("elements", elements)
    u.add_TopologyAttr("resnames", ["LIG"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    return u


def _identity_mapping(mol: Chem.Mol) -> dict:
    """1:1 RDKit→MDA identity mapping."""
    return {i: i for i in range(mol.GetNumAtoms())}


def _flipped_mapping(mol: Chem.Mol) -> dict:
    """1:1 MDA→RDKit identity mapping (values are RDKit indices)."""
    return {i: i for i in range(mol.GetNumAtoms())}


def _offset_mapping(mol: Chem.Mol, offset: int = 100) -> tuple:
    """RDKit→MDA mapping where MDA indices are offset (non-overlapping ranges)."""
    n = mol.GetNumAtoms()
    u = Mda.Universe.empty(
        n_atoms=n + offset,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * (n + offset),
        residue_segindex=[0],
    )
    names = [f"X{i}" for i in range(n + offset)]
    elements = ["C"] * (n + offset)
    for i, atom in enumerate(mol.GetAtoms()):
        elements[i + offset] = atom.GetSymbol()
        names[i + offset] = atom.GetSymbol() + str(i)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("elements", elements)
    u.add_TopologyAttr("resnames", ["LIG"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    mapping = {i: i + offset for i in range(n)}
    return u, mapping


@pytest.fixture
def benzene_mol():
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def benzene_universe(benzene_mol):
    return _universe_from_mol(benzene_mol)


@pytest.fixture
def benzene_mapping(benzene_mol):
    return _identity_mapping(benzene_mol)


@pytest.fixture
def toluene_mol():
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def toluene_universe(toluene_mol):
    return _universe_from_mol(toluene_mol)


@pytest.fixture
def toluene_mapping(toluene_mol):
    return _identity_mapping(toluene_mol)


@pytest.fixture
def ethanol_mol():
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def ethanol_universe(ethanol_mol):
    return _universe_from_mol(ethanol_mol)


@pytest.fixture
def ethanol_mapping(ethanol_mol):
    return _identity_mapping(ethanol_mol)


@pytest.fixture
def chlorobenzene_mol():
    mol = Chem.MolFromSmiles("Clc1ccccc1")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def chlorobenzene_universe(chlorobenzene_mol):
    return _universe_from_mol(chlorobenzene_mol)


@pytest.fixture
def chlorobenzene_mapping(chlorobenzene_mol):
    return _identity_mapping(chlorobenzene_mol)


@pytest.fixture
def acetic_acid_mol():
    """Deprotonated acetate [CH3][C](=O)[O-] for negative charge detection."""
    mol = Chem.MolFromSmiles("CC(=O)[O-]")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def acetic_acid_universe(acetic_acid_mol):
    return _universe_from_mol(acetic_acid_mol)


@pytest.fixture
def acetic_acid_mapping(acetic_acid_mol):
    return _identity_mapping(acetic_acid_mol)


@pytest.fixture
def methylamine_mol():
    """Protonated methylamine [CH3][NH3+] for positive charge detection."""
    mol = Chem.MolFromSmiles("C[NH3+]")
    mol = Chem.AddHs(mol)
    return mol


@pytest.fixture
def methylamine_universe(methylamine_mol):
    return _universe_from_mol(methylamine_mol)


@pytest.fixture
def methylamine_mapping(methylamine_mol):
    return _identity_mapping(methylamine_mol)



class TestDetectAromaticAtoms:

    def test_benzene_returns_one_ring(self, benzene_universe, benzene_mol, benzene_mapping):
        rings = _detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        assert isinstance(rings, list)
        assert len(rings) == 1
        assert isinstance(rings[0], AtomGroup)
        assert len(rings[0]) == 6

    def test_benzene_ring_elements_all_carbon(self, benzene_universe, benzene_mol, benzene_mapping):
        rings = _detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        for atom in rings[0]:
            assert atom.element == "C"

    def test_toluene_returns_one_ring(self, toluene_universe, toluene_mol, toluene_mapping):
        rings = _detect_aromatic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        assert len(rings) == 1
        assert len(rings[0]) == 6

    def test_no_aromatic_returns_empty(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        rings = _detect_aromatic_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(rings, list)
        assert len(rings) == 0

    def test_naphthalene_returns_two_rings(self):
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        rings = _detect_aromatic_atoms(u, mol, mapping)
        assert len(rings) >= 2

    def test_flipped_mapping_direction(self, benzene_mol):
        u = _universe_from_mol(benzene_mol)
        mapping = _flipped_mapping(benzene_mol)
        rings = _detect_aromatic_atoms(u, benzene_mol, mapping)
        assert len(rings) == 1
        assert len(rings[0]) == 6

    def test_offset_mapping(self, benzene_mol):
        u, mapping = _offset_mapping(benzene_mol, offset=50)
        rings = _detect_aromatic_atoms(u, benzene_mol, mapping)
        assert len(rings) == 1
        assert len(rings[0]) == 6

    def test_invalid_mapping_raises(self, benzene_universe, benzene_mol):
        bad_mapping = {999: 999, 1000: 1000}
        with pytest.raises(RuntimeError, match="mapping"):
            _detect_aromatic_atoms(benzene_universe, benzene_mol, bad_mapping)

    def test_empty_mapping_no_rings(self, benzene_universe, benzene_mol):
        rings = _detect_aromatic_atoms(benzene_universe, benzene_mol, {})
        assert isinstance(rings, list)

    def test_partial_mapping_skips_unmapped(self, benzene_mol):
        u = _universe_from_mol(benzene_mol)
        n = benzene_mol.GetNumAtoms()
        partial = {i: i for i in range(n // 2)}
        rings = _detect_aromatic_atoms(u, benzene_mol, partial)
        assert isinstance(rings, list)

    def test_duplicate_rings_deduplicated(self, benzene_universe, benzene_mol, benzene_mapping):
        rings = _detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        ring_keys = [tuple(sorted(int(a.index) for a in r)) for r in rings]
        assert len(ring_keys) == len(set(ring_keys))

    def test_out_of_bounds_mapping_raises(self, benzene_mol):
        u = _universe_from_mol(benzene_mol)
        n = benzene_mol.GetNumAtoms()
        bad_mapping = {i: i + 9999 for i in range(n)}
        with pytest.raises(RuntimeError, match="out of bounds"):
            _detect_aromatic_atoms(u, benzene_mol, bad_mapping)

    def test_returns_list_of_atomgroups(self, benzene_universe, benzene_mol, benzene_mapping):
        result = _detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        assert isinstance(result, list)
        for item in result:
            assert isinstance(item, AtomGroup)


class TestDetectMetalAtoms:

    def test_no_metals_returns_empty(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_metal_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)
        assert len(ag) == 0

    def test_single_zinc(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = _detect_metal_atoms(u, mol, mapping)
        assert len(ag) == 1
        assert ag[0].element == "Zn"

    def test_multiple_metals(self):
        mol = Chem.MolFromSmiles("[Zn].[Fe].[Ca]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = _detect_metal_atoms(u, mol, mapping)
        assert len(ag) == 3

    def test_metal_with_organic(self):
        mol = Chem.MolFromSmiles("[Zn].CCO")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = _detect_metal_atoms(u, mol, mapping)
        assert len(ag) == 1

    def test_flipped_mapping(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        mapping = _flipped_mapping(mol)
        ag = _detect_metal_atoms(u, mol, mapping)
        assert len(ag) == 1

    def test_offset_mapping(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u, mapping = _offset_mapping(mol, offset=10)
        ag = _detect_metal_atoms(u, mol, mapping)
        assert len(ag) == 1

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        bad_mapping = {999: 999}
        with pytest.raises(RuntimeError, match="mapping"):
            _detect_metal_atoms(ethanol_universe, ethanol_mol, bad_mapping)

    def test_empty_mapping_returns_empty(self, ethanol_universe, ethanol_mol):
        ag = _detect_metal_atoms(ethanol_universe, ethanol_mol, {})
        assert isinstance(ag, AtomGroup)
        assert len(ag) == 0

    def test_out_of_bounds_mapping_raises(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        bad_mapping = {0: 9999}
        with pytest.raises(RuntimeError, match="out of bounds"):
            _detect_metal_atoms(u, mol, bad_mapping)

    def test_unmapped_metal_still_returns(self):
        mol = Chem.MolFromSmiles("[Zn].[Fe]")
        u = _universe_from_mol(mol)
        partial = {0: 0}
        ag = _detect_metal_atoms(u, mol, partial)
        assert len(ag) == 1

    def test_all_metal_elements_recognized(self):
        for sym in METAL_ELEMENTS:
            sym_title = sym.title()
            smiles = f"[{sym_title}]"
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            u = _universe_from_mol(mol)
            mapping = _identity_mapping(mol)
            ag = _detect_metal_atoms(u, mol, mapping)
            assert len(ag) >= 1, f"Metal {sym} not detected"

    def test_returns_atomgroup(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        result = _detect_metal_atoms(u, mol, mapping)
        assert isinstance(result, AtomGroup)


class TestDetectAtomsBySmarts:

    def test_simple_pattern_match(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#8]"],
            allowed_elements=None,
            label="Oxygen",
        )
        assert isinstance(ag, AtomGroup)
        assert len(ag) == 1
        assert ag[0].element == "O"

    def test_no_match_returns_empty(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#35]"],
            allowed_elements=None,
            label="Bromine",
        )
        assert len(ag) == 0

    def test_allowed_elements_filter(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6,#8]"],
            allowed_elements=["O"],
            label="OxygenOnly",
        )
        for atom in ag:
            assert atom.element.upper() == "O"

    def test_allowed_elements_none_no_filter(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6,#8]"],
            allowed_elements=None,
            label="CandO",
        )
        assert len(ag) >= 2

    def test_invalid_smarts_raises(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        with pytest.raises(RuntimeError, match="Invalid SMARTS"):
            _detect_atoms_by_smarts(
                ethanol_universe, ethanol_mol, ethanol_mapping,
                patterns=["[INVALID_SMARTS!!!"],
                allowed_elements=None,
                label="Bad",
            )

    def test_mix_valid_and_invalid_smarts_raises(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        with pytest.raises(RuntimeError, match="Invalid SMARTS"):
            _detect_atoms_by_smarts(
                ethanol_universe, ethanol_mol, ethanol_mapping,
                patterns=["[#8]", "[NOT_VALID!!!"],
                allowed_elements=None,
                label="Mixed",
            )

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError, match="mapping"):
            _detect_atoms_by_smarts(
                ethanol_universe, ethanol_mol, {999: 999},
                patterns=["[#8]"],
                allowed_elements=None,
                label="BadMap",
            )

    def test_empty_mapping_returns_empty(self, ethanol_universe, ethanol_mol):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, {},
            patterns=["[#8]"],
            allowed_elements=None,
            label="EmptyMap",
        )
        assert len(ag) == 0

    def test_empty_patterns_returns_empty(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=[],
            allowed_elements=None,
            label="Empty",
        )
        assert len(ag) == 0

    def test_multiple_patterns_union(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6]", "[#8]"],
            allowed_elements=None,
            label="C_and_O",
        )
        elements = {atom.element for atom in ag}
        assert "C" in elements
        assert "O" in elements

    def test_flipped_mapping(self, ethanol_mol):
        u = _universe_from_mol(ethanol_mol)
        mapping = _flipped_mapping(ethanol_mol)
        ag = _detect_atoms_by_smarts(
            u, ethanol_mol, mapping,
            patterns=["[#8]"],
            allowed_elements=None,
            label="Flipped",
        )
        assert len(ag) == 1

    def test_offset_mapping(self, ethanol_mol):
        u, mapping = _offset_mapping(ethanol_mol, offset=20)
        ag = _detect_atoms_by_smarts(
            u, ethanol_mol, mapping,
            patterns=["[#8]"],
            allowed_elements=None,
            label="Offset",
        )
        assert len(ag) == 1

    def test_out_of_bounds_mapping_raises(self, ethanol_mol):
        u = _universe_from_mol(ethanol_mol)
        n = ethanol_mol.GetNumAtoms()
        bad_mapping = {i: i + 9999 for i in range(n)}
        with pytest.raises((RuntimeError, IndexError)):
            _detect_atoms_by_smarts(
                u, ethanol_mol, bad_mapping,
                patterns=["[#8]"],
                allowed_elements=None,
                label="OOB",
            )

    def test_allowed_elements_case_insensitive(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag_lower = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6,#8]"],
            allowed_elements=["o"],
            label="Lower",
        )
        ag_upper = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6,#8]"],
            allowed_elements=["O"],
            label="Upper",
        )
        assert len(ag_lower) == len(ag_upper)

    def test_allowed_elements_empty_list_filters_all(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = _detect_atoms_by_smarts(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            patterns=["[#6,#8]"],
            allowed_elements=[],
            label="EmptyFilter",
        )
        assert len(ag) == 0

    def test_custom_label_in_error(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError, match="MyLabel"):
            _detect_atoms_by_smarts(
                ethanol_universe, ethanol_mol, {999: 999},
                patterns=["[#8]"],
                allowed_elements=None,
                label="MyLabel",
            )

    def test_unique_indices_no_duplicates(self, benzene_universe, benzene_mol, benzene_mapping):
        ag = _detect_atoms_by_smarts(
            benzene_universe, benzene_mol, benzene_mapping,
            patterns=["[#6]", "[c]"],
            allowed_elements=None,
            label="Dedup",
        )
        indices = [int(a.index) for a in ag]
        assert len(indices) == len(set(indices))


class TestLogCount:

    def test_log_count_does_not_raise(self):
        _log_count("Test label", 42)

    def test_log_count_zero(self):
        _log_count("Zero count", 0)

    def test_log_count_negative(self):
        _log_count("Negative", -1)

    def test_log_count_large(self):
        _log_count("Large", 999999)


class TestDetectHydrophobicAtoms:

    def test_returns_atomgroup(self, toluene_universe, toluene_mol, toluene_mapping):
        ag = detect_hydrophobic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        assert isinstance(ag, AtomGroup)

    def test_toluene_has_hydrophobic(self, toluene_universe, toluene_mol, toluene_mapping):
        ag = detect_hydrophobic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        assert len(ag) > 0

    def test_custom_label(self, toluene_universe, toluene_mol, toluene_mapping):
        ag = detect_hydrophobic_atoms(
            toluene_universe, toluene_mol, toluene_mapping,
            label="Custom Hydrophobic",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, toluene_universe, toluene_mol):
        with pytest.raises(RuntimeError):
            detect_hydrophobic_atoms(toluene_universe, toluene_mol, {999: 999})

    def test_elements_in_allowed_set(self, toluene_universe, toluene_mol, toluene_mapping):
        ag = detect_hydrophobic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        allowed = {e.upper() for e in HYDROPHOBIC_ELEMENTS}
        for atom in ag:
            assert atom.element.upper() in allowed


class TestDetectHydrogenBondAcceptorAtoms:

    def test_returns_atomgroup(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_acceptor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)

    def test_ethanol_oxygen_is_acceptor(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_acceptor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        elements = {atom.element for atom in ag}
        assert "O" in elements

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_acceptor_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom HBA",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_hydrogen_bond_acceptor_atoms(ethanol_universe, ethanol_mol, {999: 999})

    def test_elements_in_allowed_set(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_acceptor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        allowed = {e.upper() for e in HYDROGEN_BOND_ACCEPTOR_ELEMENTS}
        for atom in ag:
            assert atom.element.upper() in allowed


class TestDetectHydrogenBondDonorAtoms:

    def test_returns_atomgroup(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_donor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)

    def test_ethanol_has_donor(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_donor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(ag) > 0

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_donor_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom HBD",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_hydrogen_bond_donor_atoms(ethanol_universe, ethanol_mol, {999: 999})

    def test_elements_in_allowed_set(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_hydrogen_bond_donor_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        allowed = {e.upper() for e in HYDROGEN_BOND_DONOR_ELEMENTS}
        for atom in ag:
            assert atom.element.upper() in allowed


class TestDetectPositiveChargeAtoms:

    def test_returns_atomgroup(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_positive_charge_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)

    def test_no_positive_in_ethanol(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_positive_charge_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(ag) == 0

    def test_methylamine_has_positive(self, methylamine_universe, methylamine_mol, methylamine_mapping):
        ag = detect_positive_charge_atoms(methylamine_universe, methylamine_mol, methylamine_mapping)
        assert len(ag) >= 1

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_positive_charge_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom Positive",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_positive_charge_atoms(ethanol_universe, ethanol_mol, {999: 999})


class TestDetectNegativeChargeAtoms:

    def test_returns_atomgroup(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_negative_charge_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)

    def test_no_negative_in_ethanol(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_negative_charge_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(ag) == 0

    def test_acetate_has_negative(self, acetic_acid_universe, acetic_acid_mol, acetic_acid_mapping):
        ag = detect_negative_charge_atoms(acetic_acid_universe, acetic_acid_mol, acetic_acid_mapping)
        assert len(ag) >= 1

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_negative_charge_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom Negative",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_negative_charge_atoms(ethanol_universe, ethanol_mol, {999: 999})


class TestDetectHalogenAtoms:

    def test_returns_atomgroup(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_halogen_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert isinstance(ag, AtomGroup)

    def test_no_halogens_in_ethanol(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_halogen_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(ag) == 0

    def test_chlorobenzene_has_halogen(self, chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping):
        ag = detect_halogen_atoms(chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping)
        assert len(ag) >= 1

    def test_chlorobenzene_halogen_is_cl(self, chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping):
        ag = detect_halogen_atoms(chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping)
        elements = {atom.element.upper() for atom in ag}
        assert "CL" in elements

    def test_fluorine_detected(self):
        mol = Chem.MolFromSmiles("FC")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = detect_halogen_atoms(u, mol, mapping)
        assert len(ag) >= 1

    def test_bromine_detected(self):
        mol = Chem.MolFromSmiles("BrC")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = detect_halogen_atoms(u, mol, mapping)
        assert len(ag) >= 1

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        ag = detect_halogen_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom Halogen",
        )
        assert isinstance(ag, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_halogen_atoms(ethanol_universe, ethanol_mol, {999: 999})

    def test_elements_in_allowed_set(self, chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping):
        ag = detect_halogen_atoms(chlorobenzene_universe, chlorobenzene_mol, chlorobenzene_mapping)
        allowed = {e.upper() for e in HALOGEN_ELEMENTS}
        for atom in ag:
            assert atom.element.upper() in allowed


class TestDetectAromaticAtomsPublic:

    def test_returns_list(self, benzene_universe, benzene_mol, benzene_mapping):
        result = detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        assert isinstance(result, list)

    def test_benzene_one_ring(self, benzene_universe, benzene_mol, benzene_mapping):
        result = detect_aromatic_atoms(benzene_universe, benzene_mol, benzene_mapping)
        assert len(result) == 1
        assert len(result[0]) == 6

    def test_no_aromatic(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        result = detect_aromatic_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(result) == 0

    def test_custom_label(self, benzene_universe, benzene_mol, benzene_mapping):
        result = detect_aromatic_atoms(
            benzene_universe, benzene_mol, benzene_mapping,
            label="Custom Aromatic",
        )
        assert isinstance(result, list)

    def test_invalid_mapping_raises(self, benzene_universe, benzene_mol):
        with pytest.raises(RuntimeError):
            detect_aromatic_atoms(benzene_universe, benzene_mol, {999: 999})


class TestDetectMetalAtomsPublic:

    def test_returns_atomgroup(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        result = detect_metal_atoms(u, mol, mapping)
        assert isinstance(result, AtomGroup)

    def test_zinc_detected(self):
        mol = Chem.MolFromSmiles("[Zn]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        result = detect_metal_atoms(u, mol, mapping)
        assert len(result) == 1

    def test_no_metals(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        result = detect_metal_atoms(ethanol_universe, ethanol_mol, ethanol_mapping)
        assert len(result) == 0

    def test_custom_label(self, ethanol_universe, ethanol_mol, ethanol_mapping):
        result = detect_metal_atoms(
            ethanol_universe, ethanol_mol, ethanol_mapping,
            label="Custom Metal",
        )
        assert isinstance(result, AtomGroup)

    def test_invalid_mapping_raises(self, ethanol_universe, ethanol_mol):
        with pytest.raises(RuntimeError):
            detect_metal_atoms(ethanol_universe, ethanol_mol, {999: 999})


class TestMappingNormalization:

    def test_rdkit_to_mda_identity(self, benzene_universe, benzene_mol):
        mapping = {i: i for i in range(benzene_mol.GetNumAtoms())}
        rings = _detect_aromatic_atoms(benzene_universe, benzene_mol, mapping)
        assert len(rings) == 1

    def test_mda_to_rdkit_flipped(self, benzene_mol):
        n = benzene_mol.GetNumAtoms()
        u, mapping = _offset_mapping(benzene_mol, offset=50)
        flipped = {v: k for k, v in mapping.items()}
        rings = _detect_aromatic_atoms(u, benzene_mol, flipped)
        assert len(rings) == 1

    def test_string_keys_raise(self, benzene_universe, benzene_mol):
        bad = {"a": 0, "b": 1}
        with pytest.raises((RuntimeError, ValueError)):
            _detect_aromatic_atoms(benzene_universe, benzene_mol, bad)

    def test_negative_index_keys_accepted_as_flipped(self, benzene_universe, benzene_mol):
        bad = {-1: 0, -2: 1}
        result = _detect_aromatic_atoms(benzene_universe, benzene_mol, bad)
        assert isinstance(result, list)


class TestEdgeCases:

    def test_single_atom_molecule_no_aromatic(self):
        mol = Chem.MolFromSmiles("[He]")
        if mol is None:
            pytest.skip("RDKit cannot parse [He]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        rings = _detect_aromatic_atoms(u, mol, mapping)
        assert rings == []

    def test_single_atom_no_hydrophobic(self):
        mol = Chem.MolFromSmiles("[He]")
        if mol is None:
            pytest.skip("RDKit cannot parse [He]")
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        ag = detect_hydrophobic_atoms(u, mol, mapping)
        assert len(ag) == 0

    def test_multi_ring_fused(self):
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        rings = detect_aromatic_atoms(u, mol, mapping)
        assert len(rings) >= 2

    def test_charged_and_aromatic_combined(self):
        mol = Chem.MolFromSmiles("c1cc[nH+]cc1")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        rings = detect_aromatic_atoms(u, mol, mapping)
        pos = detect_positive_charge_atoms(u, mol, mapping)
        assert len(rings) >= 1
        assert len(pos) >= 1

    def test_halogen_and_aromatic(self):
        mol = Chem.MolFromSmiles("Clc1ccccc1")
        mol = Chem.AddHs(mol)
        u = _universe_from_mol(mol)
        mapping = _identity_mapping(mol)
        rings = detect_aromatic_atoms(u, mol, mapping)
        halogens = detect_halogen_atoms(u, mol, mapping)
        assert len(rings) >= 1
        assert len(halogens) >= 1

    def test_all_detectors_on_same_molecule(self, toluene_universe, toluene_mol, toluene_mapping):
        aromatic = detect_aromatic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        metal = detect_metal_atoms(toluene_universe, toluene_mol, toluene_mapping)
        hydrophobic = detect_hydrophobic_atoms(toluene_universe, toluene_mol, toluene_mapping)
        hba = detect_hydrogen_bond_acceptor_atoms(toluene_universe, toluene_mol, toluene_mapping)
        hbd = detect_hydrogen_bond_donor_atoms(toluene_universe, toluene_mol, toluene_mapping)
        pos = detect_positive_charge_atoms(toluene_universe, toluene_mol, toluene_mapping)
        neg = detect_negative_charge_atoms(toluene_universe, toluene_mol, toluene_mapping)
        halogen = detect_halogen_atoms(toluene_universe, toluene_mol, toluene_mapping)

        assert len(aromatic) >= 1
        assert len(metal) == 0
        assert len(hydrophobic) > 0
        assert len(hba) == 0
        assert len(hbd) == 0
        assert len(pos) == 0
        assert len(neg) == 0
        assert len(halogen) == 0
