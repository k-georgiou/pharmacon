"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.analyzer.properties — RDKit-based molecular descriptors.

All tests use benzene (SMILES: c1ccccc1) as the reference molecule because its
properties are well-known and stable across RDKit versions.
"""
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacon.analyzer.properties import (
    get_element_counts,
    count_total_atoms,
    get_fragment_counts,
    get_molecular_volume,
    get_molecular_weight,
    get_num_rotatable_bonds,
    get_number_of_rings,
    get_tpsa,
    get_number_of_aromatic_rings,
    generate_morgan_fingerprint,
    generate_topological_torsion_fingerprint,
    generate_atom_pair_fingerprint,
    generate_maccs_keys,
    get_stereo_centers,
    get_net_charge,
    get_logp,
    fingerprint_to_string,
)


@pytest.fixture(scope="module")
def benzene():
    """Benzene with explicit hydrogens — c1ccccc1, MW ≈ 78.11."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    return Chem.AddHs(mol)


@pytest.fixture(scope="module")
def ethanol():
    """Ethanol (CCO) — one OH, one rotatable bond."""
    mol = Chem.MolFromSmiles("CCO")
    return Chem.AddHs(mol)


@pytest.fixture(scope="module")
def alanine():
    """L-Alanine — one stereocentre, zwitterion ignored (neutral form)."""
    mol = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    return Chem.AddHs(mol)


class TestElementCounts:
    def test_benzene_has_only_C_and_H(self, benzene):
        counts = get_element_counts(benzene)
        assert set(counts.keys()) == {"C", "H"}

    def test_benzene_has_6_carbons(self, benzene):
        assert get_element_counts(benzene)["C"] == 6

    def test_benzene_has_6_hydrogens(self, benzene):
        assert get_element_counts(benzene)["H"] == 6

    def test_ethanol_contains_oxygen(self, ethanol):
        counts = get_element_counts(ethanol)
        assert "O" in counts

    def test_returns_dict(self, benzene):
        assert isinstance(get_element_counts(benzene), dict)


class TestCountTotalAtoms:
    def test_benzene_12_atoms(self, benzene):
        assert count_total_atoms(benzene) == 12

    def test_returns_int(self, benzene):
        assert isinstance(count_total_atoms(benzene), int)


class TestMolecularWeight:
    def test_benzene_weight_approx_78(self, benzene):
        mw = get_molecular_weight(benzene)
        assert abs(mw - 78.11) < 0.5

    def test_returns_float(self, benzene):
        assert isinstance(get_molecular_weight(benzene), float)

    def test_ethanol_heavier_than_water(self, ethanol):
        assert get_molecular_weight(ethanol) > 18.0


class TestLogP:
    def test_benzene_logp_positive(self, benzene):
        assert get_logp(benzene) > 0

    def test_ethanol_logp_less_than_benzene(self, benzene, ethanol):
        assert get_logp(ethanol) < get_logp(benzene)

    def test_returns_float(self, benzene):
        assert isinstance(get_logp(benzene), float)


class TestTPSA:
    def test_benzene_tpsa_zero(self, benzene):
        assert get_tpsa(benzene) == 0.0

    def test_ethanol_tpsa_positive(self, ethanol):
        assert get_tpsa(ethanol) > 0

    def test_returns_float(self, benzene):
        assert isinstance(get_tpsa(benzene), float)


class TestRings:
    def test_benzene_one_ring(self, benzene):
        assert get_number_of_rings(benzene) == 1

    def test_ethanol_no_rings(self, ethanol):
        assert get_number_of_rings(ethanol) == 0

    def test_benzene_one_aromatic_ring(self, benzene):
        assert get_number_of_aromatic_rings(benzene) == 1

    def test_ethanol_no_aromatic_rings(self, ethanol):
        assert get_number_of_aromatic_rings(ethanol) == 0


class TestRotatableBonds:
    def test_benzene_no_rotatable_bonds(self, benzene):
        assert get_num_rotatable_bonds(benzene) == 0

    def test_ethanol_one_rotatable_bond(self, ethanol):
        assert get_num_rotatable_bonds(ethanol) >= 1

    def test_returns_int(self, benzene):
        assert isinstance(get_num_rotatable_bonds(benzene), int)


class TestStereoAndCharge:
    def test_benzene_no_stereocentres(self, benzene):
        assert get_stereo_centers(benzene) == 0

    def test_alanine_one_stereocentre(self, alanine):
        assert get_stereo_centers(alanine) == 1

    def test_benzene_net_charge_zero(self, benzene):
        assert get_net_charge(benzene) == 0

    def test_returns_int_for_stereo(self, benzene):
        assert isinstance(get_stereo_centers(benzene), int)


class TestMolecularVolume:
    def test_benzene_volume_positive(self, benzene):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        vol = get_molecular_volume(mol)
        assert vol > 0

    def test_returns_float(self, benzene):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        assert isinstance(get_molecular_volume(mol), float)


class TestFingerprints:
    def test_morgan_fingerprint_not_none(self, benzene):
        fp = generate_morgan_fingerprint(benzene)
        assert fp is not None

    def test_topological_torsion_not_none(self, benzene):
        fp = generate_topological_torsion_fingerprint(benzene)
        assert fp is not None

    def test_atom_pair_not_none(self, benzene):
        fp = generate_atom_pair_fingerprint(benzene)
        assert fp is not None

    def test_maccs_keys_not_none(self, benzene):
        fp = generate_maccs_keys(benzene)
        assert fp is not None

    def test_fingerprint_to_string_returns_string(self, benzene):
        fp = generate_morgan_fingerprint(benzene)
        s = fingerprint_to_string(fp)
        assert isinstance(s, str)

    def test_fingerprint_to_string_non_empty(self, benzene):
        fp = generate_morgan_fingerprint(benzene)
        s = fingerprint_to_string(fp)
        assert len(s) > 0

    def test_different_molecules_different_morgan_fps(self, benzene, ethanol):
        fp_b = fingerprint_to_string(generate_morgan_fingerprint(benzene))
        fp_e = fingerprint_to_string(generate_morgan_fingerprint(ethanol))
        assert fp_b != fp_e


class TestFragmentCounts:
    def test_returns_dict(self, benzene):
        result = get_fragment_counts(benzene)
        assert isinstance(result, dict)

    def test_values_are_ints(self, benzene):
        for v in get_fragment_counts(benzene).values():
            assert isinstance(v, int)
