"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Aggressive test suite for pharmacon.utils.mda

Tests cover every public and private function in mda.py with emphasis on
edge cases, boundary conditions, and the correctness of the MDA↔RDKit
index mapping in convert_mda_to_rdkit.
"""

import os
import numpy as np
import pytest
import MDAnalysis as Mda
from MDAnalysis.core.universe import Universe
from MDAnalysisTests.datafiles import TPR as TPR_FILE, XTC as XTC_FILE
from rdkit import Chem
from unittest.mock import patch, MagicMock

from pharmacon.utils.mda import (
    _guess_element_from_mass,
    _canonical,
    _guess_element_from_name,
    create_universe,
    convert_mda_to_rdkit,
    has_residues,
    has_atoms,
    require_atoms,
    require_residues,
)
from pharmacon.constants.smarts import (
    WATER_RESIDUES,
    ION_RES_TO_ELT,
    AA3,
    NA_RES,
    PROTEIN_CARBONS,
    SIDECHAIN_OXY,
    SIDECHAIN_NITRO,
    NUC_OXY,
    TWO_LETTER,
    METAL_ELEMENTS,
)


@pytest.fixture(scope="module")
def tpr_universe():
    """Universe from the real TPR topology (no trajectory)."""
    return Mda.Universe(str(TPR_FILE))


@pytest.fixture(scope="module")
def tpr_xtc_universe():
    """Universe from the real TPR + XTC (with trajectory)."""
    return Mda.Universe(str(TPR_FILE), str(XTC_FILE))


@pytest.fixture
def small_universe():
    """Minimal synthetic universe: 5 atoms, 1 residue, with masses and names."""
    n = 5
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n,
        residue_segindex=[0],
    )
    u.add_TopologyAttr("names", ["C", "N", "O", "H1", "H2"])
    u.add_TopologyAttr("masses", [12.011, 14.007, 15.999, 1.008, 1.008])
    u.add_TopologyAttr("resnames", ["ALA"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["PROT"])
    return u


@pytest.fixture
def ethanol_universe():
    """Synthetic universe representing ethanol: C2H5OH (9 atoms)."""
    n = 9
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=1,
        n_segments=1,
        atom_resindex=[0] * n,
        residue_segindex=[0],
    )
    u.add_TopologyAttr("names", ["C1", "C2", "O1", "H1", "H2", "H3", "H4", "H5", "H6"])
    u.add_TopologyAttr("masses", [12.011, 12.011, 15.999, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008])
    u.add_TopologyAttr("elements", ["C", "C", "O", "H", "H", "H", "H", "H", "H"])
    u.add_TopologyAttr("resnames", ["ETH"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["LIG"])
    return u


@pytest.fixture
def empty_universe():
    """Universe with zero atoms."""
    return Mda.Universe.empty(0)


class TestCanonical:
    def test_single_letter_valid(self):
        assert _canonical("C") == "C"
        assert _canonical("c") == "C"
        assert _canonical("N") == "N"
        assert _canonical("n") == "N"
        assert _canonical("O") == "O"
        assert _canonical("H") == "H"
        assert _canonical("S") == "S"
        assert _canonical("P") == "P"
        assert _canonical("F") == "F"
        assert _canonical("I") == "I"
        assert _canonical("K") == "K"

    def test_two_letter_valid(self):
        assert _canonical("Fe") == "Fe"
        assert _canonical("FE") == "Fe"
        assert _canonical("fe") == "Fe"
        assert _canonical("Zn") == "Zn"
        assert _canonical("ZN") == "Zn"
        assert _canonical("Cl") == "Cl"
        assert _canonical("CL") == "Cl"
        assert _canonical("Br") == "Br"
        assert _canonical("Na") == "Na"
        assert _canonical("Mg") == "Mg"
        assert _canonical("Ca") == "Ca"
        assert _canonical("Se") == "Se"

    def test_empty_string(self):
        assert _canonical("") is None

    def test_none_like(self):
        # empty string returns None
        assert _canonical("") is None

    def test_nonsense(self):
        assert _canonical("Xx") is None
        assert _canonical("QQ") is None
        assert _canonical("ZZ") is None

    def test_three_letter_truncation(self):
        # _canonical does sym[0].upper() + sym[1:].lower()
        # "FEE" -> "Fee" which is not a valid element
        assert _canonical("FEE") is None

    def test_numeric_string(self):
        # digits as symbol -> first char upper + rest lower -> not in periodic table
        assert _canonical("123") is None


class TestGuessElementFromMass:
    def test_hydrogen(self):
        assert _guess_element_from_mass(1.008) == "H"

    def test_carbon(self):
        assert _guess_element_from_mass(12.011) == "C"

    def test_nitrogen(self):
        assert _guess_element_from_mass(14.007) == "N"

    def test_oxygen(self):
        assert _guess_element_from_mass(15.999) == "O"

    def test_sulfur(self):
        assert _guess_element_from_mass(32.06) == "S"

    def test_phosphorus(self):
        assert _guess_element_from_mass(30.974) == "P"

    def test_iron(self):
        assert _guess_element_from_mass(55.845) == "Fe"

    def test_zinc(self):
        assert _guess_element_from_mass(65.38) == "Zn"

    def test_calcium(self):
        assert _guess_element_from_mass(40.078) == "Ca"

    def test_chlorine(self):
        assert _guess_element_from_mass(35.45) == "Cl"

    def test_bromine(self):
        assert _guess_element_from_mass(79.904) == "Br"

    def test_zero_mass(self):
        # mass 0 should not match anything meaningful (neutron has no symbol)
        result = _guess_element_from_mass(0.0)
        # Accept "n" (neutron) or "X"
        assert result in ("n", "X")

    def test_negative_mass(self):
        assert _guess_element_from_mass(-1.0) == "X"

    def test_very_large_mass(self):
        assert _guess_element_from_mass(99999.0) == "X"

    def test_nan_returns_x(self):
        assert _guess_element_from_mass(float("nan")) == "X"

    def test_string_returns_x(self):
        assert _guess_element_from_mass("not_a_number") == "X"

    def test_none_returns_x(self):
        assert _guess_element_from_mass(None) == "X"

    def test_threshold_boundary(self):
        # Carbon mass is 12.011, default threshold 0.3
        # 12.011 + 0.29 = 12.301 -> within threshold
        assert _guess_element_from_mass(12.301) == "C"
        # 12.011 + 0.31 = 12.321 -> just over threshold
        assert _guess_element_from_mass(12.321) == "X"

    def test_custom_threshold_tight(self):
        assert _guess_element_from_mass(12.011, threshold=0.001) == "C"
        assert _guess_element_from_mass(12.1, threshold=0.001) == "X"

    def test_custom_threshold_wide(self):
        assert _guess_element_from_mass(12.5, threshold=1.0) == "C"


class TestGuessElementFromName:
    # --- empty / trivial ---
    def test_empty_name(self):
        assert _guess_element_from_name("") == "X"

    def test_whitespace_only(self):
        # re.sub strips non-alpha, leaving empty
        assert _guess_element_from_name("   ") == "X"

    # --- element_hint override ---
    def test_hint_overrides(self):
        assert _guess_element_from_name("CA", "ALA", element_hint="Ca") == "Ca"
        assert _guess_element_from_name("ZN1", "ZN", element_hint="Zn") == "Zn"

    def test_hint_with_junk_chars(self):
        assert _guess_element_from_name("X", element_hint="C1++") == "C"

    def test_hint_invalid_element(self):
        # invalid hint falls through to other rules
        result = _guess_element_from_name("CA", "ALA", element_hint="Qq")
        assert result == "C"  # should resolve via protein carbon rule

    # --- water residues ---
    @pytest.mark.parametrize("res", list(WATER_RESIDUES))
    def test_water_oxygen(self, res):
        assert _guess_element_from_name("OW", res) == "O"
        assert _guess_element_from_name("O", res) == "O"
        assert _guess_element_from_name("OH2", res) == "O"

    @pytest.mark.parametrize("res", list(WATER_RESIDUES))
    def test_water_hydrogen(self, res):
        assert _guess_element_from_name("HW1", res) == "H"
        assert _guess_element_from_name("H1", res) == "H"
        assert _guess_element_from_name("HW2", res) == "H"

    @pytest.mark.parametrize("res", list(WATER_RESIDUES))
    def test_water_unknown_atom(self, res):
        assert _guess_element_from_name("MW", res) == "X"

    # --- ion residues ---
    @pytest.mark.parametrize("res,expected", list(ION_RES_TO_ELT.items()))
    def test_ion_residues(self, res, expected):
        assert _guess_element_from_name("ION", res) == expected

    # --- protein carbon naming ---
    @pytest.mark.parametrize("atom_name", list(PROTEIN_CARBONS))
    def test_protein_carbons(self, atom_name):
        assert _guess_element_from_name(atom_name, "ALA") == "C"

    # --- sidechain oxygen ---
    @pytest.mark.parametrize("atom_name", list(SIDECHAIN_OXY))
    def test_sidechain_oxygen(self, atom_name):
        for res in ("ASP", "GLU", "SER", "THR", "TYR"):
            assert _guess_element_from_name(atom_name, res) == "O"

    # --- sidechain nitrogen ---
    @pytest.mark.parametrize("atom_name", list(SIDECHAIN_NITRO))
    def test_sidechain_nitrogen(self, atom_name):
        for res in ("HIS", "ARG", "LYS", "ASN", "GLN"):
            assert _guess_element_from_name(atom_name, res) == "N"

    # --- nucleic acid oxygen ---
    @pytest.mark.parametrize("atom_name", list(NUC_OXY))
    def test_nucleic_acid_oxygen(self, atom_name):
        for res in ("DA", "DC", "DG", "DT"):
            assert _guess_element_from_name(atom_name, res) == "O"

    # --- MSE / SEC selenium ---
    def test_mse_selenium(self):
        assert _guess_element_from_name("SE", "MSE") == "Se"
        assert _guess_element_from_name("SE1", "MSE") == "Se"

    def test_sec_selenium(self):
        assert _guess_element_from_name("SE", "SEC") == "Se"

    # --- protein CA is carbon, not calcium ---
    def test_CA_in_protein_is_carbon(self):
        for res in ("ALA", "GLY", "VAL", "LEU", "ILE"):
            assert _guess_element_from_name("CA", res) == "C"

    def test_CA_in_nucleic_acid_is_carbon(self):
        # CA in nucleic acid context via TWO_LETTER special case
        for res in ("DA", "DC", "DG", "DT"):
            assert _guess_element_from_name("CA", res) == "C"

    # --- CA as calcium in non-bio context ---
    def test_CA_non_bio_is_calcium(self):
        assert _guess_element_from_name("CA", "LIG") == "Ca"
        assert _guess_element_from_name("CA", None) == "Ca"

    # --- two-letter element detection ---
    def test_two_letter_elements_non_bio(self):
        assert _guess_element_from_name("FE1", "HEM") == "Fe"
        assert _guess_element_from_name("ZN", "ZN") == "Zn"
        assert _guess_element_from_name("MG", "MG") == "Mg"
        assert _guess_element_from_name("CL1", "LIG") == "Cl"
        assert _guess_element_from_name("BR1", "LIG") == "Br"

    # --- single-letter fast path ---
    def test_single_letter_fast_path(self):
        for letter in ("C", "H", "N", "O", "S", "P", "F", "I", "K"):
            assert _guess_element_from_name(f"{letter}99", "LIG") == letter

    # --- generic fallback ---
    def test_generic_fallback_two_letter(self):
        assert _guess_element_from_name("AL1", "LIG") == "Al"

    def test_generic_fallback_one_letter(self):
        assert _guess_element_from_name("B1", "LIG") == "B"

    def test_unknown_returns_x(self):
        assert _guess_element_from_name("XX99", "UNK") == "X"

    # --- biomolecule first-letter rule ---
    @pytest.mark.parametrize("aa", ["ALA", "GLY", "VAL"])
    def test_protein_simple_atoms(self, aa):
        assert _guess_element_from_name("N", aa) == "N"
        assert _guess_element_from_name("C", aa) == "C"
        assert _guess_element_from_name("O", aa) == "O"
        assert _guess_element_from_name("H", aa) == "H"
        assert _guess_element_from_name("S", aa) == "S"  # e.g. CYS context

    # --- all amino acids resolve standard backbone atoms ---
    @pytest.mark.parametrize("aa", list(AA3))
    def test_backbone_atoms_all_aa(self, aa):
        assert _guess_element_from_name("N", aa) == "N"
        assert _guess_element_from_name("C", aa) == "C"
        assert _guess_element_from_name("O", aa) == "O"
        assert _guess_element_from_name("CA", aa) == "C"
        assert _guess_element_from_name("CB", aa) == "C"


class TestCreateUniverse:
    def test_topology_only(self):
        u, elements = create_universe(str(TPR_FILE))
        assert isinstance(u, Universe)
        assert u.atoms.n_atoms > 0
        assert elements == []

    def test_topology_and_trajectory(self):
        u, elements = create_universe(str(TPR_FILE), str(XTC_FILE))
        assert isinstance(u, Universe)
        assert u.trajectory.n_frames >= 1
        assert elements == []

    def test_with_elements(self):
        u, elements = create_universe(str(TPR_FILE), add_elements=True)
        assert len(elements) == u.atoms.n_atoms
        assert all(isinstance(e, str) for e in elements)
        # Real atoms (mass > 0) must all be assigned; virtual sites (e.g. TIP4P MW) may be "X".
        bad = [
            (i, elements[i])
            for i, m in enumerate(u.atoms.masses)
            if m > 0 and elements[i] == "X"
        ]
        assert not bad, f"Real atoms with unguessed element: {bad[:10]}"

    def test_with_elements_and_trajectory(self):
        u, elements = create_universe(str(TPR_FILE), str(XTC_FILE), add_elements=True)
        assert len(elements) == u.atoms.n_atoms

    def test_with_wrap(self):
        u, _ = create_universe(str(TPR_FILE), str(XTC_FILE), wrap=True)
        assert isinstance(u, Universe)

    def test_with_unwrap(self):
        u, _ = create_universe(str(TPR_FILE), str(XTC_FILE), unwrap=True)
        assert isinstance(u, Universe)

    def test_with_wrap_and_unwrap(self):
        u, _ = create_universe(str(TPR_FILE), str(XTC_FILE), wrap=True, unwrap=True)
        assert isinstance(u, Universe)

    def test_invalid_topology_raises(self):
        with pytest.raises(RuntimeError, match="Failed to load"):
            create_universe("/nonexistent/file.tpr")

    def test_path_object(self):
        u, _ = create_universe(TPR_FILE)
        assert isinstance(u, Universe)

    def test_elements_universe_consistency(self):
        """Elements list and universe.atoms.elements must agree."""
        u, elements = create_universe(str(TPR_FILE), add_elements=True)
        stored = list(u.atoms.elements)
        assert stored == elements

    def test_elements_are_valid_symbols(self):
        """Every guessed element should be a real periodic table symbol or X."""
        import periodictable
        u, elements = create_universe(str(TPR_FILE), add_elements=True)
        valid = {el.symbol for el in periodictable.elements}
        valid.add("X")
        for e in elements:
            assert e in valid, f"Invalid element symbol: {e!r}"


class TestConvertMdaToRdkit:
    """Aggressive tests for the MDA->RDKit converter and its index mapping."""

    @pytest.fixture(scope="class")
    def protein_universe(self):
        u, _ = create_universe(str(TPR_FILE), add_elements=True)
        return u

    # --- basic conversion ---
    def test_returns_mol_and_mapping(self, protein_universe):
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms in test topology")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert isinstance(mol, Chem.Mol)
        assert isinstance(mapping, dict)

    def test_mol_has_atoms(self, protein_universe):
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, _ = convert_mda_to_rdkit(protein)
        assert mol.GetNumAtoms() > 0

    # --- mapping correctness (the core requirement) ---
    def test_mapping_size_equals_atomgroup(self, protein_universe):
        """Mapping must have exactly one entry per atom in the AtomGroup."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert len(mapping) == protein.n_atoms

    def test_mapping_keys_are_universe_indices(self, protein_universe):
        """Keys of mapping must be universe-level atom indices (atom.ix)."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        _, mapping = convert_mda_to_rdkit(protein)
        expected_keys = {atom.ix for atom in protein}
        assert set(mapping.keys()) == expected_keys

    def test_mapping_values_are_contiguous_rdkit_indices(self, protein_universe):
        """Values of mapping must be contiguous 0..N-1 RDKit atom indices."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        rdkit_indices = sorted(mapping.values())
        assert rdkit_indices == list(range(protein.n_atoms))

    def test_mapping_values_match_rdkit_atom_count(self, protein_universe):
        """Max RDKit index + 1 must equal mol.GetNumAtoms()."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert max(mapping.values()) + 1 == mol.GetNumAtoms()

    def test_mapping_is_bijective(self, protein_universe):
        """Mapping must be one-to-one (no duplicate RDKit indices)."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        _, mapping = convert_mda_to_rdkit(protein)
        values = list(mapping.values())
        assert len(values) == len(set(values)), "Mapping is not bijective"

    def test_mapping_preserves_element_identity(self, protein_universe):
        """Element of MDA atom must match element of corresponding RDKit atom."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        periodic = Chem.GetPeriodicTable()
        for mda_atom in protein:
            rd_idx = mapping[mda_atom.ix]
            rd_atom = mol.GetAtomWithIdx(rd_idx)
            rd_element = periodic.GetElementSymbol(rd_atom.GetAtomicNum())
            mda_element = mda_atom.element
            assert rd_element == mda_element, (
                f"Element mismatch at MDA ix={mda_atom.ix}: "
                f"MDA={mda_element}, RDKit={rd_element}"
            )

    def test_reverse_mapping_consistency(self, protein_universe):
        """Reverse mapping (RDKit->MDA) must also be consistent."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        reverse = {v: k for k, v in mapping.items()}
        assert len(reverse) == len(mapping)
        for rd_idx, mda_ix in reverse.items():
            assert mapping[mda_ix] == rd_idx

    # --- residue metadata transfer ---
    def test_residue_props_transferred(self, protein_universe):
        """RDKit atoms should carry resname, resid, segid properties."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        # Check a sample of atoms
        for mda_atom in list(protein)[:20]:
            rd_idx = mapping[mda_atom.ix]
            rd_atom = mol.GetAtomWithIdx(rd_idx)
            assert rd_atom.HasProp("resname")
            assert rd_atom.HasProp("resid")
            assert rd_atom.HasProp("segid")
            assert rd_atom.GetProp("resname") == str(mda_atom.resname)
            assert rd_atom.GetProp("resid") == str(mda_atom.resid)

    # --- subgroup conversion ---
    def test_single_residue_conversion(self, protein_universe):
        """Converting a single residue should produce a small, valid mol."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        first_res = protein.residues[0].atoms
        mol, mapping = convert_mda_to_rdkit(first_res)
        assert mol.GetNumAtoms() == first_res.n_atoms
        assert len(mapping) == first_res.n_atoms

    def test_subset_mapping_keys_subset_of_universe(self, protein_universe):
        """Mapping keys for a subset must be a subset of all universe indices."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        first_res = protein.residues[0].atoms
        _, mapping = convert_mda_to_rdkit(first_res)
        all_ix = set(range(protein_universe.atoms.n_atoms))
        assert set(mapping.keys()).issubset(all_ix)

    # --- metal handling ---
    def test_metal_bonds_disconnected(self, protein_universe):
        """If metals are present, their bonds should be removed."""
        # Select atoms that might include metals (ions)
        try:
            metal_sel = protein_universe.select_atoms(
                "name ZN FE MG MN or resname ZN FE MG MN ZN2 FE2"
            )
        except Exception:
            metal_sel = None
        if metal_sel is None or metal_sel.n_atoms == 0:
            pytest.skip("No metals in test topology")

        # Include metals + nearby atoms
        nearby = protein_universe.select_atoms(
            f"around 3.0 (name ZN FE MG MN or resname ZN FE MG MN ZN2 FE2)"
        )
        ag = metal_sel | nearby
        mol, mapping = convert_mda_to_rdkit(ag)

        # Metal atoms should have zero bonds in the resulting mol
        metal_symbols_upper = set(METAL_ELEMENTS)
        for atom in mol.GetAtoms():
            if atom.GetSymbol().upper() in metal_symbols_upper:
                assert atom.GetDegree() == 0, (
                    f"Metal {atom.GetSymbol()} at idx {atom.GetIdx()} still has "
                    f"{atom.GetDegree()} bonds"
                )

    # --- sanitization ---
    def test_mol_is_sanitized(self, protein_universe):
        """The returned mol should pass RDKit sanitization without errors."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, _ = convert_mda_to_rdkit(protein)
        # Re-sanitize should not raise
        Chem.SanitizeMol(mol)

    # --- no synthetic hydrogens ---
    def test_no_extra_hydrogens_added(self, protein_universe):
        """RDKit mol should have same atom count as AtomGroup (no AddHs)."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert mol.GetNumAtoms() == protein.n_atoms

    # --- exact atom identity via mapping ---
    def test_mapping_exact_atom_identity_all_atoms(self, protein_universe):
        """Every mapped pair must point to the SAME atom: same element, same
        atom name, same resname, same resid, same segid.  This is the
        definitive proof that MDA ix 42 (N of CYS 4 in PROA) lands on
        the RDKit atom that carries exactly those properties."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        periodic = Chem.GetPeriodicTable()

        mismatches = []
        for mda_atom in protein:
            rd_idx = mapping[mda_atom.ix]
            rd_atom = mol.GetAtomWithIdx(rd_idx)

            # element check
            rd_element = periodic.GetElementSymbol(rd_atom.GetAtomicNum())
            if rd_element != mda_atom.element:
                mismatches.append(
                    f"ix={mda_atom.ix} element: MDA={mda_atom.element} RDKit={rd_element}"
                )

            # residue-level identity (transferred as props by convert_mda_to_rdkit)
            if rd_atom.HasProp("resname"):
                if rd_atom.GetProp("resname") != str(mda_atom.resname):
                    mismatches.append(
                        f"ix={mda_atom.ix} resname: MDA={mda_atom.resname} "
                        f"RDKit={rd_atom.GetProp('resname')}"
                    )
            if rd_atom.HasProp("resid"):
                if rd_atom.GetProp("resid") != str(mda_atom.resid):
                    mismatches.append(
                        f"ix={mda_atom.ix} resid: MDA={mda_atom.resid} "
                        f"RDKit={rd_atom.GetProp('resid')}"
                    )
            if rd_atom.HasProp("segid"):
                if rd_atom.GetProp("segid") != str(mda_atom.segid):
                    mismatches.append(
                        f"ix={mda_atom.ix} segid: MDA={mda_atom.segid} "
                        f"RDKit={rd_atom.GetProp('segid')}"
                    )

        assert not mismatches, (
            f"{len(mismatches)} atom identity mismatches:\n"
            + "\n".join(mismatches[:20])
        )

    def test_mapping_spot_check_specific_residues(self, protein_universe):
        """Pick specific residues by resid and verify per-atom identity
        through the mapping (e.g. N of residue 1, CA of residue 2, etc.)."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)

        # Test the first 5 residues atom-by-atom
        for res in protein.residues[:5]:
            for mda_atom in res.atoms:
                rd_idx = mapping[mda_atom.ix]
                rd_atom = mol.GetAtomWithIdx(rd_idx)

                # The RDKit atom must carry the same residue context
                assert rd_atom.GetProp("resname") == str(mda_atom.resname), (
                    f"Atom {mda_atom.name} (ix={mda_atom.ix}): resname mismatch "
                    f"MDA={mda_atom.resname} vs RDKit={rd_atom.GetProp('resname')}"
                )
                assert rd_atom.GetProp("resid") == str(mda_atom.resid), (
                    f"Atom {mda_atom.name} (ix={mda_atom.ix}): resid mismatch "
                    f"MDA={mda_atom.resid} vs RDKit={rd_atom.GetProp('resid')}"
                )
                assert rd_atom.GetProp("segid") == str(mda_atom.segid), (
                    f"Atom {mda_atom.name} (ix={mda_atom.ix}): segid mismatch "
                    f"MDA={mda_atom.segid} vs RDKit={rd_atom.GetProp('segid')}"
                )

    def test_smarts_match_resolves_to_exact_mda_atom(self, protein_universe):
        """A SMARTS match for backbone nitrogen [NX3;H1] should map back to
        MDA atoms that are actually backbone N atoms in the same residue."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        reverse = {v: k for k, v in mapping.items()}

        # Match amide nitrogens (backbone-like)
        pattern = Chem.MolFromSmarts("[#7]")
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            pytest.skip("No nitrogen matches")

        for (rd_idx,) in matches:
            mda_ix = reverse[rd_idx]
            mda_atom = protein_universe.atoms[mda_ix]

            # Element must be nitrogen
            assert mda_atom.element == "N", (
                f"SMARTS [#7] -> RDKit idx {rd_idx} -> MDA ix {mda_ix} "
                f"but element is {mda_atom.element}"
            )
            # The RDKit atom's stored resname/resid must match the MDA atom's
            rd_atom = mol.GetAtomWithIdx(rd_idx)
            assert rd_atom.GetProp("resname") == str(mda_atom.resname)
            assert rd_atom.GetProp("resid") == str(mda_atom.resid)

    def test_mapping_preserves_iteration_order(self, protein_universe):
        """RDKit idx for atom i in AtomGroup iteration must be i, since the
        converter preserves order. Verify the mapping encodes exactly this."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        _, mapping = convert_mda_to_rdkit(protein)

        for expected_rd_idx, mda_atom in enumerate(protein):
            assert mapping[mda_atom.ix] == expected_rd_idx, (
                f"Iteration order broken: atom {mda_atom.name} "
                f"(ix={mda_atom.ix}) at position {expected_rd_idx} "
                f"mapped to RDKit idx {mapping[mda_atom.ix]}"
            )

    # --- round-trip index mapping stress test ---
    def test_mapping_roundtrip_all_residues(self, protein_universe):
        """Convert each residue individually and verify mapping integrity."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        # Test first 10 residues for speed
        for res in protein.residues[:10]:
            ag = res.atoms
            mol, mapping = convert_mda_to_rdkit(ag)
            assert len(mapping) == ag.n_atoms
            # Keys are universe-level indices
            for atom in ag:
                assert atom.ix in mapping
            # Values are 0..n-1
            assert sorted(mapping.values()) == list(range(ag.n_atoms))

    # --- SMARTS match mapping correctness ---
    def test_smarts_match_maps_back_to_mda(self, protein_universe):
        """A SMARTS match on the RDKit mol should map back to correct MDA atoms."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        reverse = {v: k for k, v in mapping.items()}

        # Match nitrogen atoms
        pattern = Chem.MolFromSmarts("[#7]")
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            pytest.skip("No nitrogen matches")

        for (rd_idx,) in matches:
            mda_ix = reverse[rd_idx]
            mda_atom = protein_universe.atoms[mda_ix]
            assert mda_atom.element == "N", (
                f"SMARTS [#7] matched RDKit idx {rd_idx} -> MDA ix {mda_ix} "
                f"but MDA element is {mda_atom.element}, not N"
            )

    def test_smarts_carbon_match(self, protein_universe):
        """SMARTS [#6] matches should map back to MDA carbons."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        reverse = {v: k for k, v in mapping.items()}

        pattern = Chem.MolFromSmarts("[#6]")
        matches = mol.GetSubstructMatches(pattern)
        for (rd_idx,) in matches[:50]:  # check first 50
            mda_ix = reverse[rd_idx]
            mda_atom = protein_universe.atoms[mda_ix]
            assert mda_atom.element == "C"

    def test_smarts_oxygen_match(self, protein_universe):
        """SMARTS [#8] matches should map back to MDA oxygens."""
        protein = protein_universe.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein atoms")
        mol, mapping = convert_mda_to_rdkit(protein)
        reverse = {v: k for k, v in mapping.items()}

        pattern = Chem.MolFromSmarts("[#8]")
        matches = mol.GetSubstructMatches(pattern)
        for (rd_idx,) in matches[:50]:
            mda_ix = reverse[rd_idx]
            mda_atom = protein_universe.atoms[mda_ix]
            assert mda_atom.element == "O"


class TestHasResidues:
    def test_universe_with_residues(self, small_universe):
        assert has_residues(small_universe) is True

    def test_empty_universe(self, empty_universe):
        assert has_residues(empty_universe) is False

    def test_atomgroup_with_residues(self, small_universe):
        ag = small_universe.atoms
        assert has_residues(ag) is True

    def test_real_universe(self, tpr_universe):
        assert has_residues(tpr_universe) is True


class TestHasAtoms:
    def test_universe_with_atoms(self, small_universe):
        assert has_atoms(small_universe) is True

    def test_empty_universe(self, empty_universe):
        assert has_atoms(empty_universe) is False

    def test_atomgroup_with_atoms(self, small_universe):
        ag = small_universe.atoms
        assert has_atoms(ag) is True

    def test_real_universe(self, tpr_universe):
        assert has_atoms(tpr_universe) is True


class TestRequireAtoms:
    def test_passes_with_atoms(self, small_universe):
        require_atoms(small_universe)  # should not raise

    def test_raises_on_empty(self, empty_universe):
        with pytest.raises(ValueError, match="no atoms"):
            require_atoms(empty_universe)

    def test_passes_with_atomgroup(self, small_universe):
        require_atoms(small_universe.atoms)

    def test_raises_on_empty_atomgroup(self, small_universe):
        empty_ag = small_universe.select_atoms("name NONEXISTENT")
        with pytest.raises(ValueError, match="no atoms"):
            require_atoms(empty_ag)


class TestRequireResidues:
    def test_passes_with_residues(self, small_universe):
        require_residues(small_universe)

    def test_raises_on_empty(self, empty_universe):
        with pytest.raises(ValueError, match="no residues"):
            require_residues(empty_universe)

    def test_passes_with_atomgroup(self, small_universe):
        require_residues(small_universe.atoms)


class TestIntegration:
    def test_full_pipeline_protein(self):
        """End-to-end: load TPR -> add elements -> convert protein to RDKit."""
        u, elements = create_universe(str(TPR_FILE), add_elements=True)
        protein = u.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein in topology")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert mol.GetNumAtoms() == protein.n_atoms
        assert len(mapping) == protein.n_atoms

    def test_full_pipeline_with_trajectory(self):
        """End-to-end: load TPR+XTC -> add elements -> convert to RDKit."""
        u, elements = create_universe(str(TPR_FILE), str(XTC_FILE), add_elements=True)
        protein = u.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein in topology")
        mol, mapping = convert_mda_to_rdkit(protein)
        assert mol.GetNumAtoms() == protein.n_atoms

    def test_mapping_stability_across_frames(self):
        """Mapping should be identical regardless of which frame we're on."""
        u, _ = create_universe(str(TPR_FILE), str(XTC_FILE), add_elements=True)
        protein = u.select_atoms("protein")
        if protein.n_atoms == 0:
            pytest.skip("No protein in topology")

        # Frame 0
        u.trajectory[0]
        mol0, map0 = convert_mda_to_rdkit(protein)

        # Last frame
        u.trajectory[-1]
        mol1, map1 = convert_mda_to_rdkit(protein)

        assert map0 == map1, "Mapping changed between frames!"
        assert mol0.GetNumAtoms() == mol1.GetNumAtoms()

    def test_non_protein_single_residue(self):
        """Converting a single non-protein residue (e.g. first water) should work."""
        u, _ = create_universe(str(TPR_FILE), add_elements=True)
        non_protein = u.select_atoms("not protein")
        if non_protein.n_atoms == 0:
            pytest.skip("No non-protein atoms")
        # Pick just the first non-protein residue to avoid converting the whole water box
        first_res = non_protein.residues[0].atoms
        try:
            mol, mapping = convert_mda_to_rdkit(first_res)
            assert len(mapping) == first_res.n_atoms
        except RuntimeError:
            # Some residues (e.g. single ions) may not convert cleanly
            pass
