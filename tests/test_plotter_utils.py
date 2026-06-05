"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pure-utility functions in pharmacon.plotter.interactions.
"""
from __future__ import annotations

import pytest

from pharmacon.plotter.interactions import (
    parse_range_dict,
    parse_frame_interaction_record,
    _PLI_STACK_ORDER,
    _PLI_STACK_RANK,
)


class TestParseRangeDict:
    def test_dict_input_passthrough(self):
        result = parse_range_dict({"A": (1, 10), "B": (20, 30)})
        assert result == {"A": (1, 10), "B": (20, 30)}

    def test_string_literal_parsed(self):
        result = parse_range_dict("{'A': (1, 10)}")
        assert result == {"A": (1, 10)}

    def test_empty_string_returns_empty_dict(self):
        assert parse_range_dict("") == {}
        assert parse_range_dict("   ") == {}

    def test_list_values_coerced_to_tuples(self):
        result = parse_range_dict({"A": [1, 5]})
        assert result["A"] == (1, 5)

    def test_invalid_syntax_raises(self):
        with pytest.raises(ValueError, match="Invalid range mapping syntax"):
            parse_range_dict("not-a-dict")

    def test_non_dict_input_raises(self):
        with pytest.raises(ValueError, match="must be dict or string"):
            parse_range_dict(42)

    def test_non_dict_after_eval_raises(self):
        with pytest.raises(ValueError, match="must be a dict"):
            parse_range_dict("[1, 2, 3]")

    def test_non_string_key_raises(self):
        with pytest.raises(ValueError, match="Invalid key type"):
            parse_range_dict({1: (1, 2)})

    def test_value_not_pair_raises(self):
        with pytest.raises(ValueError, match="tuple of length 2"):
            parse_range_dict({"A": (1, 2, 3)})

    def test_negative_value_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            parse_range_dict({"A": (-1, 2)})

    def test_lower_greater_than_upper_raises(self):
        with pytest.raises(ValueError, match="Lower bound > upper bound"):
            parse_range_dict({"A": (10, 5)})

    def test_equal_bounds_accepted(self):
        result = parse_range_dict({"A": (5, 5)})
        assert result["A"] == (5, 5)


# Reference record layouts (frame_number=0 is dropped before this function):
#   Atom-atom: [label, a1_idx, a1_name, a1_id, a1_type, a1_element,
#               a1_resn, a1_resid, a1_chain, a1_segid,
#               a2_idx, a2_name, ...]
# parse_frame_interaction_record expects the record WITHOUT frame_number.

class TestParseFrameInteractionRecord:
    def _atom_atom_record(self, label: str) -> list:
        return [
            label,
            1, "CA", 2, "C", "BB",  # atom1 metadata
            "PHE", 45, "A", "PROA",  # atom1 residue: rec[6..9]
            10, "CB", 11, "C", "SC",  # atom2 metadata
            "LIG", 1, "A", "LIGA",   # atom2 residue: rec[15..18]
        ]

    def test_hydrophobic_extracts_protein_and_ligand(self):
        rec = self._atom_atom_record("HYDROPHOBIC")
        # ligand_atom is at index 11 (atom2 name "CB")
        out = parse_frame_interaction_record(rec)
        assert out["interaction"] == "HYDROPHOBIC"
        assert out["protein"] == {"resname": "PHE", "resid": 45,
                                   "chainid": "A", "segid": "PROA"}
        assert out["ligand_atoms"] == ["CB"]

    def test_hydrogen_bond_uppercased(self):
        rec = self._atom_atom_record("hydrogen-bond")
        out = parse_frame_interaction_record(rec)
        assert out["interaction"] == "HYDROGEN-BOND"
        assert out["protein"]["resid"] == 45

    def test_pi_stacking_splits_comma_separated_ring(self):
        # PI-STACKING: protein at indices 6..9 are comma-separated (one
        # entry per ring atom). The function reads ligand atoms from
        # rec[12] verbatim — we mirror that contract here.
        rec = [
            "PI-STACKING",
            "474,475,477,479,488",          # 1
            "CG,CD1,NE1,CE2,CD2",           # 2
            "475,476,478,480,489",          # 3
            "C,C,N,C,C",                    # 4
            "SC,SC,SC,SC,SC",               # 5
            "TRP,TRP,TRP,TRP,TRP",          # 6 → resname
            "32,32,32,32,32",               # 7 → resid
            "A,A,A,A,A",                    # 8 → chainid
            "PROA,PROA,PROA,PROA,PROA",     # 9 → segid
            "100,101,102,103,104,105",      # 10
            "ring2_names_unused",           # 11
            "C1,C2,C3,C4,C5,C6",            # 12 → ligand_atoms per function
            "C,C,C,C,C,C",                  # 13
            "SC,SC,SC,SC,SC,SC",            # 14
            "LIG,LIG,LIG,LIG,LIG,LIG",      # 15
            "1,1,1,1,1,1",                  # 16
            "A,A,A,A,A,A",                  # 17
            "LIGA,LIGA,LIGA,LIGA,LIGA,LIGA",  # 18
        ]
        out = parse_frame_interaction_record(rec)
        assert out["interaction"] == "PI-STACKING"
        assert out["protein"] == {"resname": "TRP", "resid": 32,
                                   "chainid": "A", "segid": "PROA"}
        assert out["ligand_atoms"] == ["C1", "C2", "C3", "C4", "C5", "C6"]

    def test_pi_cation_layout(self):
        # PI-CATION: g1_role at index 1; protein info at indices 7..10
        # (single ring residue, NOT comma-separated for the cation side).
        rec = [
            "PI-CATION",
            "ring",  # g1_role
            "682,683,685,687,689,691",  # ring atom indices
            "CG,CD1,CE1,CZ,CE2,CD2",
            "683,684,686,688,690,692",
            "C,C,C,C,C,C",
            "SC,SC,SC,SC,SC,SC",
            "PHE",   # 7 → protein resname
            45,      # 8 → resid
            "A",     # 9 → chain
            "PROA",  # 10 → segid
            "NH1",   # cation atom name
            "LIG, LIG, LIG",  # 12 → ligand atoms (comma-joined)
        ]
        out = parse_frame_interaction_record(rec)
        assert out["interaction"] == "PI-CATION"
        assert out["protein"] == {"resname": "PHE", "resid": 45,
                                   "chainid": "A", "segid": "PROA"}
        assert out["ligand_atoms"] == ["LIG", "LIG", "LIG"]

    def test_normal_with_no_ligand_field(self):
        # Trim record so index 11 is missing → ligand_atoms must be empty.
        rec = self._atom_atom_record("IONIC")[:11]
        out = parse_frame_interaction_record(rec)
        assert out["interaction"] == "IONIC"
        assert out["ligand_atoms"] == []


class TestPliStackOrder:
    def test_constants_length_match(self):
        assert len(_PLI_STACK_ORDER) == len(_PLI_STACK_RANK)
        assert set(_PLI_STACK_ORDER) == set(_PLI_STACK_RANK)

    def test_hydrophobic_is_first(self):
        # Hydrophobic should sit at the bottom of every PLI stack.
        assert _PLI_STACK_ORDER[0] == "hydrophobic"
        assert _PLI_STACK_RANK["hydrophobic"] == 0

    def test_chaotic_order_sorts_to_canonical(self):
        chaos = list(reversed(_PLI_STACK_ORDER))
        chaos.sort(key=lambda k: (_PLI_STACK_RANK.get(k, len(_PLI_STACK_ORDER)), k))
        assert tuple(chaos) == _PLI_STACK_ORDER

    def test_unknown_labels_go_to_end(self):
        items = ["zzz_unknown", "hydrophobic", "aaa_unknown", "metal_contact"]
        items.sort(key=lambda k: (_PLI_STACK_RANK.get(k, len(_PLI_STACK_ORDER)), k))
        # Known labels first (in canonical order), then unknowns alphabetically.
        assert items.index("hydrophobic") < items.index("metal_contact")
        assert items.index("metal_contact") < items.index("aaa_unknown")
        assert items.index("aaa_unknown") < items.index("zzz_unknown")
