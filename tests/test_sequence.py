"""
Tests for pharmacon.analyzer.sequence — get_sequence and sequence_dict_to_fasta
"""
import pytest
import numpy as np
import MDAnalysis as Mda

from pharmacon.analyzer.sequence import get_sequence, sequence_dict_to_fasta


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_protein_universe(aa3_list, chainids=None, resids=None):
    """Build a minimal Universe with one CA per residue (protein-like)."""
    n = len(aa3_list)
    u = Mda.Universe.empty(
        n_atoms=n,
        n_residues=n,
        n_segments=1,
        atom_resindex=list(range(n)),
        residue_segindex=[0] * n,
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["CA"] * n)
    u.add_TopologyAttr("resnames", list(aa3_list))
    u.add_TopologyAttr("resids", list(range(1, n + 1)) if resids is None else list(resids))
    u.add_TopologyAttr("segids", ["SYS"])
    cids = list(chainids) if chainids is not None else [""] * n
    u.add_TopologyAttr("chainIDs", cids)
    u.atoms.positions = np.zeros((n, 3), dtype=np.float32)
    return u


# ---------------------------------------------------------------------------
# get_sequence — no chain IDs
# ---------------------------------------------------------------------------

class TestGetSequenceNoChain:
    def test_result_key_is_no_chainid(self):
        u = _make_protein_universe(["ALA", "GLY"])
        result = get_sequence(u)
        assert "NO_CHAINID" in result

    def test_aa1_seq_correct(self):
        u = _make_protein_universe(["ALA", "GLY", "PHE"])
        result = get_sequence(u)
        assert result["NO_CHAINID"]["aa1_seq"] == "AGF"

    def test_aa3_list_correct(self):
        u = _make_protein_universe(["ALA", "GLY"])
        result = get_sequence(u)
        assert result["NO_CHAINID"]["aa3_list"] == ["ALA", "GLY"]

    def test_resids_topology_dtype_int(self):
        u = _make_protein_universe(["ALA", "GLY"])
        result = get_sequence(u)
        assert result["NO_CHAINID"]["resids_topology"].dtype == int

    def test_resid_seq_starts_at_one(self):
        u = _make_protein_universe(["ALA", "GLY", "VAL"])
        result = get_sequence(u)
        np.testing.assert_array_equal(
            result["NO_CHAINID"]["resid_seq"], np.array([1, 2, 3])
        )

    def test_resids_topology_matches_input(self):
        u = _make_protein_universe(["ALA", "GLY"], resids=[10, 20])
        result = get_sequence(u)
        np.testing.assert_array_equal(
            result["NO_CHAINID"]["resids_topology"], np.array([10, 20])
        )

    def test_no_protein_atoms_raises(self):
        u = Mda.Universe.empty(
            n_atoms=2, n_residues=1, n_segments=1,
            atom_resindex=[0, 0], residue_segindex=[0], trajectory=True,
        )
        u.add_TopologyAttr("names", ["N", "C"])
        u.add_TopologyAttr("resnames", ["WATER"])
        u.add_TopologyAttr("resids", [1])
        u.add_TopologyAttr("segids", ["SYS"])
        u.add_TopologyAttr("chainIDs", ["", ""])
        u.atoms.positions = np.zeros((2, 3), dtype=np.float32)
        with pytest.raises(RuntimeError, match="No protein CA atoms"):
            get_sequence(u)

    def test_single_residue(self):
        u = _make_protein_universe(["ALA"])
        result = get_sequence(u)
        assert result["NO_CHAINID"]["aa1_seq"] == "A"
        assert len(result["NO_CHAINID"]["aa3_list"]) == 1


# ---------------------------------------------------------------------------
# get_sequence — with chain IDs
# ---------------------------------------------------------------------------

class TestGetSequenceWithChains:
    def test_two_chains_produce_two_keys(self):
        u = _make_protein_universe(
            ["ALA", "GLY", "PHE", "LYS"],
            chainids=["A", "A", "B", "B"],
        )
        result = get_sequence(u)
        assert set(result.keys()) == {"A", "B"}

    def test_chain_a_sequence(self):
        u = _make_protein_universe(
            ["ALA", "GLY", "PHE", "LYS"],
            chainids=["A", "A", "B", "B"],
        )
        result = get_sequence(u)
        assert result["A"]["aa1_seq"] == "AG"

    def test_chain_b_sequence(self):
        u = _make_protein_universe(
            ["ALA", "GLY", "PHE", "LYS"],
            chainids=["A", "A", "B", "B"],
        )
        result = get_sequence(u)
        assert result["B"]["aa1_seq"] == "FK"

    def test_single_chain_explicit_id(self):
        u = _make_protein_universe(
            ["ALA", "GLY"],
            chainids=["X", "X"],
        )
        result = get_sequence(u)
        assert "X" in result
        assert "NO_CHAINID" not in result


# ---------------------------------------------------------------------------
# sequence_dict_to_fasta
# ---------------------------------------------------------------------------

class TestSequenceDictToFasta:
    @pytest.fixture
    def simple_data(self):
        return {
            "A": {
                "aa1_seq": "AGFK",
                "resids_topology": np.array([1, 2, 3, 4]),
                "aa3_list": ["ALA", "GLY", "PHE", "LYS"],
                "resid_seq": np.array([1, 2, 3, 4]),
            }
        }

    def test_header_format(self, simple_data):
        fasta = sequence_dict_to_fasta(simple_data)
        assert fasta.startswith(">chain:A")

    def test_header_contains_resid_range(self, simple_data):
        fasta = sequence_dict_to_fasta(simple_data)
        assert "resid=1-4" in fasta

    def test_sequence_in_output(self, simple_data):
        fasta = sequence_dict_to_fasta(simple_data)
        assert "AGFK" in fasta

    def test_custom_header_prefix(self, simple_data):
        fasta = sequence_dict_to_fasta(simple_data, header_prefix="seq")
        assert fasta.startswith(">seq:A")

    def test_wrapping_at_80(self):
        long_seq = "A" * 160
        data = {
            "A": {
                "aa1_seq": long_seq,
                "resids_topology": np.array([1, 160]),
                "aa3_list": [],
                "resid_seq": np.arange(1, 161),
            }
        }
        fasta = sequence_dict_to_fasta(data, wrap=80)
        lines = fasta.split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        assert all(len(l) <= 80 for l in seq_lines)
        assert len(seq_lines) == 2  # 160 chars / 80 = 2 lines

    def test_no_wrapping_when_wrap_zero(self):
        seq = "A" * 200
        data = {
            "A": {
                "aa1_seq": seq,
                "resids_topology": np.array([1, 200]),
                "aa3_list": [],
                "resid_seq": np.arange(1, 201),
            }
        }
        fasta = sequence_dict_to_fasta(data, wrap=0)
        lines = fasta.split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        assert len(seq_lines) == 1
        assert len(seq_lines[0]) == 200

    def test_two_chains_sorted_by_key(self):
        data = {
            "B": {"aa1_seq": "GG", "resids_topology": np.array([5, 6]), "aa3_list": [], "resid_seq": np.array([1, 2])},
            "A": {"aa1_seq": "AA", "resids_topology": np.array([1, 2]), "aa3_list": [], "resid_seq": np.array([1, 2])},
        }
        fasta = sequence_dict_to_fasta(data)
        a_pos = fasta.index(">chain:A")
        b_pos = fasta.index(">chain:B")
        assert a_pos < b_pos  # A sorted before B

    def test_empty_aa1_seq_skipped(self):
        data = {
            "A": {"aa1_seq": "", "resids_topology": np.array([]), "aa3_list": [], "resid_seq": np.array([])},
            "B": {"aa1_seq": "AG", "resids_topology": np.array([1, 2]), "aa3_list": [], "resid_seq": np.array([1, 2])},
        }
        fasta = sequence_dict_to_fasta(data)
        assert ">chain:A" not in fasta
        assert ">chain:B" in fasta

    def test_no_resids_no_range_in_header(self):
        data = {
            "A": {"aa1_seq": "AG", "resids_topology": np.array([]), "aa3_list": [], "resid_seq": np.array([])},
        }
        fasta = sequence_dict_to_fasta(data)
        assert "resid=" not in fasta
