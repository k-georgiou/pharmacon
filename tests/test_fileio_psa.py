"""
Tests for pharmacon.fileio.psa.PharmaconPSAFile
"""
import pytest
import numpy as np
from pathlib import Path

from pharmacon.fileio.psa import PharmaconPSAFile


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open(path, mode="a"):
    return PharmaconPSAFile(path, mode=mode)


def _make_sequences():
    """Minimal two-chain sequence dict matching the write_sequence schema."""
    return {
        "A": {
            "aa1_seq": "ACDE",
            "aa3_list": ["ALA", "CYS", "ASP", "GLU"],
            "resids_topology": np.array([1, 2, 3, 4], dtype=int),
            "resid_seq": np.array([1, 2, 3, 4], dtype=int),
        },
        "B": {
            "aa1_seq": "FG",
            "aa3_list": ["PHE", "GLY"],
            "resids_topology": np.array([1, 2], dtype=int),
            "resid_seq": np.array([101, 102], dtype=int),
        },
    }


# ---------------------------------------------------------------------------
# Schema constants
# ---------------------------------------------------------------------------

class TestPSASchemaConstants:
    def test_scalar_property_columns_non_empty(self):
        assert len(PharmaconPSAFile.SCALAR_PROPERTY_COLUMNS) > 0

    def test_scalar_property_columns_contains_expected(self):
        cols = PharmaconPSAFile.SCALAR_PROPERTY_COLUMNS
        assert "molecular_weight" in cols
        assert "logP" in cols
        assert "tpsa" in cols

    def test_int_scalar_columns_subset_of_scalar(self):
        assert PharmaconPSAFile.INT_SCALAR_COLUMNS.issubset(
            set(PharmaconPSAFile.SCALAR_PROPERTY_COLUMNS)
        )

    def test_fingerprint_fields_has_four_types(self):
        assert len(PharmaconPSAFile.FINGERPRINT_FIELDS) == 4
        assert "morgan" in PharmaconPSAFile.FINGERPRINT_FIELDS


# ---------------------------------------------------------------------------
# write_sequence / read_sequence roundtrip
# ---------------------------------------------------------------------------

class TestSequenceRoundtrip:
    def test_write_and_read_returns_both_chains(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            result = f.read_sequence()
        assert set(result.keys()) == {"A", "B"}

    def test_aa1_seq_roundtrip(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            result = f.read_sequence()
        assert result["A"]["aa1_seq"] == "ACDE"
        assert result["B"]["aa1_seq"] == "FG"

    def test_aa3_list_roundtrip(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            result = f.read_sequence()
        assert result["A"]["aa3_list"] == ["ALA", "CYS", "ASP", "GLU"]

    def test_resid_seq_roundtrip(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            result = f.read_sequence()
        np.testing.assert_array_equal(result["A"]["resid_seq"], [1, 2, 3, 4])

    def test_resids_topology_roundtrip(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            result = f.read_sequence()
        np.testing.assert_array_equal(result["A"]["resids_topology"], [1, 2, 3, 4])

    def test_duplicate_write_without_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            with pytest.raises(ValueError):
                f.write_sequence(sequences=seqs, overwrite=False)

    def test_duplicate_write_with_overwrite_succeeds(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        new_seqs = {
            "X": {
                "aa1_seq": "W",
                "aa3_list": ["TRP"],
                "resids_topology": np.array([1], dtype=int),
                "resid_seq": np.array([1], dtype=int),
            }
        }
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence(sequences=new_seqs, overwrite=True)
            result = f.read_sequence()
        assert "X" in result
        assert "A" not in result

    def test_read_missing_sequence_group_raises(self, tmp_path):
        path = tmp_path / "test.psa"
        with _open(path, mode="w") as f:
            with pytest.raises(KeyError):
                f.read_sequence(group_name="nonexistent")


# ---------------------------------------------------------------------------
# write_sequence_fasta
# ---------------------------------------------------------------------------

class TestSequenceFASTA:
    def test_fasta_file_created(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        fasta_path = tmp_path / "out.fasta"
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_fasta(fasta_path)
        assert fasta_path.exists()

    def test_fasta_contains_chain_headers(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        fasta_path = tmp_path / "out.fasta"
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_fasta(fasta_path)
        text = fasta_path.read_text()
        assert ">A" in text
        assert ">B" in text

    def test_fasta_contains_sequences(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        fasta_path = tmp_path / "out.fasta"
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_fasta(fasta_path)
        text = fasta_path.read_text()
        assert "ACDE" in text
        assert "FG" in text

    def test_fasta_no_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.psa"
        seqs = _make_sequences()
        fasta_path = tmp_path / "out.fasta"
        fasta_path.write_text("existing")
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            with pytest.raises(FileExistsError):
                f.write_sequence_fasta(fasta_path, overwrite=False)


# ---------------------------------------------------------------------------
# write_sequence_to_csv
# ---------------------------------------------------------------------------

class TestSequenceCSV:
    def test_csv_file_created(self, tmp_path):
        path = tmp_path / "test.psa"
        csv_path = tmp_path / "out.csv"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_to_csv(csv_path)
        assert csv_path.exists()

    def test_csv_has_header_row(self, tmp_path):
        path = tmp_path / "test.psa"
        csv_path = tmp_path / "out.csv"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_to_csv(csv_path)
        lines = csv_path.read_text().splitlines()
        assert "chain_id" in lines[0]

    def test_csv_row_count_matches_total_residues(self, tmp_path):
        path = tmp_path / "test.psa"
        csv_path = tmp_path / "out.csv"
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            f.write_sequence_to_csv(csv_path)
        lines = csv_path.read_text().splitlines()
        # header + 4 (chain A) + 2 (chain B) = 7
        assert len(lines) == 7

    def test_csv_no_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.psa"
        csv_path = tmp_path / "out.csv"
        csv_path.write_text("old")
        seqs = _make_sequences()
        with _open(path, mode="w") as f:
            f.write_sequence(sequences=seqs)
            with pytest.raises(FileExistsError):
                f.write_sequence_to_csv(csv_path, overwrite=False)
