"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for the RMSF-specific PTAFile API:
- write_rmsf_data / read_rmsf_data round-trip
- _iter_selections helper
- write_rmsf_statistics (incl. argmax)
- CSV / TSV exporters (non-merged + merged schema)
- Reader fallback to merged-record schema (mean/std/n)
- Error paths: empty input, duplicate labels, mismatched array lengths,
  overwrite semantics
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import (  # noqa: E402
    DEFAULT_RMSF_SELECTIONS,
    build_rmsf_pta,
    build_rmsf_merged_pta,
)

from pharmacon.fileio.pta import PharmaconPTAFile  # noqa: E402


class TestRoundTrip:
    def test_read_returns_dict_keyed_by_label(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        assert set(data) == {"calpha", "backbone"}

    def test_per_atom_values_round_trip(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        assert data["calpha"]["rmsf"].tolist() == [0.5, 1.2, 0.8]
        assert data["calpha"]["resids"].tolist() == [1, 2, 3]
        assert data["calpha"]["resnames"].tolist() == ["ALA", "GLY", "SER"]
        assert data["calpha"]["atom_names"].tolist() == ["CA", "CA", "CA"]
        assert data["calpha"]["atom_indices"].tolist() == [0, 4, 9]

    def test_selection_string_preserved(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        assert data["calpha"]["selection_string"] == "name CA"
        assert data["backbone"]["selection_string"] == "backbone"

    def test_numeric_fields_keep_their_dtypes(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        assert np.issubdtype(data["calpha"]["atom_indices"].dtype, np.integer)
        assert np.issubdtype(data["calpha"]["resids"].dtype, np.integer)
        assert np.issubdtype(data["calpha"]["rmsf"].dtype, np.floating)


class TestIterSelections:
    def test_yields_sorted_labels(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            labels = list(pta._iter_selections("rmsf"))
        # alphabetical
        assert labels == ["backbone", "calpha"]

    def test_skips_statistics_subgroup(self, tmp_path):
        # build_rmsf_pta also writes /rmsf/statistics — must not appear here
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            labels = list(pta._iter_selections("rmsf"))
        assert "statistics" not in labels


class TestStatistics:
    def test_statistics_table_exists(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            assert "rmsf/statistics/table" in pta.file

    def test_mean_min_max_correct(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            dset = pta.file["rmsf/statistics/table"]
            recs = {json.loads(row)["label"]: json.loads(row) for row in dset}

        calpha_vals = [0.5, 1.2, 0.8]
        r = recs["calpha"]
        assert r["n_atoms"] == 3
        assert pytest.approx(r["mean"]) == sum(calpha_vals) / 3
        assert r["min"] == min(calpha_vals)
        assert r["max"] == max(calpha_vals)
        assert pytest.approx(r["median"]) == sorted(calpha_vals)[1]

    def test_argmax_identifies_most_flexible_atom(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            dset = pta.file["rmsf/statistics/table"]
            recs = {json.loads(row)["label"]: json.loads(row) for row in dset}

        # Both selections have GLY:CA = 1.2 as max
        for label in ("calpha", "backbone"):
            assert recs[label]["argmax_resname"] == "GLY"
            assert recs[label]["argmax_atom_name"] == "CA"
            assert recs[label]["argmax_resid"] == 2

    def test_fallback_to_mean_on_merged_records(self, tmp_path):
        # write_rmsf_statistics should not crash on a merged file — it
        # falls back to rec["mean"] instead of rec["rmsf"]. Stats are
        # computed inside the build session because PharmaconPTAFile
        # refuses to reopen an existing file in mode='a'.
        pta_path = build_rmsf_merged_pta(
            tmp_path / "merged.pta", compute_statistics=True,
        )
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            dset = pta.file["rmsf/statistics/table"]
            recs = [json.loads(row) for row in dset]
        # Same shape as non-merged stats
        assert {r["label"] for r in recs} == {"calpha", "backbone"}


class TestMergedReadback:
    def test_read_returns_mean_std_n_for_merged_files(self, tmp_path):
        pta_path = build_rmsf_merged_pta(tmp_path / "merged.pta", std_value=0.15)
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        entry = data["calpha"]
        assert "rmsf" not in entry
        assert set(entry) >= {"mean", "std", "n", "atom_indices", "resids",
                              "resnames", "atom_names", "selection_string"}
        assert entry["mean"].tolist() == [0.5, 1.2, 0.8]
        assert entry["std"].tolist() == [0.15, 0.15, 0.15]
        assert entry["n"].tolist() == [2, 2, 2]


class TestExporters:
    def test_csv_non_merged_columns_and_row_count(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        out = tmp_path / "data.csv"
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            pta.write_rmsf_data_to_csv(out, group_name="rmsf", is_merged=False)
        with open(out) as fh:
            rows = list(csv.reader(fh))
        header = rows[0]
        assert header == ["selection", "atom_index", "resid", "resname",
                          "atom_name", "rmsf"]
        # Total atoms across both selections = 3 + 6 = 9
        assert len(rows) - 1 == 9

    def test_csv_merged_columns(self, tmp_path):
        pta_path = build_rmsf_merged_pta(tmp_path / "merged.pta")
        out = tmp_path / "data.csv"
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            pta.write_rmsf_data_to_csv(out, group_name="rmsf", is_merged=True)
        with open(out) as fh:
            header = next(csv.reader(fh))
        assert header == ["selection", "atom_index", "resid", "resname",
                          "atom_name", "mean", "std"]

    def test_tsv_uses_tabs_and_same_schema(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        out = tmp_path / "data.tsv"
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            pta.write_rmsf_data_to_tsv(out, group_name="rmsf", is_merged=False)
        text = out.read_text()
        assert "\t" in text
        header = text.splitlines()[0].split("\t")
        assert header == ["selection", "atom_index", "resid", "resname",
                          "atom_name", "rmsf"]

    def test_statistics_csv_one_row_per_selection(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        out = tmp_path / "stats.csv"
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            pta.write_rmsf_statistics_to_csv(out, group_name="rmsf")
        with open(out) as fh:
            rows = list(csv.reader(fh))
        # header + one row per selection
        assert len(rows) == 1 + 2
        assert rows[0][0] == "selection"
        labels = {row[0] for row in rows[1:]}
        assert labels == {"calpha", "backbone"}

    def test_atom_indices_round_trip_as_ints_in_csv(self, tmp_path):
        # Regression: merged-engine output also preserved as int strings,
        # not booleans / floats / quoted strings.
        pta_path = build_rmsf_merged_pta(tmp_path / "merged.pta")
        out = tmp_path / "data.csv"
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            pta.write_rmsf_data_to_csv(out, group_name="rmsf", is_merged=True)
        with open(out) as fh:
            rows = list(csv.reader(fh))[1:]
        for r in rows:
            # atom_index column is index 1; must look like an int (no quotes,
            # no decimal, no "True"/"False")
            int(r[1])  # raises if not


class TestErrors:
    def _open(self, path):
        return PharmaconPTAFile(path, overwrite=True, mode="a",
                                command="Trajectory Analysis", subcommand="rmsf")

    def test_empty_selections_raises(self, tmp_path):
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            with pytest.raises(ValueError, match="No selections"):
                pta.write_rmsf_data(
                    selections=[],
                    group_name="rmsf",
                    frame_begin=0, frame_end=10, frame_step=1,
                    fitting_group="name CA",
                )

    def test_duplicate_labels_raises(self, tmp_path):
        sel = {
            "label": "dup",
            "selection_string": "name CA",
            "atom_indices": np.array([0], dtype=int),
            "resids":       np.array([1], dtype=int),
            "resnames":     np.array(["ALA"]),
            "atom_names":   np.array(["CA"]),
            "rmsf":         np.array([0.5]),
        }
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            with pytest.raises(ValueError, match="Duplicate"):
                pta.write_rmsf_data(
                    selections=[sel, dict(sel)],
                    group_name="rmsf",
                    frame_begin=0, frame_end=10, frame_step=1,
                    fitting_group="name CA",
                )

    def test_mismatched_array_lengths_raises(self, tmp_path):
        bad = {
            "label": "x",
            "selection_string": "name CA",
            "atom_indices": np.array([0, 1, 2], dtype=int),
            "resids":       np.array([1, 2], dtype=int),   # too short
            "resnames":     np.array(["ALA", "GLY", "SER"]),
            "atom_names":   np.array(["CA", "CA", "CA"]),
            "rmsf":         np.array([0.5, 0.6, 0.7]),
        }
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            with pytest.raises(ValueError, match="resids"):
                pta.write_rmsf_data(
                    selections=[bad],
                    group_name="rmsf",
                    frame_begin=0, frame_end=10, frame_step=1,
                    fitting_group="name CA",
                )

    def test_missing_required_key_raises(self, tmp_path):
        bad = {
            "label": "x",
            "selection_string": "name CA",
            "atom_indices": np.array([0], dtype=int),
            # missing resids, resnames, atom_names, rmsf
        }
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            with pytest.raises(ValueError, match="missing"):
                pta.write_rmsf_data(
                    selections=[bad],
                    group_name="rmsf",
                    frame_begin=0, frame_end=10, frame_step=1,
                    fitting_group="name CA",
                )

    def test_overwrite_false_raises_on_existing_selection(self, tmp_path):
        sel = {
            "label": "calpha",
            "selection_string": "name CA",
            "atom_indices": np.array([0], dtype=int),
            "resids":       np.array([1], dtype=int),
            "resnames":     np.array(["ALA"]),
            "atom_names":   np.array(["CA"]),
            "rmsf":         np.array([0.5]),
        }
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            pta.write_rmsf_data(
                selections=[sel],
                group_name="rmsf",
                frame_begin=0, frame_end=10, frame_step=1,
                fitting_group="name CA",
            )
            with pytest.raises(FileExistsError):
                pta.write_rmsf_data(
                    selections=[sel],
                    group_name="rmsf",
                    frame_begin=0, frame_end=10, frame_step=1,
                    fitting_group="name CA",
                    overwrite=False,
                )

    def test_overwrite_true_replaces_selection(self, tmp_path):
        sel1 = {
            "label": "calpha", "selection_string": "name CA",
            "atom_indices": np.array([0], dtype=int),
            "resids":       np.array([1], dtype=int),
            "resnames":     np.array(["ALA"]),
            "atom_names":   np.array(["CA"]),
            "rmsf":         np.array([0.5]),
        }
        sel2 = dict(sel1, rmsf=np.array([9.99]))  # different value
        with self._open(tmp_path / "rmsf.pta") as pta:
            pta.create_group("rmsf")
            pta.write_rmsf_data(
                selections=[sel1], group_name="rmsf",
                frame_begin=0, frame_end=10, frame_step=1,
                fitting_group="name CA",
            )
            pta.write_rmsf_data(
                selections=[sel2], group_name="rmsf",
                frame_begin=0, frame_end=10, frame_step=1,
                fitting_group="name CA",
                overwrite=True,
            )
            data = pta.read_rmsf_data(group_name="rmsf")
        assert data["calpha"]["rmsf"].tolist() == [9.99]

    def test_read_on_missing_group_raises(self, tmp_path):
        pta_path = build_rmsf_pta(tmp_path / "rmsf.pta")
        with PharmaconPTAFile(pta_path, mode="r") as pta:
            with pytest.raises(ValueError, match="does not exist"):
                pta.read_rmsf_data(group_name="nope")
