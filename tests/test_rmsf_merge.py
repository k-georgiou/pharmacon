"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for the RMSF merge pipeline (`merge results` subcommand).

Drives `merge_results.run()` directly, bypassing `validate()` so the
test doesn't need to forge a valid artifact token. The Namespace passed
in mirrors what the CLI argparse layer would build.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import (  # noqa: E402
    DEFAULT_RMSF_SELECTIONS,
    build_rmsf_pta,
)

from pharmacon.fileio.pta import PharmaconPTAFile  # noqa: E402
from pharmacon.command_line.merge import results as merge_results  # noqa: E402


def _run_merge(inputs, output, *, tmp_path, max_warnings=100,
               file_level="WARNING", terminal_level="WARNING"):
    """Build the Namespace `merge_results.run` expects and execute it."""
    ns = argparse.Namespace(
        input=list(inputs),
        output=output,
        overwrite=True,
        max_warnings=max_warnings,
        warn_only=False,
        log=tmp_path / "merge.log",
        file_logging_level=file_level,
        terminal_logging_level=terminal_level,
    )
    merge_results.run(ns)
    return output


class TestHappyPath:
    def test_merged_file_is_created(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        assert out.is_file() and out.stat().st_size > 0

    def test_merged_file_has_rmsf_group_with_selections(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            assert "rmsf" in pta.file
            labels = list(pta._iter_selections("rmsf"))
        assert set(labels) == {"calpha", "backbone"}

    def test_per_atom_records_carry_mean_std_n(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        assert "mean" in data["calpha"]
        assert "std" in data["calpha"]
        assert "n" in data["calpha"]
        # Both inputs identical → mean equals the original rmsf, std = 0
        assert data["calpha"]["mean"].tolist() == [0.5, 1.2, 0.8]
        np.testing.assert_allclose(data["calpha"]["std"], 0.0, atol=1e-12)
        assert data["calpha"]["n"].tolist() == [2, 2, 2]


class TestArithmetic:
    def test_mean_std_across_two_distinct_inputs(self, tmp_path):
        # Same atom keys in both inputs but different RMSF values.
        sel_a = (("calpha", "name CA", [
            (0, 1, "ALA", "CA", 0.4),
            (4, 2, "GLY", "CA", 1.0),
        ]),)
        sel_b = (("calpha", "name CA", [
            (0, 1, "ALA", "CA", 0.6),
            (4, 2, "GLY", "CA", 1.4),
        ]),)
        a = build_rmsf_pta(tmp_path / "a.pta", selections=sel_a)
        b = build_rmsf_pta(tmp_path / "b.pta", selections=sel_b)
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")

        # Mean = (a + b) / 2; std with ddof=0
        np.testing.assert_allclose(data["calpha"]["mean"], [0.5, 1.2])
        np.testing.assert_allclose(data["calpha"]["std"], [0.1, 0.2])
        assert data["calpha"]["n"].tolist() == [2, 2]


class TestTypePreservation:
    def test_atom_index_and_resid_remain_ints(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            dset = pta.file["rmsf/selection_calpha/atoms"]
            rec = json.loads(dset[0])
        assert isinstance(rec["atom_index"], int)
        assert isinstance(rec["resid"], int)
        assert isinstance(rec["resname"], str)
        assert isinstance(rec["atom_name"], str)


class TestAttrPropagation:
    def test_fitting_group_and_frame_range_carry_over(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta",
                           fitting_group="protein and name CA",
                           frame_begin=10, frame_end=100, frame_step=5)
        b = build_rmsf_pta(tmp_path / "b.pta",
                           fitting_group="protein and name CA",
                           frame_begin=10, frame_end=100, frame_step=5)
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            attrs = dict(pta.file["rmsf/selection_calpha/atoms"].attrs)
        assert attrs["fitting_group"] == "protein and name CA"
        assert attrs["frame_begin"] == "10"
        assert attrs["frame_end"] == "100"
        assert attrs["frame_step"] == "5"
        assert attrs["merged"] == "True"

    def test_n_atoms_reflects_merged_count(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            attrs = dict(pta.file["rmsf/selection_calpha/atoms"].attrs)
        # All 3 calpha atoms are common to both inputs → n_atoms stays 3
        assert attrs["n_atoms"] == "3"


class TestKeyIntersection:
    def test_atom_only_in_one_input_is_dropped(self, tmp_path):
        # B has an extra atom that A lacks → must be dropped
        sel_a = (("calpha", "name CA", [
            (0, 1, "ALA", "CA", 0.5),
            (4, 2, "GLY", "CA", 1.2),
        ]),)
        sel_b = (("calpha", "name CA", [
            (0, 1, "ALA", "CA", 0.6),
            (4, 2, "GLY", "CA", 1.0),
            (9, 3, "SER", "CA", 0.8),  # extra
        ]),)
        a = build_rmsf_pta(tmp_path / "a.pta", selections=sel_a)
        b = build_rmsf_pta(tmp_path / "b.pta", selections=sel_b)
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            data = pta.read_rmsf_data(group_name="rmsf")
        # Only the 2 shared atoms survive
        assert len(data["calpha"]["atom_indices"]) == 2
        assert set(data["calpha"]["atom_indices"].tolist()) == {0, 4}


class TestFileMetadata:
    def test_merged_marker_set(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            attrs = dict(pta.file.attrs)
        assert attrs.get("is_merged") == "True"
        assert attrs.get("n_inputs") == "2"

    def test_per_input_reference_groups_written(self, tmp_path):
        a = build_rmsf_pta(tmp_path / "a.pta")
        b = build_rmsf_pta(tmp_path / "b.pta")
        out = tmp_path / "merged.pta"
        _run_merge([a, b], out, tmp_path=tmp_path)
        with PharmaconPTAFile(out, mode="r") as pta:
            assert "file1_metadata" in pta.file
            assert "file2_metadata" in pta.file
