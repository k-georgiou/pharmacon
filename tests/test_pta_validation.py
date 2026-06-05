"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for the shared validator: pharmacon.utils.pta_validation.

Covers happy paths and all named failure modes for both PTA and PSA
files, plus the merge-specific `allow_merged=False` guard.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_pli_pta

from pharmacon.utils.pta_validation import validate_pharmacon_file
from pharmacon.command_line.exceptions import ValidationError


class TestHappyPaths:
    def test_pta_returns_attrs_dict(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        attrs = validate_pharmacon_file(path, expected_format="pta")
        assert attrs["command"] == "Trajectory Analysis"
        assert attrs["subcommand"] == "pl-interactions"
        assert attrs["pharmacon_version"] == "1.0.0"
        assert attrs["is_merged"] == "False"

    def test_pta_old_version_allowed(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli_old.pta", n_frames=10,
                             pharmacon_version="0.1.0")
        attrs = validate_pharmacon_file(path, expected_format="pta")
        assert attrs["pharmacon_version"] == "0.1.0"


class TestBadFormatKw:
    def test_unknown_format_raises_value_error(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with pytest.raises(ValueError, match="expected_format"):
            validate_pharmacon_file(path, expected_format="xyz")


def _open_attrs(path: Path):
    return h5py.File(path, "a")


class TestFailurePaths:
    def test_unopenable_file_raises(self, tmp_path):
        with pytest.raises(ValidationError, match="Cannot open PTA file"):
            validate_pharmacon_file(tmp_path / "nope.pta", expected_format="pta")

    def test_artifact_status_failed_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            f.attrs["artifact_status"] = "FAILED"
        with pytest.raises(ValidationError, match="artifact_status is 'FAILED'"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_missing_blueprint_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            del f.attrs["blueprint"]
        with pytest.raises(ValidationError, match="missing 'blueprint'"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_bad_token_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            f.attrs["artifact_token"] = "INVALID-TOKEN"
        with pytest.raises(ValidationError, match="artifact_token does not match"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_missing_signature_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            del f.attrs["signature"]
        with pytest.raises(ValidationError, match="'signature' and/or 'fingerprint'"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_missing_version_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            del f.attrs["pharmacon_version"]
        with pytest.raises(ValidationError, match="missing 'pharmacon_version'"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_future_version_raises_with_clear_message(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10,
                             pharmacon_version="9.9.9")
        with pytest.raises(ValidationError, match=r"requires Pharmacon >= 9\.9\.9"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_tampered_signature_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            f.attrs["signature"] = "PHARMACON::TAMPERED::WITH::DEAD-BEEF"
        with pytest.raises(ValidationError, match="signature mismatch"):
            validate_pharmacon_file(path, expected_format="pta")

    def test_incomplete_group_raises(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            f["pl_interactions"].attrs["completed"] = "False"
        with pytest.raises(ValidationError, match="does not have completed=True"):
            validate_pharmacon_file(path, expected_format="pta")


class TestAllowMergedGuard:
    def test_merged_file_passes_when_allowed(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        # Flip is_merged and re-sign nothing (is_merged is not in the signed
        # payload so the integrity check should still pass).
        with _open_attrs(path) as f:
            f.attrs["is_merged"] = "True"
        attrs = validate_pharmacon_file(
            path, expected_format="pta", allow_merged=True,
        )
        assert attrs["is_merged"] == "True"

    def test_merged_file_rejected_when_disallowed(self, tmp_path):
        path = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(path) as f:
            f.attrs["is_merged"] = "True"
        with pytest.raises(ValidationError, match="already a merged file"):
            validate_pharmacon_file(
                path, expected_format="pta", allow_merged=False,
            )
