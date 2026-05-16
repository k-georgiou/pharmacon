"""Tests for command_line.plot.pta.validate() against a mock PTA.

Verifies the happy path and the five tamper/missing-metadata failure modes.
"""

from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path

import h5py
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_pli_pta

from pharmacon.command_line.plot.pta import validate
from pharmacon.command_line.exceptions import ValidationError


def _make_args(*, input_pta: Path, output_dir: Path, log_path: Path,
               config: Path | None = None, overwrite: bool = False) -> Namespace:
    return Namespace(
        input=str(input_pta),
        output=str(output_dir),
        log=str(log_path),
        config=str(config) if config else None,
        overwrite=overwrite,
        maxwarnings=0,
        file_logging_level="TRACE",
        terminal_logging_level="INFO",
    )


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestHappyPath:
    def test_clean_mock_pta_passes(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=50)
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        validate(args)  # should not raise
        # Side-effects: args.input/output coerced to Path; config_overrides set
        assert isinstance(args.input, Path)
        assert isinstance(args.output, Path)
        assert args.config is None
        assert args.config_overrides == {}

    def test_with_example_config(self, tmp_path):
        repo_root = Path(__file__).resolve().parent.parent
        cfg = repo_root / "examples" / "plot_ini" / "pli_stacked_column_1.ini"
        assert cfg.exists()

        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=50)
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
            config=cfg,
        )
        validate(args)
        assert isinstance(args.config_overrides, dict)
        # At least one section parsed from the INI
        assert len(args.config_overrides) >= 1


# ---------------------------------------------------------------------------
# Failure paths — corrupt the PTA in distinct ways
# ---------------------------------------------------------------------------

def _open_attrs(path: Path):
    return h5py.File(path, "a")


class TestValidationFailures:
    def test_signature_mismatch_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            f.attrs["signature"] = "PHARMACON::TAMPERED::WITH::DEAD-BEEF"
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="signature mismatch"):
            validate(args)

    def test_missing_blueprint_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            del f.attrs["blueprint"]
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="missing 'blueprint'"):
            validate(args)

    def test_bad_artifact_token_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            f.attrs["artifact_token"] = "AAAAA-DEFINITELY-NOT-VALID"
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="artifact_token does not match"):
            validate(args)

    def test_missing_version_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            del f.attrs["pharmacon_version"]
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="missing 'pharmacon_version'"):
            validate(args)

    def test_incomplete_group_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            f["pl_interactions"].attrs["completed"] = "False"
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="does not have completed=True"):
            validate(args)

    def test_artifact_status_not_success_raises(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            f.attrs["artifact_status"] = "FAILED"
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="artifact_status is 'FAILED'"):
            validate(args)


# ---------------------------------------------------------------------------
# Version-ordering: file < runtime passes; file > runtime fails clearly
# ---------------------------------------------------------------------------

class TestVersionOrdering:
    def test_old_file_validates_against_newer_runtime(self, tmp_path):
        """A PTA written by a previous Pharmacon release must still be
        readable by the current (newer) runtime."""
        pta = build_pli_pta(
            tmp_path / "pli_old.pta", n_frames=10,
            pharmacon_version="0.1.0",
        )
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        validate(args)  # must not raise

    def test_newer_file_rejected_with_clear_message(self, tmp_path):
        """A PTA written by a future Pharmacon must fail with a clear
        'upgrade Pharmacon' error, not a confusing signature mismatch."""
        pta = build_pli_pta(
            tmp_path / "pli_new.pta", n_frames=10,
            pharmacon_version="9.9.9",
        )
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match=r"requires Pharmacon >= 9\.9\.9"):
            validate(args)

    def test_tampered_version_breaks_signature(self, tmp_path):
        """If somebody edits pharmacon_version without re-signing the
        file, the signature regenerated under the new (claimed) version
        will not match the stored one."""
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=10)
        with _open_attrs(pta) as f:
            # Lower-than-runtime so we don't trip the version gate first
            f.attrs["pharmacon_version"] = "0.5.0"
        args = _make_args(
            input_pta=pta,
            output_dir=tmp_path / "out",
            log_path=tmp_path / "plot.log",
        )
        with pytest.raises(ValidationError, match="signature mismatch"):
            validate(args)
