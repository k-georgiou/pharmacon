"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Shared integrity validation for Pharmacon PTA/PSA files.

This module centralises the on-disk metadata checks every consumer of a
Pharmacon analysis file performs before reading data — `plot pta`,
`dump pta`/`dump psa`, `export pta`/`export psa`, and `merge`. Pulling the
logic here keeps the checks consistent across all call sites and prevents
the drift that occurred when each command-line module carried its own copy.

The validator opens the file via h5py (not via :class:`PharmaconHDF5File`)
to avoid triggering the same checks while the file is being read and to
stay decoupled from the file-handle lifecycle of the caller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final, Dict

from pharmacon.command_line.exceptions import ValidationError


__all__ = [
    "validate_pharmacon_file",
]


# Per-format token secret used to verify artifact tokens.
_FORMAT_TOKEN_SECRET: Final[Dict[str, str]] = {
    "pta": "trajectory_analysis",
    "psa": "structure_analysis",
}


def _parse_version(s: str) -> tuple[int, ...]:
    """Parse "X.Y.Z" style version strings as int tuples; empty tuple on failure."""
    try:
        return tuple(int(p) for p in s.split("."))
    except (ValueError, AttributeError):
        return ()


def validate_pharmacon_file(
    path: Path,
    *,
    expected_format: str,
    allow_merged: bool = True,
) -> Dict[str, str]:
    """Open *path* and run the full Pharmacon-file integrity check.

    Performs the following in order; any failure raises
    :class:`ValidationError` with a message that names the offending file
    and field:

    1. File is openable via h5py.
    2. ``artifact_status`` attr is present and equals ``"SUCCESS"``.
    3. ``artifact_token`` and ``blueprint`` attrs are present, and the
       token matches the blueprint under the format's secret.
    4. ``signature`` and ``fingerprint`` attrs are present.
    5. ``pharmacon_version`` is present and **not greater than** the
       current runtime version (older files are accepted).
    6. The stored ``signature``/``fingerprint`` match what
       :func:`create_pharmacon_signature` produces under the file's own
       version (this stays stable across runtime upgrades and only fails
       on tampering).
    7. Every top-level group in the file has ``completed=True``.
    8. If *allow_merged* is False, ``is_merged`` must be ``False``
       (merge refuses to re-merge already-merged inputs).

    :param path: Path to the file.
    :param expected_format: ``"pta"`` or ``"psa"``. Determines the token
        secret and the human-readable prefix in error messages.
    :param allow_merged: If False, an ``is_merged=True`` file raises.
    :return: A dict of the file's validated attrs as ``str`` values, so
        callers can use values like ``pharmacon_version``, ``command``,
        ``subcommand``, ``blueprint``, ``is_merged`` without re-reading.
    :raises ValidationError: On any check failure.
    """

    import h5py

    from pharmacon.constants import __version__
    from pharmacon.utils.fingerprint import create_pharmacon_signature
    from pharmacon.utils.identifiers import validate_mda_artifact_token

    fmt = expected_format.strip().lower()
    if fmt not in _FORMAT_TOKEN_SECRET:
        raise ValueError(
            f"expected_format must be one of {sorted(_FORMAT_TOKEN_SECRET)}, "
            f"got {expected_format!r}"
        )
    label = fmt.upper()  # "PTA" / "PSA"
    secret = _FORMAT_TOKEN_SECRET[fmt]

    try:
        f = h5py.File(path, "r")
    except Exception as exc:
        raise ValidationError(f"Cannot open {label} file '{path}': {exc}") from exc

    try:
        attrs: Dict[str, str] = {
            str(k): str(v) for k, v in f.attrs.items()
        }

        # (2) artifact_status
        artifact_status = attrs.get("artifact_status", "").strip().upper()
        if not artifact_status:
            raise ValidationError(
                f"{label} file '{path}' is missing 'artifact_status' metadata — "
                "the file may be incomplete or corrupted."
            )
        if artifact_status != "SUCCESS":
            raise ValidationError(
                f"{label} file '{path}' artifact_status is "
                f"'{artifact_status}' (expected 'SUCCESS') — the file may be "
                "corrupted or the analysis did not complete."
            )

        # (3) artifact_token + blueprint
        artifact_token = attrs.get("artifact_token", "").strip()
        if not artifact_token:
            raise ValidationError(
                f"{label} file '{path}' is missing 'artifact_token' metadata."
            )
        blueprint = attrs.get("blueprint", "").strip()
        if not blueprint:
            raise ValidationError(
                f"{label} file '{path}' is missing 'blueprint' metadata."
            )
        token_valid = validate_mda_artifact_token(
            artifact_token=artifact_token,
            blueprint=blueprint,
            secret=secret,
            namespace="pharmacon",
        )
        if not token_valid:
            raise ValidationError(
                f"{label} file '{path}' artifact_token does not match the "
                "blueprint — the file may have been tampered with or corrupted."
            )

        # (4) signature + fingerprint present
        signature = attrs.get("signature", "").strip()
        fingerprint = attrs.get("fingerprint", "").strip()
        if not signature or not fingerprint:
            raise ValidationError(
                f"{label} file '{path}' is missing 'signature' and/or "
                "'fingerprint' metadata."
            )

        # (5) pharmacon_version present + ordering
        file_version = attrs.get("pharmacon_version", "").strip()
        if not file_version:
            raise ValidationError(
                f"{label} file '{path}' is missing 'pharmacon_version' metadata."
            )
        runtime_version = str(__version__)
        f_parts = _parse_version(file_version)
        r_parts = _parse_version(runtime_version)
        if f_parts and r_parts and f_parts > r_parts:
            raise ValidationError(
                f"{label} file '{path}' requires Pharmacon >= {file_version}; "
                f"runtime is {runtime_version}. Please upgrade Pharmacon."
            )

        # (6) signature/fingerprint match under file's own version
        command = attrs.get("command", "").strip()
        subcommand = attrs.get("subcommand", "").strip()
        if command and subcommand:
            expected_sig = create_pharmacon_signature(
                format_name=fmt,
                command=command,
                subcommand=subcommand,
                version=file_version,
            )
            if expected_sig.signature != signature:
                raise ValidationError(
                    f"{label} file '{path}' signature mismatch.\n"
                    f"  Expected : {expected_sig.signature}\n"
                    f"  Found    : {signature}\n"
                    "The file may have been tampered with or corrupted."
                )
            if expected_sig.fingerprint != fingerprint:
                raise ValidationError(
                    f"{label} file '{path}' fingerprint mismatch.\n"
                    f"  Expected : {expected_sig.fingerprint}\n"
                    f"  Found    : {fingerprint}\n"
                    "The file may have been tampered with or corrupted."
                )

        # (7) every top-level group must have completed=True
        for key in f:
            obj = f[key]
            if not hasattr(obj, "attrs") or not hasattr(obj, "keys"):
                continue
            completed = str(obj.attrs.get("completed", "")).strip().lower()
            if completed != "true":
                raise ValidationError(
                    f"{label} file '{path}': group '{key}' does not have "
                    "completed=True — the analysis for this group may not "
                    "have finished."
                )

        # (8) is_merged guard (merge-specific)
        is_merged = attrs.get("is_merged", "False").strip().lower() == "true"
        if (not allow_merged) and is_merged:
            raise ValidationError(
                f"{label} file '{path}' is already a merged file "
                "(is_merged=True); refusing to re-merge."
            )

        attrs["is_merged"] = "True" if is_merged else "False"
        return attrs

    finally:
        f.close()
