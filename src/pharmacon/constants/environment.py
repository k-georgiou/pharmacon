"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Runtime environment helpers for Pharmacon.

Provides safe accessors for environment-controlled paths used by the
job infrastructure (scratch/tmp directories, per-job working directories).

Environment variables
---------------------
PHARMACON_TMPDIR
    Root scratch directory for all Pharmacon job working directories.

    If unset or empty, a safe per-user fallback directory is created under
    the system temporary directory and exported automatically, e.g.::

        /tmp/pharmacon/<username>

    You may still override it explicitly, for example::

        export PHARMACON_TMPDIR=/scratch/$USER/pharmacon

Public API
----------
safe_get_pharmacon_tmpdir() -> Path
    Resolve, initialize if needed, and validate PHARMACON_TMPDIR.

safe_create_jobdir(prefix, username) -> Path
    Atomically create a uniquely-named job directory under PHARMACON_TMPDIR.
"""

from __future__ import annotations

import getpass
import os
import stat
import tempfile
from pathlib import Path



__all__ = [
    "EnvironmentConfigError",
    "safe_get_pharmacon_tmpdir",
    "safe_create_jobdir",
]


_ENV_VAR = "PHARMACON_TMPDIR"


# Exceptions


class EnvironmentConfigError(RuntimeError):
    """Raised when a required environment variable is missing or invalid."""


# Internal helpers


def _safe_default_tmpdir() -> Path:
    """
    Constructs a temporary directory path specific to the current user. The path
    is built using the system's temporary directory combined with a predefined
    subdirectory named "pharmacon" and the username of the current user. If the
    username cannot be retrieved, a default value of "unknown" is used instead.

    :return: A `Path` object representing the constructed temporary directory
        path.
    :rtype: Path
    """
    try:
        username = getpass.getuser().strip()
    except Exception:
        username = ""

    if not username:
        username = "unknown"

    return (Path(tempfile.gettempdir()) / "pharmacon" / username).expanduser()


def _ensure_directory(path: Path) -> None:
    """
    Ensures the provided directory path exists, is a directory, and has appropriate
    permissions applied. If the directory does not exist, it will be created with the
    specified mode. Permissions are set on a best-effort basis and do not guarantee
    writability but are checked subsequently.

    :param path: The directory path to ensure exists and is correctly configured.
    :type path: Path
    :raises EnvironmentConfigError: If the directory cannot be created, doesn't exist
        post-creation, or is not a directory.
    """
    try:
        path.mkdir(mode=0o700, parents=True, exist_ok=True)
    except OSError as exc:
        raise EnvironmentConfigError(
            f"Failed to create directory {path}: {exc}"
        ) from exc

    if not path.exists():
        raise EnvironmentConfigError(f"Directory does not exist after creation: {path}")

    if not path.is_dir():
        raise EnvironmentConfigError(f"Path is not a directory: {path}")

    try:
        path.chmod(stat.S_IRWXU)
    except OSError:
        # Best-effort only; validation below still checks writability.
        pass


# Public API
def safe_get_pharmacon_tmpdir() -> Path:
    """
    Retrieves a temporary directory path for the application, ensuring it is valid, writable,
    and exists. It tries to fetch the directory path from an environment variable, falling
    back to a default location if the variable is not set. Directory permissions and existence
    are verified, and necessary directories will be created if they do not already exist.

    If ``PHARMACON_TMPDIR`` is unset or empty, a safe per-user fallback directory
    under the system temp directory is created and exported into the process
    environment automatically.

    :raises EnvironmentConfigError: If the directory path is not valid, doesn't exist and
        cannot be created, is not a directory, or is not writable by the current user.
    :return: A `Path` object representing the resolved temporary directory.
    :rtype: Path
    """
    raw = os.environ.get(_ENV_VAR, "").strip()

    if not raw:
        fallback = _safe_default_tmpdir()
        _ensure_directory(fallback)
        os.environ[_ENV_VAR] = str(fallback)
        raw = str(fallback)

    try:
        tmpdir = Path(raw).expanduser().resolve(strict=False)
    except Exception as exc:
        raise EnvironmentConfigError(
            f"{_ENV_VAR}={raw!r} is not a valid path: {exc}"
        ) from exc

    if not tmpdir.exists():
        try:
            tmpdir.mkdir(mode=0o700, parents=True, exist_ok=True)
        except OSError as exc:
            raise EnvironmentConfigError(
                f"{_ENV_VAR} directory does not exist and could not be created: {tmpdir}\n"
                f"Create it first or choose another directory.\n"
                f"Original error: {exc}"
            ) from exc

    if not tmpdir.is_dir():
        raise EnvironmentConfigError(
            f"{_ENV_VAR} path is not a directory: {tmpdir}"
        )

    if not os.access(tmpdir, os.W_OK | os.X_OK):
        raise EnvironmentConfigError(
            f"{_ENV_VAR} directory is not writable by the current user: {tmpdir}"
        )

    return tmpdir


def safe_create_jobdir(prefix: str, username: str) -> Path:
    """
    Creates a temporary job directory with a specific prefix and username.

    This function securely generates a unique directory for a job in the temporary
    directory of the system. The prefix and username values are validated to ensure
    they meet the criteria for naming. It also ensures proper permission settings
    for the created directory.

    The directory name follows the pattern::

        <PHARMACON_TMPDIR>/<username>_<prefix>_<random-suffix>/

    The random suffix is generated by :func:`tempfile.mkdtemp` which is
    guaranteed to be unique and creates the directory atomically, preventing
    race conditions even between concurrent MPI ranks that call this function
    simultaneously.

    Permissions are set to ``0o700`` (owner read/write/execute only).

    :param prefix: A non-empty alphanumeric string (hyphens/underscores are allowed)
        used as part of the directory name.
    :param username: A non-empty string representing the username to be included
        in the directory name.
    :return: The path to the newly created job directory as a ``Path`` object.
    :rtype: Path
    :raises ValueError: If the `prefix` is not a valid alphanumeric string or
        if `username` is an empty string.
    """
    if not prefix or not prefix.replace("-", "").replace("_", "").isalnum():
        raise ValueError(
            f"prefix must be a non-empty alphanumeric string (hyphens/underscores allowed); "
            f"got {prefix!r}."
        )
    if not username:
        raise ValueError("username must be a non-empty string.")

    tmpdir = safe_get_pharmacon_tmpdir()

    # mkdtemp creates the directory atomically.
    dir_prefix = f"{username}_{prefix}_"
    job_path = Path(tempfile.mkdtemp(prefix=dir_prefix, dir=tmpdir))

    # Ensure permissions even if umask was permissive.
    job_path.chmod(stat.S_IRWXU)

    return job_path
