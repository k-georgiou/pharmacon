"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

MPI-safe file handler for structured plain-text log output.

Design
------
* ``FileFormatter`` produces timestamped, structured plain-text lines (no ANSI).
* ``MPISafeFileHandler`` supports two strategies for concurrent MPI processes:

  **per-rank files** (default, ``per_rank=True``)
      Each MPI rank appends ``.rank_N`` to the log filename:
      ``pharmacon.log`` → ``pharmacon.rank_0.log``, ``pharmacon.rank_1.log``, …
      No locking is required because each rank owns its file exclusively.

  **shared file with POSIX locking** (``per_rank=False``)
      All ranks write to one file.  Every ``emit()`` acquires an exclusive
      ``fcntl.flock`` lock, writes, flushes, then releases the lock.  This is
      safe on most POSIX shared filesystems (Lustre, GPFS, NFS with proper
      lock-manager setup).  On Windows a ``threading.Lock`` is used instead.

MPI rank detection
------------------
The rank is discovered from common runtime environment variables without
requiring ``mpi4py``.  ``mpi4py`` is used as a fallback when available.

Log line format::

    2024-01-15 10:30:15.123 | RANK 0 | INFO     | >> message text   | module:lineno

Usage::

    from pharmacon.logger.file import make_file_handler

    handler = make_file_handler("run.log", level=logging.DEBUG)
    logging.getLogger("pharmacon").addHandler(handler)
"""

from __future__ import annotations

import fcntl
import logging
import os
import threading
from pathlib import Path
from typing import IO



__all__ = [
    "get_mpi_rank",
    "get_mpi_size",
    "FileFormatter",
    "MPISafeFileHandler",
    "make_file_handler",
]


_PREFIX = ">> "

# MPI rank detection

_MPI_RANK_VARS = (
    "OMPI_COMM_WORLD_RANK",   # Open MPI
    "PMI_RANK",               # MPICH / MVAPICH
    "SLURM_PROCID",           # Slurm
    "MPI_LOCALRANKID",        # IBM Spectrum MPI
    "ALPS_APP_PE",            # Cray ALPS
)


def get_mpi_rank() -> int | None:
    """
    Determines the MPI rank of the current process.

    This function searches for MPI rank information from the environment variables
    defined in `_MPI_RANK_VARS`. If no valid rank is found in the environment, it
    attempts to utilize the `mpi4py` library to retrieve the MPI rank of the current
    process. If neither approach succeeds, it will return `None`.

    :raises ValueError: If an environment variable within `_MPI_RANK_VARS` is set
        but cannot be converted to an integer.
    :raises ImportError: If the `mpi4py` library is not installed and required.

    :return: The MPI rank as an integer if found, otherwise `None`.
    :rtype: int | None
    """
    for var in _MPI_RANK_VARS:
        val = os.environ.get(var)
        if val is not None:
            try:
                return int(val)
            except ValueError:
                continue
    try:
        from mpi4py import MPI  # type: ignore[import]
        return MPI.COMM_WORLD.Get_rank()
    except ImportError:
        pass
    return None


def get_mpi_size() -> int | None:
    """
    Retrieves the size of the MPI (Message Passing Interface) world, if available.

    This utility function checks various environment variables that might be set in different MPI
    implementations to determine the number of processes in the MPI communicator. If no relevant
    environment variable is found or is invalid, it attempts to use the `mpi4py` library to get the
    MPI communicator size. If `mpi4py` is not available and no environment variable provides
    the information, it returns `None`.

    :raises ValueError: If an environment variable is found but its value cannot be converted to an integer.

    :return: The size of the MPI world as an integer, or `None` if it cannot be determined.
    :rtype: int | None
    """
    for var in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "SLURM_NTASKS"):
        val = os.environ.get(var)
        if val is not None:
            try:
                return int(val)
            except ValueError:
                continue
    try:
        from mpi4py import MPI  # type: ignore[import]
        return MPI.COMM_WORLD.Get_size()
    except ImportError:
        pass
    return None


# Formatter


class FileFormatter(logging.Formatter):
    """
    Custom logging formatter class for formatting log messages with a specific structure.

    This class enhances log messages by adding components such as timestamps,
    log levels, rank information (if applicable), and exception details when
    present. The formatter is particularly useful for creating a consistent and
    readable log output.

        Produces lines of the form::

        2024-01-15 10:30:15.123 | RANK 0 | INFO     | >> message   | name:lineno

    :ivar _rank: Optional rank identifier that is included in the log messages if provided.
    :type _rank: int | None
    """

    _DATEFMT = "%Y-%m-%d %H:%M:%S"
    _LEVEL_WIDTH = 8   # "CRITICAL" is 8 chars

    def __init__(self, rank: int | None = None) -> None:
        """
        Initializes an instance with optional ranking information.

        This initializer sets up the rank tag based on the rank provided. If no rank
        is given, the tag will remain an empty string.

        :param rank: An optional integer representing the rank. If provided, this
            will be used to create a formatted rank tag. If None, no rank tag is
            set.
        """
        super().__init__(datefmt=self._DATEFMT)
        self._rank = rank
        self._rank_tag = f"RANK {rank:>3d}" if rank is not None else ""

    def format(self, record: logging.LogRecord) -> str:
        """
        Formats a logging record object into a string representation using a custom
        format, including timestamp, logging level, message, and other metadata.
        Optionally appends exception information if available in the record.

        :param record: The logging record to be formatted.
        :type record: logging.LogRecord
        :return: A formatted string representing the log record.
        :rtype: str
        """
        # Timestamp with milliseconds
        ts = self.formatTime(record, self._DATEFMT)
        ms = int(record.msecs)
        timestamp = f"{ts}.{ms:03d}"

        level = record.levelname.ljust(self._LEVEL_WIDTH)
        msg = record.getMessage()

        parts = [timestamp]
        if self._rank_tag:
            parts.append(self._rank_tag)
        parts.append(level)
        parts.append(f"{_PREFIX}{msg}")
        parts.append(f"{record.name}:{record.lineno}")

        line = " | ".join(parts)

        # Append exception info if present
        if record.exc_info and not record.exc_text:
            record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            line = f"{line}\n{record.exc_text}"

        return line


# MPI-safe handler
class MPISafeFileHandler(logging.FileHandler):
    """
    A custom file handler for logging that integrates with MPI environments.

    This logging handler manages files either shared across all ranks or independent for
    each rank in an MPI environment. It ensures both intra-process and inter-process
    safety while writing logs. For shared files, it employs POSIX file locks, while
    individual files rely on rank-specific naming to avoid conflicts. The class also
    handles directory creation for specified log file paths.

    :ivar _rank: The MPI rank of the process using this handler.
    :type _rank: int
    :ivar _per_rank: Indicates whether each rank should write to its own file.
    :type _per_rank: bool
    :ivar _flock_available: Indicates whether POSIX file lock `flock` is available.
    :type _flock_available: bool
    :ivar _thread_lock: A threading lock ensuring intra-process safety while writing logs.
    :type _thread_lock: threading.Lock
    """

    def __init__(self,
                 filename: str | Path,
                 mode: str = "a",
                 encoding: str | None = "utf-8",
                 delay: bool = False,
                 *,
                 per_rank: bool = True,
                 rank: int | None = None) -> None:
        """
        Initializes the class.

        This constructor sets up a file handler with support for rank-specific files
        and intra-process locking. It resolves the file path, creates necessary
        directories, and initializes the base class with parameters suited for file
        handling.

        :param filename: The target file's name or path.
        :param mode: The mode in which the file is opened, e.g., "a" for append.
        :param encoding: The encoding used for the file. Defaults to "utf-8".
            If None, the system default is applied.
        :param delay: If set to True, defers the file opening until the first
            call to emit().
        :param per_rank: When True, appends rank-specific identifiers to the
            filename for distributed setups.
        :param rank: The rank identifier for file naming. If None, it defaults
            to the MPI rank.
        """
        self._rank = rank if rank is not None else get_mpi_rank()
        self._per_rank = per_rank
        self._flock_available = hasattr(fcntl, "flock")
        self._thread_lock = threading.Lock()  # intra-process safety

        resolved = _resolve_path(filename, self._rank if per_rank else None)
        resolved.parent.mkdir(parents=True, exist_ok=True)

        super().__init__(str(resolved), mode=mode, encoding=encoding, delay=delay)

    # emit

    def emit(self, record: logging.LogRecord) -> None:
        """
        Emit a log record, handling thread safety and inter-process locking
        based on the `_per_rank` configuration.

        When `_per_rank` is enabled, each rank writes to its own file without
        requiring inter-process locking, relying only on thread safety. When
        disabled, inter-process locks ensure exclusive write access to the
        shared file.

        :param record: The log record to be processed.
        :type record: logging.LogRecord
        :return: None
        :rtype: None
        """
        if self._per_rank:
            # Each rank owns its file → no inter-process locking needed.
            with self._thread_lock:
                super().emit(record)
        else:
            # Shared file → exclusive inter-process lock per write.
            with self._thread_lock:
                self._emit_locked(record)

    def _emit_locked(self, record: logging.LogRecord) -> None:
        """
        Logs a record to the stream, ensuring thread or process-level locking
        mechanisms based on the current platform's capabilities.

        This method attempts to acquire an exclusive file lock when supported.
        Otherwise, it relies on a thread-only locking mechanism to prevent
        concurrent log writes into the same stream. After writing to the stream,
        this method flushes the stream to guarantee immediate log output delivery.

        :param record: The log record to emit.
        :type record: logging.LogRecord
        """
        stream: IO[str] = self.stream  # type: ignore[assignment]
        if stream is None:
            stream = self._open()

        if self._flock_available:
            try:
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
                logging.StreamHandler.emit(self, record)
                stream.flush()
            finally:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        else:
            # Windows fallback — intra-process thread lock only.
            logging.StreamHandler.emit(self, record)
            stream.flush()


def _resolve_path(filename: str | Path, rank: int | None) -> Path:
    """
    Resolve the file path based on the provided filename and rank. If a rank is specified,
    generates a new file name including the rank as part of the file stem; otherwise,
    returns the original file path. This function does not perform any validation for
    the existence of the file or directory.

    :param filename: File path as a string or Path object. This represents the
        base path of the file to resolve.
    :param rank: Integer rank to append to the file name. If None, the function
        returns the original path without modification.

    :return: A Path object representing the resolved file path, either with the rank
        appended to the file stem or unchanged if rank is None.
    :rtype: Path
    """
    p = Path(filename)
    if rank is None:
        return p
    return p.with_name(f"{p.stem}.rank_{rank}{p.suffix}")


# Factory


def make_file_handler(path: str | Path,
                      level: int = logging.DEBUG,
                      *,
                      per_rank: bool = True,
                      rank: int | None = None,
                      mode: str = "w",
                      encoding: str = "utf-8") -> MPISafeFileHandler:
    """
    Creates an instance of MPISafeFileHandler with custom configuration parameters.
    This function configures a file handler to log messages with specific attributes
    defined by the logging framework, which is beneficial for working in distributed
    systems (like MPI environments). It assigns a format to the log output, sets the
    desired logging level, and includes support for writing separate logs per rank.

    :param path: The path to the file where the log messages will be written.
    :type path: str | Path
    :param level: The minimum logging level that this handler will process (default is DEBUG).
    :type level: int, optional
    :param per_rank: A flag to determine if logs should maintain separate files per rank.
    :type per_rank: bool, optional
    :param rank: The rank identifier for specifying the log file. If not provided, it is fetched automatically.
    :type rank: int | None, optional
    :param mode: The mode in which the log file will be opened (default is 'w').
    :type mode: str, optional
    :param encoding: The file encoding for the log file (default is 'utf-8').
    :type encoding: str, optional
    :return: An instance of MPISafeFileHandler configured with the provided parameters.
    :rtype: MPISafeFileHandler
    """
    resolved_rank = rank if rank is not None else get_mpi_rank()
    handler = MPISafeFileHandler(
        path, mode=mode, encoding=encoding, per_rank=per_rank, rank=resolved_rank
    )
    handler.setLevel(level)
    handler.setFormatter(FileFormatter(rank=resolved_rank))
    return handler
