"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Pharmacon logger package.

Public API
----------

Setup (call once, early in ``main()``)::

    from Pharmacon.logger import setup_logger

    setup_logger(
        terminal_level="INFO",       # or logging.INFO
        file_level="DEBUG",          # independent of terminal
        log_file="run.log",          # omit to disable file logging
        terminal=True,
        file=True,
        per_rank=True,               # MPI-safe: separate file per rank
    )

Acquire loggers (anywhere in the codebase)::

    from Pharmacon.logger import get_logger

    log = get_logger(__name__)

    log.trace("Very detailed trace")
    log.debug("Debug information")
    log.info("Normal progress message")
    log.warning("Something unexpected")
    log.error("Recoverable error")
    log.critical("Fatal condition")

Section headers::

    from Pharmacon.logger import header, subheader

    header("Protein Preparation")
    subheader("Fixing missing atoms")

Architecture
------------
* The ``Pharmacon`` root logger owns all handlers.
* ``get_logger(__name__)`` returns a child logger (e.g. ``Pharmacon.prepare.protein``);
  records propagate up to the root and are handled there.
* Terminal and file handlers are independent: each has its own level,
  formatter, and lock strategy.
* ``setup_logger()`` is idempotent: calling it twice with the same arguments
  has no effect; calling it with different arguments replaces the handlers.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Union

from pharmacon.logger.levels import TRACE, PharmaconLogger, register as _register_levels
from pharmacon.logger.headers import header, subheader

__all__ = [
    "setup_logger",
    "get_logger",
    "header",
    "subheader",
    "TRACE",
]

# Register PharmaconLogger globally so every getLogger() returns one.
_register_levels()

_LEVEL_MAP: dict[str, int] = {
    "trace":    TRACE,
    "debug":    logging.DEBUG,
    "info":     logging.INFO,
    "warning":  logging.WARNING,
    "warn":     logging.WARNING,
    "error":    logging.ERROR,
    "critical": logging.CRITICAL,
}

_ROOT_NAME = "pharmacon"


# Internal helpers


def _parse_level(level: Union[int, str]) -> int:
    """
    Parses and resolves the given logging level to its integer representation.

    The function accepts a logging level, which can either be an integer or a string
    representing a known logging level name. If the input is an integer, it is
    returned as is. If the input is a string, the method attempts to resolve it by
    matching it to a predefined mapping of level names. An exception is raised if
    the string does not match any known levels.

    :param level: The logging level to resolve. Can be an integer directly
        representing the level or a string corresponding to a level name.
    :type level: Union[int, str]
    :return: The integer representation of the resolved logging level.
    :rtype: int
    :raises ValueError: If the input string does not match any known logging
        level names.
    """
    if isinstance(level, int):
        return level
    resolved = _LEVEL_MAP.get(level.lower())
    if resolved is None:
        raise ValueError(
            f"Unknown log level {level!r}. "
            f"Valid names: {list(_LEVEL_MAP)}."
        )
    return resolved


def _clear_handlers(logger: logging.Logger) -> None:
    """
    Clear all handlers from the specified logger and close them.

    Removes every handler associated with the provided logger instance and ensures
    that each handler is properly closed.

    :param logger: The logger instance from which handlers will be removed and
                   closed.
    :type logger: logging.Logger
    :return: None
    """
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()


# Public API
def setup_logger(*,
                 terminal: bool = True,
                 terminal_level: Union[int, str] = logging.INFO,
                 file: bool = False,
                 file_level: Union[int, str] = logging.DEBUG,
                 log_file: Union[str, Path, None] = None,
                 log_queue=None,
                 queue_level: Union[int, str] = logging.DEBUG,
                 per_rank: bool = True,
                 mpi_rank: int | None = None,
                 replace: bool = True) -> PharmaconLogger:
    """
    Sets up a logger for the application, enabling logging to the terminal, file, or both.
    It allows for customization of log levels, handling of multi-rank scenarios in MPI setups,
    and replacement of existing handlers in the root logger.

    :param terminal: Flag to enable or disable logging to the terminal.
    :param terminal_level: Logging level for the terminal (e.g., INFO, DEBUG).
    :param file: Flag to enable or disable logging to a file.
    :param file_level: Logging level for the file (e.g., INFO, DEBUG).
    :param log_file: Path to the log file to write logs when `file` is True.
    :param per_rank: Determines if separate log files should be created per rank in MPI.
    :param mpi_rank: Rank of the current process in an MPI environment (if applicable).
    :param replace: If True, replaces all existing handlers in the root logger.
    :return: Configured PharmaconLogger instance.
    """
    t_level = _parse_level(terminal_level)
    f_level = _parse_level(file_level)
    q_level = _parse_level(queue_level)

    root: PharmaconLogger = logging.getLogger(_ROOT_NAME)  # type: ignore[assignment]

    if replace:
        _clear_handlers(root)

    # The root logger's own level must be ≤ the minimum of all handlers so
    # that records are not filtered before reaching any handler.
    effective_min = t_level if terminal else logging.CRITICAL
    if file:
        effective_min = min(effective_min, f_level)
    if log_queue is not None:
        effective_min = min(effective_min, q_level)
    root.setLevel(effective_min)

    # Terminal handler
    if terminal:
        from pharmacon.logger.terminal import make_terminal_handler

        t_handler = make_terminal_handler(level=t_level)
        root.addHandler(t_handler)

    # File handler
    if file:
        if log_file is None:
            raise ValueError("log_file must be provided when file=True.")
        from pharmacon.logger.file import make_file_handler

        f_handler = make_file_handler(
            log_file,
            level=f_level,
            per_rank=per_rank,
            rank=mpi_rank,
        )
        root.addHandler(f_handler)

    # Queue handler — ships LogRecords to another process (or thread) that
    # owns the real handlers. Used by multi-process workers so the main
    # process can write to a single shared log file in real time.
    if log_queue is not None:
        from logging.handlers import QueueHandler

        q_handler = QueueHandler(log_queue)
        q_handler.setLevel(q_level)
        root.addHandler(q_handler)

    # Prevent records from leaking into the root Python logger.
    root.propagate = False

    return root


def get_logger(name: str | None = None) -> PharmaconLogger:
    """
    Retrieve and configure a logger instance for the specified logger name. This function ensures that the
    logger hierarchy is properly maintained under the root logger of the library. If a name is not provided,
    defaults to the root logger of the "Pharmacon" framework. Ensures compatibility with modules outside
    the "Pharmacon" package namespace by appropriately qualifying the logger name.

    :param name: Optional; the name of the logger to retrieve. If None or corresponds to the root logger name,
        retrieves the root logger. If the provided name does not align with the package namespace,
        it is automatically qualified under the "Pharmacon" namespace.
    :return: An instance of `PharmaconLogger` corresponding to the specified logger name.
    """
    if name is None or name == _ROOT_NAME:
        return logging.getLogger(_ROOT_NAME)  # type: ignore[return-value]

    # Ensure the child is rooted under "Pharmacon" even when called from a
    # module outside the package namespace.
    if not name.startswith(_ROOT_NAME + ".") and name != _ROOT_NAME:
        qualified = f"{_ROOT_NAME}.{name}"
    else:
        qualified = name

    return logging.getLogger(qualified)  # type: ignore[return-value]
