"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Custom log levels and the PharmaconLogger class.

Extends the standard hierarchy with TRACE (5), below DEBUG (10).

    TRACE    =  5
    DEBUG    = 10
    INFO     = 20
    WARNING  = 30
    ERROR    = 40
    CRITICAL = 50

Register before any logger is created by calling ``register()``.
"""

from __future__ import annotations

import logging

# ── Level constant ─────────────────────────────────────────────────────────────



__all__ = [
    "TRACE",
    "PharmaconLogger",
    "register",
]


TRACE: int = 5
logging.addLevelName(TRACE, "TRACE")


# ── Logger subclass ────────────────────────────────────────────────────────────


class PharmaconLogger(logging.Logger):
    """Logger subclass that adds a :meth:`trace` convenience method."""

    def trace(self, msg: object, *args: object, **kwargs: object) -> None:
        """Log at TRACE level (verbosity below DEBUG)."""
        if self.isEnabledFor(TRACE):
            self._log(TRACE, msg, args, **kwargs)


# ── Registration ───────────────────────────────────────────────────────────────


def register() -> None:
    """Register PharmaconLogger as the default logger class.

    Must be called once, before any ``logging.getLogger()`` call, to ensure
    all subsequently created loggers are ``PharmaconLogger`` instances.
    """
    logging.setLoggerClass(PharmaconLogger)


# Backward-compatible alias used by the logger __init__ module.
