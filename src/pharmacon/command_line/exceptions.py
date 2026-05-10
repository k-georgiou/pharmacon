"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

CLI-specific exception hierarchy.
"""

from __future__ import annotations




__all__ = [
    "PharmaconError",
    "ValidationError",
    "DiscoveryError",
    "DispatchError",
]


class PharmaconError(Exception):
    """Base exception for all Pharmacon CLI errors."""


class ValidationError(PharmaconError):
    """Raised when argument validation fails (missing file, bad value, …)."""


class DiscoveryError(PharmaconError):
    """Raised when dynamic command/subcommand discovery fails."""


class DispatchError(PharmaconError):
    """Raised when the dispatcher cannot route the request."""
