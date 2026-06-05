"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

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
