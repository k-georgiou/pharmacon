"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Shared parent parsers for common argument groups.

Usage in a subcommand module::

    from pharmacon.command_line.common_args import VERBOSE_PARENT, LOGGING_PARENT

    def build_parser(subparsers, parents=None):
        parents = [VERBOSE_PARENT, LOGGING_PARENT, *(parents or [])]
        parser = subparsers.add_parser(..., parents=parents, ...)
        ...

Each factory returns a *new* ArgumentParser so argparse can safely merge them
(parents are consumed during add_parser).  The module-level singletons are
convenience aliases for the most common patterns — callers that need a fresh
instance should call the factory directly.
"""

from __future__ import annotations

import argparse


# Factories




__all__ = [
    "make_verbose_parser",
    "make_logging_parser",
    "make_trajectory_input_parser",
    "make_time_range_parser",
    "PARENTS",
]


def make_verbose_parser() -> argparse.ArgumentParser:
    """
    Creates and returns an ArgumentParser with a predefined argument for enabling
    verbose or debug output.

    The argument '-v' or '--verbose' can be used to toggle verbose/debug mode.

    :return: ArgumentParser object with the verbose argument preset.
    """
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose / debug output.",
    )
    return p


def make_logging_parser() -> argparse.ArgumentParser:
    """
    Creates an argument parser for configuring logging options.

    This function sets up a CLI argument parser with options for managing
    logging behavior, such as verbosity level and log output file.

    :return: A configured instance of `argparse.ArgumentParser`
             for managing logging parameters.
    :rtype: argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(add_help=False)
    grp = p.add_argument_group("Logging Options")
    grp.add_argument(
        "--log-level",
        metavar="LEVEL",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    grp.add_argument(
        "--log-file",
        metavar="FILE",
        default=None,
        help="Write log output to FILE instead of stderr.",
    )
    return p


def make_trajectory_input_parser() -> argparse.ArgumentParser:
    """
    Creates and returns an argument parser for processing trajectory input options. The parser
    contains argument definitions for specifying topology/structure file and trajectory files.
    It is designed to handle repeated input for multiple trajectory files.

    :return: An :class:`argparse.ArgumentParser` instance configured with input options.
    :rtype: argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(add_help=False)
    grp = p.add_argument_group("Input Options")
    grp.add_argument(
        "-s", "--structure",
        required=True,
        metavar="FILE",
        help="Topology / structure file (.tpr, .gro, .pdb).",
    )
    grp.add_argument(
        "-f", "--trajectory",
        required=True,
        action="append",
        metavar="FILE",
        dest="trajectory",
        help="Trajectory file(s); may be repeated: -f a.xtc -f b.xtc.",
    )
    return p


def make_time_range_parser() -> argparse.ArgumentParser:
    """
    Creates an ArgumentParser instance pre-configured with options for specifying
    a time range.

    The parser provides options to define the start and end time in picoseconds,
    as well as a step value for selecting every Nth frame. These options are
    grouped under a "Time Range Options" argument group.

    :return: An ArgumentParser instance with time range arguments configured.
    :rtype: argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(add_help=False)
    grp = p.add_argument_group("Time Range Options")
    grp.add_argument(
        "-b", "--begin",
        metavar="PS",
        type=float,
        default=0.0,
        help="Start time in ps (default: 0.0).",
    )
    grp.add_argument(
        "-e", "--end",
        metavar="PS",
        type=float,
        default=-1.0,
        help="End time in ps; -1 means last frame (default: -1.0).",
    )
    grp.add_argument(
        "--step",
        metavar="N",
        type=int,
        default=1,
        help="Process every Nth frame (default: 1).",
    )
    return p


# Module-level singletons (convenience)
# These are re-created each time so argparse does not exhaust them.

def _fresh(factory):  # type: ignore[no-untyped-def]
    """Descriptor that calls factory() on every attribute access."""
    class _Prop:
        def __get__(self, obj, cls):  # noqa: ANN001
            return factory()
    return _Prop()


class _Parents:
    """Namespace of always-fresh parent parsers."""
    verbose = _fresh(make_verbose_parser)
    logging = _fresh(make_logging_parser)
    trajectory_input = _fresh(make_trajectory_input_parser)
    time_range = _fresh(make_time_range_parser)


PARENTS = _Parents()
