"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Dispatcher: routes the raw ``sys.argv[1:]`` list to the right subcommand.

Flow
----
1. No args / ``-h`` / ``--help``  →  top-level help
2. Unknown command                →  error + top-level help
3. ``<command>`` only or ``-h``  →  command help
4. Unknown subcommand             →  error + command help
5. ``<command> <subcommand> …``   →  build parser, validate, run
   (if ``-h`` / ``--help`` appears in the tail argparse prints its own help)
"""

from __future__ import annotations

import argparse
import sys

from pharmacon.command_line import formatter as fmt
from pharmacon.command_line.exceptions import ValidationError
from pharmacon.command_line.registry import CommandRegistry




__all__ = [
    "dispatch",
]


_HELP_FLAGS: frozenset[str] = frozenset({"-h", "--help", "help"})


def dispatch(argv: list[str], registry: CommandRegistry, exe_name: str) -> None:
    """
    Dispatches commands and subcommands by parsing arguments and running the appropriate
    command or subcommand logic. Handles unknown commands, unknown subcommands, validation,
    and execution of commands, along with appropriate error messaging and help outputs.

    :param argv: List of command-line arguments to parse and dispatch.
    :type argv: list[str]
    :param registry: A registry object containing the available commands and their specifications.
    :type registry: CommandRegistry
    :param exe_name: The name of the executable or application being run.
    :type exe_name: str
    :return: None
    """

    # top-level help
    if not argv or argv[0] in _HELP_FLAGS:
        fmt.print_top_level_help(exe_name, registry)
        return

    command_name = argv[0]

    # unknown command
    if command_name not in registry.commands:
        fmt.print_error(
            f"Unknown command [bold]{command_name!r}[/bold].",
            hint=f"Run '{exe_name} -h' to see available commands.",
        )
        sys.exit(1)

    command_spec = registry.commands[command_name]
    rest = argv[1:]

    # command-level help
    if not rest or rest[0] in _HELP_FLAGS:
        fmt.print_command_help(exe_name, command_spec)
        return

    subcommand_name = rest[0]

    # unknown subcommand
    if subcommand_name not in command_spec.subcommands:
        fmt.print_error(
            f"Unknown subcommand [bold]{subcommand_name!r}[/bold] "
            f"for command [bold]{command_name!r}[/bold].",
            hint=f"Run '{exe_name} {command_name} -h' to see available subcommands.",
        )
        sys.exit(1)

    subcommand_spec = command_spec.subcommands[subcommand_name]
    sub_args = rest[1:]

    # build parser and parse
    prog = f"{exe_name} {command_name}"
    root_parser = argparse.ArgumentParser(prog=prog, add_help=False)
    subparsers = root_parser.add_subparsers(dest="subcommand")

    subcommand_spec.build_parser_fn(subparsers, parents=None)

    # If the user asked for subcommand help, argparse will print it and exit.
    # The banner was already printed by main(), so nothing else is needed here.
    try:
        args, unknown = root_parser.parse_known_args([subcommand_name, *sub_args])
    except SystemExit as exc:
        sys.exit(exc.code)

    if unknown:
        fmt.print_error(
            f"Unrecognised argument(s): {' '.join(unknown)}",
            hint=f"Run '{exe_name} {command_name} {subcommand_name} -h' for usage.",
        )
        sys.exit(1)

    # validate
    try:
        subcommand_spec.validate_fn(args)
    except ValidationError as exc:
        fmt.print_error(str(exc))
        sys.exit(1)

    # execute
    try:
        subcommand_spec.run_fn(args)
    except KeyboardInterrupt:
        fmt.console.print("\n[bold yellow]Interrupted.[/bold yellow]")
        sys.exit(130)
    except Exception as exc:  # noqa: BLE001
        fmt.print_error(f"Runtime error: {exc}")
        sys.exit(1)
