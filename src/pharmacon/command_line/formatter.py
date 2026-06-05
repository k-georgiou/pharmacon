"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Rich-powered help formatter for the pharmacon CLI.

Responsibilities
----------------
* Banner  — pyfiglet title + version / author / manuscript / github / cwd / command
* Top-level help   — table of available commands
* Command help     — table of available subcommands
* Subcommand help  — delegated to argparse + RichHelpFormatter
* Error display    — colourised error / hint messages
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pyfiglet
from rich.console import Console
from rich.rule import Rule
from rich.table import Table
from rich.text import Text
from rich_argparse import RichHelpFormatter

from pharmacon.constants.colors import COLORS, make_rich_argparse_styles



__all__ = [
    "console",
    "print_banner",
    "print_top_level_help",
    "print_command_help",
    "make_formatter_class",
    "print_error",
]


if TYPE_CHECKING:
    from pharmacon.command_line.base import CommandSpec
    from pharmacon.command_line.registry import CommandRegistry

console = Console(highlight=False)

_FIGLET_FONT = "ansi_shadow"

# Banner

def print_banner(exe_name: str,
                 version: str,
                 author: str,
                 manuscript: str,
                 github: str,
                 cwd: str,
                 argv: list[str]) -> None:
    """
    Prints a banner with application details.

    This function displays an ASCII art title and various metadata rows,
    providing information about the application's version, author, manuscript,
    GitHub repository, current working directory, and the command used to
    invoke the program.

    :param exe_name: Name of the executable or application.
    :type exe_name: str
    :param version: Current version of the application.
    :type version: str
    :param author: Name of the author of the application.
    :type author: str
    :param manuscript: Path or reference to the manuscript related
        to the application.
    :type manuscript: str
    :param github: URL of the GitHub repository associated with the
        application.
    :type github: str
    :param cwd: Current working directory when executing the
        application.
    :type cwd: str
    :param argv: List of command-line arguments provided to the
        application.
    :type argv: list[str]
    :return: This function does not return any value.
    :rtype: None
    """
    # ASCII art title
    ascii_title = pyfiglet.figlet_format("PHARMACON", font=_FIGLET_FONT)
    console.print(Text(ascii_title.rstrip(), style=COLORS.title))

    # Metadata rows
    cli_cmd = exe_name + (" " + " ".join(argv) if argv else "")

    _meta_row("version",    f"v{version}",  COLORS.version)
    _meta_row("author",     author,         COLORS.description)
    _meta_row("manuscript", manuscript,     COLORS.info)
    _meta_row("github",     github,         COLORS.info)
    _meta_row("cwd",        cwd,            COLORS.path)
    _meta_row("command",    cli_cmd,        COLORS.command)
    console.print()


def _meta_row(label: str, value: str, value_style: str) -> None:
    """
    Generates and prints a formatted line containing a label and its corresponding
    value, styled appropriately for console output.

    :param label: The label to display, aligned to the left within a fixed width.
    :param value: The value associated with the label.
    :param value_style: The style to apply to the value for console output.
    :return: None
    """
    line = Text()
    line.append(f"  {label:<12}", style=COLORS.dim)
    line.append(": ", style=COLORS.dim)
    line.append(value, style=value_style)
    console.print(line)

# Top-level help:  exe -h

def print_top_level_help(exe_name: str, registry: CommandRegistry) -> None:
    """
    Generates and prints the top-level help information to the console for a
    command-line interface. This includes usage instructions, available
    commands, and guidance on listing subcommands using the provided
    command registry.

    :param exe_name: The name of the executable or primary command-line program.
    :type exe_name: str
    :param registry: The registry containing all registered commands and
        their metadata.
    :type registry: CommandRegistry
    :return: None
    """
    console.print()
    console.print(
        f"  [bold]Usage:[/bold]  "
        f"[{COLORS.command}]{exe_name}[/]  "
        f"[{COLORS.subcommand}]<command>[/]  "
        f"[{COLORS.subcommand}]<subcommand>[/]  "
        f"[{COLORS.metavar}]\\[flags...][/]"
    )
    console.print()
    console.print(Rule("Available Commands", style=COLORS.separator))
    console.print()

    tbl = Table(box=None, padding=(0, 4), show_header=False)
    tbl.add_column("cmd", style=COLORS.command, no_wrap=True, min_width=14)
    tbl.add_column("desc", style=COLORS.description)

    for cmd in registry.commands.values():
        tbl.add_row(cmd.name, cmd.summary)

    console.print(tbl)
    console.print()
    console.print(
        f"  Run  [{COLORS.command}]{exe_name}[/]  "
        f"[{COLORS.subcommand}]<command>[/]  "
        f"[{COLORS.metavar}]-h[/]  "
        f"[{COLORS.dim}]to list subcommands.[/]"
    )
    console.print()


# Command help:  exe prepare -h

def print_command_help(exe_name: str, command_spec: CommandSpec) -> None:
    """
    Provides a detailed description of how to print a command's help message,
    including its usage and associated subcommands, using formatted console output.

    :param exe_name: The name of the executable or command-line tool.
    :type exe_name: str
    :param command_spec: The specification of the command including its name,
        summary, and subcommands.
    :type command_spec: CommandSpec
    :return: This function does not return any value; it outputs the help
        content directly to the console.
    :rtype: None
    """
    console.print()
    console.print(
        f"  [bold]Usage:[/bold]  "
        f"[{COLORS.command}]{exe_name}[/]  "
        f"[{COLORS.command}]{command_spec.name}[/]  "
        f"[{COLORS.subcommand}]<subcommand>[/]  "
        f"[{COLORS.metavar}]\\[flags...][/]"
    )
    if command_spec.summary:
        console.print(f"  [{COLORS.dim}]{command_spec.summary}[/]")
    console.print()
    console.print(Rule(f"Subcommands  —  {command_spec.name}", style=COLORS.separator))
    console.print()

    tbl = Table(box=None, padding=(0, 4), show_header=False)
    tbl.add_column("sub", style=COLORS.subcommand, no_wrap=True, min_width=18)
    tbl.add_column("desc", style=COLORS.description)

    for sc in command_spec.subcommands.values():
        tbl.add_row(sc.name, sc.summary)

    console.print(tbl)
    console.print()
    console.print(
        f"  Run  [{COLORS.command}]{exe_name}[/]  "
        f"[{COLORS.command}]{command_spec.name}[/]  "
        f"[{COLORS.subcommand}]<subcommand>[/]  "
        f"[{COLORS.metavar}]-h[/]  "
        f"[{COLORS.dim}]for detailed flag help.[/]"
    )
    console.print()


# Subcommand formatter class (rich-argparse)

def make_formatter_class() -> type:
    """
    Creates a custom argparse help formatter class with rich styling and highlighting
    capabilities for use in command-line argument parsing interfaces. The created
    formatter class provides enhanced presentation for CLI help messages with rich
    text styles and syntax highlighting based on predefined patterns.

    :return: The custom help formatter class.
    :rtype: type
    """

    _styles = make_rich_argparse_styles()

    class PharmaconHelpFormatter(RichHelpFormatter):
        """Custom argparse help formatter class with rich styling and highlighting."""
        styles = _styles  # type: ignore[assignment]
        highlights = [
            r"(?P<args>-{1,2}[\w-]+)",        # CLI flags
            r"(?P<metavar>\b[A-Z][A-Z0-9_]+\b)",  # METAVAR tokens
        ]

    return PharmaconHelpFormatter


# Error display

def print_error(message: str, hint: str | None = None) -> None:
    """
    Prints an error message to the console with an optional hint.

    This function formats and prints a given error message in a styled manner using a
    console output library. If a hint is provided, it is printed below the error message
    in a distinct style. A blank line is added after the error message for better readability.

    :param message: The main error message to be displayed.
    :type message: str
    :param hint: Optional hint providing additional context or suggestion. If None, it will be omitted.
    :type hint: str | None
    :return: None
    """
    console.print(f"[{COLORS.error}]Error:[/]  {message}")
    if hint:
        console.print(f"  [{COLORS.dim}]Hint:  {hint}[/]")
    console.print()
