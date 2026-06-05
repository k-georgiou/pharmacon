"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Pharmacon CLI entrypoint.

Invoked as:
    pharmacon [command] [subcommand] [flags...]

Registered in pyproject.toml as:
    [project.scripts]
    pharmacon = "pharmacon.main:main"
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pyjokes

from pharmacon.command_line import formatter as fmt
from pharmacon.command_line.dispatcher import dispatch
from pharmacon.command_line.registry import get_registry
from pharmacon.constants import __version__, __author__, __github__, __manuscript__




__all__ = [
    "main",
]


def _exe_name() -> str:
    """
    Determines the executable name from command-line arguments or defaults to
    a predefined value.

    If a script is being executed in development mode (e.g., using "__main__.py",
    "main.py", or "-c"), the function will return a default executable name
    "pharmacon". Otherwise, it returns the name of the current script.

    :return: The normalized or default executable name.
    :rtype: str
    """
    name = Path(sys.argv[0]).name if sys.argv else "pharmacon"
    # Normalise common dev-mode launcher names
    if name in ("__main__.py", "main.py", "-c"):
        return "pharmacon"
    return name


def main() -> None:
    """
    This is the main entry point for the script execution. It performs the
    following tasks:

    1. Retrieves the executable name and current working directory.
    2. Processes command-line arguments provided to the script.
    3. Displays a formatted banner that includes information like
       executable name, version, author, manuscript reference, repository URL,
       working directory, and command-line arguments used.
    4. Obtains a command registry.
    5. Dispatches the provided command-line arguments for execution within
       the command registry.
    6. Handles interruption signals (like keyboard interrupts) gracefully
       and terminates with appropriate messaging.
    7. Catches all other unexpected exceptions, optionally re-raises them
       in debug mode, and logs an error when not debugging.
    8. Displays a random programming joke in the console upon successful
       script execution for user delight.

    :param exe: The name of the executable.
    :type exe: str
    :param cwd: The current working directory from which the script is
        being executed.
    :type cwd: str
    :param argv: List of command-line arguments provided to the script.
    :type argv: list
    :return: None
    """
    exe = _exe_name()
    cwd = os.getcwd()
    argv = sys.argv[1:]



    fmt.print_banner(
        exe_name=exe,
        version=__version__,
        author=__author__,
        manuscript=__manuscript__,
        github=__github__,
        cwd=cwd,
        argv=argv)


    registry = get_registry()

    try:
        dispatch(argv, registry, exe)
    except KeyboardInterrupt:
        fmt.console.print("\n[bold yellow]Interrupted.[/bold yellow]")
        sys.exit(130)
    except Exception as exc:  # noqa: BLE001
        fmt.print_error(f"Unexpected error: {exc}")
        if os.getenv("PHARMACON_DEBUG") or "--debug" in sys.argv:
            raise
        sys.exit(1)

    # Success joke
    joke = pyjokes.get_joke()
    fmt.console.print(f"\n[bold italic #98FB98]{joke}[/]")


if __name__ == "__main__":
    main()
