"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Dynamic command and subcommand discovery.

The registry walks ``pharmacon/command_line/`` at import time, treats every
sub-directory as a *command*, and every ``.py`` file inside it (except
``__init__`` and ``common_args``) as a *subcommand*.

A module is accepted as a subcommand if and only if it exports:
    SUBCOMMAND_NAME : str
    SUMMARY         : str
    run(args)       : callable
    build_parser(subparsers, parents) : callable
    validate(args)  : callable

Adding a new subcommand requires nothing more than dropping a new ``.py``
file in the appropriate command directory.
"""

from __future__ import annotations

import importlib
import logging
from pathlib import Path

from pharmacon.command_line.base import CommandSpec, SubcommandSpec



__all__ = [
    "logger",
    "CommandRegistry",
    "get_registry",
]


logger = logging.getLogger(__name__)

_SKIP_NAMES: frozenset[str] = frozenset({"__init__", "common_args"})
_REQUIRED: tuple[str, ...] = ("SUBCOMMAND_NAME", "SUMMARY", "run", "build_parser", "validate")

# Registry class

class CommandRegistry:
    """Discovered snapshot of all available commands and subcommands."""

    def __init__(self) -> None:
        """
        Initializes an instance of the class.

        The constructor sets up the initial state of the object by initializing a
        dictionary to hold commands and invoking a method to populate it with data.
        This ensures that the object is configured properly and ready for use.

        Attributes:
            _commands (dict[str, CommandSpec]): Holds the mapping of command identifiers to
                their respective specifications.

        """
        self._commands: dict[str, CommandSpec] = {}
        self._discover()

    # public API

    @property
    def commands(self) -> dict[str, CommandSpec]:
        """
        Provides access to the commands dictionary. This property returns a mapping
        of command names to their specifications (`CommandSpec`). It is typically
        used to retrieve or manage the commands associated with this instance.

        :return: The dictionary of commands mapped to their specifications.
        :rtype: dict[str, CommandSpec]
        """
        return self._commands

    def get_command(self, name: str) -> CommandSpec | None:
        """
        Retrieves a command specification by its name.

        This method searches for a command specification corresponding
        to the provided name from the internal command registry. If no
        command with the specified name exists, the method returns `None`.

        :param name: The name of the command to retrieve.
        :type name: str
        :return: The command specification if found, otherwise `None`.
        :rtype: CommandSpec | None
        """
        return self._commands.get(name)

    def get_subcommand(self, command: str, subcommand: str) -> SubcommandSpec | None:
        """
        Gets a specific subcommand from a given command if it exists.

        This method attempts to retrieve a subcommand by its identifier
        if the parent command is available. If the parent command is not found,
        or the subcommand does not exist within it, the method will return None.

        :param command: The identifier of the parent command.
        :type command: str
        :param subcommand: The identifier of the subcommand to retrieve.
        :type subcommand: str
        :return: An instance of SubcommandSpec if the subcommand exists, otherwise None.
        :rtype: SubcommandSpec | None
        """
        cmd = self.get_command(command)
        return None if cmd is None else cmd.subcommands.get(subcommand)

    # discovery

    def _discover(self) -> None:
        """
        Discovers and registers commands from the `pharmacon.command_line` package.

        This method scans the directory of the `pharmacon.command_line` package to identify
        available commands. Each command directory is expected to represent a command, and
        its name is processed to generate a command name. Hidden directories or those starting
        with an underscore are skipped. Summary information and subcommands are extracted for
        each command and stored in the `_commands` attribute.

        The method outputs debug logs for tracking the discovered commands and their associated
        subcommands.

        :raises OSError: If there are issues accessing the command directory.
        """
        import pharmacon.command_line as _cl_pkg

        cl_dir = Path(_cl_pkg.__file__).parent  # type: ignore[arg-type]

        for item in sorted(cl_dir.iterdir()):
            if not item.is_dir() or item.name.startswith("_"):
                continue

            command_name = item.name.replace("_", "-")
            summary = self._command_summary(item)
            subcommands = self._discover_subcommands(item, command_name)

            self._commands[command_name] = CommandSpec(
                name=command_name,
                summary=summary,
                subcommands=subcommands,
            )
            logger.debug(
                "command '%s'  →  %d subcommand(s)", command_name, len(subcommands)
            )

    def _command_summary(self, cmd_dir: Path) -> str:
        """
        Generates a summary of a command based on the provided directory by attempting
        to dynamically load a corresponding module and retrieving its `COMMAND_SUMMARY`
        attribute. If the module does not exist or fails to load, an empty string is
        returned.

        :param cmd_dir: The directory representing the command whose summary should be
            retrieved.
        :type cmd_dir: Path
        :return: A string containing the command summary if available, otherwise an
            empty string.
        :rtype: str
        """
        pkg = f"pharmacon.command_line.{cmd_dir.name}"
        try:
            mod = importlib.import_module(pkg)
            return getattr(mod, "COMMAND_SUMMARY", "")
        except ImportError:
            return ""

    def _discover_subcommands(self, cmd_dir: Path, cmd_name: str) -> dict[str, SubcommandSpec]:
        """
        Discover and validate subcommands within a specified directory.

        This method scans the provided directory for Python files to identify
        and validate subcommands based on certain required attributes. It imports
        each module and confirms required attributes are present before registering
        them as discovered subcommands. The canonical name for each subcommand,
        its module path, and associated functions are stored in a dictionary.

        :param cmd_dir: The directory containing potential subcommand files.
        :type cmd_dir: Path
        :param cmd_name: The name of the parent command for which subcommands
            are being discovered.
        :type cmd_name: str
        :return: A dictionary mapping subcommand names to their respective
            specifications.
        :rtype: dict[str, SubcommandSpec]
        """
        pkg = f"pharmacon.command_line.{cmd_dir.name}"
        found: dict[str, SubcommandSpec] = {}

        for py_file in sorted(cmd_dir.glob("*.py")):
            stem = py_file.stem
            if stem in _SKIP_NAMES or stem.startswith("_"):
                continue

            full = f"{pkg}.{stem}"
            try:
                mod = importlib.import_module(full)
            except ImportError as exc:
                logger.warning("Cannot import '%s': %s", full, exc)
                continue

            missing = [a for a in _REQUIRED if not hasattr(mod, a)]
            if missing:
                logger.debug("Skipping '%s': missing %s", full, missing)
                continue

            sc_name: str = mod.SUBCOMMAND_NAME  # canonical name (may have hyphens)
            found[sc_name] = SubcommandSpec(
                name=sc_name,
                summary=mod.SUMMARY,
                module_path=full,
                run_fn=mod.run,
                build_parser_fn=mod.build_parser,
                validate_fn=mod.validate,
            )
            logger.debug("  subcommand '%s %s'", cmd_name, sc_name)

        return found


# Module-level singleton

_registry: CommandRegistry | None = None


def get_registry() -> CommandRegistry:
    """Return the shared registry, creating it on first call."""
    global _registry
    if _registry is None:
        _registry = CommandRegistry()
    return _registry
