"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Core data-classes and structural interfaces for the CLI framework.

Hierarchy
---------
CommandSpec
    └─ SubcommandSpec   (one per .py module in a command directory)

ParseResult             carries the fully-parsed Namespace after dispatch.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Any, Callable


# Argument specification


__all__ = [
    "ArgumentSpec",
    "SubcommandSpec",
    "CommandSpec",
    "ParseResult",
]


@dataclass(frozen=True)
class ArgumentSpec:
    """
    A class that defines the specifications for command-line arguments.

    The ArgumentSpec class is designed to describe various properties and behaviors of
    command-line arguments, making it easier to define, organize, and handle them
    consistently in a program. This includes options like required arguments,
    default values, flags, choices, and more. It is particularly useful for
    generating argument parsers programmatically and capturing all related
    information for an argument in a structured way.

    :ivar name: The long name for the argument, without dashes (e.g. "input").
    :type name: str
    :ivar short_name: An optional single-character representation of the argument, without a dash
        (e.g. "i").
    :type short_name: str | None
    :ivar description: An optional descriptive text for the argument.
    :type description: str
    :ivar required: A flag indicating whether the argument is required. Defaults to False.
    :type required: bool
    :ivar default: The default value of the argument if not explicitly provided.
    :type default: Any
    :ivar is_flag: A flag indicating whether the argument functions as a boolean switch.
        When set to True, it implies a "store_true" behavior for the argument. Defaults
        to False.
    :type is_flag: bool
    :ivar multi_value: A flag indicating whether the argument accepts multiple
        values and can be repeated (e.g. action="append"). Defaults to False.
    :type multi_value: bool
    :ivar nargs: Defines the number of arguments expected (e.g., "+", "*", or an integer).
    :type nargs: str | int | None
    :ivar metavar: Provides a name for the additional values expected for display
        in parser help messages.
    :type metavar: str | None
    :ivar choices: A list of valid choices for the argument. If provided, the parser
        enforces these choices.
    :type choices: list[str] | None
    :ivar type: The expected data type of the argument. Defaults to str.
    :type type: type
    :ivar group: The organizational group to which the argument belongs, primarily used
        for grouping options in parser help. Defaults to "Options".
    :type group: str
    """

    # Names
    name: str                        # long name without dashes  (e.g. "input")
    short_name: str | None = None    # single-char without dash   (e.g. "i")

    # Behaviour
    description: str = ""
    required: bool = False
    default: Any = None
    is_flag: bool = False            # boolean switch → store_true
    multi_value: bool = False        # can be repeated → action="append"
    nargs: str | int | None = None   # for inline multi-value ("+" / "*" / N)
    metavar: str | None = None
    choices: list[str] | None = None
    type: type = str

    # Organisation
    group: str = "Options"

    # derived helpers

    @property
    def flags(self) -> list[str]:
        """
        Generates command-line flags based on the attributes of the instance.

        The method dynamically constructs flags for a command-line tool
        based on the `name` and `short_name` attributes of the object. The
        result will include a short version flag (if `short_name` is
        specified) and a long version flag derived from the `name`.

        :return: A list of strings representing the generated command-line
                 flags.
        :rtype: list[str]
        """
        parts: list[str] = []
        if self.short_name:
            parts.append(f"-{self.short_name}")
        parts.append(f"--{self.name.replace('_', '-')}")
        return parts

    def to_argparse_kwargs(self) -> dict[str, Any]:
        """
        Converts instance attributes to a dictionary of arguments suitable for defining
        a command-line argument in argparse.

        This method constructs a dictionary of keyword arguments based on the instance
        attributes. It supports different modes of argument configuration such as
        flags, multi-value arguments, and standard arguments. The resulting dictionary
        can directly be used with argparse's `add_argument` method.

        :raises AttributeError: If required attributes for argparse keyword generation
            are missing from the instance.

        :return: A dictionary containing keyword arguments for argparse's
            `add_argument` method.
        :rtype: dict[str, Any]
        """
        kwargs: dict[str, Any] = {"help": self.description}

        if self.is_flag:
            kwargs["action"] = "store_true"
            kwargs["default"] = False if self.default is None else self.default
        elif self.multi_value:
            kwargs["action"] = "append"
            kwargs["default"] = self.default
            kwargs["metavar"] = self.metavar or self.name.upper()
            kwargs["type"] = self.type
        else:
            kwargs["required"] = self.required
            kwargs["default"] = self.default
            if self.metavar:
                kwargs["metavar"] = self.metavar
            if self.choices:
                kwargs["choices"] = self.choices
            if self.nargs is not None:
                kwargs["nargs"] = self.nargs
            kwargs["type"] = self.type

        return kwargs


# Registry entries


@dataclass
class SubcommandSpec:
    """
    Representation of a subcommand specification for a command-line utility.

    This dataclass defines the structure for storing information related to a
    specific subcommand in a command-line application. It includes the subcommand's
    name, summary description, the module it belongs to, and callable functions for
    running the subcommand, building its argument parser, and validating its usage.

    :ivar name: Canonical name of the subcommand with hyphens (e.g. "water-clusters").
    :type name: str
    :ivar summary: A brief, one-line description displayed in help output.
    :type summary: str
    :ivar module_path: The dotted import path to the subcommand's related module
        (e.g. "pharmacon.command_line.analysis.rmsd").
    :type module_path: str
    """

    name: str                  # canonical name with hyphens  (e.g. "water-clusters")
    summary: str               # one-line description shown in command help
    module_path: str           # dotted import path  (e.g. "pharmacon.command_line.analysis.rmsd")

    # Callables loaded from the module at discovery time
    run_fn: Callable[[argparse.Namespace], None]
    build_parser_fn: Callable[
        [argparse._SubParsersAction, list[argparse.ArgumentParser] | None],
        argparse.ArgumentParser,
    ]
    validate_fn: Callable[[argparse.Namespace], None]


@dataclass
class CommandSpec:
    """
    Represents the specification for a command, including its name, summary, and
    any associated subcommands.

    This class is designed to encapsulate the structure of a command in a command-line
    interface application. It includes details about the command's name, a short
    summary explaining the command’s functionality, and a dictionary of subcommands
    associated with the command. Each subcommand is represented using instances of
    the `SubcommandSpec` class.

    :ivar name: The name of the command. Can include hyphens and serves as an identifier
        for the command in the application.
    :type name: str
    :ivar summary: A brief description or summary of the command's functionality. Retrieved
        from the `COMMAND_SUMMARY` attribute during initialization.
    :type summary: str
    :ivar subcommands: A dictionary mapping subcommand names (as strings) to their
        corresponding `SubcommandSpec` instances. Represents subcommands associated
        with this command.
    :type subcommands: dict[str, SubcommandSpec]
    """
    name: str                                   # directory name, hyphens allowed
    summary: str                                # from __init__.COMMAND_SUMMARY
    subcommands: dict[str, SubcommandSpec] = field(default_factory=dict)


# Parse result

@dataclass
class ParseResult:
    """
    Represents the result of parsing a command-line interface operation.

    This class is used to encapsulate the outcome of parsing a command-line
    input. It includes the main command, the subcommand, the parsed arguments,
    and the specification for the subcommand. It provides a standardized way
    to handle and store parsed CLI data.

    :ivar command: The main command issued by the user.
    :type command: str
    :ivar subcommand: The specific subcommand invoked within the main command.
    :type subcommand: str
    :ivar args: Parsed arguments associated with the command and subcommand.
    :type args: argparse.Namespace
    :ivar subcommand_spec: The specification detailing the behavior and structure
        of the invoked subcommand.
    :type subcommand_spec: SubcommandSpec
    """

    command: str
    subcommand: str
    args: argparse.Namespace
    subcommand_spec: SubcommandSpec
