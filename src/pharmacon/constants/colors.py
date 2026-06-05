"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Centralised color and style tokens for terminal formatting.

All CLI code imports from here so that the visual theme can be changed in one place.
Rich style strings are accepted by Console.print() markup, Text.stylize(), and
RichHelpFormatter.styles.
"""

from __future__ import annotations

from dataclasses import dataclass




__all__ = [
    "COLORS",
    "make_rich_argparse_styles",
]


@dataclass(frozen=True)
class _ColorPalette:
    # Palette anchors derived from rich-argparse defaults:
    #   #008080  cyan        → argparse.args   (flags / args)
    #   #ff8700  dark_orange → argparse.groups  (group headers / commands)
    #   #00af87  dark_cyan   → argparse.metavar (metavar / values)
    #   #808080  grey50      → argparse.prog    (dim / secondary chrome)

    # Branding
    title: str = "bold #008B8B"          # dark_cyan bold
    subtitle: str = "#ff8700"            # dark_orange
    version: str = "bold #00af87"        # dark_cyan bold

    # Paths / values
    path: str = "#00af87"                # dark_cyan
    value: str = "#00af87"              # dark_cyan

    # Status
    error: str = "bold red"              # keep standard red
    warning: str = "bold #ff8700"        # dark_orange bold
    info: str = "#008080"               # cyan
    success: str = "bold #00af87"        # dark_cyan bold

    # CLI structure
    command: str = "bold #ff8700"        # dark_orange bold  (≈ argparse.groups)
    subcommand: str = "#ff8700"          # dark_orange       (≈ argparse.groups)
    flag: str = "#008080"               # cyan              (≈ argparse.args)
    metavar: str = "#00af87"            # dark_cyan         (≈ argparse.metavar)
    positional: str = "#008080"          # cyan              (≈ argparse.args)

    # Text / chrome
    description: str = "default"         # terminal default  (≈ argparse.help)
    example_code: str = "#008080"        # cyan
    example_desc: str = "italic #808080" # grey50 italic     (≈ argparse.default)
    header: str = "bold default"         # bold default
    group_header: str = "bold #ff8700"   # dark_orange bold  (≈ argparse.groups)
    separator: str = "#808080"           # grey50            (≈ argparse.prog)
    dim: str = "#808080"                 # grey50            (≈ argparse.prog)

    # rich-argparse style keys — exact defaults
    argparse_prog: str = "#808080"       # grey50   (exact default)
    argparse_args: str = "#008080"       # cyan     (exact default)
    argparse_groups: str = "#ff8700"     # dark_orange (exact default)
    argparse_help: str = "default"       # default  (exact default)
    argparse_metavar: str = "#00af87"    # dark_cyan (exact default)
    argparse_syntax: str = "bold"        # bold     (exact default)
    argparse_default: str = "italic"     # italic   (exact default)
    argparse_text: str = "default"       # default  (exact default)


COLORS = _ColorPalette()


def make_rich_argparse_styles() -> dict[str, str]:
    """
    Creates and returns a dictionary of styles used by `argparse` with rich formatting
    styles. The returned dictionary maps specific `argparse` elements to their
    corresponding styles for improved visual appearance in terminal outputs.

    :return: A dictionary where the keys are `argparse` element identifiers (e.g.
        'argparse.prog', 'argparse.args') and the values are the corresponding
        rich style names used for formatting these elements.
    :rtype: dict[str, str]
    """
    return {
        "argparse.prog": COLORS.argparse_prog,
        "argparse.args": COLORS.argparse_args,
        "argparse.groups": COLORS.argparse_groups,
        "argparse.help": COLORS.argparse_help,
        "argparse.metavar": COLORS.argparse_metavar,
        "argparse.syntax": COLORS.argparse_syntax,
        "argparse.default": COLORS.argparse_default,
        "argparse.text": COLORS.argparse_text,
    }
