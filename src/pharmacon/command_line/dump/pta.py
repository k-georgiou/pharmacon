"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

dump pta — Dumps to STDOUT information about a Pharmacon Trajectory Analysis file.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

from pharmacon.command_line.exceptions import ValidationError
from pharmacon.command_line.formatter import make_formatter_class
from pharmacon.utils.validation import (
    validate_bool_flag,
    validate_existing_input_file,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "pta"
SUMMARY = "Dumps to STDOUT information about a Pharmacon Trajectory Analysis file."

_EPILOG = """\
Examples:

  pharmacon dump pta -i input_file.pta

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures the parser for handling a specific subcommand.

    This function adds a new subparser to the provided `subparsers` object for managing
    Pharmacon Trajectory Analysis (PTA) files. It includes options for specifying input
    files, toggling features such as metadata output, and managing verbosity.

    :param subparsers: A subparsers action object to which the new parser will be added.
    :param parents: A list of argparse ArgumentParser objects to inherit arguments from,
        or `None` to avoid inheritance.
    :return: Configured argparse.ArgumentParser instance for the subcommand.
    """
    parser: argparse.ArgumentParser = subparsers.add_parser(
        SUBCOMMAND_NAME,
        help=SUMMARY,
        description=SUMMARY,
        epilog=_EPILOG,
        formatter_class=make_formatter_class(),
        parents=parents or [],
        add_help=True,
    )

    # Input Options
    inp = parser.add_argument_group("Input Options")
    inp.add_argument(
        "-i", "--input",
        required=True,
        metavar="FILE",
        help="Input Pharmacon Trajectory Analysis (pta) file. [REQUIRED]",
    )

    # Toggle Options
    tog = parser.add_argument_group("Toggle Options")
    tog.add_argument(
        "--disable-file-metadata",
        action="store_true",
        default=False,
        help="Disable dumping of file metadata to the STDOUT. Default='False'. Type='bool' [OPTIONAL]",
        required=False,
    )
    tog.add_argument(
        "--disable-group-meta",
        action="store_true",
        help="Disable dumping group metadata of PTA file. Default='False'. Type=bool [OPTIONAL]")

    tog.add_argument(
        "--disable-dataset-meta",
        action="store_true",
        help="Disable dump dataset metadata of PTA file. Default='False'. Type=bool [OPTIONAL]")

    tog.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose dump to STDOUT. Default='False'. Type=bool [OPTIONAL]")

    return parser


def validate(args: argparse.Namespace) -> None:
    """
    Validates the arguments provided in an `argparse.Namespace` object for correctness
    and ensures that the input PTA file meets specific metadata integrity requirements.

    This function performs a series of checks and validations on the provided input,
    including ensuring the existence and type correctness of input files and flags.
    It also verifies metadata properties of the PTA file to ensure it is complete,
    unchanged, and valid according to pre-defined integrity constraints.

    :param args: Arguments parsed and encapsulated in an `argparse.Namespace` object.
                 The function modifies this object in-place with processed and
                 validated values.
    :type args: argparse.Namespace

    :raises ValidationError: If any validation criteria are not met, including:
        - File existence checks
        - PTA metadata completeness and correctness
        - Consistency in artifact token, status, blueprint, signature, fingerprint, etc.
    """
    data: dict[str, object] = vars(args).copy()

    # Input file
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        [".pta"],
    )
    data["input"] = input_file

    # Boolean flags
    disable_file_metadata: bool = validate_bool_flag(
        data.get("disable_file_metadata", False),
        "disable_file_metadata",
    )
    data["disable_file_metadata"] = disable_file_metadata

    disable_group_meta: bool = validate_bool_flag(
        data.get("disable_group_meta", False),
        "disable_group_meta",
    )
    data["disable_group_meta"] = disable_group_meta

    disable_dataset_meta: bool = validate_bool_flag(
        data.get("disable_dataset_meta", False),
        "disable_dataset_meta",
    )
    data["disable_dataset_meta"] = disable_dataset_meta

    verbose: bool = validate_bool_flag(
        data.get("verbose", False),
        "verbose",
    )
    data["verbose"] = verbose

    # PTA file metadata integrity checks
    import h5py
    from pharmacon.utils.fingerprint import create_pharmacon_signature
    from pharmacon.utils.identifiers import validate_mda_artifact_token

    try:
        f = h5py.File(input_file, "r")
    except Exception as exc:
        raise ValidationError(f"Cannot open PTA file: {exc}") from exc

    try:
        attrs = dict(f.attrs)

        # artifact_status must be SUCCESS
        artifact_status = str(attrs.get("artifact_status", "")).strip().upper()
        if not artifact_status:
            raise ValidationError(
                "PTA file is missing 'artifact_status' metadata — "
                "the file may be incomplete or corrupted."
            )
        if artifact_status != "SUCCESS":
            raise ValidationError(
                f"PTA file artifact_status is '{artifact_status}' (expected 'SUCCESS') — "
                "the file may be corrupted or the analysis did not complete."
            )

        # artifact_token must be present and valid
        artifact_token = str(attrs.get("artifact_token", "")).strip()
        if not artifact_token:
            raise ValidationError(
                "PTA file is missing 'artifact_token' metadata — "
                "the file may be corrupted."
            )

        blueprint = str(attrs.get("blueprint", "")).strip()
        if not blueprint:
            raise ValidationError(
                "PTA file is missing 'blueprint' metadata — "
                "the file may be corrupted."
            )

        token_valid = validate_mda_artifact_token(
            artifact_token=artifact_token,
            blueprint=blueprint,
            secret="trajectory_analysis",
            namespace="pharmacon",
        )
        if not token_valid:
            raise ValidationError(
                "PTA file artifact_token does not match the blueprint — "
                "the file may have been tampered with or corrupted."
            )

        # signature / fingerprint must be present and consistent
        signature = str(attrs.get("signature", "")).strip()
        fingerprint = str(attrs.get("fingerprint", "")).strip()
        command = str(attrs.get("command", "")).strip()
        subcommand = str(attrs.get("subcommand", "")).strip()

        if not signature or not fingerprint:
            raise ValidationError(
                "PTA file is missing 'signature' and/or 'fingerprint' metadata — "
                "the file may be corrupted."
            )

        if command and subcommand:
            expected_sig = create_pharmacon_signature(
                format_name="pta",
                command=command,
                subcommand=subcommand,
            )
            if expected_sig.signature != signature:
                raise ValidationError(
                    f"PTA file signature mismatch.\n"
                    f"  Expected : {expected_sig.signature}\n"
                    f"  Found    : {signature}\n"
                    "The file may have been tampered with or corrupted."
                )
            if expected_sig.fingerprint != fingerprint:
                raise ValidationError(
                    f"PTA file fingerprint mismatch.\n"
                    f"  Expected : {expected_sig.fingerprint}\n"
                    f"  Found    : {fingerprint}\n"
                    "The file may have been tampered with or corrupted."
                )

        # pharmacon_version must be present
        version = str(attrs.get("pharmacon_version", "")).strip()
        if not version:
            raise ValidationError(
                "PTA file is missing 'pharmacon_version' metadata — "
                "the file may be corrupted."
            )

    finally:
        f.close()

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes a command-line utility function to parse and display metadata and structural
    information about a given Pharmacon-compatible HDF5 file. The function inspects file
    validity, versioning, and provides configurable options for displaying metadata and
    tree structure.

    :param args: An argparse.Namespace object containing the following attributes:
                 - input: Path to the input HDF5 file.
                 - verbose: Whether verbose output is enabled.
                 - disable_file_metadata: Disable display of file-level metadata.
                 - disable_group_meta: Disable display of group metadata.
                 - disable_dataset_meta: Disable display of dataset metadata.
    :return: None
    """
    from rich.console import Console
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.fileio.base import PharmaconHDF5File
    from pharmacon.logger import header, subheader

    console = Console()

    # Header
    header("Dump PTA")

    # Run configuration
    cfg = Table(
        show_header=False,
        box=SIMPLE,
        pad_edge=False,
        padding=(0, 1),
    )
    cfg.add_column("Key", style="bold cyan", no_wrap=True, min_width=20)
    cfg.add_column("Value", style="white")

    cfg.add_row("Input", str(args.input))
    cfg.add_row("Verbose", str(args.verbose))
    cfg.add_row("File metadata", "disabled" if args.disable_file_metadata else "enabled")
    cfg.add_row("Group metadata", "disabled" if args.disable_group_meta else "enabled")
    cfg.add_row("Dataset metadata", "disabled" if args.disable_dataset_meta else "enabled")
    console.print(cfg)

    # Open the PTA file in read-only mode
    with PharmaconHDF5File(args.input, mode="r") as pta:

        # File validity status
        subheader("File Validity")
        pta.print_file_validity_status()

        # Version check
        pta.validate_version(strict=False)

        # ile metadata
        if not args.disable_file_metadata:
            subheader("File Metadata")
            pta.print_file_metadata()

        # Tree structure
        subheader("File Structure")
        pta.print_tree(
            group_meta=not args.disable_group_meta,
            dataset_meta=not args.disable_dataset_meta,
            compact=not args.verbose,
            max_items=5,
            file_attrs=args.disable_file_metadata,
        )

        # Per-input metadata for merged files (verbose only)
        if args.verbose:
            _print_input_file_metadata_tables(pta, console)

    # Done
    header("Done")


def _print_input_file_metadata_tables(pta, console) -> None:
    """
    Prints the metadata tables of input files from the provided file object and
    displays them via the console.

    This function scans the given file-like object for metadata groups matching
    a specific pattern, retrieves relevant key-value pairs, and formats them
    as rich tables for display. If no matching metadata groups are found,
    the function does nothing.

    :param pta: A file object (e.g., HDF5 file handler) containing groups and attributes.
    :type pta: Any
    :param console: A console object from the `rich` library used for rendering tables.
    :type console: Console
    :return: This function does not return anything; outputs are printed to the console.
    :rtype: None
    """
    import re
    import h5py
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.logger import subheader

    pattern = re.compile(r"^file(\d+)_metadata$")

    matches: list[tuple[int, str]] = []
    for key in pta.file:
        m = pattern.match(key)
        if m is None:
            continue
        obj = pta.file[key]
        if not isinstance(obj, h5py.Group):
            continue
        matches.append((int(m.group(1)), key))

    if not matches:
        return

    matches.sort(key=lambda kv: kv[0])

    subheader("Input Files Metadata")

    for _, group_name in matches:
        grp = pta.file[group_name]

        table = Table(
            title=group_name,
            show_header=True,
            header_style="bold cyan",
            box=SIMPLE,
            pad_edge=False,
            padding=(0, 1),
        )
        table.add_column("Key", style="bold cyan", no_wrap=True, min_width=20)
        table.add_column("Value", style="white", overflow="fold")

        for key in sorted(grp.attrs):
            table.add_row(str(key), str(grp.attrs[key]))

        console.print(table)
