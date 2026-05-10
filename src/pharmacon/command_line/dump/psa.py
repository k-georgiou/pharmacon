"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

dump psa — Dumps to STDOUT information about a Pharmacon Structure Analysis file.
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


SUBCOMMAND_NAME = "psa"
SUMMARY = "Dumps to STDOUT information about a Pharmacon Structure Analysis file."

_EPILOG = """\
Examples:

  pharmacon dump psa -i input_file.psa

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and returns an argparse.ArgumentParser configured for handling command-line
    arguments specific to the PSA file analysis. The parser includes different groups
    of options for input arguments and toggle switches that control how data is
    processed and output.

    :param subparsers: Subparsers action object where this parser will be added.
    :type subparsers: argparse._SubParsersAction
    :param parents: List of argparse.ArgumentParser objects that contain options
        to be inherited by this parser. Defaults to None.
    :type parents: List[argparse.ArgumentParser] | None
    :return: Configured argparse.ArgumentParser instance for the PSA command-line
        options.
    :rtype: argparse.ArgumentParser
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

    # ── Input Options ──────────────────────────────────────────────────────────
    inp = parser.add_argument_group("Input Options")
    inp.add_argument(
        "-i", "--input",
        required=True,
        metavar="FILE",
        help="Input Pharmacon Structure Analysis (psa) file. [REQUIRED]",
    )

    # ── Toggle Options ──────────────────────────────────────────────────────────
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
        help="Disable dumping group metadata of PSA file. Default='False'. Type=bool [OPTIONAL]")

    tog.add_argument(
        "--disable-dataset-meta",
        action="store_true",
        help="Disable dump dataset metadata of PSA file. Default='False'. Type=bool [OPTIONAL]")

    tog.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose dump to STDOUT. Default='False'. Type=bool [OPTIONAL]")

    return parser


def validate(args: argparse.Namespace) -> None:
    """
    Validates the command-line arguments passed to the program. The function ensures
    that the input file exists, the boolean flags are correctly parsed, and performs
    metadata integrity checks for PSA files. The function modifies the provided
    ``argparse.Namespace`` object in-place to include validated attributes.

    :param args: Namespace object containing the command-line arguments to validate.
    :type args: argparse.Namespace
    :raises ValidationError: If the input file is invalid, necessary metadata is
        missing or corrupted, or integrity checks fail.
    :return: None
    """
    data: dict[str, object] = vars(args).copy()

    # Input file
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        [".psa"],
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

    # PSA file metadata integrity checks
    import h5py
    from pharmacon.utils.fingerprint import create_pharmacon_signature
    from pharmacon.utils.identifiers import validate_mda_artifact_token

    try:
        f = h5py.File(input_file, "r")
    except Exception as exc:
        raise ValidationError(f"Cannot open PSA file: {exc}") from exc

    try:
        attrs = dict(f.attrs)

        # -- artifact_status must be SUCCESS ----------------------------------
        artifact_status = str(attrs.get("artifact_status", "")).strip().upper()
        if not artifact_status:
            raise ValidationError(
                "PSA file is missing 'artifact_status' metadata — "
                "the file may be incomplete or corrupted."
            )
        if artifact_status != "SUCCESS":
            raise ValidationError(
                f"PSA file artifact_status is '{artifact_status}' (expected 'SUCCESS') — "
                "the file may be corrupted or the analysis did not complete."
            )

        # -- artifact_token must be present and valid -------------------------
        artifact_token = str(attrs.get("artifact_token", "")).strip()
        if not artifact_token:
            raise ValidationError(
                "PSA file is missing 'artifact_token' metadata — "
                "the file may be corrupted."
            )

        blueprint = str(attrs.get("blueprint", "")).strip()
        if not blueprint:
            raise ValidationError(
                "PSA file is missing 'blueprint' metadata — "
                "the file may be corrupted."
            )

        token_valid = validate_mda_artifact_token(
            artifact_token=artifact_token,
            blueprint=blueprint,
            secret="structure_analysis",
            namespace="pharmacon",
        )
        if not token_valid:
            raise ValidationError(
                "PSA file artifact_token does not match the blueprint — "
                "the file may have been tampered with or corrupted."
            )

        # -- signature / fingerprint must be present and consistent -----------
        signature = str(attrs.get("signature", "")).strip()
        fingerprint = str(attrs.get("fingerprint", "")).strip()
        command = str(attrs.get("command", "")).strip()
        subcommand = str(attrs.get("subcommand", "")).strip()

        if not signature or not fingerprint:
            raise ValidationError(
                "PSA file is missing 'signature' and/or 'fingerprint' metadata — "
                "the file may be corrupted."
            )

        if command and subcommand:
            expected_sig = create_pharmacon_signature(
                format_name="psa",
                command=command,
                subcommand=subcommand,
            )
            if expected_sig.signature != signature:
                raise ValidationError(
                    f"PSA file signature mismatch.\n"
                    f"  Expected : {expected_sig.signature}\n"
                    f"  Found    : {signature}\n"
                    "The file may have been tampered with or corrupted."
                )
            if expected_sig.fingerprint != fingerprint:
                raise ValidationError(
                    f"PSA file fingerprint mismatch.\n"
                    f"  Expected : {expected_sig.fingerprint}\n"
                    f"  Found    : {fingerprint}\n"
                    "The file may have been tampered with or corrupted."
                )

        # -- pharmacon_version must be present --------------------------------
        version = str(attrs.get("pharmacon_version", "")).strip()
        if not version:
            raise ValidationError(
                "PSA file is missing 'pharmacon_version' metadata — "
                "the file may be corrupted."
            )

    finally:
        f.close()

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Runs the process of printing PSA file content and metadata in a structured
    manner using rich formatting. It includes displaying configuration, validity
    status, file metadata, and tree structure of the PSA file.

    :param args: Parsed command-line arguments containing the configuration for
        the execution, including file input and metadata display preferences.
        This argument must be of type argparse.Namespace.
    :return: None

    """
    from rich.console import Console
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.fileio.base import PharmaconHDF5File
    from pharmacon.logger import header, subheader

    console = Console()

    # Header
    header("Dump PSA")

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

    # Open the PSA file in read-only mode
    with PharmaconHDF5File(args.input, mode="r") as psa:

        # File validity status
        subheader("File Validity")
        psa.print_file_validity_status()

        # Version check
        psa.validate_version(strict=False)

        # File metadata
        if not args.disable_file_metadata:
            subheader("File Metadata")
            psa.print_file_metadata()

        # Tree structure
        subheader("File Structure")
        _print_psa_tree(
            psa,
            group_meta=not args.disable_group_meta,
            dataset_meta=not args.disable_dataset_meta,
            max_items=None if args.verbose else 5,
            file_attrs=args.disable_file_metadata,
        )

        # Per-input metadata for merged files (verbose only)
        if args.verbose:
            _print_input_file_metadata_tables(psa, console)

    # Done
    header("Done")


def _print_psa_tree(psa, *, group_meta: bool, dataset_meta: bool,
                    file_attrs: bool, max_items: int | None = 5) -> None:
    """
    Prints a hierarchical tree structure of an HDF5 file, showing attributes,
    groups, and datasets within the specified HDF5 group/file. The tree is
    visually rich with color-coded elements for better readability.

    :param psa: The HDF5 group or file to be printed.
    :type psa: h5py.Group | h5py.File
    :param group_meta: Whether to include group-level attributes in the tree.
    :type group_meta: bool
    :param dataset_meta: Whether to include dataset-level attributes in the tree.
    :type dataset_meta: bool
    :param file_attrs: Whether to include file-level attributes at the root.
    :type file_attrs: bool
    :param max_items: The maximum number of items to display per group.
        If None, all items are shown.
    :type max_items: int | None
    :return: This function does not return any value.
    :rtype: None
    """
    import h5py
    from rich.console import Console
    from rich.tree import Tree

    console = Console()

    root = Tree(
        f"[bold magenta]{psa.path.name}[/bold magenta]  "
        f"[dim italic]{psa.path.parent}[/dim italic]",
        guide_style="bold bright_blue",
    )

    def add_attrs(tree, attrs) -> None:
        """
        Add attributes to a tree in a formatted string representation.

        This function iterates through the provided attributes and adds each key-value
        pair to the tree as a formatted string. Each attribute is displayed in a
        consistent style. If the string representation of a value exceeds 80 characters,
        it is truncated to 77 characters with an ellipsis.

        :param tree: The tree structure to which attributes will be added.
        :param attrs: A dictionary containing attributes as key-value pairs.
        :return: None
        """
        for k, v in sorted(attrs.items()):
            val_str = str(v)
            if len(val_str) > 80:
                val_str = val_str[:77] + "..."
            tree.add(f"[dim cyan]@{k}[/dim cyan] [dim]=[/dim] [dim white]{val_str}[/dim white]")

    def walk(h5group, tree, depth: int = 0, is_root: bool = False) -> None:
        """
        Recursively walks through an HDF5 group structure, adding its hierarchical representation
        to a tree object. This function explores groups and datasets, formats their metadata,
        and integrates their attributes, if applicable, into the provided tree structure.

        :param h5group: The HDF5 group to explore.
        :param tree: The tree object where the hierarchical representation will be added.
        :param depth: The current depth level of the recursive walk. Defaults to 0.
        :param is_root: A flag indicating whether the current group is the root. Defaults to False.

        :return: None
        """
        show_attrs = (file_attrs or not is_root) and group_meta and h5group.attrs
        if show_attrs:
            add_attrs(tree, h5group.attrs)

        total = len(h5group)
        shown = 0

        for name in h5group:
            if max_items is not None and shown >= max_items:
                break
            obj = h5group[name]

            if isinstance(obj, h5py.Group):
                n_children = len(obj)
                branch = tree.add(
                    f"[bold green]{name}[/bold green]  "
                    f"[dim bright_black]{n_children} "
                    f"item{'s' if n_children != 1 else ''}[/dim bright_black]"
                )
                walk(obj, branch, depth + 1, is_root=False)

            elif isinstance(obj, h5py.Dataset):
                shape_str = "x".join(str(d) for d in obj.shape) if obj.shape else "scalar"
                dnode = tree.add(
                    f"[yellow]{name}[/yellow]  "
                    f"[dim bright_black]{shape_str}[/dim bright_black] "
                    f"[dim]{obj.dtype}[/dim]"
                )
                if dataset_meta and obj.attrs:
                    add_attrs(dnode, obj.attrs)

            shown += 1

        if total > shown:
            remaining = total - shown
            tree.add(
                f"[dim italic]+ {remaining} more "
                f"item{'s' if remaining != 1 else ''} ...[/dim italic]"
            )

    walk(psa.file, root, is_root=True)
    console.print()
    console.print(root)
    console.print()


def _print_input_file_metadata_tables(psa, console) -> None:
    """
    Prints metadata tables for input files found in the given psa object to the provided console.

    This function searches for metadata groups in the psa object that match a specific naming
    pattern, extracts their attributes, and displays them in a formatted table using the rich
    library. Metadata groups that do not match the expected naming convention or do not belong to
    the correct type are ignored.

    :param psa: An object containing an HDF5 file structure with metadata groups.
    :param console: A rich console object used for printing formatted tables.
    :return: None
    """
    import re
    import h5py
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.logger import subheader

    pattern = re.compile(r"^file(\d+)_metadata$")

    matches: list[tuple[int, str]] = []
    for key in psa.file:
        m = pattern.match(key)
        if m is None:
            continue
        obj = psa.file[key]
        if not isinstance(obj, h5py.Group):
            continue
        matches.append((int(m.group(1)), key))

    if not matches:
        return

    matches.sort(key=lambda kv: kv[0])

    subheader("Input Files Metadata")

    for _, group_name in matches:
        grp = psa.file[group_name]

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
