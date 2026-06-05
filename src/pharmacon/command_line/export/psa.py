"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

export psa — Export Pharmacon Structure Analysis data to CSV, TSV, or FASTA.
"""

from __future__ import annotations

import argparse
import h5py
from pathlib import Path
from typing import List

from pharmacon.command_line.exceptions import ValidationError
from pharmacon.command_line.formatter import make_formatter_class
from pharmacon.utils.validation import validate_bool_flag, validate_existing_input_file



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "psa"
SUMMARY = "Export Pharmacon Structure Analysis (.psa) data to CSV, TSV, or FASTA."

_EPILOG = """\
Examples:

  pharmacon export psa -i properties.psa -f csv   -o ./results/

"""

_SUPPORTED_FORMATS: frozenset[str] = frozenset({"csv", "tsv", "fasta"})

# Analysis group detection
# Maps known PSA top-level group names to their analysis type. Each entry is
# dispatched to a dedicated writer inside :func:`run`.

_PROPERTIES_GROUPS = frozenset({
    "properties",
})

_SEQUENCE_GROUPS = frozenset({
    "sequence",
})


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures a command-line argument parser for subcommands. This parser
    specifically handles input and output options related to Pharmacon Structure
    Analysis files, allowing the user to specify required input file, output format,
    and output directory.

    :param subparsers: Subparsers action object to which the new subcommand parser
        is added.
    :param parents: Optional; List of parent parsers that provide default arguments
        for this parser. Defaults to None, indicating no parent parsers.
    :return: A fully configured argparse.ArgumentParser instance ready for parsing
        arguments for this subcommand.
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
        help="Input Pharmacon Structure Analysis (.psa) file. [REQUIRED]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-f", "--format",
        required=True,
        choices=["csv", "tsv", "fasta"],
        help="Output format: csv, tsv, or fasta (sequence only). [REQUIRED]",
    )
    out.add_argument(
        "-o", "--output",
        required=True,
        metavar="DIR",
        help="Output directory for exported files. [REQUIRED]",
    )
    out.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output files. Default='False'. Type='bool' [OPTIONAL]",
    )

    return parser


def validate(args: argparse.Namespace) -> None:
    """
    Validates the given argparse namespace for correctness and sets validated
    values into the namespace. This function ensures that the required input,
    output, and format-related options are validated to ensure proper processing.

    Additionally, it performs comprehensive checks on the metadata integrity
    for PSA (Pharmacodynamics Structural Analysis) files and validates whether
    the export format is compatible with the data contained in the input file.

    :param args: The namespace object containing parsed command-line arguments.
                 This namespace will be modified in-place with validated and
                 required values.
    :type args: argparse.Namespace
    :raises ValidationError: If any validation checks fail, such as missing
                             required input data, invalid paths, unsupported
                             formats, or metadata inconsistencies in PSA files.
    """
    data: dict[str, object] = vars(args).copy()

    # Input file
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        [".psa"],
    )
    data["input"] = input_file

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Format
    fmt = str(data.get("format", "")).strip().lower()
    if fmt not in _SUPPORTED_FORMATS:
        raise ValidationError(
            f"Unsupported export format: '{fmt}'. "
            f"Supported: {', '.join(sorted(_SUPPORTED_FORMATS))}."
        )
    data["format"] = fmt

    # Output directory
    import shutil

    output_dir = Path(str(data.get("output", ""))).expanduser().resolve()

    if output_dir.exists():
        if output_dir.is_dir():
            if overwrite:
                shutil.rmtree(output_dir)
                output_dir.mkdir(parents=True)
        else:
            raise ValidationError(
                f"Output path is not a directory: '{output_dir}'."
            )
    else:
        output_dir.mkdir(parents=True)

    data["output"] = output_dir

    # PSA file metadata integrity checks (shared validator)
    from pharmacon.utils.pta_validation import validate_pharmacon_file

    attrs = validate_pharmacon_file(input_file, expected_format="psa")
    data["is_merged"] = (attrs.get("is_merged", "False") == "True")

    # FASTA is sequence-only — reject if the file contains properties groups.
    if fmt == "fasta":
        with h5py.File(input_file, "r") as fh:
            top_level_groups = {
                key for key in fh if isinstance(fh[key], h5py.Group)
            }
        properties_present = top_level_groups & _PROPERTIES_GROUPS
        if properties_present:
            raise ValidationError(
                f"Format 'fasta' is not supported for PSA files containing "
                f"properties groups (found: {sorted(properties_present)}). "
                f"FASTA is a sequence-only format — use 'csv' or 'tsv' to "
                f"export properties, or point '-i' at a sequence-only PSA file."
            )

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the main export process for the PSA file based on the provided arguments.

    This function orchestrates the process of exporting data from a PSA file,
    including validating the file, writing metadata, discovering analysis groups,
    and exporting the content to the specified formats and directory. It provides
    a structured output of steps in the process for a clear and comprehensive
    understanding of the operations performed.

    :param args: The Namespace object from argparse containing the parsed command-line arguments.
    :type args: argparse.Namespace
    """
    from rich.console import Console
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.fileio.psa import PharmaconPSAFile
    from pharmacon.logger import header, subheader

    console = Console()
    fmt: str = args.format
    delimiter: str = "\t" if fmt == "tsv" else ","
    output_dir: Path = args.output
    is_merged: bool = args.is_merged

    # Header
    header("Export PSA")

    cfg = Table(show_header=False, box=SIMPLE, pad_edge=False, padding=(0, 1))
    cfg.add_column("Key", style="bold cyan", no_wrap=True, min_width=20)
    cfg.add_column("Value", style="white")
    cfg.add_row("Input", str(args.input))
    cfg.add_row("Format", fmt.upper())
    cfg.add_row("Output directory", str(output_dir))
    cfg.add_row("Overwrite", str(args.overwrite))
    cfg.add_row("Merged file", str(is_merged))
    console.print(cfg)

    # Open the PSA file
    with PharmaconPSAFile(args.input, mode="r") as psa:

        # File validity
        subheader("File Validity")
        psa.print_file_validity_status()
        psa.validate_version(strict=False)

        # Write metadata.json
        subheader("Metadata")
        _write_metadata_json(psa, output_dir)
        console.print(f"  [green]Written:[/green] {output_dir / 'metadata.json'}")

        # Discover analysis groups
        subheader("Discovering Analysis Groups")

        groups = psa.get_groups()
        if not groups:
            console.print("[yellow]No analysis groups found in PSA file.[/yellow]")
            header("Done")
            return

        exported: list[str] = []

        # Properties groups
        for group_name in sorted(groups):
            if group_name not in _PROPERTIES_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            # Pre-flight overwrite checks so we fail before writing anything
            base_path = output_dir / group_name
            stem = base_path.name
            parent = base_path.parent
            ext = "tsv" if fmt == "tsv" else "csv"

            scalars_path = parent / f"{stem}_scalars.{ext}"
            meta_path = parent / f"{stem}.meta.json"
            fp_paths = [
                parent / f"{stem}_fp_{kind}.parquet"
                for kind in ("morgan", "maccs", "atom_pair", "topological_torsion")
            ]
            for p in [scalars_path, meta_path, *fp_paths]:
                _check_overwrite(p, args.overwrite)

            outputs = psa.write_ml_ready_export(
                base_path=base_path,
                group_name=group_name,
                overwrite=args.overwrite,
                include_fingerprints=True,
                fingerprint_kinds=("morgan", "maccs", "atom_pair", "topological_torsion"),
                delimiter=delimiter,
            )
            for label, path in outputs.items():
                console.print(f"  [green]Written:[/green] ({label}) {path}")
                exported.append(str(path))

        # Sequence groups
        for group_name in sorted(groups):
            if group_name not in _SEQUENCE_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            if fmt == "fasta":
                out_path = output_dir / f"{group_name}.fasta"
                _check_overwrite(out_path, args.overwrite)
                psa.write_sequence_fasta(
                    out_path,
                    group_name=group_name,
                    overwrite=args.overwrite,
                )
                console.print(f"  [green]Written:[/green] {out_path}")
                exported.append(str(out_path))
            else:
                ext = "tsv" if fmt == "tsv" else "csv"
                out_path = output_dir / f"{group_name}.{ext}"
                _check_overwrite(out_path, args.overwrite)
                if fmt == "csv":
                    psa.write_sequence_to_csv(
                        out_path,
                        group_name=group_name,
                        overwrite=args.overwrite,
                    )
                else:
                    psa.write_sequence_to_tsv(
                        out_path,
                        group_name=group_name,
                        overwrite=args.overwrite,
                    )
                console.print(f"  [green]Written:[/green] {out_path}")
                exported.append(str(out_path))

    # Summary
    subheader("Summary")
    if exported:
        console.print(f"  Exported [bold]{len(exported)}[/bold] file(s) to [bold]{output_dir}[/bold]")
    else:
        console.print("[yellow]  No exportable analysis groups found in PSA file.[/yellow]")

    # Done
    header("Done")


# Helpers

def _check_overwrite(path: Path, overwrite: bool) -> None:
    """
    Validates whether an existing file can be overwritten.

    This function checks if the file at the specified path already exists. If
    the file exists and the 'overwrite' flag is set to False, a ValidationError
    is raised, notifying the user that the file cannot be overwritten without
    explicitly setting the 'overwrite' flag to True.

    :param path: The path to the file that needs to be checked.
    :type path: Path
    :param overwrite: A flag indicating if overwriting the file is permitted.
    :type overwrite: bool
    :return: None
    """
    if path.exists() and not overwrite:
        raise ValidationError(
            f"Output file already exists: '{path}'. "
            "Use --overwrite to replace existing files."
        )


def _write_metadata_json(psa, output_dir: Path) -> None:
    """
    Writes metadata from an HDF5 file to a JSON file.

    This function extracts metadata from an HDF5 file, including file-level attributes
    and group-level metadata, and serializes it into a JSON file in the specified output directory.

    :param psa: HDF5 processing object containing references to the file and its metadata.
    :type psa: h5py.File
    :param output_dir: Directory where the JSON metadata file will be written.
    :type output_dir: Path
    :return: None
    """
    import json
    import h5py

    metadata: dict = {}

    # File-level attributes
    metadata["file"] = {
        "path": str(psa.path),
        "name": psa.path.name,
    }
    metadata["file_metadata"] = {
        str(k): str(v) for k, v in sorted(psa.file.attrs.items())
    }

    # Group-level metadata
    groups_meta: dict = {}

    def _collect_group_meta(h5group, prefix: str = "") -> None:
        for name, obj in sorted(h5group.items()):
            full_path = f"{prefix}/{name}" if prefix else name
            if isinstance(obj, h5py.Group):
                entry: dict = {}
                if obj.attrs:
                    entry["attributes"] = {
                        str(k): str(v) for k, v in sorted(obj.attrs.items())
                    }
                # count datasets and subgroups
                n_datasets = sum(1 for v in obj.values() if isinstance(v, h5py.Dataset))
                n_subgroups = sum(1 for v in obj.values() if isinstance(v, h5py.Group))
                entry["n_datasets"] = n_datasets
                entry["n_subgroups"] = n_subgroups
                groups_meta[full_path] = entry
                _collect_group_meta(obj, full_path)

    _collect_group_meta(psa.file)
    metadata["groups"] = groups_meta

    # Write
    out_path = output_dir / "metadata.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
