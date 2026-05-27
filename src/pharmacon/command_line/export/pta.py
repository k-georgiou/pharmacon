"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

export pta — Export Pharmacon Trajectory Analysis data to CSV or TSV.
"""

from __future__ import annotations

import argparse
import h5py
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
SUMMARY = "Export Pharmacon Trajectory Analysis (.pta) data to CSV or TSV."

_EPILOG = """\
Examples:

  pharmacon export pta -i analysis.pta -f csv -o ./results/

"""

# Analysis group detection
# Maps known HDF5 top-level group prefixes to their analysis type.
# Each entry: (type, dataset_kind)
#   dataset_kind is used to select the right csv/tsv writer.

_INTERACTION_GROUPS = frozenset({
    "pl_interactions",
    "pp_interactions",
    "hbonds",
})

_DISTANCE_GROUPS = frozenset({
    "distances",
})

_ANGLE_GROUPS = frozenset({
    "angles",
})

_RMSD_GROUPS = frozenset({
    "rmsd",
})

_RMSF_GROUPS = frozenset({
    "rmsf",
})

_PCA_GROUPS = frozenset({
    "pca",
})


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and returns an argparse.ArgumentParser object configured for handling
    the specified subcommand.

    The parser is created as a subparser from the given subparsers object and is
    configured with command-line arguments for input and output options.

    :param subparsers: An argparse._SubParsersAction object to which a new parser
        will be added as a subcommand.
    :param parents: A list of argparse.ArgumentParser objects that serve as
        parent parsers for the new parser. Defaults to None.
    :return: An argparse.ArgumentParser instance configured with the subcommand's
        argument structure.
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
        help="Input Pharmacon Trajectory Analysis (.pta) file. [REQUIRED]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-f", "--format",
        required=True,
        choices=["csv", "tsv"],
        help="Output format: csv or tsv. [REQUIRED]",
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
    Validates command-line arguments and performs various integrity checks and preprocessing to
    ensure the input arguments and associated files meet defined requirements for further operations.

    :param args: Namespace object containing parsed command-line arguments.
    :type args: argparse.Namespace

    :raises ValidationError: If any input validation, metadata integrity, or preprocessing checks fail.
    """
    data: dict[str, object] = vars(args).copy()

    # ── Input file ───────────────────────────────────────────────────────────
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        [".pta"],
    )
    data["input"] = input_file

    # ── Overwrite flag ───────────────────────────────────────────────────────
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # ── Format ───────────────────────────────────────────────────────────────
    fmt = str(data.get("format", "")).strip().lower()
    if fmt not in {"csv", "tsv"}:
        raise ValidationError(
            f"Unsupported export format: '{fmt}'. Supported: csv, tsv."
        )
    data["format"] = fmt

    # ── Output directory ─────────────────────────────────────────────────────
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

    # ── PTA file metadata integrity checks (shared validator) ───────────────
    from pharmacon.utils.pta_validation import validate_pharmacon_file

    attrs = validate_pharmacon_file(input_file, expected_format="pta")
    data["is_merged"] = (attrs.get("is_merged", "False") == "True")

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the primary functionality for exporting PTA analysis data based on the
    provided configuration.

    This function processes analysis data contained in a PTA file and exports it
    into the specified format (CSV or TSV). It includes features for validating the
    PTA file, writing metadata, discovering groups, and exporting specific types
    of data (interactions, distances, angles, etc.). The output is organized into
    relevant files in the specified output directory.

    :param args: Namespace object containing configuration and arguments for
        the function execution, including input file path, output directory,
        export format, overwrite flag, and merging options.
    :type args: argparse.Namespace
    :return: None
    """
    from rich.console import Console
    from rich.table import Table
    from rich.box import SIMPLE

    from pharmacon.fileio.pta import PharmaconPTAFile
    from pharmacon.logger import header, subheader

    console = Console()
    fmt: str = args.format
    ext: str = f".{fmt}"
    output_dir: Path = args.output
    is_merged: bool = args.is_merged

    # Header
    header("Export PTA")

    cfg = Table(show_header=False, box=SIMPLE, pad_edge=False, padding=(0, 1))
    cfg.add_column("Key", style="bold cyan", no_wrap=True, min_width=20)
    cfg.add_column("Value", style="white")
    cfg.add_row("Input", str(args.input))
    cfg.add_row("Format", fmt.upper())
    cfg.add_row("Output directory", str(output_dir))
    cfg.add_row("Overwrite", str(args.overwrite))
    cfg.add_row("Merged file", str(is_merged))
    console.print(cfg)

    # Open the PTA file
    with PharmaconPTAFile(args.input, mode="r") as pta:

        # File validity
        subheader("File Validity")
        pta.print_file_validity_status()
        pta.validate_version(strict=False)

        # Write metadata.json
        subheader("Metadata")
        _write_metadata_json(pta, output_dir)
        console.print(f"  [green]Written:[/green] {output_dir / 'metadata.json'}")

        # Discover analysis groups
        subheader("Discovering Analysis Groups")

        groups = pta.get_groups()
        if not groups:
            console.print("[yellow]No analysis groups found in PTA file.[/yellow]")
            header("Done")
            return

        exported: list[str] = []

        # Interaction groups
        for group_name in sorted(groups):
            if group_name not in _INTERACTION_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            if not is_merged:
                # per-frame interactions
                out_path = output_dir / f"{group_name}_interactions{ext}"
                _check_overwrite(out_path, args.overwrite)

                if fmt == "csv":
                    pta.write_interactions_to_csv(out_path, group_name=group_name)
                else:
                    pta.write_interactions_to_tsv(out_path, group_name=group_name)
                console.print(f"  [green]Written:[/green] {out_path}")
                exported.append(str(out_path))

                # interaction modes
                modes_root = f"{group_name}/modes"
                if modes_root in pta:
                    mode_names = [k for k in pta[modes_root] if k.startswith("mode")]
                    if mode_names:
                        mode_out_files = [
                            output_dir / f"{group_name}_{m}{ext}"
                            for m in sorted(mode_names)
                        ]
                        for p in mode_out_files:
                            _check_overwrite(p, args.overwrite)

                        mode_flags = {
                            "mode1": "mode1" in mode_names,
                            "mode2": "mode2" in mode_names,
                            "mode3": "mode3" in mode_names,
                        }
                        if fmt == "csv":
                            pta.write_interaction_modes_to_csv(
                                mode_out_files, group_name=group_name, **mode_flags,
                            )
                        else:
                            pta.write_interaction_modes_to_tsv(
                                mode_out_files, group_name=group_name, **mode_flags,
                            )
                        for p in mode_out_files:
                            console.print(f"  [green]Written:[/green] {p}")
                            exported.append(str(p))
            else:
                # merged interaction modes (only source present in merged files)
                modes_merged_root = f"{group_name}/modes_merged"
                if modes_merged_root in pta:
                    for mode_name in sorted(pta[modes_merged_root]):
                        out_path = output_dir / f"{group_name}_merged_{mode_name}{ext}"
                        _check_overwrite(out_path, args.overwrite)

                        if fmt == "csv":
                            pta.write_merged_interaction_modes_to_csv(
                                out_path, group_name=group_name, mode_name=mode_name,
                            )
                        else:
                            pta.write_merged_interaction_modes_to_tsv(
                                out_path, group_name=group_name, mode_name=mode_name,
                            )
                        console.print(f"  [green]Written:[/green] {out_path}")
                        exported.append(str(out_path))

        # Distance groups
        for group_name in sorted(groups):
            if group_name not in _DISTANCE_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            out_path = output_dir / f"{group_name}{ext}"
            _check_overwrite(out_path, args.overwrite)

            if fmt == "csv":
                pta.write_distances_to_csv(out_path, group_name=group_name, is_merged=is_merged)
            else:
                pta.write_distances_to_tsv(out_path, group_name=group_name, is_merged=is_merged)
            console.print(f"  [green]Written:[/green] {out_path}")
            exported.append(str(out_path))

            # distance statistics
            stats_path = f"{group_name}/statistics"
            if stats_path in pta:
                out_stats = output_dir / f"{group_name}_statistics{ext}"
                _check_overwrite(out_stats, args.overwrite)

                if fmt == "csv":
                    pta.write_distance_statistics_to_csv(out_stats, group_name=group_name)
                else:
                    pta.write_distance_statistics_to_tsv(out_stats, group_name=group_name)
                console.print(f"  [green]Written:[/green] {out_stats}")
                exported.append(str(out_stats))

        # Angle groups
        for group_name in sorted(groups):
            if group_name not in _ANGLE_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            out_path = output_dir / f"{group_name}{ext}"
            _check_overwrite(out_path, args.overwrite)

            if fmt == "csv":
                pta.write_angles_to_csv(out_path, group_name=group_name, is_merged=is_merged)
            else:
                pta.write_angles_to_tsv(out_path, group_name=group_name, is_merged=is_merged)
            console.print(f"  [green]Written:[/green] {out_path}")
            exported.append(str(out_path))

            # angle statistics
            stats_path = f"{group_name}/statistics"
            if stats_path in pta:
                out_stats = output_dir / f"{group_name}_statistics{ext}"
                _check_overwrite(out_stats, args.overwrite)

                if fmt == "csv":
                    pta.write_angle_statistics_to_csv(out_stats, group_name=group_name)
                else:
                    pta.write_angle_statistics_to_tsv(out_stats, group_name=group_name)
                console.print(f"  [green]Written:[/green] {out_stats}")
                exported.append(str(out_stats))

        # RMSD groups
        for group_name in sorted(groups):
            if group_name not in _RMSD_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            out_path = output_dir / f"{group_name}{ext}"
            _check_overwrite(out_path, args.overwrite)

            if fmt == "csv":
                pta.write_rmsd_data_to_csv(out_path, group_name=group_name, is_merged=is_merged)
            else:
                pta.write_rmsd_data_to_tsv(out_path, group_name=group_name, is_merged=is_merged)
            console.print(f"  [green]Written:[/green] {out_path}")
            exported.append(str(out_path))

            # RMSD statistics
            stats_path = f"{group_name}/statistics"
            if stats_path in pta:
                out_stats = output_dir / f"{group_name}_statistics{ext}"
                _check_overwrite(out_stats, args.overwrite)

                if fmt == "csv":
                    pta.write_rmsd_statistics_to_csv(out_stats, group_name=group_name)
                else:
                    pta.write_rmsd_statistics_to_tsv(out_stats, group_name=group_name)
                console.print(f"  [green]Written:[/green] {out_stats}")
                exported.append(str(out_stats))

        # RMSF groups
        for group_name in sorted(groups):
            if group_name not in _RMSF_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            out_path = output_dir / f"{group_name}{ext}"
            _check_overwrite(out_path, args.overwrite)

            if fmt == "csv":
                pta.write_rmsf_data_to_csv(out_path, group_name=group_name, is_merged=is_merged)
            else:
                pta.write_rmsf_data_to_tsv(out_path, group_name=group_name, is_merged=is_merged)
            console.print(f"  [green]Written:[/green] {out_path}")
            exported.append(str(out_path))

            # RMSF statistics (absent in merged files — decision 2a)
            stats_path = f"{group_name}/statistics"
            if stats_path in pta:
                out_stats = output_dir / f"{group_name}_statistics{ext}"
                _check_overwrite(out_stats, args.overwrite)

                if fmt == "csv":
                    pta.write_rmsf_statistics_to_csv(out_stats, group_name=group_name)
                else:
                    pta.write_rmsf_statistics_to_tsv(out_stats, group_name=group_name)
                console.print(f"  [green]Written:[/green] {out_stats}")
                exported.append(str(out_stats))

        # PCA groups
        for group_name in sorted(groups):
            if group_name not in _PCA_GROUPS:
                continue

            subheader(f"Exporting: {group_name}")

            out_path = output_dir / f"{group_name}{ext}"
            _check_overwrite(out_path, args.overwrite)

            if fmt == "csv":
                pta.write_pca_to_csv(out_path, group_name=group_name)
            else:
                pta.write_pca_to_tsv(out_path, group_name=group_name)
            console.print(f"  [green]Written:[/green] {out_path}")
            exported.append(str(out_path))

    # Summary
    subheader("Summary")
    if exported:
        console.print(f"  Exported [bold]{len(exported)}[/bold] file(s) to [bold]{output_dir}[/bold]")
    else:
        console.print("[yellow]  No exportable analysis groups found in PTA file.[/yellow]")

    # Done
    header("Done")


# Helpers

def _check_overwrite(path: Path, overwrite: bool) -> None:
    """
    Checks whether the specified file path exists and raises an exception if it does
    and overwrite is not enabled.

    :param path: Path to the file that needs to be checked
    :type path: Path
    :param overwrite: Flag indicating whether to overwrite the file if it exists
    :type overwrite: bool
    :raises ValidationError: If the file exists and overwrite is not enabled
    """
    if path.exists() and not overwrite:
        raise ValidationError(
            f"Output file already exists: '{path}'. "
            "Use --overwrite to replace existing files."
        )


def _write_metadata_json(pta, output_dir: Path) -> None:
    """
    Writes metadata from an HDF5 file to a JSON file.

    This function extracts metadata from an HDF5 file, including file-level attributes,
    group-level metadata such as attributes, and counts of datasets and subgroups. The metadata
    is serialized into a JSON file and saved to the specified output directory.

    :param pta: An object containing an HDF5 file (`h5py.File`). It must have a `path` attribute
                representing its file path and a `file` attribute which is an `h5py.File` object.
    :param output_dir: The directory where the metadata JSON file will be written.
    :type output_dir: Path
    :return: This function does not return a value.
    """
    import json
    import h5py

    metadata: dict = {}

    # File-level attributes
    metadata["file"] = {
        "path": str(pta.path),
        "name": pta.path.name,
    }
    metadata["file_metadata"] = {
        str(k): str(v) for k, v in sorted(pta.file.attrs.items())
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

    _collect_group_meta(pta.file)
    metadata["groups"] = groups_meta

    # Write
    out_path = output_dir / "metadata.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
