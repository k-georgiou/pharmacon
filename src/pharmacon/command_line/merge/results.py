"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

merge results — Merge multiple pharmacon trajectory analysis files into one.
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
    validate_logging_level,
    validate_output_file,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "results"
SUMMARY = "Merge multiple pharmacon trajectory analysis files into one."

_EPILOG = """\
Examples:

  pharmacon merge results -i run1.pta run2.pta -o merged.pta

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures an argument parser for the 'merge' subcommand.

    The function is responsible for creating and adding an argument parser
    to the specified subparsers collection for handling the 'merge'
    subcommand. It configures various argument groups including input, output,
    warning, and logging options.

    :param subparsers: A collection of subparser actions to which the 'merge'
                       subcommand parser will be added.
    :type subparsers: argparse._SubParsersAction
    :param parents: A list of argument parsers whose arguments should
                    be included in the created parser. Defaults to None.
    :type parents: List[argparse.ArgumentParser] | None
    :return: Configured argument parser for the 'merge' subcommand.
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

    # Input Options
    inp = parser.add_argument_group("Input Options")
    inp.add_argument(
        "-i", "--input",
        required=True,
        metavar="FILE",
        nargs="+",
        help="Input Pharmacon Trajectory Analysis file(s) to merge. [REQUIRED]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-o", "--output",
        required=False,
        metavar="FILE",
        default="merged.pta",
        help="Output Pharmacon Trajectory Analysis file (default=merged.pta). [OPTIONAL]",
    )
    out.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output file (default=False). [OPTIONAL]",
    )

    # Warning Options
    warn = parser.add_argument_group("Warning Options")
    warn.add_argument(
        "-mw", "--max-warnings",
        required=False,
        metavar="N",
        type=int,
        default=0,
        help="Maximum number of warnings to tolerate before aborting (default: 0). [OPTIONAL]",
    )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="merge_results.log",
                     help="File to log to (default: merge_results.log). [OPTIONAL]", )
    log.add_argument("-fl", "--file-logging-level",
                     required=False,
                     metavar="LEVEL",
                     default="DEBUG",
                     choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
                              "trace", "debug", "info", "warning", "error", "critical"],
                     help="Logging level for file logging (default: DEBUG). [OPTIONAL]", )
    log.add_argument("-tl", "--terminal-logging-level",
                     required=False,
                     metavar="LEVEL",
                     default="INFO",
                     choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
                              "trace", "debug", "info", "warning", "error", "critical"],
                     help="Logging level for terminal logging (default: INFO). [OPTIONAL]", )

    return parser


def _validate_pta_inputs(paths: list[Path]) -> None:
    """
    Validates the Pharmacon Trajectory Analysis (PTA) input files.

    Delegates per-file integrity checks (artifact status, token, signature,
    fingerprint, version ordering, group completion) to the shared
    :func:`validate_pharmacon_file` helper, then enforces merge-specific
    cross-input invariants: every input must share the same command and
    subcommand, and PCA analyses are not mergeable.

    :param paths: List of paths to the Pharmacon Trajectory Analysis input
                  files to validate.
    :raises ValidationError: On any per-file or cross-input failure.
    """
    from pharmacon.utils.pta_validation import validate_pharmacon_file

    first_command: str | None = None
    first_subcommand: str | None = None
    first_path: Path | None = None

    for path in paths:
        # Per-file checks: shared validator (refuses already-merged inputs).
        attrs = validate_pharmacon_file(
            path, expected_format="pta", allow_merged=False,
        )

        command = attrs.get("command", "")
        subcommand = attrs.get("subcommand", "")

        if subcommand.strip().lower() == "pca":
            raise ValidationError(
                f"Input '{path}' is a PCA analysis; merging PCA results is not supported."
            )

        if first_path is None:
            first_command = command
            first_subcommand = subcommand
            first_path = path
            continue

        if command != first_command:
            raise ValidationError(
                f"Input '{path}' has command={command!r}, but '{first_path}' has "
                f"command={first_command!r}. All inputs must share the same command."
            )
        if subcommand != first_subcommand:
            raise ValidationError(
                f"Input '{path}' has subcommand={subcommand!r}, but '{first_path}' has "
                f"subcommand={first_subcommand!r}. All inputs must share the same subcommand."
            )


def validate(args: argparse.Namespace) -> None:
    """
    Validates the provided command-line arguments and ensures that the required
    conditions are met for execution. The validation process includes checking
    file paths, duplicate entries, value ranges, and data types. Updated values
    are written back to the input `argparse.Namespace` object in-place.

    :param args: The parsed command-line arguments.
    :type args: argparse.Namespace
    :raises ValidationError: If the provided arguments fail to meet the validation
        criteria, including file existence, file duplication, invalid flag values,
        and unsupported value types or ranges.
    """
    data: dict[str, object] = vars(args).copy()

    # Input files
    raw_inputs = data.get("input")
    if not isinstance(raw_inputs, (list, tuple)) or len(raw_inputs) == 0:
        raise ValidationError("At least one '--input' file must be provided.")

    inputs: list[Path] = []
    for item in raw_inputs:
        path: Path = validate_existing_input_file(
            item,
            "input",
            {".pta"},
        )
        inputs.append(path)

    # Duplicate input guard
    seen: set[Path] = set()
    for path in inputs:
        if path in seen:
            raise ValidationError(f"Duplicate input file: '{path}'.")
        seen.add(path)
    data["input"] = inputs

    # Per-file PTA content validation (hard errors).
    _validate_pta_inputs(inputs)

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Output file
    output: Path = validate_output_file(
        data.get("output", "merged.pta"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=[".pta"],
    )
    data["output"] = output

    for path in inputs:
        if output == path:
            raise ValidationError(
                f"Output file must not be the same as an input file: '{path}'."
            )

    # Max warnings
    max_warnings_raw = data.get("max_warnings", 0)
    if not isinstance(max_warnings_raw, int) or isinstance(max_warnings_raw, bool):
        raise ValidationError(
            f"Argument 'max_warnings' must be an integer, got {max_warnings_raw!r}."
        )
    if max_warnings_raw < 0:
        raise ValidationError(
            f"Argument 'max_warnings' must be non-negative, got {max_warnings_raw}."
        )
    data["max_warnings"] = max_warnings_raw

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "merge_results.log"),
        "log",
        overwrite=True,
        allowed_suffixes={".log"},
    )
    file_logging_level: str = validate_logging_level(
        data.get("file_logging_level", "TRACE"),
        "file_logging_level",
    )
    terminal_logging_level: str = validate_logging_level(
        data.get("terminal_logging_level", "INFO"),
        "terminal_logging_level",
    )

    data["log"] = log_path
    data["file_logging_level"] = file_logging_level
    data["terminal_logging_level"] = terminal_logging_level

    for path in inputs:
        if log_path == path:
            raise ValidationError(
                f"Log file must not be the same as an input file: '{path}'."
            )
    if log_path == output:
        raise ValidationError("Log file must not be the same as the output file.")

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Execute the merging of PTA (Pharmacon Trajectory Analysis) files based on
    provided arguments. The function processes multiple input PTA files, merges
    data based on specific subcommands, and writes the merged result to an
    output file. It ensures consistency of metadata across input files and
    performs checks to maintain data integrity during the merge operation.

    :param args: Parsed `argparse.Namespace` object containing the command-line
        arguments passed to the script. The expected arguments include:
        - input (List[str]): List of paths to the input PTA files.
        - output (str): Path to the output PTA file.
        - log (str): Path to the log file for capturing runtime logs.
        - max_warnings (int): Maximum number of warnings allowed before aborting
          the merging process.
        - terminal_logging_level (str): Logging level for terminal output
          (e.g., DEBUG, INFO).
        - file_logging_level (str): Logging level for the file output
          (e.g., DEBUG, INFO).

    :raises ValueError: Raised if an unsupported command or subcommand is detected.
    :raises RuntimeError: Raised if the number of warnings exceeds the allowed
        threshold defined by `args.max_warnings`.

    :return: None
    """
    import json
    from collections import defaultdict

    import h5py
    import numpy as np

    from pharmacon.fileio import PTAFile
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.utils.identifiers import create_mda_artifact_token

    # Configure logging for this subcommand
    setup_logger(
        terminal=True,
        terminal_level=args.terminal_logging_level,
        file=True,
        file_level=args.file_logging_level,
        log_file=args.log,
    )
    log = get_logger(__name__)

    # Header
    header("Merge Results")

    for idx, path in enumerate(args.input, start=1):
        log.info("Input #%-3d    : %s", idx, path)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)
    log.info("Max warnings   : %d", args.max_warnings)

    # Warning budget
    warning_state: dict[str, int] = {"total": 0, "max": int(args.max_warnings)}

    def emit_warning(message: str) -> None:
        warning_state["total"] += 1
        log.warning(message)
        if warning_state["total"] > warning_state["max"]:
            raise RuntimeError(
                "Too many warnings; aborting merge. "
                "Use --max-warnings to increase the threshold if you know what you're doing."
            )

    # Read per-file metadata
    # Consistency of command/subcommand is already enforced by validate();
    # blueprint equality is treated as a soft warning.
    commands: List[str] = []
    subcommands: List[str] = []
    blueprints: List[str] = []

    log.debug("Reading per-file metadata from %d input file(s)...", len(args.input))
    for idx, path in enumerate(args.input, start=1):
        with PTAFile(path, mode="r") as pta_in:
            meta = pta_in.file.attrs
            commands.append(str(meta.get("command", "")).strip())
            subcommands.append(str(meta.get("subcommand", "")).strip())
            blueprints.append(str(meta.get("blueprint", "")).strip())
        log.debug(
            "  input #%d: command=%r, subcommand=%r, blueprint=%r",
            idx, commands[-1], subcommands[-1], blueprints[-1],
        )

    command = commands[0]
    subcommand = subcommands[0]
    blueprint = blueprints[0]

    if len(set(blueprints)) != 1:
        emit_warning("Blueprints differ between inputs.")
    else:
        log.debug("All inputs share the same blueprint.")

    normalized_command = command.strip().lower()
    normalized_subcommand = subcommand.strip().lower()

    # Map subcommand → HDF5 group name
    SUBCOMMAND_TO_GROUP: dict[str, str] = {
        "rmsd": "rmsd",
        "rmsf": "rmsf",
        "angles": "angles",
        "distances": "distances",
        "pp-interactions": "pp_interactions",
        "pl-interactions": "pl_interactions",
        "h-bonds": "hbonds",
    }

    if normalized_command != "trajectory analysis":
        raise ValueError(f"Unsupported command for merge: {command}")
    if normalized_subcommand == "pca":
        raise ValueError("PCA merging is not supported.")
    if normalized_subcommand not in SUBCOMMAND_TO_GROUP:
        raise ValueError(f"Unsupported subcommand for merge: {subcommand}")

    group_name = SUBCOMMAND_TO_GROUP[normalized_subcommand]
    log.debug(
        "Resolved dispatch: command=%r subcommand=%r -> group=%r",
        normalized_command, normalized_subcommand, group_name,
    )

    subheader("Merging PTA Files")
    log.info("Command        : %s", command)
    log.info("Subcommand     : %s", subcommand)
    log.info("Group          : %s", group_name)
    log.info("Inputs         : %d file(s)", len(args.input))

    # Helpers
    STANDARD_TIME_KEYS: tuple[str, ...] = ("time_ps", "time_ns", "time_us")

    def merge_per_frame_simple(*,
                               inputs_: list,
                               pta_out_,
                               group: str,
                               dataset: str,
                               key_fields: tuple[str, ...],
                               value_field: str,
                               schema: str,
                               units: str | None = None) -> None:
        """
        Merges data per frame across multiple input datasets while retaining only
        common keys shared across all inputs, and computes statistical summaries such
        as mean and standard deviation for common records. The merged output datasets
        are stored in a specified output group.

        :param inputs_: A list of input datasets containing data to be merged.
        :type inputs_: list
        :param pta_out_: The output dataset handle where merged results will be stored.
        :param group: The group name within the datasets that contains the data frames.
        :type group: str
        :param dataset: The name of the dataset to merge within the group.
        :type dataset: str
        :param key_fields: Tuple of field names used as keys to identify unique records.
        :type key_fields: tuple[str, ...]
        :param value_field: The field name representing the value to compute statistics for.
        :type value_field: str
        :param schema: The schema identifier, used for metadata in the output.
        :type schema: str
        :param units: Optional string representing the units of the data, if applicable.
        :type units: str | None
        :return: None
        """
        frames = list(inputs_[0]._iter_frames(group))
        pta_out_.create_group(group)
        log.debug(
            "merge_per_frame_simple: group=%r dataset=%r frames=%d inputs=%d "
            "key_fields=%s value_field=%r",
            group, dataset, len(frames), len(inputs_),
            list(key_fields), value_field,
        )

        for frame in frames:
            file1_dset_path = f"{group}/frame_{frame}/{dataset}"
            if file1_dset_path not in inputs_[0].file:
                continue

            file1_dset = inputs_[0].file[file1_dset_path]
            file1_attrs = dict(file1_dset.attrs)

            # Build file1's canonical record table for this frame: key -> value.
            file1_records: dict[tuple, float] = {}
            for row in file1_dset:
                rec = json.loads(row)
                key = tuple(str(rec[f]) for f in key_fields)
                file1_records[key] = float(rec[value_field])
            file1_keys = set(file1_records.keys())
            log.trace(
                "  frame %d: file1 has %d record(s)",
                frame, len(file1_records),
            )

            # Each file contributes a {key -> value} map for this frame.
            # file1 goes first so its key order is preserved in the output.
            per_file_records: list[dict[tuple, float]] = [file1_records]
            frame_ok = True

            for idx in range(1, len(inputs_)):
                pta_in = inputs_[idx]
                dset_path = f"{group}/frame_{frame}/{dataset}"
                if dset_path not in pta_in.file:
                    emit_warning(
                        f"{group}: frame {frame} missing in input #{idx + 1}; "
                        f"skipping frame."
                    )
                    frame_ok = False
                    break

                dset = pta_in.file[dset_path]

                for k in STANDARD_TIME_KEYS:
                    v1 = str(file1_attrs.get(k, "")).strip()
                    vi = str(dset.attrs.get(k, "")).strip()
                    if v1 != vi:
                        emit_warning(
                            f"{group}: frame {frame}, input #{idx + 1} "
                            f"has {k}={vi!r} but file1 has {v1!r}."
                        )
                        break

                file_records: dict[tuple, float] = {}
                for row in dset:
                    rec = json.loads(row)
                    key = tuple(str(rec[f]) for f in key_fields)
                    file_records[key] = float(rec[value_field])
                per_file_records.append(file_records)
                log.trace(
                    "  frame %d: input #%d has %d record(s)",
                    frame, idx + 1, len(file_records),
                )

            if not frame_ok:
                continue

            # Strict intersection: a key must exist in every file.
            common_keys = set(file1_keys)
            for file_records in per_file_records[1:]:
                common_keys &= set(file_records.keys())

            dropped_from_file1 = file1_keys - common_keys
            if dropped_from_file1:
                emit_warning(
                    f"{group}: frame {frame}, {len(dropped_from_file1)} key(s) "
                    f"from file1 (e.g. {sorted(dropped_from_file1)[0]}) not "
                    f"present in every input; dropped from merged output."
                )
            for idx in range(1, len(per_file_records)):
                extra = set(per_file_records[idx].keys()) - file1_keys
                if extra:
                    emit_warning(
                        f"{group}: frame {frame}, input #{idx + 1} has "
                        f"{len(extra)} key(s) not in file1 "
                        f"(e.g. {sorted(extra)[0]}); dropped from merged output."
                    )

            log.trace(
                "  frame %d: %d key(s) common across all inputs",
                frame, len(common_keys),
            )

            if not common_keys:
                log.trace("  frame %d: no common keys, skipping write", frame)
                continue

            records = []
            for key in file1_records:
                if key not in common_keys:
                    continue
                vals = [fr[key] for fr in per_file_records]
                arr = np.asarray(vals, dtype=float)
                rec_out: dict = {f: k for f, k in zip(key_fields, key)}
                rec_out["mean"] = float(arr.mean())
                rec_out["std"] = float(arr.std(ddof=0))
                rec_out["n"] = len(vals)
                records.append(rec_out)

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

            frame_group = f"{group}/frame_{frame}"
            pta_out_.create_group(frame_group)

            meta_out: dict[str, str] = {
                "frame_index": str(frame),
                "merged": "True",
                "n_inputs": str(len(inputs_)),
                "schema": schema,
                "format": "json-per-row",
            }
            for k in STANDARD_TIME_KEYS:
                if k in file1_attrs:
                    meta_out[k] = str(file1_attrs[k])
            if units is not None:
                meta_out["units"] = units

            pta_out_.create_dataset(
                group_name=frame_group,
                dataset_name=dataset,
                data=data,
                metadata=meta_out,
            )
            log.trace(
                "  frame %d: wrote %d merged record(s) to %s/%s",
                frame, len(records), frame_group, dataset,
            )

        log.debug(
            "merge_per_frame_simple: finished group=%r (%d frames)",
            group, len(frames),
        )

    def merge_per_selection_simple(*,
                                   inputs_: list,
                                   pta_out_,
                                   group: str,
                                   dataset: str,
                                   key_fields: tuple[str, ...],
                                   value_field: str,
                                   schema: str) -> None:
        """
        Merges per-atom data across multiple input datasets that are organised
        per *selection* rather than per *frame*. Used by analyses with no time
        axis (e.g. RMSF).

        For each selection label discovered in ``inputs_[0]``, this:
        - reads ``/<group>/selection_<label>/<dataset>`` from every input;
        - keys records by the supplied ``key_fields`` (full identity tuple);
        - takes the strict intersection of keys across all inputs (records
          present in only some inputs are dropped with a warning);
        - for each surviving key, computes mean / std / n over the
          ``value_field`` and writes a merged record carrying the original
          ``key_fields`` (numeric types preserved from file1's row).

        Per-selection dataset attrs (``selection_string``, ``fitting_group``,
        ``frame_begin``, ``frame_end``, ``frame_step``) are copied from
        file1; mismatches across inputs are warned but not blocking, matching
        the policy used by :func:`merge_per_frame_simple` for time keys.

        :param inputs_: List of opened PTA input handles.
        :param pta_out_: Output PTA handle to write into.
        :param group: Parent HDF5 group (e.g. ``"rmsf"``).
        :param dataset: Dataset name within each selection (e.g. ``"atoms"``).
        :param key_fields: Tuple of field names identifying an atom uniquely
            across files.
        :param value_field: Field whose value is aggregated (e.g. ``"rmsf"``).
        :param schema: Schema string written to the merged dataset attrs.
        :return: None
        """
        labels = list(inputs_[0]._iter_selections(group))
        pta_out_.create_group(group)
        log.debug(
            "merge_per_selection_simple: group=%r dataset=%r selections=%d "
            "inputs=%d key_fields=%s value_field=%r",
            group, dataset, len(labels), len(inputs_),
            list(key_fields), value_field,
        )

        PROPAGATE_ATTRS: tuple[str, ...] = (
            "selection_string", "fitting_group",
            "frame_begin", "frame_end", "frame_step",
        )

        for label in labels:
            file1_path = f"{group}/selection_{label}/{dataset}"
            if file1_path not in inputs_[0].file:
                continue

            file1_dset = inputs_[0].file[file1_path]
            file1_attrs = dict(file1_dset.attrs)

            # Parse file1's records, keeping the full typed record so we can
            # preserve numeric key types (atom_index, resid) on output.
            file1_records: dict[tuple, float] = {}
            file1_typed: dict[tuple, dict] = {}
            for row in file1_dset:
                rec = json.loads(row)
                key = tuple(str(rec[f]) for f in key_fields)
                file1_records[key] = float(rec[value_field])
                file1_typed[key] = rec
            file1_keys = set(file1_records.keys())
            log.trace(
                "  selection %r: file1 has %d record(s)",
                label, len(file1_records),
            )

            per_file_records: list[dict[tuple, float]] = [file1_records]
            sel_ok = True

            for idx in range(1, len(inputs_)):
                pta_in = inputs_[idx]
                dset_path = f"{group}/selection_{label}/{dataset}"
                if dset_path not in pta_in.file:
                    emit_warning(
                        f"{group}: selection {label!r} missing in input "
                        f"#{idx + 1}; skipping selection."
                    )
                    sel_ok = False
                    break

                dset = pta_in.file[dset_path]

                # Sanity-check propagated attrs (warn, do not block).
                for k in PROPAGATE_ATTRS:
                    v1 = str(file1_attrs.get(k, "")).strip()
                    vi = str(dset.attrs.get(k, "")).strip()
                    if v1 != vi:
                        emit_warning(
                            f"{group}: selection {label!r}, input "
                            f"#{idx + 1} has {k}={vi!r} but file1 has {v1!r}."
                        )

                file_records: dict[tuple, float] = {}
                for row in dset:
                    rec = json.loads(row)
                    key = tuple(str(rec[f]) for f in key_fields)
                    file_records[key] = float(rec[value_field])
                per_file_records.append(file_records)
                log.trace(
                    "  selection %r: input #%d has %d record(s)",
                    label, idx + 1, len(file_records),
                )

            if not sel_ok:
                continue

            # Strict intersection: a key must exist in every file.
            common_keys = set(file1_keys)
            for file_records in per_file_records[1:]:
                common_keys &= set(file_records.keys())

            dropped_from_file1 = file1_keys - common_keys
            if dropped_from_file1:
                emit_warning(
                    f"{group}: selection {label!r}, "
                    f"{len(dropped_from_file1)} key(s) from file1 "
                    f"(e.g. {sorted(dropped_from_file1)[0]}) not present in "
                    f"every input; dropped from merged output."
                )
            for idx in range(1, len(per_file_records)):
                extra = set(per_file_records[idx].keys()) - file1_keys
                if extra:
                    emit_warning(
                        f"{group}: selection {label!r}, input #{idx + 1} "
                        f"has {len(extra)} key(s) not in file1 "
                        f"(e.g. {sorted(extra)[0]}); dropped from merged output."
                    )

            log.trace(
                "  selection %r: %d key(s) common across all inputs",
                label, len(common_keys),
            )

            if not common_keys:
                log.trace("  selection %r: no common keys, skipping write", label)
                continue

            records = []
            for key in file1_records:
                if key not in common_keys:
                    continue
                vals = [fr[key] for fr in per_file_records]
                arr = np.asarray(vals, dtype=float)
                # Recover original typed key fields from file1's parsed row.
                rec_out: dict = {f: file1_typed[key][f] for f in key_fields}
                rec_out["mean"] = float(arr.mean())
                rec_out["std"] = float(arr.std(ddof=0))
                rec_out["n"] = len(vals)
                records.append(rec_out)

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

            sel_group_out = f"{group}/selection_{label}"
            pta_out_.create_group(sel_group_out)

            meta_out: dict[str, str] = {
                "label": label,
                "merged": "True",
                "n_inputs": str(len(inputs_)),
                "n_atoms": str(len(records)),
                "schema": schema,
                "format": "json-per-row",
            }
            for k in PROPAGATE_ATTRS:
                if k in file1_attrs:
                    meta_out[k] = str(file1_attrs[k])

            pta_out_.create_dataset(
                group_name=sel_group_out,
                dataset_name=dataset,
                data=data,
                metadata=meta_out,
            )
            log.trace(
                "  selection %r: wrote %d merged record(s) to %s/%s",
                label, len(records), sel_group_out, dataset,
            )

        log.debug(
            "merge_per_selection_simple: finished group=%r (%d selections)",
            group, len(labels),
        )

    def merge_interaction_modes(*,
                                inputs_: list,
                                pta_out_,
                                group: str,
                                mode_names: tuple[str, ...]) -> None:
        """
        Merges interaction modes from multiple input datasets into a single output dataset.
        The merged dataset is stored in the specified group within the `pta_out_`
        object, creating a union of keys across all inputs and computing aggregated
        statistics for each key.

        The function handles missing data by treating missing files or tables as
        empty, and will emit warnings for unexpected or malformed keys encountered
        in the input datasets.

        For each requested mode this:

        - reads ``/<group>/modes/<mode>/table`` from every input and indexes
          it by ``((r1 tuple), (r2 tuple), label)`` — the native key format
          emitted by :meth:`PharmaconPTAFile.build_interaction_modes`;
        - unions the keys across all inputs (missing entries count as 0.0
          frequency);
        - computes avg and std of ``frequency`` per key;
        - writes the merged table under ``/<group>/modes_merged/<mode>/table``
          with schema ``{key: [[r1], [r2], label], mean_frequency,
          std_frequency, n_files}`` — the exact schema consumed by
          :meth:`PharmaconPTAFile.read_merged_interactions` and the merged
          interaction export writers.

        :param inputs_: List of input datasets to be merged. Each entry is expected
            to provide interaction tables in the specified group and modes.
        :type inputs_: list
        :param pta_out_: Output dataset where the merged interaction modes will be
            stored.
        :param group: Name of the group where interaction modes will be located.
        :type group: str
        :param mode_names: Tuple of mode names to be merged under the specified group.
        :type mode_names: tuple[str, ...]
        :return: This function does not return any value.
        :rtype: None
        """
        pta_out_.create_group(group)
        pta_out_.create_group(f"{group}/modes_merged")

        n_inputs = len(inputs_)
        log.debug(
            "merge_interaction_modes: group=%r modes=%s inputs=%d",
            group, list(mode_names), n_inputs,
        )

        for mode_name in mode_names:
            log.debug("  merging mode %r...", mode_name)
            per_file_tables: list[dict[tuple, float]] = []

            for idx, pta_in in enumerate(inputs_):
                table_path = f"{group}/modes/{mode_name}/table"
                table: dict[tuple, float] = {}

                if table_path not in pta_in.file:
                    emit_warning(
                        f"{group}: input #{idx + 1} is missing '{table_path}'; treating as empty."
                    )
                    per_file_tables.append(table)
                    continue

                dset = pta_in.file[table_path]
                for row in dset:
                    rec = json.loads(row)
                    key_obj = rec.get("key")
                    if (
                        not isinstance(key_obj, (list, tuple))
                        or len(key_obj) != 3
                    ):
                        emit_warning(
                            f"{group}/{mode_name}: input #{idx + 1} has an "
                            f"unexpected interaction key {key_obj!r}; skipping."
                        )
                        continue
                    r1_raw, r2_raw, label = key_obj
                    r1 = tuple(r1_raw) if isinstance(r1_raw, (list, tuple)) else (str(r1_raw),)
                    r2 = tuple(r2_raw) if isinstance(r2_raw, (list, tuple)) else (str(r2_raw),)
                    key = (r1, r2, str(label))
                    frequency = float(rec.get("frequency", 0.0))
                    table[key] = frequency

                per_file_tables.append(table)
                log.trace(
                    "    mode %r: input #%d has %d row(s)",
                    mode_name, idx + 1, len(table),
                )

            all_keys: set[tuple] = set()
            for table in per_file_tables:
                all_keys.update(table.keys())
            log.debug(
                "    mode %r: union across inputs = %d key(s)",
                mode_name, len(all_keys),
            )

            records = []
            for key in sorted(all_keys):
                freqs = [table.get(key, 0.0) for table in per_file_tables]
                freqs_arr = np.asarray(freqs, dtype=float)

                r1, r2, label = key
                records.append({
                    "key": [list(r1), list(r2), label],
                    "mean_frequency": float(freqs_arr.mean()),
                    "std_frequency": float(freqs_arr.std(ddof=0)),
                    "n_files": n_inputs,
                })

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

            mode_group = f"{group}/modes_merged/{mode_name}"
            pta_out_.create_group(mode_group)
            pta_out_.create_dataset(
                group_name=mode_group,
                dataset_name="table",
                data=data,
                metadata={
                    "n_files": str(n_inputs),
                    "n_rows": str(len(data)),
                    "merged": "True",
                    "missing_policy": "zero-fill",
                    "schema": (
                        "{key: [[r1], [r2], label], "
                        "mean_frequency, std_frequency, n_files}"
                    ),
                    "format": "json-per-row",
                },
            )
            log.debug(
                "    mode %r: wrote %d merged row(s) to %s/table",
                mode_name, len(records), mode_group,
            )

        log.debug(
            "merge_interaction_modes: finished group=%r",
            group,
        )

    # Merge
    log.debug("Opening output PTA file %s (overwrite=%s)", args.output, args.overwrite)
    with PTAFile(args.output, overwrite=args.overwrite,
                 command=command, subcommand=subcommand) as pta_out:

        pta_out.add_file_metadata(
            metadata={
                "description": "Merged trajectory analysis",
                "blueprint": blueprint,
                "is_merged": "True",
                "merge_strategy": "frame-wise avg+std",
                "merge_missing_policy": "zero-fill",
                "n_inputs": str(len(args.input)),
                "input_files": json.dumps([str(p) for p in args.input]),
            },
            overwrite=True,
        )

        log.debug("Opening %d input PTA file(s) for merge...", len(args.input))
        inputs = [PTAFile(f, mode="r") for f in args.input]
        try:
            log.debug("Dispatching merge for subcommand %r", normalized_subcommand)
            match normalized_subcommand:
                case "rmsd":
                    merge_per_frame_simple(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        dataset="rmsd",
                        key_fields=("label",),
                        value_field="value",
                        schema="{label, mean, std, n}",
                    )
                case "rmsf":
                    merge_per_selection_simple(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        dataset="atoms",
                        key_fields=("atom_index", "resid", "resname", "atom_name"),
                        value_field="rmsf",
                        schema="{atom_index, resid, resname, atom_name, mean, std, n}",
                    )
                case "angles":
                    merge_per_frame_simple(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        dataset="angles",
                        key_fields=("label", "kind"),
                        value_field="value",
                        schema="{label, kind, mean, std, n}",
                        units="degrees",
                    )
                case "distances":
                    merge_per_frame_simple(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        dataset="distances",
                        key_fields=("label", "method"),
                        value_field="distance",
                        schema="{label, method, mean, std, n}",
                    )
                case "pp-interactions" | "pl-interactions":
                    merge_interaction_modes(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        mode_names=("mode1", "mode2", "mode3"),
                    )
                case "h-bonds":
                    merge_interaction_modes(
                        inputs_=inputs,
                        pta_out_=pta_out,
                        group=group_name,
                        mode_names=("mode1", "mode2"),
                    )
                case _:
                    raise ValueError(f"Unsupported subcommand for merge: {subcommand}")
        finally:
            log.debug("Closing %d input PTA file(s)", len(inputs))
            for pta_in in inputs:
                pta_in.close()

        # Per-input file metadata (hidden reference groups)
        # Each input's original file-level attrs are preserved under
        # ``/file{N}_metadata`` so downstream consumers (or users) can trace
        # back to the individual runs. These groups carry ``completed=True``
        # to satisfy the export validator and do not match any known analysis
        # group name, so the export module ignores them. The dump module
        # surfaces them as tables when run with ``--verbose``.
        for idx, path in enumerate(args.input, start=1):
            meta_group = f"file{idx}_metadata"
            pta_out.create_group(meta_group)
            with PTAFile(path, mode="r") as pta_in:
                source_attrs: dict[str, str] = {
                    str(k): str(v) for k, v in pta_in.file.attrs.items()
                }
            source_attrs["source_index"] = str(idx)
            source_attrs["source_path"] = str(path)
            source_attrs["completed"] = "True"
            pta_out.add_group_metadata(
                group_name=meta_group,
                metadata=source_attrs,
                overwrite=True,
            )
            log.debug(
                "Wrote hidden group %s (%d attrs) for input #%d: %s",
                meta_group, len(source_attrs), idx, path,
            )

        # Final success metadata + artifact token (MUST BE LAST)
        log.debug("Writing final success metadata to output file")
        pta_out.add_file_metadata(
            metadata={
                "artifact_status": "SUCCESS",
                "artifact_status_code": 0,
                "completed": "True",
            },
            overwrite=True,
        )
        log.debug("Generating artifact token for merged output")
        token = create_mda_artifact_token(
            blueprint=blueprint,
            secret="trajectory_analysis",
            namespace="pharmacon",
        )
        pta_out.add_file_metadata(
            {"artifact_token": token, "artifact_token_version": "1"},
            overwrite=True,
        )
        log.debug("Artifact token written; output file finalized")
        pta_out.add_group_metadata(
            group_name=group_name,
            metadata={"completed": "True"},
        )

    header("Done")
    log.info("Merge complete. %d warning(s).", warning_state["total"])
    log.info("Results written to: %s", args.output)
