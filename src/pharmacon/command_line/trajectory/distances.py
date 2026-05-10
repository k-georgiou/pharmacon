"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory distances — Calculate distances between atom groups through a trajectory.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

from pharmacon.command_line.exceptions import ValidationError
from pharmacon.command_line.formatter import make_formatter_class
from pharmacon.constants.formats import (supported_trajectory_analysis_topology_formats,
                                         supported_trajectory_analysis_trajectory_formats)
from pharmacon.utils.validation import (
    validate_bool_flag,
    validate_existing_input_file,
    validate_frame_range,
    validate_logging_level,
    validate_output_file,
    validate_string_list,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "distances"
SUMMARY = "Calculate distances between atom groups through a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory distances -p topol.tpr -x traj.xtc -o distances.pta -sel1 "resname LIG" -sel2 "protein and name CA" -m COM -n lig_to_ca

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures an argument parser for the specified subcommand.

    The function creates an argument parser for handling the input, output, selection,
    time range, periodic boundary condition (PBC), and logging options used in the
    Pharmacon Trajectory Analysis process. This parser is built with detailed
    configuration to parse command-line arguments pertaining to the execution of the
    analysis.

    :param subparsers: A subparsers action object for attaching the subcommand parser.
    :param parents: A list of parent parsers, if any, to inherit arguments from (default: None).
    :return: The configured argument parser for the defined subcommand.
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
        "-p", "--topology",
        required=True,
        metavar="FILE",
        help="Input topology file. [REQUIRED]",
    )
    inp.add_argument(
        "-x", "--trajectory",
        required=True,
        metavar="FILE",
        help="Input trajectory file. [REQUIRED]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-o", "--output",
        required=False,
        metavar="FILE",
        default="distances.pta",
        help="Output Pharmacon Trajectory Analysis file (default=distances.pta). [OPTIONAL]",
    )

    out.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output file (default=False). [OPTIONAL]",
    )

    # Selection Options
    sel = parser.add_argument_group("Selection Options")
    sel.add_argument(
        "-sel1", "--selection1",
        metavar="SEL",
        default=None,
        required=True,
        nargs="+",
        help="First selection group(s) for distance calculation. [REQUIRED]",
    )

    sel.add_argument(
        "-sel2", "--selection2",
        metavar="SEL",
        default=None,
        required=True,
        nargs="+",
        help="Second selection group(s) for distance calculation. [REQUIRED]",
    )

    sel.add_argument(
        "-m", "--method",
        metavar="METHOD",
        default=None,
        required=True,
        nargs="+",
        choices=["MIN", "MAX", "COM", "COG",
                 "min", "max", "com", "cog"],
        help="Distance method per pair: MIN (minimum), MAX (maximum), COM (center of mass), COG (center of geometry). [REQUIRED]",
    )

    sel.add_argument(
        "-n", "--names",
        metavar="NAME",
        default=None,
        required=True,
        nargs="+",
        help="Name(s) to assign to each distance measurement. [REQUIRED]",
    )

    # Time Range Options
    tr = parser.add_argument_group("Time Range Options")
    tr.add_argument(
        "-b", "--begin",
        metavar="FRAME",
        type=int,
        default=0,
        required=False,
        help="Start frame (default: 0). [OPTIONAL]")

    tr.add_argument(
        "-e", "--end",
        metavar="FRAME",
        type=int,
        default=None,
        required=False,
        help="End frame; None means last frame (default: None). [OPTIONAL]",
    )

    tr.add_argument(
        "-s", "--step",
        metavar="N",
        type=int,
        default=1,
        help="Iterate every Nth frame (default: 1). [OPTIONAL]")

    # PBC Options
    pbc = parser.add_argument_group("PBC Options")
    pbc.add_argument("-at", "--add-transformations",
                     required=False,
                     action="store_true",
                     help="Add transformations to MDAnalysis Universe (default: False). [OPTIONAL]. See documentation for more details.", )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="distances.log",
                     help="File to log to (default: distances.log). [OPTIONAL]", )
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


def validate(args: argparse.Namespace) -> None:
    """
    Validates the input arguments provided via an argparse.Namespace object for consistency, correctness,
    and compliance with expected formats. This function verifies various inputs such as file paths,
    boolean flags, selection lists, methods, logging levels, and frame range, ensuring cohesive and
    valid configurations are maintained. It also throws validation errors for improper combinations
    or formats and writes validated values back to the given namespace.

    :param args: Input arguments encapsulated in an argparse.Namespace object. Arguments include:
        - `topology`: Input topology file path, validated for existence and supported formats.
        - `trajectory`: Input trajectory file path, validated for existence and supported formats.
        - `output`: Output file path, validated for overwrite flag and supported suffixes.
        - `overwrite`: Boolean flag indicating whether existing files can be overwritten.
        - `selection1`: List of selection strings for the first group.
        - `selection2`: List of selection strings for the second group.
        - `method`: List containing distance calculation methods, validated against allowed values.
        - `names`: List of names corresponding to pairwise distance calculations.
        - Other optional logging configurations and transformations.
    :type args: argparse.Namespace

    :raises ValidationError: Raised in cases of invalid inputs, mismatched argument combinations, unsupported
        methods, or inconsistent list lengths among related arguments.

    :return: This function modifies the given argparse.Namespace object in-place and does not return any value.
    """
    data: dict[str, object] = vars(args).copy()

    # Input files
    topology: Path = validate_existing_input_file(
        data.get("topology"),
        "topology",
        supported_trajectory_analysis_topology_formats,
    )
    trajectory: Path = validate_existing_input_file(
        data.get("trajectory"),
        "trajectory",
        supported_trajectory_analysis_trajectory_formats,
    )

    if topology == trajectory:
        raise ValidationError(
            "Topology and trajectory must be different files."
        )

    data["topology"] = topology
    data["trajectory"] = trajectory

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Output file
    output: Path = validate_output_file(
        data.get("output", "distances.pta"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=[".pta"],
    )
    data["output"] = output

    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selections, methods, and names
    selection1: list[str] = validate_string_list(data.get("selection1"), "selection1")
    selection2: list[str] = validate_string_list(data.get("selection2"), "selection2")
    names: list[str] = validate_string_list(data.get("names"), "names")

    raw_method = data.get("method")
    if not isinstance(raw_method, list) or not raw_method:
        raise ValidationError(
            "Argument 'method' requires at least one value."
        )
    valid_methods = {"MIN", "MAX", "COM", "COG"}
    method: list[str] = []
    for i, m in enumerate(raw_method):
        if not isinstance(m, str):
            raise ValidationError(
                f"Method at position {i + 1} is not a string."
            )
        upper = m.strip().upper()
        if upper not in valid_methods:
            raise ValidationError(
                f"Invalid method at position {i + 1}: '{m}'. "
                f"Valid methods: {', '.join(sorted(valid_methods))}."
            )
        method.append(upper)

    data["selection1"] = selection1
    data["selection2"] = selection2
    data["method"] = method
    data["names"] = names

    lengths = {
        "selection1": len(selection1),
        "selection2": len(selection2),
        "method": len(method),
        "names": len(names),
    }
    unique_lengths = set(lengths.values())
    if len(unique_lengths) != 1:
        parts = [f"{k}={v}" for k, v in lengths.items()]
        raise ValidationError(
            f"selection1, selection2, method, and names must all have the same "
            f"number of entries. Got: {', '.join(parts)}."
        )

    # Frame range
    begin, end, step = validate_frame_range(data)

    # Boolean flags
    add_transformations: bool = validate_bool_flag(
        data.get("add_transformations", False),
        "add_transformations",
    )
    data["add_transformations"] = add_transformations

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "distances.log"),
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

    if log_path == topology:
        raise ValidationError("Log file must not be the same as the topology file.")
    if log_path == trajectory:
        raise ValidationError("Log file must not be the same as the trajectory file.")
    if log_path == output:
        raise ValidationError("Log file must not be the same as the output file.")

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the trajectory distance analysis subcommand.

    This function reads molecular dynamics trajectory and topology files, calculates
    pairwise distances between specified atomic selections across the trajectory frames
    for given distance definitions, and writes the calculated distances along with
    statistics into an output file.

    The function also validates input selections, processes trajectory frames, applies
    optional transformations, and logs relevant steps and results. Distance calculations
    are performed using the specified method for each distance definition.

    :param args: Parsed command-line arguments containing input/output file paths,
                 logging levels, distance definitions, frame range, and other options.
    :type args: argparse.Namespace
    :return: None
    """
    import time
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
    from pharmacon.utils.mda import create_universe
    from pharmacon.fileio import PTAFile
    from pharmacon.utils.identifiers import generate_mda_blueprint, create_mda_artifact_token
    from pharmacon.analyzer.distances import calculate_frame_distances

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
    header("Distance Analysis")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Distance definitions
    subheader("Distance Definitions")
    for name, s1, s2, m in zip(args.names, args.selection1, args.selection2, args.method):
        log.info("  %-14s : %s  <-->  %s  [%s]", name, s1, s2, m)

    # Frame range
    subheader("Frame Range")
    log.info("Begin          : %d", args.begin)
    log.info("End            : %s", args.end if args.end is not None else "last")
    log.info("Step           : %d", args.step)

    if args.add_transformations:
        log.info("PBC transformations enabled")

    # Run
    subheader("Running Analysis")
    log.debug("Setting up MDAnalysis Universe ...")

    begin_time = time.time()
    log.debug("Creating the pharmacon workspace...")
    ws = Workspace(is_tmp_dir_needed=False)
    log.debug("Pharmacon workspace created successfully.")

    log.info("Reading input topology and trajectory files...")

    if args.add_transformations:
        wrapping: bool = True
        unwrapping: bool = True
    else:
        wrapping: bool = False
        unwrapping: bool = False

    u, elements = create_universe(topology=args.topology, trajectory=args.trajectory,
                                  add_elements=True, wrap=wrapping, unwrap=unwrapping,
                                  compound="atoms")

    n_frames = int(u.trajectory.n_frames)

    begin = args.begin
    if args.end is None:
        end = n_frames
    else:
        end = int(args.end)

    if end > n_frames:
        log.error("End frame index %d is greater than the total number of frames %d.", end, n_frames)
        log.warning("Setting end frame index to the last frame.")
        end = n_frames

    step = int(args.step)

    # Validate selections
    for name, s1, s2 in zip(args.names, args.selection1, args.selection2):
        for label, sel in [("selection1", s1), ("selection2", s2)]:
            log.trace("Validating %s for '%s': '%s'", label, name, sel)
            try:
                ag = u.select_atoms(sel)
            except Exception as exc:
                log.critical("Invalid MDAnalysis selection for %s '%s': '%s' — %s", label, name, sel, exc)
                raise RuntimeError(f"Invalid {label} selection for '{name}': '{sel}' — {exc}") from exc
            if ag.n_atoms == 0:
                log.critical("Selection for %s '%s' returned 0 atoms — selection='%s'", label, name, sel)
                raise RuntimeError(f"{label} selection for '{name}' has 0 atoms - selection='{sel}'")
            log.info("  %s '%s' — %d atoms selected.", label, name, ag.n_atoms)

    log.info("Calculating distances...")
    calc_start_time = time.time()

    # Use updating=True so coordinate-dependent selections are re-evaluated each frame
    atom_groups1 = [u.select_atoms(sel, updating=True) for sel in args.selection1]
    atom_groups2 = [u.select_atoms(sel, updating=True) for sel in args.selection2]

    with PTAFile(args.output, overwrite=args.overwrite,
                 command="Trajectory Analysis", subcommand="distances") as pta:
        pta.add_file_metadata(metadata={"description": "Trajectory Distance Analysis",
                                        "topology_file": str(args.topology),
                                        "trajectory_file": str(args.trajectory),
                                        "selection1": str(args.selection1),
                                        "selection2": str(args.selection2),
                                        "method": str(args.method),
                                        "names": str(args.names),
                                        "begin": str(begin), "end": str(end), "step": str(step),
                                        "is_merged": "False",
                                        "total_frames": str(n_frames),
                                        })
        pta.create_group("distances")

        for ts in u.trajectory[begin:end:step]:
            box = ts.dimensions
            frame_distances = calculate_frame_distances(atom_groups1=atom_groups1,
                                                        atom_groups2=atom_groups2,
                                                        methods=args.method,
                                                        labels=args.names,
                                                        box=box)
            pta.write_frame_distances(frame_index=ts.frame,
                                      group_name="distances",
                                      distances=frame_distances,
                                      overwrite=True,
                                      time_ps=ts.time)
            log.trace("Processed frame #%d", ts.frame)

        calc_end_time = time.time()
        calc_time = calc_end_time - calc_start_time
        log.info("Calculations completed. Total time taken: %.2f seconds.", calc_time)

        log.info("Building distance statistics...")
        pta.build_distance_statistics(group_name="distances",
                                      begin=begin, end=end, step=step)
        log.debug("Distance statistics built successfully.")

        log.debug("Finalizing PTA file data")
        blueprint = str(generate_mda_blueprint(u=u,
                                               command="Trajectory Analysis",
                                               subcommand="distances",
                                               selection1=str(args.selection1),
                                               selection2=str(args.selection2),
                                               method=str(args.method),
                                               names=str(args.names),
                                               begin=begin, end=end, step=step,
                                               ))
        meta = {
            "completed": "True",
            "artifact_status": "SUCCESS",
            "artifact_status_code": 0,
            "blueprint": blueprint,
        }
        token = create_mda_artifact_token(blueprint=blueprint,
                                          secret="trajectory_analysis",
                                          namespace="pharmacon",)
        pta.add_file_metadata({**meta, "artifact_token": token, "artifact_token_version": "1"})
        pta.add_group_metadata(group_name="distances", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Distance analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Trajectory Distance Analysis completed successfully. Time taken: %.2f seconds.", total_time)
