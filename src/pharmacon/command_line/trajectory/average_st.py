"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory average-st — Compute the average structure from a trajectory.
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
    normalize_selection,
    validate_bool_flag,
    validate_existing_input_file,
    validate_frame_range,
    validate_logging_level,
    validate_non_negative_int,
    validate_output_file,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "average-st"
SUMMARY = "Compute the average structure from a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory average-st -p topol.tpr -x traj.xtc -o average.pdb -sel "protein and name CA" -r 0

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds a parser for handling specific command-line arguments to configure the subcommand.

    This function creates an `argparse.ArgumentParser` instance with pre-defined
    argument groups and accepted inputs for configuration. The parser supports input,
    output, selection, time range, periodic boundary condition (PBC), and logging
    options. Each group enables users to fine-tune the behavior of the program by
    specifying detailed parameters and flags. The function also integrates optional
    parent parsers for added flexibility and reusability.

    :param subparsers: Sub-parser object used to add a specific subcommand to the parent
        parser.
    :type subparsers: argparse._SubParsersAction
    :param parents: List of parent parsers to inherit the arguments from. Can be `None`
        if no parent parsers are required.
    :type parents: List[argparse.ArgumentParser] | None
    :return: The configured `ArgumentParser` instance for the subcommand.
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
        default="average_st.pdb",
        help="Output structure file (.pdb, .crd, .gro) (default=average_st.pdb). [OPTIONAL]",
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
        "-sel", "--selection",
        metavar="SEL",
        default=None,
        required=True,
        help="MDAnalysis selection string for atoms to average. Must not be a dynamic MDAnalysis selection string. [REQUIRED]",
    )

    sel.add_argument(
        "-ws", "--write-selection",
        metavar="SEL",
        default=None,
        required=False,
        help="MDAnalysis selection string for atoms to write to the output file. "
             "If not provided, defaults to --selection. Use this to e.g. average "
             "on 'protein and name CA' but write all 'protein' atoms. [OPTIONAL]",
    )

    sel.add_argument(
        "-r", "--reference-frame",
        metavar="FRAME",
        type=int,
        default=0,
        required=False,
        help="Reference frame index for alignment (default: 0). [OPTIONAL]",
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
                     default="average_st.log",
                     help="File to log to (default: average_st.log). [OPTIONAL]", )
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
    Validates input arguments for a molecular dynamics analysis task, ensuring consistency,
    type correctness, and that essential requirements are met for further processing.

    The function validates and processes input arguments provided through an
    `argparse.Namespace`, performing the following actions:

    - Verifies the existence and format of required input files for topology and trajectory.
    - Ensures specific constraints, such as topology and trajectory being different files.
    - Validates boolean flags and their consistency with user-provided values.
    - Checks and processes output file settings with conditions on overwriting permissions
      and disallowed output file conflicts with input files.
    - Validates and normalizes the selection and write selection strings for molecular
      object selections.
    - Manages reference frame and frame range settings, ensuring values are acceptable.
    - Sets up logging attributes, including log file path, file logging level, and terminal
      logging level, while avoiding conflicts with other files.
    - Writes all validated and normalized argument values back into the parsed namespace
      to overwrite the original inputs in-place.

    :param args: The arguments namespace containing values parsed by the argparse module.
    :type args: argparse.Namespace
    :return: None
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
        data.get("output", "average_st.pdb"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=[".pdb", ".crd", ".gro"],
    )
    data["output"] = output

    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selection string
    selection: str | None = normalize_selection(
        data.get("selection"),
        "selection",
        required=True,
        default=None,
    )
    assert selection is not None
    data["selection"] = selection

    # Write selection (defaults to --selection if not provided)
    write_selection: str | None = normalize_selection(
        data.get("write_selection"),
        "write_selection",
        required=False,
        default=None,
    )
    if write_selection is None:
        write_selection = selection
    data["write_selection"] = write_selection

    # Reference frame
    reference_frame: int = validate_non_negative_int(
        data.get("reference_frame", 0),
        "reference_frame",
    )
    data["reference_frame"] = reference_frame

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
        data.get("log", "average_st.log"),
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
    Executes the average structure analysis based on the provided input arguments.
    The analysis involves reading in a molecular dynamics simulation trajectory,
    validating selections, computing an average structure, and identifying the most
    representative frame based on RMSD to the average.

    This function operates within the Pharmacon framework, utilizing MDAnalysis for
    trajectory processing and custom utilities for workspace management and logging.

    :param args: Command-line arguments encapsulated in an argparse.Namespace object.
        Expected parameters include input topology and trajectory file paths, selection
        strings, reference frame, output file paths, logging configurations, and options
        for simulation transformations like periodic boundary conditions (PBC) wrapping/unwrapping.
    :type args: argparse.Namespace

    :return: None. Outputs include logs of the analysis process, a saved file of
        the most representative structure, and additional computed data like RMSD if
        required.
    """
    import time
    import numpy as np
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
    from pharmacon.utils.mda import create_universe
    from pharmacon.analyzer.average_structure import avg_st_process_trajectory, extract_trajectory_frame

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
    header("Average Structure")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selection
    subheader("Selection")
    log.info("Selection      : %s", args.selection)
    log.info("Write selection: %s", args.write_selection)
    log.info("Reference frame: %d", args.reference_frame)

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
        end = n_frames - 1
    else:
        end = int(args.end)

    if end >= n_frames:
        log.error("End frame index %d is >= the total number of frames %d.", end, n_frames)
        log.warning("Setting end frame index to the last frame.")
        end = n_frames - 1

    step = int(args.step)
    ref = args.reference_frame

    if ref >= n_frames:
        raise RuntimeError(
            f"reference-frame ({ref}) must be < total frames ({n_frames})"
        )

    # Validate selection
    log.trace("Selecting atoms for group...")
    try:
        selection = u.select_atoms(args.selection)
    except Exception as exc:
        log.critical("Invalid MDAnalysis selection: '%s' — %s", args.selection, exc)
        raise RuntimeError(f"Invalid selection: '{args.selection}' — {exc}") from exc
    if selection.n_atoms == 0:
        log.critical("Selection returned 0 atoms — selection='%s'", args.selection)
        raise RuntimeError(f"Selection has 0 atoms - selection='{args.selection}'")
    log.info("Selected %d atoms for averaging.", selection.n_atoms)

    # Validate write selection
    log.trace("Selecting atoms for output writing...")
    try:
        write_atoms = u.select_atoms(args.write_selection)
    except Exception as exc:
        log.critical("Invalid MDAnalysis write selection: '%s' — %s", args.write_selection, exc)
        raise RuntimeError(f"Invalid write selection: '{args.write_selection}' — {exc}") from exc
    if write_atoms.n_atoms == 0:
        log.critical("Write selection returned 0 atoms — selection='%s'", args.write_selection)
        raise RuntimeError(f"Write selection has 0 atoms - selection='{args.write_selection}'")
    log.info("Selected %d atoms for output writing.", write_atoms.n_atoms)

    # Compute average structure
    log.info("Computing average structure...")
    calc_start_time = time.time()

    avg_positions, rmsd_to_avg, rmsd_to_ref = avg_st_process_trajectory(
        u=u,
        selection=selection,
        start=begin,
        stop=end,
        step=step,
        reference_frame=ref,
        memory_efficient=True,
    )

    calc_end_time = time.time()
    calc_time = calc_end_time - calc_start_time
    log.info("Average structure computation completed. Time taken: %.2f seconds.", calc_time)

    # Find most representative frame
    closest_rel = int(np.argmin(rmsd_to_avg))
    closest_frame = begin + closest_rel * step
    log.info("Most representative frame: %d (RMSD to average: %.4f)", closest_frame, rmsd_to_avg[closest_rel])

    # Save representative structure
    output_format: str = str(args.output.suffix[1:]).lower()
    log.info("Saving representative structure (format: %s)...", output_format.upper())

    extract_trajectory_frame(
        frame_idx=closest_frame,
        u=u,
        ag=write_atoms,
        output_file_path=args.output,
        output_format=output_format,
        title="Most representative structure saved by Pharmacon",
        selection=args.write_selection,
    )
    log.info("File saved at %s", args.output)

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Average structure computation complete.")
    log.info("Results written to: %s", args.output)
    log.info("Trajectory Average Structure Analysis completed successfully. Time taken: %.2f seconds.", total_time)
