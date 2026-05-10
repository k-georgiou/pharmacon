"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory angles — Calculate angles (3-atom, vector, dihedral) through a trajectory.
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


SUBCOMMAND_NAME = "angles"
SUMMARY = "Calculate angles (3-atom, vector, dihedral) through a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory angles -p topol.tpr -x traj.xtc -o angles.pta -sel "index 4 or index 5 -> x-axis" -n vector1

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures a subcommand parser for the specified subparsers action.

    This function adds specific argument groups and parameters related to input, output,
    selection, time range, periodic boundary conditions (PBC), and logging options. It
    uses these options to define the behavior and arguments required for the subcommand
    parser. The returned parser can be used to handle user inputs and execute the
    intended operations of the subcommand.

    :param subparsers: The argparse `_SubParsersAction` object where the parser
        for this subcommand should be registered.
    :param parents: An optional list of parent `argparse.ArgumentParser`s whose
        arguments should be inherited by this parser. Defaults to None.

    :return: An `argparse.ArgumentParser` configured with argument groups and
        parameters for the specific subcommand.
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
        default="angles.pta",
        help="Output Pharmacon Trajectory Analysis file (default=angles.pta). [OPTIONAL]",
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
        "-sel", "--selections",
        metavar="SEL",
        default=None,
        required=True,
        nargs="+",
        help="Selection(s) string(s) to define angles (default: None). [REQUIRED]",
    )

    sel.add_argument(
        "-n", "--names",
        metavar="NAME",
        default=None,
        required=True,
        nargs="+",
        help="Name(s) to assign to each angle (default: None). [REQUIRED]",
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
                     default="angles.log",
                     help="File to log to (default: angles.log). [OPTIONAL]", )
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
    Validates the given arguments and ensures the input configuration is consistent
    with expected parameters and constraints. This includes validation of input files,
    output files, and specific parameter relationships such as selection and names
    lists or file-specific restrictions. Any inconsistencies or violations of
    constraints will raise a `ValidationError`.

    :param args: The parsed command-line arguments namespace to validate. It contains
        user-defined parameters and file paths.
    :type args: argparse.Namespace

    :raises ValidationError: If any validation step fails, such as mismatching
        file types, constraints violations, or invalid parameter definitions.

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
        data.get("output", "angles.pta"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=[".pta"],
    )
    data["output"] = output

    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selections and names
    selections: list[str] = validate_string_list(data.get("selections"), "selections")
    names: list[str] = validate_string_list(data.get("names"), "names")
    data["selections"] = selections
    data["names"] = names

    if len(selections) != len(names):
        raise ValidationError(
            f"Number of selections ({len(selections)}) must match "
            f"the number of names ({len(names)})."
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
        data.get("log", "angles.log"),
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
    Executes the trajectory angle analysis sub-command. This method reads molecular dynamics
    topology and trajectory files, calculates specified angle metrics, and outputs the results
    to a designated file. It handles necessary configurations, frame range slicing, and unit
    transformations during analysis.

    :param args: Parsed command-line arguments containing attributes such as topology file path,
        trajectory file path, output path, frame range, and angle definitions. These arguments
        define the inputs and parameters for the angle analysis process.
    :type args: argparse.Namespace

    :raises RuntimeError: If invalid angle selections are provided in the input arguments.
    """
    import time
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
    from pharmacon.utils.mda import create_universe
    from pharmacon.fileio import PTAFile
    from pharmacon.utils.identifiers import generate_mda_blueprint, create_mda_artifact_token
    from pharmacon.analyzer.angles import create_mda_angle_selections, calculate_frame_angles

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
    header("Angle Analysis")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selections
    subheader("Angle Definitions")
    for name, sel in zip(args.names, args.selections):
        log.info("  %-14s : %s", name, sel)

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

    # Resolve angle specifications
    log.info("Resolving angle specifications...")
    angle_specs = []
    for name, sel in zip(args.names, args.selections):
        log.trace("Resolving angle spec for '%s': '%s'", name, sel)
        try:
            spec = create_mda_angle_selections(u, sel)
        except Exception as exc:
            log.critical("Invalid angle selection for '%s': '%s' — %s", name, sel, exc)
            raise RuntimeError(f"Invalid angle selection for '{name}': '{sel}' — {exc}") from exc
        angle_specs.append(spec)

        kind = spec["type"]
        if kind == "vector-angle":
            log.info("  %-14s : VECTOR–VECTOR", name)
            left_type, left_obj = spec["left"]
            right_type, right_obj = spec["right"]

            def _vec_to_str(vtype, obj):
                if vtype == "axis":
                    return "axis (%.1f, %.1f, %.1f)" % (obj[0], obj[1], obj[2])
                a, b = obj.atoms
                return "atoms (%d -> %d)" % (int(a.index), int(b.index))

            log.debug("    left : %s", _vec_to_str(left_type, left_obj))
            log.debug("    right: %s", _vec_to_str(right_type, right_obj))
        elif kind == "angle":
            atoms = spec["atoms"].atoms
            indices = " - ".join(str(int(a.index)) for a in atoms)
            log.info("  %-14s : ANGLE (3 atoms: %s)", name, indices)
        elif kind == "dihedral":
            atoms = spec["atoms"].atoms
            indices = " - ".join(str(int(a.index)) for a in atoms)
            log.info("  %-14s : DIHEDRAL (4 atoms: %s)", name, indices)
        else:
            log.info("  %-14s : %s", name, kind)

    log.info("Calculating angles...")
    calc_start_time = time.time()

    with PTAFile(args.output, overwrite=args.overwrite,
                 command="Trajectory Analysis", subcommand="angles") as pta:
        pta.add_file_metadata(metadata={"description": "Trajectory Angle Analysis",
                                        "topology_file": str(args.topology),
                                        "trajectory_file": str(args.trajectory),
                                        "selections": str(args.selections),
                                        "names": str(args.names),
                                        "begin": str(begin), "end": str(end), "step": str(step),
                                        "is_merged": "False",
                                        "total_frames": str(n_frames),
                                        })
        pta.create_group("angles")

        for ts in u.trajectory[begin:end:step]:
            box = ts.dimensions
            frame_angles = calculate_frame_angles(specs=angle_specs,
                                                  labels=args.names,
                                                  box=box)
            pta.write_frame_angles(frame_index=ts.frame,
                                   group_name="angles",
                                   angles=frame_angles,
                                   overwrite=True,
                                   time_ps=ts.time)
            log.trace("Processed frame #%d", ts.frame)

        calc_end_time = time.time()
        calc_time = calc_end_time - calc_start_time
        log.info("Calculations completed. Total time taken: %.2f seconds.", calc_time)

        log.info("Building angle statistics...")
        pta.build_angle_statistics(group_name="angles",
                                   begin=begin, end=end, step=step)
        log.debug("Angle statistics built successfully.")

        log.debug("Finalizing PTA file data")
        blueprint = str(generate_mda_blueprint(u=u,
                                               command="Trajectory Analysis",
                                               subcommand="angles",
                                               selections=str(args.selections),
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
        pta.add_group_metadata(group_name="angles", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Angle analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Trajectory Angle Analysis completed successfully. Time taken: %.2f seconds.", total_time)
