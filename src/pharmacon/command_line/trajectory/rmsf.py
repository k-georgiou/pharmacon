"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

trajectory rmsf — Calculate per-atom RMSF over a trajectory.
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
    validate_string_list,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "rmsf"
SUMMARY = "Calculate per-atom RMSF over a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory rmsf -p topol.tpr -x traj.xtc -o rmsf.pta -sel "protein and name CA" "backbone" -f "protein and name CA" -n calpha bb

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Constructs and configures an argparse subparser for a command-line tool.

    This function defines and adds a subparser for a specific subcommand, allowing
    customization of input, output, selection options, time range, and logging
    preferences. It facilitates the parsing of arguments used in the Pharmacon
    Trajectory Analysis process. Each parameter and option is categorized into
    distinct groups for better organization and usability.

    :param subparsers: Subparsers action object from argparse where this subparser
        will be added.
    :param parents: A list of argparse.ArgumentParser instances that serve as parent
        parsers for this subparser. Can be None if no parents are used.
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
        default="rmsf.pta",
        help="Output Pharmacon Trajectory Analysis file (default=rmsf.pta). [OPTIONAL]",
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
        help="MDAnalysis selection string(s) for per-atom RMSF calculation. [REQUIRED]",
    )

    sel.add_argument(
        "-n", "--names",
        metavar="NAME",
        default=None,
        required=True,
        nargs="+",
        help="Name(s) to assign to each selection. [REQUIRED]",
    )

    sel.add_argument(
        "-f", "--fitting-group",
        metavar="SEL",
        default=None,
        required=True,
        help=("MDAnalysis selection string for the alignment group. The trajectory "
              "is aligned in-memory on this selection before RMSF is computed. "
              "[REQUIRED]"),
    )

    sel.add_argument(
        "-r", "--reference-frame",
        metavar="FRAME",
        type=int,
        default=None,
        required=False,
        help=("Frame index used as the alignment reference. Omit to align to "
              "the average structure (two-pass alignment, recommended). "
              "[OPTIONAL]"),
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

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="rmsf.log",
                     help="File to log to (default: rmsf.log). [OPTIONAL]", )
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
    Validates the input arguments and updates the `argparse.Namespace` object in-place
    with processed and validated data. The function ensures that various files, flags,
    and configuration parameters meet the defined expectations and conform to specific
    validation rules, and raises validation errors if checks fail.

    :param args: The arguments namespace containing user-provided inputs for
        trajectory and topology analysis.
    :type args: argparse.Namespace
    :raises ValidationError: If provided arguments are invalid or contradictory.
    """
    data: dict[str, object] = vars(args).copy()

    # Input files
    topology: Path = validate_existing_input_file(
        data.get("topology"), "topology",
        supported_trajectory_analysis_topology_formats,
    )
    trajectory: Path = validate_existing_input_file(
        data.get("trajectory"), "trajectory",
        supported_trajectory_analysis_trajectory_formats,
    )
    if topology == trajectory:
        raise ValidationError("Topology and trajectory must be different files.")
    data["topology"] = topology
    data["trajectory"] = trajectory

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Output file
    output: Path = validate_output_file(
        data.get("output", "rmsf.pta"), "output",
        overwrite=overwrite, allowed_suffixes=[".pta"],
    )
    data["output"] = output
    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selections and names
    selections: list[str] = validate_string_list(data.get("selections"), "selections")
    names: list[str] = validate_string_list(data.get("names"), "names")
    if len(selections) != len(names):
        raise ValidationError(
            f"Number of selections ({len(selections)}) must match "
            f"the number of names ({len(names)})."
        )
    data["selections"] = selections
    data["names"] = names

    # Fitting group
    fitting_group: str | None = normalize_selection(
        data.get("fitting_group"), "fitting_group", required=True, default=None,
    )
    assert fitting_group is not None
    data["fitting_group"] = fitting_group

    # Reference frame: None = align to average structure (two-pass).
    raw_ref = data.get("reference_frame", None)
    if raw_ref is None:
        data["reference_frame"] = None
    else:
        data["reference_frame"] = validate_non_negative_int(raw_ref, "reference_frame")

    # Frame range
    begin, end, step = validate_frame_range(data)

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "rmsf.log"), "log",
        overwrite=True, allowed_suffixes={".log"},
    )
    file_logging_level: str = validate_logging_level(
        data.get("file_logging_level", "TRACE"), "file_logging_level",
    )
    terminal_logging_level: str = validate_logging_level(
        data.get("terminal_logging_level", "INFO"), "terminal_logging_level",
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
    Executes per-atom RMSF (Root Mean Square Fluctuation) analysis over a trajectory.

    Steps:
    - Configure logging.
    - Load topology and trajectory via MDAnalysis.
    - Validate the alignment selection and each per-selection AtomGroup.
    - Align the trajectory in-memory: first against the initial frame, then
      against the average structure computed on the alignment selection
      (canonical MDAnalysis recipe for RMSF).
    - For each selection, compute per-atom RMSF via `rms.RMSF`.
    - Write per-selection RMSF data + aggregate statistics to a PTA file.

    :param args: Parsed command-line arguments (`topology`, `trajectory`,
        `fitting_group`, `selections`, `names`, frame range, logging).
    :type args: argparse.Namespace
    :return: None
    """
    import time
    import numpy as np
    from MDAnalysis.analysis import rms, align
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
    from pharmacon.utils.mda import create_universe
    from pharmacon.fileio import PTAFile
    from pharmacon.utils.identifiers import generate_mda_blueprint, create_mda_artifact_token

    # Configure logging for this subcommand
    setup_logger(
        terminal=True,
        terminal_level=args.terminal_logging_level,
        file=True,
        file_level=args.file_logging_level,
        log_file=args.log,
    )
    log = get_logger(__name__)

    header("RMSF Analysis")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    subheader("Selections")
    for name, sel in zip(args.names, args.selections):
        log.info("  %-14s : %s", name, sel)

    subheader("RMSF Parameters")
    log.info("Fitting group  : %s", args.fitting_group)
    if args.reference_frame is None:
        log.info("Alignment      : two-pass (initial frame, then average structure)")
    else:
        log.info("Alignment      : single-pass to reference frame %d", args.reference_frame)

    subheader("Frame Range")
    log.info("Begin          : %d", args.begin)
    log.info("End            : %s", args.end if args.end is not None else "last")
    log.info("Step           : %d", args.step)

    subheader("Running Analysis")
    log.debug("Setting up MDAnalysis Universe ...")

    begin_time = time.time()
    log.debug("Creating the pharmacon workspace...")
    Workspace(is_tmp_dir_needed=False)
    log.debug("Pharmacon workspace created successfully.")

    log.info("Reading input topology and trajectory files...")
    u, elements = create_universe(topology=args.topology,
                                  trajectory=args.trajectory,
                                  add_elements=True,
                                  wrap=False,
                                  unwrap=False,
                                  compound="atoms")

    n_frames = int(u.trajectory.n_frames)
    begin = int(args.begin)
    end = n_frames if args.end is None else int(args.end)

    if end > n_frames:
        log.error("End frame index %d is greater than the total number of frames %d.", end, n_frames)
        log.warning("Setting end frame index to the last frame.")
        end = n_frames

    step = int(args.step)

    ref = args.reference_frame
    if ref is not None and ref >= n_frames:
        raise RuntimeError(
            f"reference-frame ({ref}) must be < total frames ({n_frames})"
        )

    log.trace("Selecting atoms for fitting group...")
    try:
        fitting_atoms = u.select_atoms(args.fitting_group)
    except Exception as exc:
        log.critical("Invalid MDAnalysis selection for fitting group: '%s' — %s", args.fitting_group, exc)
        raise RuntimeError(f"Invalid fitting group selection: '{args.fitting_group}' — {exc}") from exc
    if fitting_atoms.n_atoms == 0:
        log.critical("Fitting group selection returned 0 atoms — selection='%s'", args.fitting_group)
        raise RuntimeError(f"Fitting group selection has 0 atoms - selection='{args.fitting_group}'")
    log.info("Selected %d atoms for fitting group.", fitting_atoms.n_atoms)

    for name, sel in zip(args.names, args.selections):
        log.trace("Selecting atoms for calculation group '%s'...", name)
        try:
            atoms = u.select_atoms(sel)
        except Exception as exc:
            log.critical("Invalid MDAnalysis selection for group '%s': '%s' — %s", name, sel, exc)
            raise RuntimeError(f"Invalid selection for group '{name}': '{sel}' — {exc}") from exc
        if atoms.n_atoms == 0:
            log.critical("Selection for group '%s' returned 0 atoms — selection='%s'", name, sel)
            raise RuntimeError(f"Group selection has 0 atoms - selection='{sel}'")
        log.info("Selected %d atoms for calculation group '%s'.", atoms.n_atoms, name)

    # rms.RMSF performs no fitting itself, so the trajectory must be
    # pre-aligned. Two-pass (None): align to frame 0, compute the average
    # structure, then re-align to the average — tighter mean position.
    # Single-pass (int): align all frames to the user-chosen reference frame.
    if ref is None:
        log.info("Aligning trajectory to initial frame (in-memory)...")
        align.AlignTraj(u, u, select=args.fitting_group,
                        ref_frame=0, in_memory=True).run()

        log.info("Computing average structure on the alignment selection...")
        avg = align.AverageStructure(u, u, select=args.fitting_group,
                                     ref_frame=0).run()

        log.info("Re-aligning trajectory to the average structure...")
        align.AlignTraj(u, avg.results.universe,
                        select=args.fitting_group, in_memory=True).run()
    else:
        log.info("Aligning trajectory to reference frame %d (in-memory)...", ref)
        align.AlignTraj(u, u, select=args.fitting_group,
                        ref_frame=ref, in_memory=True).run()
    log.info("Alignment complete.")

    log.info("Calculating RMSF...")
    calc_start_time = time.time()
    selections_data: List[dict] = []
    for name, sel in zip(args.names, args.selections):
        atoms = u.select_atoms(sel)
        log.debug("Running RMSF on group '%s' (%d atoms)...", name, atoms.n_atoms)
        r = rms.RMSF(atoms).run(start=begin, stop=end, step=step)
        selections_data.append({
            "label": name,
            "selection_string": sel,
            "atom_indices": np.asarray(atoms.indices),
            "resids":       np.asarray(atoms.resids),
            "resnames":     np.asarray(atoms.resnames),
            "atom_names":   np.asarray(atoms.names),
            "rmsf":         np.asarray(r.results.rmsf),
        })
    calc_end_time = time.time()
    log.info("RMSF calculation complete. Total time taken: %.2f seconds.",
             calc_end_time - calc_start_time)

    with PTAFile(args.output, overwrite=args.overwrite,
                 command="Trajectory Analysis", subcommand="rmsf") as pta:
        pta.add_file_metadata(metadata={"description": "Trajectory RMSF Analysis",
                                        "topology_file": str(args.topology),
                                        "trajectory_file": str(args.trajectory),
                                        "fitting_group": str(args.fitting_group),
                                        "selections": str(args.selections),
                                        "names": str(args.names),
                                        "reference_frame": "average" if ref is None else str(ref),
                                        "begin": str(begin), "end": str(end), "step": str(step),
                                        "is_merged": "False",
                                        "total_frames": str(n_frames),
                                        })
        pta.create_group("rmsf")

        log.info("Writing RMSF data...")
        pta.write_rmsf_data(
            selections=selections_data,
            group_name="rmsf",
            frame_begin=begin,
            frame_end=end,
            frame_step=step,
            fitting_group=args.fitting_group,
            overwrite=args.overwrite,
        )
        log.debug("RMSF data written successfully.")

        log.info("Writing RMSF statistics...")
        pta.write_rmsf_statistics(group_name="rmsf")
        log.debug("RMSF statistics written successfully.")

        log.debug("Finalizing PTA file data")
        blueprint = str(generate_mda_blueprint(u=u,
                                               command="Trajectory Analysis",
                                               subcommand="rmsf",
                                               fitting_group=args.fitting_group,
                                               selections=str(args.selections),
                                               names=str(args.names),
                                               reference_frame="average" if ref is None else ref,
                                               begin=begin,
                                               end=end,
                                               step=step,
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
        pta.add_file_metadata({
            **meta,
            "artifact_token": token,
            "artifact_token_version": "1",
        })

        pta.add_group_metadata(group_name="rmsf", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    header("Done")
    log.info("RMSF analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Trajectory RMSF Analysis completed successfully. Time taken: %.2f seconds.", total_time)
