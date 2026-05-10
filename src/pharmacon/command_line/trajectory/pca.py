"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory pca — Principal Component Analysis on a trajectory.
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
    validate_positive_int,
    validate_string_list,
)



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "pca"
SUMMARY = "Principal Component Analysis on a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory pca -p topol.tpr -x traj.xtc -o pca.pta -sel "protein and name CA" -c 5

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures a command-line argument parser for the specified subcommand.

    The parser supports various argument groups, organized by functionality such as
    input options, output options, PCA options, and more. It enables users to
    configure the PCA analysis in a highly customizable manner directly through the
    command-line interface.

    :param subparsers: A subparsers action to which this parser will be added.
    :param parents: A list of parent argument parsers to inherit options from. Defaults to None.
    :return: A configured argparse.ArgumentParser instance for the subcommand.
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
        default="pca.pta",
        help="Output Pharmacon Trajectory Analysis file (default=pca.pta). [OPTIONAL]",
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
        nargs="+",
        help="MDAnalysis selection string(s) for PCA. [REQUIRED]",
    )

    sel.add_argument(
        "-n", "--names",
        metavar="NAME",
        default=None,
        required=True,
        nargs="+",
        help="Name(s) to assign to each selection. [REQUIRED]",
    )

    # PCA Options
    pca_grp = parser.add_argument_group("PCA Options")
    pca_grp.add_argument(
        "-c", "--components",
        metavar="N",
        type=int,
        default=3,
        required=False,
        help="Number of principal components to compute (default: 3). [OPTIONAL]",
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

    # Parallel Options
    par = parser.add_argument_group("Parallel Options")
    par.add_argument("--parallel",
                     required=False,
                     action="store_true",
                     default=False,
                     help="Run PCA for each selection in parallel using separate universe copies (default: False). [OPTIONAL]",)

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="pca.log",
                     help="File to log to (default: pca.log). [OPTIONAL]", )
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
    Validates the input arguments passed as an argparse.Namespace object. This function ensures that
    the required file paths, flags, ranges, and other parameters are valid and conform to the
    expected data types and constraints. If validation fails, an exception is raised.

    :param args: Namespace object containing all the program's arguments to be validated.
    :type args: argparse.Namespace

    :raises ValidationError: If any of the validation checks fail due to incorrectly provided arguments,
        missing required files, conflicting parameter values, or invalid data types.

    :return: None. Mutates the provided args object in place with validated and transformed values.
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
        data.get("output", "pca.pta"), "output",
        overwrite=overwrite, allowed_suffixes=[".pta"],
    )
    data["output"] = output
    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selection and names
    selection: list[str] = validate_string_list(data.get("selection"), "selection")
    names: list[str] = validate_string_list(data.get("names"), "names")
    if len(selection) != len(names):
        raise ValidationError(
            f"Number of selections ({len(selection)}) must match "
            f"the number of names ({len(names)})."
        )
    data["selection"] = selection
    data["names"] = names

    # Components
    components: int = validate_positive_int(data.get("components", 3), "components")
    data["components"] = components

    # Frame range
    begin, end, step = validate_frame_range(data)

    # Boolean flags
    add_transformations: bool = validate_bool_flag(
        data.get("add_transformations", False), "add_transformations",
    )
    parallel: bool = validate_bool_flag(
        data.get("parallel", False), "parallel",
    )
    data["add_transformations"] = add_transformations
    data["parallel"] = parallel

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "pca.log"), "log",
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


def _run_pca_for_selection(*,
                           name: str,
                           sel_str: str,
                           topology: Path,
                           trajectory: Path,
                           n_components: int,
                           begin: int,
                           end: int,
                           step: int,
                           wrapping: bool,
                           unwrapping: bool,
                           log,
                           log_raw) -> dict:
    """
    Executes Principal Component Analysis (PCA) on a specified atomic selection from
    a molecular dynamics simulation. This process includes universe creation,
    selection verification, PCA computation, and data transformation.

    Raises exceptions for invalid selections or selections returning zero atoms.
    Logs warnings when selections contain hydrogens, as they may introduce noise
    or yield less meaningful results.

    :param name: Identifier for the PCA process
    :param sel_str: Atom selection string (MDAnalysis compatible)
    :param topology: Path to the topology file
    :param trajectory: Path to the trajectory file
    :param n_components: Number of principal components to compute
    :param begin: Starting frame for PCA calculation
    :param end: Ending frame for PCA calculation
    :param step: Frame step size for PCA calculation
    :param wrapping: Boolean indicating whether to process with periodic wrapping
    :param unwrapping: Boolean indicating whether to process with unwrapping for trajectory
    :param log: Logger object for informational messages
    :param log_raw: Logger object for frame-level debug information
    :return: Dictionary with keys:
        - "name": The provided identifier for the PCA process
        - "sel_str": The atom selection string
        - "pc": PCA object from MDAnalysis
        - "transformed": Transformed PCA data array
        - "universe": MDAnalysis Universe object
    """
    import time
    import numpy as np
    from MDAnalysis.analysis import pca
    from pharmacon.utils.mda import create_universe
    from pharmacon.analyzer.frame_logger import log_every_frame

    log.info("Running PCA for selection '%s': '%s'", name, sel_str)

    u, _ = create_universe(topology=topology, trajectory=trajectory,
                           add_elements=True, wrap=wrapping, unwrap=unwrapping,
                           compound="atoms")

    try:
        selection = u.select_atoms(sel_str)
    except Exception as exc:
        log.critical("Invalid MDAnalysis selection for '%s': '%s' — %s", name, sel_str, exc)
        raise RuntimeError(f"Invalid selection for '{name}': '{sel_str}' — {exc}") from exc
    if selection.n_atoms == 0:
        log.critical("Selection for '%s' returned 0 atoms — selection='%s'", name, sel_str)
        raise RuntimeError(f"Selection has 0 atoms for '{name}' - selection='{sel_str}'")
    log.info("Selected %d atoms for '%s'.", selection.n_atoms, name)

    # Hydrogen warning
    atom_names = np.asarray(selection.atoms.names, dtype=str)
    is_h = np.char.startswith(atom_names, "H")
    n_total = selection.n_atoms
    n_h = int(is_h.sum())
    n_heavy = n_total - n_h

    if n_h == n_total:
        log.warning(
            "PCA selection '%s' contains only hydrogen atoms. "
            "PCA on hydrogens is generally not meaningful due to high-frequency motions.",
            name,
        )
    elif n_h > 0 and n_heavy > 0:
        log.warning(
            "PCA selection '%s' contains both heavy and hydrogen atoms. "
            "Including hydrogens typically adds noise. Consider heavy atoms only or Ca atoms.",
            name,
        )

    log.info("Calculating PCA for '%s'...", name)
    timer1 = time.time()
    with log_every_frame(u=u, logger=log_raw, start=begin, stop=end, step=step,
                         prefix=f"PCA [{name}] frame"):
        pc = pca.PCA(
            universe=u,
            select=sel_str,
            align=True,
            mean=None,
            n_components=n_components,
            verbose=False,
        ).run(start=begin, stop=end, step=step)
    timer2 = time.time()
    log.info("PCA calculation for '%s' complete. Time taken: %.2f seconds.", name, timer2 - timer1)

    log.info("Transforming data for '%s'...", name)
    timer1 = time.time()
    with log_every_frame(u=u, logger=log_raw, start=begin, stop=end, step=step,
                         every=1, prefix=f"Transform [{name}] frame"):
        transformed = pc.transform(
            selection,
            n_components=n_components,
            start=begin, stop=end, step=step,
        )  # (n_frames, n_components)
    timer2 = time.time()
    log.info("Transformation for '%s' complete. Time taken: %.2f seconds.", name, timer2 - timer1)

    return {
        "name": name,
        "sel_str": sel_str,
        "pc": pc,
        "transformed": transformed,
        "universe": u,
    }


def run(args: argparse.Namespace) -> None:
    """
    Runs the Principal Component Analysis (PCA) on molecular dynamics data. This function configures
    logging, sets up the analysis workspace, validates user-defined atomic selections, and performs
    the PCA computation for each specified selection. It supports both parallel and sequential execution.

    :param args: Namespace object containing the arguments required for PCA analysis.
        The `Namespace` must include the following fields:
            topology (str): Path to the topology file.
            trajectory (str): Path to the trajectory file.
            output (str): Path where output files will be saved.
            log (str): Path to the log file.
            terminal_logging_level (str): Logging level for terminal output.
            file_logging_level (str): Logging level for file-based logging.
            components (int): Number of PCA components to calculate.
            parallel (bool): Specifies whether to perform PCA in parallel.
            begin (int): Frame index to start analysis.
            end (Optional[int]): Frame index to end analysis.
            step (int): Step interval for frame analysis.
            add_transformations (bool): Toggle for periodic boundary condition transformations.
            names (list[str]): Names for the selections.
            selection (list[str]): Atom selection strings for PCA.

    :raises RuntimeError: If any atomic selection is invalid or contains zero atoms.
    :raises ValueError: If the end frame index exceeds the total frame count.
    """
    import time
    import logging
    from concurrent.futures import ThreadPoolExecutor, as_completed
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
    log_raw = logging.getLogger(__name__)

    # Header
    header("Principal Component Analysis")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selections
    subheader("Selections")
    for name, sel in zip(args.names, args.selection):
        log.info("  %-14s : %s", name, sel)

    # PCA parameters
    subheader("PCA Parameters")
    log.info("Components     : %d", args.components)
    log.info("Parallel       : %s", args.parallel)

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

    # Create a primary universe to determine n_frames
    u, elements = create_universe(topology=args.topology, trajectory=args.trajectory,
                                  add_elements=True, wrap=wrapping, unwrap=unwrapping,
                                  compound="atoms")

    n_frames = int(u.trajectory.n_frames)
    begin = int(args.begin)
    if args.end is None:
        end = n_frames
    else:
        end = int(args.end)

    if end > n_frames:
        log.error("End frame index %d is greater than the total number of frames %d.", end, n_frames)
        log.warning("Setting end frame index to the last frame.")
        end = n_frames

    step = int(args.step)

    # Validate all selections against the primary universe
    for name, sel in zip(args.names, args.selection):
        try:
            atoms = u.select_atoms(sel)
        except Exception as exc:
            log.critical("Invalid MDAnalysis selection for '%s': '%s' — %s", name, sel, exc)
            raise RuntimeError(f"Invalid selection for '{name}': '{sel}' — {exc}") from exc
        if atoms.n_atoms == 0:
            log.critical("Selection for '%s' returned 0 atoms — selection='%s'", name, sel)
            raise RuntimeError(f"Selection has 0 atoms for '{name}' - selection='{sel}'")
        log.info("Validated selection '%s': %d atoms.", name, atoms.n_atoms)

    # Common kwargs for _run_pca_for_selection
    common_kwargs = dict(
        topology=args.topology,
        trajectory=args.trajectory,
        n_components=args.components,
        begin=begin, end=end, step=step,
        wrapping=wrapping, unwrapping=unwrapping,
        log=log, log_raw=log_raw,
    )

    # Memory check for parallel mode
    if args.parallel and len(args.selection) > 1:
        import os

        n_atoms = u.atoms.n_atoms
        n_analysis_frames = len(range(begin, end, step))
        # Estimate per-universe: positions (n_atoms * 3 * 8 bytes) per frame
        # plus PCA covariance matrix, transformed array, and MDAnalysis overhead
        per_frame_bytes = n_atoms * 3 * 8  # float64 positions
        universe_base_bytes = per_frame_bytes * 10  # topology + buffers overhead
        pca_work_bytes = n_analysis_frames * n_atoms * 3 * 8  # coordinates array
        covariance_bytes = (n_atoms * 3) ** 2 * 8  # covariance matrix
        per_worker_bytes = universe_base_bytes + pca_work_bytes + covariance_bytes
        total_estimated = per_worker_bytes * len(args.selection)

        # Get available RAM — try psutil first, fall back to OS-level methods
        available_ram: int | None = None
        try:
            import psutil
            available_ram = psutil.virtual_memory().available
        except ImportError:
            try:
                # macOS / BSD
                page_size = os.sysconf("SC_PAGE_SIZE")
                avail_pages = os.sysconf("SC_AVPHYS_PAGES")
                if page_size > 0 and avail_pages > 0:
                    available_ram = page_size * avail_pages
            except (ValueError, OSError, AttributeError):
                pass
            # If os.sysconf didn't work, available_ram stays None

        if available_ram is not None and available_ram > 0:
            ram_ratio = total_estimated / available_ram

            log.debug("Memory estimate: %.1f MB per worker x %d workers = %.1f MB total",
                      per_worker_bytes / 1e6, len(args.selection), total_estimated / 1e6)
            log.debug("Available RAM: %.1f MB (usage ratio: %.0f%%)",
                      available_ram / 1e6, ram_ratio * 100)

            if ram_ratio > 0.7:
                log.warning(
                    "Parallel PCA estimated to use ~%.1f GB but only %.1f GB RAM available (%.0f%% usage). "
                    "Consider running without --parallel to avoid memory pressure.",
                    total_estimated / 1e9, available_ram / 1e9, ram_ratio * 100,
                )
            elif ram_ratio > 0.4:
                log.warning(
                    "Parallel PCA estimated to use ~%.1f GB of %.1f GB available RAM (%.0f%%). "
                    "Monitor memory usage during the run.",
                    total_estimated / 1e9, available_ram / 1e9, ram_ratio * 100,
                )
        else:
            log.debug("Memory estimate: %.1f MB per worker x %d workers = %.1f MB total (available RAM unknown)",
                      per_worker_bytes / 1e6, len(args.selection), total_estimated / 1e6)

    # Run PCA per selection: parallel or sequential
    log.info("Calculating PCA...")
    calc_start_time = time.time()
    results: list[dict] = []

    if args.parallel and len(args.selection) > 1:
        log.info("Running %d PCA selections in parallel...", len(args.selection))
        with ThreadPoolExecutor(max_workers=len(args.selection)) as executor:
            futures = {
                executor.submit(
                    _run_pca_for_selection,
                    name=name, sel_str=sel, **common_kwargs,
                ): name
                for name, sel in zip(args.names, args.selection)
            }
            for future in as_completed(futures):
                name = futures[future]
                result = future.result()
                results.append(result)
                log.info("PCA for '%s' completed.", name)
        # Re-order results to match original selection order
        order = {name: i for i, name in enumerate(args.names)}
        results.sort(key=lambda r: order[r["name"]])
    else:
        if len(args.selection) > 1:
            log.info("Running %d PCA selections sequentially...", len(args.selection))
        for name, sel in zip(args.names, args.selection):
            result = _run_pca_for_selection(name=name, sel_str=sel, **common_kwargs)
            results.append(result)

    calc_end_time = time.time()
    calc_time = calc_end_time - calc_start_time
    log.info("All PCA calculations completed. Total time taken: %.2f seconds.", calc_time)

    # Write results to PTA file
    log.info("Writing PCA results to PTA file...")
    with PTAFile(args.output, overwrite=args.overwrite,
                 command="Trajectory Analysis", subcommand="pca") as pta:
        pta.add_file_metadata(metadata={"description": "Trajectory Principal Component Analysis (PCA)",
                                        "topology_file": str(args.topology),
                                        "trajectory_file": str(args.trajectory),
                                        "selection": str(args.selection),
                                        "names": str(args.names),
                                        "n_components": str(args.components),
                                        "begin": str(begin), "end": str(end), "step": str(step),
                                        "parallel": str(args.parallel),
                                        "is_merged": "False",
                                        "total_frames": str(n_frames),
                                        })

        for res in results:
            group_name = f"pca_{res['name']}"
            pta.create_group(group_name)

            pc = res["pc"]
            transformed = res["transformed"]
            res_u = res["universe"]

            variances = pc.variance
            variance_ratios = variances / variances.sum()

            for i, frame in enumerate(range(begin, end, step)):
                pta.write_pca_data(
                    frame_index=frame,
                    pc_values=transformed[i],
                    variances=variances,
                    variance_ratios=variance_ratios,
                    group_name=group_name,
                    time_ps=res_u.trajectory[frame].time,
                    overwrite=True,
                )

            pta.add_group_metadata(group_name=group_name, metadata={"completed": "True"})
            log.debug("PCA data written for '%s'.", res["name"])

        log.debug("Finalizing PTA file data")
        blueprint = str(generate_mda_blueprint(u=u,
                                               command="Trajectory Analysis",
                                               subcommand="pca",
                                               selection=str(args.selection),
                                               names=str(args.names),
                                               n_components=args.components,
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

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("PCA analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Trajectory PCA completed successfully. Time taken: %.2f seconds.", total_time)
