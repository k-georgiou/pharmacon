"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

trajectory h-bonds — Calculate hydrogen bonds through a trajectory.
"""

from __future__ import annotations

import time
import argparse
import itertools
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
    validate_output_file,
    validate_positive_int,
)
from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
from pharmacon.utils.mda import create_universe, convert_mda_to_rdkit
from pharmacon.analyzer.detector import (detect_hydrogen_bond_donor_atoms,
                                         detect_hydrogen_bond_acceptor_atoms)
from pharmacon.fileio import PTAFile
from pharmacon.utils.identifiers import generate_mda_blueprint, create_mda_artifact_token
from pharmacon.analyzer.interactions import interactions_process_frame, deduplicate_interactions



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "h-bonds"
SUMMARY = "Calculate hydrogen bonds through a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory h-bonds -p topol.tpr -x traj.xtc -o hbonds.pta -sel "protein" --workers 4

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Creates and configures a subparser for the specified subcommand. The parser defines various
    options and argument groups for controlling input, output, selection, time range, PBC transformations,
    parallelization, and logging settings. This function returns the configured parser to be used
    within a larger argparse-based command-line interface.

    :param subparsers: The subparsers object from an argparse.ArgumentParser where the subcommand parser
                       will be added.
    :param parents: A list of parent argument parsers to inherit global arguments from, or None if
                    there are no parents.
    :return: The configured argparse.ArgumentParser for the specific subcommand.
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
        default="hbonds.pta",
        help="Output Pharmacon Trajectory Analysis file (default=hbonds.pta). [OPTIONAL]",
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
        help="MDAnalysis selection string for hydrogen bond analysis. [REQUIRED]",
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
                     help="Add transformations to MDAnalysis Universe (default: False). [OPTIONAL]. "
                          "See documentation for more details.",)

    # Parallel Options
    par = parser.add_argument_group("Parallel Options")
    par.add_argument(
        "--workers",
        metavar="N",
        type=int,
        default=1,
        required=False,
        help="Number of parallel workers for frame processing (default: 1). [OPTIONAL]",
    )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="hbonds.log",
                     help="File to log to (default: hbonds.log). [OPTIONAL]",)
    log.add_argument("-fl", "--file-logging-level",
                     required=False,
                     metavar="LEVEL",
                     default="DEBUG",
                     choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
                              "trace", "debug", "info", "warning", "error", "critical"],
                     help="Logging level for file logging (default: DEBUG). [OPTIONAL]",)
    log.add_argument("-tl", "--terminal-logging-level",
                     required=False,
                     metavar="LEVEL",
                     default="INFO",
                     choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
                              "trace", "debug", "info", "warning", "error", "critical"],
                     help="Logging level for terminal logging (default: INFO). [OPTIONAL]",)

    return parser


def validate(args: argparse.Namespace) -> None:
    """
    Validate and normalize input arguments for a molecular dynamics analysis task.

    This function ensures that all provided input arguments are valid, well-formed, and
    meet the required specifications. It validates input files, boolean flags, selections,
    frame ranges, and logging configuration. Critical checks are performed to guarantee
    that specific inputs do not reference the same files (e.g., topology, trajectory,
    output, and log files). Additionally, the function prepares the argument namespace
    with updated validated values, written back in place for downstream use.

    :param args: Namespace object containing arguments to validate.
    :type args: argparse.Namespace

    :raises ValidationError: If validation fails due to invalid or conflicting input values.
    :raises AssertionError: If critical assumptions about derived values are violated.

    :return: None. The provided `args` namespace is modified directly to include validated values.
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
        data.get("output", "hbonds.pta"), "output",
        overwrite=overwrite, allowed_suffixes=[".pta"],
    )
    data["output"] = output
    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selection
    selection: str | None = normalize_selection(
        data.get("selection"), "selection", required=True, default=None,
    )
    assert selection is not None
    data["selection"] = selection

    # Frame range
    begin, end, step = validate_frame_range(data)

    # Boolean flags
    add_transformations: bool = validate_bool_flag(
        data.get("add_transformations", False), "add_transformations",
    )
    data["add_transformations"] = add_transformations

    # Workers
    workers: int = validate_positive_int(data.get("workers", 1), "workers")
    data["workers"] = workers

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "hbonds.log"), "log",
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


def _detect_hbond_atoms(u, group_atoms, label):
    """
    Detects hydrogen bond donor and acceptor atoms from a molecular dynamics universe and specific group
    of atoms. This function makes use of converted RDKit molecule representations and mappings between
    the MDAnalysis universe and RDKit representation. The detected donor and acceptor atoms are grouped
    and returned in a dictionary format, organized by their respective categories.

    :param u: The molecular dynamics universe to analyze.
    :type u: MDAnalysis.core.universe.Universe
    :param group_atoms: A group of atoms from the MDAnalysis universe to detect hydrogen bond donors
                        and acceptors for.
    :type group_atoms: MDAnalysis.core.groups.AtomGroup
    :param label: A descriptive label for the analyzed group, for grouping and identification purposes.
                  This enhances any reporting or debugging.
    :type label: str
    :return: A dictionary categorizing detected hydrogen bond donor and acceptor atoms for the given
             atomic group. Contains keys for both donor and acceptor atoms across group 1 and group 2.
    :rtype: dict
    """
    rdkit_mol, mapping_dict = convert_mda_to_rdkit(group_atoms)

    donors = detect_hydrogen_bond_donor_atoms(
        u, rdkit_mol, mapping_dict, label=f"Hydrogen Bond Donor atoms in {label}")
    acceptors = detect_hydrogen_bond_acceptor_atoms(
        u, rdkit_mol, mapping_dict, label=f"Hydrogen Bond Acceptor atoms in {label}")

    return dict(
        hydrogen_bond_donor_atoms_group1=donors,
        hydrogen_bond_donor_atoms_group2=donors,
        hydrogen_bond_acceptor_atoms_group1=acceptors,
        hydrogen_bond_acceptor_atoms_group2=acceptors,
    )


def _extract_atom_indices(atom_groups: dict) -> dict:
    """
    Extracts atom indices from a dictionary of atom groups.

    This function processes a dictionary where each value contains an object with
    an `indices` attribute. It converts these `indices` attributes, which are
    expected to be arrays, into lists and constructs a new dictionary preserving
    the original keys mapped to the converted lists.

    :param atom_groups: A dictionary with keys of any type and values that are
        objects containing an `indices` attribute.
    :type atom_groups: dict
    :return: A dictionary where the keys are the same as `atom_groups`, and the
        values are lists originating from the `indices` attributes of the
        corresponding values from `atom_groups`.
    :rtype: dict
    """
    return {key: ag.indices.tolist() for key, ag in atom_groups.items()}


def _rebuild_atom_groups(u, atom_indices: dict) -> dict:
    """
    Rebuild atom groups by mapping provided atom indices to corresponding atoms
    in the given universe.

    :param u: The molecular universe containing atom information.
    :type u: Any
    :param atom_indices: Dictionary mapping group names to indices of atoms
        in the molecular universe.
    :type atom_indices: dict
    :return: A dictionary mapping group names to updated atom groups based
        on the provided indices.
    :rtype: dict
    """
    return {key: u.atoms[indices] for key, indices in atom_indices.items()}


def _worker_process_frames(*,
                           worker_id: int,
                           topology,
                           trajectory,
                           wrapping: bool,
                           unwrapping: bool,
                           frame_indices: list,
                           disable_flags: dict,
                           atom_indices: dict,
                           db_path,
                           log_file=None,
                           log_queue=None,
                           file_logging_level: str = "TRACE"
                           ):
    """
    Processes frames for interaction analysis in a worker thread.

    This function is intended to be executed in a worker process, enabling parallel
    processing of molecular dynamics trajectory frames. Each worker is responsible
    for processing a subset of the frames provided in `frame_indices`. The function
    also manages logging specific to the worker, handles atom group reconstruction,
    and writes processed results to a SQLite database.

    :param worker_id: The unique identifier of the worker processing frames.
    :type worker_id: int
    :param topology: The topology file required to construct the molecular simulation
        Universe.
    :param trajectory: The trajectory file containing molecular simulation data.
    :param wrapping: Flag indicating whether to apply per-frame periodic boundary
        corrections to all coordinates (wrapping).
    :type wrapping: bool
    :param unwrapping: Flag indicating whether to apply trajectory-wide periodic boundary
        corrections to all coordinates (unwrapping).
    :type unwrapping: bool
    :param frame_indices: List of indices of frames assigned to the worker for
        processing.
    :type frame_indices: list
    :param disable_flags: Dictionary of boolean flags to disable certain processing
        features (e.g., certain interaction types).
    :type disable_flags: dict
    :param atom_indices: Dictionary containing pre-computed atom group indices,
        used for reconstructing AtomGroups in the molecular simulation Universe.
    :type atom_indices: dict
    :param db_path: Path to the SQLite database where the processed interaction data
        will be stored.
    :param log_file: Optional path to a log file to enable file-based logging for this
        worker. If ``None``, no logging will occur to a file.
    :param file_logging_level: String specifying the logging level for the log file
        if `log_file` is provided. Defaults to "TRACE".
    :type file_logging_level: str

    :return: A tuple containing the worker ID and the number of frames successfully
        processed.
    :rtype: tuple
    """
    import json
    import sqlite3
    from pharmacon.logger import setup_logger, get_logger

    # Set up per-worker logger. If a log_queue was supplied by the parent,
    # records are streamed live into the parent's main log via the
    # QueueListener running there (single unified log file). Otherwise
    # fall back to a per-rank file (still safe; was the prior behaviour).
    if log_queue is not None:
        setup_logger(
            terminal=False,
            file=False,
            log_queue=log_queue,
            queue_level=file_logging_level,
            replace=False,
        )
    elif log_file is not None:
        setup_logger(
            terminal=False,
            file=True,
            file_level=file_logging_level,
            log_file=log_file,
            per_rank=True,
            mpi_rank=worker_id,
            replace=False,
        )
    log = get_logger(f"{__name__}.worker_{worker_id}")

    log.trace("Worker %d starting: %d frames assigned.", worker_id, len(frame_indices))

    u, _ = create_universe(topology=topology, trajectory=trajectory,
                           add_elements=True, wrap=wrapping, unwrap=unwrapping,
                           compound="atoms")

    # Rebuild AtomGroups on this worker's universe from pre-computed indices
    atom_groups = _rebuild_atom_groups(u, atom_indices)

    frame_kwargs = dict(
        **disable_flags,
        **atom_groups,
        u=u,
    )

    conn = sqlite3.connect(str(db_path), timeout=60)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    for fidx in frame_indices:
        ts = u.trajectory[fidx]
        frame_kwargs["box"] = ts.dimensions
        frame_contacts = interactions_process_frame(**frame_kwargs)
        flat = itertools.chain.from_iterable(frame_contacts)
        flat = deduplicate_interactions(flat)

        rows = [(fidx, json.dumps(rec)) for rec in flat]
        if rows:
            conn.executemany(
                "INSERT INTO interactions (frame_index, data) VALUES (?, ?)",
                rows,
            )
            conn.commit()

        log.trace("Worker %d processed frame #%d", worker_id, fidx)

    conn.close()
    log.trace("Worker %d finished: %d frames.", worker_id, len(frame_indices))
    return worker_id, len(frame_indices)


def _estimate_memory_and_warn(*, n_atoms, n_workers, n_analysis_frames, log):
    """
    Estimates the memory usage for a parallel operation and logs warnings or debug
    information based on the available system RAM compared to the estimated
    requirements. The function calculates memory usage per worker, the total memory
    usage, and emits logging messages if memory consumption might exceed safe
    limits or needs attention. The method attempts multiple approaches to determine
    the available RAM, handling cases where such information is unavailable.

    :param n_atoms: The number of atoms involved in the computation.
    :type n_atoms: int
    :param n_workers: The number of parallel workers to perform computations.
    :type n_workers: int
    :param n_analysis_frames: The number of frames considered in the analysis.
    :type n_analysis_frames: int
    :param log: The logging object used to record debug and warning messages.
    :type log: logging.Logger
    :return: None
    """
    import os

    per_frame_bytes = n_atoms * 3 * 8
    universe_base_bytes = per_frame_bytes * 10
    work_bytes = n_analysis_frames * n_atoms * 3 * 8
    per_worker_bytes = universe_base_bytes + work_bytes
    total_estimated = per_worker_bytes * n_workers

    available_ram = None
    try:
        import psutil
        available_ram = psutil.virtual_memory().available
    except ImportError:
        try:
            page_size = os.sysconf("SC_PAGE_SIZE")
            avail_pages = os.sysconf("SC_AVPHYS_PAGES")
            if page_size > 0 and avail_pages > 0:
                available_ram = page_size * avail_pages
        except (ValueError, OSError, AttributeError):
            pass

    if available_ram is not None and available_ram > 0:
        ram_ratio = total_estimated / available_ram

        log.debug("Memory estimate: %.1f MB per worker x %d workers = %.1f MB total",
                  per_worker_bytes / 1e6, n_workers, total_estimated / 1e6)
        log.debug("Available RAM: %.1f MB (usage ratio: %.0f%%)",
                  available_ram / 1e6, ram_ratio * 100)

        if ram_ratio > 0.7:
            log.warning(
                "Parallel mode estimated to use ~%.1f GB but only %.1f GB RAM available (%.0f%%). "
                "Consider reducing --workers to avoid memory pressure.",
                total_estimated / 1e9, available_ram / 1e9, ram_ratio * 100,
            )
        elif ram_ratio > 0.4:
            log.warning(
                "Parallel mode estimated to use ~%.1f GB of %.1f GB available RAM (%.0f%%). "
                "Monitor memory usage during the run.",
                total_estimated / 1e9, available_ram / 1e9, ram_ratio * 100,
            )
    else:
        log.debug("Memory estimate: %.1f MB per worker x %d workers = %.1f MB total (available RAM unknown)",
                  per_worker_bytes / 1e6, n_workers, total_estimated / 1e6)


def run(args: argparse.Namespace) -> None:
    """
    Executes the hydrogen bond analysis based on the provided arguments. The function configures
    logging, validates the inputs, sets up the working environment, and performs the analysis. It
    also handles frame selection, workspace management, atom selection, interaction detection, and
    result storage to the specified output.

    :param args: Namespace object containing the command-line arguments and configurations for
                 the analysis. Includes options for input files, logging levels, output, frame
                 range, selection, and computational flags.
    :type args: argparse.Namespace
    :raises RuntimeError: Raised when invalid selections are provided, no atoms are selected, no
                          frames are available to process, or critical errors occur during analysis.
    :return: None
    """
    from pharmacon.logger import setup_logger, get_logger, header, subheader

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
    header("Hydrogen Bond Analysis")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selection
    subheader("Selection")
    log.info("Selection      : %s", args.selection)

    # Frame range
    subheader("Frame Range")
    log.info("Begin          : %d", args.begin)
    log.info("End            : %s", args.end if args.end is not None else "last")
    log.info("Step           : %d", args.step)

    if args.add_transformations:
        log.info("PBC transformations enabled")

    log.info("Workers        : %d", args.workers)

    # Run
    subheader("Running Analysis")
    log.debug("Setting up MDAnalysis Universe ...")

    begin_time = time.time()
    need_tmp = args.workers > 1
    log.debug("Creating the pharmacon workspace...")
    # cleanup_on_exit=True ensures the temp directory is removed on normal
    # exit, unhandled exceptions, and KeyboardInterrupt (via atexit).
    ws = Workspace(is_tmp_dir_needed=need_tmp, cleanup_on_exit=True)
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

    # Validate selection
    log.trace("Selecting atoms for selection group...")
    try:
        group_atoms = u.select_atoms(args.selection)
    except Exception as exc:
        log.critical("Invalid MDAnalysis selection: '%s' — %s", args.selection, exc)
        raise RuntimeError(f"Invalid selection: '{args.selection}' — {exc}") from exc
    if group_atoms.n_atoms == 0:
        log.critical("Selection returned 0 atoms — selection='%s'", args.selection)
        raise RuntimeError(f"Selection has 0 atoms - selection='{args.selection}'")
    log.info("Selected %d atoms for hydrogen bond analysis.", group_atoms.n_atoms)

    # Build the full frame list
    all_frames = list(range(begin, end, step))
    n_analysis_frames = len(all_frames)

    if n_analysis_frames == 0:
        raise RuntimeError("No frames to process. Check --begin, --end, --step.")

    n_workers = min(args.workers, n_analysis_frames)

    # Disable flags: only hbonds enabled
    disable_flags = dict(
        disable_hydrophobic=True,
        disable_hbonds=False,
        disable_water_bridges=True,
        disable_pi_stacking=True,
        disable_pi_cation=True,
        disable_ionic=True,
        disable_metal_complexation=True,
        disable_halogen=True,
    )

    log.info("Calculating hydrogen bonds...")
    calc_start_time = time.time()

    if n_workers <= 1:
        # ── Single-worker: original sequential path (no SQLite overhead) ─────
        log.info("Detecting H-bond donors and acceptors...")
        atom_groups = _detect_hbond_atoms(u, group_atoms, "Selection")

        frame_kwargs = dict(**disable_flags, **atom_groups, u=u)

        with PTAFile(args.output, overwrite=args.overwrite,
                     command="Trajectory Analysis", subcommand="h-bonds") as pta:
            pta.add_file_metadata(metadata={"description": "Hydrogen Bond Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "selection": str(args.selection),
                                            "begin": str(begin), "end": str(end), "step": str(step),
                                            "is_merged": "False",
                                            "total_frames": str(n_frames),
                                            "workers": "1",
                                            })
            pta.create_group("hbonds")
            for ts in u.trajectory[begin:end:step]:
                frame_kwargs["box"] = ts.dimensions
                frame_contacts = interactions_process_frame(**frame_kwargs)
                flat = itertools.chain.from_iterable(frame_contacts)
                flat = deduplicate_interactions(flat)
                pta.write_frame_interactions(frame_index=ts.frame,
                                             group_name="hbonds",
                                             interactions=flat,
                                             overwrite=True)
                log.trace("Processed frame #%d", ts.frame)

            calc_end_time = time.time()
            calc_time = calc_end_time - calc_start_time
            log.info("Calculations completed. Total time taken: %.2f seconds.", calc_time)

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="hbonds",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="h-bonds",
                                                   selection=args.selection,
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
            pta.add_group_metadata(group_name="hbonds", metadata={"completed": "True"})

    else:
        # Multi-worker: SQLite temp storage
        import json
        import sqlite3
        import math
        import logging
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from logging.handlers import QueueListener

        log.info("Parallel mode: %d workers for %d frames.", n_workers, n_analysis_frames)

        # Memory warning
        _estimate_memory_and_warn(
            n_atoms=u.atoms.n_atoms,
            n_workers=n_workers,
            n_analysis_frames=n_analysis_frames,
            log=log,
        )

        # Detect H-bond atoms ONCE on the main universe, then extract indices
        log.info("Detecting H-bond donors and acceptors (once, on main universe)...")
        atom_groups = _detect_hbond_atoms(u, group_atoms, "Selection")
        atom_indices = _extract_atom_indices(atom_groups)

        # Create SQLite database in workspace temp directory
        db_path = ws.temp_directory / "interactions.db"
        conn = sqlite3.connect(str(db_path))
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute(
            "CREATE TABLE interactions ("
            "  frame_index INTEGER NOT NULL,"
            "  data TEXT NOT NULL"
            ")"
        )
        conn.execute("CREATE INDEX idx_frame ON interactions(frame_index)")
        conn.commit()
        conn.close()

        # Split frames into chunks for each worker
        chunk_size = math.ceil(n_analysis_frames / n_workers)
        frame_chunks = [
            all_frames[i * chunk_size : (i + 1) * chunk_size]
            for i in range(n_workers)
        ]
        # Remove any empty trailing chunks
        frame_chunks = [ch for ch in frame_chunks if ch]

        log.info("Dispatching %d workers...", len(frame_chunks))
        for i, chunk in enumerate(frame_chunks):
            log.debug("  Worker %d: frames %d–%d (%d frames)", i, chunk[0], chunk[-1], len(chunk))

        # Live-merge worker logs into the main log via a multiprocessing queue.
        log_manager = multiprocessing.Manager()
        log_queue = log_manager.Queue()
        listener_handlers = list(logging.getLogger("pharmacon").handlers)
        log_listener = QueueListener(log_queue, *listener_handlers, respect_handler_level=True)
        log_listener.start()

        try:
            with ProcessPoolExecutor(max_workers=len(frame_chunks)) as executor:
                futures = {
                    executor.submit(
                        _worker_process_frames,
                        worker_id=i,
                        topology=args.topology,
                        trajectory=args.trajectory,
                        wrapping=wrapping,
                        unwrapping=unwrapping,
                        frame_indices=chunk,
                        disable_flags=disable_flags,
                        atom_indices=atom_indices,
                        db_path=db_path,
                        log_queue=log_queue,
                        file_logging_level=args.file_logging_level,
                    ): i
                    for i, chunk in enumerate(frame_chunks)
                }
                for future in as_completed(futures):
                    wid = futures[future]
                    _, n_done = future.result()
                    log.info("Worker %d finished (%d frames processed).", wid, n_done)
        finally:
            log_listener.stop()
            log_manager.shutdown()

        calc_end_time = time.time()
        calc_time = calc_end_time - calc_start_time
        log.info("All workers completed. Total calculation time: %.2f seconds.", calc_time)

        # Read SQLite and write to PTA frame-by-frame
        log.info("Writing interactions from temp database to PTA file...")
        conn = sqlite3.connect(str(db_path))

        with PTAFile(args.output, overwrite=args.overwrite,
                     command="Trajectory Analysis", subcommand="h-bonds") as pta:
            pta.add_file_metadata(metadata={"description": "Hydrogen Bond Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "selection": str(args.selection),
                                            "begin": str(begin), "end": str(end), "step": str(step),
                                            "is_merged": "False",
                                            "total_frames": str(n_frames),
                                            "workers": str(n_workers),
                                            })
            pta.create_group("hbonds")

            for fidx in all_frames:
                cursor = conn.execute(
                    "SELECT data FROM interactions WHERE frame_index = ? ORDER BY rowid",
                    (fidx,),
                )
                interactions = [tuple(json.loads(row[0])) for row in cursor]
                pta.write_frame_interactions(frame_index=fidx,
                                             group_name="hbonds",
                                             interactions=interactions,
                                             overwrite=True)
                log.debug("Wrote frame %d (%d interactions)", fidx, len(interactions))

            conn.close()

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="hbonds",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="h-bonds",
                                                   selection=args.selection,
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
            pta.add_group_metadata(group_name="hbonds", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Hydrogen bond analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Hydrogen Bond Analysis completed successfully. Time taken: %.2f seconds.", total_time)
