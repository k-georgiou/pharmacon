"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory pp-interactions — Calculate protein-protein interactions.
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
from pharmacon.utils.validation import (normalize_selection, validate_bool_flag, validate_existing_input_file,
                                        validate_frame_range, validate_logging_level, validate_output_file,
                                        validate_positive_int)
from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
from pharmacon.utils.mda import create_universe, convert_mda_to_rdkit
from pharmacon.analyzer.detector import (detect_hydrophobic_atoms, detect_hydrogen_bond_donor_atoms, detect_metal_atoms,
                                         detect_halogen_atoms, detect_aromatic_atoms, detect_negative_charge_atoms,
                                         detect_positive_charge_atoms, detect_hydrogen_bond_acceptor_atoms)
from pharmacon.fileio import PTAFile
from pharmacon.utils.identifiers import generate_mda_blueprint, create_mda_artifact_token
from pharmacon.analyzer.interactions import interactions_process_frame, deduplicate_interactions
from pharmacon.constants.smarts import WATER_RESIDUES_MDA_SELECTION_STR



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "pp-interactions"
SUMMARY = "Calculate protein-protein interactions through a trajectory."

_EPILOG = """\
Examples:

  pharmacon trajectory pp-interactions -p topol.tpr -x traj.xtc -o pp_interactions.pta --prt1 "chainid A" --prt2 "chainid B"

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures an argument parser for the given subcommand.

    This function creates an argparse parser with various groups of arguments
    to ensure flexibility and completeness in the command-line interface. The
    parser supports input, output, selection, time range, parallelization,
    and logging options, among others. It also allows toggling of specific
    interaction types (e.g., hydrophobic, hydrogen bonds) and enables optional
    transformations in periodic boundary conditions (PBCs). The function configures
    logging levels and file handling options for better debugging and operational
    transparency.

    :param subparsers: subparsers object to attach the new argument parser.
    :type subparsers: argparse._SubParsersAction
    :param parents: optionally include parent parsers for shared argument parsing.
    :type parents: List[argparse.ArgumentParser] | None
    :return: fully configured argument parser for the subcommand.
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
        default="pp_interactions.pta",
        help="Output Pharmacon Trajectory Analysis file (default=pp_interactions.pta). [OPTIONAL]",
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
        "-prt1", "--protein1",
        metavar="SEL",
        default=None,
        required=True,
        help="MDAnalysis selection string for protein 1 atoms. [REQUIRED]",
    )

    sel.add_argument(
        "-prt2", "--protein2",
        metavar="SEL",
        default=None,
        required=True,
        help="MDAnalysis selection string for protein 2 atoms. [REQUIRED]",
    )

    sel.add_argument(
        "-w", "--water",
        metavar="SEL",
        default=None,
        required=False,
        help="MDAnalysis selection string for water atoms (default: None). [OPTIONAL]",
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
        help="Irritate every Nth frame (default: 1). [OPTIONAL]")

    # PBC Options
    pbc = parser.add_argument_group("PBC Options")
    pbc.add_argument("-at", "--add-transformations",
                     required=False,
                     action="store_true",
                     help="Add transformations to MDAnalysis Universe (default: False). [OPTIONAL]. See documentation for more details.", )

    # Parallel Options
    par = parser.add_argument_group("Parallel Options")
    par.add_argument(
        "--workers",
        metavar="N",
        type=int,
        default=1,
        required=False,
        help="Number of parallel workers. Each worker processes a chunk of frames "
             "using its own universe copy. Results are streamed through a temporary "
             "SQLite database to avoid memory pressure (default: 1). [OPTIONAL]",
    )

    # Toggle Options
    toggle = parser.add_argument_group("Interaction Toggle Options")
    toggle.add_argument("--disable-hydrophobic",
                        required=False,
                        action="store_true",
                        help="Disable hydrophobic interactions (default: False). [OPTIONAL]",
                        )
    toggle.add_argument("--disable-hbonds",
                        required=False,
                        action="store_true",
                        help="Disable hydrogen bonding interactions (default: False). [OPTIONAL]",
                        )
    toggle.add_argument("--disable-pi-stacking",
                        required=False,
                        action="store_true",
                        help="Disable pi-stacking interactions (default: False). [OPTIONAL]", )
    toggle.add_argument("--disable-pi-cation",
                        required=False,
                        action="store_true",
                        help="Disable pi-cation interactions (default: False). [OPTIONAL]", )
    toggle.add_argument("--disable-ionic",
                        required=False,
                        action="store_true",
                        help="Disable ionic interactions (default: False). [OPTIONAL]", )
    toggle.add_argument("--disable-water-bridges",
                        required=False,
                        action="store_true",
                        help="Disable water bridging interactions (default: False). [OPTIONAL]", )
    toggle.add_argument("--disable-halogen",
                        required=False,
                        action="store_true",
                        help="Disable halogen bonding interactions (default: False). [OPTIONAL]", )
    toggle.add_argument("--disable-metal",
                        required=False,
                        action="store_true",
                        help="Disable metal ion interactions (default: False). [OPTIONAL]", )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="pp_interactions.log",
                     help="File to log to (default: pp_interactions.log). [OPTIONAL]", )
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
    Validates and processes the command-line arguments for a trajectory analysis task.
    The function ensures all input and output files are properly specified, various boolean
    flags are appropriately set, and frame ranges are validated. Additionally, interactions
    and logging configurations are verified. It also enforces required constraints, such as
    unique inputs and the presence of valid configurations.

    :param args: Namespace containing all parsed arguments and their values.
    :type args: argparse.Namespace

    :raises ValidationError:
        - If topology and trajectory are not different files.
        - If output file is the same as either topology or trajectory file.
        - If Protein 1 and Protein 2 selections are identical.
        - If all interaction types are disabled.
        - If water bridge interactions are enabled without a water selection.
        - If log file is the same as topology, trajectory, or output file.
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
        data.get("output", "pp_interactions.pta"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=[".pta"],
    )
    data["output"] = output

    if output == topology:
        raise ValidationError("Output file must not be the same as the topology file.")
    if output == trajectory:
        raise ValidationError("Output file must not be the same as the trajectory file.")

    # Selection strings
    protein1: str | None = normalize_selection(
        data.get("protein1"),
        "protein1",
        required=True,
        default=None,
    )
    protein2: str | None = normalize_selection(
        data.get("protein2"),
        "protein2",
        required=True,
        default=None,
    )
    water: str | None = normalize_selection(
        data.get("water"),
        "water",
        required=False,
        default=None,
    )

    assert protein1 is not None
    assert protein2 is not None

    if protein1 == protein2:
        raise ValidationError(
            "Protein 1 and Protein 2 selections must be different."
        )

    data["protein1"] = protein1
    data["protein2"] = protein2
    data["water"] = water

    # Frame range
    begin, end, step = validate_frame_range(data)

    # Boolean flags
    add_transformations: bool = validate_bool_flag(
        data.get("add_transformations", False),
        "add_transformations",
    )
    disable_hydrophobic: bool = validate_bool_flag(
        data.get("disable_hydrophobic", False),
        "disable_hydrophobic",
    )
    disable_hbonds: bool = validate_bool_flag(
        data.get("disable_hbonds", False),
        "disable_hbonds",
    )
    disable_pi_stacking: bool = validate_bool_flag(
        data.get("disable_pi_stacking", False),
        "disable_pi_stacking",
    )
    disable_pi_cation: bool = validate_bool_flag(
        data.get("disable_pi_cation", False),
        "disable_pi_cation",
    )
    disable_ionic: bool = validate_bool_flag(
        data.get("disable_ionic", False),
        "disable_ionic",
    )
    disable_water_bridges: bool = validate_bool_flag(
        data.get("disable_water_bridges", False),
        "disable_water_bridges",
    )
    disable_halogen: bool = validate_bool_flag(
        data.get("disable_halogen", False),
        "disable_halogen",
    )
    disable_metal: bool = validate_bool_flag(
        data.get("disable_metal", False),
        "disable_metal",
    )

    data["add_transformations"] = add_transformations
    data["disable_hydrophobic"] = disable_hydrophobic
    data["disable_hbonds"] = disable_hbonds
    data["disable_pi_stacking"] = disable_pi_stacking
    data["disable_pi_cation"] = disable_pi_cation
    data["disable_ionic"] = disable_ionic
    data["disable_water_bridges"] = disable_water_bridges
    data["disable_halogen"] = disable_halogen
    data["disable_metal"] = disable_metal

    enabled_interactions: dict[str, bool] = {
        "hydrophobic": not disable_hydrophobic,
        "hbonds": not disable_hbonds,
        "pi_stacking": not disable_pi_stacking,
        "pi_cation": not disable_pi_cation,
        "ionic": not disable_ionic,
        "water_bridges": not disable_water_bridges,
        "halogen": not disable_halogen,
        "metal": not disable_metal,
    }

    if not any(enabled_interactions.values()):
        raise ValidationError(
            "All interaction types are disabled. "
            "Enable at least one interaction calculation."
        )

    if enabled_interactions["water_bridges"] and water is None:
        raise ValidationError(
            "Water bridge interactions are enabled, but no water selection was "
            "provided. Use '--water <SEL>' or disable water bridges with "
            "'--disable-water-bridges'."
        )

    data["enabled_interactions"] = tuple(
        name for name, enabled in enabled_interactions.items() if enabled
    )

    # Workers
    workers: int = validate_positive_int(data.get("workers", 1), "workers")
    data["workers"] = workers

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "pp_interactions.log"),
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
    Executes the main workflow for analyzing protein-protein interactions in molecular
    dynamics simulations, leveraging user-provided arguments and configurations. Handles
    logging, trajectory reading, atom selection, and calculation of interactions.

    :param args: Parsed arguments from the command-line interface, containing user-configured
                 options such as paths for input/output files, atom group selections,
                 interaction flags, and computational settings.
    :type args: argparse.Namespace
    :return: This function does not return a value; it performs the analysis workflow
             and writes results to the specified output.
    :rtype: None
    :raises RuntimeError: If any essential operation such as atom selection, trajectory
                          processing, or frame selection fails due to invalid or insufficient
                          input data.
    """
    from pharmacon.logger import setup_logger, get_logger, header, subheader
    from pharmacon.command_line.trajectory.pl_interactions import (
        _detect_atoms, _extract_atom_indices, _worker_process_frames, _estimate_memory_and_warn,
    )

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
    header("Protein–Protein Interactions")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selections
    subheader("Selections")
    log.info("Protein 1      : %s", args.protein1)
    log.info("Protein 2      : %s", args.protein2)
    log.info("Water          : %s", args.water or "(none)")

    # Frame range
    subheader("Frame Range")
    log.info("Begin          : %d", args.begin)
    log.info("End            : %s", args.end if args.end is not None else "last")
    log.info("Step           : %d", args.step)

    # Enabled interactions
    subheader("Enabled Interactions")
    for interaction in args.enabled_interactions:
        log.info("  + %s", interaction)

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

    log.trace("Selecting atoms for Protein 1...")
    group1_atoms = u.select_atoms(args.protein1)
    if group1_atoms.n_atoms == 0:
        raise RuntimeError(f"Protein 1 selection has 0 atoms - selection='{args.protein1}'")
    log.info("Selected %d atoms for Protein 1.", group1_atoms.n_atoms)

    log.trace("Selecting atoms for Protein 2...")
    group2_atoms = u.select_atoms(args.protein2)
    if group2_atoms.n_atoms == 0:
        raise RuntimeError(f"Protein 2 selection has 0 atoms - selection='{args.protein2}'")
    log.info("Selected %d atoms for Protein 2.", group2_atoms.n_atoms)

    log.trace("Selecting water molecules...")
    water_selection = args.water
    if water_selection is None:
        log.debug("No water selection specified. Selecting water molecules based on internal logic...")
        water_selection = WATER_RESIDUES_MDA_SELECTION_STR
        log.debug("Selected water molecules based on internal logic. Selection='%s'", water_selection)
    else:
        _ = u.select_atoms(water_selection)

    # Build the full frame list
    all_frames = list(range(begin, end, step))
    n_analysis_frames = len(all_frames)

    if n_analysis_frames == 0:
        raise RuntimeError("No frames to process. Check --begin, --end, --step.")

    n_workers = min(args.workers, n_analysis_frames)

    # Disable flags dict (shared across workers)
    disable_flags = dict(
        disable_hydrophobic=args.disable_hydrophobic,
        disable_hbonds=args.disable_hbonds,
        disable_water_bridges=args.disable_water_bridges,
        disable_pi_stacking=args.disable_pi_stacking,
        disable_pi_cation=args.disable_pi_cation,
        disable_ionic=args.disable_ionic,
        disable_metal_complexation=args.disable_metal,
        disable_halogen=args.disable_halogen,
    )

    log.info("Calculating interactions...")
    calc_start_time = time.time()

    if n_workers <= 1:
        # Single-worker: original sequential path (no SQLite overhead)
        log.trace("Converting Protein 1 atoms to rdkit object...")
        log.trace("Converting Protein 2 atoms to rdkit object...")
        log.info("Detecting atoms...")
        atom_groups = _detect_atoms(u, group1_atoms, group2_atoms, "Protein 1", "Protein 2")

        frame_kwargs = dict(**disable_flags, **atom_groups, u=u)

        with PTAFile(args.output, overwrite=args.overwrite,
                     command="Trajectory Analysis", subcommand="pp-interactions") as pta:
            pta.add_file_metadata(metadata={"description": "Protein-Protein Interaction Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "protein1": str(args.protein1),
                                            "protein2": str(args.protein2),
                                            "water": str(args.water),
                                            "begin": str(begin), "end": str(end), "step": str(step),
                                            "disable_hydrophobic": str(args.disable_hydrophobic),
                                            "disable_hbonds": str(args.disable_hbonds),
                                            "disable_ionic": str(args.disable_ionic),
                                            "disable_halogen": str(args.disable_halogen),
                                            "disable_water_bridges": str(args.disable_water_bridges),
                                            "disable_pi_stacking": str(args.disable_pi_stacking),
                                            "disable_pi_cation": str(args.disable_pi_cation),
                                            "disable_metal": str(args.disable_metal),
                                            "is_merged": "False",
                                            "total_frames": str(n_frames),
                                            "workers": "1",
                                            })
            pta.create_group("pp_interactions")
            for ts in u.trajectory[begin:end:step]:
                frame_kwargs["box"] = ts.dimensions
                frame_contacts = interactions_process_frame(**frame_kwargs)
                flat = itertools.chain.from_iterable(frame_contacts)
                flat = deduplicate_interactions(flat)
                pta.write_frame_interactions(frame_index=ts.frame,
                                             group_name="pp_interactions",
                                             interactions=flat,
                                             overwrite=True)
                log.trace("Processed frame #%d", ts.frame)

            calc_end_time = time.time()
            calc_time = calc_end_time - calc_start_time
            log.info("Calculations completed. Total time taken: %.2f seconds.", calc_time)

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="pp_interactions",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="pp-interactions",
                                                   protein1=args.protein1,
                                                   protein2=args.protein2,
                                                   water=str(args.water),
                                                   begin=begin, end=end, step=step,
                                                   disable_hydrophobic=args.disable_hydrophobic,
                                                   disable_hbonds=args.disable_hbonds,
                                                   disable_ionic=args.disable_ionic,
                                                   disable_halogen=args.disable_halogen,
                                                   disable_water_bridges=args.disable_water_bridges,
                                                   disable_pi_stacking=args.disable_pi_stacking,
                                                   disable_pi_cation=args.disable_pi_cation,
                                                   disable_metal=args.disable_metal,
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
            pta.add_group_metadata(group_name="pp_interactions", metadata={"completed": "True"})

    else:
        # Multi-worker: SQLite temp storage
        import json
        import sqlite3
        import math
        from concurrent.futures import ProcessPoolExecutor, as_completed

        log.info("Parallel mode: %d workers for %d frames.", n_workers, n_analysis_frames)

        # Memory warning
        _estimate_memory_and_warn(
            n_atoms=u.atoms.n_atoms,
            n_workers=n_workers,
            n_analysis_frames=n_analysis_frames,
            log=log,
        )

        # Detect atoms ONCE on the main universe, then extract indices
        log.info("Detecting atoms (once, on main universe)...")
        atom_groups = _detect_atoms(u, group1_atoms, group2_atoms, "Protein 1", "Protein 2")
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
        frame_chunks = [ch for ch in frame_chunks if ch]

        log.info("Dispatching %d workers...", len(frame_chunks))
        for i, chunk in enumerate(frame_chunks):
            log.debug("  Worker %d: frames %d–%d (%d frames)", i, chunk[0], chunk[-1], len(chunk))

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
                    log_file=args.log,
                    file_logging_level=args.file_logging_level,
                ): i
                for i, chunk in enumerate(frame_chunks)
            }
            for future in as_completed(futures):
                wid = futures[future]
                _, n_done = future.result()
                log.info("Worker %d finished (%d frames processed).", wid, n_done)

        calc_end_time = time.time()
        calc_time = calc_end_time - calc_start_time
        log.info("All workers completed. Total calculation time: %.2f seconds.", calc_time)

        # Read SQLite and write to PTA frame-by-frame
        log.info("Writing interactions from temp database to PTA file...")
        conn = sqlite3.connect(str(db_path))

        with PTAFile(args.output, overwrite=args.overwrite,
                     command="Trajectory Analysis", subcommand="pp-interactions") as pta:
            pta.add_file_metadata(metadata={"description": "Protein-Protein Interaction Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "protein1": str(args.protein1),
                                            "protein2": str(args.protein2),
                                            "water": str(args.water),
                                            "begin": str(begin), "end": str(end), "step": str(step),
                                            "disable_hydrophobic": str(args.disable_hydrophobic),
                                            "disable_hbonds": str(args.disable_hbonds),
                                            "disable_ionic": str(args.disable_ionic),
                                            "disable_halogen": str(args.disable_halogen),
                                            "disable_water_bridges": str(args.disable_water_bridges),
                                            "disable_pi_stacking": str(args.disable_pi_stacking),
                                            "disable_pi_cation": str(args.disable_pi_cation),
                                            "disable_metal": str(args.disable_metal),
                                            "is_merged": "False",
                                            "total_frames": str(n_frames),
                                            "workers": str(n_workers),
                                            })
            pta.create_group("pp_interactions")

            for fidx in all_frames:
                cursor = conn.execute(
                    "SELECT data FROM interactions WHERE frame_index = ? ORDER BY rowid",
                    (fidx,),
                )
                interactions = [tuple(json.loads(row[0])) for row in cursor]
                pta.write_frame_interactions(frame_index=fidx,
                                             group_name="pp_interactions",
                                             interactions=interactions,
                                             overwrite=True)
                log.debug("Wrote frame %d (%d interactions)", fidx, len(interactions))

            conn.close()

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="pp_interactions",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="pp-interactions",
                                                   protein1=args.protein1,
                                                   protein2=args.protein2,
                                                   water=str(args.water),
                                                   begin=begin, end=end, step=step,
                                                   disable_hydrophobic=args.disable_hydrophobic,
                                                   disable_hbonds=args.disable_hbonds,
                                                   disable_ionic=args.disable_ionic,
                                                   disable_halogen=args.disable_halogen,
                                                   disable_water_bridges=args.disable_water_bridges,
                                                   disable_pi_stacking=args.disable_pi_stacking,
                                                   disable_pi_cation=args.disable_pi_cation,
                                                   disable_metal=args.disable_metal,
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
            pta.add_group_metadata(group_name="pp_interactions", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Protein-protein interaction analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Protein-Protein Interaction Analysis completed successfully. Time taken: %.2f seconds.", total_time)
