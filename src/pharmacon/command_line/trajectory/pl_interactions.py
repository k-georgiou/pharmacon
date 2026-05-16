"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

trajectory pl-interactions — Calculate protein-ligand interactions.
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


SUBCOMMAND_NAME = "pl-interactions"
SUMMARY = "Calculate protein-ligand interactions through a trajectory."


_EPILOG = """\
Examples:

  pharmacon trajectory pl-interactions -p topol.tpr -x traj.xtc -o pl_interactions.pta -prt "protein" -lig "resname LIG -w "resname WAT"
  
"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and adds a command-line subparser for the 'pl_interactions' program to an existing
    argparse._SubParsersAction instance. The parser includes groups of options for input and
    output file handling, selection of molecular elements, time range constraints, periodic
    boundary condition (PBC) transformations, parallelism settings, interaction toggles, and
    logging customization.

    The constructed parser is returned for further use.

    :param subparsers: The subparsers object from argparse to which the 'pl_interactions'
        command parser is added.
    :type subparsers: argparse._SubParsersAction
    :param parents: List of parent argparse.ArgumentParser objects whose arguments are
        inherited by the created parser. Defaults to None.
    :type parents: List[argparse.ArgumentParser], optional
    :return: The configured argparse.ArgumentParser instance added to the subparsers object.
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
        default="pl_interactions.pta",
        help="Output Pharmacon Trajectory Analysis file (default=pl_interactions.pta). [OPTIONAL]",
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
        "-prt", "--protein",
        metavar="SEL",
        default=None,
        required=False,
        help="MDAnalysis selection string for protein atoms (default: None). [REQUIRED]",
    )

    sel.add_argument(
        "-lig", "--ligand",
        metavar="SEL",
        default=None,
        required=True,
        help="MDAnalysis selection string for ligand atoms (default: None). [REQUIRED]",
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
                        help="Disable pi-stacking interactions (default: False). [OPTIONAL]",)
    toggle.add_argument("--disable-pi-cation",
                        required=False,
                        action="store_true",
                        help="Disable pi-cation interactions (default: False). [OPTIONAL]",)
    toggle.add_argument("--disable-ionic",
                        required=False,
                        action="store_true",
                        help="Disable ionic interactions (default: False). [OPTIONAL]",)
    toggle.add_argument("--disable-water-bridges",
                        required=False,
                        action="store_true",
                        help="Disable water bridging interactions (default: False). [OPTIONAL]",)
    toggle.add_argument("--disable-halogen",
                        required=False,
                        action="store_true",
                        help="Disable halogen bonding interactions (default: False). [OPTIONAL]",)
    toggle.add_argument("--disable-metal",
                        required=False,
                        action="store_true",
                        help="Disable metal ion interactions (default: False). [OPTIONAL]",)

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="pl_interactions.log",
                     help="File to log to (default: pl_interactions.log). [OPTIONAL]", )
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
                     help="Logging level for terminal logging (default: INFO). [OPTIONAL]",)

    return parser


def validate(args: argparse.Namespace) -> None:
    """
    Validates input arguments provided via an `argparse.Namespace`, ensuring
    appropriate file paths, argument consistency, and configurations for a
    trajectory analysis process.

    The validation process encompasses several tasks:
    - Ensures input files (topology and trajectory) are valid, exist, and are
      different from each other.
    - Verifies output files comply with necessary constraints, including overwrite
      flags, valid suffixes, and uniqueness compared to input files.
    - Normalizes selection strings for entities like protein, ligand, and water
      while handling optional and required configurations accordingly.
    - Validates frame ranges, boolean flags (e.g., disabling specific interactions),
      and worker configurations.
    - Generates and validates logging configurations, ensuring paths, levels, and
      uniqueness of file destinations are compliant.

    Raises `ValidationError` when input arguments violate constraints or result in
    logical inconsistencies.

    :param args: An `argparse.Namespace` object containing the arguments to be
        validated.
    :type args: argparse.Namespace
    :return: None. The function updates the `args` object in-place with validated
        and normalized values.
    :rtype: NoneType
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
        data.get("output", "pl_interactions.pta"),
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
    protein: str | None = normalize_selection(
        data.get("protein"),
        "protein",
        required=False,
        default=None,
    )
    ligand: str | None = normalize_selection(
        data.get("ligand"),
        "ligand",
        required=True,
        default=None,
    )
    water: str | None = normalize_selection(
        data.get("water"),
        "water",
        required=False,
        default=None,
    )

    assert ligand is not None

    data["protein"] = protein
    data["ligand"] = ligand
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
        data.get("log", "pl_interactions.log"),
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


def _detect_atoms(u, group1_atoms, group2_atoms, group1_label, group2_label):
    """
    Detects various chemical groupings and atom types in two defined sets of atoms within a Universe `u`.

    This function processes two groups of atoms, converts them into RDKit molecular objects, and analyzes
    the atoms based on several chemical characteristics such as hydrophobicity, charge, aromaticity, and
    metallicity. The result is a dictionary that maps each detected attribute to the corresponding atom
    group in each of the two atom sets.

    :param u: Molecular universe object containing atomic data.
    :param group1_atoms: Atoms corresponding to the first group.
    :param group2_atoms: Atoms corresponding to the second group.
    :param group1_label: Descriptive label for the first atom group.
    :param group2_label: Descriptive label for the second atom group.
    :return: A dictionary with detected atom types and groups within each atom set.
    :rtype: dict
    """
    rdkit_mol_group1, mapping_dict_group1 = convert_mda_to_rdkit(group1_atoms)
    rdkit_mol_group2, mapping_dict_group2 = convert_mda_to_rdkit(group2_atoms)

    return dict(
        hydrophobic_atoms_group1=detect_hydrophobic_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Hydrophobic atoms in {group1_label}"),
        hydrophobic_atoms_group2=detect_hydrophobic_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Hydrophobic atoms in {group2_label}"),
        hydrogen_bond_acceptor_atoms_group1=detect_hydrogen_bond_acceptor_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Hydrogen Bond Acceptor atoms in {group1_label}"),
        hydrogen_bond_acceptor_atoms_group2=detect_hydrogen_bond_acceptor_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Hydrogen Bond Acceptor atoms in {group2_label}"),
        hydrogen_bond_donor_atoms_group1=detect_hydrogen_bond_donor_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Hydrogen Bond Donor atoms in {group1_label}"),
        hydrogen_bond_donor_atoms_group2=detect_hydrogen_bond_donor_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Hydrogen Bond Donor atoms in {group2_label}"),
        aromatic_atoms_group1=detect_aromatic_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Aromatic rings in {group1_label}"),
        aromatic_atoms_group2=detect_aromatic_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Aromatic rings in {group2_label}"),
        positive_charged_atoms_group1=detect_positive_charge_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Positive Charge atoms in {group1_label}"),
        positive_charged_atoms_group2=detect_positive_charge_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Positive Charge atoms in {group2_label}"),
        negative_charged_atoms_group1=detect_negative_charge_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Negative Charge atoms in {group1_label}"),
        negative_charged_atoms_group2=detect_negative_charge_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Negative Charge atoms in {group2_label}"),
        halogen_atoms_group1=detect_halogen_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Halogen atoms in {group1_label}"),
        halogen_atoms_group2=detect_halogen_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Halogen atoms in {group2_label}"),
        metal_atoms_group1=detect_metal_atoms(
            u, rdkit_mol_group1, mapping_dict_group1, label=f"Metal atoms in {group1_label}"),
        metal_atoms_group2=detect_metal_atoms(
            u, rdkit_mol_group2, mapping_dict_group2, label=f"Metal atoms in {group2_label}"),
    )


_AROMATIC_KEYS: frozenset[str] = frozenset({
    "aromatic_atoms_group1",
    "aromatic_atoms_group2",
})


def _extract_atom_indices(atom_groups: dict) -> dict:
    """
    Extract atom indices from the provided atom groups dictionary.

    This function processes a dictionary of atom groups, classifying them based on
    specific keys and extracting their indices. For keys identified as aromatic,
    it processes each ring in the group to extract indices. Non-aromatic keys are
    treated differently by directly extracting the indices.

    :param atom_groups: A dictionary where the keys represent atom group categories
        and the values are objects (e.g., aromatic rings or other atom groups) that
        contain atom indices information.
    :type atom_groups: dict

    :return: A dictionary where each key corresponds to the input keys of the atom
        groups and the values are lists of indices extracted from the respective
        atom group objects.
    :rtype: dict
    """
    out: dict = {}
    for key, ag in atom_groups.items():
        if key in _AROMATIC_KEYS:
            out[key] = [ring.indices.tolist() for ring in ag]
        else:
            out[key] = ag.indices.tolist()
    return out


def _rebuild_atom_groups(u, atom_indices: dict) -> dict:
    """
    Rebuilds atom groups based on provided atom indices and their respective groups, distinguishing between
    aromatic and non-aromatic keys.

    Given a dictionary of atom indices grouped under specific keys, this function iterates through the dictionary
    and constructs atom groups either as a list of atoms for non-aromatic keys or as a list of atom subgroups for
    aromatic keys. The distinction between aromatic and non-aromatic keys is determined based on a predefined set
    of aromatic keys.

    :param u: Atom universe or object containing atom data.
    :param atom_indices: Dictionary where keys correspond to atom group categories, and values are lists of atom
        indices or indices arranged by their relevant groups (e.g., aromatic rings).
    :return: Dictionary with the same keys as the input atom_indices, where the values correspond to lists of
        atom objects (or subgroups for aromatic keys).
    :rtype: dict
    """
    out: dict = {}
    for key, indices in atom_indices.items():
        if key in _AROMATIC_KEYS:
            out[key] = [u.atoms[ring] for ring in indices]
        else:
            out[key] = u.atoms[indices]
    return out


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
                           file_logging_level: str = "TRACE"):
    """
    Processes specific frames of a molecular trajectory using a worker process, handling
    interactions between atoms and storing results in a database.

    :param worker_id: Identifier for the worker process.
    :type worker_id: int
    :param topology: Topology information of the molecular system.
    :param trajectory: Trajectory data of the molecular system.
    :param wrapping: Specifies whether to apply periodic boundary condition wrapping.
    :type wrapping: bool
    :param unwrapping: Specifies whether to unwrap the trajectory for processing.
    :type unwrapping: bool
    :param frame_indices: List of frame indices to process.
    :type frame_indices: list
    :param disable_flags: Dictionary of flags to disable specific calculation features.
    :type disable_flags: dict
    :param atom_indices: Dictionary containing pre-computed atomic group indices.
    :type atom_indices: dict
    :param db_path: Path to the SQLite database file where the interaction results will
        be stored.
    :param log_file: Optional log file path for the worker-specific logger.
    :type log_file: str or None
    :param file_logging_level: Logging level for the log file output. The default is
        "TRACE".
    :type file_logging_level: str
    :return: A tuple containing the worker's identifier and the number of frames it
        processed successfully.
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
    # (avoids expensive RDKit conversion + SMARTS matching per worker)
    atom_groups = _rebuild_atom_groups(u, atom_indices)

    frame_kwargs = dict(
        **disable_flags,
        **atom_groups,
        u=u,
    )

    conn = sqlite3.connect(str(db_path), timeout=60)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    for _n, fidx in enumerate(frame_indices, 1):
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

        log.debug("Worker %d: frame %d/%d (idx %d)", worker_id, _n, len(frame_indices), fidx)
        log.trace("Worker %d processed frame #%d", worker_id, fidx)

    conn.close()
    log.trace("Worker %d finished: %d frames.", worker_id, len(frame_indices))
    return worker_id, len(frame_indices)


def _estimate_memory_and_warn(*, n_atoms, n_workers, n_analysis_frames, log):
    """
    Estimates the memory usage for a parallel operation, providing warnings if the estimated
    memory requirement exceeds certain thresholds relative to the available system RAM.
    The function logs memory estimates, available RAM, and recommendations to adjust
    parameters when necessary.

    :param n_atoms: Number of atoms to be processed in the analysis.
    :type n_atoms: int
    :param n_workers: Number of workers used for the parallel operation.
    :type n_workers: int
    :param n_analysis_frames: Total number of analysis frames to process.
    :type n_analysis_frames: int
    :param log: Logger instance to record debug and warning messages.
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
    Executes the primary function of analyzing protein-ligand interactions within a given molecular
    dynamics trajectory using user-specified parameters and options. The function sets up a logger,
    reads input trajectory and topology files, performs atomic selections, processes frame ranges, and
    calculates specified interaction types. Results are written to an output file.

    :param args: Parsed command-line arguments containing user-specified options for analysis.
    :type args: argparse.Namespace
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
    header("Protein–Ligand Interactions")

    log.info("Topology       : %s", args.topology)
    log.info("Trajectory     : %s", args.trajectory)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Selections
    subheader("Selections")
    log.info("Protein        : %s", args.protein or "(auto-detect)")
    log.info("Ligand         : %s", args.ligand)
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
    ws = Workspace(is_tmp_dir_needed=need_tmp, cleanup_on_exit=True)
    log.debug("Pharmacon workspace created successfully.")
    # NOTE: Workspace registers atexit cleanup. Temp dir is removed on
    # normal exit, unhandled exceptions, and KeyboardInterrupt.
    # Only SIGKILL (kill -9) can bypass cleanup.

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

    log.trace("Selecting atoms for Protein...")
    group1_atoms = u.select_atoms(args.protein)
    if group1_atoms.n_atoms == 0:
        raise RuntimeError(f"Protein selection has 0 atoms - selection='{args.protein}'")
    log.info("Selected %d atoms for protein.", group1_atoms.n_atoms)

    log.trace("Selecting atoms for Ligand...")
    group2_atoms = u.select_atoms(args.ligand)
    if group2_atoms.n_atoms == 0:
        raise RuntimeError(f"Ligand selection has 0 atoms - selection='{args.ligand}'")
    log.info("Selected %d atoms for ligand.", group2_atoms.n_atoms)

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
        # ── Single-worker: original sequential path (no SQLite overhead) ─────
        log.trace("Converting protein atoms to rdkit object...")
        log.trace("Converting ligand atoms to rdkit object...")
        log.info("Detecting atoms...")
        atom_groups = _detect_atoms(u, group1_atoms, group2_atoms, "Protein", "Ligand")

        frame_kwargs = dict(**disable_flags, **atom_groups, u=u)

        with PTAFile(args.output, overwrite=args.overwrite,
                     command="Trajectory Analysis", subcommand="pl-interactions") as pta:
            pta.add_file_metadata(metadata={"description": "Protein-Ligand Interaction Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "protein": str(args.protein),
                                            "ligand": str(args.ligand),
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
            pta.create_group("pl_interactions")
            for _n, ts in enumerate(u.trajectory[begin:end:step], 1):
                frame_kwargs["box"] = ts.dimensions
                frame_contacts = interactions_process_frame(**frame_kwargs)
                flat = itertools.chain.from_iterable(frame_contacts)
                flat = deduplicate_interactions(flat)
                pta.write_frame_interactions(frame_index=ts.frame,
                                             group_name="pl_interactions",
                                             interactions=flat,
                                             overwrite=True)
                log.debug("Frame %d/%d (idx %d)", _n, n_analysis_frames, ts.frame)
                log.trace("Processed frame #%d", ts.frame)

            calc_end_time = time.time()
            calc_time = calc_end_time - calc_start_time
            log.info("Calculations completed. Total time taken: %.2f seconds.", calc_time)

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="pl_interactions",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="pl-interactions",
                                                   protein=args.protein,
                                                   ligand=args.ligand,
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
            pta.add_group_metadata(group_name="pl_interactions", metadata={"completed": "True"})

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

        # Detect atoms ONCE on the main universe, then extract indices
        log.info("Detecting atoms (once, on main universe)...")
        atom_groups = _detect_atoms(u, group1_atoms, group2_atoms, "Protein", "Ligand")
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
        # The Manager-backed Queue is picklable across spawn/fork contexts, so it
        # passes safely as a kwarg to ProcessPoolExecutor workers.
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
                     command="Trajectory Analysis", subcommand="pl-interactions") as pta:
            pta.add_file_metadata(metadata={"description": "Protein-Ligand Interaction Analysis",
                                            "topology_file": str(args.topology),
                                            "trajectory_file": str(args.trajectory),
                                            "protein": str(args.protein),
                                            "ligand": str(args.ligand),
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
            pta.create_group("pl_interactions")

            for fidx in all_frames:
                cursor = conn.execute(
                    "SELECT data FROM interactions WHERE frame_index = ? ORDER BY rowid",
                    (fidx,),
                )
                interactions = [tuple(json.loads(row[0])) for row in cursor]
                pta.write_frame_interactions(frame_index=fidx,
                                             group_name="pl_interactions",
                                             interactions=interactions,
                                             overwrite=True)
                log.debug("Wrote frame %d (%d interactions)", fidx, len(interactions))

            conn.close()

            log.info("Building interaction modes...")
            pta.build_interaction_modes(group_name="pl_interactions",
                                        begin=begin, end=end, step=step)
            log.debug("Interaction modes built successfully.")
            log.debug("Finalizing PTA file data")

            blueprint = str(generate_mda_blueprint(u=u,
                                                   command="Trajectory Analysis",
                                                   subcommand="pl-interactions",
                                                   protein=args.protein,
                                                   ligand=args.ligand,
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
            pta.add_group_metadata(group_name="pl_interactions", metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Protein-ligand interaction analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Protein-Ligand Interaction Analysis completed successfully. Time taken: %.2f seconds.", total_time)
