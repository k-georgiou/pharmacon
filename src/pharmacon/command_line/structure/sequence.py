"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

structure sequence — Extract amino-acid sequence(s) from a structure topology.
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


SUBCOMMAND_NAME = "sequence"
SUMMARY = "Extract amino-acid sequence(s) from a structure topology."

_INPUT_SUFFIXES: frozenset[str] = frozenset({".pdb", ".gro", ".crd"})
_OUTPUT_SUFFIXES: frozenset[str] = frozenset({".psa"})

_EPILOG = """\
Examples:

  pharmacon structure sequence -p protein.pdb  -o sequence.psa

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and configures a subparser for a specific subcommand.

    This function adds a subparser to the specified argument parsers group. It
    configures the input options, output options, and logging options required
    for the operation of the subcommand. Each group is tailored to define both
    required and optional arguments for the user.

    :param subparsers: A subparsers action object to which this command's
        parser will be added.
    :param parents: Optionally, a list of argparse.ArgumentParser instances that
        serve as parents for the created subparser. Defaults to None.
    :return: Configured instance of argparse.ArgumentParser for the subcommand.
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
        help="Input topology file (.pdb, .gro, .crd). [REQUIRED]",
    )
    inp.add_argument(
        "-sel", "--selection",
        required=False,
        metavar="STR",
        default="protein",
        help="MDAnalysis selection string identifying the protein atoms from "
             "which to extract the sequence (default: 'protein'). [OPTIONAL]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-o", "--output",
        required=False,
        metavar="FILE",
        default="sequence.psa",
        help="Output Pharmacon Structure Analysis file (default=sequence.psa). [OPTIONAL]",
    )
    out.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output file (default=False). [OPTIONAL]",
    )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="sequence.log",
                     help="File to log to (default: sequence.log). [OPTIONAL]", )
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
    Validates the provided `argparse.Namespace` object. This function performs multiple
    validation and sanitization checks on the input data, including validating file paths,
    parsing and verifying selection strings for MDAnalysis, and ensuring logical consistency
    among file-related inputs. All modifications are applied directly to the provided
    `Namespace`.

    :param args:
        An `argparse.Namespace` containing configuration data to be validated and updated.
        Expected keys and their values include:
        - `topology` (str | Path): The path to the topology file.
        - `selection` (str): A non-empty MDAnalysis selection string.
        - `overwrite` (bool): A flag indicating whether to overwrite output files.
        - `output` (str | Path): The path to the output file for results.
        - `log` (str | Path): The path to the log file for logging information.
        - `file_logging_level` (str): The logging level for file-based logging.
        - `terminal_logging_level` (str): The logging level for terminal-based logging.

    :raises ValidationError:
        If any validation check fails during processing. This includes:
        - Missing or invalid topology file.
        - Improper selection string for MDAnalysis.
        - Conflicting or illegal file paths for topology, output, and log files.
        - Invalid logging level values.

    :raises ImportError:
        If the required dependency `MDAnalysis` is not available in the environment.

    :return:
        The function does not return any value. The `args` object is updated in place
        with validated and parsed data.
    """
    data: dict[str, object] = vars(args).copy()

    # Topology file
    topology_file: Path = validate_existing_input_file(
        data.get("topology"),
        "topology",
        _INPUT_SUFFIXES,
    )
    data["topology"] = topology_file

    # Selection string
    selection_raw = data.get("selection", "")
    if not isinstance(selection_raw, str) or not selection_raw.strip():
        raise ValidationError(
            "Selection string must be a non-empty MDAnalysis selection."
        )
    selection: str = selection_raw.strip()

    try:
        import MDAnalysis as Mda
    except ImportError as exc:
        raise ValidationError(
            "MDAnalysis is required to validate the selection string. "
            "Install with: pip install MDAnalysis"
        ) from exc

    try:
        u = Mda.Universe(str(topology_file))
    except Exception as exc:
        raise ValidationError(
            f"Cannot load topology '{topology_file}' with MDAnalysis: {exc}"
        ) from exc

    try:
        atoms = u.select_atoms(selection)
    except Exception as exc:
        raise ValidationError(
            f"Invalid MDAnalysis selection '{selection}': {exc}"
        ) from exc

    if atoms.n_atoms == 0:
        raise ValidationError(
            f"Selection '{selection}' resolved to 0 atoms in topology "
            f"'{topology_file}'."
        )

    data["selection"] = selection

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Output file
    output: Path = validate_output_file(
        data.get("output", "sequence.psa"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=_OUTPUT_SUFFIXES,
    )
    data["output"] = output

    if output == topology_file:
        raise ValidationError("Output file must not be the same as the topology file.")

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "sequence.log"),
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

    if log_path == topology_file:
        raise ValidationError("Log file must not be the same as the topology file.")
    if log_path == output:
        raise ValidationError("Log file must not be the same as the output file.")

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the structure sequence analysis subcommand.

    The `run` function processes molecular structure data to extract sequence
    information for chains present in the topology file. It utilizes the MDAnalysis
    library for specifying atom selections and extracts sequence details based
    on user-provided configurations. The extracted sequences and metadata are
    written to a PSA file for further usage or analysis. Additionally, the function
    handles detailed logging for both debugging and runtime information.

    :param args: Parsed namespace containing the arguments passed by the user.
    :type args: argparse.Namespace
    :raises RuntimeError: If the provided atom selection resolves to 0 atoms
        in the specified topology.
    :return: None
    """
    import time
    import MDAnalysis as Mda

    from pharmacon.analyzer.sequence import get_sequence, sequence_dict_to_fasta
    from pharmacon.fileio import PSAFile
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
    header("Structure Sequence")

    log.info("Topology       : %s", args.topology)
    log.info("Selection      : %s", args.selection)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Run
    subheader("Running Analysis")
    begin_time = time.time()

    log.debug("Loading topology with MDAnalysis...")
    u = Mda.Universe(str(args.topology))
    log.debug("Loaded topology with %d atoms.", u.atoms.n_atoms)

    log.debug("Applying selection: %s", args.selection)
    selection_atoms = u.select_atoms(args.selection)
    if selection_atoms.n_atoms == 0:
        raise RuntimeError(
            f"Selection '{args.selection}' resolved to 0 atoms in topology "
            f"'{args.topology}'."
        )
    log.info("Selection matched %d atoms.", selection_atoms.n_atoms)

    # Build a restricted Universe containing only the selected atoms so that
    # analyzer.get_sequence() operates on the user's region of interest while
    # preserving its internal `protein and name CA` filter.
    log.debug("Building restricted universe from selection...")
    restricted_u = Mda.Merge(selection_atoms)

    log.info("Extracting per-chain sequences...")
    calc_start_time = time.time()
    sequences = get_sequence(restricted_u)
    calc_time = time.time() - calc_start_time
    log.info("Extracted %d chain(s) in %.6f seconds.", len(sequences), calc_time)

    for chainid, data in sorted(sequences.items()):
        aa3_list = data.get("aa3_list") or []
        first_aa = aa3_list[0] if aa3_list else "-"
        last_aa = aa3_list[-1] if aa3_list else "-"
        log.info(
            "  chain %-8s : %4d residues  (%s ... %s)",
            chainid, len(data["aa1_seq"]), first_aa, last_aa,
        )

    # FASTA preview at debug level
    for line in sequence_dict_to_fasta(sequences).splitlines():
        log.debug("  %s", line)

    # Write PSA file
    log.info("Writing PSA file data...")

    with PSAFile(args.output, overwrite=args.overwrite,
                 command="Structure Analysis", subcommand="sequence") as psa:
        psa.add_file_metadata(metadata={
            "description": "Structure Sequence Analysis",
            "topology_file": str(args.topology),
            "selection": args.selection,
            "is_merged": "False",
            "n_chains": str(len(sequences)),
        })

        psa.write_sequence(
            sequences=sequences,
            group_name="sequence",
            overwrite=True,
        )

        blueprint = "None"
        meta = {
            "completed": "True",
            "artifact_status": "SUCCESS",
            "artifact_status_code": 0,
            "blueprint": blueprint,
        }
        token = create_mda_artifact_token(
            blueprint=blueprint,
            secret="structure_analysis",
            namespace="pharmacon",
        )
        psa.add_file_metadata({
            **meta,
            "artifact_token": token,
            "artifact_token_version": "1",
        })
        psa.add_group_metadata(group_name="sequence",
                               metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Structure sequence analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Structure Sequence Analysis completed successfully. Time taken: %.2f seconds.", total_time)
