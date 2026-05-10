"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

structure properties — Compute structural/chemical properties of a molecule.
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

from pharmacon.analyzer.properties import (read_smi_file, read_mol2_file, generate_morgan_fingerprint,
                                           generate_topological_torsion_fingerprint, generate_maccs_keys,
                                           generate_atom_pair_fingerprint, count_total_atoms, get_molecular_weight,
                                           get_logp, get_net_charge, get_molecular_volume, get_num_rotatable_bonds,
                                           get_tpsa, get_stereo_centers, get_number_of_rings, get_number_of_aromatic_rings,
                                           get_element_counts, get_fragment_counts, fingerprint_to_string
                                           )



__all__ = [
    "SUBCOMMAND_NAME",
    "SUMMARY",
    "build_parser",
    "validate",
    "run",
]


SUBCOMMAND_NAME = "properties"
SUMMARY = "Compute structural/chemical properties of a molecule."

_INPUT_SUFFIXES: frozenset[str] = frozenset({".mol2", ".sdf", ".smi"})
_OUTPUT_SUFFIXES: frozenset[str] = frozenset({".psa"})

_EPILOG = """\
Examples:

  pharmacon structure properties -i ligands.smi  -o ligand.psa

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Builds and returns a subcommand parser.

    This function adds a subcommand parser to the given subparsers object.
    The parser enables handling command-line arguments for input options, output
    options, and logging features. It allows specifying the required input file,
    optional output file, overwrite behavior, and logging levels. You can further
    customize or extend this parser as needed.

    :param subparsers: An instance of `argparse._SubParsersAction` used to define
        subcommand parsers.
    :param parents: A list of parent parsers that provide shared arguments.
        Defaults to None.
    :return: The created `argparse.ArgumentParser` instance configured with
        specific options and arguments.
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
        help="Input structure file (.mol2, .sdf, .smi). [REQUIRED]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-o", "--output",
        required=False,
        metavar="FILE",
        default="properties.psa",
        help="Output Pharmacon Structure Analysis file (default=properties.psa). [OPTIONAL]",
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
                     default="properties.log",
                     help="File to log to (default: properties.log). [OPTIONAL]", )
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
    Validates the provided argparse.Namespace object, ensuring values meet specified criteria
    and updates the namespace in place with validated data. Validates input and output files,
    logging configuration, and overwrite flag. The function raises errors if critical conditions
    such as file duplication or invalid flags are encountered.

    :param args: The argparse.Namespace object containing command-line arguments to validate.
    :type args: argparse.Namespace
    :raises ValidationError: If certain validation rules are violated, such as input and output
        files being the same or invalid file paths/flags.
    :raises TypeError: If argument types do not meet expected types, such as invalid type
        for overwrite flag or logging level.
    :return: None
    """
    data: dict[str, object] = vars(args).copy()

    # Input file
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        _INPUT_SUFFIXES,
    )
    data["input"] = input_file

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Output file
    output: Path = validate_output_file(
        data.get("output", "properties.psa"),
        "output",
        overwrite=overwrite,
        allowed_suffixes=_OUTPUT_SUFFIXES,
    )
    data["output"] = output

    if output == input_file:
        raise ValidationError("Output file must not be the same as the input file.")

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "properties.log"),
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

    if log_path == input_file:
        raise ValidationError("Log file must not be the same as the input file.")
    if log_path == output:
        raise ValidationError("Log file must not be the same as the output file.")

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the structure properties analysis process, which involves reading molecule
    data from the specified input file, processing the molecular properties, generating
    fingerprints, and writing the analyzed data to the output file in a PSA format.

    :param args: Command-line arguments specifying input/output files, logging
        levels, and other configurations.
        - `args.input`: Path to the input file containing molecule specifications
          (.smi or .mol2 format).
        - `args.output`: Path to the output PSA file for storing the analysis.
        - `args.log`: Path to the log file for logging execution details.
        - `args.overwrite`: Flag indicating whether the output file may be
          overwritten if it already exists.
        - `args.terminal_logging_level`: Logging level for terminal output.
        - `args.file_logging_level`: Logging level for logging to file.

    :return: None
    """
    import json
    import time
    from typing import Dict

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
    header("Structure Properties")

    log.info("Input          : %s", args.input)
    log.info("Output         : %s", args.output)
    log.info("Log file       : %s", args.log)

    # Run
    subheader("Running Analysis")
    begin_time = time.time()

    master_data_dict: Dict[str, Dict[str, Dict[str, str]]] = {}
    calc_start_time = time.time()

    input_file: Path = args.input
    log.info("Reading input file: %s", input_file)
    extension = input_file.suffix.lower()

    if extension == ".smi":
        mol_names, mol_list = read_smi_file(input_file)
    elif extension == ".mol2":
        mol_names, mol_list = read_mol2_file(input_file)
    elif extension == ".sdf":
        # TODO: wire up an .sdf reader (e.g. RDKit Chem.SDMolSupplier).
        raise NotImplementedError(
            "SDF input is not yet supported; use .smi or .mol2 for now."
        )
    else:
        raise ValueError(f"Unsupported input file extension: {extension}")

    log.info("Loaded %d molecule(s) from %s", len(mol_list), input_file.name)

    file_data: Dict[str, Dict[str, str]] = {}

    for mol_name, mol in zip(mol_names, mol_list):
        log.debug("Processing molecule: %s", mol_name)

        morgan_fp = generate_morgan_fingerprint(mol)
        topol_fp = generate_topological_torsion_fingerprint(mol)
        maccs_fp = generate_maccs_keys(mol)
        atom_pair_fp = generate_atom_pair_fingerprint(mol)

        file_data[mol_name] = {
            "total_atoms": str(count_total_atoms(mol)),
            "molecular_weight": str(get_molecular_weight(mol)),
            "logP": str(get_logp(mol)),
            "net_charge": str(get_net_charge(mol)),
            "volume": str(get_molecular_volume(mol)),
            "rotatable_bonds": str(get_num_rotatable_bonds(mol)),
            "tpsa": str(get_tpsa(mol)),
            "stereo_centers": str(get_stereo_centers(mol)),
            "rings": str(get_number_of_rings(mol)),
            "aromatic_rings": str(get_number_of_aromatic_rings(mol)),
            "elements_dictionary": json.dumps(get_element_counts(mol), sort_keys=True),
            "fragments_dictionary": json.dumps(get_fragment_counts(mol), sort_keys=True),
            "morgan_fingerprint": fingerprint_to_string(morgan_fp),
            "topological_torsion_fingerprint": fingerprint_to_string(topol_fp),
            "maccs_keys": fingerprint_to_string(maccs_fp),
            "atom_pair_fingerprint": fingerprint_to_string(atom_pair_fp),
        }

    master_data_dict[input_file.name] = file_data

    calc_time = time.time() - calc_start_time
    log.info("Calculations completed. Total time taken: %.6f seconds.", calc_time)

    # Write PSA file
    log.info("Writing PSA file data...")

    with PSAFile(args.output, overwrite=args.overwrite,
                 command="Structure Analysis", subcommand="properties") as psa:
        psa.add_file_metadata(metadata={
            "description": "Structure Properties Analysis",
            "input_file": str(input_file),
            "is_merged": "False",
            "n_molecules": str(len(file_data)),
        })

        psa.create_group("properties")
        for source_file, molecules in master_data_dict.items():
            src_group = f"properties/{source_file}"
            psa.create_group(src_group)
            for mol_name, props in molecules.items():
                mol_group = f"{src_group}/{mol_name}"
                psa.create_group(mol_group)
                psa.add_group_metadata(group_name=mol_group,
                                       metadata=props, overwrite=True)
                log.trace("Wrote properties for molecule: %s", mol_name)

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
        psa.add_group_metadata(group_name="properties",
                               metadata={"completed": "True"})

    total_time = time.time() - begin_time

    # Done
    header("Done")
    log.info("Structure properties analysis complete.")
    log.info("Results written to: %s", args.output)
    log.info("Structure Properties Analysis completed successfully. Time taken: %.2f seconds.", total_time)
