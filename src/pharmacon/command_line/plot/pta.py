"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

plot pta — Render plots from a Pharmacon Trajectory Analysis file.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, Tuple, Type

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


SUBCOMMAND_NAME = "pta"
SUMMARY = "Render plots from a Pharmacon Trajectory Analysis (.pta) file."

# Analysis group detection
_INTERACTION_GROUPS = frozenset({"pl_interactions", "pp_interactions", "hbonds"})
_DISTANCE_GROUPS = frozenset({"distances"})
_ANGLE_GROUPS = frozenset({"angles"})
_RMSD_GROUPS = frozenset({"rmsd"})
_PCA_GROUPS = frozenset({"pca"})

_INI_SUFFIXES = frozenset({".ini", ".inp", ".conf"})

_EPILOG = """\
Examples:

  pharmacon plot pta -i merged.pta   -o ./plots/ -c custom.ini -mw 5

"""


def build_parser(subparsers: argparse._SubParsersAction,  # type: ignore[type-arg]
                 parents: List[argparse.ArgumentParser] | None = None) -> argparse.ArgumentParser:
    """
    Build an argparse parser for handling the command-line arguments.

    This function creates a parser and configures various groups of command-line
    options for input, output, and logging configurations. The configured
    argparse parser is prepared to handle command-line arguments required for
    processing Pharmacon Trajectory Analysis (.pta) files, setting plot options
    via configuration files, and managing output and logging behaviors.

    :param subparsers: Sub-command parsers action object allowing addition of new
        sub-command parsers.
    :type subparsers: argparse._SubParsersAction
    :param parents: A list of parent argument parsers from which arguments will be
        inherited. Default=None.
    :type parents: List[argparse.ArgumentParser] | None
    :return: Configured argument parser for the functionality.
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
        help="Input Pharmacon Trajectory Analysis (.pta) file. [REQUIRED]",
    )
    inp.add_argument(
        "-c", "--config",
        required=False,
        metavar="FILE",
        default=None,
        help="Pharmacon INI config file with plot settings overrides. [OPTIONAL]",
    )

    # Output Options
    out = parser.add_argument_group("Output Options")
    out.add_argument(
        "-o", "--output",
        required=True,
        metavar="DIR",
        help="Output directory where plots will be written. [REQUIRED]",
    )
    out.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output directory. Default='False'. Type='bool' [OPTIONAL]",
    )
    out.add_argument(
        "-mw", "--maxwarnings",
        required=False,
        type=int,
        default=0,
        metavar="N",
        help="Maximum number of plot-settings coercion warnings tolerated per "
             "plot before aborting. Default=0. Type='int' [OPTIONAL]",
    )

    # Logging Options
    log = parser.add_argument_group("Logging Options")
    log.add_argument("-l", "--log",
                     required=False,
                     metavar="FILE",
                     default="plot_pta.log",
                     help="File to log to (default: plot_pta.log). [OPTIONAL]", )
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
    Validates and prepares the input arguments for execution, ensuring that provided
    parameters such as input/output files, logging settings, and PTA metadata are correctly
    set, consistent, and meet the expected requirements. The function performs transformations
    on argument values and raises errors if invalid data or configurations are detected.

    :param args: Parsed command-line arguments encapsulated in an argparse.Namespace instance.
    :type args: argparse.Namespace

    :raises ValidationError: If any input, configuration, or metadata parameter is
        invalid, inconsistent, or results in errors during validation.
    :return: None
    """
    import h5py

    from pharmacon.utils.fingerprint import create_pharmacon_signature
    from pharmacon.utils.identifiers import validate_mda_artifact_token

    data: dict[str, object] = vars(args).copy()

    # Input PTA file
    input_file: Path = validate_existing_input_file(
        data.get("input"),
        "input",
        [".pta"],
    )
    data["input"] = input_file

    # Overwrite flag
    overwrite: bool = validate_bool_flag(data.get("overwrite", False), "overwrite")
    data["overwrite"] = overwrite

    # Max warnings
    try:
        maxwarnings = int(data.get("maxwarnings", 0))
    except (TypeError, ValueError) as exc:
        raise ValidationError(
            f"--maxwarnings must be an integer, got '{data.get('maxwarnings')}'."
        ) from exc
    if maxwarnings < 0:
        raise ValidationError(
            f"--maxwarnings must be >= 0, got {maxwarnings}."
        )
    data["maxwarnings"] = maxwarnings

    # Output directory
    import shutil

    output_dir = Path(str(data.get("output", ""))).expanduser().resolve()

    if output_dir.exists():
        if output_dir.is_dir():
            if overwrite:
                shutil.rmtree(output_dir)
                output_dir.mkdir(parents=True)
            elif any(output_dir.iterdir()):
                raise ValidationError(
                    f"Output directory is not empty: '{output_dir}'. "
                    "Use --overwrite to replace it."
                )
        else:
            raise ValidationError(
                f"Output path is not a directory: '{output_dir}'."
            )
    else:
        output_dir.mkdir(parents=True)

    data["output"] = output_dir

    # Logging
    log_path: Path = validate_output_file(
        data.get("log", "plot_pta.log"), "log",
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

    if log_path == input_file:
        raise ValidationError("Log file must not be the same as the input file.")

    # Config file (optional)
    config_overrides: Dict[str, Dict[str, Any]] = {}
    config_value = data.get("config")
    if config_value is not None:
        config_file: Path = validate_existing_input_file(
            config_value,
            "config",
            _INI_SUFFIXES,
        )
        data["config"] = config_file

        from pharmacon.utils.ini import namespace_to_dict, read_ini

        try:
            ini_ns: SimpleNamespace = read_ini(config_file)
        except (ImportError, FileNotFoundError, ValueError) as exc:
            raise ValidationError(f"Failed to read config file: {exc}") from exc

        ini_dict = namespace_to_dict(ini_ns)
        for section_name, section_body in ini_dict.items():
            if not isinstance(section_body, dict):
                # Top-level scalar keys are ignored — only aliased sections
                # apply plot settings.
                continue
            config_overrides[section_name.upper()] = _flatten_section(section_body)
    else:
        data["config"] = None

    data["config_overrides"] = config_overrides

    # PTA file metadata integrity checks
    try:
        f = h5py.File(input_file, "r")
    except Exception as exc:
        raise ValidationError(f"Cannot open PTA file: {exc}") from exc

    try:
        attrs = dict(f.attrs)

        artifact_status = str(attrs.get("artifact_status", "")).strip().upper()
        if not artifact_status:
            raise ValidationError(
                "PTA file is missing 'artifact_status' metadata — "
                "the file may be incomplete or corrupted."
            )
        if artifact_status != "SUCCESS":
            raise ValidationError(
                f"PTA file artifact_status is '{artifact_status}' (expected 'SUCCESS') — "
                "the file may be corrupted or the analysis did not complete."
            )

        artifact_token = str(attrs.get("artifact_token", "")).strip()
        if not artifact_token:
            raise ValidationError(
                "PTA file is missing 'artifact_token' metadata — "
                "the file may be corrupted."
            )

        blueprint = str(attrs.get("blueprint", "")).strip()
        if not blueprint:
            raise ValidationError(
                "PTA file is missing 'blueprint' metadata — "
                "the file may be corrupted."
            )

        token_valid = validate_mda_artifact_token(
            artifact_token=artifact_token,
            blueprint=blueprint,
            secret="trajectory_analysis",
            namespace="pharmacon",
        )
        if not token_valid:
            raise ValidationError(
                "PTA file artifact_token does not match the blueprint — "
                "the file may have been tampered with or corrupted."
            )

        signature = str(attrs.get("signature", "")).strip()
        fingerprint = str(attrs.get("fingerprint", "")).strip()
        command = str(attrs.get("command", "")).strip()
        subcommand = str(attrs.get("subcommand", "")).strip()

        if not signature or not fingerprint:
            raise ValidationError(
                "PTA file is missing 'signature' and/or 'fingerprint' metadata — "
                "the file may be corrupted."
            )

        if command and subcommand:
            expected_sig = create_pharmacon_signature(
                format_name="pta",
                command=command,
                subcommand=subcommand,
            )
            if expected_sig.signature != signature:
                raise ValidationError(
                    f"PTA file signature mismatch.\n"
                    f"  Expected : {expected_sig.signature}\n"
                    f"  Found    : {signature}\n"
                    "The file may have been tampered with or corrupted."
                )
            if expected_sig.fingerprint != fingerprint:
                raise ValidationError(
                    f"PTA file fingerprint mismatch.\n"
                    f"  Expected : {expected_sig.fingerprint}\n"
                    f"  Found    : {fingerprint}\n"
                    "The file may have been tampered with or corrupted."
                )

        version = str(attrs.get("pharmacon_version", "")).strip()
        if not version:
            raise ValidationError(
                "PTA file is missing 'pharmacon_version' metadata — "
                "the file may be corrupted."
            )

        from pharmacon.constants import __version__
        runtime_version = str(__version__)
        try:
            f_parts = tuple(map(int, version.split(".")))
            r_parts = tuple(map(int, runtime_version.split(".")))
            if f_parts > r_parts:
                import warnings
                warnings.warn(
                    f"PTA file was created with Pharmacon v{version} "
                    f"but the current runtime is v{runtime_version}. "
                    "Some features may not be supported.",
                    UserWarning,
                    stacklevel=2,
                )
        except (ValueError, AttributeError):
            pass

        is_merged = str(attrs.get("is_merged", "False")).strip().lower() == "true"
        data["is_merged"] = is_merged

        for key in f:
            obj = f[key]
            if not hasattr(obj, "attrs") or not hasattr(obj, "keys"):
                continue
            completed = str(obj.attrs.get("completed", "")).strip().lower()
            if completed != "true":
                raise ValidationError(
                    f"Group '{key}' does not have completed=True — "
                    f"the analysis for this group may not have finished."
                )

    finally:
        f.close()

    # Write validated values back into the namespace in-place.
    for key, value in data.items():
        setattr(args, key, value)


def run(args: argparse.Namespace) -> None:
    """
    Executes the given arguments to generate plots and manage the PTA outputs and configurations.

    This function sets up logger configurations, processes input arguments, reads PTA files,
    and generates specified plots for analysis groups present in the PTA file. It handles
    universal timeseries plots as per the specified configurations and produces PCA (Principal Component
    Analysis) plots where applicable.

    :param args: Parsed arguments containing file and plot settings.
    :type args: argparse.Namespace
    """
    from rich.box import SIMPLE
    from rich.console import Console
    from rich.table import Table

    from pharmacon.constants.plots import (PCAPlotFESHeatmapSettings,
                                           PCAPlotProbabilityHeatmapSettings,
                                           PCAPlotScatterSettings,
                                           PCAPlotTimeSeriesSettings,
                                           PCAPlotVarianceRatioSettings,
                                           PlotSettingsBase,
                                           PlotUniversalSettings,
                                           ProteinLigandInteractionsHeatmap1PlotSettings,
                                           ProteinLigandInteractionsHeatmap2PlotSettings,
                                           ProteinLigandInteractionsLigandMonitorSettings,
                                           ProteinLigandInteractionsPieCharts1PlotSettings,
                                           ProteinLigandInteractionsStackedColumn1PlotSettings,
                                           ProteinLigandInteractionsStackedColumn2PlotSettings,
                                           ProteinProteinInteractionsHeatmap,
                                           ProteinProteinInteractionsStackedColumnSettings,
                                           ProteinProteinInteractionsTimelinePairs
                                           )
    from pharmacon.fileio.pta import PharmaconPTAFile
    from pharmacon.logger import get_logger, header, setup_logger, subheader
    from pharmacon.utils.workspace import PharmaconWorkspace as Workspace
    from pharmacon.plotter.interactions import (plot_protein_ligand_interactions_heatmap_1_from_file,
                                                plot_protein_ligand_interactions_heatmap_2_from_file,
                                                plot_protein_ligand_interactions_ligand_monitor_from_file,
                                                plot_protein_ligand_interactions_pie_charts_from_file,
                                                plot_protein_ligand_interactions_stacked_column_1_from_file,
                                                plot_protein_ligand_interactions_stacked_column_2_from_file,
                                                plot_protein_protein_heatmap_freq_from_file,
                                                plot_protein_protein_interactions_stacked_column_from_file,
                                                plot_protein_protein_timeline_pairs_from_file,
                                                )
    from pharmacon.plotter.universal import (plot_pca_fes_from_file,
                                             plot_pca_probability_from_file,
                                             plot_pca_scatter_from_file,
                                             plot_pca_timeseries_from_file,
                                             plot_pca_variance_ratio_from_file,
                                             plot_pta_timeseries_from_file,
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

    console = Console()

    output_dir: Path = args.output
    is_merged: bool = args.is_merged
    maxwarnings: int = args.maxwarnings
    overrides_by_alias: Dict[str, Dict[str, Any]] = args.config_overrides

    # Header
    header("Plot PTA")

    log.info("Input          : %s", args.input)
    log.info("Output dir     : %s", output_dir)
    log.info("Config         : %s", args.config if args.config else "(none)")
    log.info("Overwrite      : %s", args.overwrite)
    log.info("Max warnings   : %d", maxwarnings)
    log.info("Merged file    : %s", is_merged)
    log.info("Log file       : %s", args.log)

    cfg = Table(show_header=False, box=SIMPLE, pad_edge=False, padding=(0, 1))
    cfg.add_column("Key", style="bold cyan", no_wrap=True, min_width=20)
    cfg.add_column("Value", style="white")
    cfg.add_row("Input", str(args.input))
    cfg.add_row("Output directory", str(output_dir))
    cfg.add_row("Config", str(args.config) if args.config else "(none)")
    cfg.add_row("Overwrite", str(args.overwrite))
    cfg.add_row("Max warnings", str(maxwarnings))
    cfg.add_row("Merged file", str(is_merged))
    cfg.add_row("Log file", str(args.log))
    console.print(cfg)

    # Workspace
    log.debug("Creating the pharmacon workspace...")
    ws = Workspace(is_tmp_dir_needed=False)
    log.debug("Pharmacon workspace created at: %s", ws.working_directory)

    rendered: list[str] = []
    skipped: list[tuple[str, str]] = []

    def _settings_for(cls: Type[PlotSettingsBase]) -> PlotSettingsBase | None:
        """
        Construct and return an instance of the specified PlotSettingsBase subclass with
        its configuration derived from provided overrides. If settings construction
        fails or exceeds the allowed number of coercion warnings, None is returned.

        :param cls: The class type inheriting from PlotSettingsBase for which settings
            are to be created.
        :return: An instance of the PlotSettingsBase subclass initialized with overridden
            settings or None if the configuration fails or exceeds the allowed warnings.
        """
        overrides = _lookup_overrides(cls, overrides_by_alias)
        try:
            instance = cls.from_dict(overrides)
        except Exception as exc:
            skipped.append((cls.__name__, f"settings construction failed: {exc}"))
            return None
        n_warnings = int(getattr(instance, "_current_warnings", 0) or 0)
        if n_warnings > maxwarnings:
            skipped.append((
                cls.__name__,
                f"{n_warnings} coercion warning(s) exceed --maxwarnings {maxwarnings}",
            ))
            return None
        return instance

    def _invoke(label: str, fn: Callable[..., None], **kwargs: Any) -> None:
        """
        Invoke a callable function and handle any exceptions that occur during its execution.
        This function attempts to execute the provided callable with the given keyword arguments.
        In case of an exception, it logs the error message along with the label and ensures the
        process continues without interruption.

        :param label: A string representing the label associated with the function being invoked.
        :param fn: A callable object which will be executed with the provided keyword arguments.
        :param kwargs: Additional keyword arguments passed to the callable function.
        :return: This function does not return a value.
        """
        try:
            fn(**kwargs)
        except Exception as exc:
            skipped.append((label, f"{type(exc).__name__}: {exc}"))
            console.print(f"  [red]Failed:[/red] {label} — {exc}")
            return
        rendered.append(label)
        console.print(f"  [green]Rendered:[/green] {label}")

    # Open the PTA file
    with PharmaconPTAFile(args.input, mode="r") as pta:

        subheader("File Validity")
        pta.print_file_validity_status()
        pta.validate_version(strict=False)

        subheader("Discovering Analysis Groups")
        groups = pta.get_groups()
        if not groups:
            console.print("[yellow]No analysis groups found in PTA file.[/yellow]")
            header("Done")
            return

        console.print(f"  Found [bold]{len(groups)}[/bold] top-level group(s): "
                      f"{', '.join(sorted(groups))}")

        # rmsd / angles / distances (unified timeseries)
        for group_name in sorted(groups):
            if group_name not in (_RMSD_GROUPS | _ANGLE_GROUPS | _DISTANCE_GROUPS):
                continue

            subheader(f"Plotting: {group_name}")

            settings = _settings_for(PlotUniversalSettings)
            if settings is None:
                continue

            _invoke(
                f"{group_name}/timeseries",
                plot_pta_timeseries_from_file,
                pta_file=pta,
                group_name=group_name,
                settings=settings,
                out_dir=output_dir,
                is_merged=is_merged,
            )

        # PCA
        for group_name in sorted(groups):
            if group_name not in _PCA_GROUPS:
                continue

            subheader(f"Plotting: {group_name}")

            if is_merged:
                console.print(
                    "[yellow]  PCA plots are not supported on merged PTA files — "
                    "skipping.[/yellow]"
                )
                skipped.append((group_name, "PCA not supported on merged files"))
                continue

            pca_plots: Tuple[Tuple[str, Type[PlotSettingsBase], Callable[..., None]], ...] = (
                ("pca/timeseries",       PCAPlotTimeSeriesSettings,        plot_pca_timeseries_from_file),
                ("pca/scatter",          PCAPlotScatterSettings,           plot_pca_scatter_from_file),
                ("pca/variance_ratio",   PCAPlotVarianceRatioSettings,     plot_pca_variance_ratio_from_file),
                ("pca/probability",      PCAPlotProbabilityHeatmapSettings, plot_pca_probability_from_file),
                ("pca/fes",              PCAPlotFESHeatmapSettings,        plot_pca_fes_from_file),
            )
            for label, cls, fn in pca_plots:
                settings = _settings_for(cls)
                if settings is None:
                    continue
                _invoke(
                    label,
                    fn,
                    pta_file=pta,
                    group_name=group_name,
                    settings=settings,
                    out_dir=output_dir,
                )

        # Protein–Ligand interactions
        for group_name in sorted(groups):
            if group_name != "pl_interactions":
                continue

            subheader(f"Plotting: {group_name}")

            mode_names = _discover_modes(pta, group_name, is_merged)

            # Per-mode plots that exist in both merged and non-merged files.
            for mode_name in mode_names:
                sc1 = _settings_for(ProteinLigandInteractionsStackedColumn1PlotSettings)
                if sc1 is not None:
                    _invoke(
                        f"{group_name}/{mode_name}/stacked_column_1",
                        plot_protein_ligand_interactions_stacked_column_1_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        mode_name=mode_name,
                        settings=sc1,
                        out_dir=output_dir,
                        attach_to_name=mode_name,
                        is_merged=is_merged,
                    )

            # Non-merged-only plots.
            if not is_merged:
                heatmap1 = _settings_for(ProteinLigandInteractionsHeatmap1PlotSettings)
                if heatmap1 is not None:
                    _invoke(
                        f"{group_name}/heatmap_1",
                        plot_protein_ligand_interactions_heatmap_1_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=heatmap1,
                        out_dir=output_dir,
                        is_merged=False,
                    )

                heatmap2 = _settings_for(ProteinLigandInteractionsHeatmap2PlotSettings)
                if heatmap2 is not None:
                    _invoke(
                        f"{group_name}/heatmap_2",
                        plot_protein_ligand_interactions_heatmap_2_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=heatmap2,
                        out_dir=output_dir,
                        is_merged=False,
                    )

                pie = _settings_for(ProteinLigandInteractionsPieCharts1PlotSettings)
                if pie is not None:
                    _invoke(
                        f"{group_name}/pie_charts",
                        plot_protein_ligand_interactions_pie_charts_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=pie,
                        out_dir=output_dir,
                        is_merged=False,
                    )

                sc2 = _settings_for(ProteinLigandInteractionsStackedColumn2PlotSettings)
                if sc2 is not None:
                    _invoke(
                        f"{group_name}/stacked_column_2",
                        plot_protein_ligand_interactions_stacked_column_2_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=sc2,
                        out_dir=output_dir,
                        is_merged=False,
                    )

                monitor = _settings_for(ProteinLigandInteractionsLigandMonitorSettings)
                if monitor is not None:
                    _invoke(
                        f"{group_name}/ligand_monitor",
                        plot_protein_ligand_interactions_ligand_monitor_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=monitor,
                        out_dir=output_dir,
                        is_merged=False,
                    )

        # Protein–Protein interactions
        for group_name in sorted(groups):
            if group_name != "pp_interactions":
                continue

            subheader(f"Plotting: {group_name}")

            mode_names = _discover_modes(pta, group_name, is_merged)

            for mode_name in mode_names:
                ppi_heatmap = _settings_for(ProteinProteinInteractionsHeatmap)
                if ppi_heatmap is not None:
                    _invoke(
                        f"{group_name}/{mode_name}/heatmap",
                        plot_protein_protein_heatmap_freq_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        mode_name=mode_name,
                        settings=ppi_heatmap,
                        out_dir=output_dir,
                        attach_to_name=mode_name,
                        is_merged=is_merged,
                    )

                ppi_stacked = _settings_for(ProteinProteinInteractionsStackedColumnSettings)
                if ppi_stacked is not None:
                    _invoke(
                        f"{group_name}/{mode_name}/stacked_column",
                        plot_protein_protein_interactions_stacked_column_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        mode_name=mode_name,
                        settings=ppi_stacked,
                        out_dir=output_dir,
                        attach_to_name=mode_name,
                        is_merged=is_merged,
                    )

            if not is_merged:
                timeline = _settings_for(ProteinProteinInteractionsTimelinePairs)
                if timeline is not None:
                    _invoke(
                        f"{group_name}/timeline_pairs",
                        plot_protein_protein_timeline_pairs_from_file,
                        pta_file=pta,
                        group_name=group_name,
                        settings=timeline,
                        out_dir=output_dir,
                        is_merged=False,
                    )

    # Summary
    subheader("Summary")
    console.print(f"  Rendered [bold]{len(rendered)}[/bold] plot(s) to [bold]{output_dir}[/bold]")
    if skipped:
        console.print(f"  Skipped  [bold]{len(skipped)}[/bold] plot(s):")
        for label, reason in skipped:
            console.print(f"    [yellow]- {label}[/yellow]: {reason}")

    # Unknown config sections warning
    unknown_sections = _unknown_config_sections(overrides_by_alias)
    if unknown_sections:
        console.print(
            f"  [yellow]Warning:[/yellow] config file contained "
            f"{len(unknown_sections)} section(s) that did not match any "
            f"known plot-settings alias: {', '.join(sorted(unknown_sections))}"
        )

    header("Done")


# Helpers

def _flatten_section(section: Dict[str, Any]) -> Dict[str, Any]:
    """
    Flatten a given section dictionary by removing nested dictionaries.

    This function processes a dictionary and extracts all key-value pairs where
    values are not dictionaries. It is useful for simplifying a structured dictionary
    by discarding nested subsections while retaining their key-value pairs.

    :param section: The dictionary section to be flattened.
    :type section: Dict[str, Any]
    :return: A flattened dictionary containing non-dictionary values only.
    :rtype: Dict[str, Any]
    """
    return {k: v for k, v in section.items() if not isinstance(v, dict)}


def _lookup_overrides(cls: Any, overrides_by_alias: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Looks up and retrieves override settings for a specified class based on its aliases. The method
    inspects the 'alias' attribute of the given class to search for matching override mappings in the
    provided dictionary. If a match is found, the corresponding dictionary is returned. If no
    match exists, an empty dictionary is returned.

    :param cls: The class to inspect for aliases.
    :type cls: Any
    :param overrides_by_alias: A dictionary where keys are alias names in uppercase, and values
        are dictionaries that contain override settings for the corresponding alias.
    :type overrides_by_alias: Dict[str, Dict[str, Any]]
    :return: A dictionary of override settings corresponding to the matched alias,
        or an empty dictionary if no match is found.
    :rtype: Dict[str, Any]
    """
    for alias in getattr(cls, "alias", ()):
        key = alias.upper()
        if key in overrides_by_alias:
            return overrides_by_alias[key]
    return {}


def _unknown_config_sections(overrides_by_alias: Dict[str, Dict[str, Any]]) -> List[str]:
    """
    Identifies and returns configuration sections that are not recognized from the given
    dictionary of overrides defined by aliases. The recognition is based on a predefined
    set of known aliases.

    :param overrides_by_alias: A dictionary where keys are section aliases as strings,
        and values are dictionaries representing configuration overrides for those aliases.
    :return: A list of string keys representing sections in `overrides_by_alias` that are
        not part of the predefined set of recognized aliases.
    """
    from pharmacon.constants import PlotSettingsBase

    known = {a.upper() for a in PlotSettingsBase.get_all_aliases()}
    return [section for section in overrides_by_alias if section not in known]


def _discover_modes(pta: Any, group_name: str, is_merged: bool) -> List[str]:
    """
    Discovers and returns a list of mode names from the specified PTA group.

    This function checks the specified PTA group in the given file for modes
    associated with the provided `group_name`. If the `is_merged` flag is True,
    it searches for merged modes. Otherwise, it searches for unmerged modes.

    :param pta: An object containing a file-like attribute that is searched for modes.
    :param group_name: The name of the group within the PTA file to search for modes.
    :param is_merged: A boolean flag indicating whether to look for merged modes.
    :return: A sorted list of discovered mode names. Returns an empty list if no
        modes are found.
    :rtype: List[str]
    """
    modes_parent = f"{group_name}/modes_merged" if is_merged else f"{group_name}/modes"
    if modes_parent not in pta.file:
        return []
    parent = pta.file[modes_parent]
    return sorted(name for name in parent if name.startswith("mode"))
