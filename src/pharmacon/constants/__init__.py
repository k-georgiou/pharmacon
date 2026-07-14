"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

"""

from typing import Final
from .formats import (supported_trajectory_analysis_topology_formats, supported_structure_analysis_topology_formats,
                      supported_trajectory_analysis_trajectory_formats, supported_pharmacon_ini_formats
                      )
from .smarts import (AA3, AA3_to_AA1, NA_RES, WATER_RESIDUES, ION_RES_TO_ELT,
                     PROTEIN_CARBONS, SIDECHAIN_OXY, SIDECHAIN_NITRO, NUC_OXY, TWO_LETTER,
                     BB_RE, HYDROPHOBIC_PATTERNS, HYDROPHOBIC_ELEMENTS,
                     HYDROGEN_BOND_DONOR_PATTERNS, HYDROGEN_BOND_ACCEPTOR_PATTERNS,
                     HYDROGEN_BOND_DONOR_ELEMENTS, HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
                     NEGATIVELY_CHARGED_PATTERNS, POSITIVELY_CHARGED_PATTERNS,
                     METAL_ELEMENTS, HALOGEN_PATTERNS, HALOGEN_ELEMENTS
                     )

# Plot setting classes are exposed lazily via PEP 562 __getattr__ so that
# importing pharmacon.constants (for __version__, BASE_PHARMACON_META, etc.)
# does not pull in the full plots subpackage.
# Mapping: public name on pharmacon.constants -> real attribute on pharmacon.constants.plots
_LAZY_PLOTS: Final[dict] = {
    "PLIStackedSettings1": "ProteinLigandInteractionsStackedColumn1PlotSettings",
    "PLIStackedSettings2": "ProteinLigandInteractionsStackedColumn2PlotSettings",
    "PLIHeatmapSettings1": "ProteinLigandInteractionsHeatmap1PlotSettings",
    "PLIHeatmapSettings2": "ProteinLigandInteractionsHeatmap2PlotSettings",
    "PLIPieChartsSettings1": "ProteinLigandInteractionsPieCharts1PlotSettings",
    "PLILigandMonitorSettings": "ProteinLigandInteractionsLigandMonitorSettings",
    "PPITimelinePairsSettings": "ProteinProteinInteractionsTimelinePairs",
    "PPIHeatmapSettings": "ProteinProteinInteractionsHeatmap",
    "PPIStackedColumnSettings": "ProteinProteinInteractionsStackedColumnSettings",
    "PlotUniversalSettings": "PlotUniversalSettings",
    "PCAPlotTimeSeriesSettings": "PCAPlotTimeSeriesSettings",
    "PCAPlotScatterSettings": "PCAPlotScatterSettings",
    "PCAPlotVarianceRatioSettings": "PCAPlotVarianceRatioSettings",
    "PCAPlotFESHeatmapSettings": "PCAPlotFESHeatmapSettings",
    "PCAPlotProbabilityHeatmapSettings": "PCAPlotProbabilityHeatmapSettings",
    "PlotSettingsBase": "PlotSettingsBase",
}


def __getattr__(name):
    if name in _LAZY_PLOTS:
        from . import plots as _plots
        value = getattr(_plots, _LAZY_PLOTS[name])
        globals()[name] = value  # cache for subsequent lookups
        return value
    raise AttributeError(f"module 'pharmacon.constants' has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_PLOTS.keys()))


__version__: Final[str]     = "1.0.1"
__author__: Final[str]      = "Kyriakos Georgiou"
__manuscript__: Final[str]  = "https://doi.org/10.1021/acs.jcim.6c00837"
__github__: Final[str]      = "https://github.com/k-georgiou/pharmacon"

# On-disk HDF5 layout version. Bump when the PTA/PSA schema changes in a way
# that older readers would misinterpret. Independent of __version__.
__schema_version__: Final[str] = "1.0"


BASE_PHARMACON_META: Final[dict] = {
    "pharmacon_version": __version__,
    "schema_version": __schema_version__,
    "file_type": "",
    "created_at": "",
    "command": "",
    "subcommand": "",
    "signature": "",
    "fingerprint": "",
}


__all__ = ["__version__", "__schema_version__", "__author__", "__manuscript__", "__github__",
           "BASE_PHARMACON_META",
           "supported_trajectory_analysis_topology_formats",
           "supported_structure_analysis_topology_formats",
           "supported_trajectory_analysis_trajectory_formats",
           "supported_pharmacon_ini_formats",
           "PLIStackedSettings1", "PLIStackedSettings2",
           "PLIHeatmapSettings1", "PLIHeatmapSettings2",
           "PLIPieChartsSettings1", "PLILigandMonitorSettings",
           "PPITimelinePairsSettings", "PPIHeatmapSettings",
           "PPIStackedColumnSettings",
           "PlotUniversalSettings", "PCAPlotTimeSeriesSettings",
           "PCAPlotScatterSettings", "PCAPlotVarianceRatioSettings",
           "PCAPlotFESHeatmapSettings", "PCAPlotProbabilityHeatmapSettings",
           "PlotSettingsBase"]

