"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Pharmacon: A Molecular Dynamics Suite.
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
from .plots import (ProteinLigandInteractionsStackedColumn1PlotSettings as PLIStackedSettings1,
                    ProteinLigandInteractionsStackedColumn2PlotSettings as PLIStackedSettings2,
                    ProteinLigandInteractionsHeatmap1PlotSettings as PLIHeatmapSettings1,
                    ProteinLigandInteractionsHeatmap2PlotSettings as PLIHeatmapSettings2,
                    ProteinLigandInteractionsPieCharts1PlotSettings as PLIPieChartsSettings1,
                    ProteinLigandInteractionsLigandMonitorSettings as PLILigandMonitorSettings,
                    ProteinProteinInteractionsTimelinePairs as PPITimelinePairsSettings,
                    ProteinProteinInteractionsHeatmap as PPIHeatmapSettings,
                    ProteinProteinInteractionsStackedColumnSettings as PPIStackedColumnSettings,
                    PlotUniversalSettings,
                    PCAPlotTimeSeriesSettings,
                    PCAPlotScatterSettings,
                    PCAPlotVarianceRatioSettings,
                    PCAPlotFESHeatmapSettings,
                    PCAPlotProbabilityHeatmapSettings,
                    PlotSettingsBase,
                    )


__version__: Final[str]     = "1.0.0"
__author__: Final[str]      = "Kyriakos Georgiou"
__manuscript__: Final[str]  = "https://arxiv.org/abs/2302.04242"
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

