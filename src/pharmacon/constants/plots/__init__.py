"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Plot settings package.

Historically lived in a single ``constants/plots.py`` module; split into one
submodule per plot family (``pli``, ``ppi``, ``universal``, ``pca``) sharing
``_base``. All public names are re-exported from this package so that
existing call sites ``from pharmacon.constants.plots import X`` keep working.
"""
from ._base import (
    VALID_EXTENSIONS,
    VALID_FONT_WEIGHTS,
    AVAILABLE_FONTS,
    VALID_LINE_STYLES,
    VALID_LEGEND_LOCS,
    PlotSettingsBase,
    logger,
)
from .pli import (
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinLigandInteractionsStackedColumn2PlotSettings,
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsPieCharts1PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
)
from .ppi import (
    ProteinProteinInteractionsTimelinePairs,
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsStackedColumnSettings,
)
from .universal import PlotUniversalSettings
from .rmsf import RMSFPlotSettings
from .pca import (
    PCAPlotTimeSeriesSettings,
    PCAPlotScatterSettings,
    PCAPlotVarianceRatioSettings,
    PCAPlotFESHeatmapSettings,
    PCAPlotProbabilityHeatmapSettings,
)


__all__ = [
    "VALID_EXTENSIONS",
    "VALID_FONT_WEIGHTS",
    "AVAILABLE_FONTS",
    "VALID_LINE_STYLES",
    "VALID_LEGEND_LOCS",
    "PlotSettingsBase",
    "logger",
    "ProteinLigandInteractionsStackedColumn1PlotSettings",
    "ProteinLigandInteractionsStackedColumn2PlotSettings",
    "ProteinLigandInteractionsHeatmap1PlotSettings",
    "ProteinLigandInteractionsHeatmap2PlotSettings",
    "ProteinLigandInteractionsPieCharts1PlotSettings",
    "ProteinLigandInteractionsLigandMonitorSettings",
    "ProteinProteinInteractionsTimelinePairs",
    "ProteinProteinInteractionsHeatmap",
    "ProteinProteinInteractionsStackedColumnSettings",
    "PlotUniversalSettings",
    "RMSFPlotSettings",
    "PCAPlotTimeSeriesSettings",
    "PCAPlotScatterSettings",
    "PCAPlotVarianceRatioSettings",
    "PCAPlotFESHeatmapSettings",
    "PCAPlotProbabilityHeatmapSettings",
]
