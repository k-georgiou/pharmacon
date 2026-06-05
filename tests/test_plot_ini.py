"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for the plot INI examples under examples/plot_ini/.

Verifies:
  - Every example INI parses without error and produces a non-empty dict.
  - Scalar coercion (bool/int/float) works on representative keys.
  - Alias resolution: all aliases for a settings class are registered
    in PlotSettingsBase._registry and map to that class.
  - The master `all_plots.ini` contains a section recognized by every
    PLI/PPI/PCA settings class.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pharmacon.utils.ini import read_ini, namespace_to_dict
from pharmacon.constants.plots import (
    PlotSettingsBase,
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinLigandInteractionsStackedColumn2PlotSettings,
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsPieCharts1PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
    ProteinProteinInteractionsTimelinePairs,
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsStackedColumnSettings,
    PlotUniversalSettings,
    PCAPlotTimeSeriesSettings,
    PCAPlotScatterSettings,
    PCAPlotVarianceRatioSettings,
    PCAPlotFESHeatmapSettings,
    PCAPlotProbabilityHeatmapSettings,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = REPO_ROOT / "examples" / "plot_ini"

ALL_EXAMPLE_INIS = sorted(EXAMPLES_DIR.glob("*.ini"))


class TestExampleInisParse:
    def test_examples_dir_exists(self):
        assert EXAMPLES_DIR.is_dir(), f"missing examples directory: {EXAMPLES_DIR}"

    def test_at_least_one_example(self):
        assert ALL_EXAMPLE_INIS, "expected at least one example INI"

    @pytest.mark.parametrize("ini_path", ALL_EXAMPLE_INIS,
                             ids=lambda p: p.name)
    def test_round_trip(self, ini_path):
        ns = read_ini(ini_path)
        d = namespace_to_dict(ns)
        assert isinstance(d, dict)
        # Every example INI defines at least one section
        sections = {k: v for k, v in d.items() if isinstance(v, dict)}
        assert sections, f"no sections in {ini_path.name}"


class TestScalarCoercion:
    def test_pli_stacked_column_1_coerced(self):
        ini = EXAMPLES_DIR / "pli_stacked_column_1.ini"
        ns = read_ini(ini)
        d = namespace_to_dict(ns)
        # Section name is one of the aliases — find a section that has fig_dpi.
        section = None
        for k, v in d.items():
            if isinstance(v, dict) and "fig_dpi" in v:
                section = v
                break
        assert section is not None, f"no section with fig_dpi in {ini}"
        # fig_dpi should be an int after coercion
        assert isinstance(section["fig_dpi"], int)
        # fig_transparent should be coerced to bool
        if "fig_transparent" in section:
            assert isinstance(section["fig_transparent"], bool)


CLASSES_WITH_ALIASES = [
    ProteinLigandInteractionsStackedColumn1PlotSettings,
    ProteinLigandInteractionsStackedColumn2PlotSettings,
    ProteinLigandInteractionsHeatmap1PlotSettings,
    ProteinLigandInteractionsHeatmap2PlotSettings,
    ProteinLigandInteractionsPieCharts1PlotSettings,
    ProteinLigandInteractionsLigandMonitorSettings,
    ProteinProteinInteractionsTimelinePairs,
    ProteinProteinInteractionsHeatmap,
    ProteinProteinInteractionsStackedColumnSettings,
    PlotUniversalSettings,
    PCAPlotTimeSeriesSettings,
    PCAPlotScatterSettings,
    PCAPlotVarianceRatioSettings,
    PCAPlotFESHeatmapSettings,
    PCAPlotProbabilityHeatmapSettings,
]


class TestAliasResolution:
    @pytest.mark.parametrize("cls", CLASSES_WITH_ALIASES,
                             ids=lambda c: c.__name__)
    def test_class_has_aliases(self, cls):
        aliases = getattr(cls, "alias", ())
        assert aliases, f"{cls.__name__} declares no aliases"

    @pytest.mark.parametrize("cls", CLASSES_WITH_ALIASES,
                             ids=lambda c: c.__name__)
    def test_aliases_registered(self, cls):
        registry = PlotSettingsBase._registry
        registered = {a.upper() for a in PlotSettingsBase.get_all_aliases()}
        for a in cls.alias:
            assert a.upper() in registered, (
                f"alias {a!r} of {cls.__name__} not in registry"
            )

    @pytest.mark.parametrize("cls", CLASSES_WITH_ALIASES,
                             ids=lambda c: c.__name__)
    def test_aliases_map_to_correct_class(self, cls):
        for a in cls.alias:
            # Registry keys/lookup are case-insensitive; we exercise
            # both the declared form and uppercased.
            mapped = PlotSettingsBase._registry.get(a)
            if mapped is None:
                mapped = PlotSettingsBase._registry.get(a.upper())
            assert mapped is cls, (
                f"alias {a!r} resolved to {mapped} (expected {cls.__name__})"
            )


class TestAllPlotsIni:
    @pytest.fixture(scope="class")
    def all_sections(self):
        path = EXAMPLES_DIR / "all_plots.ini"
        ns = read_ini(path)
        return namespace_to_dict(ns)

    @pytest.mark.parametrize("cls", CLASSES_WITH_ALIASES,
                             ids=lambda c: c.__name__)
    def test_every_settings_class_has_a_section(self, all_sections, cls):
        # all_plots.ini sections may use any alias; assert at least one of
        # the class's aliases appears as a key (case-insensitive).
        upper_keys = {str(k).upper() for k in all_sections.keys()}
        aliases_upper = {a.upper() for a in cls.alias}
        common = upper_keys & aliases_upper
        assert common, (
            f"no section in all_plots.ini for {cls.__name__}; "
            f"expected one of {sorted(aliases_upper)}"
        )


class TestFromDictTolerance:
    def test_unknown_key_ignored(self):
        # Adding a key that doesn't exist on the dataclass should not raise.
        cls = ProteinLigandInteractionsStackedColumn1PlotSettings
        overrides = {"this_key_does_not_exist": 42, "fig_dpi": 150}
        instance = cls.from_dict(overrides)
        # Real key still applied
        assert int(instance.fig_dpi) == 150
        # Unknown key did not survive
        assert not hasattr(instance, "this_key_does_not_exist")
