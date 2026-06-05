"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for plot-data builders on mock PTAs.

Coverage:
  - build_pli_normal_data        (modes/<m>/table; non-merged)
  - build_pli_merged_stacked_data (modes_merged/<m>/table; merged)
  - PPI: same builders with group_name='pp_interactions'
  - PCA: _get_all_pca_components_from_file, _normalize_pca_component_pairs,
         _collect_pca_xy_from_file
  - Universal (rmsd/angles/distances) reads through plot_pta_timeseries_from_file
    are covered indirectly via the smoke layer; here we just verify the
    mock fixture writes records the readers parse without error.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import (
    build_pli_pta,
    build_pli_merged_pta,
    build_ppi_pta,
    build_ppi_merged_pta,
    build_pca_pta,
    build_universal_pta,
    DEFAULT_PLI_ROWS,
)

from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.plotter.interactions import (
    build_pli_normal_data,
    build_pli_merged_stacked_data,
    _PLI_STACK_ORDER,
)
from pharmacon.plotter.universal import (
    _get_all_pca_components_from_file,
    _normalize_pca_component_pairs,
    _collect_pca_xy_from_file,
)


class TestBuildPliNormalData:
    def test_basic_shape(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, interactions, values = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False, renumber_int=None,
            )
        assert values.shape == (len(residues), len(interactions))
        # Every (residue, interaction) cell must be a finite float >= 0.
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0)

    def test_residues_sorted_by_resid(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, _, _ = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False, renumber_int=None,
            )
        resids = [r[1] for r in residues]
        assert resids == sorted(resids)

    def test_interactions_in_canonical_order(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            _, interactions, _ = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False, renumber_int=None,
            )
        # Every label in `interactions` must appear in canonical order
        # (matches _PLI_STACK_ORDER for the subset that is present).
        present_canonical = [k for k in _PLI_STACK_ORDER if k in interactions]
        assert interactions == present_canonical

    def test_threshold_drops_low_frequencies(self, tmp_path):
        # Above threshold=0.5 only HYDROPHOBIC (0.80), HYDROGEN-BOND (0.55),
        # PI-STACKING (0.66) should remain among the defaults.
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, interactions, _ = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.5, aa3_to_aa1=False, renumber=False, renumber_int=None,
            )
        labels = set(interactions)
        assert labels.issubset({"hydrophobic", "hydrogen_bonds", "pi_stacking"})
        # And at least one residue survived
        assert len(residues) > 0

    def test_aa3_to_aa1_conversion(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, _, _ = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=True, renumber=False, renumber_int=None,
            )
        # Residue names should be single-letter after conversion.
        for resname, _, _, _ in residues:
            assert len(resname) == 1, f"expected 1-letter, got {resname!r}"

    def test_renumber_starts_at_supplied_int(self, tmp_path):
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, _, _ = build_pli_normal_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=True, renumber_int=1000,
            )
        resids = [r[1] for r in residues]
        # All renumbered ids must be >= 1000 and increase contiguously.
        assert min(resids) == 1000
        assert resids == sorted(resids)


class TestBuildPliMergedStackedData:
    def test_basic_shape(self, tmp_path):
        pta = build_pli_merged_pta(tmp_path / "pli_merged.pta", n_files=4)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, interactions, values, errors = build_pli_merged_stacked_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False,
                renumber_int=None, debug=False,
            )
        assert values.shape == errors.shape == (len(residues), len(interactions))
        assert np.all(np.isfinite(values))
        assert np.all(np.isfinite(errors))
        # std should be reasonable (we set std = mean * 0.1 in the mock)
        assert np.all(errors <= values + 1e-9)

    def test_canonical_interaction_order(self, tmp_path):
        pta = build_pli_merged_pta(tmp_path / "pli_merged.pta", n_files=3)
        with PharmaconPTAFile(pta, mode="r") as f:
            _, interactions, _, _ = build_pli_merged_stacked_data(
                f, group_name="pl_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False,
                renumber_int=None, debug=False,
            )
        present_canonical = [k for k in _PLI_STACK_ORDER if k in interactions]
        assert interactions == present_canonical

    def test_missing_modes_merged_raises(self, tmp_path):
        # build_pli_pta writes /modes/, not /modes_merged/ — expect RuntimeError.
        pta = build_pli_pta(tmp_path / "pli.pta", n_frames=100)
        with PharmaconPTAFile(pta, mode="r") as f:
            with pytest.raises(RuntimeError, match="Dataset not found"):
                build_pli_merged_stacked_data(
                    f, group_name="pl_interactions", mode_name="mode1",
                    threshold=0.0, aa3_to_aa1=False, renumber=False,
                    renumber_int=None, debug=False,
                )


class TestBuilderForPpi:
    def test_normal_data_works_for_ppi(self, tmp_path):
        pta = build_ppi_pta(tmp_path / "ppi.pta", n_frames=50)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, interactions, values = build_pli_normal_data(
                f, group_name="pp_interactions", mode_name="mode2",
                threshold=0.0, aa3_to_aa1=False, renumber=False, renumber_int=None,
            )
        assert values.shape == (len(residues), len(interactions))
        assert len(interactions) > 0

    def test_merged_data_works_for_ppi(self, tmp_path):
        pta = build_ppi_merged_pta(tmp_path / "ppi_merged.pta", n_files=2)
        with PharmaconPTAFile(pta, mode="r") as f:
            residues, interactions, values, errors = build_pli_merged_stacked_data(
                f, group_name="pp_interactions", mode_name="mode1",
                threshold=0.0, aa3_to_aa1=False, renumber=False,
                renumber_int=None, debug=False,
            )
        assert values.shape == errors.shape == (len(residues), len(interactions))


class TestPcaHelpers:
    def test_get_all_components_sorted_unique(self, tmp_path):
        pta = build_pca_pta(tmp_path / "pca.pta", n_frames=5, n_components=4)
        with PharmaconPTAFile(pta, mode="r") as f:
            pcs = _get_all_pca_components_from_file(f, group_name="pca")
        assert pcs == [1, 2, 3, 4]

    def test_get_all_components_too_few_raises(self, tmp_path):
        pta = build_pca_pta(tmp_path / "pca.pta", n_frames=3, n_components=1,
                            variance_ratios=[1.0])
        with PharmaconPTAFile(pta, mode="r") as f:
            with pytest.raises(ValueError, match="At least two PCA components"):
                _get_all_pca_components_from_file(f, group_name="pca")

    def test_normalize_component_pairs_single(self):
        pairs = _normalize_pca_component_pairs(
            (1, 2), allow_multiple=False, all_components=[1, 2, 3])
        assert pairs == [(1, 2)]

    def test_normalize_component_pairs_list(self):
        pairs = _normalize_pca_component_pairs(
            [1, 2, 3], allow_multiple=True, all_components=[1, 2, 3, 4])
        assert pairs == [(1, 2), (1, 3), (2, 3)]

    def test_normalize_component_pairs_none_with_multiple(self):
        pairs = _normalize_pca_component_pairs(
            None, allow_multiple=True, all_components=[1, 2, 3])
        assert pairs == [(1, 2), (1, 3), (2, 3)]

    def test_normalize_component_pairs_none_single(self):
        pairs = _normalize_pca_component_pairs(
            None, allow_multiple=False, all_components=[1, 2, 3])
        assert pairs == [(1, 2)]

    def test_normalize_component_pairs_invalid(self):
        with pytest.raises(ValueError):
            _normalize_pca_component_pairs(
                "garbage", allow_multiple=True, all_components=[1, 2])

    def test_collect_xy_returns_per_frame_values(self, tmp_path):
        pta = build_pca_pta(tmp_path / "pca.pta", n_frames=10, n_components=3)
        with PharmaconPTAFile(pta, mode="r") as f:
            x, y = _collect_pca_xy_from_file(f, group_name="pca", pc_x=1, pc_y=2)
        assert x.shape == y.shape == (10,)
        assert np.all(np.isfinite(x))
        assert np.all(np.isfinite(y))


class TestUniversalFixturesAreReadable:
    """The actual plot function is hit in the smoke layer. Here we just
    verify the mock writes records the reader will parse without error."""

    def test_rmsd_mock_has_frames(self, tmp_path):
        pta = build_universal_pta(tmp_path / "rmsd.pta", group="rmsd",
                                  series={"backbone": [1.0, 1.1, 1.2, 1.3, 1.4]})
        with PharmaconPTAFile(pta, mode="r") as f:
            frames = list(f._iter_frames("rmsd"))
        assert frames == [0, 1, 2, 3, 4]

    def test_angles_mock_has_nested_records(self, tmp_path):
        pta = build_universal_pta(
            tmp_path / "ang.pta", group="angles",
            series={"phi1": {"backbone": [60.0, 61.0, 62.0]}})
        import json
        with PharmaconPTAFile(pta, mode="r") as f:
            dset = f.file["angles/frame_0/angles"]
            rec = json.loads(dset[0])
        assert rec["label"] == "phi1"
        assert rec["kind"] == "backbone"
        assert rec["value"] == pytest.approx(60.0)

    def test_distances_mock_records(self, tmp_path):
        pta = build_universal_pta(
            tmp_path / "dis.pta", group="distances",
            series={"ca-ca": {"min_dist": [3.5, 3.6, 3.7]}})
        import json
        with PharmaconPTAFile(pta, mode="r") as f:
            dset = f.file["distances/frame_2/distances"]
            rec = json.loads(dset[0])
        assert rec["label"] == "ca-ca"
        assert rec["method"] == "min_dist"
        assert rec["distance"] == pytest.approx(3.7)
