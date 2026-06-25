"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for the PPI residue×residue heatmap matrix builder.

Previously the heatmap reused the PLI per-residue builder
(`build_pli_normal_data`), which discards the partner residue
(`rec["key"][1]`), and then placed per-residue totals on the DIAGONAL only.
The result showed none of the actual residue-residue contacts, dropped every
residue that only appeared as the second partner, and made `symmetric` a
permanent no-op.

`build_ppi_heatmap_matrix` now uses BOTH residues of each key, populating the
off-diagonal pair frequencies.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make tests/helpers importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from helpers.mock_pta import build_ppi_pta

from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.plotter.interactions import (
    build_ppi_heatmap_matrix,
    plot_protein_protein_heatmap_freq_from_file,
)
from pharmacon.constants.plots import ProteinProteinInteractionsHeatmap as H

# The mock PPI rows are 4 distinct cross-chain residue pairs:
#   A:PHE45  ↔ B:VAL200   (0.40)
#   A:ASP102 ↔ B:LYS215   (0.55)
#   A:ARG150 ↔ B:GLU220   (0.22)
#   A:TYR77  ↔ B:TRP230   (0.18)


def _matrix(tmp_path, **overrides):
    pta = build_ppi_pta(tmp_path / "ppi.pta", n_frames=100)
    settings = H.from_dict(overrides)
    with PharmaconPTAFile(pta, mode="r") as f:
        return build_ppi_heatmap_matrix(
            f, group_name="pp_interactions", mode_name="mode1",
            settings=settings, is_merged=False,
        )


class TestPairAwareMatrix:
    def test_offdiagonal_populated_not_diagonal(self, tmp_path):
        # The core regression: contacts land OFF the diagonal, not on it.
        labels, m = _matrix(tmp_path)
        diag_nonzero = int(np.count_nonzero(np.diag(m)))
        off_nonzero = int((m > 0).sum()) - diag_nonzero
        assert off_nonzero == 4      # 4 residue-residue pairs
        assert diag_nonzero == 0     # no self-contacts on the diagonal

    def test_both_partners_appear_on_axis(self, tmp_path):
        # Residues that only ever appear as the 2nd partner (chain B) must
        # still appear — the old builder dropped them entirely.
        labels, _ = _matrix(tmp_path)
        assert len(labels) == 8
        for lbl in ("A:PHE45", "B:VAL200", "A:ASP102", "B:LYS215",
                    "A:ARG150", "B:GLU220", "A:TYR77", "B:TRP230"):
            assert lbl in labels

    def test_specific_pair_frequency(self, tmp_path):
        labels, m = _matrix(tmp_path)
        idx = {l: i for i, l in enumerate(labels)}
        assert m[idx["A:PHE45"], idx["B:VAL200"]] == pytest.approx(0.40)
        assert m[idx["A:ASP102"], idx["B:LYS215"]] == pytest.approx(0.55)
        # un-symmetrised: the mirror cell is still zero before symmetrise
        assert m[idx["B:VAL200"], idx["A:PHE45"]] == 0.0

    def test_matrix_asymmetric_so_symmetric_flag_matters(self, tmp_path):
        labels, m = _matrix(tmp_path)
        assert not np.allclose(m, m.T)                 # was a no-op before
        sym = m + m.T - np.diag(np.diag(m))
        assert np.allclose(sym, sym.T)
        # symmetrising mirrors every (i, j) into (j, i): doubles the nonzeros
        assert int((sym > 0).sum()) == 2 * int((m > 0).sum())

    def test_threshold_filters_pairs(self, tmp_path):
        # threshold above 0.22 / 0.18 keeps only the 0.40 and 0.55 pairs
        labels, m = _matrix(tmp_path, threshold=0.30)
        assert int((m > 0).sum()) == 2

    def test_aa3_to_aa1_labels(self, tmp_path):
        labels, _ = _matrix(tmp_path, aa3_to_aa1=True)
        assert "A:F45" in labels and "B:V200" in labels


@pytest.mark.slow
def test_render_ppi_heatmap_produces_file(tmp_path):
    pta = build_ppi_pta(tmp_path / "ppi.pta", n_frames=100)
    out = tmp_path / "out"
    out.mkdir()
    settings = H.from_dict({"fig_dpi": 80, "fig_size_width": 4, "fig_size_height": 4})
    with PharmaconPTAFile(pta, mode="r") as f:
        plot_protein_protein_heatmap_freq_from_file(
            f, group_name="pp_interactions", mode_name="mode1",
            settings=settings, out_dir=out, attach_to_name="mode1",
            is_merged=False,
        )
    files = [p for p in out.iterdir() if p.is_file()]
    assert files and all(p.stat().st_size > 0 for p in files)
