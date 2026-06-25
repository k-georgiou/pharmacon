"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for selective interaction-mode building.

h-bonds analysis produces only HYDROGEN-BOND contacts, so mode3 (which only
differs from mode1 in how it dedups HYDROPHOBIC contacts) is byte-identical to
mode1 — a redundant duplicate. The h-bonds command now builds only modes 1 & 2
(`mode3=False`), and ``build_interaction_modes`` writes a ``mode_definition``
metadata string that describes only the modes actually built.
"""

from __future__ import annotations

import h5py

from pharmacon.fileio.pta import PharmaconPTAFile


def _hbond_rec(resname1, resid1, resname2, resid2):
    # Atom-atom interaction record layout consumed by extract_mode_key.
    return ["HYDROGEN-BOND",
            1, "N", 1, "N", "BB", resname1, resid1, "A", "SYSTEM",
            2, "O", 2, "O", "BB", resname2, resid2, "A", "SYSTEM",
            3.0, True]


def _build(tmp_path, **mode_flags):
    path = tmp_path / "hb.pta"
    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand="h-bonds") as pta:
        pta.add_file_metadata({"is_merged": "False"})
        pta.create_group("hbonds")
        # Two frames; a residue pair recurs so mode1 (count all) > mode2 (dedup).
        for frame in (0, 1):
            pta.write_frame_interactions(
                frame_index=frame, group_name="hbonds",
                interactions=[_hbond_rec("ARG", 4, "THR", 118),
                              _hbond_rec("ARG", 4, "THR", 118),
                              _hbond_rec("ALA", 6, "ALA", 120)],
                overwrite=True,
            )
        pta.build_interaction_modes(group_name="hbonds", begin=0, end=1, step=1,
                                    **mode_flags)
    return path


def _modes_and_def(path):
    f = h5py.File(path, "r")
    modes = sorted(f["hbonds/modes"].keys())
    first = modes[0]
    mdef = f["hbonds/modes"][first]["table"].attrs.get("mode_definition")
    return modes, mdef


def test_mode3_false_builds_only_two_modes(tmp_path):
    modes, mdef = _modes_and_def(_build(tmp_path, mode3=False))
    assert modes == ["mode1", "mode2"]
    assert "mode3" not in mdef
    assert "mode1=count all rows" in mdef
    assert "mode2=deduplicate residue-key per frame" in mdef


def test_default_builds_all_three_modes(tmp_path):
    modes, mdef = _modes_and_def(_build(tmp_path))
    assert modes == ["mode1", "mode2", "mode3"]
    assert "mode3=hydrophobic dedup per frame, others count all rows" in mdef


def test_mode_definition_lists_only_built_modes(tmp_path):
    _, mdef = _modes_and_def(_build(tmp_path, mode2=False, mode3=False))
    assert mdef == "mode1=count all rows"


def test_hbonds_mode1_and_mode2_differ(tmp_path):
    # Sanity: the two kept modes are genuinely different (mode1 counts the
    # duplicate ARG4-THR118 pair twice per frame, mode2 once).
    path = _build(tmp_path, mode3=False)
    f = h5py.File(path, "r")
    m1 = {bytes(x).decode() for x in f["hbonds/modes/mode1/table"][:]}
    m2 = {bytes(x).decode() for x in f["hbonds/modes/mode2/table"][:]}
    assert m1 != m2
