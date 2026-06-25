"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Mock PTA builders for the plot-module test suite.

Each builder produces a fully-valid PTA file (signed, tokenized, completed
groups) that passes `command_line/plot/pta.py::validate()`. Builders write
to a caller-supplied path so tests can use `tmp_path` and let pytest clean
up automatically.

Schemas mirrored:

    PLI :  /pl_interactions/modes/mode{1,2,3}/table        JSON-per-row
                                                           {key:[r1,r2,label], count, frequency}
           /pl_interactions/frame_0/interactions           empty schema-stable dataset
    PPI :  /pp_interactions/...                            same shape, different group
    PCA :  /pca/frame_<i>/pca                              {pc, value, variance_ratio}
    RMSD:  /rmsd/frame_<i>/rmsd                            {label, value}
    ANG :  /angles/frame_<i>/angles                        {label, kind, value}
    DIST:  /distances/frame_<i>/distances                  {label, method, distance}
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List, Mapping, Sequence

import h5py
import numpy as np

from pharmacon.fileio.pta import PharmaconPTAFile
from pharmacon.utils.identifiers import create_mda_artifact_token


__all__ = [
    "build_pli_pta",
    "build_pli_frames_pta",
    "build_pli_merged_pta",
    "build_ppi_pta",
    "build_ppi_frames_pta",
    "build_ppi_merged_pta",
    "build_pca_pta",
    "build_universal_pta",
    "build_rmsf_pta",
    "build_rmsf_merged_pta",
    "DEFAULT_PLI_ROWS",
    "DEFAULT_PPI_ROWS",
    "DEFAULT_PLI_FRAME_SPECS",
    "DEFAULT_PPI_FRAME_SPECS",
    "DEFAULT_RMSF_SELECTIONS",
]


# Per-frame interaction-record specs. Each spec describes one
# protein↔partner contact and the fraction of frames it is present in.
# A record is the 20-field, frame_number-less row that
# ``parse_frame_interaction_record`` consumes:
#
#   [label,
#    a1_idx, a1_name, a1_id, a1_element, a1_type,   a1_resn, a1_resid, a1_chain, a1_seg,
#    a2_idx, a2_name, a2_id, a2_element, a2_type,   a2_resn, a2_resid, a2_chain, a2_seg,
#    details]
#
# atom1 is the protein (a1_type "BB"/"SC" drives the backbone/side-chain
# split); atom2 is the ligand (PLI) or the partner protein (PPI), and its
# name is what the ligand-atom monitor buckets on.
#
# Spec tuple:
#   (label, a1_type, a1_resn, a1_resid, a1_chain, a1_seg,
#    a2_name, a2_type, a2_resn, a2_resid, a2_chain, a2_seg, frequency)
DEFAULT_PLI_FRAME_SPECS: tuple = (
    ("HYDROPHOBIC",  "SC", "PHE", 45,  "A", "PROA", "C1", "SC", "LIG", 1, "A", "LIGA", 0.80),
    ("HYDROPHOBIC",  "BB", "ILE", 175, "A", "PROA", "C2", "SC", "LIG", 1, "A", "LIGA", 0.50),
    ("HYDROGEN-BOND", "SC", "ASP", 102, "A", "PROA", "N4", "SC", "LIG", 1, "A", "LIGA", 0.60),
    ("HYDROGEN-BOND", "BB", "ASN", 184, "A", "PROA", "O2", "SC", "LIG", 1, "A", "LIGA", 0.40),
    ("IONIC",        "SC", "ARG", 150, "A", "PROA", "N1", "SC", "LIG", 1, "A", "LIGA", 0.30),
)

DEFAULT_PPI_FRAME_SPECS: tuple = (
    ("HYDROPHOBIC",  "SC", "PHE", 45,  "A", "PROA", "CA", "SC", "VAL", 200, "B", "PROB", 0.50),
    ("HYDROGEN-BOND", "SC", "ASP", 102, "A", "PROA", "CA", "SC", "LYS", 215, "B", "PROB", 0.60),
    ("IONIC",        "SC", "ARG", 150, "A", "PROA", "CA", "SC", "GLU", 220, "B", "PROB", 0.30),
)


def _spec_to_record(spec: Sequence) -> list:
    """Turn a frame spec into a frame_number-less interaction record."""
    (label, a1_type, a1_rn, a1_rid, a1_ch, a1_seg,
     a2_name, a2_type, a2_rn, a2_rid, a2_ch, a2_seg, _freq) = spec
    return [
        str(label),
        1, "CA", 1, "C", str(a1_type),
        str(a1_rn), int(a1_rid), str(a1_ch), str(a1_seg),
        2, str(a2_name), 2, "C", str(a2_type),
        str(a2_rn), int(a2_rid), str(a2_ch), str(a2_seg),
        {"distance": 3.5},
    ]


# A small, deterministic PLI dataset used by tests. Each tuple is
# (interaction, protein_resname, protein_resid, protein_chainid,
#  ligand_resname, ligand_resid, ligand_chainid, frequency).
DEFAULT_PLI_ROWS: tuple = (
    ("HYDROPHOBIC",     "PHE", 45,  "A", "LIG", 1, "A", 0.80),
    ("HYDROPHOBIC",     "ILE", 175, "A", "LIG", 1, "A", 0.27),
    ("HYDROGEN-BOND",   "ASP", 102, "A", "LIG", 1, "A", 0.55),
    ("HYDROGEN-BOND",   "ASN", 184, "A", "LIG", 1, "A", 0.36),
    ("IONIC",           "ARG", 150, "A", "LIG", 1, "A", 0.41),
    ("HALOGEN-BOND",    "TYR", 77,  "A", "LIG", 1, "A", 0.12),
    ("PI-STACKING",     "TRP", 32,  "A", "LIG", 1, "A", 0.66),
    ("PI-CATION",       "GLU", 68,  "A", "LIG", 1, "A", 0.28),
    ("METAL-CONTACT",   "HIS", 90,  "A", "LIG", 1, "A", 0.07),
    ("WATER-BRIDGE-1",  "GLN", 92,  "A", "LIG", 1, "A", 0.30),
)


# A small PPI dataset; protein↔protein pairs. Same column layout.
DEFAULT_PPI_ROWS: tuple = (
    ("HYDROPHOBIC",     "PHE", 45,  "A", "VAL", 200, "B", 0.40),
    ("HYDROGEN-BOND",   "ASP", 102, "A", "LYS", 215, "B", 0.55),
    ("IONIC",           "ARG", 150, "A", "GLU", 220, "B", 0.22),
    ("PI-STACKING",     "TYR", 77,  "A", "TRP", 230, "B", 0.18),
)


def _override_version(path: Path, version: str) -> None:
    """Re-sign a freshly-built PTA with a non-default pharmacon_version.

    PharmaconPTAFile always stamps the current runtime version at creation;
    use this only when a test explicitly wants to simulate a file produced
    by another release. Updates the version attr AND regenerates the
    signature/fingerprint under the new version so validation can succeed
    on the cooked file.
    """
    import h5py
    from pharmacon.utils.fingerprint import create_pharmacon_signature

    fmt = path.suffix.lstrip(".").lower() or "pta"
    with h5py.File(path, "a") as fh:
        attrs = fh.attrs
        attrs["pharmacon_version"] = version
        command = str(attrs.get("command", "")).strip()
        subcommand = str(attrs.get("subcommand", "")).strip()
        if command and subcommand:
            sig = create_pharmacon_signature(
                format_name=fmt,
                command=command, subcommand=subcommand, version=version,
            )
            attrs["signature"] = sig.signature
            attrs["fingerprint"] = sig.fingerprint


def _write_interaction_pta(*, path: Path, group: str, subcommand: str,
                           rows: Sequence[Sequence], n_frames: int,
                           description: str,
                           frame_specs: Sequence | None = None,
                           pharmacon_version: str | None = None) -> Path:
    """Internal: build a PLI- or PPI-shaped PTA. Used by both builders.

    When ``frame_specs`` is given, real per-frame ``frame_<N>/interactions``
    records are written (one row per active spec per frame) so plots that
    walk per-frame data (heatmaps, pie charts, ligand monitor, PPI timeline)
    have something to render. Otherwise a single empty frame_0 is written.
    """

    records: list[dict] = []
    for r in rows:
        label, p_rn, p_rid, p_ch, l_rn, l_rid, l_ch, freq = r
        key = [
            [str(p_rn), int(p_rid), str(p_ch or ""), ""],
            [str(l_rn), int(l_rid), str(l_ch or ""), ""],
            str(label),
        ]
        count = int(round(float(freq) * n_frames))
        records.append({"key": key, "count": count, "frequency": float(freq)})

    begin, end, step = 0, max(0, n_frames - 1), 1
    blueprint = f"mock-{subcommand}::{path.name}::{n_frames}::{len(records)}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand=subcommand) as pta:
        pta.add_file_metadata({
            "description": description,
            "topology_file": "mock://topology",
            "trajectory_file": "mock://trajectory",
            "begin": str(begin), "end": str(end), "step": str(step),
            "total_frames": str(n_frames),
            "is_merged": "False",
            "workers": "1",
            "labels": subcommand,
            "reference_frame": "0",
        })

        pta.create_group(group)
        pta.add_group_metadata(group_name=group, metadata={"completed": "True"})

        if frame_specs:
            # Real per-frame records: each spec is present in the first
            # round(frequency * n_frames) frames.
            prepared = [
                (_spec_to_record(s), int(round(float(s[-1]) * n_frames)))
                for s in frame_specs
            ]
            for fidx in range(max(1, n_frames)):
                active = [rec for rec, until in prepared if fidx < until]
                pta.write_frame_interactions(
                    frame_index=fidx, interactions=active,
                    group_name=group, overwrite=True,
                )
        else:
            # Schema-stable empty frame_0/interactions so consumers that iterate
            # frame_* under the group find at least one entry.
            pta.write_frame_interactions(
                frame_index=0, interactions=[], group_name=group, overwrite=True,
            )

        modes_root = f"{group}/modes"
        pta.create_group(modes_root)

        for mode_name, description in (
            ("mode1", "Count frequency interactions analysis data (residue-level, count all occurrences)"),
            ("mode2", "Residue-level once-per-frame interactions analysis data (deduplicate within frame)"),
            ("mode3", "Hybrid interactions analysis data (hydrophobic once-per-frame; others count all occurrences)"),
        ):
            if records:
                data = np.array(
                    [json.dumps(rec) for rec in records],
                    dtype=h5py.string_dtype("utf-8"),
                )
            else:
                data = np.empty((0,), dtype=h5py.string_dtype("utf-8"))

            mode_group = f"{modes_root}/{mode_name}"
            pta.create_group(mode_group)
            pta.create_dataset(
                group_name=mode_group, dataset_name="table", data=data,
                metadata={
                    "n_frames": str(n_frames),
                    "n_rows": str(len(data)),
                    "frame_begin": str(begin),
                    "frame_end": str(end),
                    "frame_step": str(step),
                    "end_inclusive": "True",
                    "normalization": "frequency = count / n_frames",
                    "format": "json-per-row",
                    "mode_definition": (
                        "mode1=count all rows; "
                        "mode2=deduplicate residue-key per frame; "
                        "mode3=hydrophobic dedup per frame, others count all rows"
                    ),
                },
            )
            pta.add_group_metadata(
                group_name=mode_group,
                metadata={"description": description, "mode": mode_name, "completed": "True"},
                overwrite=True,
            )

        pta.add_group_metadata(
            group_name=modes_root,
            metadata={
                "completed": "True",
                "n_frames": str(n_frames),
                "frame_begin": str(begin),
                "frame_end": str(end),
                "frame_step": str(step),
                "end_inclusive": "True",
            },
            overwrite=True,
        )

        token = create_mda_artifact_token(
            blueprint=blueprint, secret="trajectory_analysis", namespace="pharmacon",
        )
        pta.add_file_metadata(
            {
                "completed": "True",
                "artifact_status": "SUCCESS",
                "artifact_status_code": "0",
                "blueprint": blueprint,
                "artifact_token": token,
                "artifact_token_version": "1",
            },
            overwrite=True,
        )

    if pharmacon_version is not None:
        _override_version(path, pharmacon_version)

    return path


def build_pli_pta(path: Path, *, rows: Sequence = DEFAULT_PLI_ROWS,
                  n_frames: int = 100,
                  pharmacon_version: str | None = None) -> Path:
    """Build a mock PL-interactions PTA at `path` and return the path."""
    return _write_interaction_pta(
        path=path, group="pl_interactions", subcommand="pl-interactions",
        rows=rows, n_frames=n_frames,
        description="Protein-Ligand Interaction Analysis (mock)",
        pharmacon_version=pharmacon_version,
    )


def build_ppi_pta(path: Path, *, rows: Sequence = DEFAULT_PPI_ROWS,
                  n_frames: int = 100,
                  pharmacon_version: str | None = None) -> Path:
    """Build a mock PP-interactions PTA at `path` and return the path."""
    return _write_interaction_pta(
        path=path, group="pp_interactions", subcommand="pp-interactions",
        rows=rows, n_frames=n_frames,
        description="Protein-Protein Interaction Analysis (mock)",
        pharmacon_version=pharmacon_version,
    )


def build_pli_frames_pta(path: Path, *,
                         frame_specs: Sequence = DEFAULT_PLI_FRAME_SPECS,
                         rows: Sequence = DEFAULT_PLI_ROWS,
                         n_frames: int = 100) -> Path:
    """Build a PL-interactions PTA with real per-frame interaction records.

    Suitable for the plots that walk frame_<N>/interactions: stacked-column-2
    (backbone/side-chain), heatmaps 1 & 2, pie charts, and the ligand-atom
    monitor. Also carries the modes tables, so stacked-column-1 still works.
    """
    return _write_interaction_pta(
        path=path, group="pl_interactions", subcommand="pl-interactions",
        rows=rows, n_frames=n_frames,
        description="Protein-Ligand Interaction Analysis, per-frame (mock)",
        frame_specs=frame_specs,
    )


def build_ppi_frames_pta(path: Path, *,
                         frame_specs: Sequence = DEFAULT_PPI_FRAME_SPECS,
                         rows: Sequence = DEFAULT_PPI_ROWS,
                         n_frames: int = 100) -> Path:
    """Build a PP-interactions PTA with real per-frame interaction records.

    Needed by the PPI timeline-pairs plot, which builds a residue-pair ×
    frame matrix from frame_<N>/interactions.
    """
    return _write_interaction_pta(
        path=path, group="pp_interactions", subcommand="pp-interactions",
        rows=rows, n_frames=n_frames,
        description="Protein-Protein Interaction Analysis, per-frame (mock)",
        frame_specs=frame_specs,
    )


def build_hbonds_frames_pta(path: Path, *,
                            frame_specs: Sequence = DEFAULT_PPI_FRAME_SPECS,
                            rows: Sequence = DEFAULT_PPI_ROWS,
                            n_frames: int = 100) -> Path:
    """Build an h-bonds PTA (group ``hbonds``) with per-frame interaction
    records + modes, for the h-bond plots (heatmap, timeline, count, occupancy).
    """
    return _write_interaction_pta(
        path=path, group="hbonds", subcommand="h-bonds",
        rows=rows, n_frames=n_frames,
        description="Hydrogen Bond Analysis, per-frame (mock)",
        frame_specs=frame_specs,
    )


def _write_merged_interaction_pta(*, path: Path, group: str, subcommand: str,
                                  rows: Sequence[Sequence], n_files: int,
                                  description: str) -> Path:
    """Internal: build a merged-shape PTA (modes_merged/...)."""

    records: list[dict] = []
    for r in rows:
        label, p_rn, p_rid, p_ch, l_rn, l_rid, l_ch, mean_freq, std_freq = r
        key = [
            [str(p_rn), int(p_rid), str(p_ch or ""), ""],
            [str(l_rn), int(l_rid), str(l_ch or ""), ""],
            str(label),
        ]
        records.append({
            "key": key,
            "mean_frequency": float(mean_freq),
            "std_frequency": float(std_freq),
            "n_files": int(n_files),
        })

    blueprint = f"mock-merged-{subcommand}::{path.name}::{n_files}::{len(records)}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Merge Results", subcommand=subcommand) as pta:
        pta.add_file_metadata({
            "description": description,
            "topology_file": "mock://topology",
            "trajectory_file": "mock://trajectory",
            "is_merged": "True",
            "n_files": str(n_files),
            "labels": subcommand,
            "reference_frame": "0",
        })

        pta.create_group(group)
        pta.add_group_metadata(group_name=group, metadata={"completed": "True"})

        merged_root = f"{group}/modes_merged"
        pta.create_group(merged_root)

        for mode_name in ("mode1", "mode2", "mode3"):
            if records:
                data = np.array(
                    [json.dumps(rec) for rec in records],
                    dtype=h5py.string_dtype("utf-8"),
                )
            else:
                data = np.empty((0,), dtype=h5py.string_dtype("utf-8"))

            mode_group = f"{merged_root}/{mode_name}"
            pta.create_group(mode_group)
            pta.create_dataset(
                group_name=mode_group, dataset_name="table", data=data,
                metadata={
                    "n_files": str(n_files),
                    "n_rows": str(len(data)),
                    "normalization": "mean_frequency, std_frequency across files",
                    "format": "json-per-row",
                },
            )
            pta.add_group_metadata(
                group_name=mode_group,
                metadata={"mode": mode_name, "completed": "True"},
                overwrite=True,
            )

        pta.add_group_metadata(
            group_name=merged_root,
            metadata={"completed": "True"},
            overwrite=True,
        )

        _finalize_metadata(pta, blueprint=blueprint)

    return path


def _default_merged_rows(rows: Sequence) -> list:
    """Convert (..., frequency) rows to (..., mean_frequency, std_frequency)."""
    out = []
    for r in rows:
        label, p_rn, p_rid, p_ch, l_rn, l_rid, l_ch, freq = r
        out.append((label, p_rn, p_rid, p_ch, l_rn, l_rid, l_ch,
                    float(freq), float(freq) * 0.1))
    return out


def build_pli_merged_pta(path: Path, *, rows: Sequence | None = None,
                         n_files: int = 3) -> Path:
    """Build a merged-mode PL-interactions PTA (writes modes_merged/...)."""
    if rows is None:
        rows = _default_merged_rows(DEFAULT_PLI_ROWS)
    return _write_merged_interaction_pta(
        path=path, group="pl_interactions", subcommand="pl-interactions",
        rows=rows, n_files=n_files,
        description="Merged Protein-Ligand Interaction Analysis (mock)",
    )


def build_ppi_merged_pta(path: Path, *, rows: Sequence | None = None,
                         n_files: int = 3) -> Path:
    """Build a merged-mode PP-interactions PTA (writes modes_merged/...)."""
    if rows is None:
        rows = _default_merged_rows(DEFAULT_PPI_ROWS)
    return _write_merged_interaction_pta(
        path=path, group="pp_interactions", subcommand="pp-interactions",
        rows=rows, n_files=n_files,
        description="Merged Protein-Protein Interaction Analysis (mock)",
    )


def _finalize_metadata(pta: PharmaconPTAFile, *, blueprint: str) -> None:
    token = create_mda_artifact_token(
        blueprint=blueprint, secret="trajectory_analysis", namespace="pharmacon",
    )
    pta.add_file_metadata(
        {
            "completed": "True",
            "artifact_status": "SUCCESS",
            "artifact_status_code": "0",
            "blueprint": blueprint,
            "artifact_token": token,
            "artifact_token_version": "1",
        },
        overwrite=True,
    )


def build_pca_pta(path: Path, *, n_frames: int = 10, n_components: int = 4,
                  variance_ratios: Sequence[float] | None = None,
                  group: str = "pca") -> Path:
    """Build a mock PCA PTA at `path`. Each frame stores n_components rows.

    ``group`` defaults to ``"pca"`` but the real analyzer writes one group per
    selection named ``pca_<selection>`` (e.g. ``pca_protein``); pass that to
    exercise the production group-naming convention.
    """

    if variance_ratios is None:
        # geometric decay so the variance ratio array sums close to 1
        raw = [0.5 ** i for i in range(n_components)]
        s = sum(raw)
        variance_ratios = [v / s for v in raw]
    assert len(variance_ratios) == n_components, "variance_ratios length must match n_components"

    rng = np.random.default_rng(0)
    blueprint = f"mock-pca::{path.name}::{n_frames}::{n_components}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand="pca") as pta:
        pta.add_file_metadata({
            "description": "Principal Component Analysis (mock)",
            "topology_file": "mock://topology",
            "trajectory_file": "mock://trajectory",
            "begin": "0", "end": str(max(0, n_frames - 1)), "step": "1",
            "total_frames": str(n_frames),
            "is_merged": "False",
            "labels": "pca",
            "reference_frame": "0",
        })

        pta.create_group(group)
        pta.add_group_metadata(group_name=group, metadata={"completed": "True"})

        for frame in range(n_frames):
            frame_group = f"{group}/frame_{frame}"
            pta.create_group(frame_group)
            # Reproducible-but-varied projections
            samples = rng.normal(size=n_components).tolist()
            rows = [
                {"pc": pc, "value": float(samples[pc - 1]),
                 "variance_ratio": float(variance_ratios[pc - 1])}
                for pc in range(1, n_components + 1)
            ]
            data = np.array(
                [json.dumps(r) for r in rows],
                dtype=h5py.string_dtype("utf-8"),
            )
            pta.create_dataset(
                group_name=frame_group, dataset_name="pca", data=data,
                metadata={"frame_index": str(frame),
                          "n_rows": str(len(data)),
                          "time_ps": str(float(frame) * 10.0),
                          "format": "json-per-row"},
            )

        _finalize_metadata(pta, blueprint=blueprint)

    return path


def build_universal_pta(path: Path, *, group: str,
                        series: Mapping[str, Sequence[float] | Mapping[str, Sequence[float]]],
                        n_frames: int | None = None,
                        time_ps_step: float = 10.0,
                        units: str | None = None) -> Path:
    """Build a mock universal-timeseries PTA (rmsd / angles / distances).

    `series` shape depends on `group`:
        rmsd      : {label: [values per frame]}
        angles    : {label: {kind: [values per frame]}}
        distances : {label: {method: [values per frame]}}
    """

    valid = {"rmsd", "angles", "distances"}
    if group not in valid:
        raise ValueError(f"group must be one of {valid}, got {group!r}")

    # Determine n_frames from any series if not provided.
    if n_frames is None:
        try:
            if group == "rmsd":
                n_frames = len(next(iter(series.values())))
            else:
                # nested
                first_label = next(iter(series.values()))
                n_frames = len(next(iter(first_label.values())))
        except (StopIteration, TypeError):
            raise ValueError("series is empty; specify n_frames explicitly")

    subcommand = {"rmsd": "rmsd", "angles": "angles", "distances": "distances"}[group]
    blueprint = f"mock-{subcommand}::{path.name}::{n_frames}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand=subcommand) as pta:
        pta.add_file_metadata({
            "description": f"{group} timeseries (mock)",
            "topology_file": "mock://topology",
            "trajectory_file": "mock://trajectory",
            "begin": "0", "end": str(max(0, n_frames - 1)), "step": "1",
            "total_frames": str(n_frames),
            "is_merged": "False",
            "labels": subcommand,
            "reference_frame": "0",
        })

        pta.create_group(group)
        pta.add_group_metadata(group_name=group, metadata={"completed": "True"})

        for frame in range(n_frames):
            frame_group = f"{group}/frame_{frame}"
            pta.create_group(frame_group)

            rows: list[dict] = []
            for label, payload in series.items():
                if group == "rmsd":
                    rows.append({"label": str(label), "value": float(payload[frame])})
                elif group == "angles":
                    for kind, values in payload.items():
                        rows.append({"label": str(label), "kind": str(kind),
                                     "value": float(values[frame])})
                else:  # distances
                    for method, values in payload.items():
                        rows.append({"label": str(label), "method": str(method),
                                     "distance": float(values[frame])})

            data = np.array(
                [json.dumps(r) for r in rows],
                dtype=h5py.string_dtype("utf-8"),
            ) if rows else np.empty((0,), dtype=h5py.string_dtype("utf-8"))

            frame_metadata = {"frame_index": str(frame),
                              "n_rows": str(len(data)),
                              "time_ps": str(float(frame) * time_ps_step),
                              "format": "json-per-row"}
            if units is not None:
                frame_metadata["units"] = str(units)

            pta.create_dataset(
                group_name=frame_group, dataset_name=group, data=data,
                metadata=frame_metadata,
            )

        _finalize_metadata(pta, blueprint=blueprint)

    return path


# A small, deterministic RMSF dataset shared by the RMSF tests.
# Each entry is (label, selection_string, atoms) where atoms is a list of
# (atom_index, resid, resname, atom_name, rmsf) tuples.
DEFAULT_RMSF_SELECTIONS: tuple = (
    ("calpha", "name CA", [
        (0, 1, "ALA", "CA", 0.50),
        (4, 2, "GLY", "CA", 1.20),
        (9, 3, "SER", "CA", 0.80),
    ]),
    ("backbone", "backbone", [
        (0, 1, "ALA", "N",  0.30),
        (1, 1, "ALA", "CA", 0.50),
        (2, 1, "ALA", "C",  0.40),
        (3, 2, "GLY", "N",  0.90),
        (4, 2, "GLY", "CA", 1.20),
        (5, 2, "GLY", "C",  1.10),
    ]),
)


def _selections_to_dicts(selections: Sequence) -> list[dict]:
    """Convert (label, sel_str, [(i, resid, resname, atom_name, rmsf), ...])
    tuples into the dict-list shape `write_rmsf_data` consumes."""
    out: list[dict] = []
    for label, sel_str, atoms in selections:
        out.append({
            "label": label,
            "selection_string": sel_str,
            "atom_indices": np.array([a[0] for a in atoms], dtype=int),
            "resids":       np.array([a[1] for a in atoms], dtype=int),
            "resnames":     np.array([a[2] for a in atoms]),
            "atom_names":   np.array([a[3] for a in atoms]),
            "rmsf":         np.array([a[4] for a in atoms], dtype=float),
        })
    return out


def build_rmsf_pta(
    path: Path,
    *,
    selections: Sequence | None = None,
    fitting_group: str = "protein and name CA",
    frame_begin: int = 0,
    frame_end: int = 100,
    frame_step: int = 1,
) -> Path:
    """Build a mock non-merged RMSF PTA at `path`.

    Uses the production `write_rmsf_data` + `write_rmsf_statistics` writers
    so the file shape matches what `trajectory rmsf` actually produces.
    """
    if selections is None:
        selections = DEFAULT_RMSF_SELECTIONS

    sel_dicts = _selections_to_dicts(selections)
    blueprint = f"mock-rmsf::{path.name}::{len(sel_dicts)}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand="rmsf") as pta:
        pta.add_file_metadata({
            "description": "RMSF (mock)",
            "topology_file": "mock://topology",
            "trajectory_file": "mock://trajectory",
            "fitting_group": fitting_group,
            "selections": str([s[1] for s in selections]),
            "names": str([s[0] for s in selections]),
            "reference_frame": "average",
            "begin": str(frame_begin),
            "end": str(frame_end),
            "step": str(frame_step),
            "is_merged": "False",
            "total_frames": str(frame_end),
        })
        pta.create_group("rmsf")
        pta.write_rmsf_data(
            selections=sel_dicts,
            group_name="rmsf",
            frame_begin=frame_begin,
            frame_end=frame_end,
            frame_step=frame_step,
            fitting_group=fitting_group,
        )
        pta.write_rmsf_statistics(group_name="rmsf")
        pta.add_group_metadata(group_name="rmsf", metadata={"completed": "True"})
        _finalize_metadata(pta, blueprint=blueprint)

    return path


def build_rmsf_merged_pta(
    path: Path,
    *,
    selections: Sequence | None = None,
    n_inputs: int = 2,
    std_value: float = 0.1,
    fitting_group: str = "protein and name CA",
    frame_begin: int = 0,
    frame_end: int = 100,
    frame_step: int = 1,
    compute_statistics: bool = False,
) -> Path:
    """Build a mock merged RMSF PTA. Per-atom records carry
    {atom_index, resid, resname, atom_name, mean, std, n} instead of the
    non-merged {..., rmsf} schema.

    Writes per-selection subgroups directly via `create_dataset` to avoid
    going through the production writer (which expects "rmsf"). When
    `compute_statistics=True`, also calls `write_rmsf_statistics` inside
    the same write session (the real `merge results` does NOT, per
    decision 2a; this knob exists for tests that exercise the merged-
    record fallback in the statistics builder).
    """
    if selections is None:
        selections = DEFAULT_RMSF_SELECTIONS

    blueprint = f"mock-rmsf-merged::{path.name}::{len(selections)}"

    with PharmaconPTAFile(path, overwrite=True, mode="a",
                          command="Trajectory Analysis", subcommand="rmsf") as pta:
        pta.add_file_metadata({
            "description": "Merged RMSF (mock)",
            "fitting_group": fitting_group,
            "is_merged": "True",
            "n_inputs": str(n_inputs),
            "merge_strategy": "per-selection avg+std",
        })
        pta.create_group("rmsf")

        for label, sel_str, atoms in selections:
            sel_group = f"rmsf/selection_{label}"
            pta.create_group(sel_group)

            records = [
                {
                    "atom_index": int(i),
                    "resid":      int(resid),
                    "resname":    str(resname),
                    "atom_name":  str(atom_name),
                    "mean":       float(rmsf),
                    "std":        float(std_value),
                    "n":          int(n_inputs),
                }
                for (i, resid, resname, atom_name, rmsf) in atoms
            ]
            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
            pta.create_dataset(
                group_name=sel_group,
                dataset_name="atoms",
                data=data,
                metadata={
                    "label": label,
                    "selection_string": sel_str,
                    "n_atoms": str(len(atoms)),
                    "fitting_group": fitting_group,
                    "frame_begin": str(frame_begin),
                    "frame_end":   str(frame_end),
                    "frame_step":  str(frame_step),
                    "merged":      "True",
                    "n_inputs":    str(n_inputs),
                    "schema":      "{atom_index, resid, resname, atom_name, mean, std, n}",
                    "format":      "json-per-row",
                },
            )

        if compute_statistics:
            pta.write_rmsf_statistics(group_name="rmsf")

        pta.add_group_metadata(group_name="rmsf", metadata={"completed": "True"})
        _finalize_metadata(pta, blueprint=blueprint)

    return path
