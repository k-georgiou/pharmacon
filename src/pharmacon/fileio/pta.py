"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Handles PharmaconPTA specific HDF5 file operations.

This module contains the `PharmaconPTAFile` class, which is a specialized
implementation for handling PharmaconPTA-related HDF5 file operations. It is
designed to handle the specific data formats and functionalities required for
PharmaconPTA, while leveraging the base `PharmaconHDF5File` class. The class
supports initialization with overwrite options and specific file modes.

Classes:
    PharmaconPTAFile: Manages PharmaconPTA-specific HDF5 file operations.
"""


import csv
import json
import h5py
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Union, Iterable, Tuple, Final, List, Dict, Counter

from .base import PharmaconHDF5File
from pharmacon.logger import get_logger, PharmaconLogger



__all__ = [
    "logger",
    "PharmaconPTAFile",
]


logger: PharmaconLogger = get_logger(__name__)


class PharmaconPTAFile(PharmaconHDF5File):
    """
    Handles PharmaconPTA specific HDF5 file operations.

    This class is a specialized implementation for handling PharmaconPTA-related
    HDF5 file operations. It extends the base `PharmaconHDF5File`, introducing
    modifications or specific functionality required to work with PharmaconPTA
    data formats. It allows initialization with various options, including
    overwrite capability and HDF5 file modes.

    :ivar path: The file path of the HDF5 file.
    :type path: Union[str, Path]
    :ivar overwrite: Whether to overwrite the file if it exists.
    :type overwrite: bool
    :ivar mode: The file mode for accessing the HDF5 file. Common options are
        read, write, or append.
    :type mode: str
    """

    INTERACTION_COLUMNS: Final[List[str]] = [
        "frame_number",
        "interaction",

        # atom 1
        "atom1_index",
        "atom1_name",
        "atom1_id",
        "atom1_type",
        "atom1_element",
        "atom1_resname",
        "atom1_resid",
        "atom1_chainid",
        "atom1_segid",

        # atom 2
        "atom2_index",
        "atom2_name",
        "atom2_id",
        "atom2_type",
        "atom2_element",
        "atom2_resname",
        "atom2_resid",
        "atom2_chainid",
        "atom2_segid",

        "details",
    ]
    INTERACTION_DETAIL_SCHEMA: Final[Dict[str, Tuple[str, ...]]] = {

        "HYDROPHOBIC": ("distance",
                        "is_hydrogen",),

        "HYDROGEN-BOND": ("distance",
            "angle_dha",
            "angle_hax",
            "orientation",),

        "IONIC": ("distance",
                  "orientation",),

        "HALOGEN-BOND": ("distance",
                         "angle_cxa",
                         "angle_xay",
                         "orientation",),

        "METAL-CONTACT": ("distance",
                          "role",
                          "orientation"),

        "WATER-BRIDGE-1": ("water_index",
                           "water_name",
                           "water_id",
                           "water_element",
                           "water_type",
                           "water_resname",
                           "water_resid",
                           "water_chainid",
                           "water_segid",
                           "d_g1_water",
                           "d_g2_water",
                           "angle_1_dha",
                           "angle_2_dha",
                           "angle_1_hax",
                           "angle_2_hax",
                           "orientation"),

        "WATER-BRIDGE-2": ("water1_index",
                           "water1_name",
                           "water1_id",
                           "water1_element",
                           "water1_type",
                           "water1_resname",
                           "water1_resid",
                           "water1_chainid",
                           "water1_segid",
                           "water2_index",
                           "water2_name",
                           "water2_id",
                           "water2_element",
                           "water2_type",
                           "water2_resname",
                           "water2_resid",
                           "water2_chainid",
                           "water2_segid",
                           "d_g1_w1",
                           "d_w1_w2",
                           "d_w2_g2",
                           "angle_1_dha",
                           "angle_ww_dha",
                           "angle_2_dha",
                           "angle_1_hax",
                           "angle_2_hax",
                           "orientation"),

        "PI-CATION": ("distance",
                      "theta",
                      "orientation"),

        "PI-STACKING": ("distance",
                        "angle",
                        "stacking_type",
                        "orientation"),
    }

    def __init__(self, path: Union[str, Path], *, overwrite: bool = False, mode: str = "a",
                 command: str = "", subcommand: str = "") -> None:
        super().__init__(path=path, overwrite=overwrite, mode=mode,
                         command=command, subcommand=subcommand)


    def write_frame_interactions(self, *, frame_index: int, interactions: Iterable[Tuple],
                                 group_name: str, overwrite: bool = False) -> None:
        """
        Writes interaction data for a specific frame into a hierarchical data structure.

        Interactions are stored as JSON-per-row inside an HDF5 dataset:

            /<group_name>/frame_<frame_index>/interactions

        This method is **frame-atomic**: each call fully defines one frame.

        This method organizes interaction records for a given frame into a
        dedicated group, serialized in a stable JSON-per-row format. If the
        specified group for the frame already exists, you can optionally
        overwrite it, provided the `overwrite` flag is set to True.

        :param frame_index: The numerical index of the frame. Must be a non-negative integer.
        :type frame_index: int
        :param interactions: Iterable collection of tuples representing interaction records.
        :param group_name: Name of the parent group under which the frame interactions will be stored.
        :type group_name: str
        :param overwrite: Boolean flag specifying whether to overwrite existing frame data.
        :type overwrite: bool
        :return: None
        :rtype: NoneType

        :raises ValueError: If the provided `frame_index` is not a non-negative integer.
        :raises FileExistsError: If the frame group already exists under the specified parent group and `overwrite` is False.
        """
        # --- validate frame index ---
        if not isinstance(frame_index, int) or frame_index < 0:
            raise ValueError("frame_index must be a non-negative integer")

        frame_group = f"{group_name}/frame_{frame_index}"

        # --- overwrite semantics ---
        if self.group_exists(frame_group):
            if not overwrite:
                raise FileExistsError(
                    f"Frame {frame_index} already exists under '{group_name}'"
                )
            self.delete_group(frame_group)

        self.create_group(frame_group)

        # --- materialize interactions (small, per-frame) ---
        records = list(interactions)

        # --- serialize JSON-per-row ---
        if records:
            data = np.array(
                [json.dumps(rec) for rec in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
        else:
            # empty but schema-stable dataset
            data = np.empty(
                (0,),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

        # --- write dataset + metadata ---
        self.create_dataset(
            group_name=frame_group,
            dataset_name="interactions",
            data=data,
            metadata={
                "frame_index": str(frame_index),
                "n_interactions": str(len(data)),
                "format": "json-per-row",
            },
        )

    def read_frame_interactions(self, *, frame_index: int, group_name: str) -> Tuple[Tuple, ...]:
        """
        Read all interaction records for a specific frame.

        Interactions are read from:

            /<group_name>/frame_<frame_index>/interactions

        Records are deserialized from JSON-per-row and returned as tuples.

        :param frame_index: Frame index to read.
        :type frame_index: int
        :param group_name: Parent group name.
        :type group_name: str
        :return: Tuple of interaction records.
        :rtype: Tuple[Tuple, ...]
        :raises ValueError: If frame_index is invalid.
        :raises KeyError: If the frame group or dataset does not exist.
        """

        if not isinstance(frame_index, int) or frame_index < 0:
            raise ValueError("frame_index must be a non-negative integer")

        dset_path = f"{group_name}/frame_{frame_index}/interactions"

        if dset_path not in self.file:
            raise KeyError(f"No interactions found for frame {frame_index}")

        dset = self.file[dset_path]

        if dset.size == 0:
            return ()

        # rows are:
        # ['PI-CATION', 'ring', '682,683,685,687,689,691', 'CG,CD1,CE1,CZ,CE2,CD2', '683,684,686,688,690,692',
        # 'C,C,C,C,C,C', 'SC,SC,SC,SC,SC,SC', 'PHE,PHE,PHE,PHE,PHE,PHE', '45,45,45,45,45,45', '',
        # 'SYSTEM,SYSTEM,SYSTEM,SYSTEM,SYSTEM,SYSTEM', 1602, 'NH1', 1603, 'N', 'SC', 'ARG', 105, '',
        # 'SYSTEM', 5.770456248419717, 24.10333551636779, 'cation-above-face']
        #
        # ['PI-STACKING', '474,475,477,479,488', 'CG,CD1,NE1,CE2,CD2', '475,476,478,480,489', 'C,C,N,C,C',
        # 'SC,SC,SC,SC,SC', 'TRP,TRP,TRP,TRP,TRP', '32,32,32,32,32', '', 'SYSTEM,SYSTEM,SYSTEM,SYSTEM,SYSTEM',
        # '480,479,488,486,484,482', 'CZ2,CE2,CD2,CE3,CZ3,CH2', '481,480,489,487,485,483', 'C,C,C,C,C,C',
        # 'SC,SC,SC,SC,SC,SC', 'TRP,TRP,TRP,TRP,TRP,TRP', '32,32,32,32,32,32', '',
        # 'SYSTEM,SYSTEM,SYSTEM,SYSTEM,SYSTEM,SYSTEM', 2.2205118417955965, 8.563199790120406,
        # 'parallel', '(G1)-ring···ring-(G2)']
        #
        # ['WATER-BRIDGE-1', 1051, 'O', 1052, 'O', 'BB', 'ILE', 69, '', 'SYSTEM', 1105, 'O', 1106, 'O', 'BB',
        # 'PRO', 73, '', 'SYSTEM', 63700, 'O', 63701, 'O', 'W', 'WAT', 9123, '', 'SYSTEM', 1.6196940030310187,
        # 1.8355580146241308, 166.068708622592, 172.6066731028374, 135.60529128877832,
        # 152.10048811999542, '(G1)-A···H–O–H···A-(G2)']
        # ['HYDROGEN-BOND', 4027, 'N', 4028, 'N', 'BB', 'HIE', 251, '', 'SYSTEM', 3974, 'O', 3975, 'O', 'BB',
        # 'TRP', 247, '', 'SYSTEM', 2.528356928860922, 140.6208709567424, 130.64733970961, '(G1)-D–H···A-(G2)']
        # etc
        #

        return tuple(json.loads(row) for row in dset)

    def _interaction_to_row(self, frame_number: int, rec: Tuple) -> List:
        label = rec[0]

        # PI–CATION
        if label == "PI-CATION":
            g1_role = rec[1]  # "ring" or "cation"
            g1 = rec[2:11]  # ALWAYS first
            g2 = rec[11:20]  # ALWAYS second

            dist = rec[20]
            theta = rec[21]
            kind = rec[22] if len(rec) > 22 else ""

            details = {
                "distance": dist,
                "theta": theta,
                "orientation": kind,
                "g1_role": g1_role,
            }

            # Swap element[3] ↔ type[4] to match INTERACTION_COLUMNS order
            return [
                frame_number,
                "PI-CATION",

                # atom / pseudo-atom 1 (type before element)
                g1[0], g1[1], g1[2], g1[4], g1[3], g1[5], g1[6], g1[7], g1[8],

                # atom / pseudo-atom 2 (type before element)
                g2[0], g2[1], g2[2], g2[4], g2[3], g2[5], g2[6], g2[7], g2[8],

                str(details),
            ]

        elif label == "PI-STACKING":
            g1 = rec[1:10]
            g2 = rec[10:19]

            dist = rec[19]
            angle = rec[20]

            stacking_type = rec[21] if len(rec) > 21 else ""
            orient = rec[22] if len(rec) > 22 else ""

            def ring9_to_row_fields(r9):
                indices, names, bonded, elems, atypes, resn, resid, chain, segid = r9
                return (
                    indices,
                    names,
                    bonded,
                    atypes,
                    elems,
                    resn,
                    resid,
                    chain or "",
                    segid or "",

                )

            a1 = ring9_to_row_fields(g1)
            a2 = ring9_to_row_fields(g2)

            details = {
                "distance": dist,
                "angle": angle,
                "stacking_type": stacking_type,
                "orientation": orient,
            }

            return [
                frame_number,
                "PI-STACKING",
                a1[0], a1[1], a1[2], a1[3], a1[4], a1[5], a1[6], a1[7], a1[8],
                a2[0], a2[1], a2[2], a2[3], a2[4], a2[5], a2[6], a2[7], a2[8],
                str(details),
            ]

        # Default: atom–atom & water-bridge interactions
        row = [
            frame_number,
            label,

            # atom 1
            rec[1],  # atom1_index
            rec[2],  # atom1_name
            rec[3],  # atom1_id
            rec[5],  # atom1_type
            rec[4],  # atom1_element
            rec[6],  # atom1_resname
            rec[7],  # atom1_resid
            rec[8],  # atom1_chainid
            rec[9],  # atom1_segid

            # atom 2
            rec[10],  # atom2_index
            rec[11],  # atom2_name
            rec[12],  # atom2_id
            rec[14],  # atom2_type
            rec[13],  # atom2_element
            rec[15],  # atom2_resname
            rec[16],  # atom2_resid
            rec[17],  # atom2_chainid
            rec[18],  # atom2_segid
        ]

        # everything else is opaque details
        details_raw = rec[19:] if len(rec) > 19 else ()

        if details_raw:
            schema = self.INTERACTION_DETAIL_SCHEMA.get(label)
            if schema:
                details = {k: v for k, v in zip(schema, details_raw)}
            else:
                # fallback for unknown interaction types
                details = {"values": details_raw}
        else:
            details = {}

        row.append(str(details))
        return row

    def _details_to_dict(self, label: str, tail: Tuple) -> Dict:
        """
        Converts interaction details described by a label and a tuple of data
        into a dictionary representation. The keys for the dictionary are
        retrieved based on the label from a pre-defined schema. If the label
        does not correspond to any schema keys, the raw data is returned as
        a dictionary.

        :param label: The label identifying the schema for interpreting the
            interaction details.
        :type label: str
        :param tail: The tuple containing the interaction data to be mapped
            to the schema keys.
        :type tail: Tuple
        :return: A dictionary where keys are the schema-defined identifiers
            and values are the corresponding elements of the input tuple. If
            the label is not found in the schema, the tuple is returned as
            raw data in a dictionary with the key 'raw'.
        :rtype: Dict
        """
        keys = self.INTERACTION_DETAIL_SCHEMA.get(label)

        if not keys:
            # fallback: preserve raw data
            return {"raw": tail}

        out = {}
        for k, v in zip(keys, tail):
            out[k] = v

        return out

    def _iter_frames(self, group_name: str):
        """
        Iterates over frames within a specified group in an HDF5 file. The method
        analyzes group contents, identifies keys starting with "frame_", and sorts
        them numerically before yielding each frame index.

        :param group_name: The group name within the HDF5 file to inspect.
        :type group_name: str
        :return: Generator yielding sorted frame indices within the specified group.
        :rtype: Iterator[int]
        """
        root = self.file[group_name]
        frames = []

        for name in root.keys():
            if name.startswith("frame_"):
                frames.append(int(name.split("_")[1]))

        for frame in sorted(frames):
            yield frame


    def build_interaction_modes(self, *, group_name: str, begin: int, end: int, step: int,
                                mode1: bool = True, mode2: bool = True, mode3: bool = True,
                                debug: bool = False, debug_frames: int = 1) -> None:
        """
        Build interaction modes from per-frame JSON-per-row interaction datasets.

        Definitions (residue-level keys):
          - Mode 1: Count EVERY row occurrence (atom-level multiplicity preserved)
          - Mode 2: Count ONCE per (residue-key) per frame
          - Mode 3: HYDROPHOBIC counted ONCE per (residue-key) per frame; all others counted like Mode 1

        Storage:
          /<group_name>/modes/<modeN>/table  (json-per-row)
          Each row: {"key": (<r1>, <r2>, <label>), "count": int, "frequency": float}

        Notes:
          - `extract_mode_key(tuple(rec))` must return (r1, r2, label) at residue level.
          - `rec[0]` must be the raw interaction label (e.g. "HYDROPHOBIC", "PI-CATION", ...).
        """

        DESCRIPTION_MODE1: Final[str] = (
            "Count frequency interactions analysis data (residue-level, count all occurrences)"
        )
        DESCRIPTION_MODE2: Final[str] = (
            "Residue-level once-per-frame interactions analysis data (deduplicate within frame)"
        )
        DESCRIPTION_MODE3: Final[str] = (
            "Hybrid interactions analysis data (hydrophobic once-per-frame; others count all occurrences)"
        )

        # Lazy import: keep MDAnalysis-heavy analyzer module out of fileio's import path
        from pharmacon.analyzer.interactions import extract_mode_key

        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        # Discover frames (ONLY under group_name)
        frame_indices = list(self._iter_frames(group_name))

        # Apply begin/end/step (end is inclusive)
        frames = [
            f for f in frame_indices
            if f >= begin and f <= end and ((f - begin) % step == 0)
        ]


        if not frames:
            raise ValueError("No frames selected for interaction mode construction")

        n_frames = float(len(frames))  # normalization base

        mode1_counts: dict = defaultdict(int)
        mode2_counts: dict = defaultdict(int)
        mode3_counts: dict = defaultdict(int)

        dbg_frames_left = max(0, int(debug_frames)) if debug else 0

        for frame in frames:
            dset_path = f"{group_name}/frame_{frame}/interactions"
            if dset_path not in self.file:
                # If frames list is derived from _iter_frames, this shouldn't happen.
                # But keep it safe.
                continue

            dset = self.file[dset_path]

            # Mode2: once per residue-key per frame
            seen_mode2: set = set()

            # Mode3: hydrophobic once per residue-key per frame
            seen_mode3_hydro: set = set()

            # Optional debug: verify multiplicities per frame at residue-level
            if dbg_frames_left > 0:
                logger.debug("DEBUG FRAME: %d", frame)
                logger.debug("Dataset path: %s", dset_path)
                logger.debug("Dataset size: %d", int(getattr(dset, "size", 0)))
                # Count residue-keys per label (to prove whether mode1!=mode2 should happen)
                per_label_counts: dict[str, Counter] = defaultdict(Counter)

                # Peek the first row structure
                for i_row, row in enumerate(dset):
                    rec = json.loads(row)
                    if i_row == 0:
                        logger.debug("----- RAW STORED ROW -----")
                        logger.debug("Length: %d", len(rec))
                        logger.debug("Label: %s", rec[0] if rec else None)
                        for i, v in enumerate(rec):
                            logger.debug("[%d] -> %s (%s)", i, v, type(v))
                        logger.debug("--------------------------")

                    key = extract_mode_key(tuple(rec))
                    if key is None:
                        continue
                    label = str(rec[0])
                    per_label_counts[label][key] += 1

                # Log only labels where there is multiplicity > 1 for any key
                for lbl, ctr in per_label_counts.items():
                    multi = [(k, c) for k, c in ctr.items() if c > 1]
                    if multi:
                        logger.debug("Multiplicity detected for label=%s:", lbl)
                        for k, c in sorted(multi, key=lambda x: -x[1])[:25]:
                            logger.debug("  %s -> %d", k, c)

                dbg_frames_left -= 1

            # Main counting
            for row in dset:
                rec = json.loads(row)
                if not rec:
                    continue

                key = extract_mode_key(tuple(rec))
                if key is None:
                    continue

                label = str(rec[0])

                # MODE 1: count every occurrence (aggregated by residue-key)
                if mode1:
                    mode1_counts[key] += 1

                # MODE 2: count once per residue-key per frame
                if mode2:
                    if key not in seen_mode2:
                        mode2_counts[key] += 1
                        seen_mode2.add(key)

                # MODE 3: hydrophobic once-per-frame; others count all occurrences
                if mode3:
                    if label == "HYDROPHOBIC":
                        if key not in seen_mode3_hydro:
                            mode3_counts[key] += 1
                            seen_mode3_hydro.add(key)
                    else:
                        mode3_counts[key] += 1

        # Write modes to HDF5
        modes_root = f"{group_name}/modes"

        if self.group_exists(modes_root):
            self.delete_group(modes_root)
        self.create_group(modes_root)

        def _write_mode(name: str, counts: dict, description: str) -> None:
            records = [
                {
                    "key": key,
                    "count": int(count),
                    "frequency": float(count) / n_frames,
                }
                for key, count in counts.items()
            ]

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype("utf-8"),
            )

            mode_group = f"{modes_root}/{name}"
            self.create_group(mode_group)
            self.create_dataset(
                group_name=mode_group,
                dataset_name="table",
                data=data,
                metadata={
                    "n_frames": str(int(n_frames)),
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

            self.add_group_metadata(
                group_name=mode_group,
                metadata={
                    "description": description,
                    "mode": name,
                    "completed": "True",
                },
                overwrite=True,
            )

        if mode1:
            _write_mode("mode1", mode1_counts, DESCRIPTION_MODE1)
        if mode2:
            _write_mode("mode2", mode2_counts, DESCRIPTION_MODE2)
        if mode3:
            _write_mode("mode3", mode3_counts, DESCRIPTION_MODE3)

        self.add_group_metadata(
            group_name=modes_root,
            metadata={
                "completed": "True",
                "n_frames": str(int(n_frames)),
                "frame_begin": str(begin),
                "frame_end": str(end),
                "frame_step": str(step),
                "end_inclusive": "True",
            },
            overwrite=True,
        )

    def write_interactions_to_csv(self, output_file: str | Path, *, group_name: str) -> None:
        """
        Writes interaction data to a CSV file.

        This method processes interaction data from a specified interaction group
        (e.g. pp_interactions, pl_interactions, hbonds) and writes it to a CSV file.

        :param output_file: The path to the output CSV file.
        :param group_name: Name of the interaction group to export.
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        if not self.group_exists(group_name):
            raise KeyError(f"Interaction group does not exist: {group_name}")

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(self.INTERACTION_COLUMNS)

            root = self.file[group_name]

            for frame in sorted(
                    int(k.split("_")[1])
                    for k in root.keys()
                    if k.startswith("frame_")
            ):
                interactions = self.read_frame_interactions(
                    frame_index=frame,
                    group_name=group_name,
                )

                for rec in interactions:
                    writer.writerow(self._interaction_to_row(frame, rec))

    def write_interaction_modes_to_csv(self, output_files: List[str | Path], *, group_name: str = "pp_interactions",
                                       mode1: bool = True, mode2: bool = True,  mode3: bool = True) -> None:
        """
        Writes interaction mode data from the file to specified CSV output files. Each mode is
        written to a separate CSV file specified by `output_files`. The function processes the
        `mode1`, `mode2`, and `mode3` flags to determine which interaction modes to consider.
        The output CSV files contain information about the interactions along with unique
        statistics such as count and frequency.

        :param output_files: A list of paths where the interaction mode data will be
            saved as CSV files. The number of paths must match the number of enabled modes.
        :param group_name: Name of the group within the dataset from which interaction mode
            data is fetched. Defaults to "pp_interactions".
        :param mode1: A boolean indicating whether to process and save data for "mode1".
            Defaults to True.
        :param mode2: A boolean indicating whether to process and save data for "mode2".
            Defaults to True.
        :param mode3: A boolean indicating whether to process and save data for "mode3".
            Defaults to True.
        :return: This function does not return a value but writes CSV files.
        :raises ValueError: If the number of output files does not match the number of enabled modes.
        :raises KeyError: If the specified mode dataset is not found in the file.
        """
        modes: List[str] = []
        if mode1:
            modes.append("mode1")
        if mode2:
            modes.append("mode2")
        if mode3:
            modes.append("mode3")

        if len(output_files) != len(modes):
            raise ValueError(
                f"Expected {len(modes)} output files, got {len(output_files)}"
            )

        modes_root = f"{group_name}/modes"

        for mode_name, out_path in zip(modes, output_files):
            out_path = Path(out_path).expanduser().resolve()

            dset_path = f"{modes_root}/{mode_name}/table"
            if dset_path not in self.file:
                raise KeyError(f"Mode dataset not found: {dset_path}")

            dset = self.file[dset_path]

            with open(out_path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)

                # NOTE: no "mode" column
                writer.writerow([
                    "interaction",
                    "res1_resname", "res1_resid", "res1_chainid", "res1_segid",
                    "res2_resname", "res2_resid", "res2_chainid", "res2_segid",
                    "count",
                    "frequency",
                ])

                for row in dset:
                    rec = json.loads(row)

                    (r1, r2, label) = rec["key"]
                    count = rec["count"]
                    freq = rec["frequency"]

                    writer.writerow([
                        label,
                        r1[0], r1[1], r1[2], r1[3],
                        r2[0], r2[1], r2[2], r2[3],
                        count,
                        f"{freq:.3f}",  # <<< 3 decimal places
                    ])

    def read_merged_interactions(self, *, group_name: str, mode_name: str) -> List[Dict]:
        """
        Read merged interaction mode records from a PTA file.

        Data are read from:
            /<group_name>/modes_merged/<mode_name>/table

        Each row is JSON-per-row and decoded into a Python dict.
        The interaction key is normalized into tuples.

        :param group_name: Interaction group name (e.g. "pl-interactions")
        :param mode_name: Mode name ("mode1", "mode2", "mode3")
        :return: List of merged interaction records
        :raises KeyError: If merged mode table does not exist
        :raises ValueError: If dataset is empty or malformed
        """

        dset_path = f"{group_name}/modes_merged/{mode_name}/table"

        if dset_path not in self.file:
            raise KeyError(f"Merged interaction mode not found: {dset_path}")

        dset = self.file[dset_path]

        if dset.size == 0:
            raise ValueError(f"Merged interaction table is empty: {dset_path}")

        records: list[dict] = []

        for i, row in enumerate(dset):
            rec = json.loads(row)

            # --- validate required fields ---
            for field in ("key", "mean_frequency", "std_frequency", "n_files"):
                if field not in rec:
                    raise ValueError(
                        f"Missing field '{field}' in merged record {i}"
                    )

            # --- normalize key ---
            key = rec["key"]
            if not isinstance(key, (list, tuple)) or len(key) != 3:
                raise ValueError(f"Invalid key in record {i}: {key}")

            r1, r2, label = key

            rec["key"] = (
                tuple(r1),
                tuple(r2),
                str(label),
            )

            rec["mean_frequency"] = float(rec["mean_frequency"])
            rec["std_frequency"] = float(rec["std_frequency"])
            rec["n_files"] = int(rec["n_files"])

            records.append(rec)

        return records

    def write_interactions_to_tsv(self, output_file: str | Path, *, group_name: str) -> None:
        """
        Writes interaction data to a TSV file.

        This method processes interaction data from a specified interaction group
        (e.g. pp_interactions, pl_interactions, hbonds) and writes it to a TSV file.

        :param output_file: The target file path for saving the interaction data.
        :param group_name: Name of the interaction group to export.
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        if not self.group_exists(group_name):
            raise KeyError(f"Interaction group does not exist: {group_name}")

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")

            # header
            writer.writerow(self.INTERACTION_COLUMNS)

            root = self.file[group_name]

            for frame in sorted(
                    int(k.split("_")[1])
                    for k in root.keys()
                    if k.startswith("frame_")
            ):
                interactions = self.read_frame_interactions(
                    frame_index=frame,
                    group_name=group_name,
                )

                for rec in interactions:
                    writer.writerow(self._interaction_to_row(frame, rec))

    def write_interaction_modes_to_tsv(self, output_files: List[str | Path], *, group_name: str = "pp_interactions",
                                       mode1: bool = True, mode2: bool = True, mode3: bool = True) -> None:
        """
        Writes specified interaction modes data into individual TSV files.

        This method processes interaction modes based on the provided arguments and writes each mode's data into
        respective TSV files. TSV files maintain tabular structure with a set of predefined headers and include
        information on pairs of interacting residues along with interaction frequency and count. The number of
        output files provided must match the number of selected modes.

        :param output_files: List of file paths where the data for each interaction mode will be written.
        :param group_name: Base group name in the dataset hierarchy where the modes data exists. Defaults to "pp_interactions".
        :param mode1: Boolean to include data for mode1. Defaults to True.
        :param mode2: Boolean to include data for mode2. Defaults to True.
        :param mode3: Boolean to include data for mode3. Defaults to True.
        :return: None
        :raises ValueError: If the count of output_files is not equal to the number of selected modes.
        :raises KeyError: If the required mode dataset is not found in the current file.
        """

        modes: List[str] = []
        if mode1:
            modes.append("mode1")
        if mode2:
            modes.append("mode2")
        if mode3:
            modes.append("mode3")

        if len(output_files) != len(modes):
            raise ValueError(
                f"Expected {len(modes)} output files, got {len(output_files)}"
            )

        modes_root = f"{group_name}/modes"

        for mode_name, out_path in zip(modes, output_files):
            out_path = Path(out_path).expanduser().resolve()

            dset_path = f"{modes_root}/{mode_name}/table"
            if dset_path not in self.file:
                raise KeyError(f"Mode dataset not found: {dset_path}")

            dset = self.file[dset_path]

            with open(out_path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh, delimiter="\t")

                # header (NO mode column)
                writer.writerow([
                    "interaction",
                    "res1_resname", "res1_resid", "res1_chainid", "res1_segid",
                    "res2_resname", "res2_resid", "res2_chainid", "res2_segid",
                    "count",
                    "frequency",
                ])

                for row in dset:
                    rec = json.loads(row)

                    (r1, r2, label) = rec["key"]
                    count = rec["count"]
                    freq = rec["frequency"]

                    writer.writerow([
                        label,
                        r1[0], r1[1], r1[2], r1[3],
                        r2[0], r2[1], r2[2], r2[3],
                        count,
                        f"{freq:.3f}",
                    ])

    def write_frame_distances(self, *, frame_index: int, distances: Iterable[Tuple[str, str, float]],
                              group_name: str, time_ps: float | None = None,  overwrite: bool = False) -> None:
        """
        Writes distance data for a specific frame into the HDF5 file.

        Distances are stored as JSON-per-row inside:

            /<group_name>/frame_<frame_index>/distances

        Each record has the schema:
            {
                "label": <str>,
                "method": <str>,
                "distance": <float>
            }

        Frame-level time information is stored as dataset metadata:
            - time_ps
            - time_ns
            - time_us

        This method is frame-atomic: each call fully defines one frame.

        :param frame_index: Frame index (must be >= 0).
        :param distances: Iterable of (label, method, distance).
        :param group_name: Parent group name.
        :param time_ps: Simulation time in picoseconds (optional).
        :param overwrite: Whether to overwrite existing frame data.
        :raises ValueError: If frame_index is invalid.
        :raises FileExistsError: If frame exists and overwrite=False.
        """

        # --- validate frame index ---
        if not isinstance(frame_index, int) or frame_index < 0:
            raise ValueError("frame_index must be a non-negative integer")

        frame_group = f"{group_name}/frame_{frame_index}"

        # --- overwrite semantics ---
        if self.group_exists(frame_group):
            if not overwrite:
                raise FileExistsError(
                    f"Frame {frame_index} already exists under '{group_name}'"
                )
            self.delete_group(frame_group)

        self.create_group(frame_group)

        # --- materialize records (small, per-frame) ---
        records = [
            {
                "label": str(label),
                "method": str(method),
                "distance": float(distance),
            }
            for label, method, distance in distances
        ]

        # --- serialize JSON-per-row ---
        if records:
            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
        else:
            data = np.empty(
                (0,),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

        # --- time metadata ---
        meta = {
            "frame_index": str(frame_index),
            "n_distances": str(len(data)),
            "format": "json-per-row",
            "schema": "{label, method, distance}",
        }

        if time_ps is not None:
            time_ps = float(time_ps)
            meta.update({
                "time_ps": f"{time_ps:.6f}",
                "time_ns": f"{time_ps / 1_000.0:.6f}",
                "time_us": f"{time_ps / 1_000_000.0:.6f}",
            })

        # --- write dataset + metadata ---
        self.create_dataset(
            group_name=frame_group,
            dataset_name="distances",
            data=data,
            metadata=meta,
        )

    def read_distance_data(self, *, group_name: str = "distances") -> np.ndarray:
        """
        Reads and processes distance data from a specified group within the file. The function collects
        information from datasets present in the specified group and organizes it into a NumPy array
        containing structured distance data. If the group does not exist or if no valid data is found,
        an exception is raised.

        :param group_name: The name of the group from which to read the distance data. The default value is "distances".
        :type group_name: str
        :raises ValueError: If the specified group does not exist.
        :raises ValueError: If no valid distance data is found within the specified group.
        :return: A NumPy array containing structured distance data. Each row of the array contains
            the frame index, time in picoseconds, label, method, and distance as extracted from the datasets.
        :rtype: np.ndarray
        """
        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        rows = []

        for frame in self._iter_frames(group_name):
            dset_path = f"{group_name}/frame_{frame}/distances"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]
            attrs = dset.attrs

            time_ps = float(attrs.get("time_ps", 0.0))

            for row in dset:
                rec = json.loads(row)
                rows.append([
                    frame,
                    time_ps,
                    rec["label"],
                    rec["method"],
                    float(rec["distance"]),
                ])

        if not rows:
            raise ValueError("No distance data found")

        return np.asarray(rows, dtype=object)

    def build_distance_statistics(self, *, group_name: str, begin: int, end: int, step: int,) -> None:
        """
        Build distance statistics across frames for each (label, method) pair.

        Reads per-frame distance records from:

            /<group_name>/frame_<frame>/distances

        and produces aggregated statistics stored at:

            /<group_name>/statistics/table

        Statistics computed per (label, method):
            - count
            - mean
            - std
            - min
            - max

        Data are written as JSON-per-row and processed in a memory-safe streaming manner.

        :param group_name: Parent group containing frame distance data.
        :param begin: First frame index (inclusive).
        :param end: Last frame index (inclusive).
        :param step: Frame stride.
        :raises ValueError: If no frames are selected or group does not exist.
        """

        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        # discover frames
        frame_indices = list(self._iter_frames(group_name))

        frames = [
            f for f in frame_indices
            if f >= begin and f <= end and ((f - begin) % step == 0)
        ]

        if not frames:
            raise ValueError("No frames selected for distance statistics")

        # accumulators
        sums = defaultdict(float)
        sums_sq = defaultdict(float)
        mins = {}
        maxs = {}
        counts = defaultdict(int)

        # stream frames
        for frame in frames:
            dset_path = f"{group_name}/frame_{frame}/distances"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]

            for row in dset:
                rec = json.loads(row)
                key = (rec["label"], rec["method"])
                val = float(rec["distance"])

                counts[key] += 1
                sums[key] += val
                sums_sq[key] += val * val

                if key not in mins:
                    mins[key] = val
                    maxs[key] = val
                else:
                    mins[key] = min(mins[key], val)
                    maxs[key] = max(maxs[key], val)

        # build output records
        records = []
        for (label, method), n in counts.items():
            mean = sums[(label, method)] / n
            var = (sums_sq[(label, method)] / n) - mean * mean
            std = var ** 0.5 if var > 0 else 0.0

            records.append(
                {
                    "label": label,
                    "method": method,
                    "count": n,
                    "mean": mean,
                    "std": std,
                    "min": mins[(label, method)],
                    "max": maxs[(label, method)],
                }
            )

        # write HDF5
        stats_group = f"{group_name}/statistics"

        if self.group_exists(stats_group):
            self.delete_group(stats_group)

        self.create_group(stats_group)

        data = np.array(
            [json.dumps(r) for r in records],
            dtype=h5py.string_dtype("utf-8"),
        )

        self.create_dataset(
            group_name=stats_group,
            dataset_name="table",
            data=data,
            metadata={
                "n_rows": str(len(records)),
                "frame_begin": str(begin),
                "frame_end": str(end),
                "frame_step": str(step),
                "end_inclusive": "True",
                "schema": "{label, method, count, mean, std, min, max}",
                "format": "json-per-row",
            },
        )

        self.add_group_metadata(
            group_name=stats_group,
            metadata={
                "completed": "True",
                "n_frames": str(len(frames)),
                "description": "Aggregated distance statistics per (label, method)",
            },
            overwrite=True,
        )

    def write_distances_to_csv(self, output_file: str | Path, *, group_name: str = "distances",
                            is_merged: bool = False) -> None:
        """
        Writes distances data to a CSV file. The method processes frames from the specified
        group and writes either detailed or merged summary information to the output file.
        The function creates and writes to the file in CSV format, ensuring appropriate
        headers are added depending on the `is_merged` flag.

        :param output_file: The target file path to create or overwrite with the distances
            data. Accepts a string or a Path object.
        :param group_name: The name of the group containing the data to process.
            Defaults to "distances".
        :param is_merged: Flag indicating whether the output data should be written in
            merged form (mean and standard deviation included) or detailed form.
            Defaults to False.
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)

            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "method",
                    "distance",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "method",
                    "mean",
                    "std",
                ])

            for frame in self._iter_frames(group_name):
                dset = self.file[f"{group_name}/frame_{frame}/distances"]
                attrs = dset.attrs

                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["method"],
                            rec["distance"],
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["method"],
                            rec["mean"],
                            rec["std"],
                        ])

    def write_distances_to_tsv(self, output_file: str | Path, *, group_name: str = "distances",
                               is_merged: bool = False) -> None:
        """
        Writes distances data to a TSV file. The method processes frames from the specified
        group and writes either detailed or merged summary information to the output file.
        The function creates and writes to the file in TSV format, ensuring appropriate
        headers are added depending on the `is_merged` flag.

        :param output_file: The target file path to create or overwrite with the distances
            data. Accepts a string or a Path object.
        :param group_name: The name of the group containing the data to process.
            Defaults to "distances".
        :param is_merged: Flag indicating whether the output data should be written in
            merged form (mean and standard deviation included) or detailed form.
            Defaults to False.
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")

            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "method",
                    "distance",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "method",
                    "mean",
                    "std",
                ])

            for frame in self._iter_frames(group_name):
                dset = self.file[f"{group_name}/frame_{frame}/distances"]
                attrs = dset.attrs

                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["method"],
                            rec["distance"],
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["method"],
                            rec["mean"],
                            rec["std"],
                        ])

    def write_distance_statistics_to_csv(self, output_file: str | Path, *, group_name: str = "distances") -> None:
        """
        Writes distance statistics to a CSV file. The statistics are fetched from the
        specified group in the input file and include details such as label, method, count,
        mean, standard deviation, minimum, and maximum values for each group.

        :param output_file: The path to the output CSV file where the statistics will be written.
        :param group_name: The group name in the file from which statistics are retrieved.
                           Defaults to "distances".
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "label",
                "method",
                "count",
                "mean",
                "std",
                "min",
                "max",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["method"],
                    rec["count"],
                    f"{rec['mean']:.6f}",
                    f"{rec['std']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                ])

    def write_distance_statistics_to_tsv(self, output_file: str | Path, *, group_name: str = "distances") -> None:
        """
        Writes distance statistics from an HDF5 dataset to a tab-separated value (TSV) file.

        This function extracts statistical data related to distances from a specific group
        in an HDF5 file and writes it to a TSV file. Each row in the resulting file
        represents a set of statistics (e.g., count, mean, standard deviation, minimum, and
        maximum values) associated with a particular label and method.

        :param output_file: Path to the output TSV file where the statistics will be written.
                           The path is expanded to include user directories and resolved to
                           an absolute path.
        :type output_file: str | Path
        :param group_name: Name of the group from which to retrieve the statistics table in
                           the HDF5 file. Defaults to "distances".
        :type group_name: str
        :return: This function does not return a value as it writes to a file.
        :rtype: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "label",
                "method",
                "count",
                "mean",
                "std",
                "min",
                "max",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["method"],
                    rec["count"],
                    f"{rec['mean']:.6f}",
                    f"{rec['std']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                ])

    def write_frame_angles(self, *,  frame_index: int, angles: Iterable[Tuple[str, str, float]],
                           group_name: str, time_ps: float | None = None, overwrite: bool = False) -> None:
        """
        This method writes frame angles to a specific group within a data container. It creates a new group
        for the frame or overwrites the existing one based on the provided parameters. The method serializes
        the angles data into JSON records and supports additional metadata like time in different resolutions.

        /<group_name>/frame_<frame_index>/angles

        :param frame_index: Index of the frame to write.
            Must be a non-negative integer.
        :param angles: Collection of angle records, where each record is a tuple containing label, kind, and value.
        :param group_name: Name of the group under which the frame angles are stored.
        :param time_ps: The timestamp in picoseconds (optional).
        :param overwrite: Boolean flag indicating whether to overwrite existing data for the same frame.
        :return: None
        """

        # --- validate frame index ---
        if not isinstance(frame_index, int) or frame_index < 0:
            raise ValueError("frame_index must be a non-negative integer")

        frame_group = f"{group_name}/frame_{frame_index}"

        # --- overwrite semantics ---
        if self.group_exists(frame_group):
            if not overwrite:
                raise FileExistsError(
                    f"Frame {frame_index} already exists under '{group_name}'"
                )
            self.delete_group(frame_group)

        self.create_group(frame_group)

        # --- materialize records (small, per-frame) ---
        records = [
            {
                "label": str(label),
                "kind": str(kind),
                "value": float(value),
            }
            for label, kind, value in angles
        ]

        # --- serialize JSON-per-row ---
        if records:
            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
        else:
            data = np.empty(
                (0,),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

        # --- metadata ---
        meta = {
            "frame_index": str(frame_index),
            "n_angles": str(len(data)),
            "format": "json-per-row",
            "schema": "{label, kind, value}",
            "units": "degrees",
        }

        if time_ps is not None:
            time_ps = float(time_ps)
            meta.update({
                "time_ps": f"{time_ps:.6f}",
                "time_ns": f"{time_ps / 1_000.0:.6f}",
                "time_us": f"{time_ps / 1_000_000.0:.6f}",
            })

        # --- write dataset ---
        self.create_dataset(
            group_name=frame_group,
            dataset_name="angles",
            data=data,
            metadata=meta,
        )

    def read_angle_data(self, *, group_name: str = "angles") -> np.ndarray:
        """
        Reads and retrieves angle data from the specified group in the dataset.

        This method iterates through the frames in the given group, extracts angle
        information, and compiles it into a numpy array. Each record includes the
        frame number, the timestamp in picoseconds, the label of the angle, its
        type (kind), and its value.

        :param group_name: The name of the group from which to read the angle data.
            Default is "angles".
        :type group_name: str
        :return: A numpy array containing the angle data. Each entry in the array
            represents a record with the frame number, timestamp, label, kind, and
            value of the angle in the specified group.
        :rtype: numpy.ndarray
        :raises ValueError: If the specified group does not exist or if no angle
            data is found in the specified group.
        """

        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        rows = []

        for frame in self._iter_frames(group_name):
            dset_path = f"{group_name}/frame_{frame}/angles"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]
            attrs = dset.attrs

            time_ps = float(attrs.get("time_ps", 0.0))

            for row in dset:
                rec = json.loads(row)
                rows.append([
                    frame,
                    time_ps,
                    rec["label"],
                    rec["kind"],
                    float(rec["value"]),
                ])

        if not rows:
            raise ValueError("No angle data found")

        return np.asarray(rows, dtype=object)

    def build_angle_statistics(self, *, group_name: str, begin: int, end: int, step: int) -> None:
        """
        Build angle statistics across frames for each (label, kind) pair.
        Reads per-frame angle records from:

            /<group_name>/frame_<frame>/angles
        and produces aggregated statistics stored at:
            /<group_name>/statistics/table

        Statistics computed per (label, kind):
            - count
            - mean
            - std
            - min
            - max

        Data are written as JSON-per-row and processed in a memory-safe streaming manner.

        :param group_name: Parent group containing frame angle data (e.g. "angles").
        :param begin: First frame index (inclusive).
        :param end: Last frame index (inclusive).
        :param step: Frame stride.
        :raises ValueError: If no frames are selected or group does not exist.
        """

        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        # --- discover frames ---
        frame_indices = list(self._iter_frames(group_name))

        frames = [
            f for f in frame_indices
            if f >= begin and f <= end and ((f - begin) % step == 0)
        ]

        if not frames:
            raise ValueError("No frames selected for angle statistics")

        # --- accumulators ---
        sums = defaultdict(float)
        sums_sq = defaultdict(float)
        mins: Dict[tuple, float] = {}
        maxs: Dict[tuple, float] = {}
        counts = defaultdict(int)

        # --- stream frames ---
        for frame in frames:
            dset_path = f"{group_name}/frame_{frame}/angles"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]

            for row in dset:
                rec = json.loads(row)

                key = (rec["label"], rec["kind"])
                val = float(rec["value"])

                counts[key] += 1
                sums[key] += val
                sums_sq[key] += val * val

                if key not in mins:
                    mins[key] = val
                    maxs[key] = val
                else:
                    mins[key] = min(mins[key], val)
                    maxs[key] = max(maxs[key], val)

        # --- build output records ---
        records = []
        for (label, kind), n in counts.items():
            mean = sums[(label, kind)] / n
            var = (sums_sq[(label, kind)] / n) - mean * mean
            std = var ** 0.5 if var > 0 else 0.0

            records.append(
                {
                    "label": label,
                    "kind": kind,
                    "count": n,
                    "mean": mean,
                    "std": std,
                    "min": mins[(label, kind)],
                    "max": maxs[(label, kind)],
                }
            )

        # --- write HDF5 ---
        stats_group = f"{group_name}/statistics"

        if self.group_exists(stats_group):
            self.delete_group(stats_group)

        self.create_group(stats_group)

        data = np.array(
            [json.dumps(r) for r in records],
            dtype=h5py.string_dtype("utf-8"),
        )

        self.create_dataset(
            group_name=stats_group,
            dataset_name="table",
            data=data,
            metadata={
                "n_rows": str(len(records)),
                "frame_begin": str(begin),
                "frame_end": str(end),
                "frame_step": str(step),
                "end_inclusive": "True",
                "schema": "{label, kind, count, mean, std, min, max}",
                "units": "degrees",
                "format": "json-per-row",
            },
        )

        self.add_group_metadata(
            group_name=stats_group,
            metadata={
                "completed": "True",
                "n_frames": str(len(frames)),
                "description": "Aggregated angle statistics per (label, kind)",
            },
            overwrite=True,
        )

    def write_angles_to_csv(self, output_file: str | Path, *, group_name: str = "angles",
                            is_merged: bool = False) -> None:
        """
        Writes angle data to a CSV file. The method processes angle datasets for
        each frame under the specified group and writes the information to the
        given `output_file`. The output format depends on the `is_merged`
        parameter. Non-merged output includes individual values for angles,
        whereas merged output aggregates data with mean and standard deviation.

        :param output_file: The path to the CSV file where the data will be written.
                            This path is expanded to the user's home directory and
                            resolved to its absolute path.
        :type output_file: str | Path
        :param group_name: The name of the group in the HDF5 file to be processed.
                           Default is "angles".
        :type group_name: str
        :param is_merged: Specifies whether the angles data is aggregated (merged)
                          or written directly as individual records. If True, mean
                          and standard deviation are written instead of individual
                          values. Default is False.
        :type is_merged: bool
        :returns: None
        """
        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)

            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "kind",
                    "value",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "kind",
                    "mean",
                    "std",
                ])

            for frame in self._iter_frames(group_name):
                dset = self.file[f"{group_name}/frame_{frame}/angles"]
                attrs = dset.attrs

                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["kind"],
                            rec["value"],
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["kind"],
                            rec["mean"],
                            rec["std"],
                        ])

    def write_angles_to_tsv(self, output_file: str | Path, *, group_name: str = "angles", is_merged: bool = False) -> None:
        """
        Writes angle information to a TSV file. The method exports data in either a detailed or
        summarized format based on the `is_merged` parameter. Angles data are read from frames
        within a specified group and written to the specified file in tabular form.

        :param output_file: The path to the output file where the TSV data will be written.
        :type output_file: str | Path
        :param group_name: The name of the group in the file from which angle data will
                           be extracted. Defaults to "angles".
        :type group_name: str
        :param is_merged: A flag indicating whether the data should be exported in merged format
                          (e.g., with statistical mean and standard deviation for values). Defaults
                          to False.
        :type is_merged: bool
        :return: None
        """

        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")

            # --- header ---
            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "kind",
                    "value",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "kind",
                    "mean",
                    "std",
                ])

            # --- data ---
            for frame in self._iter_frames(group_name):
                dset = self.file[f"{group_name}/frame_{frame}/angles"]
                attrs = dset.attrs

                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["kind"],
                            rec["value"],
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            rec["kind"],
                            rec["mean"],
                            rec["std"],
                        ])

    def write_angle_statistics_to_csv(self, output_file: str | Path, *, group_name: str = "angles") -> None:
        """
        Writes angle statistics from a specified group in the HDF5 file to a CSV file.

        This function extracts statistical data from the given group in the HDF5 file
        and writes it to a specified CSV file. The resulting file will contain rows of
        data with headers describing the statistics for each angle, such as 'label',
        'count', 'mean', and others. This is particularly useful for exporting preprocessing
        or analysis statistics to a portable format.

        label, kind, count, mean, std, min, max

        :param output_file: Path to the output CSV file where angle statistics data will
            be written. Can be provided as a string or a Path instance.
        :param group_name: Name of the HDF5 group containing the statistics data.
            Defaults to "angles".
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "label",
                "kind",
                "count",
                "mean",
                "std",
                "min",
                "max",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["kind"],
                    rec["count"],
                    f"{rec['mean']:.6f}",
                    f"{rec['std']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                ])

    def write_angle_statistics_to_tsv(self, output_file: str | Path, *, group_name: str = "angles") -> None:
        """
        Writes angle statistics to a TSV file. This involves reading data from a dataset,
        processing the information, and writing it as tab-separated values into an output file.
        The output file includes details such as label, kind, count, mean, standard deviation,
        minimum, and maximum values.

        :param output_file: Path to the output TSV file where the angle statistics will be saved.
        :type output_file: str | Path
        :param group_name: Name of the group in the file to access the statistics data.
            Defaults to "angles".
        :type group_name: str
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "label",
                "kind",
                "count",
                "mean",
                "std",
                "min",
                "max",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["kind"],
                    rec["count"],
                    f"{rec['mean']:.6f}",
                    f"{rec['std']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                ])

    def write_pca_data(self, *, frame_index: int, pc_values: Iterable[float], variances: Iterable[float],
                       variance_ratios: Iterable[float], group_name: str = "pca", time_ps: float | None = None,
                       overwrite: bool = False) -> None:
        """
        Writes PCA projection data for a single frame.

        Data are stored as JSON-per-row at:

            /<group_name>/frame_<frame_index>/pca

        Each row has schema:
            {
                "pc": <int>,
                "value": <float>,
                "variance": <float>,
                "variance_ratio": <float>
            }

        Frame-level time metadata (ps/ns/us) are stored as dataset attributes.

        This method is frame-atomic.

        :param frame_index: Frame index (must be >= 0).
        :param pc_values: Iterable of PCA projection values for this frame.
        :param variances: Iterable of variances for each principal component.
        :param variance_ratios: Iterable of variance ratios for each component.
        :param group_name: Parent group name (default: "pca").
        :param time_ps: Simulation time in picoseconds (optional).
        :param overwrite: Whether to overwrite existing frame data.
        """

        # --- validate frame index ---
        if not isinstance(frame_index, int) or frame_index < 0:
            raise ValueError("frame_index must be a non-negative integer")

        pc_values = list(pc_values)
        variances = list(variances)
        variance_ratios = list(variance_ratios)

        if not (len(pc_values) == len(variances) == len(variance_ratios)):
            raise ValueError("pc_values, variances, and variance_ratios must have equal length")

        frame_group = f"{group_name}/frame_{frame_index}"

        # --- overwrite semantics ---
        if self.group_exists(frame_group):
            if not overwrite:
                raise FileExistsError(
                    f"Frame {frame_index} already exists under '{group_name}'"
                )
            self.delete_group(frame_group)

        self.create_group(frame_group)

        # --- materialize records ---
        records = [
            {
                "pc": i + 1,
                "value": float(val),
                "variance": float(var),
                "variance_ratio": float(ratio),
            }
            for i, (val, var, ratio) in enumerate(
                zip(pc_values, variances, variance_ratios)
            )
        ]

        # --- serialize JSON-per-row ---
        data = np.array(
            [json.dumps(r) for r in records],
            dtype=h5py.string_dtype(encoding="utf-8"),
        )

        # --- metadata ---
        meta = {
            "frame_index": str(frame_index),
            "n_components": str(len(records)),
            "format": "json-per-row",
            "schema": "{pc, value, variance, variance_ratio}",
        }

        if time_ps is not None:
            time_ps = float(time_ps)
            meta.update({
                "time_ps": f"{time_ps:.6f}",
                "time_ns": f"{time_ps / 1_000.0:.6f}",
                "time_us": f"{time_ps / 1_000_000.0:.6f}",
            })

        # --- write dataset ---
        self.create_dataset(
            group_name=frame_group,
            dataset_name="pca",
            data=data,
            metadata=meta,
        )

    def write_pca_to_csv(self, output_file: str | Path, *, group_name: str = "pca") -> None:
        """
        Writes principal component analysis (PCA) data to a specified CSV file. This function
        retrieves PCA data from a hierarchical data structure and writes it into a CSV file with
        relevant metadata and PCA-related metrics for further analysis or sharing.

        :param output_file: The target file where PCA data will be written. This should be a
            string or a Path object representing the file path.
        :param group_name: The name of the group in the hierarchical data structure from
            which PCA data is extracted. Defaults to "pca".
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "frame_index",
                "time_ps",
                "time_ns",
                "time_us",
                "pc",
                "value",
                "variance",
                "variance_ratio",
            ])

            for frame in self._iter_frames(group_name):
                dset_path = f"{group_name}/frame_{frame}/pca"
                dset = self.file[dset_path]

                attrs = dset.attrs
                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)
                    writer.writerow([
                        frame,
                        time_ps,
                        time_ns,
                        time_us,
                        rec["pc"],
                        f"{rec['value']:.6f}",
                        f"{rec['variance']:.6f}",
                        f"{rec['variance_ratio']:.6f}",
                    ])

    def write_pca_to_tsv(self, output_file: str | Path, *, group_name: str = "pca") -> None:
        """
        Writes PCA (Principal Component Analysis) data to a TSV (Tab-Separated Values) file.
        The data is retrieved from the specified group in the current HDF5 file, processed,
        and written frame by frame to the output TSV file. Each row in the output corresponds
        to a PCA record with information about the principal component, its value, variance,
        and variance ratio.

        :param output_file: Path to the output file where the TSV data will be written. The path
            will be resolved to an absolute path before use.
        :type output_file: str | Path
        :param group_name: The name of the group in the HDF5 file from which PCA data should be
            retrieved. Default is "pca".
        :type group_name: str
        :return: None
        """

        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "frame_index",
                "time_ps",
                "time_ns",
                "time_us",
                "pc",
                "value",
                "variance",
                "variance_ratio",
            ])

            for frame in self._iter_frames(group_name):
                dset_path = f"{group_name}/frame_{frame}/pca"
                dset = self.file[dset_path]

                attrs = dset.attrs
                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)
                    writer.writerow([
                        frame,
                        time_ps,
                        time_ns,
                        time_us,
                        rec["pc"],
                        f"{rec['value']:.6f}",
                        f"{rec['variance']:.6f}",
                        f"{rec['variance_ratio']:.6f}",
                    ])

    def write_rmsd_data(self, *, rmsd_array: np.ndarray, labels: List[str],
                        group_name: str = "rmsd", overwrite: bool = False) -> None:
        """
        Writes RMSD data to the specified group in the storage system. This method processes
        a 2D RMSD array, verifies its structure, and writes it into hierarchical groups.
        Each frame's data is stored in a separate subgroup under the provided group name.
        Optional reference RMSD columns and overwriting behavior are handled by this method.

        Stored at:
            /<group_name>/frame_<frame>/rmsd

        :param rmsd_array: A 2D numpy array where rows represent frames and columns include
                           frame index, time in picoseconds, and RMSD values for multiple groups.
        :param labels: A list of strings representing the labels for individual groups' RMSD values.
        :param group_name: A string representing the base name of the group to store the data.
                           Default is "rmsd".
        :param overwrite: A boolean flag indicating whether existing data for a frame should be
                          overwritten if it exists under the same group name. Default is False.
        :return: None
        """

        rmsd_array = np.asarray(rmsd_array)

        if rmsd_array.ndim != 2 or rmsd_array.size == 0:
            raise ValueError("Invalid or empty RMSD array")

        ncols = rmsd_array.shape[1]
        expected = 2 + len(labels)  # frame + time + groups

        # --- handle optional reference RMSD column ---
        if ncols == expected + 1:
            # legacy behavior: drop column index 2
            rmsd_array = np.delete(rmsd_array, 2, axis=1)
            ncols -= 1

        if ncols != expected:
            raise ValueError(
                f"Unexpected RMSD array shape {rmsd_array.shape}, "
                f"expected {expected} columns"
            )

        # columns now are:
        # frame_index, time_ps, rmsd_group1, rmsd_group2, ...
        n_groups = len(labels)

        for row in rmsd_array:
            frame = int(row[0])
            time_ps = float(row[1])

            frame_group = f"{group_name}/frame_{frame}"

            # overwrite semantics
            if self.group_exists(frame_group):
                if not overwrite:
                    raise FileExistsError(
                        f"Frame {frame} already exists under '{group_name}'"
                    )
                self.delete_group(frame_group)

            self.create_group(frame_group)

            # --- per-group RMSD records ---
            records = [
                {
                    "label": labels[i],
                    "value": float(row[i + 2]),
                }
                for i in range(n_groups)
            ]

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

            # --- metadata (legacy-compatible) ---
            meta = {
                "frame_index": str(frame),
                "n_groups": str(n_groups),
                "format": "json-per-row",
                "schema": "{label, value}",
                "time_ps": f"{time_ps:.6f}",
                "time_ns": f"{time_ps / 1_000.0:.6f}",
                "time_us": f"{time_ps / 1_000_000.0:.6f}",
            }

            self.create_dataset(
                group_name=frame_group,
                dataset_name="rmsd",
                data=data,
                metadata=meta,
            )

    def read_rmsd_data(self, *, group_name: str = "rmsd") -> Tuple[np.ndarray, List[str]]:
        """
        Reads RMSD (Root Mean Square Deviation) data from a specified group within a dataset.

        The function retrieves all RMSD frames from the provided group, verifies that the
        frames exist and are consistent, and reconstructs a dataset containing RMSD values
        along with relevant metadata. The function ensures that all frame labels match across
        frames and raises errors when data is missing or inconsistent.

        rmsd_array: np.ndarray of shape (n_frames, 2 + n_labels)
                Columns:
                    frame_index, time_ps, rmsd(label1), rmsd(label2), ...
            labels: List[str] in column order

        :param group_name: The name of the group containing RMSD data. Defaults to "rmsd".
        :return: A tuple containing:
                 - An ndarray with RMSD data including frame indices, time, and RMSD values.
                 - A list of strings representing the labels for the RMSD values.
        :rtype: Tuple[np.ndarray, List[str]]
        :raises ValueError: If the group does not exist, if no frames are found, if labels
                            across different frames are inconsistent, or if data reconstruction
                            fails.
        """
        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        frames = list(self._iter_frames(group_name))
        if not frames:
            raise ValueError("No RMSD frames found")

        frames = sorted(frames)

        rows = []
        labels: List[str] | None = None

        for frame in frames:
            dset_path = f"{group_name}/frame_{frame}/rmsd"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]
            attrs = dset.attrs

            time_ps = float(attrs.get("time_ps", 0.0))

            records = [json.loads(row) for row in dset]

            if not records:
                continue

            # establish label order ONCE (from first frame)
            if labels is None:
                labels = [rec["label"] for rec in records]
            else:
                frame_labels = [rec["label"] for rec in records]
                if frame_labels != labels:
                    raise ValueError(
                        f"Inconsistent RMSD labels in frame {frame}: "
                        f"{frame_labels} != {labels}"
                    )

            values = [float(rec["value"]) for rec in records]

            rows.append([frame, time_ps, *values])

        if not rows or labels is None:
            raise ValueError("No RMSD data could be reconstructed")

        rmsd_array = np.asarray(rows, dtype=float)

        return rmsd_array, labels

    def write_rmsd_statistics(self, *, group_name: str = "rmsd", begin: int | None = None, end: int | None = None,
                              step: int = 1) -> None:
        """
        Builds RMSD summary statistics across frames for each label.

        Statistics per label:
            - mean
            - min
            - max
            - median
            - n

        Stored at:
            /<group_name>/statistics/table

        Data format:
            JSON-per-row

        :param group_name: RMSD parent group (default: "rmsd")
        :param begin: First frame index (inclusive). None = auto-detect.
        :param end: Last frame index (inclusive). None = auto-detect.
        :param step: Frame stride.
        """

        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        # --- discover frames ---
        frames = list(self._iter_frames(group_name))
        if not frames:
            raise ValueError("No RMSD frames found")

        if begin is None:
            begin = frames[0]
        if end is None:
            end = frames[-1]

        selected_frames = [
            f for f in frames
            if f >= begin and f <= end and ((f - begin) % step == 0)
        ]

        if not selected_frames:
            raise ValueError("No frames selected for RMSD statistics")

        # --- accumulate values per label ---
        values: Dict[str, list] = defaultdict(list)

        for frame in selected_frames:
            dset_path = f"{group_name}/frame_{frame}/rmsd"
            if dset_path not in self.file:
                continue

            dset = self.file[dset_path]

            for row in dset:
                rec = json.loads(row)
                label = rec["label"]
                val = float(rec["value"])
                values[label].append(val)

        if not values:
            raise ValueError("No RMSD values collected")

        # --- build statistics ---
        records = []
        for label, vals in values.items():
            arr = np.asarray(vals, dtype=float)

            records.append(
                {
                    "label": label,
                    "mean": float(arr.mean()),
                    "min": float(arr.min()),
                    "max": float(arr.max()),
                    "median": float(np.median(arr)),
                    "n": int(arr.size),
                }
            )

        # --- write HDF5 ---
        stats_group = f"{group_name}/statistics"

        if self.group_exists(stats_group):
            self.delete_group(stats_group)

        self.create_group(stats_group)

        data = np.array(
            [json.dumps(r) for r in records],
            dtype=h5py.string_dtype("utf-8"),
        )

        self.create_dataset(
            group_name=stats_group,
            dataset_name="table",
            data=data,
            metadata={
                "n_rows": str(len(records)),
                "frame_begin": str(begin),
                "frame_end": str(end),
                "frame_step": str(step),
                "end_inclusive": "True",
                "schema": "{label, mean, min, max, median, n}",
                "format": "json-per-row",
            },
        )

        self.add_group_metadata(
            group_name=stats_group,
            metadata={
                "completed": "True",
                "n_frames": str(len(selected_frames)),
                "description": "Aggregated RMSD statistics per label",
            },
            overwrite=True,
        )

    def write_rmsd_data_to_csv(self, output_file: str | Path, *, group_name: str = "rmsd", is_merged: bool = False) -> None:
        """
        Writes Root Mean Square Deviation (RMSD) data to a specified CSV file. The function
        retrieves frame-wise RMSD data for a given group, formats it, and writes it into
        the CSV file. It supports writing either individual RMSD data or merged RMSD statistics
        (mean and standard deviation) depending on the given option.

        Non-merged schema:
            frame_index, time_ps, time_ns, time_us, label, value

        Merged schema:
            frame_index, time_ps, time_ns, time_us, label, mean, std

        :param output_file: Path to the output CSV file where the RMSD data will be stored.
        :type output_file: str | Path
        :param group_name: Name of the group in the file containing the RMSD data to be processed.
        :type group_name: str
        :param is_merged: If True, writes merged RMSD statistics (mean and standard
            deviation) to the CSV file. If False, writes individual RMSD data per frame.
        :type is_merged: bool
        :return: None
        """

        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)

            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "value",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "mean",
                    "std",
                ])

            for frame in self._iter_frames(group_name):
                dset_path = f"{group_name}/frame_{frame}/rmsd"
                dset = self.file[dset_path]

                attrs = dset.attrs
                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            f"{rec['value']:.6f}",
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            f"{rec['mean']:.6f}",
                            f"{rec['std']:.6f}",
                        ])

    def write_rmsd_data_to_tsv(self, output_file: str | Path, *, group_name: str = "rmsd", is_merged: bool = False) -> None:
        """
        Writes RMSD (Root Mean Square Deviation) data into a TSV (Tab-Separated Values) file using
        data from an HDF5 file. The output format may vary depending on whether the data is merged
        or not. Each row in the generated TSV file corresponds to a frame and its associated data.

        Non-merged schema:
            frame_index, time_ps, time_ns, time_us, label, value

        Merged schema:
            frame_index, time_ps, time_ns, time_us, label, mean, std

        This method retrieves RMSD data for a given group based on a group name and processes it
        for writing into the specified output file. If the RMSD data is merged, additional fields
        such as the standard deviation are included in the output.

        :param output_file: The target file path where TSV data will be written. This can
            be a string or a Path object.
        :param group_name: The name of the group in the HDF5 file to extract RMSD data
            from. Defaults to "rmsd".
        :param is_merged: A boolean indicating whether the RMSD data is merged. If True,
            additional columns for "mean" and "std" are included in the output. Defaults
            to False.
        :return: None
        """

        output_file = Path(output_file).expanduser().resolve()

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")

            if not is_merged:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "value",
                ])
            else:
                writer.writerow([
                    "frame_index",
                    "time_ps",
                    "time_ns",
                    "time_us",
                    "label",
                    "mean",
                    "std",
                ])

            for frame in self._iter_frames(group_name):
                dset_path = f"{group_name}/frame_{frame}/rmsd"
                dset = self.file[dset_path]

                attrs = dset.attrs
                time_ps = attrs.get("time_ps", "")
                if time_ps != "" and "time_ns" not in attrs:
                    time_ns = float(time_ps) / 1_000.0
                    time_us = float(time_ps) / 1_000_000.0
                else:
                    time_ns = attrs.get("time_ns", "")
                    time_us = attrs.get("time_us", "")

                for row in dset:
                    rec = json.loads(row)

                    if not is_merged:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            f"{rec['value']:.6f}",
                        ])
                    else:
                        writer.writerow([
                            frame,
                            time_ps,
                            time_ns,
                            time_us,
                            rec["label"],
                            f"{rec['mean']:.6f}",
                            f"{rec['std']:.6f}",
                        ])

    def write_rmsd_statistics_to_csv(self, output_file: str | Path, *, group_name: str = "rmsd") -> None:
        """
        Writes RMSD statistics stored in a dataset to a specified CSV file.

        This method reads data from a specified group in the file object and writes
        the statistics into a CSV file, including detailed columns such as label,
        mean, min, max, median, and n. Each row in the dataset is processed as JSON,
        and the extracted values are formatted and written into the CSV.

        :param output_file: The output file's path where the CSV data will be saved.
        :type output_file: str | Path
        :param group_name: The name of the group in the file containing the statistics
            table dataset. Defaults to "rmsd".
        :type group_name: str
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "label",
                "mean",
                "min",
                "max",
                "median",
                "n",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    f"{rec['mean']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                    f"{rec['median']:.6f}",
                    rec["n"],
                ])

    def write_rmsd_statistics_to_tsv(self, output_file: str | Path, *, group_name: str = "rmsd") -> None:
        """
        Writes RMSD (Root Mean Square Deviation) statistics to a TSV (Tab-Separated Values) file based on data
        read from a predefined group in an HDF5 file. The function retrieves data from the specified group under
        the `statistics/table` key, processes it, and writes it to the output file.

        :param output_file: Path to the output TSV file where the RMSD statistics will be written.
        :type output_file: str | Path
        :param group_name: The name of the group in the HDF5 file that contains the `statistics/table` dataset.
                           Defaults to "rmsd".
        :type group_name: str
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "label",
                "mean",
                "min",
                "max",
                "median",
                "n",
            ])

            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    f"{rec['mean']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                    f"{rec['median']:.6f}",
                    rec["n"],
                ])

    def write_rmsf_data(self, *, selections: List[Dict], group_name: str = "rmsf",
                        frame_begin: int, frame_end: int, frame_step: int = 1,
                        fitting_group: str | None = None,
                        overwrite: bool = False) -> None:
        """
        Writes per-atom RMSF data, one HDF5 subgroup per selection.

        Stored at:
            /<group_name>/selection_<label>/atoms

        Each entry in `selections` is a dict with keys:
            - "label":            str             — short name (e.g. "calpha")
            - "selection_string": str             — MDAnalysis selection
            - "atom_indices":     np.ndarray[int] — absolute topology indices
            - "resids":           np.ndarray[int]
            - "resnames":         np.ndarray[str]
            - "atom_names":       np.ndarray[str]
            - "rmsf":             np.ndarray[float] — per-atom RMSF

        All arrays within one selection must share the same length (n_atoms).

        :param selections: List of per-selection dicts.
        :param group_name: Parent HDF5 group (default "rmsf").
        :param frame_begin: First frame (inclusive) used to compute RMSF.
        :param frame_end: Last frame (inclusive) used to compute RMSF.
        :param frame_step: Frame stride used to compute RMSF.
        :param fitting_group: MDAnalysis selection used for the AlignTraj pre-pass.
        :param overwrite: If True, replace any existing per-selection subgroup.
        :raises ValueError: If `selections` is empty, labels are duplicated, or
            array lengths within a selection are inconsistent.
        :raises FileExistsError: If a selection subgroup already exists and
            `overwrite=False`.
        """
        if not selections:
            raise ValueError("No selections provided")

        if not self.group_exists(group_name):
            self.create_group(group_name)

        seen_labels: set = set()
        required_keys = ("label", "selection_string", "atom_indices",
                         "resids", "resnames", "atom_names", "rmsf")

        for sel in selections:
            missing = [k for k in required_keys if k not in sel]
            if missing:
                raise ValueError(
                    f"Selection entry missing required keys: {missing}"
                )

            label = str(sel["label"])
            if label in seen_labels:
                raise ValueError(f"Duplicate selection label: '{label}'")
            seen_labels.add(label)

            atom_indices = np.asarray(sel["atom_indices"])
            resids = np.asarray(sel["resids"])
            resnames = np.asarray(sel["resnames"])
            atom_names = np.asarray(sel["atom_names"])
            rmsf_values = np.asarray(sel["rmsf"], dtype=float)

            n_atoms = atom_indices.size
            if n_atoms == 0:
                raise ValueError(f"Selection '{label}' has 0 atoms")
            for arr_name, arr in (("resids", resids),
                                  ("resnames", resnames),
                                  ("atom_names", atom_names),
                                  ("rmsf", rmsf_values)):
                if arr.size != n_atoms:
                    raise ValueError(
                        f"Selection '{label}': '{arr_name}' length ({arr.size}) "
                        f"does not match atom_indices length ({n_atoms})"
                    )

            sel_group = f"{group_name}/selection_{label}"
            if self.group_exists(sel_group):
                if not overwrite:
                    raise FileExistsError(
                        f"Selection '{label}' already exists under '{group_name}'"
                    )
                self.delete_group(sel_group)
            self.create_group(sel_group)

            records = [
                {
                    "atom_index": int(atom_indices[i]),
                    "resid": int(resids[i]),
                    "resname": str(resnames[i]),
                    "atom_name": str(atom_names[i]),
                    "rmsf": float(rmsf_values[i]),
                }
                for i in range(n_atoms)
            ]

            data = np.array(
                [json.dumps(r) for r in records],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

            meta = {
                "label": label,
                "selection_string": str(sel["selection_string"]),
                "n_atoms": str(n_atoms),
                "fitting_group": "" if fitting_group is None else str(fitting_group),
                "frame_begin": str(int(frame_begin)),
                "frame_end": str(int(frame_end)),
                "frame_step": str(int(frame_step)),
                "mean": f"{float(rmsf_values.mean()):.6f}",
                "min": f"{float(rmsf_values.min()):.6f}",
                "max": f"{float(rmsf_values.max()):.6f}",
                "median": f"{float(np.median(rmsf_values)):.6f}",
                "schema": "{atom_index, resid, resname, atom_name, rmsf}",
                "format": "json-per-row",
            }

            self.create_dataset(
                group_name=sel_group,
                dataset_name="atoms",
                data=data,
                metadata=meta,
            )

    def _iter_selections(self, group_name: str):
        """
        Iterates over selection labels within a specified group. Identifies
        keys starting with "selection_" and yields each label (suffix) in
        sorted order. Reserved subgroups (e.g. "statistics") are skipped
        because they do not carry the "selection_" prefix.

        :param group_name: The parent group (e.g. "rmsf").
        :type group_name: str
        :return: Generator yielding selection labels in sorted order.
        :rtype: Iterator[str]
        """
        root = self.file[group_name]
        labels = []
        for name in root.keys():
            if name.startswith("selection_"):
                labels.append(name[len("selection_"):])
        for label in sorted(labels):
            yield label

    def read_rmsf_data(self, *, group_name: str = "rmsf") -> Dict[str, Dict]:
        """
        Reads per-atom RMSF data for every selection under `group_name`.

        Returns a dict-of-dicts (not a single packed ndarray) because each
        selection can have a different atom count.

        Returns:
            {
                label: {
                    "selection_string": str,
                    "atom_indices":     np.ndarray[int],
                    "resids":           np.ndarray[int],
                    "resnames":         np.ndarray[str],
                    "atom_names":       np.ndarray[str],
                    "rmsf":             np.ndarray[float],
                },
                ...
            }

        :param group_name: Parent group name (default "rmsf").
        :raises ValueError: If the group is missing or no selections are found.
        """
        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        out: Dict[str, Dict] = {}
        for label in self._iter_selections(group_name):
            dset_path = f"{group_name}/selection_{label}/atoms"
            if dset_path not in self.file:
                continue
            dset = self.file[dset_path]
            attrs = dset.attrs
            records = [json.loads(row) for row in dset]
            if not records:
                continue

            entry: Dict = {
                "selection_string": str(attrs.get("selection_string", "")),
                "atom_indices": np.asarray([r["atom_index"] for r in records], dtype=int),
                "resids":       np.asarray([r["resid"] for r in records], dtype=int),
                "resnames":     np.asarray([r["resname"] for r in records], dtype=object),
                "atom_names":   np.asarray([r["atom_name"] for r in records], dtype=object),
            }

            # Merged datasets carry mean/std/n in place of a single rmsf value.
            if "rmsf" in records[0]:
                entry["rmsf"] = np.asarray([r["rmsf"] for r in records], dtype=float)
            else:
                entry["mean"] = np.asarray([r["mean"] for r in records], dtype=float)
                entry["std"]  = np.asarray([r["std"]  for r in records], dtype=float)
                entry["n"]    = np.asarray([r["n"]    for r in records], dtype=int)
            out[label] = entry

        if not out:
            raise ValueError("No RMSF selections found")
        return out

    def write_rmsf_statistics(self, *, group_name: str = "rmsf") -> None:
        """
        Builds RMSF summary statistics across atoms for each selection.

        Frame-range provenance (begin/end/step) is read from each selection's
        `atoms` dataset attrs (they were locked at write time), so no
        begin/end/step args are needed here.

        Per-selection record schema:
            {label, n_atoms, mean, min, max, median,
             argmax_atom_index, argmax_resid, argmax_resname, argmax_atom_name}

        Stored at:
            /<group_name>/statistics/table

        Data format:
            JSON-per-row
        """
        if not self.group_exists(group_name):
            raise ValueError(f"Group '{group_name}' does not exist")

        labels = list(self._iter_selections(group_name))
        if not labels:
            raise ValueError("No RMSF selections found")

        records = []
        frame_begin: int | None = None
        frame_end: int | None = None
        frame_step: int | None = None

        for label in labels:
            dset_path = f"{group_name}/selection_{label}/atoms"
            if dset_path not in self.file:
                continue
            dset = self.file[dset_path]
            attrs = dset.attrs

            if frame_begin is None:
                frame_begin = int(attrs.get("frame_begin", 0))
                frame_end = int(attrs.get("frame_end", 0))
                frame_step = int(attrs.get("frame_step", 1))

            rows = [json.loads(r) for r in dset]
            if not rows:
                continue

            # Merged datasets carry "mean" instead of "rmsf"; use whichever is present.
            rmsf_arr = np.asarray(
                [r.get("rmsf", r.get("mean")) for r in rows], dtype=float,
            )
            argmax = int(rmsf_arr.argmax())

            records.append({
                "label": label,
                "n_atoms": int(rmsf_arr.size),
                "mean": float(rmsf_arr.mean()),
                "min": float(rmsf_arr.min()),
                "max": float(rmsf_arr.max()),
                "median": float(np.median(rmsf_arr)),
                "argmax_atom_index": int(rows[argmax]["atom_index"]),
                "argmax_resid": int(rows[argmax]["resid"]),
                "argmax_resname": str(rows[argmax]["resname"]),
                "argmax_atom_name": str(rows[argmax]["atom_name"]),
            })

        if not records:
            raise ValueError("No RMSF values collected")

        stats_group = f"{group_name}/statistics"
        if self.group_exists(stats_group):
            self.delete_group(stats_group)
        self.create_group(stats_group)

        data = np.array(
            [json.dumps(r) for r in records],
            dtype=h5py.string_dtype("utf-8"),
        )

        self.create_dataset(
            group_name=stats_group,
            dataset_name="table",
            data=data,
            metadata={
                "n_rows": str(len(records)),
                "frame_begin": str(frame_begin if frame_begin is not None else 0),
                "frame_end": str(frame_end if frame_end is not None else 0),
                "frame_step": str(frame_step if frame_step is not None else 1),
                "end_inclusive": "True",
                "schema": "{label, n_atoms, mean, min, max, median, "
                          "argmax_atom_index, argmax_resid, argmax_resname, "
                          "argmax_atom_name}",
                "format": "json-per-row",
            },
        )

        self.add_group_metadata(
            group_name=stats_group,
            metadata={
                "completed": "True",
                "n_selections": str(len(records)),
                "description": "Aggregated RMSF statistics per selection",
            },
            overwrite=True,
        )

    def write_rmsf_data_to_csv(self, output_file: str | Path, *,
                               group_name: str = "rmsf",
                               is_merged: bool = False) -> None:
        """
        Writes per-atom RMSF data to a CSV file.

        Non-merged schema:
            selection, atom_index, resid, resname, atom_name, rmsf
        Merged schema:
            selection, atom_index, resid, resname, atom_name, mean, std

        :param output_file: Path of the CSV file to write.
        :param group_name: Parent group containing the RMSF data (default "rmsf").
        :param is_merged: If True, expects merged per-atom records that carry
            "mean" and "std" keys instead of "rmsf".
        """
        output_file = Path(output_file).expanduser().resolve()
        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            if not is_merged:
                writer.writerow([
                    "selection", "atom_index", "resid", "resname",
                    "atom_name", "rmsf",
                ])
            else:
                writer.writerow([
                    "selection", "atom_index", "resid", "resname",
                    "atom_name", "mean", "std",
                ])

            for label in self._iter_selections(group_name):
                dset_path = f"{group_name}/selection_{label}/atoms"
                if dset_path not in self.file:
                    continue
                dset = self.file[dset_path]
                for row in dset:
                    rec = json.loads(row)
                    if not is_merged:
                        writer.writerow([
                            label,
                            rec["atom_index"],
                            rec["resid"],
                            rec["resname"],
                            rec["atom_name"],
                            f"{rec['rmsf']:.6f}",
                        ])
                    else:
                        writer.writerow([
                            label,
                            rec["atom_index"],
                            rec["resid"],
                            rec["resname"],
                            rec["atom_name"],
                            f"{rec['mean']:.6f}",
                            f"{rec['std']:.6f}",
                        ])

    def write_rmsf_data_to_tsv(self, output_file: str | Path, *,
                               group_name: str = "rmsf",
                               is_merged: bool = False) -> None:
        """
        Writes per-atom RMSF data to a TSV file. Schema matches
        write_rmsf_data_to_csv but separated by tabs.
        """
        output_file = Path(output_file).expanduser().resolve()
        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            if not is_merged:
                writer.writerow([
                    "selection", "atom_index", "resid", "resname",
                    "atom_name", "rmsf",
                ])
            else:
                writer.writerow([
                    "selection", "atom_index", "resid", "resname",
                    "atom_name", "mean", "std",
                ])

            for label in self._iter_selections(group_name):
                dset_path = f"{group_name}/selection_{label}/atoms"
                if dset_path not in self.file:
                    continue
                dset = self.file[dset_path]
                for row in dset:
                    rec = json.loads(row)
                    if not is_merged:
                        writer.writerow([
                            label,
                            rec["atom_index"],
                            rec["resid"],
                            rec["resname"],
                            rec["atom_name"],
                            f"{rec['rmsf']:.6f}",
                        ])
                    else:
                        writer.writerow([
                            label,
                            rec["atom_index"],
                            rec["resid"],
                            rec["resname"],
                            rec["atom_name"],
                            f"{rec['mean']:.6f}",
                            f"{rec['std']:.6f}",
                        ])

    def write_rmsf_statistics_to_csv(self, output_file: str | Path, *,
                                     group_name: str = "rmsf") -> None:
        """
        Writes RMSF per-selection statistics to a CSV file.

        Schema:
            selection, n_atoms, mean, min, max, median,
            argmax_atom_index, argmax_resid, argmax_resname, argmax_atom_name
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "selection", "n_atoms", "mean", "min", "max", "median",
                "argmax_atom_index", "argmax_resid", "argmax_resname",
                "argmax_atom_name",
            ])
            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["n_atoms"],
                    f"{rec['mean']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                    f"{rec['median']:.6f}",
                    rec["argmax_atom_index"],
                    rec["argmax_resid"],
                    rec["argmax_resname"],
                    rec["argmax_atom_name"],
                ])

    def write_rmsf_statistics_to_tsv(self, output_file: str | Path, *,
                                     group_name: str = "rmsf") -> None:
        """
        Writes RMSF per-selection statistics to a TSV file. Schema matches
        write_rmsf_statistics_to_csv but separated by tabs.
        """
        output_file = Path(output_file).expanduser().resolve()
        dset = self.file[f"{group_name}/statistics/table"]

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "selection", "n_atoms", "mean", "min", "max", "median",
                "argmax_atom_index", "argmax_resid", "argmax_resname",
                "argmax_atom_name",
            ])
            for row in dset:
                rec = json.loads(row)
                writer.writerow([
                    rec["label"],
                    rec["n_atoms"],
                    f"{rec['mean']:.6f}",
                    f"{rec['min']:.6f}",
                    f"{rec['max']:.6f}",
                    f"{rec['median']:.6f}",
                    rec["argmax_atom_index"],
                    rec["argmax_resid"],
                    rec["argmax_resname"],
                    rec["argmax_atom_name"],
                ])

    def write_merged_interaction_modes_to_csv(self, output_file: str | Path, *, group_name: str, mode_name: str) -> None:
        """
        Writes merged interaction mode data to a CSV file. The data is retrieved based on
        the specified group and mode names, and it is saved in a structured CSV format containing
        interaction details and statistical frequency information.

        Reads from:
            /<group_name>/modes_merged/<mode_name>/table

        :param output_file: The file path where the CSV output will be saved. It can be a string
            or a Path-like object.
        :param group_name: The name of the interaction group from which to retrieve the data.
        :param mode_name: The specific mode name within the group for which the data will
            be written to the CSV.
        :return: This function does not return any value. All processed data will be written to
            the specified CSV file.
        """

        output_file = Path(output_file).expanduser().resolve()

        records = self.read_merged_interactions(
            group_name=group_name,
            mode_name=mode_name,
        )

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)

            writer.writerow([
                "interaction",
                "res1_resname", "res1_resid", "res1_chainid", "res1_segid",
                "res2_resname", "res2_resid", "res2_chainid", "res2_segid",
                "mean_frequency",
                "std_frequency",
                "n_files",
            ])

            for rec in records:
                (r1, r2, label) = rec["key"]

                writer.writerow([
                    label,
                    r1[0], r1[1], r1[2], r1[3],
                    r2[0], r2[1], r2[2], r2[3],
                    f"{rec['mean_frequency']:.6f}",
                    f"{rec['std_frequency']:.6f}",
                    rec["n_files"],
                ])

    def write_merged_interaction_modes_to_tsv(self, output_file: str | Path, *, group_name: str, mode_name: str) -> None:
        """
        Writes merged interaction modes to a TSV file. The method retrieves interaction
        records for a specified group and mode name, and writes detailed interaction
        data into a tab-separated file, including information about interacting residues,
        frequency statistics, and other interaction details.

        :param output_file: The path where the TSV file will be created.
        :type output_file: str | Path
        :param group_name: The name of the interaction group to be processed.
        :type group_name: str
        :param mode_name: The specific interaction mode within the group to be processed.
        :type mode_name: str
        :return: None
        """
        output_file = Path(output_file).expanduser().resolve()

        records = self.read_merged_interactions(
            group_name=group_name,
            mode_name=mode_name,
        )

        with open(output_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")

            writer.writerow([
                "interaction",
                "res1_resname", "res1_resid", "res1_chainid", "res1_segid",
                "res2_resname", "res2_resid", "res2_chainid", "res2_segid",
                "mean_frequency",
                "std_frequency",
                "n_files",
            ])

            for rec in records:
                (r1, r2, label) = rec["key"]

                writer.writerow([
                    label,
                    r1[0], r1[1], r1[2], r1[3],
                    r2[0], r2[1], r2[2], r2[3],
                    f"{rec['mean_frequency']:.6f}",
                    f"{rec['std_frequency']:.6f}",
                    rec["n_files"],
                ])
