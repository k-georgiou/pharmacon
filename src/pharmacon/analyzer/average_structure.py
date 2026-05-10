"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

A module for processing molecular dynamics trajectories, including operations
such as structure averaging, alignment, and frame extraction.

This module defines essential functionalities for analyzing molecular dynamics
simulations performed using the MDAnalysis library. It includes methods for
trajectory processing, averaging atom positions, calculating RMSD metrics, and
extracting specific frames from simulations for further analysis or visualization.
"""


import os
import warnings
import tempfile
import numpy as np
from pathlib import Path
from datetime import datetime
from MDAnalysis import AtomGroup, Universe
import MDAnalysis.transformations as trans
from MDAnalysis.analysis.align import rotation_matrix


from pharmacon.constants import __version__
from pharmacon.logger import get_logger, PharmaconLogger




__all__ = [
    "logger",
    "avg_st_process_trajectory",
    "extract_trajectory_frame",
]


warnings.filterwarnings("ignore")

logger: PharmaconLogger = get_logger(__name__)


def avg_st_process_trajectory(u: Universe, selection: AtomGroup,
                              start: int, stop: int, step: int = 1,
                              reference_frame: int = 0, memory_efficient: bool = True):
    """
    This function processes a trajectory from a molecular dynamics universe and calculates
    an averaged structure and Root Mean Square Deviation (RMSD) metrics. The selection of
    atoms is aligned to a reference frame and averaged over the specified range of frames.
    Two modes of operation are available based on the `memory_efficient` parameter. If
    set to `True`, the function processes the trajectory incrementally to reduce memory
    usage. If `False`, it constructs a full aligned structure array in memory.

    :param u: Universe object containing the molecular trajectory to analyze.
    :type u: Universe
    :param selection: Subset of atoms to be processed within the Universe object.
    :type selection: AtomGroup
    :param start: Starting frame index (inclusive) for trajectory processing.
    :type start: int
    :param stop: Ending frame index (inclusive) for trajectory processing.
    :type stop: int
    :param step: Step interval for iterating through frames, default is 1.
    :type step: int
    :param reference_frame: Frame index used as the initial reference for alignment,
        must be within the range of start and stop frames.
    :type reference_frame: int
    :param memory_efficient: Flag to indicate whether trajectory processing is conducted
        in a memory-efficient manner. Defaults to True.
    :type memory_efficient: bool
    :return: A tuple containing:
        - `avg_positions` (numpy.ndarray): Averaged 3D coordinates of the selected atoms.
        - `rmsd_to_avg` (numpy.ndarray): RMSD values to the averaged positions for each
          considered frame.
        - `rmsd_to_ref` (numpy.ndarray): RMSD values to the reference frame for each
          considered frame.
    :rtype: tuple
    :raises ValueError: If `reference_frame` is out of the specified range or exceeds
        the total number of trajectory frames, or if `step` is non-positive.
    :raises RuntimeError: If there are issues accessing the trajectory or during trajectory
        processing.
    """
    try:
        n_frames = len(u.trajectory)
    except Exception as e:
        raise RuntimeError(f"Failed to access trajectory: {e}")

    # reference frame validation
    if reference_frame < 0:
        raise ValueError("reference_frame must be >= 0")

    if reference_frame < start or reference_frame > stop:
        raise ValueError(
            f"reference_frame ({reference_frame}) must be between start ({start}) and stop ({stop})"
        )

    if reference_frame >= n_frames:
        raise ValueError(
            f"reference_frame ({reference_frame}) exceeds total frames ({n_frames})"
        )

    if step <= 0:
        raise ValueError("step must be > 0")

    u.trajectory[reference_frame]
    ref = selection.positions.copy()
    ref_com = ref.mean(axis=0)

    try:
        if memory_efficient:
            count = 0
            mean = np.zeros_like(ref, dtype=float)
            rmsd_to_ref = []

            for i, fr in enumerate(range(start, stop + 1, step)):
                u.trajectory[fr]
                mob = selection.positions.copy()
                mob_com = mob.mean(axis=0)

                R, rmsd_val = rotation_matrix(mob - mob_com, ref - ref_com)
                mob_aligned = (mob - mob_com) @ R + ref_com

                count += 1
                delta = mob_aligned - mean
                mean += delta / count

                rmsd_to_ref.append(float(rmsd_val))
                logger.debug(f"Averaging coordinates processed frame: {fr}")

            avg_positions = mean
            logger.info("Averaging coordinates completed.")

            rmsd_to_avg = []
            for i, fr in enumerate(range(start, stop + 1, step)):
                u.trajectory[fr]
                mob = selection.positions.copy()
                mob_com = mob.mean(axis=0)

                R, _ = rotation_matrix(mob - mob_com, ref - ref_com)
                mob_aligned = (mob - mob_com) @ R + ref_com

                diff = mob_aligned - avg_positions
                rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
                rmsd_to_avg.append(float(rmsd))
                logger.debug(f"RMSD to avg processed frame: {fr}")

            return avg_positions, np.asarray(rmsd_to_avg), np.asarray(rmsd_to_ref)

        else:
            aligned_positions_per_frame = []
            rmsd_to_ref = []

            for i, fr in enumerate(range(start, stop + 1, step)):
                u.trajectory[fr]
                mob = selection.positions.copy()
                mob_com = mob.mean(axis=0)

                R, rmsd_val = rotation_matrix(mob - mob_com, ref - ref_com)
                mob_aligned = (mob - mob_com) @ R + ref_com

                aligned_positions_per_frame.append(mob_aligned)
                rmsd_to_ref.append(float(rmsd_val))
                logger.debug(f"Averaging coordinates processed frame: {fr}")

            all_pos = np.stack(aligned_positions_per_frame, axis=0)
            avg_positions = all_pos.mean(axis=0)

            diff = all_pos - avg_positions[None, :, :]
            rmsd_to_avg = np.sqrt(np.mean(np.sum(diff * diff, axis=2), axis=1))

            return avg_positions, rmsd_to_avg, np.asarray(rmsd_to_ref, dtype=float)

    except Exception as e:
        raise RuntimeError(f"Failed during trajectory averaging: {e}")


def extract_trajectory_frame(u: Universe, ag: AtomGroup, selection: str, frame_idx: int, *,
                             output_file_path: str | Path, output_format: str,
                             temporary_directory: str | Path | None = None,
                             title: str = "Most representative structure saved by Pharmacon") -> None:
    """
    Extracts a specific frame from a molecular dynamics trajectory and saves it
    to a specified file format. This operation involves frame validation, optional
    trajectory transformations, atom selection, and saving the output with additional
    metadata.

    :param u: MDAnalysis Universe object representing the molecular dynamics system
        and trajectory.
    :param ag: AtomGroup object representing the subset of atoms to be processed.
    :param selection: A selection string for selecting specific atoms in the system.
    :param frame_idx: Index of the trajectory frame to extract. Must be greater
        than or equal to 0 and less than the number of trajectory frames.
    :param output_file_path: Path to the output file where the extracted frame
        will be written.
    :param output_format: Desired output file format (e.g., pdb, gro, crd).
    :param temporary_directory: Optional temporary directory for creating
        intermediate files.
    :param title: Title or description to include in output file headers or remarks.
    :return: None
    """

    logger.debug(f"Extracting frame {frame_idx} from trajectory")
    today = datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    # --- frame index validation ---
    n_frames = len(u.trajectory)
    if frame_idx < 0:
        raise ValueError("frame_idx must be >= 0")
    if frame_idx >= n_frames:
        raise ValueError(f"frame_idx ({frame_idx}) exceeds trajectory length ({n_frames})")

    try:
        u.trajectory[frame_idx]

        workflow = [
            trans.unwrap(u.atoms),
            trans.center_in_box(ag, wrap=True),
        ]
        u.trajectory.add_transformations(*workflow)

        sel = u.select_atoms(selection)

    except Exception as e:
        raise RuntimeError(f"Failed preparing trajectory for frame extraction: {e}")

    try:
        with tempfile.NamedTemporaryFile(
                dir=temporary_directory,
                delete=False,
                suffix=f".{output_format}",
        ) as tmp:
            tmp_name = tmp.name

        sel.write(tmp_name)

        with open(output_file_path, "w") as out_f:
            if output_format == "pdb":
                out_f.write(f"TITLE     {title}\n")
                out_f.write(f"REMARK    Version: {__version__}\n")
                out_f.write(f"REMARK    DateTime {today}\n")
                out_f.write(f"REMARK    Frame {frame_idx}\n")
            elif output_format == "gro":
                out_f.write(f"{title}-Frame{frame_idx}-V_{__version__}\n")
            elif output_format == "crd":
                out_f.write(f"* {title}-Frame{frame_idx}-V_{__version__}\n")

            with open(tmp_name, "r") as tmp_f:
                for i, line in enumerate(tmp_f):
                    if output_format == "pdb" and (
                            line.startswith("TITLE") or line.startswith("REMARK")
                    ):
                        continue
                    elif output_format == "gro" and i == 0:
                        continue
                    out_f.write(line)

        os.remove(tmp_name)

    except Exception as e:
        raise RuntimeError(f"Failed to extract and write frame {frame_idx}: {e}")
