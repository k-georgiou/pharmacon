"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Functions to calculate various types of distances between given pairs of
atom groups using specified methods.

Provides functionality to compute distances such as center of mass,
center of geometry, minimum, and maximum inter-atomic distances for
provided atom groups, with optional support for periodic boundary
conditions.
"""


import numpy as np
from typing import List, Tuple

from MDAnalysis import AtomGroup
from MDAnalysis.lib.distances import distance_array





__all__ = [
    "calculate_frame_distances",
]


def calculate_frame_distances(atom_groups1: List[AtomGroup],
                              atom_groups2: List[AtomGroup],
                              methods: List[str],
                              labels: List[str],
                              box: np.ndarray | None) -> List[Tuple[str, str, float]]:
    """
    Calculates distances between pairs of atom groups using specified methods.

    Supported methods:
      - "com" : center of mass distance
      - "cog" : center of geometry distance
      - "min" : minimum inter-atomic distance
      - "max" : maximum inter-atomic distance

    :param atom_groups1: First AtomGroup list.
    :param atom_groups2: Second AtomGroup list.
    :param methods: Distance methods ("com", "cog", "min", "max").
    :param labels: Semantic labels for each distance (method-independent).
    :param box: Unit cell dimensions for PBC, or None.
    :return: List of (label, method, distance) tuples.
    """

    if not (
        len(atom_groups1)
        == len(atom_groups2)
        == len(methods)
        == len(labels)
    ):
        raise ValueError(
            "atom_groups1, atom_groups2, methods, and labels must have equal length"
        )

    out: List[Tuple[str, str, float]] = []

    for ag1, ag2, method, label in zip(atom_groups1, atom_groups2, methods, labels):
        if ag1.n_atoms == 0:
            raise RuntimeError(
                f"[EMPTY GROUP1] label='{label}' | method='{method}' | "
            )

        if ag2.n_atoms == 0:
            raise RuntimeError(
                f"[EMPTY GROUP2] label='{label}' | method='{method}' | "
            )

        method = method.lower().strip()

        if method == "com":
            p1 = ag1.center_of_mass(wrap=False)
            p2 = ag2.center_of_mass(wrap=False)
            d = distance_array(
                p1.reshape(1, 3),
                p2.reshape(1, 3),
                box=box,
            )[0, 0]

        elif method == "cog":
            p1 = ag1.center_of_geometry(wrap=False)
            p2 = ag2.center_of_geometry(wrap=False)
            d = distance_array(
                p1.reshape(1, 3),
                p2.reshape(1, 3),
                box=box,
            )[0, 0]

        elif method in {"min", "max"}:
            A = ag1.positions
            B = ag2.positions
            D = distance_array(A, B, box=box)
            d = float(np.min(D) if method == "min" else np.max(D))

        else:
            raise RuntimeError(f"Unknown distance method '{method}'")

        out.append((label, method, float(d)))

    return out

