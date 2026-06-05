"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Geometric computations and selection utilities for molecular dynamics analysis.

This module provides methods for calculating angles, dihedral angles, and parsing
selection specifications for molecular dynamics analyses using MDAnalysis. It
integrates utility functions and handles edge cases with custom checks for atom
selection counts and zero-length vectors. The module also supports periodic
boundary conditions as required for molecular simulations.
"""


import re
import numpy as np
from MDAnalysis.core.groups import Atom
from MDAnalysis import AtomGroup, Universe
from typing_extensions import Final, Dict, List, Tuple, Any
from MDAnalysis.lib.distances import calc_angles, calc_dihedrals, minimize_vectors




__all__ = [
    "GeometrySpecError",
    "AXIS_VECTORS",
    "calculate_angle_between_vectors",
    "calculate_angle_3_atoms",
    "calculate_dihedral_angle",
    "create_mda_angle_selections",
    "calculate_frame_angles",
]


class GeometrySpecError(ValueError):
    """
    Represents an error specific to geometry specifications.

    This exception is derived from the ValueError class. It is used to
    indicate issues or inconsistencies within geometric specifications
    or configurations that deviate from expected or required standards.

    """
    pass


AXIS_VECTORS: Final[Dict[str, np.ndarray]] = {
    "x-axis": np.array([1.0, 0.0, 0.0], dtype=float),
    "y-axis": np.array([0.0, 1.0, 0.0], dtype=float),
    "z-axis": np.array([0.0, 0.0, 1.0], dtype=float),
}

_ZERO_TOL: Final[float] = 1e-12

_NAME_TOKEN_RE = re.compile(r"\bname\s+([A-Za-z0-9_]+)\b", re.IGNORECASE)


def _displacement(a1: Atom, a2: Atom, box: np.ndarray | None) -> np.ndarray:
    """
    Calculates the displacement vector between two atoms, considering periodic
    boundary conditions if a box is provided.

    The function computes the direct vector difference between the positions
    of two atoms. If a periodic boundary condition is defined by a box,
    the vector is adjusted to account for the minimum image convention.

    :param a1: The first atom, represented as an Atom object.
    :param a2: The second atom, represented as an Atom object.
    :param box: A 2D numpy array representing the periodic boundary
        conditions, or None if no boundary conditions are applied.
    :return: A numpy array representing the displacement vector.
    :rtype: np.ndarray
    """
    v = (a2.position - a1.position).astype(float, copy=False)
    if box is not None:
        v = minimize_vectors(v[None, :], box)[0]
    return v


def calculate_angle_between_vectors(a1: Atom, a2: Atom, a3: Atom, a4: Atom, box: np.ndarray | None = None) -> float:
    """
    Calculates the angle (in degrees) between two vectors defined by four atomic positions.

    This function computes the displacement vectors between two pairs of atoms and then
    determines the angle between these displacement vectors. The angle is calculated in
    degrees and considers periodic boundary conditions if a simulation box is provided.

    :param a1: The first atom in the first pair defining the first vector.
    :type a1: Atom
    :param a2: The second atom in the first pair defining the first vector.
    :type a2: Atom
    :param a3: The first atom in the second pair defining the second vector.
    :type a3: Atom
    :param a4: The second atom in the second pair defining the second vector.
    :type a4: Atom
    :param box: A 3x3 matrix representing the periodic simulation box (optional).
    :type box: numpy.ndarray | None
    :return: The angle between the two vectors in degrees.
    :rtype: float
    """
    v1 = _displacement(a1, a2, box)
    v2 = _displacement(a3, a4, box)

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < _ZERO_TOL or n2 < _ZERO_TOL:
        raise ValueError("Zero-length vector encountered")

    cos_theta = float(np.dot(v1, v2) / (n1 * n2))
    cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_theta)))


def calculate_angle_3_atoms(a1: Atom, a2: Atom, a3: Atom, box: np.ndarray | None = None) -> float:
    """
    Calculate the angle formed by three atoms.

    This function computes the angle, in degrees, between three atoms (a1, a2, and a3). The atoms
    are represented as instances of the `Atom` class. The optional parameter `box` allows for
    calculations considering periodic boundary conditions if applicable.

    :param a1: The first atom in the angle calculation.
    :type a1: Atom
    :param a2: The second atom in the angle calculation, serving as the vertex of the angle.
    :type a2: Atom
    :param a3: The third atom in the angle calculation.
    :type a3: Atom
    :param box: A 3x3 numpy array representing the periodic boundary conditions box. If None, no
        periodic boundary conditions are applied.
    :type box: np.ndarray or None
    :return: The angle formed by the three atoms in degrees.
    :rtype: float
    """
    angle_rad = calc_angles(
        a1.position[None, :],
        a2.position[None, :],
        a3.position[None, :],
        box=box,
    )[0]
    return float(np.degrees(angle_rad))


def calculate_dihedral_angle(a1: Atom, a2: Atom, a3: Atom, a4: Atom, box: np.ndarray | None = None) -> float:
    """
    Calculates the dihedral angle between four atoms in space. The dihedral angle measures the
    angle between the two planes formed by the first three atoms and the last three atoms.
    The calculation is performed in degrees and can optionally account for periodic boundary
    conditions using a provided box.

    :param a1: The first atom in the sequence, used in calculating the dihedral angle.
    :type a1: Atom
    :param a2: The second atom in the sequence, used in calculating the dihedral angle.
    :type a2: Atom
    :param a3: The third atom in the sequence, used in calculating the dihedral angle.
    :type a3: Atom
    :param a4: The fourth atom in the sequence, used in calculating the dihedral angle.
    :type a4: Atom
    :param box: Optional. An array defining the periodic boundary conditions. If None, no
        boundary conditions are applied.
    :type box: np.ndarray | None
    :return: The dihedral angle in degrees as a floating-point number.
    :rtype: float
    """
    dih_rad = calc_dihedrals(
        a1.position[None, :],
        a2.position[None, :],
        a3.position[None, :],
        a4.position[None, :],
        box=box,
    )[0]
    return float(np.degrees(dih_rad))


def _parse_vector_spec(u: Universe | AtomGroup, spec: str) -> Tuple[str, Any]:
    """
    Parses a vector specification and determines its type, either as a predefined
    axis vector or based on a selection of atoms from the universe. The function
    validates the given specification and returns a tuple representing the type of
    vector and its associated data.

    :param u: An instance of Universe or AtomGroup from which the atom selection
        will be made.
    :param spec: The string specification representing the vector. This can be a
        predefined axis vector keyword or an atom selection string.
    :return: A tuple containing the type of the vector as a string ("axis" or
        "atoms") and the corresponding data (axis vector or selected atom group).
    :rtype: Tuple[str, Any]
    :raises GeometrySpecError: If the atom selection does not resolve to exactly
        two atoms.
    """
    spec_norm = spec.strip().lower()

    # axis vector
    if spec_norm in AXIS_VECTORS:
        return ("axis", AXIS_VECTORS[spec_norm])

    # MDAnalysis selection → must be exactly 2 atoms
    ag = u.select_atoms(spec)
    if ag.n_atoms != 2:
        raise GeometrySpecError(
            f"Vector selection must resolve to exactly 2 atoms, got {ag.n_atoms}: '{spec}'"
        )
    return ("atoms", ag)


def create_mda_angle_selections(u: Universe | AtomGroup, selection_str: str) -> Dict[str, Any]:
    """
    Analyze and parse a molecular dynamic selection string into angle or dihedral
    geometry components, or evaluate a vector-to-vector angle interpretation.

    :param u: MDAnalysis Universe or AtomGroup, which provides the atomic data for
        parsing the selection string.
    :param selection_str: A formatted string specifying the geometry selection,
        such as atoms involved in an angle or dihedral, or vectors for angles.
    :return: A dictionary containing details of the parsed geometry component,
        including its type ('vector-angle', 'angle', 'dihedral'), and the associated
        atoms or vectors as applicable.
    :raises GeometrySpecError: If the selection string contains more than one '->'
        symbol, or if the atom selection fails to resolve to 3 atoms (angle) or 4
        atoms (dihedral).
    """
    arrow_count = selection_str.count("->")
    if arrow_count > 1:
        raise GeometrySpecError("Only one '->' is allowed.")

    # VECTOR–VECTOR ANGLE
    if arrow_count == 1:
        left, right = map(str.strip, selection_str.split("->", 1))
        left_type, left_obj = _parse_vector_spec(u, left)
        right_type, right_obj = _parse_vector_spec(u, right)
        return {"type": "vector-angle", "left": (left_type, left_obj), "right": (right_type, right_obj)}

    # ATOM-BASED ANGLE / DIHEDRAL
    ag = u.select_atoms(selection_str)

    if ag.n_atoms == 3:
        return {"type": "angle", "atoms": ag, "spec_str": selection_str}

    if ag.n_atoms == 4:
        return {"type": "dihedral", "atoms": ag, "spec_str": selection_str}

    raise GeometrySpecError(
        f"Selection must resolve to 3 atoms (angle) or 4 atoms (dihedral), got {ag.n_atoms}"
    )


def _ordered_names_from_selection_string(selection_str: str) -> List[str]:
    """
    Extracts and orders names from a selection string based on unique occurrences. The returned
    list contains uppercase names in the order they first appear in the selection string.

    :param selection_str: The input string containing name tokens to extract and order.
    :type selection_str: str
    :return: A list of unique, uppercase names extracted in their first occurrence order.
    :rtype: List[str]

    Example:
      "(resname UNA and name C14) or (resname UNA and name C09) or ..."
    -> ["C14","C09","N08","N07"]  (unique, in order)

    """
    names = [m.group(1).upper() for m in _NAME_TOKEN_RE.finditer(selection_str)]
    out, seen = [], set()
    for n in names:
        if n not in seen:
            out.append(n)
            seen.add(n)
    return out


def _reorder_ag_by_names(ag: AtomGroup, ordered_names: List[str]) -> AtomGroup:
    """
    Reorders the atoms of an AtomGroup based on the specified order of atom names.

    This function checks if the number of atoms in the selection matches the
    length of the ordered name list. If the counts do not match, or if duplicate
    atom names or missing atoms are found, an exception will be raised. The atoms
    are then reordered based on the provided ordered names.

    :param ag: The AtomGroup whose atoms need to be reordered.
    :type ag: AtomGroup
    :param ordered_names: A list of atom names defining the desired order for the
        atoms in `ag`.
    :type ordered_names: List[str]
    :return: A new AtomGroup instance with its atoms reordered based on
        `ordered_names`.
    :rtype: AtomGroup
    :raises GeometrySpecError: If the number of atoms in the AtomGroup does not
        match the length of the `ordered_names` list, if duplicate names are found
        in the selection, or if some atom names in `ordered_names` are missing
        in the AtomGroup.
    """
    if ag.n_atoms != len(ordered_names):
        raise GeometrySpecError(
            f"Cannot reorder: selection has {ag.n_atoms} atoms but inferred order has {len(ordered_names)}."
        )

    sel_names = [a.name.upper() for a in ag.atoms]
    if len(set(sel_names)) != len(sel_names):
        raise GeometrySpecError(
            f"Ambiguous dihedral selection: duplicate atom names in selection {sel_names}. "
            f"Constrain selection to a single residue/ligand (add resid/segid/chainid)."
        )

    by_name = {a.name.upper(): a for a in ag.atoms}
    missing = [n for n in ordered_names if n not in by_name]
    if missing:
        raise GeometrySpecError(
            f"Cannot reorder dihedral: missing atoms {missing}. Got {sel_names}."
        )

    ordered_atoms = [by_name[n] for n in ordered_names]
    return ag.universe.atoms[[a.ix for a in ordered_atoms]]


def calculate_frame_angles(specs: List[Dict[str, Any]],
                           labels: List[str],
                           box: np.ndarray | None) -> List[Tuple[str, str, float]]:
    """
    Calculate angles based on provided specifications and labels. This function computes vector–vector angles,
    3-atom angles, or 4-atom dihedral angles depending on the configuration defined in the input specifications.
    The results are returned as a list of tuples containing the label, angle type, and the computed angle value.

    :param specs: List of dictionaries defining the calculation type and necessary data such as atoms, vectors,
        or selection strings.
    :type specs: List[Dict[str, Any]]

    :param labels: List of string labels corresponding to each specification for identifying the computed result.
    :type labels: List[str]

    :param box: Optional array or None representing periodic boundary conditions for calculating displacement vectors
        between atoms.
    :type box: numpy.ndarray | None

    :return: A list of tuples where each tuple contains a label, the angle type (e.g., "vector-angle", "angle", "dihedral"),
        and the calculated angle in degrees.
    :rtype: List[Tuple[str, str, float]]
    """
    if len(specs) != len(labels):
        raise ValueError("specs and labels must have the same length")

    results = []
    for spec, label in zip(specs, labels):
        kind = spec["type"]

        # VECTOR–VECTOR ANGLE
        if kind == "vector-angle":
            (l_type, l_obj) = spec["left"]
            (r_type, r_obj) = spec["right"]

            # left vector
            if l_type == "axis":
                v1 = np.asarray(l_obj, dtype=float)
            else:
                a1, a2 = l_obj.atoms
                v1 = _displacement(a1, a2, box)

            # right vector
            if r_type == "axis":
                v2 = np.asarray(r_obj, dtype=float)
            else:
                a3, a4 = r_obj.atoms
                v2 = _displacement(a3, a4, box)

            n1 = np.linalg.norm(v1)
            n2 = np.linalg.norm(v2)
            if n1 < _ZERO_TOL or n2 < _ZERO_TOL:
                raise ValueError(f"Zero-length vector for label '{label}'")

            cos_theta = float(np.dot(v1, v2) / (n1 * n2))
            cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
            value = float(np.degrees(np.arccos(cos_theta)))

            results.append((label, "vector-angle", value))

        # 3-ATOM ANGLE
        elif kind == "angle":
            ag: AtomGroup = spec["atoms"]

            # Dynamic ordering based on selection string
            spec_str = spec.get("spec_str", "")
            inferred = _ordered_names_from_selection_string(spec_str)

            # Only reorder if we inferred exactly 3 unique names
            if len(inferred) == 3:
                ag = _reorder_ag_by_names(ag, inferred)

            # Now guaranteed deterministic order
            a1, a2, a3 = ag.atoms
            value = calculate_angle_3_atoms(a1, a2, a3, box=box)

            results.append((label, "angle", value))

        # 4-ATOM DIHEDRAL
        elif kind == "dihedral":
            ag: AtomGroup = spec["atoms"]

            # Dynamic ordering based on the selection string, when possible
            spec_str = spec.get("spec_str", "")
            inferred = _ordered_names_from_selection_string(spec_str)

            # Only attempt reorder if we inferred exactly 4 unique names
            if len(inferred) == 4:
                ag = _reorder_ag_by_names(ag, inferred)

            # If we couldn't infer 4 names, we fall back to topology order (old behavior)
            # You can also choose to raise instead if you want strictness.

            a1, a2, a3, a4 = ag.atoms
            value = calculate_dihedral_angle(a1, a2, a3, a4, box=box)
            results.append((label, "dihedral", value))

        else:
            raise RuntimeError(f"Unknown geometry type '{kind}'")

    return results

