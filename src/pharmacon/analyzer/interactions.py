"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Module for analyzing atomic interactions, including hydrophobic and hydrogen bond
contacts, between molecular structures utilizing atomistic data.

This module provides utilities to classify atom types, retrieve chain or segment
information, and detect specific interactions like hydrophobic contacts between
molecular groups. The interaction analysis applies distance criteria and optionally
accounts for periodic boundary conditions (PBC).

It supports efficient serialization of detected interactions into a schema-stable,
flat tuple format for downstream analysis or storage.
"""

import numpy as np
import MDAnalysis as Mda
from typing import Tuple, List
from MDAnalysis import AtomGroup
from collections import defaultdict
from MDAnalysis.core.groups import Atom
from MDAnalysis.lib.distances import distance_array, calc_bonds, calc_angles


from pharmacon.constants.smarts import WATER_RESIDUES_MDA_SELECTION_STR, AA3, BB_RE, WATER_RESIDUES
from pharmacon.logger import get_logger, PharmaconLogger




__all__ = [
    "logger",
    "calculate_hydrophobic_contacts",
    "calculate_hydrogen_bond_contacts",
    "calculate_first_degree_water_bridge_contacts",
    "calculate_ionic_contacts",
    "calculate_halogen_bond_contacts",
    "calculate_pi_stacking_contacts",
    "calculate_pi_cation_contacts",
    "calculate_metal_contacts",
    "interactions_process_frame",
    "deduplicate_interactions",
    "extract_mode_key",
    "hbonds_process_frame",
]


logger: PharmaconLogger = get_logger(__name__)


def _get_chain(atom_group: AtomGroup, idx: int) -> str:
    """
    Fetches the chain identifier from the given atom group at the specified index.

    The function attempts to retrieve the chain identifier corresponding to the
    provided index from the `chainIDs` attribute of the `atom_group`, if it exists.
    If the attribute is not present, or the value at the index is not set, an empty
    string is returned.

    :param atom_group:
        A data structure, expected to contain a `chainIDs` attribute,
        holding chain identifiers for the atoms.
    :param idx:
        Integer index specifying which chain identifier to retrieve.
    :return:
        The chain identifier as a string from the specified index, or an
        empty string if unavailable.
    """
    if hasattr(atom_group, "chainIDs"):
        return atom_group.chainIDs[idx] or ""
    return ""


def _get_segid(atom_group: AtomGroup, idx: int) -> str:
    """
    Retrieve the segment ID (segid) at a specified index for the given atom group. The segid is
    obtained from the `segids` attribute of the `atom_group` if it exists. If the attribute does
    not exist, or the value at the given index is None or empty, an empty string is returned.

    :param atom_group: Group of atoms that potentially contains segment IDs.
    :type atom_group: AtomGroup
    :param idx: Index of the segment ID to be retrieved.
    :type idx: int
    :return: Segment ID at the given index or an empty string if not available.
    :rtype: str
    """
    if hasattr(atom_group, "segids"):
        return atom_group.segids[idx] or ""
    return ""


def _classify_type(res_name: str, atom_name: str, *, ring_mode: bool = False) -> str:
    """
    Classifies the type of a given molecular component based on residue name and atom name.

    This function determines the classification of a molecular component as either
    water, backbone, side chain, or unknown. It considers residue names, atom names,
    and optionally a ring mode flag to differentiate between backbone and side chain
    when identifying protein residues.

    :param res_name: Residue name of the molecular component. Expected to be a string.
    :param atom_name: Atom name of the molecular component. Expected to be a string.
    :param ring_mode: Boolean flag to indicate whether ring mode classification should
        be applied for protein residues. Defaults to False.
    :return: A string representing the classification of the molecular component. Possible
        values are:
        - "W": Water
        - "BB": Backbone
        - "SC": Side chain
        - "UN": Unknown
    """
    r = (res_name or "").strip().upper()
    a = (atom_name or "").strip().upper()

    # Water first (critical!)
    if r in WATER_RESIDUES:
        return "W"

    # Protein residues
    if r in AA3:
        if ring_mode:
            return "SC"
        return "BB" if BB_RE.match(a) else "SC"

    return "UN"


def _is_hydrogen(atom) -> bool:
    """Return True if the atom's element is hydrogen."""
    return (getattr(atom, "element", "") or "").upper() == "H"


def _heavy_neighbors(atom) -> list:
    """Return non-hydrogen bonded neighbors of an atom."""
    return [nb for nb in atom.bonded_atoms if not _is_hydrogen(nb)]


def _ring_normal_cross(ring) -> np.ndarray:
    """Compute the unit normal vector of an aromatic ring via cross product."""
    coords = ring.positions.astype(np.float64)
    if coords.shape[0] < 2:
        return np.array([0.0, 0.0, 1.0], dtype=float)
    c = ring.center_of_geometry().astype(np.float64)
    v1 = None
    idx_v1 = 0
    for i in range(min(8, coords.shape[0])):
        d = coords[i] - c
        if np.linalg.norm(d) > 1e-12:
            v1 = d
            idx_v1 = i
            break
    if v1 is None:
        return np.array([0.0, 0.0, 1.0], dtype=float)
    n = None
    for j in range(idx_v1 + 1, min(idx_v1 + 8, coords.shape[0])):
        v2 = coords[j] - c
        if np.linalg.norm(v2) <= 1e-12:
            continue
        cross = np.cross(v1, v2)
        if np.linalg.norm(cross) > 1e-12:
            n = cross
            break
    if n is None:
        for a in range(coords.shape[0] - 1):
            for b in range(a + 1, coords.shape[0]):
                cross = np.cross(coords[a] - c, coords[b] - c)
                if np.linalg.norm(cross) > 1e-12:
                    n = cross
                    break
            if n is not None:
                break
    if n is None or not np.all(np.isfinite(n)):
        return np.array([0.0, 0.0, 1.0], dtype=float)
    nrm = np.linalg.norm(n)
    if nrm == 0.0 or not np.isfinite(nrm):
        return np.array([0.0, 0.0, 1.0], dtype=float)
    return n / nrm


def _ring_meta_str(ag) -> Tuple[str, str, str, str, str, str, str, str, str]:
    """Generate comma-joined metadata strings for atoms in a ring AtomGroup."""
    indices = ",".join(str(int(a.index)) for a in ag.atoms)
    names = ",".join(a.name for a in ag.atoms)
    ids = ",".join(str(int(a.id)) for a in ag.atoms)
    elems = ",".join(a.element for a in ag.atoms)
    atom_types = ",".join(
        _classify_type(a.resname, a.name, ring_mode=True) for a in ag.atoms
    )
    resn = ",".join(a.resname for a in ag.atoms)
    resid = ",".join(str(int(a.resid)) for a in ag.atoms)
    chains = ",".join(
        (a.chainID or "") for a in ag.atoms
    ) if hasattr(ag.atoms[0], "chainID") else ""
    segids = ",".join(
        (a.segid or "") for a in ag.atoms
    ) if hasattr(ag.atoms[0], "segid") else ""
    return (indices, names, ids, elems, atom_types, resn, resid, chains, segids)


def calculate_hydrophobic_contacts(group1: AtomGroup,
                                   group2: AtomGroup,
                                   box: np.ndarray | None,
                                   cutoff_heavy: float = 3.6,
                                   cutoff_hydrogen: float = 3.2) -> Tuple[Tuple, ...]:
    """
    Detect hydrophobic contacts between two atom groups using distance-based criteria
    with periodic boundary conditions (PBC).

    A hydrophobic contact is defined as an atom pair whose minimum-image distance
    satisfies a cutoff that depends on atom types:

        - heavy–heavy pairs use ``cutoff_heavy`` (Å)
        - heavy–hydrogen pairs use ``cutoff_hydrogen`` (Å)
        - hydrogen–hydrogen pairs are excluded

    Each detected contact is returned as a flat, schema-stable tuple designed for
    efficient downstream serialization and analysis.

    --------------------------------------------------------------------
    Interaction record format
    --------------------------------------------------------------------

    Each contact is represented as:

        (
            "HYDROPHOBIC",

            atom1_index, atom1_name, atom1_id, atom1_element, atom1_type,
            atom1_resname, atom1_resid, atom1_chain, atom1_segid,

            atom2_index, atom2_name, atom2_id, atom2_element, atom2_type,
            atom2_resname, atom2_resid, atom2_chain, atom2_segid,

            distance,               # float, Å
            is_hydrogen_involved    # bool
        )

    where:

        - ``is_hydrogen_involved`` is ``True`` if either atom in the pair is hydrogen
        - ``atom*_type`` is a coarse structural classification:
            * ``"BB"`` backbone atom (protein)
            * ``"SC"`` side-chain or ring atom
            * ``"UN"`` unknown / non-standard residue

    --------------------------------------------------------------------
    Parameters
    --------------------------------------------------------------------

    :param group1:
        First atom group participating in hydrophobic contact detection.
    :type group1: AtomGroup

    :param group2:
        Second atom group participating in hydrophobic contact detection.
    :type group2: AtomGroup

    :param box:
        Periodic boundary condition box used for minimum-image distance
        calculations (shape ``(3, 3)`` or equivalent).
    :type box: numpy.ndarray | None

    :param cutoff_heavy:
        Distance cutoff (Å) applied to heavy–heavy pairs.
    :type cutoff_heavy: float

    :param cutoff_hydrogen:
        Distance cutoff (Å) applied to heavy–hydrogen pairs.
    :type cutoff_hydrogen: float

    --------------------------------------------------------------------
    Returns
    --------------------------------------------------------------------

    :return:
        Tuple of hydrophobic interaction records. If no contacts are found,
        an empty tuple is returned.
    :rtype:
        Tuple[Tuple, ...]
    """

    if len(group1) == 0 or len(group2) == 0:
        return ()

    # Distances (Å) with PBC
    dists = distance_array(group1.positions, group2.positions, box=box)

    # Element flags
    g1_isH = (group1.elements == "H")
    g2_isH = (group2.elements == "H")

    # Pairwise cutoff matrix:
    #   use hydrogen cutoff if either partner is H
    #   otherwise heavy cutoff
    pairwise_cutoffs = np.where(
        g1_isH[:, None] | g2_isH[None, :],
        cutoff_hydrogen,
        cutoff_heavy
    )

    # Exclude H–H pairs and apply distance mask
    mask = (dists < pairwise_cutoffs) & ~(g1_isH[:, None] & g2_isH[None, :])

    i_indices, j_indices = np.where(mask)

    if i_indices.size == 0:
        return ()

    # Pre-extract arrays to avoid repeated attribute access in the loop
    g1_indices = group1.indices
    g1_names = group1.names
    g1_ids = group1.ids
    g1_elements = group1.elements
    g1_resnames = group1.resnames
    g1_resids = group1.resids

    g2_indices = group2.indices
    g2_names = group2.names
    g2_ids = group2.ids
    g2_elements = group2.elements
    g2_resnames = group2.resnames
    g2_resids = group2.resids

    contacts = []
    for i, j in zip(i_indices, j_indices):
        contacts.append((
            "HYDROPHOBIC",

            int(g1_indices[i]),
            g1_names[i],
            int(g1_ids[i]),
            g1_elements[i],
            _classify_type(str(g1_resnames[i]), str(g1_names[i])),
            g1_resnames[i],
            int(g1_resids[i]),
            _get_chain(group1, i),
            _get_segid(group1, i),

            int(g2_indices[j]),
            g2_names[j],
            int(g2_ids[j]),
            g2_elements[j],
            _classify_type(str(g2_resnames[j]), str(g2_names[j])),
            g2_resnames[j],
            int(g2_resids[j]),
            _get_chain(group2, j),
            _get_segid(group2, j),

            float(dists[i, j]),
            bool(g1_isH[i] or g2_isH[j])
        ))

    return tuple(contacts)


def calculate_hydrogen_bond_contacts(group1_acceptors: AtomGroup, group2_acceptors: AtomGroup,
                                     group1_donors: AtomGroup, group2_donors: AtomGroup,
                                     box: np.ndarray | None,
                                     min_angleDHA: float = 120.0,
                                     min_angleHAX: float = 90.0,
                                     cutoff_loose: float = 6.0,
                                     cutoff: float = 2.9) -> Tuple[Tuple, ...]:
    """
    Detect hydrogen-bond interactions between two molecular groups using
    distance and angular geometric criteria.

    This function evaluates **intermolecular hydrogen bonds** between two
    sets of donor and acceptor atoms, considering both interaction directions:

        1. Group 1 acceptor ← Group 2 donor
        2. Group 2 acceptor ← Group 1 donor

    For each potential donor–acceptor pair, bonded hydrogens on the donor
    are tested against the acceptor atom using strict distance and angular
    constraints. Only geometries satisfying all criteria are reported.

    --------------------------------------------------------------------
    Geometric criteria
    --------------------------------------------------------------------

    A hydrogen bond is recorded if **all** of the following are satisfied:

    - Loose donor–acceptor distance ≤ ``cutoff_loose`` (Å)
    - Hydrogen–acceptor distance ≤ ``cutoff`` (Å)
    - Donor–H–Acceptor angle (DHA) ≥ ``min_angleDHA`` (degrees)
    - Hydrogen–Acceptor–X angle (HAX) ≥ ``min_angleHAX`` (degrees),
      where X is a heavy atom bonded to the acceptor

    Periodic boundary conditions are applied to all distance and angle
    calculations.

    --------------------------------------------------------------------
    Output record format
    --------------------------------------------------------------------

    Each detected hydrogen bond is returned as a **flat tuple** with the
    following schema:

        (
            "HYDROGEN-BOND",

            atom1_index, atom1_name, atom1_id, atom1_element, atom1_type,
            atom1_resname, atom1_resid, atom1_chain, atom1_segid,

            atom2_index, atom2_name, atom2_id, atom2_element, atom2_type,
            atom2_resname, atom2_resid, atom2_chain, atom2_segid,

            distance,           # H···A distance (Å)
            angle_DHA,          # Donor–H–Acceptor angle (degrees)
            angle_HAX,          # H–Acceptor–X angle (degrees)
            direction_label     # Interaction direction descriptor
        )

    ``atom_type`` is classified using backbone/sidechain logic
    (e.g. "BB", "SC", "UN") via ``_classify_type``.

    The ordering of atom1/atom2 reflects the interaction direction and is
    **schema-stable** across all interaction types to support downstream
    deduplication and frequency analysis.

    --------------------------------------------------------------------
    Memory and performance notes
    --------------------------------------------------------------------

    - No global state is modified
    - No per-frame caching is performed
    - All results are emitted as immutable tuples
    - Designed for safe use inside tight trajectory frame loops

    --------------------------------------------------------------------
    Parameters
    --------------------------------------------------------------------

    :param group1_acceptors: Acceptor atoms from group 1
    :param group2_acceptors: Acceptor atoms from group 2
    :param group1_donors: Donor atoms from group 1 (must have bonded hydrogens)
    :param group2_donors: Donor atoms from group 2 (must have bonded hydrogens)
    :param box: Periodic simulation box (3×3 matrix)
    :param min_angleDHA: Minimum Donor–H–Acceptor angle (degrees)
    :param min_angleHAX: Minimum H–Acceptor–X angle (degrees)
    :param cutoff_loose: Loose donor–acceptor prefilter distance (Å)
    :param cutoff: Strict hydrogen–acceptor distance cutoff (Å)

    :return: Tuple of hydrogen-bond interaction records
    :rtype: Tuple[Tuple, ...]
    """

    if (len(group1_acceptors) == 0 or len(group2_donors) == 0) and \
            (len(group2_acceptors) == 0 or len(group1_donors) == 0):
        return ()

    contacts: List[Tuple] = []

    # Group1 acceptors vs Group2 donors (G1: A, G2: D)
    if len(group1_acceptors) > 0 and len(group2_donors) > 0:
        dist = distance_array(group1_acceptors.positions, group2_donors.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            acceptor = group1_acceptors[i]
            donor = group2_donors[j]

            for h in donor.bonded_atoms:
                if not _is_hydrogen(h):
                    continue

                distance = calc_bonds(acceptor.position[None, :], h.position[None, :], box=box)[0]
                if distance > cutoff:
                    continue

                dha = np.degrees(calc_angles(donor.position, h.position, acceptor.position, box=box))
                if dha < min_angleDHA:
                    continue

                for x in acceptor.bonded_atoms:
                    hax = np.degrees(calc_angles(h.position, acceptor.position, x.position, box=box))
                    if hax < min_angleHAX:
                        continue

                    contacts.append((
                        "HYDROGEN-BOND",

                        int(acceptor.index), acceptor.name, int(acceptor.id), acceptor.element,
                        _classify_type(acceptor.resname, acceptor.name),
                        acceptor.resname, int(acceptor.resid), _get_chain(group1_acceptors, i), _get_segid(group1_acceptors, i),

                        int(donor.index), donor.name, int(donor.id), donor.element,
                        _classify_type(donor.resname, donor.name),
                        donor.resname, int(donor.resid), _get_chain(group2_donors, j), _get_segid(group2_donors, j),

                        float(distance), float(dha), float(hax),
                        "(G1)-A···H–D-(G2)"
                    ))

    # Group2 acceptors vs Group1 donors (G2: A, G1: D)
    if len(group2_acceptors) > 0 and len(group1_donors) > 0:
        dist = distance_array(group2_acceptors.positions, group1_donors.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            acceptor = group2_acceptors[i]
            donor = group1_donors[j]

            for h in donor.bonded_atoms:
                if not _is_hydrogen(h):
                    continue

                distance = calc_bonds(acceptor.position[None, :], h.position[None, :], box=box)[0]

                if distance > cutoff:
                    continue

                dha = np.degrees(calc_angles(donor.position, h.position, acceptor.position, box=box))
                if dha < min_angleDHA:
                    continue

                for x in acceptor.bonded_atoms:
                    hax = np.degrees(calc_angles(h.position, acceptor.position, x.position, box=box))
                    if hax < min_angleHAX:
                        continue

                    contacts.append((
                        "HYDROGEN-BOND",

                        int(donor.index), donor.name, int(donor.id), donor.element,
                        _classify_type(donor.resname, donor.name),
                        donor.resname, int(donor.resid), _get_chain(group1_donors, j), _get_segid(group1_donors, j),

                        int(acceptor.index), acceptor.name, int(acceptor.id), acceptor.element,
                        _classify_type(acceptor.resname, acceptor.name),
                        acceptor.resname, int(acceptor.resid), _get_chain(group2_acceptors, i), _get_segid(group2_acceptors, i),

                        float(distance), float(dha), float(hax),
                        "(G1)-D–H···A-(G2)"
                    ))

    return tuple(contacts)


def calculate_first_degree_water_bridge_contacts(u: Mda.Universe,
                                                 group1_acceptors: AtomGroup, group2_acceptors: AtomGroup,
                                                 group1_donors: AtomGroup, group2_donors: AtomGroup,
                                                 box: np.ndarray | None,
                                                 water_selection: str = WATER_RESIDUES_MDA_SELECTION_STR,
                                                 min_angleDHA: float = 110.0,
                                                 min_angleHAX: float = 80.0,
                                                 cutoff_loose: float = 6.0,
                                                 cutoff: float = 2.70,
                                                 second_degree: bool = False) -> Tuple[Tuple, ...]:
    """
    Calculates first-degree water bridge contacts between two groups of donors and acceptors using specified geometric criteria.
    The method identifies water molecules that act as mediators in these interactions by forming hydrogen bonds, and evaluates their
    relationships with the specified groups, accounting for various cutoff distances and angles.

    1. **Acceptor···H–O–H···Acceptor**
           - Water donates both H atoms to acceptors in Group1 and Group2.
           - Valid if:
             - Both H···A distances ≤ `cutoff`
             - DHA angles ≥ `min_angleDHA`
             - HAX angles ≥ `min_angleHAX` for both acceptors.

    2. **Acceptor···H–O···H–Donor**
       - Water acts as both donor and acceptor.
       - H of water donates to Group1 acceptor.
       - H from Group2 donor donates to water oxygen.
       - Valid if both H···A distances and angles satisfy thresholds.

    3. **Donor–H···O–H···Acceptor**
       - Water accepts H from Group1 donor, and donates H to Group2 acceptor.
        Valid if:
            - Donor–H···O and O–H···A distances ≤ `cutoff`
            - All DHA and HAX angles satisfy the limits.

    4. **Donor–H···O···H–Donor**
        - Water acts as hydrogen bond acceptor from both donors.
        - Both donor-H···O distances ≤ `cutoff`, with valid DHA/HAX angles.

    :param u: MDAnalysis universe containing trajectory, topology, and atomic information.
    :param group1_acceptors: Atom group containing acceptors from the first group.
    :param group2_acceptors: Atom group containing acceptors from the second group.
    :param group1_donors: Atom group containing donors from the first group.
    :param group2_donors: Atom group containing donors from the second group.
    :param box: Simulation box dimensions as a NumPy array.
    :param water_selection: Selection string for identifying water molecules in the universe. Default is WATER_RESIDUES_MDA_SELECTION_STR.
    :param min_angleDHA: Minimum angle (in degrees) for donor-hydrogen-acceptor interactions for validity. Default is 110.0.
    :param min_angleHAX: Minimum angle (in degrees) for hydrogen-acceptor-extra interaction validity. Default is 80.0.
    :param cutoff_loose: Loose cutoff distance (in angstroms) for preliminary water-bridge consideration. Default is 6.0.
    :param cutoff: Strict cutoff distance (in angstroms) for valid hydrogen bonding. Default is 2.70.
    :return: A tuple of tuples, where each inner tuple describes a water-bridge interaction including atomic indices, names, and geometric properties.
    """

    water_atoms = u.select_atoms(water_selection)
    water_oxygens = water_atoms.select_atoms("element O")

    if len(water_oxygens) == 0:
        return ()

    dist_mat_g1a = distance_array(group1_acceptors.positions, water_oxygens.positions, box=box)
    dist_mat_g2a = distance_array(group2_acceptors.positions, water_oxygens.positions, box=box)
    dist_mat_g2d = distance_array(group2_donors.positions, water_oxygens.positions, box=box)
    dist_mat_g1d = distance_array(group1_donors.positions, water_oxygens.positions, box=box)

    g1acc_to_w = defaultdict(list)
    g2acc_to_w = defaultdict(list)
    g2don_to_w = defaultdict(list)
    g1don_to_w = defaultdict(list)

    for i, j in zip(*np.where(dist_mat_g1a <= cutoff_loose)):
        g1acc_to_w[j].append(i)
    for n, k in zip(*np.where(dist_mat_g2a <= cutoff_loose)):
        g2acc_to_w[k].append(n)
    for d2i, wi in zip(*np.where(dist_mat_g2d <= cutoff_loose)):
        g2don_to_w[wi].append(d2i)
    for o, p in zip(*np.where(dist_mat_g1d <= cutoff_loose)):
        g1don_to_w[p].append(o)

    contacts = []
    for w_idx in set(g1acc_to_w) | set(g2acc_to_w) | set(g2don_to_w) | set(g1don_to_w):
        wo = water_oxygens[w_idx]
        h_atoms = [h for h in wo.bonded_atoms if _is_hydrogen(h)]
        if len(h_atoms) == 0:
            continue

        # Case 1-2: Acceptor···H–O–H···Acceptor |  Acceptor···H–O···H–Donor
        for wh1 in h_atoms:
            for i in g1acc_to_w.get(w_idx, []):
                a1 = group1_acceptors[i]
                dist1 = calc_bonds(a1.position[None, :], wh1.position[None, :], box=box)[0]
                if dist1 > cutoff:
                    continue

                a1_dha = np.degrees(calc_angles(wo.position, wh1.position, a1.position, box=box))
                if a1_dha < min_angleDHA:
                    continue

                # Acceptor···H–O–H···Acceptor
                for wh2 in h_atoms:
                    if wh1.index == wh2.index:
                        continue
                    for n in g2acc_to_w.get(w_idx, []):
                        a2 = group2_acceptors[n]
                        dist2 = calc_bonds(a2.position[None, :], wh2.position[None, :], box=box)[0]
                        if dist2 > cutoff:
                            continue

                        a2_dha = np.degrees(calc_angles(wo.position, wh2.position, a2.position, box=box))
                        if a2_dha < min_angleDHA:
                            continue

                        recorded = False
                        for a1b in a1.bonded_atoms:
                            a1_hax = np.degrees(calc_angles(wh1.position, a1.position, a1b.position, box=box))
                            if a1_hax < min_angleHAX:
                                continue

                            for a2b in a2.bonded_atoms:
                                a2_hax = np.degrees(calc_angles(wh2.position, a2.position, a2b.position, box=box))
                                if a2_hax < min_angleHAX:
                                    continue

                                contacts.append(
                                    ("WATER-BRIDGE-1",

                                     int(a1.index), a1.name, int(a1.id), a1.element,
                                     _classify_type(a1.resname, a1.name),
                                     a1.resname, int(a1.resid),
                                     _get_chain(group1_acceptors, i), _get_segid(group1_acceptors, i),

                                     int(a2.index), a2.name, int(a2.id), a2.element,
                                     _classify_type(a2.resname, a2.name),
                                     a2.resname, int(a2.resid),
                                     _get_chain(group2_acceptors, n), _get_segid(group2_acceptors, n),

                                     int(wo.index), wo.name, int(wo.id), wo.element,
                                     _classify_type(wo.resname, wo.name),
                                     wo.resname, int(wo.resid),
                                     _get_chain(water_oxygens, w_idx), _get_segid(water_oxygens, w_idx),

                                     float(dist1), float(dist2),

                                     float(a1_dha), float(a2_dha),
                                     float(a1_hax), float(a2_hax),

                                     "(G1)-A···H–O–H···A-(G2)"
                                     )
                                )
                                recorded = True
                                break
                            if recorded:
                                break

                # Acceptor···H–O···H–Donor
                for d2_idx in g2don_to_w.get(w_idx, []):
                    d2 = group2_donors[d2_idx]
                    for d2h in d2.bonded_atoms:
                        if not _is_hydrogen(d2h):
                            continue
                        dist2 = calc_bonds(wo.position[None, :], d2h.position[None, :], box=box)[0]
                        if dist2 > cutoff:
                            continue

                        a2_dha = np.degrees(calc_angles(d2.position, d2h.position, wo.position, box=box))
                        if a2_dha < min_angleDHA:
                            continue

                        for wh2 in h_atoms:
                            if wh1.index == wh2.index:
                                continue

                            a2_hax = np.degrees(calc_angles(d2.position, wo.position, wh2.position, box=box))
                            if a2_hax < min_angleHAX:
                                continue

                            recorded = False
                            for a1b in a1.bonded_atoms:
                                a1_hax = np.degrees(calc_angles(wh1.position, a1.position, a1b.position, box=box))
                                if a1_hax < min_angleHAX:
                                    continue

                                contacts.append(
                                    ("WATER-BRIDGE-1",

                                     int(a1.index), a1.name, int(a1.id), a1.element,
                                     _classify_type(a1.resname, a1.name),
                                     a1.resname, int(a1.resid),
                                     _get_chain(group1_acceptors, i), _get_segid(group1_acceptors, i),

                                     int(d2.index), d2.name, int(d2.id), d2.element,
                                     _classify_type(d2.resname, d2.name),
                                     d2.resname, int(d2.resid),
                                     _get_chain(group2_donors, d2_idx), _get_segid(group2_donors, d2_idx),

                                     int(wo.index), wo.name, int(wo.id), wo.element,
                                     _classify_type(wo.resname, wo.name),
                                     wo.resname, int(wo.resid),
                                     _get_chain(water_oxygens, w_idx), _get_segid(water_oxygens, w_idx),

                                     float(dist1), float(dist2),

                                     float(a1_dha), float(a2_dha),
                                     float(a1_hax), float(a2_hax),

                                     "(G1)-A···H–O···H-D-(G2)"
                                     )
                                )
                                recorded = True
                                break
                            if recorded:
                                break

        # Case 3-4: Donor–H···O–H···Acceptor |  Donor–H···O···H–Donor
        for a1_index in g1don_to_w.get(w_idx, []):
            a1 = group1_donors[a1_index]
            for a1h in a1.bonded_atoms:
                if not _is_hydrogen(a1h):
                    continue
                dist1 = calc_bonds(a1h.position[None, :], wo.position[None, :], box=box)[0]
                if dist1 > cutoff:
                    continue

                a1_dha = np.degrees(calc_angles(a1.position, a1h.position, wo.position, box=box))
                if a1_dha < min_angleDHA:
                    continue

                for wh1 in h_atoms:

                    for n in g2acc_to_w.get(w_idx, []):
                        a2 = group2_acceptors[n]
                        dist2 = calc_bonds(wh1.position[None, :], a2.position[None, :], box=box)[0]
                        if dist2 > cutoff:
                            continue

                        a2_dha = np.degrees(calc_angles(wo.position, wh1.position, a2.position, box=box))
                        if a2_dha < min_angleDHA:
                            continue

                        a1_hax = np.degrees(calc_angles(a1h.position, wo.position, wh1.position, box=box))
                        if a1_hax < min_angleHAX:
                            continue

                        recorded = False
                        for a2b in a2.bonded_atoms:
                            a2_hax = np.degrees(calc_angles(a2b.position, a2.position, wh1.position, box=box))
                            if a2_hax < min_angleHAX:
                                continue

                            contacts.append(
                                ("WATER-BRIDGE-1",

                                 int(a1.index), a1.name, int(a1.id), a1.element,
                                 _classify_type(a1.resname, a1.name),
                                 a1.resname, int(a1.resid),
                                 _get_chain(group1_donors, a1_index), _get_segid(group1_donors, a1_index),

                                 int(a2.index), a2.name, int(a2.id), a2.element,
                                 _classify_type(a2.resname, a2.name),
                                 a2.resname, int(a2.resid),
                                 _get_chain(group2_acceptors, n), _get_segid(group2_acceptors, n),

                                 int(wo.index), wo.name, int(wo.id), wo.element,
                                 _classify_type(wo.resname, wo.name),
                                 wo.resname, int(wo.resid),
                                 _get_chain(water_oxygens, w_idx), _get_segid(water_oxygens, w_idx),

                                 float(dist1), float(dist2),

                                 float(a1_dha), float(a2_dha),
                                 float(a1_hax), float(a2_hax),

                                 "(G1)-D-H···O–H···A-(G2)"
                                 )
                            )
                            recorded = True
                            break
                        if recorded:
                            break

                # Donor–H···O···H–Donor
                for d2_idx in g2don_to_w.get(w_idx, []):
                    a2 = group2_donors[d2_idx]
                    for a2h in a2.bonded_atoms:
                        if not _is_hydrogen(a2h):
                            continue
                        dist2 = calc_bonds(a2h.position[None, :], wo.position[None, :], box=box)[0]
                        if dist2 > cutoff:
                            continue
                        a2_dha = np.degrees(calc_angles(a2.position, a2h.position, wo.position, box=box))
                        if a2_dha < min_angleDHA:
                            continue

                        recorded = False
                        for wh1 in h_atoms:
                            a1_hax = np.degrees(calc_angles(a1h.position, wo.position, wh1.position, box=box))
                            if a1_hax < min_angleHAX:
                                continue
                            for wh2 in h_atoms:
                                if wh1.index == wh2.index:
                                    continue

                                a2_hax = np.degrees(calc_angles(a2h.position, wo.position, wh2.position, box=box))
                                if a2_hax < min_angleHAX:
                                    continue

                                contacts.append(
                                    ("WATER-BRIDGE-1",

                                     int(a1.index), a1.name, int(a1.id), a1.element,
                                     _classify_type(a1.resname, a1.name),
                                     a1.resname, int(a1.resid),
                                     _get_chain(group1_donors, a1_index), _get_segid(group1_donors, a1_index),

                                     int(a2.index), a2.name, int(a2.id), a2.element,
                                     _classify_type(a2.resname, a2.name),
                                     a2.resname, int(a2.resid),
                                     _get_chain(group2_donors, d2_idx), _get_segid(group2_donors, d2_idx),

                                     int(wo.index), wo.name, int(wo.id), wo.element,
                                     _classify_type(wo.resname, wo.name),
                                     wo.resname, int(wo.resid),
                                     _get_chain(water_oxygens, w_idx), _get_segid(water_oxygens, w_idx),

                                     float(dist1), float(dist2),

                                     float(a1_dha), float(a2_dha),
                                     float(a1_hax), float(a2_hax),

                                     "(G1)-D-H···O···H-D-(G2)"
                                     )
                                )
                                recorded = True
                                break
                            if recorded:
                                break

    # ------------------------------------------------------------------
    # Second-degree water bridges:  Group1 — W1 — W2 — Group2
    #
    # Each link in the chain is a proper H-bond.  A water H atom that
    # donates in one link cannot donate in another, so we track which
    # specific H index is consumed at each stage and forbid reuse on
    # the same water molecule.
    #
    # Tuple stored per link:
    #   g-side  : (atom, local_idx, group_ref, dist, dha, hax, role,
    #              w_h_used)
    #             w_h_used = atom.index of the water-H consumed on
    #             that water, or None when the water acts as acceptor.
    #   ww-side : (dist, dha, w1_h_used, w2_h_used)
    # ------------------------------------------------------------------
    if second_degree and len(water_oxygens) >= 2:
        ww_dist = distance_array(water_oxygens.positions, water_oxygens.positions, box=box)

        ww_pairs: List[Tuple[int, int]] = []
        for wi, wj in zip(*np.where(ww_dist <= cutoff_loose)):
            if wi >= wj:
                continue
            ww_pairs.append((wi, wj))

        for w1_idx, w2_idx in ww_pairs:
            wo1 = water_oxygens[w1_idx]
            wo2 = water_oxygens[w2_idx]
            h1_atoms = [h for h in wo1.bonded_atoms if _is_hydrogen(h)]
            h2_atoms = [h for h in wo2.bonded_atoms if _is_hydrogen(h)]
            if not h1_atoms or not h2_atoms:
                continue

            # --- Collect ALL valid water–water H-bonds with H-atom tracking ---
            # Each entry: (dist, dha, w1_h_used_index|None, w2_h_used_index|None)
            ww_links: List[Tuple[float, float, int | None, int | None]] = []

            # W1 donates H to W2 oxygen  →  w1_h consumed, w2_h free
            for wh in h1_atoms:
                d = calc_bonds(wo2.position[None, :], wh.position[None, :], box=box)[0]
                if d > cutoff:
                    continue
                ang = np.degrees(calc_angles(wo1.position, wh.position, wo2.position, box=box))
                if ang >= min_angleDHA:
                    ww_links.append((float(d), float(ang), wh.index, None))

            # W2 donates H to W1 oxygen  →  w2_h consumed, w1_h free
            for wh in h2_atoms:
                d = calc_bonds(wo1.position[None, :], wh.position[None, :], box=box)[0]
                if d > cutoff:
                    continue
                ang = np.degrees(calc_angles(wo2.position, wh.position, wo1.position, box=box))
                if ang >= min_angleDHA:
                    ww_links.append((float(d), float(ang), None, wh.index))

            if not ww_links:
                continue

            # Try both assignments:  (W1↔G1, W2↔G2)  and  (W2↔G1, W1↔G2)
            for near_g1_idx, far_g2_idx, wo_near, wo_far, h_near, h_far in [
                (w1_idx, w2_idx, wo1, wo2, h1_atoms, h2_atoms),
                (w2_idx, w1_idx, wo2, wo1, h2_atoms, h1_atoms),
            ]:
                # --- G1 side links (atom, li, grp, dist, dha, hax, role, w_near_h_used) ---
                g1_links = []

                # G1 acceptor ← wo_near donates H  →  that H is consumed on wo_near
                for wh in h_near:
                    for i in g1acc_to_w.get(near_g1_idx, []):
                        a1 = group1_acceptors[i]
                        d1 = calc_bonds(a1.position[None, :], wh.position[None, :], box=box)[0]
                        if d1 > cutoff:
                            continue
                        dha1 = np.degrees(calc_angles(wo_near.position, wh.position, a1.position, box=box))
                        if dha1 < min_angleDHA:
                            continue
                        for a1b in a1.bonded_atoms:
                            hax1 = np.degrees(calc_angles(wh.position, a1.position, a1b.position, box=box))
                            if hax1 >= min_angleHAX:
                                g1_links.append((a1, i, group1_acceptors,
                                                 float(d1), float(dha1), float(hax1),
                                                 "A", wh.index))
                                break

                # G1 donor → wo_near accepts H  →  no wo_near H consumed
                for g1d_i in g1don_to_w.get(near_g1_idx, []):
                    d1_atom = group1_donors[g1d_i]
                    for d1h in d1_atom.bonded_atoms:
                        if not _is_hydrogen(d1h):
                            continue
                        d1 = calc_bonds(d1h.position[None, :], wo_near.position[None, :], box=box)[0]
                        if d1 > cutoff:
                            continue
                        dha1 = np.degrees(calc_angles(d1_atom.position, d1h.position, wo_near.position, box=box))
                        if dha1 < min_angleDHA:
                            continue
                        for wh in h_near:
                            hax1 = np.degrees(calc_angles(d1h.position, wo_near.position, wh.position, box=box))
                            if hax1 >= min_angleHAX:
                                g1_links.append((d1_atom, g1d_i, group1_donors,
                                                 float(d1), float(dha1), float(hax1),
                                                 "D", None))
                                break
                        break

                if not g1_links:
                    continue

                # --- G2 side links (atom, li, grp, dist, dha, hax, role, w_far_h_used) ---
                g2_links = []

                # G2 acceptor ← wo_far donates H  →  that H is consumed on wo_far
                for wh in h_far:
                    for j in g2acc_to_w.get(far_g2_idx, []):
                        a2 = group2_acceptors[j]
                        d2 = calc_bonds(a2.position[None, :], wh.position[None, :], box=box)[0]
                        if d2 > cutoff:
                            continue
                        dha2 = np.degrees(calc_angles(wo_far.position, wh.position, a2.position, box=box))
                        if dha2 < min_angleDHA:
                            continue
                        for a2b in a2.bonded_atoms:
                            hax2 = np.degrees(calc_angles(wh.position, a2.position, a2b.position, box=box))
                            if hax2 >= min_angleHAX:
                                g2_links.append((a2, j, group2_acceptors,
                                                 float(d2), float(dha2), float(hax2),
                                                 "A", wh.index))
                                break

                # G2 donor → wo_far accepts H  →  no wo_far H consumed
                for g2d_j in g2don_to_w.get(far_g2_idx, []):
                    d2_atom = group2_donors[g2d_j]
                    for d2h in d2_atom.bonded_atoms:
                        if not _is_hydrogen(d2h):
                            continue
                        d2 = calc_bonds(d2h.position[None, :], wo_far.position[None, :], box=box)[0]
                        if d2 > cutoff:
                            continue
                        dha2 = np.degrees(calc_angles(d2_atom.position, d2h.position, wo_far.position, box=box))
                        if dha2 < min_angleDHA:
                            continue
                        for wh in h_far:
                            hax2 = np.degrees(calc_angles(d2h.position, wo_far.position, wh.position, box=box))
                            if hax2 >= min_angleHAX:
                                g2_links.append((d2_atom, g2d_j, group2_donors,
                                                 float(d2), float(dha2), float(hax2),
                                                 "D", None))
                                break
                        break

                if not g2_links:
                    continue

                # --- Combine: check H-atom reuse constraints ---
                # wo_near is w1 or w2 depending on the orientation loop.
                # Map back to the canonical (wo1, wo2) frame so we can
                # compare against ww_link H-usage.
                near_is_w1 = (wo_near is wo1)

                for (a1, a1_li, a1_grp, dist1, dha1, hax1, role1, near_h) in g1_links:
                    for (a2, a2_li, a2_grp, dist2, dha2, hax2, role2, far_h) in g2_links:
                        for (ww_d, ww_a, ww_w1h, ww_w2h) in ww_links:
                            # Determine which canonical water-H the ww link
                            # and the g-side links consume.
                            if near_is_w1:
                                w1_h_g = near_h   # H used by G1-W1 link
                                w2_h_g = far_h    # H used by W2-G2 link
                            else:
                                w1_h_g = far_h    # W1 is on the G2 side
                                w2_h_g = near_h   # W2 is on the G1 side

                            # Reject if the same H on W1 is donated twice
                            if w1_h_g is not None and ww_w1h is not None and w1_h_g == ww_w1h:
                                continue
                            # Reject if the same H on W2 is donated twice
                            if w2_h_g is not None and ww_w2h is not None and w2_h_g == ww_w2h:
                                continue

                            g1_tag = "D" if role1 == "D" else "A"
                            g2_tag = "D" if role2 == "D" else "A"
                            bridge_label = f"(G1)-{g1_tag}···W1···W2···{g2_tag}-(G2)"

                            contacts.append(
                                ("WATER-BRIDGE-2",

                                 int(a1.index), a1.name, int(a1.id), a1.element,
                                 _classify_type(a1.resname, a1.name),
                                 a1.resname, int(a1.resid),
                                 _get_chain(a1_grp, a1_li), _get_segid(a1_grp, a1_li),

                                 int(a2.index), a2.name, int(a2.id), a2.element,
                                 _classify_type(a2.resname, a2.name),
                                 a2.resname, int(a2.resid),
                                 _get_chain(a2_grp, a2_li), _get_segid(a2_grp, a2_li),

                                 int(wo_near.index), wo_near.name, int(wo_near.id), wo_near.element,
                                 _classify_type(wo_near.resname, wo_near.name),
                                 wo_near.resname, int(wo_near.resid),
                                 _get_chain(water_oxygens, near_g1_idx), _get_segid(water_oxygens, near_g1_idx),

                                 int(wo_far.index), wo_far.name, int(wo_far.id), wo_far.element,
                                 _classify_type(wo_far.resname, wo_far.name),
                                 wo_far.resname, int(wo_far.resid),
                                 _get_chain(water_oxygens, far_g2_idx), _get_segid(water_oxygens, far_g2_idx),

                                 float(dist1), float(ww_d), float(dist2),

                                 float(dha1), float(ww_a), float(dha2),
                                 float(hax1), float(hax2),

                                 bridge_label
                                 )
                            )

    return tuple(contacts)


def calculate_ionic_contacts(group1_positive: AtomGroup, group1_negative: AtomGroup,
                             group2_positive: AtomGroup, group2_negative: AtomGroup,
                             box: np.ndarray | None,
                             cutoff: float = 5.0) -> Tuple[Tuple, ...]:
    """
    Calculate ionic contacts between two groups of atoms.

    This function identifies ionic contacts between two groups of positively and
    negatively charged atoms. It calculates distances between provided groups using
    a given cutoff and periodic boundary conditions defined by the box.

    :param group1_positive: A group of positively charged atoms in group 1.
    :type group1_positive: AtomGroup
    :param group1_negative: A group of negatively charged atoms in group 1.
    :type group1_negative: AtomGroup
    :param group2_positive: A group of positively charged atoms in group 2.
    :type group2_positive: AtomGroup
    :param group2_negative: A group of negatively charged atoms in group 2.
    :type group2_negative: AtomGroup
    :param box: The periodic boundary box as a 3x3 numpy array.
    :type box: numpy.ndarray
    :param cutoff: Distance cutoff for defining ionic contacts. Default is 5.0.
    :type cutoff: float
    :return: A tuple containing ionic contact information. Each contact is a tuple
        consisting of ionic interaction type, atom indices, atom names, atom IDs,
        elements, residue names, residue IDs, chain identifiers, segment identifiers,
        distances and contact label.
    :rtype: Tuple[Tuple, ...]
    """

    contacts: List[Tuple] = []

    def _scan(gA: AtomGroup, gB: AtomGroup, label: str) -> None:
        """
        Scans for close contacts between two groups of atoms and records the interactions.

        This function calculates the pairwise distances between atoms in two groups and,
        if the distance satisfies a cutoff criterion, it identifies the pair as having
        close contact. The details of the interacting pairs are recorded with their
        respective attributes and associated metadata. Only pairs within the given cutoff
        distance are processed.

        :param gA: The first group of atoms to be analyzed.
        :type gA: AtomGroup
        :param gB: The second group of atoms to be analyzed.
        :type gB: AtomGroup
        :param label: A label describing the context or type of interaction being analyzed.
        :type label: str
        :return: This function does not return anything. The results are appended to an
                 external global list or data structure.
        :rtype: None
        """
        if len(gA) == 0 or len(gB) == 0:
            return

        dist_mat = distance_array(gA.positions, gB.positions, box=box)
        ai, bj = np.where(dist_mat <= cutoff)

        for i, j in zip(ai, bj):
            a1 = gA[i]
            a2 = gB[j]

            contacts.append((
                "IONIC",

                int(a1.index), a1.name, int(a1.id), a1.element,
                _classify_type(a1.resname, a1.name),
                a1.resname, int(a1.resid),
                _get_chain(gA, i), _get_segid(gA, i),

                int(a2.index), a2.name, int(a2.id), a2.element,
                _classify_type(a2.resname, a2.name),
                a2.resname, int(a2.resid),
                _get_chain(gB, j), _get_segid(gB, j),

                float(dist_mat[i, j]),
                label
            ))

    # Positive (G1) — Negative (G2)
    _scan(group1_positive, group2_negative, "(G1+)···(G2-)")
    # Negative (G1) — Positive (G2)
    _scan(group1_negative, group2_positive, "(G1-)···(G2+)")

    return tuple(contacts)


def calculate_halogen_bond_contacts(group1_halogens: AtomGroup, group2_halogens: AtomGroup,
                                    group1_acceptors: AtomGroup, group2_acceptors: AtomGroup,
                                    group1_donors: AtomGroup, group2_donors: AtomGroup,
                                    box: np.ndarray | None,
                                    donor_min_angle_as_donor: float = 140.0,
                                    acceptor_min_angle_as_donor: float = 90.0,
                                    donor_min_angle_as_acceptor: float = 120.0,
                                    acceptor_min_angle_as_acceptor: float = 90.0,
                                    acceptor_max_angle: float = 170.0,
                                    cutoff_loose: float = 6.0,
                                    cutoff: float = 3.5) -> Tuple[Tuple, ...]:
    """
    Calculate and identify halogen-bond contacts between specified atom groups.

    This function evaluates potential halogen-bond contacts between two sets
    of atomic groups provided as input parameters. The contacts are determined
    by considering geometric and angular constraints that define the conditions
    for halogen bonds. The function operates bidirectionally to ensure both
    sets of atom groups are analyzed as donors and acceptors.

    1. Group 1 halogens as donors vs. Group 2 acceptors
        2. Group 2 halogens as donors vs. Group 1 acceptors
        3. Group 1 halogens as acceptors vs. Group 2 donors
        4. Group 2 halogens as acceptors vs. Group 1 donors

    The geometric definitions are:

    **Halogen acting as donor (X···A)**:
        - Loose halogen–acceptor prefilter distance: `cutoff_loose` (default 6.0 Å)
        - Strict halogen–acceptor distance cutoff: `cutoff` (default 3.5 Å)
        - Minimum X–C···A angle: `donor_min_angle_as_donor` (default 140°)
          (X = halogen donor, C = covalently bound atom to X)
        - Minimum C···A–X' angle: `acceptor_min_angle_as_donor` (default 90°)
          (A = acceptor atom, X' = heavy atom bonded to acceptor)

    **Halogen acting as acceptor (D–H···X)**:
        - Loose donor–halogen prefilter distance: `cutoff_loose` (default 6.0 Å)
        - Strict donor–halogen distance cutoff: `cutoff` (default 3.5 Å)
        - Minimum D–H···X angle: `donor_min_angle_as_acceptor` (default 120°)
          (D = heavy donor atom, H = hydrogen)
        - Minimum H···X–C angle: `acceptor_min_angle_as_acceptor` (default 90°)
          (C = atom covalently bonded to X)
        - Maximum D–X–C angle: `acceptor_max_angle` (default 170°)
          (ensures correct acceptor geometry for halogen)

    :param group1_halogens: Atom group consisting of potential halogen donors
        from group 1.
    :type group1_halogens: AtomGroup

    :param group2_halogens: Atom group consisting of potential halogen donors
        from group 2.
    :type group2_halogens: AtomGroup

    :param group1_acceptors: Atom group representing potential acceptors
        from group 1.
    :type group1_acceptors: AtomGroup

    :param group2_acceptors: Atom group representing potential acceptors
        from group 2.
    :type group2_acceptors: AtomGroup

    :param group1_donors: Atom group defining potential donors from group 1.
    :type group1_donors: AtomGroup

    :param group2_donors: Atom group defining potential donors from group 2.
    :type group2_donors: AtomGroup

    :param box: The periodic boundary conditions of the simulation box,
        expressed as a 3D numpy array.
    :type box: np.ndarray

    :param donor_min_angle_as_donor: Minimum angle (in degrees) required
        for a donor atom to be valid in a donor role.
    :type donor_min_angle_as_donor: float

    :param acceptor_min_angle_as_donor: Minimum angle (in degrees) required
        for an acceptor atom when evaluated in a donor role.
    :type acceptor_min_angle_as_donor: float

    :param donor_min_angle_as_acceptor: Minimum angle (in degrees) required
        for a donor atom in an acceptor role.
    :type donor_min_angle_as_acceptor: float

    :param acceptor_min_angle_as_acceptor: Minimum angle (in degrees) for an
        acceptor atom when acting in its default role as an acceptor.
    :type acceptor_min_angle_as_acceptor: float

    :param acceptor_max_angle: Maximum angle (in degrees) between covalent
        partners of acceptor and donor alignment to validate the acceptor geometry.
    :type acceptor_max_angle: float

    :param cutoff_loose: Distance (in angstroms) defining a loose cutoff for
        identifying nearby atom pairs for further evaluation.
    :type cutoff_loose: float

    :param cutoff: Stringent distance criterion (in angstroms) for determining
        valid halogen-bond contacts.
    :type cutoff: float

    :return: A tuple of all identified halogen-bond contacts that satisfy the
        geometric and angular constraints. Each entry in the tuple represents
        a specific contact.
    :rtype: Tuple[Tuple, ...]
    """

    contacts = []

    # ---------- Direction 1 ----------
    # Case 1a: G1 acceptors vs G2 halogen donors  (X donor: A···X–C)
    if len(group1_acceptors) > 0 and len(group2_halogens) > 0:
        dist = distance_array(group1_acceptors.positions, group2_halogens.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            A = group1_acceptors[i]  # acceptor
            X = group2_halogens[j]  # halogen donor
            C_list = _heavy_neighbors(X)  # covalent partner(s) of X
            if not C_list:
                continue

            # strict distance: X···A
            XA = calc_bonds(A.position[None, :], X.position[None, :], box=box)[0]
            if XA > cutoff:
                continue

            # angles: ∠C–X···A and ∠X···A–Y (Y heavy neighbor of A)
            Y_list = _heavy_neighbors(A)
            if not Y_list:
                continue

            valid = False
            for C in C_list:
                ang_CXA = float(np.degrees(calc_angles(C.position, X.position, A.position, box=box)))
                if ang_CXA < donor_min_angle_as_donor:
                    continue
                for Y in Y_list:
                    ang_XAY = float(np.degrees(calc_angles(X.position, A.position, Y.position, box=box)))
                    if ang_XAY < acceptor_min_angle_as_donor:
                        continue
                    # record first valid Y (avoid duplicates)
                    contacts.append((
                        "HALOGEN-BOND",

                        int(A.index), A.name, int(A.id), A.element,
                        _classify_type(A.resname, A.name),
                        A.resname, int(A.resid), _get_chain(group1_acceptors, i), _get_segid(group1_acceptors, i),

                        int(X.index), X.name, int(X.id), X.element,
                        _classify_type(X.resname, X.name),
                        X.resname, int(X.resid), _get_chain(group2_halogens, j), _get_segid(group2_halogens, j),

                        float(XA), ang_CXA, ang_XAY,
                        "(G1)-A···X–Donor-(G2)"
                    ))
                    valid = True
                    break
                if valid:
                    break

    # Case 1b: G1 halogen acceptors vs G2 donors (D–H···X)
    if len(group1_halogens) > 0 and len(group2_donors) > 0:
        dist = distance_array(group1_halogens.positions, group2_donors.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            X = group1_halogens[i]  # halogen acceptor
            D = group2_donors[j]  # donor (has H)
            C_list = _heavy_neighbors(X)  # covalent partner(s) of X
            if not C_list:
                continue

            for H in D.bonded_atoms:
                if not _is_hydrogen(H):
                    continue

                HX = calc_bonds(H.position[None, :], X.position[None, :], box=box)[0]
                if HX > cutoff:
                    continue

                ang_DHX = float(np.degrees(calc_angles(D.position, H.position, X.position, box=box)))
                if ang_DHX < donor_min_angle_as_acceptor:
                    continue

                valid = False
                for C in C_list:
                    ang_HXC = float(np.degrees(calc_angles(H.position, X.position, C.position, box=box)))
                    if ang_HXC < acceptor_min_angle_as_acceptor:
                        continue
                    # optional acceptor geometry: ∠D–X–C ≤ acceptor_max_angle
                    ang_DXC = float(np.degrees(calc_angles(D.position, X.position, C.position, box=box)))
                    if ang_DXC > acceptor_max_angle:
                        continue

                    contacts.append((
                        "HALOGEN-BOND",

                        int(X.index), X.name, int(X.id), X.element,
                        _classify_type(X.resname, X.name),
                        X.resname, int(X.resid), _get_chain(group1_halogens, i), _get_segid(group1_halogens, i),

                        int(D.index), D.name, int(D.id), D.element,
                        _classify_type(D.resname, D.name),
                        D.resname, int(D.resid), _get_chain(group2_donors, j), _get_segid(group2_donors, j),

                        float(HX), ang_DHX, ang_HXC,
                        "(G1)-X(acceptor)···H–D-(G2)"
                    ))
                    valid = True
                    break
                if valid:
                    break

    # ---------- Direction 2 ----------
    # Case 2a: G2 acceptors vs G1 halogen donors (X donor)
    if len(group2_acceptors) > 0 and len(group1_halogens) > 0:
        dist = distance_array(group2_acceptors.positions, group1_halogens.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            A = group2_acceptors[i]
            X = group1_halogens[j]
            C_list = _heavy_neighbors(X)
            if not C_list:
                continue

            XA = calc_bonds(A.position[None, :], X.position[None, :], box=box)[0]
            if XA > cutoff:
                continue

            Y_list = _heavy_neighbors(A)
            if not Y_list:
                continue

            valid = False
            for C in C_list:
                ang_CXA = float(np.degrees(calc_angles(C.position, X.position, A.position, box=box)))
                if ang_CXA < donor_min_angle_as_donor:
                    continue
                for Y in Y_list:
                    ang_XAY = float(np.degrees(calc_angles(X.position, A.position, Y.position, box=box)))
                    if ang_XAY < acceptor_min_angle_as_donor:
                        continue

                    # Mirror your hydrogen-bond ordering (donor block first in direction 2)
                    contacts.append((
                        "HALOGEN-BOND",

                        int(X.index), X.name, int(X.id), X.element,
                        _classify_type(X.resname, X.name),
                        X.resname, int(X.resid), _get_chain(group1_halogens, j), _get_segid(group1_halogens, j),

                        int(A.index), A.name, int(A.id), A.element,
                        _classify_type(A.resname, A.name),
                        A.resname, int(A.resid), _get_chain(group2_acceptors, i), _get_segid(group2_acceptors, i),

                        float(XA), ang_CXA, ang_XAY,
                        "(G1)-Donor X–A···(G2)"
                    ))
                    valid = True
                    break
                if valid:
                    break

    # Case 2b: G2 halogen acceptors vs G1 donors (D–H···X)
    if len(group2_halogens) > 0 and len(group1_donors) > 0:
        dist = distance_array(group2_halogens.positions, group1_donors.positions, box=box)
        i_indices, j_indices = np.where(dist < cutoff_loose)

        for i, j in zip(i_indices, j_indices):
            X = group2_halogens[i]
            D = group1_donors[j]
            C_list = _heavy_neighbors(X)
            if not C_list:
                continue

            for H in D.bonded_atoms:
                if not _is_hydrogen(H):
                    continue

                HX = calc_bonds(H.position[None, :], X.position[None, :], box=box)[0]
                if HX > cutoff:
                    continue

                ang_DHX = float(np.degrees(calc_angles(D.position, H.position, X.position, box=box)))
                if ang_DHX < donor_min_angle_as_acceptor:
                    continue

                valid = False
                for C in C_list:
                    ang_HXC = float(np.degrees(calc_angles(H.position, X.position, C.position, box=box)))
                    if ang_HXC < acceptor_min_angle_as_acceptor:
                        continue
                    ang_DXC = float(np.degrees(calc_angles(D.position, X.position, C.position, box=box)))
                    if ang_DXC > acceptor_max_angle:
                        continue

                    contacts.append((
                        "HALOGEN-BOND",

                        int(D.index), D.name, int(D.id), D.element,
                        _classify_type(D.resname, D.name),
                        D.resname, int(D.resid), _get_chain(group1_donors, j), _get_segid(group1_donors, j),

                        int(X.index), X.name, int(X.id), X.element,
                        _classify_type(X.resname, X.name),
                        X.resname, int(X.resid), _get_chain(group2_halogens, i), _get_segid(group2_halogens, i),

                        float(HX), ang_DHX, ang_HXC,
                        "(G1)-D–H···X(acceptor)-(G2)"
                    ))
                    valid = True
                    break
                if valid:
                    break

    return tuple(contacts)


def calculate_pi_stacking_contacts(group1_aromatic_rings: List["AtomGroup"], group2_aromatic_rings: List["AtomGroup"],
                                   box: np.ndarray | None,
                                   dmax_f2f: float = 4.4,       # face-to-face distance cutoff
                                   angmax_f2f: float = 30.0,    # face-to-face max angle (deg)
                                   dmax_f2e: float = 5.5,       # face-to-edge distance cutoff
                                   angmin_f2e: float = 60.0,    # face-to-edge min angle (deg)
                                   ) -> Tuple[Tuple, ...]:
    """
    Calculates π-stacking contacts between two groups of aromatic rings based on geometric
    criteria such as centroid distances, face-to-face angles, and face-to-edge angles.

    The function identifies and categorizes π-stacking arrangements into "parallel" (face-to-face)
    or "T-shaped" (face-to-edge) based on the provided thresholds. It computes centroid
    distances, ring normals, angles between planes, and generates metadata for the involved
    atom groups.

    π–π detection (distance + plane–plane angle only), PBC-aware distances.
    Face-to-face:  d_cc ≤ 4.4 Å and θ ≤ 30°
    Face-to-edge:  d_cc ≤ 5.5 Å and θ ≥ 60°

    :param group1_aromatic_rings: The first group of aromatic rings represented as
        a list of AtomGroup objects. Each AtomGroup contains atomic positions and
        supporting geometric methods.
    :type group1_aromatic_rings: list[AtomGroup]
    :param group2_aromatic_rings: The second group of aromatic rings represented
        similarly as the first group.
    :type group2_aromatic_rings: list[AtomGroup]
    :param box: The box vector defining periodic boundary conditions for distance
        calculations.
    :type box: numpy.ndarray
    :param dmax_f2f: Face-to-face maximum distance cutoff in Angstroms beyond which
        interactions are ignored.
    :type dmax_f2f: float
    :param angmax_f2f: Face-to-face angle cutoff in degrees. Angles below this
        threshold qualify as "parallel" stacking.
    :type angmax_f2f: float
    :param dmax_f2e: Face-to-edge maximum distance cutoff in Angstroms for "T-shaped"
        stacking.
    :type dmax_f2e: float
    :param angmin_f2e: Face-to-edge minimum angle cutoff in degrees for qualifying
        interactions as "T-shaped" stacking.
    :type angmin_f2e: float
    :return: A tuple of tuples where each inner tuple describes a detected π-stacking
        interaction. Each inner tuple includes interaction type, metadata strings for
        each aromatic group, geometric parameters, and stacking classification.
    :rtype: tuple[tuple, ...]
    """

    if not group1_aromatic_rings or not group2_aromatic_rings:
        return ()

    # Precompute centroids and normals
    g1C = np.array([r.center_of_geometry() for r in group1_aromatic_rings], dtype=np.float64)
    g2C = np.array([r.center_of_geometry() for r in group2_aromatic_rings], dtype=np.float64)
    g1N = np.array([_ring_normal_cross(r) for r in group1_aromatic_rings], dtype=np.float64)
    g2N = np.array([_ring_normal_cross(r) for r in group2_aromatic_rings], dtype=np.float64)

    # PBC-aware centroid distances
    dmat = distance_array(g1C, g2C, box=box)  # shape (n1, n2)
    # Loose prefilter: anything within the larger cutoff (f2e)
    ii, jj = np.where(dmat <= max(dmax_f2f, dmax_f2e))
    if ii.size == 0:
        return ()

    contacts: list[tuple] = []

    for a, b in zip(ii, jj):
        i, j = int(a), int(b)
        n1 = g1N[i]
        n2 = g2N[j]
        # plane-plane angle θ in [0, 90]: acos(|n1·n2|)
        cos_th = float(np.clip(abs(np.dot(n1, n2)), 0.0, 1.0))
        theta = float(np.degrees(np.arccos(cos_th)))
        dcc = float(dmat[i, j])

        is_f2f = (dcc <= dmax_f2f) and (theta <= angmax_f2f)
        is_f2e = (dcc <= dmax_f2e) and (theta >= angmin_f2e)

        if not (is_f2f or is_f2e):
            continue

        stacking_type = "parallel" if is_f2f else "T-shaped"
        direction = "(G1)-ring···ring-(G2)"

        r1_meta = _ring_meta_str(group1_aromatic_rings[i])
        r2_meta = _ring_meta_str(group2_aromatic_rings[j])

        contacts.append((
            "PI-STACKING",
            *r1_meta,
            *r2_meta,
            dcc, theta, stacking_type, direction
        ))

    return tuple(contacts)


def calculate_pi_cation_contacts(group1_aromatic_rings: List["AtomGroup"], group2_aromatic_rings: List["AtomGroup"],
                                 group1_positive: "AtomGroup", group2_positive: "AtomGroup",
                                 box: np.ndarray | None,
                                 cutoff: float = 6.6,
                                 max_angle: float = 30.0) -> Tuple[Tuple, ...]:
    """
    Calculate π-cation interactions between two groups of aromatic rings and positively charged atoms,
    considering periodic boundary conditions. The interactions are identified based on the distances
    and angles between the centroids and normals of the aromatic rings and the positively charged atoms.

    This method is commonly used in structural bioinformatics or computational chemistry to study
    π-cation interactions, which are important in molecular recognition and stability.

    :param group1_aromatic_rings: List of aromatic rings for group 1, represented as AtomGroup instances.
    :param group2_aromatic_rings: List of aromatic rings for group 2, represented as AtomGroup instances.
    :param group1_positive: AtomGroup representing positively charged atoms in group 1.
    :param group2_positive: AtomGroup representing positively charged atoms in group 2.
    :param box: Periodic boundary box represented as a NumPy ndarray.
    :param cutoff: Cutoff distance (in Ångstroms) between centroids of aromatic rings and positively
        charged atoms. Defaults to 6.6 Å.
    :param max_angle: Maximum angle (in degrees) allowed between the normal of the aromatic ring plane
        and the centroid-to-cation vector, folded within 90 degrees. Defaults to 30 degrees.
    :return: Tuple of detected contacts. Each contact is represented as a tuple containing metadata
        about the interacting groups, detailed atom information, distance, angle, and interaction type.
    """

    # Quick exits: nothing to compare
    if (not group1_aromatic_rings and not group2_aromatic_rings) or \
            ((group1_positive is None or len(group1_positive) == 0) and
             (group2_positive is None or len(group2_positive) == 0)):
        return ()

    def _precompute_rings(rings):
        """
        Precomputes the center of geometry and normal vectors for a collection of rings.

        This function calculates the center of geometry and normal vectors for each
        ring provided in the input. If no rings are provided, the function returns
        empty arrays for both centers and normals. The results are used in further
        geometry or spatial computations.

        :param rings: A collection of ring objects, each having a `center_of_geometry`
                      method for computing the geometric center and a `_ring_normal_cross`
                      function to calculate the normal vector.
        :type rings: list
        :return: A tuple where the first element is a NumPy array containing the
                 computed centers of geometry (shape: (m, 3)) and the second element
                 is a NumPy array with the normal vectors (shape: (m, 3)) for each ring.
        :rtype: tuple
        """
        if not rings:
            return (np.empty((0, 3), dtype=np.float64),
                    np.empty((0, 3), dtype=np.float64))
        m = len(rings)
        C = np.empty((m, 3), dtype=np.float64)
        N = np.empty((m, 3), dtype=np.float64)
        for idx, ring in enumerate(rings):
            C[idx] = ring.center_of_geometry()
            N[idx] = _ring_normal_cross(ring)
        return C, N

    def cation_meta(atom: Atom) -> Tuple[int, str, int, str, str, str, int, str, str]:
        """
        Extracts metadata related to a cation atom. The function processes the attributes
        of the given atom object and determines the atom's classification and metadata
        based on its properties.

        :param atom: An object representing an atom in the molecular structure. The object
            is expected to have the attributes `index`, `name`, `id`, `element`, `resname`,
            and `resid`. Optional attributes may include `chainID` and `segid`, which will be
            set to empty strings if not present.
        :type atom: object

        :return: A tuple containing the atom's metadata: index, name, id, element, atom's
            type (determined by its residue name and atom name), residue name, residue index,
            chain identifier, and segment identifier. Chain and segment identifiers default
            to empty strings if the respective attributes are not found in the input atom.
        :rtype: Tuple[int, str, int, str, str, str, int, str, str]
        """
        chain = atom.chainID if hasattr(atom, "chainID") else ""
        segid = atom.segid if hasattr(atom, "segid") else ""

        atom_type = _classify_type(atom.resname, atom.name)

        return (
            int(atom.index),
            atom.name,
            int(atom.id),
            atom.element,
            atom_type,
            atom.resname,
            int(atom.resid),
            chain,
            segid,
        )

    # Precompute ring features
    g1C, g1N = _precompute_rings(group1_aromatic_rings)
    g2C, g2N = _precompute_rings(group2_aromatic_rings)

    contacts: List[Tuple] = []

    def process(ring_group_id: int,
                rings: List["AtomGroup"], C: np.ndarray, N: np.ndarray,
                cations: "AtomGroup"):
        """
        Processes ring-cation interactions based on geometrical and distance criteria.

        This function identifies and evaluates potential interactions between rings and
        cations. It first applies a distance prefilter using periodic boundary conditions
        (PBC) to determine which cations are within a specified cutoff distance from
        ring centroids. Then, it calculates the angle between the normal direction at
        the ring centroid and the vector pointing to the cation. Only interactions
        meeting angle and distance criteria are considered valid. Outputs metadata
        describing the identified interactions.

        :param ring_group_id: Identifier for the ring group to determine the role
            assignment of group 1 and group 2 in the interaction metadata.
        :type ring_group_id: int
        :param rings: List of ring atom groups, where each group contains information
            related to its ring structure and metadata.
        :type rings: List[AtomGroup]
        :param C: Array of coordinates representing the centroids of the rings.
        :type C: numpy.ndarray
        :param N: Array of normal vectors corresponding to the rings' centroids.
        :type N: numpy.ndarray
        :param cations: Atom group containing information on cationic species and
            positions.
        :type cations: AtomGroup
        :return: None
        """

        if C.size == 0 or cations is None or len(cations) == 0:
            return

        # Distance prefilter (centroid-to-cation), PBC-aware
        dmat = distance_array(C, cations.positions, box=box)  # (n_rings, n_cations)
        ii, jj = np.where(dmat <= cutoff)
        if ii.size == 0:
            return

        # Angle at the centroid between normal direction and centroid→cation vector
        # Use pseudo point: centroid + normal (scale irrelevant for angle)
        eps = 1e-3 # Å
        P = C + N * eps
        A = P[ii]
        B = C[ii]
        C3 = cations.positions[jj]
        ang_rad = calc_angles(A, B, C3, box=box)
        ang_deg = np.degrees(ang_rad)
        ang_deg = np.minimum(ang_deg, 180.0 - ang_deg)  # fold to [0, 90]

        keep = ang_deg <= max_angle
        if not np.any(keep):
            return

        # Emit: GROUP1 FIRST (role-aware), then GROUP2
        for (i, j, theta, dist) in zip(ii[keep], jj[keep], ang_deg[keep], dmat[ii[keep], jj[keep]]):
            ring = rings[int(i)]
            cat = cations[int(j)]

            ring_fields = _ring_meta_str(ring)  # 9 fields (comma-joined strings)
            cation_fields = cation_meta(cat)  # 7 fields (atom meta)

            if ring_group_id == 1:
                # This pass: G1 supplies RING, G2 supplies CATION
                g1_role, g1_payload = "ring", ring_fields
                g2_role, g2_payload = "cation", cation_fields
            else:
                # This pass: G2 supplies RING, G1 supplies CATION
                g1_role, g1_payload = "cation", cation_fields
                g2_role, g2_payload = "ring", ring_fields


            contacts.append((
                "PI-CATION",
                g1_role,        # MUST be position 1
                *g1_payload,    # 9 fields
                *g2_payload,    # 9 fields
                float(dist),
                float(theta),
                "cation-above-face",
            ))

    # Pass 1: G1 rings vs G2 cations  (G1 is RING)
    process(1, group1_aromatic_rings, g1C, g1N, group2_positive)

    # Pass 2: G2 rings vs G1 cations  (G1 is CATION)
    process(2, group2_aromatic_rings, g2C, g2N, group1_positive, )

    return tuple(contacts)


def calculate_metal_contacts(group1_metal_atoms: AtomGroup, group2_metal_atoms: AtomGroup,
                             group1_acceptors: AtomGroup, group1_donors: AtomGroup,
                             group2_acceptors: AtomGroup, group2_donors: AtomGroup,
                             box: np.ndarray | None,
                             cutoff: float = 2.8) -> Tuple[Tuple, ...]:
    """
    Calculates and records metal contacts between two groups of metal atoms and their
    potential acceptor/donor partners within a specified cutoff distance.

    This function identifies interactions between metal atoms (from two predefined groups)
    and their partner atoms (acceptors and donors) based on spatial proximity. Each identified
    interaction is labeled with metadata, including atom indices, names, residues, and interaction
    type. The result is a tuple of interactions logged in a structured format.

    **Contact definition:**
        - Interaction is considered if metal–partner distance ≤ `cutoff` (Å).
        - Metals are provided in `group*_metal_atoms`.
        - Partners are provided in `group*_acceptors` or `group*_donors`.
        - Distance calculation is minimum-image and box-aware.
        - No angular or geometry filtering is applied—this is purely distance-based.

    **Directions evaluated:**
        1. Group 1 metals → Group 2 acceptors   (role: "acceptor")
        2. Group 1 metals → Group 2 donors      (role: "donor")
        3. Group 2 metals → Group 1 acceptors   (role: "acceptor")
        4. Group 2 metals → Group 1 donors      (role: "donor")

    :param group1_metal_atoms: Group of metal atoms in the first set to evaluate interactions.
    :type group1_metal_atoms: AtomGroup
    :param group2_metal_atoms: Group of metal atoms in the second set to evaluate interactions.
    :type group2_metal_atoms: AtomGroup
    :param group1_acceptors: Group of acceptor atoms in the first group for interactions.
    :type group1_acceptors: AtomGroup
    :param group1_donors: Group of donor atoms in the first group for interactions.
    :type group1_donors: AtomGroup
    :param group2_acceptors: Group of acceptor atoms in the second group for interactions.
    :type group2_acceptors: AtomGroup
    :param group2_donors: Group of donor atoms in the second group for interactions.
    :type group2_donors: AtomGroup
    :param box: Periodic boundary conditions to apply during distance calculation.
    :type box: numpy.ndarray
    :param cutoff: Distance cutoff in Ångstroms for identifying interactions.
    :type cutoff: float
    :return: A tuple of tuples, each corresponding to an interaction with the following
        details:
        - Interaction type (str): "METAL-CONTACT"
        - Metal atom details (tuple): Index, name, ID, element, residue name, residue ID,
          chain ID, segment ID
        - Partner atom details (tuple): Same set of details as for the metal atom
        - Distance (float): Calculated distance between the interacting pair
        - Partner role (str): "acceptor" or "donor"
        - Interaction direction label (str)
    :rtype: tuple[tuple]

    """

    def atomdata(a):
        """
        Extracts and organizes relevant atom data into a structured tuple.

        Given an atom object, this function retrieves various properties, such as its
        index, name, ID, element, type, residue name, residue ID, chain ID, and segment
        ID. The atom type is determined by a classification function based on the
        residue and atom names.

        :param a: Atom object containing data to be extracted.
        :type a: Any object with attributes `index`, `name`, `id`, `element`,
            `resname`, `resid`, `chainID`, and `segid`. Attributes `chainID`
            and `segid` are optional.
        :return: A tuple containing the following atom data:
            - Index (:class:`int`)
            - Name (:class:`str`)
            - ID (:class:`int`)
            - Element (:class:`str`)
            - Type (output of `_classify_type`)
            - Residue name (:class:`str`)
            - Residue ID (:class:`int`)
            - Chain ID (:class:`str`, optional)
            - Segment ID (:class:`str`, optional)
        :rtype: tuple
        """
        chain = a.chainID if hasattr(a, "chainID") else ""
        segid = a.segid if hasattr(a, "segid") else ""

        atom_type = _classify_type(a.resname, a.name)

        return (
            int(a.index),
            a.name,
            int(a.id),
            a.element,
            atom_type,
            a.resname,
            int(a.resid),
            chain,
            segid
        )

    contacts: List[Tuple] = []

    def record_interactions(metals: AtomGroup,
                            partners: AtomGroup,
                            role: str,
                            direction_label: str):
        """
        Records interactions between metal atoms and their partner atoms if they meet a
        specified distance cutoff, logging relevant metadata about each interaction.

        :param metals: A group of metal atoms involved in potential interactions.
        :type metals: AtomGroup
        :param partners: A group of partner atoms interacting with metals.
        :type partners: AtomGroup
        :param role: The role of the partner atoms in the interaction (e.g., "acceptor"
            or "donor").
        :type role: str
        :param direction_label: A label indicating the interaction direction or context
            (e.g., "(G1)-M···A-(G2)").
        :type direction_label: str
        :return: None. The function modifies the global `contacts` list with recorded
            interactions.
        :rtype: NoneType
        """
        if len(metals) == 0 or len(partners) == 0:
            return

        dists = distance_array(metals.positions, partners.positions, box=box)
        ii, jj = np.where(dists <= cutoff)
        for i, j in zip(ii, jj):
            m = metals[i]
            p = partners[j]
            contacts.append((
                "METAL-CONTACT",
                *atomdata(m),  # metal first
                *atomdata(p),  # partner second
                float(dists[i, j]),  # distance Å
                role,  # "acceptor" or "donor" (partner role)
                direction_label  # e.g., "(G1)-M···A-(G2)"
            ))

    # 1M – 2 Acceptors
    record_interactions(group1_metal_atoms, group2_acceptors, "acceptor", "(G1)-M···A-(G2)")
    # 1M – 2 Donors
    record_interactions(group1_metal_atoms, group2_donors, "donor", "(G1)-M···D-(G2)")
    # 2M – 1 Acceptors
    record_interactions(group2_metal_atoms, group1_acceptors, "acceptor", "(G1)-A···M-(G2)")
    # 2M – 1 Donors
    record_interactions(group2_metal_atoms, group1_donors, "donor", "(G1)-D···M-(G2)")

    return tuple(contacts)


def interactions_process_frame(**kwargs) -> List[Tuple[Tuple, ...]]:
    """
    Processes a single trajectory frame and computes intermolecular interactions
    using atom groups and control flags supplied via keyword arguments.

    This function acts as a *frame-level interaction dispatcher*: for each enabled
    interaction type, the corresponding detector is executed and its results are
    collected. No post-processing, filtering, or cross-interaction merging is
    performed at this stage.

    -----------------------------------------------------------------------
    Return structure
    -----------------------------------------------------------------------

    The function returns a list where **each element corresponds to one interaction
    class**, in a fixed order. Each element is a tuple containing zero or more
    interaction records of that type.

        frame_contacts = [
            hydrophobic_contacts,        # Tuple[Tuple, ...]
            hydrogen_bond_contacts,      # Tuple[Tuple, ...]
            water_bridge_contacts,       # Tuple[Tuple, ...]
            pi_stacking_contacts,        # Tuple[Tuple, ...]
            pi_cation_contacts,          # Tuple[Tuple, ...]
            ionic_contacts,              # Tuple[Tuple, ...]
            metal_contacts,              # Tuple[Tuple, ...]
            halogen_contacts,            # Tuple[Tuple, ...]
        ]

    If an interaction type is enabled but no contacts are detected, the
    corresponding entry is an empty tuple ``()``.

    -----------------------------------------------------------------------
    Interaction record format
    -----------------------------------------------------------------------

    Each interaction is represented as a **flat tuple** containing atom-level
    metadata and geometric descriptors. The exact fields depend on the interaction
    type, but all records begin with a string interaction label.

    Example (hydrogen bond):

        (
            "HYDROGEN-BOND",
            atom1_index, atom1_name, atom1_id, atom1_element,
            atom1_resname, atom1_resid, atom1_chain, atom1_segid,
            atom2_index, atom2_name, atom2_id, atom2_element,
            atom2_resname, atom2_resid, atom2_chain, atom2_segid,
            distance, angle1, angle2,
            direction_label
        )

    This flat, schema-stable format is designed for efficient downstream
    serialization (e.g., HDF5, NumPy, TSV) and avoids object allocation inside
    tight frame loops.

    -----------------------------------------------------------------------
    Parameters
    -----------------------------------------------------------------------

    :param kwargs:
        Dictionary supplying all atom groups, simulation box information, and
        boolean control flags required to evaluate interactions for the current
        frame.

        Expected keys include (but are not limited to):

        - Atom groups for group 1 and group 2 (e.g., hydrophobic atoms, donors,
          acceptors, aromatic rings, charged atoms, metals)
        - Periodic boundary box (`box`)
        - Boolean flags such as `disable_hydrophobic`, `disable_hbonds`,
          `disable_pi_stacking`, etc.

        All interaction-specific functions assume their required inputs are
        present if the corresponding interaction is enabled.

    :raises KeyError:
        If required keys for an enabled interaction type are missing from `kwargs`.

    -----------------------------------------------------------------------
    Returns
    -----------------------------------------------------------------------

    :return:
        A list of interaction result blocks for the current frame, where each block
        is a tuple of interaction records belonging to one interaction class.

    :rtype:
        List[Tuple[Tuple, ...]]
    """
    frame_contacts: List = []

    box = kwargs["box"]

    if not kwargs.get("disable_hydrophobic", False):
        frame_contacts.append(
            calculate_hydrophobic_contacts(
                group1=kwargs["hydrophobic_atoms_group1"],
                group2=kwargs["hydrophobic_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_hbonds", False):
        frame_contacts.append(
            calculate_hydrogen_bond_contacts(
                group1_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group1"],
                group2_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group2"],
                group1_donors=kwargs["hydrogen_bond_donor_atoms_group1"],
                group2_donors=kwargs["hydrogen_bond_donor_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_water_bridges", False):
        frame_contacts.append(
            calculate_first_degree_water_bridge_contacts(
                group1_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group1"],
                group2_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group2"],
                group1_donors=kwargs["hydrogen_bond_donor_atoms_group1"],
                group2_donors=kwargs["hydrogen_bond_donor_atoms_group2"],
                box=box,
                u=kwargs["u"],
                second_degree=kwargs.get("second_degree_water_bridges", False)
            )
        )

    if not kwargs.get("disable_pi_stacking", False):
        frame_contacts.append(
            calculate_pi_stacking_contacts(
                group1_aromatic_rings=kwargs["aromatic_atoms_group1"],
                group2_aromatic_rings=kwargs["aromatic_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_pi_cation", False):
        frame_contacts.append(
            calculate_pi_cation_contacts(
                group1_aromatic_rings=kwargs["aromatic_atoms_group1"],
                group2_aromatic_rings=kwargs["aromatic_atoms_group2"],
                group1_positive=kwargs["positive_charged_atoms_group1"],
                group2_positive=kwargs["positive_charged_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_ionic", False):
        frame_contacts.append(
            calculate_ionic_contacts(
                group1_positive=kwargs["positive_charged_atoms_group1"],
                group1_negative=kwargs["negative_charged_atoms_group1"],
                group2_positive=kwargs["positive_charged_atoms_group2"],
                group2_negative=kwargs["negative_charged_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_metal_complexation", False):
        frame_contacts.append(
            calculate_metal_contacts(
                group1_metal_atoms=kwargs["metal_atoms_group1"],
                group2_metal_atoms=kwargs["metal_atoms_group2"],
                group1_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group1"],
                group2_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group2"],
                group1_donors=kwargs["hydrogen_bond_donor_atoms_group1"],
                group2_donors=kwargs["hydrogen_bond_donor_atoms_group2"],
                box=box
            )
        )

    if not kwargs.get("disable_halogen", False):
        frame_contacts.append(
            calculate_halogen_bond_contacts(
                group1_halogens=kwargs["halogen_atoms_group1"],
                group2_halogens=kwargs["halogen_atoms_group2"],
                group1_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group1"],
                group2_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group2"],
                group1_donors=kwargs["hydrogen_bond_donor_atoms_group1"],
                group2_donors=kwargs["hydrogen_bond_donor_atoms_group2"],
                box=box
            )
        )

    return frame_contacts



def deduplicate_interactions(interactions: Tuple[Tuple, ...]) -> Tuple[Tuple, ...]:
    """
    Deduplicates molecular interaction records based on their residue or atom identity and removes
    redundant entries.

    This function processes a set of molecular interactions and eliminates redundant entries by:
    1. Removing interactions between residues or atoms that are identical.
    2. Removing reverse duplicates (e.g., interaction A-B is considered the same as B-A).

    Interactions are grouped based on their type, such as "PI-STACKING", "PI-CATION", or
    atom-to-atom interactions, and processed accordingly to identify and filter out duplicates.

    :param interactions: A tuple of tuples representing molecular interaction records. Each record is
        expected to hold metadata for a specific interaction type.
    :return: A tuple of tuples containing unique molecular interaction records.
    """

    if not interactions:
        return interactions

    seen = set()
    unique = []

    def first_token(s: str) -> str:
        """
        Extracts the first token from a given comma-separated string. If the input string
        is empty or None, the function returns an empty string.

        :param s: A comma-separated string from which the first token will be extracted.
        :type s: str
        :return: The first token extracted from the given string, or an empty string
            if the input is None or empty.
        :rtype: str
        """
        return str(s).split(",")[0] if s else ""

    for rec in interactions:
        label = rec[0]

        # PI-STACKING
        if label == "PI-STACKING":
            try:
                r1 = rec[1:10]
                r2 = rec[10:19]

                # residue identity
                res1 = (
                    first_token(r1[5]),  # resname
                    first_token(r1[6]),  # resid
                    first_token(r1[7]),  # chain
                    first_token(r1[8]),  # segid
                )
                res2 = (
                    first_token(r2[5]),
                    first_token(r2[6]),
                    first_token(r2[7]),
                    first_token(r2[8]),
                )

                # 1) Remove same residue
                if res1 == res2:
                    continue

                # Use full ring metadata string as identity
                ring1_id = "|".join(map(str, r1))
                ring2_id = "|".join(map(str, r2))

                key = (label, tuple(sorted((ring1_id, ring2_id))))

            except Exception as exc:
                logger.warning("Skipping malformed PI-STACKING record during deduplication: %s", exc)
                continue

        # PI-CATION
        elif label == "PI-CATION":
            try:
                g1_role = rec[1]
                g1 = rec[2:11]
                g2 = rec[11:20]

                if g1_role == "ring":
                    ring = g1
                    atom = g2
                else:
                    ring = g2
                    atom = g1

                ring_res = (
                    first_token(ring[5]),
                    first_token(ring[6]),
                    first_token(ring[7]),
                    first_token(ring[8]),
                )
                atom_res = (
                    first_token(atom[5]),
                    first_token(atom[6]),
                    first_token(atom[7]),
                    first_token(atom[8]),
                )

                # 1) Remove same residue
                if ring_res == atom_res:
                    continue

                ring_id = "|".join(map(str, ring))
                atom_id = "|".join(map(str, atom))

                key = (label, tuple(sorted((ring_id, atom_id))))

            except Exception as exc:
                logger.warning("Skipping malformed PI-CATION record during deduplication: %s", exc)
                continue

        # OTHER INTERACTIONS (atom–atom)
        else:
            try:
                a1 = rec[1:10]
                a2 = rec[10:19]

                res1 = (
                    first_token(a1[5]),
                    first_token(a1[6]),
                    first_token(a1[7]),
                    first_token(a1[8]),
                )
                res2 = (
                    first_token(a2[5]),
                    first_token(a2[6]),
                    first_token(a2[7]),
                    first_token(a2[8]),
                )

                # 1) Remove same residue
                if res1 == res2:
                    continue

                atom1_id = "|".join(map(str, a1))
                atom2_id = "|".join(map(str, a2))

                key = (label, tuple(sorted((atom1_id, atom2_id))))

            except Exception as exc:
                logger.warning("Skipping malformed %r record during deduplication: %s", label, exc)
                continue

        # 2) Remove reverse duplicates
        if key in seen:
            continue

        seen.add(key)
        unique.append(rec)

    return tuple(unique)


def extract_mode_key(rec: Tuple, *, debug: bool = False) -> Tuple | None:
    """
    Extracts a structured key from a given record based on interaction type.

    This function processes tuple-based interaction records and organizes
    data into a tuple format that can be used for further analysis. It
    handles various interaction types such as ATOM-ATOM interactions,
    WATER-BRIDGE, PI-STACKING, and PI-CATION. The function selectively
    processes and formats these interactions according to specific rules
    defined for each interaction type.

    :param rec: A tuple representing an interaction record. The structure
        of the tuple is assumed to align with the expected format for
        interaction records.
    :type rec: Tuple
    :param debug: A boolean indicating whether to enable debug mode. When
        set to True, debug messages will be printed to assist in tracing
        the function's logic and output.
    :type debug: bool, optional
    :return: A tuple containing residue representations and interaction
        type for the recognized interactions, or None if no valid key could
        be extracted.
    :rtype: Tuple | None
    """
    label = rec[0]

    def residue(resname, resid, chain, segid):
        """
        Creates and returns a residue representation as a tuple based on the provided input
        parameters. Each input value represents a specific property of the residue, allowing for
        flexible and consistent residue identification.

        :param resname: The residue name, typically a string identifier of the residue.
        :type resname: str
        :param resid: The residue ID, an integer or a value convertible to an integer that
            uniquely identifies the residue.
        :type resid: int
        :param chain: The chain identifier, which can be an empty string or a string value
            representing the chain.
        :type chain: str
        :param segid: The segment ID, which can be an empty string or a string value representing
            the segment.
        :type segid: str
        :return: A tuple containing the residue name, residue ID, chain, and segment ID,
            formatted for consistent residue representation.
        :rtype: tuple
        """
        return (resname, int(resid), chain or "", segid or "")

    # ATOM–ATOM interactions (preserve order)
    if label in {
        "HYDROPHOBIC",
        "HYDROGEN-BOND",
        "IONIC",
        "HALOGEN-BOND",
        "METAL-CONTACT",
    }:
        r1 = residue(rec[6], rec[7], rec[8], rec[9])
        r2 = residue(rec[15], rec[16], rec[17], rec[18])

        key = (r1, r2, label)   # ← removed sorting
        logger.trace("extract_mode_key ATOM–ATOM -> %s", key)
        return key

    # WATER BRIDGE (preserve order, keep WATER-BRIDGE-1 / WATER-BRIDGE-2 distinct)
    if label.startswith("WATER-BRIDGE"):
        r1 = residue(rec[6], rec[7], rec[8], rec[9])
        r2 = residue(rec[15], rec[16], rec[17], rec[18])

        key = (r1, r2, label)
        logger.trace("extract_mode_key %s -> %s", label, key)
        return key

    # PI–STACKING (ring–ring, preserve order)
    if label == "PI-STACKING":
        r1 = residue(
            rec[6].split(",")[0],
            rec[7].split(",")[0],
            rec[8].split(",")[0],
            rec[9].split(",")[0],
        )

        r2 = residue(
            rec[15].split(",")[0],
            rec[16].split(",")[0],
            rec[17].split(",")[0],
            rec[18].split(",")[0],
        )

        if r1 == r2:
            logger.trace("extract_mode_key PI-STACKING self-interaction dropped: %s", r1)
            return None

        key = (r1, r2, "PI-STACKING")   # ← removed sorting
        logger.trace("extract_mode_key PI-STACKING -> %s", key)
        return key

    # PI–CATION (KEEP STORED ORDER: payload1 then payload2)
    if label == "PI-CATION":
        g1_role = rec[1]

        def ring_residue(base):
            return residue(
                rec[base + 5].split(",")[0],
                rec[base + 6].split(",")[0],
                rec[base + 7].split(",")[0],
                rec[base + 8].split(",")[0],
            )

        def cation_residue(base):
            return residue(
                rec[base + 5],
                rec[base + 6],
                rec[base + 7],
                rec[base + 8],
            )

        # Parse BOTH payload residues correctly depending on what payload1 is
        # payload1 always starts at base=2, payload2 always starts at base=11
        if g1_role == "ring":
            payload1_res = ring_residue(2)
            payload2_res = cation_residue(11)
        else:
            payload1_res = cation_residue(2)
            payload2_res = ring_residue(11)

        if payload1_res == payload2_res:
            logger.trace("extract_mode_key PI-CATION self-interaction dropped: %s", payload1_res)
            return None

        # CRITICAL: keep stored order (group1 first, group2 second)
        key = (payload1_res, payload2_res, "PI-CATION")
        logger.trace("extract_mode_key PI-CATION (stored-order) -> %s", key)
        return key

    logger.trace("extract_mode_key IGNORED: %s", rec)
    return None


def hbonds_process_frame(**kwargs) -> List[Tuple[Tuple, ...]]:
    """
    Processes a single frame to determine hydrogen bonding interactions between
    specified groups of atoms within the provided simulation box. The function
    invokes a calculation utility to evaluate hydrogen bond contacts and collects
    the results accordingly.

    :param kwargs: Arbitrary keyword arguments that include:
        - 'box': Simulation box dimensions and specifications.
        - 'hydrogen_bond_acceptor_atoms_group1': Atoms that act as acceptors in group 1.
        - 'hydrogen_bond_acceptor_atoms_group2': Atoms that act as acceptors in group 2.
        - 'hydrogen_bond_donor_atoms_group1': Atoms that act as donors in group 1.
        - 'hydrogen_bond_donor_atoms_group2': Atoms that act as donors in group 2.
    :return: A list of tuples representing hydrogen bonding contacts for the frame.
    :rtype: List[Tuple[Tuple, ...]]
    """
    frame_contacts: List = []
    box = kwargs["box"]
    frame_contacts.append(
            calculate_hydrogen_bond_contacts(
                group1_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group1"],
                group2_acceptors=kwargs["hydrogen_bond_acceptor_atoms_group2"],
                group1_donors=kwargs["hydrogen_bond_donor_atoms_group1"],
                group2_donors=kwargs["hydrogen_bond_donor_atoms_group2"],
                box=box)
    )

    return frame_contacts
