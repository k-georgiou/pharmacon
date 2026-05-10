"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Utility functions for detecting specific atom groups in molecular systems using MDAnalysis and RDKit.

This module provides specialized functions for identifying aromatic atoms, metal atoms, and
atoms matching specific SMARTS patterns. The detection involves mapping atom indices between
RDKit and MDAnalysis representations, ensuring uniqueness, and creating AtomGroups for further
analysis.

Loggers and error handling are in place to ensure invalid mappings or malformed patterns are reported.

Imports:
    - MDAnalysis for atom and universe manipulations.
    - RDKit for molecular structure and SMARTS pattern matching.
    - Pharmacon constants for pre-defined patterns and elements.
    - Custom errors for domain-specific validation.

Functions:
    - _detect_aromatic_atoms: Identify aromatic rings using SMARTS patterns and create AtomGroups.
    - _detect_metal_atoms: Detect metal atoms and return them as an AtomGroup.
    - _detect_atoms_by_smarts: Match atoms using SMARTS patterns and optionally filter by element types.
"""


import warnings
from rdkit import Chem
import MDAnalysis as Mda
from MDAnalysis.core.groups import AtomGroup
from typing import Dict, List, Tuple, Optional, Sequence


from pharmacon.logger import get_logger, PharmaconLogger
from pharmacon.constants.smarts import (AROMATIC_PATTERNS, METAL_ELEMENTS, HYDROPHOBIC_PATTERNS, HYDROPHOBIC_ELEMENTS,
                                        HYDROGEN_BOND_ACCEPTOR_PATTERNS, HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
                                        HYDROGEN_BOND_DONOR_PATTERNS,HYDROGEN_BOND_DONOR_ELEMENTS, POSITIVELY_CHARGED_PATTERNS,
                                        NEGATIVELY_CHARGED_PATTERNS, HALOGEN_PATTERNS, HALOGEN_ELEMENTS)





__all__ = [
    "logger",
    "detect_hydrophobic_atoms",
    "detect_hydrogen_bond_acceptor_atoms",
    "detect_hydrogen_bond_donor_atoms",
    "detect_positive_charge_atoms",
    "detect_negative_charge_atoms",
    "detect_halogen_atoms",
    "detect_aromatic_atoms",
    "detect_metal_atoms",
]


warnings.filterwarnings("ignore")


logger: PharmaconLogger = get_logger(__name__)


def _detect_aromatic_atoms(u: Mda.Universe,
                           rdkit_mol: Chem.Mol,
                           mapping: Dict) -> List[AtomGroup]:
    """
    Detects aromatic atoms in a molecular structure by matching SMARTS patterns and mapping indices
    between RDKit and MDAnalysis (MDA).

    This function takes a molecular universe from MDA, a molecular object from RDKit, and a mapping
    of atom indices between these two representations. It then identifies aromatic rings by applying
    SMARTS patterns and ensures the matches correspond correctly between the two frameworks.
    The function validates the mapping direction (RDKit→MDA or MDA→RDKit), handles SMARTS patterns
    and unique aromatic ring detection, and builds AtomGroups corresponding to detected rings.

    :param u: Molecular universe from MDAnalysis.
    :type u: Mda.Universe
    :param rdkit_mol: RDKit molecular object for chemical structure analysis.
    :type rdkit_mol: Chem.Mol
    :param mapping: Dictionary mapping atom indices between RDKit and MDAnalysis representations.
    :type mapping: Dict
    :return: List of MDAnalysis AtomGroups corresponding to identified aromatic rings.
    :rtype: List[AtomGroup]
    """

    logger.debug("Aromatic detection: starting.")
    logger.trace(
        "Aromatic detection input summary: rdkit_atoms=%d universe_atoms=%d mapping_size=%d patterns=%d",
        rdkit_mol.GetNumAtoms(),
        u.atoms.n_atoms,
        len(mapping),
        len(AROMATIC_PATTERNS),
    )

    # Normalize mapping to RDKit -> MDA
    num_rd = rdkit_mol.GetNumAtoms()

    if all(isinstance(k, int) and 0 <= k < num_rd for k in mapping.keys()):
        rdkit_to_mda = {int(k): int(v) for k, v in mapping.items()}
        map_dir = "RDKit→MDA"
    elif all(isinstance(v, int) and 0 <= v < num_rd for v in mapping.values()):
        rdkit_to_mda = {int(v): int(k) for k, v in mapping.items()}
        map_dir = "MDA→RDKit (flipped)"
    else:
        msg = "`mapping` must be RDKit->MDA or MDA->RDKit with int indices"
        logger.error("Aromatic detection failed: invalid mapping format.")
        raise RuntimeError(msg)

    logger.debug(
        "Aromatic detection: mapping normalized (%s), rdkit_atoms=%d",
        map_dir,
        num_rd,
    )
    logger.trace(
        "Aromatic detection: normalized mapping size=%d",
        len(rdkit_to_mda),
    )

    # Match SMARTS; dedupe rings across patterns
    invalid_smarts = []
    unique_rings: Dict[Tuple[int, ...], List[int]] = {}
    total_matches = 0
    skipped_unmapped_matches = 0
    duplicate_ring_matches = 0

    for smarts in AROMATIC_PATTERNS:
        logger.trace("Aromatic detection: compiling SMARTS pattern %r", smarts)
        patt = Chem.MolFromSmarts(smarts)

        if patt is None:
            invalid_smarts.append(smarts)
            logger.warning(
                "Aromatic detection: invalid SMARTS pattern skipped: %r",
                smarts,
            )
            continue

        matches = rdkit_mol.GetSubstructMatches(patt, uniquify=True)
        total_matches += len(matches)

        logger.debug(
            "Aromatic SMARTS %r -> %d matches",
            smarts,
            len(matches),
        )

        for match_i, match in enumerate(matches, start=1):
            rd_indices = tuple(int(i) for i in match)
            logger.trace(
                "Aromatic SMARTS %r match #%d: rdkit_indices=%s",
                smarts,
                match_i,
                rd_indices,
            )

            try:
                ordered = [rdkit_to_mda[int(rd_idx)] for rd_idx in match]
                logger.trace(
                    "Aromatic SMARTS %r match #%d mapped to MDA indices=%s",
                    smarts,
                    match_i,
                    ordered,
                )
            except KeyError:
                skipped_unmapped_matches += 1
                logger.warning(
                    "Aromatic detection: skipping SMARTS %r match #%d due to unmapped RDKit index.",
                    smarts,
                    match_i,
                )
                logger.trace(
                    "Aromatic SMARTS %r match #%d contains unmapped RDKit indices=%s",
                    smarts,
                    match_i,
                    rd_indices,
                )
                continue

            key = tuple(sorted(ordered))

            if key not in unique_rings:
                unique_rings[key] = ordered
                logger.trace(
                    "Aromatic detection: registered unique ring key=%s ordered=%s",
                    key,
                    ordered,
                )
            else:
                duplicate_ring_matches += 1
                logger.trace(
                    "Aromatic detection: duplicate ring ignored key=%s ordered=%s",
                    key,
                    ordered,
                )

    if invalid_smarts:
        msg = f"Invalid aromatic SMARTS: {', '.join(invalid_smarts)}"
        logger.error(
            "Aromatic detection failed: encountered %d invalid SMARTS patterns.",
            len(invalid_smarts),
        )
        raise RuntimeError(msg)

    logger.debug(
        "Aromatic detection summary: total_matches=%d unique_rings=%d skipped_unmapped=%d duplicates=%d",
        total_matches,
        len(unique_rings),
        skipped_unmapped_matches,
        duplicate_ring_matches,
    )

    if unique_rings:
        logger.trace(
            "Aromatic detection: unique ring keys=%s",
            list(unique_rings.keys()),
        )

    # ---- bounds check & build AtomGroups ----
    n_atoms = u.atoms.n_atoms
    rings_ag: List[AtomGroup] = []

    logger.debug(
        "Unique aromatic rings detected: %d",
        len(unique_rings),
    )

    for ring_i, (key, ordered) in enumerate(unique_rings.items(), start=1):
        logger.trace(
            "Aromatic ring #%d candidate: key=%s ordered=%s",
            ring_i,
            key,
            ordered,
        )

        if any(i < 0 or i >= n_atoms for i in key):
            msg = (
                f"Aromatic detection failed: indices out of bounds in ring "
                f"{list(key)} (n_atoms={n_atoms})"
            )
            logger.error(
                "Aromatic detection failed: out-of-bounds indices in ring #%d.",
                ring_i,
            )
            raise RuntimeError(msg)

        try:
            ring_ag = u.atoms[ordered]
            logger.debug(
                "Aromatic ring #%d built: size=%d",
                ring_i,
                len(ordered),
            )
            logger.trace(
                "Aromatic ring #%d MDA indices=%s",
                ring_i,
                ordered,
            )

            if len(ordered) <= 12:
                for atom in ring_ag:
                    logger.trace(
                        "Aromatic ring #%d atom: idx=%d name=%s resname=%s resid=%s segid=%s element=%s",
                        ring_i,
                        int(atom.index),
                        getattr(atom, "name", "?"),
                        getattr(atom, "resname", "?"),
                        getattr(atom, "resid", "?"),
                        getattr(atom, "segid", "?"),
                        getattr(atom, "element", "?"),
                    )
            else:
                logger.trace(
                    "Aromatic ring #%d atom detail logging skipped (size=%d > 12).",
                    ring_i,
                    len(ordered),
                )

            rings_ag.append(ring_ag)

        except Exception as e:
            logger.error(
                "Aromatic detection: failed while creating AtomGroup for ring #%d: %r",
                ring_i,
                e,
            )
            raise RuntimeError(
                f"Aromatic detection failed while creating AtomGroup for ring {ordered}: {e}"
            ) from e

    logger.debug(
        "Aromatic detection complete: returning %d aromatic rings.",
        len(rings_ag),
    )

    return rings_ag


def _detect_metal_atoms(u: Mda.Universe,
                        rdkit_mol: Chem.Mol,
                        mapping: Dict) -> AtomGroup:
    """
    Detect metal atoms in a molecular system based on a mapping between the RDKit and MDAnalysis (MDA)
    universes. This function identifies metal atoms in the RDKit molecule and maps them to their
    corresponding indices in the MDA universe.

    :param u: Molecular dynamics universe from MDAnalysis.
    :param rdkit_mol: RDKit molecule object.
    :param mapping: Dictionary mapping between RDKit atom indices and MDAnalysis atom indices.

    :return: AtomGroup containing all identified metal atoms. If no metals are found,
             an empty AtomGroup is returned.
    """

    logger.debug("Metal detection: starting.")
    logger.trace(
        "Metal detection input summary: rdkit_atoms=%d universe_atoms=%d mapping_size=%d",
        rdkit_mol.GetNumAtoms(),
        u.atoms.n_atoms,
        len(mapping),
    )

    # ---- normalize mapping to RDKit -> MDA ----
    num_rd = rdkit_mol.GetNumAtoms()

    if all(isinstance(k, int) and 0 <= k < num_rd for k in mapping.keys()):
        rdkit_to_mda = {int(k): int(v) for k, v in mapping.items()}
        map_dir = "RDKit→MDA"

    elif all(isinstance(v, int) and 0 <= v < num_rd for v in mapping.values()):
        rdkit_to_mda = {int(v): int(k) for k, v in mapping.items()}
        map_dir = "MDA→RDKit (flipped)"

    else:
        msg = "`mapping` must be RDKit->MDA or MDA->RDKit with int indices"
        logger.error("Metal detection failed: invalid mapping format.")
        raise RuntimeError(msg)

    logger.debug(
        "Metal detection: mapping normalized (%s), rdkit_atoms=%d",
        map_dir,
        num_rd,
    )
    logger.trace(
        "Metal detection: normalized mapping size=%d",
        len(rdkit_to_mda),
    )

    # ---- collect metals ----
    metal_set = {str(z).upper() for z in METAL_ELEMENTS}

    logger.trace(
        "Metal detection: metal element set size=%d",
        len(metal_set),
    )

    metals = set()
    unmapped_metals = 0

    for atom in rdkit_mol.GetAtoms():
        rd_idx = int(atom.GetIdx())
        sym = (atom.GetSymbol() or "").strip().upper()

        logger.trace(
            "Inspecting RDKit atom idx=%d symbol=%s",
            rd_idx,
            sym,
        )

        if sym in metal_set:
            mda_idx = rdkit_to_mda.get(rd_idx)

            if mda_idx is not None:
                metals.add(int(mda_idx))

                logger.trace(
                    "Metal detected: symbol=%s rdkit_idx=%d -> mda_idx=%d",
                    sym,
                    rd_idx,
                    mda_idx,
                )

            else:
                unmapped_metals += 1
                logger.warning(
                    "Metal detected but mapping missing: symbol=%s rdkit_idx=%d",
                    sym,
                    rd_idx,
                )

    idx_sorted = sorted(metals)

    logger.debug(
        "Metal detection summary: metals_found=%d mapped=%d unmapped=%d",
        len(idx_sorted) + unmapped_metals,
        len(idx_sorted),
        unmapped_metals,
    )

    if idx_sorted:
        logger.trace(
            "Metal detection: MDA metal indices=%s",
            idx_sorted,
        )

    # ---- bounds check ----
    n_atoms = u.atoms.n_atoms
    oob = [i for i in idx_sorted if i < 0 or i >= n_atoms]

    if oob:
        msg = f"MDA indices out of bounds: {oob} (n_atoms={n_atoms})"
        logger.error("Metal detection failed: mapped indices out of bounds.")
        raise RuntimeError(msg)

    # ---- build AtomGroup ----
    try:
        ag = u.atoms[idx_sorted]

        logger.debug(
            "Metal detection: final AtomGroup size=%d",
            len(ag),
        )

        if len(ag) <= 25:
            for atom in ag:
                logger.trace(
                    "Metal atom: idx=%d name=%s resname=%s resid=%s segid=%s element=%s",
                    int(atom.index),
                    getattr(atom, "name", "?"),
                    getattr(atom, "resname", "?"),
                    getattr(atom, "resid", "?"),
                    getattr(atom, "segid", "?"),
                    getattr(atom, "element", "?"),
                )
        else:
            logger.trace(
                "Metal detection: atom detail logging skipped (count=%d > 25).",
                len(ag),
            )

        return ag

    except Exception as e:
        logger.error("Metal detection AtomGroup creation failed: %r", e)
        msg = f"while creating AtomGroup: {e}"
        raise RuntimeError(msg) from e


def _detect_atoms_by_smarts(u: Mda.Universe,
                            rdkit_mol: Chem.Mol,
                            mapping: Dict,
                            patterns: Sequence[str],
                            allowed_elements: Optional[List],
                            label: str = "SMARTS") -> AtomGroup:
    """
    Identifies atoms based on SMARTS patterns within a molecular dynamics universe
    using an RDKit molecule, providing an interface for filtering and mapping atom indices.

    :param u: Molecular dynamics universe.
    :type u: Mda.Universe
    :param rdkit_mol: RDKit molecule to be used for SMARTS matching.
    :type rdkit_mol: Chem.Mol
    :param mapping: Mapping of atom indices between RDKit and MDAnalysis (or vice versa).
    :type mapping: Dict
    :param patterns: SMARTS patterns to match against the RDKit molecule.
    :type patterns: Sequence[str]
    :param allowed_elements: Optional list of element symbols allowed for filtering matched atoms.
    :type allowed_elements: Optional[List]
    :param label: Label used for logging and error messages, defaults to "SMARTS".
    :type label: str
    :return: An MDAnalysis AtomGroup containing matched atoms.
    :rtype: AtomGroup
    """

    def _atom_log_fields(atom) -> tuple:
        """
        Extracts and returns specific attributes from an atom object as a tuple. The function
        retrieves the values for `name`, `resname`, `resid`, `chainID`, `segid`, and `element`
        from the `atom` object. If any of these attributes are missing, a default value of "?"
        is used. This facilitates safe evaluation of atom properties while debugging or logging.

        :param atom: Instance containing molecular properties. It is expected to have attributes
            named "name", "resname", "resid", "chainID" (or "chainid" as a fallback),
            "segid", and "element".
        :return: A tuple containing the retrieved property values:
            - atom name,
            - residue name,
            - residue ID,
            - chain ID,
            - segment ID,
            - element.
        """
        return (
            getattr(atom, "name", "?"),
            getattr(atom, "resname", "?"),
            getattr(atom, "resid", "?"),
            getattr(atom, "chainID", getattr(atom, "chainid", "?")),
            getattr(atom, "segid", "?"),
            getattr(atom, "element", "?"),
        )

    logger.debug("%s: starting SMARTS detection.", label)
    logger.trace(
        "%s: input summary: rdkit_atoms=%d, universe_atoms=%d, patterns=%d, allowed_elements=%r",
        label,
        rdkit_mol.GetNumAtoms(),
        u.atoms.n_atoms,
        len(patterns),
        allowed_elements,
    )

    # Normalize mapping to RDKit -> MDA
    num_rd = rdkit_mol.GetNumAtoms()

    if all(isinstance(k, int) and 0 <= k < num_rd for k in mapping.keys()):
        rdkit_to_mda = {int(k): int(v) for k, v in mapping.items()}
        map_dir = "RDKit→MDA"
    elif all(isinstance(v, int) and 0 <= v < num_rd for v in mapping.values()):
        rdkit_to_mda = {int(v): int(k) for k, v in mapping.items()}
        map_dir = "MDA→RDKit (flipped)"
    else:
        msg = (
            f"{label} detection failed: `mapping` must be RDKit->MDA "
            f"or MDA->RDKit with int indices"
        )
        logger.error("%s: invalid mapping format.", label)
        raise RuntimeError(msg)

    logger.debug(
        "%s: mapping normalized (%s), mapping_size=%d",
        label,
        map_dir,
        len(rdkit_to_mda),
    )
    logger.trace("%s: SMARTS patterns=%r", label, list(patterns))

    # match SMARTS
    mda_indices = set()
    invalid_smarts = []
    total_matches = 0
    total_mapped_atoms = 0
    total_unmapped_atoms = 0

    for smarts in patterns:
        logger.trace("%s: compiling SMARTS pattern %r", label, smarts)
        patt = Chem.MolFromSmarts(smarts)

        if patt is None:
            invalid_smarts.append(smarts)
            logger.warning("%s: invalid SMARTS pattern skipped: %r", label, smarts)
            continue

        matches = rdkit_mol.GetSubstructMatches(patt, uniquify=True)
        total_matches += len(matches)

        logger.debug(
            "%s: SMARTS %r -> %d matches",
            label,
            smarts,
            len(matches),
        )

        for match_i, match in enumerate(matches, start=1):
            rd_indices = tuple(int(i) for i in match)
            logger.trace(
                "%s: SMARTS %r match #%d: rdkit_indices=%s",
                label,
                smarts,
                match_i,
                rd_indices,
            )

            mapped_in_match = 0
            unmapped_in_match = 0

            for rd_idx in rd_indices:
                mda_idx = rdkit_to_mda.get(rd_idx)

                if mda_idx is not None:
                    mda_idx = int(mda_idx)
                    mda_indices.add(mda_idx)
                    mapped_in_match += 1
                    total_mapped_atoms += 1

                    atom = u.atoms[mda_idx]
                    atom_name, atom_resname, atom_resid, atom_chainid, atom_segid, atom_element = _atom_log_fields(atom)

                    logger.trace(
                        "%s: SMARTS %r mapped RDKit atom %d -> MDA atom %d "
                        "(name=%s resname=%s resid=%s chainID=%s segid=%s element=%s)",
                        label,
                        smarts,
                        rd_idx,
                        mda_idx,
                        atom_name,
                        atom_resname,
                        atom_resid,
                        atom_chainid,
                        atom_segid,
                        atom_element,
                    )
                else:
                    unmapped_in_match += 1
                    total_unmapped_atoms += 1
                    logger.trace(
                        "%s: SMARTS %r RDKit atom %d has no mapping to MDA",
                        label,
                        smarts,
                        rd_idx,
                    )

            logger.trace(
                "%s: SMARTS %r match #%d summary: mapped=%d unmapped=%d",
                label,
                smarts,
                match_i,
                mapped_in_match,
                unmapped_in_match,
            )

    if invalid_smarts:
        msg = f"Invalid SMARTS patterns for {label}: {', '.join(invalid_smarts)}"
        logger.error(
            "%s: encountered %d invalid SMARTS patterns.",
            label,
            len(invalid_smarts),
        )
        raise RuntimeError(msg)

    idx_list = sorted(mda_indices)

    logger.debug(
        "%s: SMARTS matching summary: total_matches=%d, unique_mda_atoms=%d, mapped_atoms=%d, unmapped_atoms=%d",
        label,
        total_matches,
        len(idx_list),
        total_mapped_atoms,
        total_unmapped_atoms,
    )

    if idx_list:
        logger.trace(
            "%s: matched MDA indices before filtering=%s",
            label,
            idx_list,
        )

        if len(idx_list) <= 25:
            for idx in idx_list:
                atom = u.atoms[idx]
                atom_name, atom_resname, atom_resid, atom_chainid, atom_segid, atom_element = _atom_log_fields(atom)
                logger.trace(
                    "%s: pre-filter atom idx=%d name=%s resname=%s resid=%s chainID=%s segid=%s element=%s",
                    label,
                    idx,
                    atom_name,
                    atom_resname,
                    atom_resid,
                    atom_chainid,
                    atom_segid,
                    atom_element,
                )
        else:
            logger.trace(
                "%s: pre-filter atom detail logging skipped because matched atom count is %d (>25).",
                label,
                len(idx_list),
            )

    # bounds check
    n_atoms = u.atoms.n_atoms
    oob = [i for i in idx_list if i < 0 or i >= n_atoms]

    if oob:
        msg = f"{label} detection failed: MDA indices out of bounds: {oob} (n_atoms={n_atoms})"
        logger.error("%s: out-of-bounds mapped indices detected.", label)
        raise RuntimeError(msg)

    # optional allowed-elements post-filter
    if allowed_elements is not None:
        allowed_set = {str(e).upper() for e in allowed_elements}
        ag_preview = u.atoms[idx_list]
        kept_idx = []

        logger.debug(
            "%s: applying allowed-elements filter with %d allowed values.",
            label,
            len(allowed_set),
        )
        logger.trace("%s: allowed element set=%s", label, sorted(allowed_set))

        for i, a in zip(idx_list, ag_preview):
            atom_name, atom_resname, atom_resid, atom_chainid, atom_segid, atom_element = _atom_log_fields(a)
            atom_element_u = (atom_element or "").upper()
            keep = atom_element_u in allowed_set

            logger.trace(
                "%s: element filter atom idx=%d name=%s resname=%s resid=%s chainID=%s segid=%s element=%s keep=%s",
                label,
                i,
                atom_name,
                atom_resname,
                atom_resid,
                atom_chainid,
                atom_segid,
                atom_element_u,
                keep,
            )

            if keep:
                kept_idx.append(i)

        removed = len(idx_list) - len(kept_idx)
        logger.debug(
            "%s: post-filter kept=%d removed=%d",
            label,
            len(kept_idx),
            removed,
        )
        idx_list = kept_idx

    # return AtomGroup
    try:
        ag = u.atoms[idx_list]
        logger.debug("%s: final matched atoms=%d", label, len(idx_list))

        if idx_list:
            logger.trace("%s: final MDA indices=%s", label, idx_list)

            if len(idx_list) <= 25:
                for atom in ag:
                    atom_idx = int(atom.index)
                    atom_name, atom_resname, atom_resid, atom_chainid, atom_segid, atom_element = _atom_log_fields(atom)

                    logger.trace(
                        "%s: final atom idx=%d name=%s resname=%s resid=%s chainID=%s segid=%s element=%s",
                        label,
                        atom_idx,
                        atom_name,
                        atom_resname,
                        atom_resid,
                        atom_chainid,
                        atom_segid,
                        atom_element,
                    )
            else:
                logger.trace(
                    "%s: final atom detail logging skipped because final atom count is %d (>25).",
                    label,
                    len(idx_list),
                )
        else:
            logger.trace("%s: no atoms matched after filtering.", label)

        return ag

    except Exception as e:
        logger.error("%s: AtomGroup creation failed: %r", label, e)
        msg = f"{label} detection failed while creating AtomGroup: {e}"
        raise RuntimeError(msg) from e


def _log_count(label: str, n: int) -> None:
    """
    Logs the count of a specific label and its corresponding value to the logger.

    This function is intended to log the provided label formatted with a uniform
    width of 40 characters alongside the given numerical value. The primary
    purpose of this method is to facilitate standardized logging for better
    readability and tracking.

    :param label: The label to be logged. It should be descriptive to indicate
        the context or type of the value being logged.
    :type label: str
    :param n: The numerical value associated with the label that needs to be
        logged.
    :type n: int
    :return: None
    """
    logger.info(f"{label:<40} : {n}")


def detect_hydrophobic_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                             label="Hydrophobic atoms in Group") -> AtomGroup:
    """
    Detect hydrophobic atoms based on the specified SMARTS patterns and allowed
    elements.

    This function identifies hydrophobic atoms within the provided Molecular
    Dynamics Analysis (MDA) universe by utilizing the given RDKit molecule,
    mapping, and defined hydrophobic elements and patterns. It also labels
    the identified atoms according to the specified or default label. The
    results include the set of hydrophobic atoms that match the defined criteria.

    :param u: The Molecular Dynamics Analysis (MDA) universe containing the
              molecular data to analyze.
    :type u: Mda.Universe
    :param rdkit_mol: The RDKit molecule object used for SMARTS-based pattern
                      matching.
    :type rdkit_mol: Chem.Mol
    :param mapping: A dictionary mapping RDKit atom indices to corresponding
                    MDA atom indices.
    :type mapping: dict
    :param label: The label to assign to the detected group of hydrophobic
                  atoms. Defaults to "Hydrophobic atoms in Group".
    :type label: str
    :return: The AtomGroup object representing the detected hydrophobic atoms.
    :rtype: AtomGroup
    """
    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=HYDROPHOBIC_PATTERNS,
        allowed_elements=list(HYDROPHOBIC_ELEMENTS),
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_hydrogen_bond_acceptor_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                                        label="Hydrogen Bond Acceptor atoms in Group") -> AtomGroup:
    """
    Detects hydrogen bond acceptor atoms in a molecular dynamics universe.

    This function identifies atoms that serve as hydrogen bond acceptors
    in a given molecular dynamics system, based on provided SMARTS patterns
    and allowed elements.

    :param u: A MDAnalysis Universe object representing the molecular system.
    :param rdkit_mol: An RDKit Mol object containing the molecular structure.
    :param mapping: Dictionary mapping between RDKit Mol and MDAnalysis
        Universe indices.
    :param label: Customized label for the group of detected hydrogen
        bond acceptor atoms. Defaults to "Hydrogen Bond Acceptor atoms in Group".
    :return: An MDAnalysis AtomGroup containing the identified hydrogen bond
        acceptor atoms.
    """
    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=HYDROGEN_BOND_ACCEPTOR_PATTERNS,
        allowed_elements=list(HYDROGEN_BOND_ACCEPTOR_ELEMENTS),
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_hydrogen_bond_donor_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                                     label="Hydrogen Bond Donor atoms in Group") -> AtomGroup:
    """
    Detects hydrogen bond donor atoms in a molecular group.

    This function identifies and labels atoms in a molecule that are potential hydrogen bond donors
    using SMARTS pattern recognition. The identified atoms are returned as an AtomGroup, suitable
    for further analysis or processing.

    :param u: A molecular dynamics Universe object containing the molecular system.
    :param rdkit_mol: An RDKit molecule object representing the molecule for pattern matching.
    :param mapping: A dictionary mapping between RDKit atom indices and MDAnalysis atom indices.
    :param label: A descriptive label for the AtomGroup being created.
    :return: An AtomGroup object containing atoms matching the hydrogen bond donor criteria.
    """
    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=HYDROGEN_BOND_DONOR_PATTERNS,
        allowed_elements=list(HYDROGEN_BOND_DONOR_ELEMENTS),
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_positive_charge_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                                 label="Positive charged atoms in Group") -> AtomGroup:
    """
    Detects positively charged atoms within the provided chemical structure.

    This function identifies atoms with positive charges from the given molecular
    data using specified SMARTS patterns. The positively charged atoms are grouped
    and labeled according to the provided parameters.

    :param u: MDAnalysis Universe object containing molecular topology and trajectory.
    :type u: mda.Universe
    :param rdkit_mol: RDKit molecule object representing the molecular structure.
    :type rdkit_mol: Chem.Mol
    :param mapping: Dictionary mapping between MDAnalysis and RDKit atom indices.
    :type mapping: Dict
    :param label: Label for grouping the positively charged atoms, with default
                  value "Positive charged atoms in Group".
    :return: Group of positively charged atoms matching the specified patterns.
    :rtype: AtomGroup
    """

    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=POSITIVELY_CHARGED_PATTERNS,
        allowed_elements=None,
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_negative_charge_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                                 label="Negative charged atoms in Group") -> AtomGroup:
    """
    Detect and return atoms with negative charge based on predefined SMARTS patterns.

    This function identifies atoms within a molecular system described by the
    MDAnalysis Universe and RDKit molecule that exhibit a negative charge. The
    detection is performed using predefined SMARTS patterns. It supports integration
    of RDKit and MDAnalysis functionality to map detected atoms and return them
    as an MDAnalysis AtomGroup for further analysis or operations.

    :param u: The MDAnalysis Universe object representing the molecular system.
    :type u: Mda.Universe
    :param rdkit_mol: The RDKit Mol object representing the molecular structure.
    :type rdkit_mol: Chem.Mol
    :param mapping: The mapping between RDKit indices and MDAnalysis indices.
    :type mapping: Dict
    :param label: Optional label used for logging detected atoms. The default label
        is "Negative charged atoms in Group".
    :type label: str, optional
    :return: The MDAnalysis AtomGroup containing atoms identified as negatively
        charged based on SMARTS patterns.
    :rtype: AtomGroup
    """

    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=NEGATIVELY_CHARGED_PATTERNS,
        allowed_elements=None,
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_halogen_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                         label="Halogen atoms in Group") -> AtomGroup:
    """
    Detects halogen atoms within a molecular system based on predefined SMARTS patterns
    and allowed elements. The function utilizes RDKit's molecular representation and provides
    a labeled AtomGroup corresponding to the detected halogen atoms.

    :param u: MDAnalysis Universe object representing the molecular system of interest.
    :param rdkit_mol: RDKit Mol object representing the chemical structure of the molecular system.
    :param mapping: A dictionary mapping atoms between the MDAnalysis and RDKit representations.
    :param label: A string label describing the halogen atom group detected, default is
        "Halogen atoms in Group".
    :return: An MDAnalysis AtomGroup containing the detected halogen atoms.
    """
    ag = _detect_atoms_by_smarts(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping,
        patterns=HALOGEN_PATTERNS,
        allowed_elements=list(HALOGEN_ELEMENTS),
        label=label,
    )
    _log_count(label, len(ag))
    return ag


def detect_aromatic_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                          label="Aromatic rings in Group") -> List[AtomGroup]:
    """
    Detects aromatic atoms in a molecular universe and maps them,
    returning a list of associated atom groups.

    This method utilizes an internal function to identify aromatic atoms
    within a molecular universe and links them back to the molecular
    structure provided by RDKit. It logs the count of detected
    aromatic atom groups with the specified label.

    :param u: The MDAnalysis Universe object representing the molecular
        universe to be analyzed.
    :type u: Mda.Universe
    :param rdkit_mol: RDKit molecular object providing the chemical
        structure to be mapped.
    :type rdkit_mol: Chem.Mol
    :param mapping: Dictionary linking MDAnalysis atoms to RDKit atoms.
    :type mapping: Dict
    :param label: Label message for logging detected aromatic atom groups
        (default is "Aromatic rings in Group").
    :type label: str
    :return: A list of MDAnalysis atom groups identified as aromatic.
    :rtype: List[AtomGroup]
    """
    ag = _detect_aromatic_atoms(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping)
    _log_count(label, len(ag))
    return ag


def detect_metal_atoms(u: Mda.Universe, rdkit_mol: Chem.Mol, mapping: Dict,
                       label="Metal atoms in Group") -> AtomGroup:
    """
    Detects and labels metal atoms in a molecular universe.

    This function identifies metal atoms within a molecular dynamics universe,
    based on a given RDKit molecule representation and a mapping between the
    universe and molecule. The detection process segregates metal atoms, allowing
    further analysis or processing.

    :param u: An MDAnalysis universe object representing the molecular dynamics
              simulation data.
    :param rdkit_mol: An RDKit molecule object defining the molecular structure
                      to be compared against the MDAnalysis universe.
    :param mapping: A dictionary mapping atom indices between the MDAnalysis
                    universe and the RDKit molecule.
    :param label: An optional label used in the log message to describe the
                  detected group of metal atoms.
    :return: An MDAnalysis AtomGroup containing the detected metal atoms.
    """
    ag = _detect_metal_atoms(
        u=u,
        rdkit_mol=rdkit_mol,
        mapping=mapping)
    _log_count(label, len(ag))
    return ag
