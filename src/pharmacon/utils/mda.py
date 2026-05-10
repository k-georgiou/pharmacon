"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

This module provides utility functions for managing chemical element determination,
atom naming conventions, and MDAnalysis Universe creation with optional trajectory
transformations.

Functions in this module primarily focus on interpreting molecular data from topology
and trajectory files, assigning chemical elements, and handling periodic boundary
conditions in molecular dynamics simulations. These utilities are designed for enhanced
compatibility with biomolecular and chemical simulation analysis.
"""

import re
import time
import numpy as np
import periodictable
from rdkit import Chem
from pathlib import Path
import MDAnalysis as Mda
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe
from typing import Optional, List, Tuple, Dict, Union
from MDAnalysis.transformations import unwrap as tf_unwrap, wrap as tf_wrap, center_in_box


from pharmacon.logger import get_logger, PharmaconLogger
from pharmacon.constants import (WATER_RESIDUES, ION_RES_TO_ELT, AA3, NA_RES, PROTEIN_CARBONS, SIDECHAIN_OXY,
                                 SIDECHAIN_NITRO, NUC_OXY, TWO_LETTER, METAL_ELEMENTS)


import warnings




__all__ = [
    "logger",
    "create_universe",
    "convert_mda_to_rdkit",
    "has_residues",
    "has_atoms",
    "require_atoms",
    "require_residues",
]


warnings.filterwarnings("ignore")


logger: PharmaconLogger = get_logger(__name__)


def _guess_element_from_mass(mass: float, threshold: float = 0.3) -> str:
    """
    Guesses the chemical element corresponding to a given mass value. This function iterates
    through the elements in the periodic table, comparing their average mass with the
    provided mass. If the absolute difference between a periodic table element's average mass
    and the provided mass is less than or equal to the given threshold and is the smallest
    difference found, that element is selected. If no match is found within the threshold,
    it returns the placeholder "X". Invalid input or errors during processing will also
    result in a return value of "X".

    :param mass: The mass value to compare against the periodic table elements.
    :type mass: float
    :param threshold: The maximum allowed absolute difference between the input mass and
        an element's average mass for a match (default is 0.3).
    :type threshold: float
    :return: The symbol of the element whose mass best matches the input mass, or "X" if no
        match is found or in case of an error.
    :rtype: str
    """

    logger.trace(
        "Guessing element from mass (mass=%.5f, threshold=%.3f)",
        mass,
        threshold,
    )

    closest_element = None
    smallest_diff = float("inf")

    try:
        m = float(mass)
    except Exception as e:
        logger.debug(
            "Mass value could not be converted to float (mass=%r); returning 'X'.",
            mass,
        )
        logger.trace("Mass conversion exception: %r", e)
        return "X"

    for element in periodictable.elements:
        try:
            avg_mass = element.mass
            if avg_mass is None:
                continue

            diff = abs(m - float(avg_mass))

            if diff <= threshold and diff < smallest_diff:
                smallest_diff = diff
                closest_element = element.symbol

        except Exception as e:
            logger.trace(
                "Skipping element during mass comparison (symbol=%s): %r",
                getattr(element, "symbol", "?"),
                e,
            )
            continue

    if closest_element:
        logger.trace(
            "Element match found for mass %.5f → %s (Δ=%.5f)",
            m,
            closest_element,
            smallest_diff,
        )
        return closest_element

    logger.debug(
        "No element match within threshold for mass %.5f (threshold=%.3f); returning 'X'.",
        m,
        threshold,
    )

    return "X"


def _canonical(sym: str) -> Optional[str]:
    """
    Returns the canonical element symbol of a given chemical symbol if it exists in
    the periodic table. The input symbol is processed to ensure the correct case
    format, where the first character is uppercase and the rest are lowercase,
    before checking its validity.

    :param sym: The input chemical symbol as a string.
    :type sym: str
    :return: The canonical chemical element symbol in proper case format if
        it exists in the periodic table; otherwise returns None.
    :rtype: Optional[str]
    """
    if not sym:
        return None
    cand = sym[0].upper() + sym[1:].lower()
    return cand if hasattr(periodictable.elements, cand) else None


def _guess_element_from_name(atom_name: str, res_name: Optional[str] = None,
                             element_hint: Optional[str] = None) -> str:
    """
    Determines the best-guess element symbol based on the provided atom name, residue name,
    and optional element hint. The function uses a series of heuristic approaches such as
    explicit hint, residue type, biomolecule conventions, and element naming rules in order
    to infer the most likely element symbol.

    :param atom_name: The name of the atom, typically provided in molecular structures.
    :type atom_name: str
    :param res_name: The name of the residue associated with the atom, often available
        in biomolecular contexts (e.g., proteins or nucleic acids). Can be None if
        not applicable.
    :type res_name: Optional[str]
    :param element_hint: An optional explicit hint about the element symbol; may override
        other heuristics.
    :type element_hint: Optional[str]
    :return: The inferred element symbol. Returns "X" if no suitable element can be determined.
    :rtype: str
    """

    # 0) trivial
    if not atom_name:
        return "X"

    # 1) explicit hint wins
    if element_hint:
        hint = re.sub(r"[^A-Za-z]", "", str(element_hint)).upper()
        can = _canonical(hint)
        if can:
            return can

    # normalize names
    res = (res_name or "").strip().upper()
    # keep letters for symbol logic; keep digits for pattern discrimination if needed
    letters = re.sub(r"[^A-Za-z]", "", atom_name).upper()

    # 2) water residues
    if res in WATER_RESIDUES:
        if letters.startswith("O"):
            return "O"
        if letters.startswith("H"):
            return "H"
        return "X"

    # 3) ion residues
    if res in ION_RES_TO_ELT:
        return ION_RES_TO_ELT[res]

    # 4) biomolecule naming conventions
    if res in AA3 or res in NA_RES or res == "MSE" or res == "SEC":
        if letters[:2] in PROTEIN_CARBONS:
            return "C"
        if letters[:2] in SIDECHAIN_OXY or letters[:2] in NUC_OXY:
            return "O"
        if letters[:2] in SIDECHAIN_NITRO:
            return "N"
        if res == "MSE" and letters.startswith("SE"):
            return "Se"
        if res == "SEC" and letters.startswith("SE"):
            return "Se"
        if letters and letters[0] in {"C", "H", "N", "O", "S", "P"}:
            return letters[0]

    # 5) two-letter true element at start, except protein CA
    if len(letters) >= 2:
        head2 = letters[:2]
        if head2 in TWO_LETTER:
            if head2 == "CA" and (res in AA3 or res in NA_RES):
                return "C"
            can = _canonical(head2)
            if can:
                return can

    # 6) single-letter fast path
    if letters and letters[0] in {"C", "H", "N", "O", "S", "P", "F", "I", "K"}:
        return letters[0]

    # 7) generic periodic-table prefix fallback
    for n in (2, 1):
        sym = letters[:n]
        can = _canonical(sym)
        if can:
            return can

    return "X"


def create_universe(topology: str | Path,
                    trajectory: Optional[str | Path] = None,
                    add_elements: bool = False,
                    wrap: bool = False,
                    unwrap: bool = False,
                    compound: str = "atoms") -> Tuple[Mda.Universe, List[str]]:
    """
    Creates and returns an MDAnalysis Universe object and optionally assigns elements to atoms.

    This function is designed to handle reading molecular dynamics topologies and trajectories
    using MDAnalysis and provides functionality for applying trajectory transformations
    (unwrap, center, wrap). It also attempts to ensure that atomic elements are properly
    guessed or provided, which is important for downstream analysis requiring chemical data.

    :param topology: Path to the topology file that defines the molecular system.
    :param trajectory: Optional path to the trajectory file containing simulation frames.
    :param add_elements: Whether to assign or guess elements for atoms in the system.
    :param wrap: Enables wrapping atoms back into the primary simulation box.
    :param unwrap: Enables unwrapping molecules across periodic boundary conditions (PBC).
    :param compound: Defines the type of wrapping compound when `wrap` is enabled
                     (e.g., "atoms", "residues", or "fragments").
    :return: A tuple containing the created MDAnalysis Universe object and the list of
             element symbols (empty if `add_elements` is False).
    """

    logger.trace(
        "create_universe called (topology=%r, trajectory=%r, add_elements=%s, wrap=%s, unwrap=%s, compound=%s)",
        str(topology),
        None if trajectory is None else str(trajectory),
        add_elements,
        wrap,
        unwrap,
        compound,
    )

    elements: List[str] = []

    # --- construct universe ---
    try:
        start_t = time.time()

        if trajectory is not None:
            logger.trace("Constructing Universe with topology and trajectory.")
            u = Mda.Universe(topology, trajectory)
            logger.info(
                "MDAnalysis Universe created (topology=%s, trajectory=%s)",
                Path(topology).name,
                Path(trajectory).name,
            )
        else:
            logger.trace("Constructing Universe with topology only.")
            u = Mda.Universe(topology)
            logger.info(
                "MDAnalysis Universe created (topology=%s, no trajectory)",
                Path(topology).name,
            )

        end_t = time.time()

        logger.debug(
            "Universe creation time: %.3f s",
            end_t - start_t,
        )

        logger.trace(
            "Universe stats: atoms=%d residues=%d segments=%d",
            getattr(u.atoms, "n_atoms", -1),
            getattr(u.residues, "n_residues", -1),
            getattr(u.segments, "n_segments", -1),
        )

    except Exception as e:
        logger.error("Universe construction failed: %r", e)
        raise RuntimeError(
            f"Failed to load topology/trajectory into MDAnalysis: {e}"
        ) from e

    # --- add optional coordinate transformations ---
    # order: unwrap -> center -> wrap
    try:
        if getattr(u, "trajectory", None) is not None and u.trajectory.n_frames >= 1:

            logger.trace(
                "Trajectory detected (frames=%d). Evaluating transformations.",
                u.trajectory.n_frames,
            )

            if not (unwrap or wrap):
                logger.debug(
                    "No trajectory transformations requested (wrap=False, unwrap=False); skipping transformation setup."
                )
                logger.trace(
                    "Centering is also skipped because no wrapping-related transformations were requested."
                )
            else:
                tfms = []

                if unwrap:
                    logger.debug(
                        "Enabling trajectory transformation: unwrap (make molecules whole across PBC)."
                    )
                    tfms.append(tf_unwrap(u.atoms))
                    logger.trace("unwrap transformation queued.")

                # select centering reference only if a transformation will be applied
                try:
                    logger.debug("Selecting reference AtomGroup for centering (protein preferred).")

                    protein = u.select_atoms("protein")
                    logger.trace("Protein selection size: %d", protein.n_atoms)

                    ref = protein if protein.n_atoms > 0 else u.atoms

                    logger.debug(
                        "Centering reference selected: %d atoms.",
                        ref.n_atoms,
                    )

                except Exception as e:
                    logger.debug(
                        "Protein selection failed; falling back to all atoms for centering."
                    )
                    logger.trace("Protein selection exception: %r", e)

                    ref = u.atoms
                    logger.debug(
                        "Centering reference selected: %d atoms.",
                        ref.n_atoms,
                    )

                tfms.append(center_in_box(ref, center="geometry"))
                logger.debug(
                    "Added trajectory transformation: center in box (geometry)."
                )
                logger.trace("center transformation queued.")

                if wrap:
                    logger.debug(
                        "Enabling trajectory transformation: wrap (compound=%s).",
                        compound,
                    )
                    tfms.append(tf_wrap(u.atoms, compound=compound))
                    logger.trace("wrap transformation queued.")

                logger.debug(
                    "Applying trajectory transformations in order: unwrap → center → wrap."
                )
                logger.trace(
                    "Number of transformations applied: %d",
                    len(tfms),
                )

                u.trajectory.add_transformations(*tfms)

                logger.debug(
                    "Trajectory transformations successfully applied."
                )

        else:
            logger.trace(
                "No trajectory frames available; skipping transformations."
            )

    except Exception as e:
        # transformations are convenience; don't fail hard
        logger.warning(
            "Trajectory transformations could not be applied (continuing without them): %s",
            e,
        )
        logger.trace("Transformation error: %r", e)

    if not add_elements:
        logger.debug(
            "Element assignment disabled; returning Universe without elements."
        )
        return u, elements

    # --- ensure/guess elements ---
    try:
        logger.debug("Checking for existing atomic element information.")

        have_attr = hasattr(u.atoms, "elements")
        current = None

        logger.trace("Atoms have elements attribute: %s", have_attr)

        if have_attr:
            try:
                current = np.asarray(u.atoms.elements, dtype=object)
                logger.debug("Existing elements attribute detected on atoms.")
                logger.trace("Elements array size: %d", current.size)
            except Exception as e:
                logger.debug("Elements attribute present but could not be read.")
                logger.trace("Reading elements failed: %r", e)

        if current is not None and current.size > 0:

            nonempty_mask = np.array(
                [bool(x) and str(x) != "None" for x in current],
                dtype=bool,
            )

            logger.trace(
                "Existing element coverage: %d/%d",
                int(nonempty_mask.sum()),
                int(current.size),
            )

            if nonempty_mask.all():
                elements = list(map(str, current.tolist()))

                logger.info(
                    "All atoms already have element assignments; skipping guessing."
                )

                return u, elements

            else:
                logger.debug(
                    "Partial element assignments detected; missing values will be guessed."
                )

                elements = list(
                    map(lambda x: "" if x is None else str(x), current.tolist())
                )

        else:
            logger.debug(
                "No usable element data found; initializing empty element array."
            )

            try:
                u.add_TopologyAttr("elements")
                logger.debug("Topology attribute 'elements' added to Universe.")
            except Exception as e:
                logger.debug(
                    "Failed to explicitly add 'elements' topology attribute; proceeding anyway."
                )
                logger.trace("add_TopologyAttr failed: %r", e)

            elements = [""] * u.atoms.n_atoms

        masses = getattr(u.atoms, "masses", None)

        can_use_mass = False
        m = None

        if masses is not None:
            m = np.asarray(masses, dtype=float)

            can_use_mass = (
                m.size == u.atoms.n_atoms
                and np.all(np.isfinite(m))
            )

            logger.trace(
                "Mass validation: size=%d finite=%s usable=%s",
                m.size,
                bool(np.all(np.isfinite(m))),
                can_use_mass,
            )

        if can_use_mass:

            logger.info(
                "Guessing atomic elements from masses (all masses present and valid)."
            )

            guessed = [
                _guess_element_from_mass(float(mm))
                for mm in m
            ]

        else:

            logger.warning(
                "Atomic masses unavailable or invalid; guessing elements from names and residue names."
            )

            names = np.asarray(u.atoms.names, dtype=object)

            try:
                resn = np.asarray(u.atoms.resnames, dtype=object)
            except Exception:
                resn = np.array([None] * len(names), dtype=object)

            guessed = [
                _guess_element_from_name(
                    str(n),
                    None if r is None else str(r),
                )
                for n, r in zip(names, resn)
            ]

        missing_before = sum(1 for x in elements if not x)
        logger.trace("Missing elements before fill: %d", missing_before)

        for i in range(u.atoms.n_atoms):
            if not elements[i]:
                elements[i] = guessed[i]

        missing_after = sum(1 for x in elements if not x)
        logger.trace("Missing elements after fill: %d", missing_after)

        u.atoms.elements = elements

        logger.info(
            "Element assignment completed for %d atoms.",
            u.atoms.n_atoms,
        )

        return u, elements

    except Exception as e:
        logger.error("Element assignment failed: %r", e)

        raise RuntimeError(
            f"Failed to guess/assign elements in Universe: {e}"
        ) from e



def convert_mda_to_rdkit(ag: Mda.core.groups.AtomGroup) -> Tuple[Chem.Mol, Dict[int, int]]:
    """
    Converts an MDAnalysis AtomGroup to an RDKit molecule while preserving a stable
    atom-index mapping between MDAnalysis and RDKit.

    This version intentionally does NOT call ``Chem.AddHs()``. The goal is to keep
    RDKit atom indices aligned with atoms that actually exist in the input
    MDAnalysis AtomGroup, so that downstream SMARTS matches can be mapped back to
    MDAnalysis atoms without introducing synthetic RDKit-only hydrogens.

    :param ag: The AtomGroup object from MDAnalysis containing molecular data
        to be converted.
    :type ag: Mda.core.groups.AtomGroup

    :return: A tuple containing:
        - ``mol``: The resulting sanitized RDKit molecule.
        - ``index_mapping``: A dictionary mapping universe-level atom indices (``atom.ix``)
          to RDKit atom indices.
    :rtype: Tuple[Chem.Mol, Dict[int, int]]

    :raises RuntimeError: If the conversion of the AtomGroup to an RDKit
        molecule fails or if the sanitization of the RDKit molecule fails.

    .. note::
        This function includes additional processing to disconnect bonds to metals
        and adjusts neighboring formal charges to ensure correctness in RDKit's
        molecular representation. Residue-level information is preserved on a
        best-effort basis using atom properties.

        Hydrogens are NOT added with ``Chem.AddHs()`` because doing so creates
        extra RDKit atoms that do not exist in MDAnalysis and therefore cannot be
        mapped back to the original AtomGroup.
    """

    logger.debug("Converting MDAnalysis AtomGroup to RDKit molecule.")
    logger.trace(
        "convert_mda_to_rdkit called for AtomGroup: n_atoms=%d, n_residues=%s, universe_atoms=%d",
        ag.n_atoms,
        getattr(ag.residues, "n_residues", "unknown"),
        ag.universe.atoms.n_atoms,
    )

    try:
        logger.trace(
            "Calling MDAnalysis RDKit converter with implicit_hydrogens=False, cache=False."
        )
        # Depending on MDAnalysis version, converters live under different namespaces.
        # If this fails, user likely needs MDAnalysis with RDKit converter installed.
        rdkit_mol = ag.convert_to("RDKIT", implicit_hydrogens=False, cache=False)  # type: ignore[attr-defined]

        logger.debug("MDAnalysis AtomGroup successfully converted to RDKit Mol.")
        logger.trace(
            "Initial RDKit Mol stats: atoms=%d, bonds=%d, hydrogens=%d",
            rdkit_mol.GetNumAtoms(),
            rdkit_mol.GetNumBonds(),
            sum(1 for a in rdkit_mol.GetAtoms() if a.GetAtomicNum() == 1),
        )
    except Exception as e:
        logger.error("MDAnalysis -> RDKit conversion failed: %r", e)
        raise RuntimeError(
            f"Failed to convert AtomGroup to RDKit Mol: {e}"
        ) from e

    # Mapping: MDA order -> RDKit order (assumes preserved order from converter)
    # noinspection PyTypeChecker
    index_mapping = {atom.ix: rd_idx for rd_idx, atom in enumerate(ag)}

    logger.debug("Index mapping created.")
    logger.trace(
        "Index mapping size=%d, atomgroup_size=%d",
        len(index_mapping),
        ag.n_atoms,
    )
    logger.trace(
        "Index mapping assumes RDKit atom order matches AtomGroup iteration order."
    )

    if len(index_mapping) != ag.n_atoms:
        logger.warning(
            "Index mapping size (%d) does not match AtomGroup size (%d).",
            len(index_mapping),
            ag.n_atoms,
        )

    # Work on an editable copy
    rwmol = Chem.RWMol(rdkit_mol)
    logger.trace(
        "Editable RDKit molecule created: atoms=%d, bonds=%d",
        rwmol.GetNumAtoms(),
        rwmol.GetNumBonds(),
    )

    # Disconnect bonds to metals and fix neighboring formal charges
    metal_symbols = set(METAL_ELEMENTS)
    metal_indices = [
        atom.GetIdx()
        for atom in rwmol.GetAtoms()
        if atom.GetSymbol().upper() in metal_symbols
    ]

    logger.debug("Scanning RDKit molecule for metal atoms.")
    logger.trace(
        "Detected %d metal atoms for bond disconnection.",
        len(metal_indices),
    )

    removed_bond_count = 0
    adjusted_charge_count = 0

    for metal_idx in metal_indices:
        metal_atom = rwmol.GetAtomWithIdx(metal_idx)
        metal_symbol = metal_atom.GetSymbol()

        neighbors = list(metal_atom.GetNeighbors())
        logger.trace(
            "Processing metal atom idx=%d symbol=%s with %d neighbors.",
            metal_idx,
            metal_symbol,
            len(neighbors),
        )

        for neighbor in neighbors:
            n_idx = neighbor.GetIdx()
            n_atom = rwmol.GetAtomWithIdx(n_idx)
            n_symbol = n_atom.GetSymbol()

            if n_symbol.upper() not in metal_symbols:
                old_charge = n_atom.GetFormalCharge()
                new_charge = old_charge - 1
                n_atom.SetFormalCharge(new_charge)
                adjusted_charge_count += 1

                logger.trace(
                    "Adjusted neighbor formal charge: metal_idx=%d neighbor_idx=%d symbol=%s charge=%d->%d",
                    metal_idx,
                    n_idx,
                    n_symbol,
                    old_charge,
                    new_charge,
                )
            else:
                logger.trace(
                    "Neighbor is also metal; skipping charge adjustment: metal_idx=%d neighbor_idx=%d symbol=%s",
                    metal_idx,
                    n_idx,
                    n_symbol,
                )

            rwmol.RemoveBond(metal_idx, n_idx)
            removed_bond_count += 1

            logger.trace(
                "Removed bond between metal atom %d (%s) and neighbor %d (%s).",
                metal_idx,
                metal_symbol,
                n_idx,
                n_symbol,
            )

    logger.debug(
        "Metal handling complete: removed %d bonds, adjusted %d neighboring formal charges.",
        removed_bond_count,
        adjusted_charge_count,
    )

    # Keep atom mapping consistent with the original AtomGroup.
    # Do NOT add synthetic RDKit hydrogens here, because they cannot be mapped back to MDA.
    mol = rwmol.GetMol()
    logger.debug(
        "Using RDKit molecule without adding extra hydrogens to preserve MDA↔RDKit mapping."
    )
    logger.trace(
        "RDKit Mol before sanitization: atoms=%d, bonds=%d, hydrogens=%d",
        mol.GetNumAtoms(),
        mol.GetNumBonds(),
        sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1),
    )

    try:
        logger.debug("Sanitizing RDKit molecule.")
        Chem.SanitizeMol(mol)
        logger.debug("RDKit molecule sanitization successful.")
    except Exception as e:
        logger.error("RDKit sanitization failed: %r", e)
        raise RuntimeError(
            f"Failed to sanitize RDKit molecule: {e}"
        ) from e

    # Preserve residue info (best-effort)
    logger.debug("Transferring residue metadata from MDAnalysis atoms to RDKit atoms.")
    rdkit_to_mda = {rd_idx: mda_idx for mda_idx, rd_idx in index_mapping.items()}
    logger.trace(
        "Reverse mapping created: size=%d",
        len(rdkit_to_mda),
    )

    transferred_props = 0
    skipped_props = 0

    for atom in mol.GetAtoms():
        rd_idx = atom.GetIdx()
        mda_idx = rdkit_to_mda.get(rd_idx)

        if mda_idx is not None:
            mda_atom = ag.universe.atoms[mda_idx]
            atom.SetProp("resname", str(mda_atom.resname))
            atom.SetProp("resid", str(mda_atom.resid))
            atom.SetProp("segid", str(mda_atom.segid))
            transferred_props += 1
        else:
            skipped_props += 1
            logger.trace(
                "RDKit atom %d has no MDA mapping during metadata transfer.",
                rd_idx,
            )

    logger.debug(
        "Residue metadata transfer complete: assigned=%d, skipped=%d.",
        transferred_props,
        skipped_props,
    )

    logger.debug("RDKit conversion complete.")
    logger.trace(
        "Final RDKit molecule stats: atoms=%d, bonds=%d, hydrogens=%d, mapped_atoms=%d",
        mol.GetNumAtoms(),
        mol.GetNumBonds(),
        sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1),
        len(index_mapping),
    )

    return mol, index_mapping


def has_residues(u: Union[Universe, AtomGroup]) -> bool:
    """
    Check if the given object has residues.

    This function evaluates whether the provided object, which can either
    be a Universe or an AtomGroup, contains any residues. Residues are a
    grouping of atoms, typically representative of molecules or molecular
    subdivisions.

    :param u: The object to check, either a Universe or AtomGroup.
    :return: A boolean indicating if the object has residues.
    :rtype: bool
    """
    return len(u.residues) > 0


def has_atoms(u: Union[Universe, AtomGroup]) -> bool:
    """
    Determines if the provided universe or atom group contains any atoms.

    This function checks whether the given `Universe` or `AtomGroup`
    object has any atoms and returns a boolean value indicating the
    presence of atoms.

    :param u: The input object, which can be either a `Universe`
        or an `AtomGroup`.
    :type u: Union[Universe, AtomGroup]
    :return: `True` if the input object contains one or more atoms;
        otherwise, `False`.
    :rtype: bool
    """

    return len(u.atoms) > 0


def require_atoms(u: Union[Universe, AtomGroup]) -> None:
    """
    Validates that the given object contains atoms.

    This function checks if the provided object has atoms. If the given universe
    or atom group does not contain atoms, a `ValueError` is raised. The type of the
    input object must either be a `Universe` or an `AtomGroup`.

    :param u: The object expected to contain atoms.
    :type u: Union[Universe, AtomGroup]
    :raises ValueError: If the object does not contain atoms.
    :return: None
    """
    if not has_atoms(u):
        raise ValueError("Universe contains no atoms")


def require_residues(u: Union[Universe, AtomGroup]) -> None:
    """
    Checks if a provided Universe or AtomGroup contains residues.

    This function validates the presence of residues in the given object.
    Residues are essential for certain workflows, and the function ensures
    that the provided data structure has them. If not, it raises an exception
    to notify the user.

    :param u: The MDAnalysis Universe or AtomGroup object to be checked.
    :type u: Union[Universe, AtomGroup]
    :raises ValueError: If the Universe or AtomGroup does not contain residues.
    :return: None
    """
    if not has_residues(u):
        raise ValueError("Universe contains no residues")

