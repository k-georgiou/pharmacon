"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

This module contains functions for reading and parsing chemical file formats such as
MOL2 and SMILES (.smi), processing molecule data into RDKit molecule objects and
managing molecule metadata such as names.

The module offers robust handling of file parsing errors, detailed logging mechanisms,
and automated fallback strategies for molecule name resolution or sanitization.
It is primarily focused on providing usable RDKit molecule objects from chemical
file formats for further computational chemistry applications.
"""


import re
import uuid
import warnings
from rdkit import Chem
from pathlib import Path
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
from rdkit.Chem import Fragments
from rdkit.Chem import Descriptors
from typing import Dict, List, Tuple
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator as FPG


from pharmacon.logger import get_logger, PharmaconLogger




__all__ = [
    "logger",
    "read_mol2_file",
    "read_smi_file",
    "get_element_counts",
    "count_total_atoms",
    "get_fragment_counts",
    "get_molecular_volume",
    "get_molecular_weight",
    "get_num_rotatable_bonds",
    "get_number_of_rings",
    "get_tpsa",
    "get_number_of_aromatic_rings",
    "generate_morgan_fingerprint",
    "generate_topological_torsion_fingerprint",
    "generate_atom_pair_fingerprint",
    "generate_maccs_keys",
    "get_stereo_centers",
    "get_net_charge",
    "get_logp",
    "fingerprint_to_string",
]


warnings.filterwarnings("ignore")


RDLogger.DisableLog("rdApp.*")


logger: PharmaconLogger = get_logger(__name__)


def read_mol2_file(topology: str | Path) -> Tuple[List[str], List[Chem.Mol]]:
    """
    Reads a MOL2 file and parses its contents into a list of molecule names and RDKit molecule objects.

    This function processes each molecule block in the MOL2 file format, attempts to construct
    the corresponding RDKit molecule, and resolves molecule names based on available metadata.
    Failed molecule conversions are logged and skipped automatically.

    :param topology: Path to the MOL2 file to be read. Can be a string or Path object.
    :type topology: str | Path

    :return: A tuple containing two lists. The first list contains the names of the
             parsed molecules, and the second list contains the corresponding RDKit
             molecule objects.
    :rtype: Tuple[List[str], List[Chem.Mol]]

    :raises RuntimeError: If the MOL2 file cannot be read due to issues such as file
                          not found or parsing errors.
    """
    names: List[str] = []
    molecules: List[Chem.Mol] = []

    def _finalize_block(block_lines: List[str], idx: int, fallback_name: str | None) -> None:
        """
        Finalize the processing of a chemical block by invoking RDKit methods to interpret the
        chemical data and handle naming resolutions for molecules.

        This method processes a MOL2 block, attempts to construct the corresponding RDKit molecule,
        and optionally resolves the molecule's name based on metadata or fallback logic.

        :param block_lines: Lines representing the chemical structure in MOL2 block format.
        :type block_lines: List[str]
        :param idx: Index of the molecule block, typically used for logging or fallback naming.
        :type idx: int
        :param fallback_name: An optional fallback name to use for the molecule if no name is
                              found in the block. Defaults to None.
        :type fallback_name: str | None

        :return: None
        """
        if not block_lines:
            return
        block = "".join(block_lines)
        mol: Chem.Mol | None = None

        # Primary attempt
        try:
            mol = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
        except Exception as e:
            logger.warning(f"MOL2[{idx}] sanitize=True failed: {e}")

        # Fallback route
        if mol is None:
            try:
                mol = Chem.MolFromMol2Block(block, sanitize=False, removeHs=False)
                if mol is not None:
                    try:
                        Chem.SanitizeMol(mol)
                        mol = Chem.RemoveHs(mol)
                    except Exception as e:
                        logger.warning(f"MOL2[{idx}] manual sanitize failed: {e}")
                        mol = None
            except Exception as e:
                logger.warning(f"MOL2[{idx}] parse sanitize=False failed: {e}")
                mol = None

        if mol is None:
            logger.warning(f"Failed to convert molecule at index {idx} to RDKit Mol; skipping.")
            return

        # Name resolution
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else (fallback_name or f"mol_{idx}")
        names.append(name)
        molecules.append(mol)

    # File reading
    block_lines: List[str] = []
    block_index: int = 0
    expect_name_line: bool = False
    pending_name: str | None = None

    try:
        with open(topology, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("@<TRIPOS>MOLECULE"):
                    if block_lines:
                        _finalize_block(block_lines, block_index, pending_name)
                        block_index += 1
                        block_lines = []
                        pending_name = None

                    block_lines.append(line)
                    expect_name_line = True
                    continue

                if expect_name_line:
                    stripped = line.strip()
                    if stripped:
                        pending_name = stripped
                        expect_name_line = False
                    block_lines.append(line)
                else:
                    block_lines.append(line)

        if block_lines:
            _finalize_block(block_lines, block_index, pending_name)

    except Exception as e:
        raise RuntimeError(f"Failed to read MOL2 file '{topology}': {e}")

    return names, molecules


def read_smi_file(topology: str | Path) -> Tuple[List[str], List[Chem.Mol]]:
    """
    Processes a .smi file to generate a list of molecule names and their corresponding
    RDKit molecule objects. The function reads the given .smi file, parses each SMILES
    string, assigns unique names to the molecules, and returns them. If any molecule
    in the file fails to parse, detailed warnings are logged, and the molecule is skipped.

    - 'SMILES  name' or 'SMILES' (name auto-generated)
      - lines starting with '#', ';', '//' are ignored.

    :param topology: The path to the .smi file to be parsed.
    :type topology: str | Path

    :return: A tuple containing a list of molecule names and a list of RDKit molecule objects.
    :rtype: Tuple[List[str], List[Chem.Mol]]
    """
    names: List[str] = []
    mols: List[Chem.Mol] = []

    def _sanitize_name(s: str | None) -> str:
        """
        Sanitizes a given name string by replacing spaces and invalid characters with underscores.
        Ensures the name contains only alphanumeric characters, underscores, periods, or hyphens.
        Stripped of leading/trailing invalid characters.

        :param s: The input string to be sanitized. Can be None.
        :type s: str | None
        :return: A sanitized string with valid characters, or an empty string if the input is invalid or empty.
        :rtype: str
        """
        s = (s or "").strip()
        if not s:
            return ""
        s = re.sub(r"\s+", "_", s)
        s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
        return s.strip("._-") or ""

    def _unique_name(base: str, seen: Dict[str, int]) -> str:
        """
        Generates a unique name based on the given base string and tracks occurrences using
        the provided dictionary. If the base string is empty, a randomly generated 8-character
        UUID hex value is used as the base. Returns a unique name which increments the count
        if the base already exists in the seen dictionary.

        :param base: The base string to derive the unique name from. If empty, a UUID-based
            string will be generated.
        :type base: str
        :param seen: A dictionary that keeps track of occurrences of unique names. Keys are
            the base names, and values are the count of how many times the base has been used.
        :type seen: Dict[str, int]
        :return: A unique name generated based on the base string and the seen dictionary.
        :rtype: str
        """
        if base == "":
            base = uuid.uuid4().hex[:8]
        if base not in seen:
            seen[base] = 0
            return base
        seen[base] += 1
        return f"{base}_{seen[base]}"

    seen_names: Dict[str, int] = {}

    try:
        with open(topology, "r", encoding="utf-8", errors="replace") as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#") or line.startswith(";") or line.startswith("//"):
                    continue

                parts = line.split(maxsplit=1)
                if not parts:
                    continue
                smiles = parts[0]
                name_in = parts[1] if len(parts) == 2 else ""

                mol = None
                try:
                    mol = Chem.MolFromSmiles(smiles, sanitize=True)
                except Exception as e:
                    logger.warning(f".smi line {lineno}: sanitize=True parse error: {e}")

                if mol is None:
                    try:
                        mol = Chem.MolFromSmiles(smiles, sanitize=False)
                        if mol is not None:
                            try:
                                Chem.SanitizeMol(mol)
                            except Exception as e:
                                logger.warning(f".smi line {lineno}: manual sanitize failed: {e}")
                                mol = None
                    except Exception as e:
                        logger.warning(f".smi line {lineno}: sanitize=False parse error: {e}")
                        mol = None

                if mol is None:
                    logger.warning(f".smi line {lineno}: failed to parse SMILES '{smiles}'; skipping.")
                    continue

                base_name = _sanitize_name(name_in)
                final_name = _unique_name(base_name, seen_names)
                try:
                    mol.SetProp("_Name", final_name)
                except Exception:
                    pass

                names.append(final_name)
                mols.append(mol)

    except Exception as e:
        raise RuntimeError(f"Failed to read SMI file '{topology}': {e}")

    return names, mols


def get_element_counts(mol: Chem.Mol) -> Dict[str, int]:
    """
    Calculates the counts of each element in a molecular structure.

    This function iterates over all the atoms in the given molecular structure and
    counts the number of occurrences of each element, using the element's symbol
    as the key in the resulting dictionary.

    :param mol: The molecular structure from which to count element occurrences.
    :type mol: Chem.Mol
    :return: A dictionary where keys are element symbols (str) and values are the
        number of occurrences of each element (int).
    :rtype: Dict[str, int]
    """
    counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
    return counts


def count_total_atoms(mol: Chem.Mol) -> int:
    """
    Counts the total number of atoms present in a molecule.

    This function takes a molecule object as input and returns the
    total number of atoms present in the molecule. It provides a simple
    interface to retrieve the atom count using properties of the molecule
    object.

    :param mol: The molecule object for which the total number of atoms
        is to be counted. It must be an instance of `Chem.Mol`.
    :return: The total number of atoms in the given molecule.
    :rtype: int
    """
    return mol.GetNumAtoms()


def get_fragment_counts(mol: Chem.Mol) -> Dict[str, int]:
    """
    Analyzes the provided molecular structure and calculates the count of specific
    fragments using predefined fragment detection functions. Each function from
    the Fragments module starting with the prefix `fr_` is applied to the molecule
    to determine the count of the corresponding fragment. Fragment names are
    retrieved by stripping the prefix `fr_`.

    :param mol: RDKit molecule object representing the molecular structure to analyze
                 for fragment counts.
    :return: Dictionary where keys are fragment names (excluding `fr_` prefix) and
             values are their respective counts in the molecule.
    :rtype: Dict[str, int]
    """
    out: Dict[str, int] = {}
    for name in dir(Fragments):
        if name.startswith("fr_"):
            func = getattr(Fragments, name)
            try:
                cnt = int(func(mol))
            except Exception:
                continue
            if cnt > 0:
                out[name[3:]] = cnt
    return out


def get_molecular_volume(mol: Chem.Mol) -> float:
    """
    Calculates the molecular volume of a given molecule. If the molecule does not
    have any conformers, a conformer is generated using embedding and then optimized
    using the Universal Force Field (UFF). The molecular volume is computed thereafter.
    If an error occurs during the computation process, a fallback value of 0.0 is
    returned.

    :param mol: The molecule for which the molecular volume is to be calculated.
    :type mol: Chem.Mol
    :return: The calculated molecular volume or 0.0 as a fallback on failure.
    :rtype: float
    """
    try:
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
            AllChem.UFFOptimizeMolecule(mol)
        vol = AllChem.ComputeMolVolume(mol)
        return float(vol)
    except Exception:
        # Best-effort fallback: 0.0
        return 0.0


def get_molecular_weight(mol: Chem.Mol) -> float:
    """
    Calculate the exact molecular weight of a given molecular structure.

    This function uses the RDKit library to calculate the exact molecular
    weight (monoisotopic weight) of a molecule.

    :param mol: The molecular structure for which the exact molecular weight
        is computed. The parameter is expected to be an instance of
        `rdkit.Chem.Mol`.

    :return: The exact molecular weight of the molecule as a float.
    """
    return float(Descriptors.rdMolDescriptors.CalcExactMolWt(mol))


def get_num_rotatable_bonds(mol: Chem.Mol) -> int:
    """
    Calculates the number of rotatable bonds in a molecular structure.

    This function utilizes the `CalcNumRotatableBonds` method from RDKit's
    `rdMolDescriptors` module to determine the number of rotatable bonds in
    the given molecular structure. A rotatable bond is generally defined as
    a single bond between two heavy (non-hydrogen) atoms that allows for
    rotation around the bond.

    :param mol: The molecular structure for which the number of rotatable
        bonds is to be calculated, represented as an RDKit `Chem.Mol` object.
    :return: The total count of rotatable bonds in the given molecular
        structure as an integer.
    """
    return int(rdMolDescriptors.CalcNumRotatableBonds(mol))


def get_number_of_rings(mol: Chem.Mol) -> int:
    """
    Calculates and returns the number of rings present in the given molecule.

    This function uses RDKit's molecular descriptor calculation to determine
    the number of rings in a molecule.

    :param mol: The molecule object of type `Chem.Mol` for which the number
        of rings needs to be calculated.
    :return: An integer representing the number of rings in the molecule.
    """
    return int(rdMolDescriptors.CalcNumRings(mol))


def get_tpsa(mol: Chem.Mol) -> float:
    """
    Calculates the Topological Polar Surface Area (TPSA) of a molecule.

    TPSA is a descriptor used in cheminformatics to predict drug absorption and
    transport properties. It is calculated based on the contributions of polar
    atoms (e.g., oxygen, nitrogen) in a molecule, excluding hydrogens.

    :param mol: The molecule for which the TPSA is to be calculated.
    :type mol: Chem.Mol
    :return: The computed value of the molecule's Topological Polar Surface Area.
    :rtype: float
    """
    return float(rdMolDescriptors.CalcTPSA(mol))


def get_number_of_aromatic_rings(mol: Chem.Mol) -> int:
    """
    Calculate the number of aromatic rings in a molecular structure.

    This function uses RDKit's `CalcNumAromaticRings` method to determine the
    total number of aromatic rings in the given molecular structure.

    :param mol: The molecular structure represented as an RDKit `Chem.Mol` object.
    :return: The number of aromatic rings in the molecular structure.
    :rtype: int
    """
    return int(rdMolDescriptors.CalcNumAromaticRings(mol))


def generate_morgan_fingerprint(mol: Chem.Mol, radius: int = 2, nBits: int = 2048):
    """
    Generate a Morgan fingerprint for a given molecule using specific parameters.

    This function leverages the RDKit library to compute the Morgan fingerprint
    (a circular fingerprint representation) for a given molecule. The fingerprint
    is determined based on a given radius and a fixed number of bits, and it includes
    chirality by default. Morgan fingerprints are often used in cheminformatics for
    encoding molecular structures for tasks such as virtual screening or similarity
    search.

    :param mol: RDKit molecule object for which the fingerprint will be generated.
    :type mol: Chem.Mol
    :param radius: Radius to consider when generating the circular fingerprint. Default is 2.
    :type radius: int
    :param nBits: Size of the fingerprint in bits. Default is 2048.
    :type nBits: int
    :return: Morgan fingerprint of the molecule as a bit vector.
    :rtype: ExplicitBitVect
    """
    gen = FPG.GetMorganGenerator(radius=radius, fpSize=nBits, includeChirality=True)
    return gen.GetFingerprint(mol)


def generate_topological_torsion_fingerprint(mol: Chem.Mol, nBits: int = 2048):
    """
    Generates a topological torsion fingerprint for the given molecule.

    This function utilizes the RDKit fingerprint generator to compute
    the topological torsion fingerprint of a molecule. It creates a
    fingerprint of the molecule using the given number of bits.

    :param mol: The molecule for which to generate the fingerprint.
    :type mol: Chem.Mol
    :param nBits: The number of bits to use for the fingerprint. Defaults to 2048.
    :type nBits: int
    :return: The computed topological torsion fingerprint for the given molecule.
    :rtype: Any
    """
    gen = FPG.GetTopologicalTorsionGenerator(fpSize=nBits)
    return gen.GetFingerprint(mol)


def generate_atom_pair_fingerprint(mol: Chem.Mol):
    """
    Generates an atom pair fingerprint for a given molecule. The atom pair fingerprint
    is a type of molecular descriptor useful in cheminformatics for representing
    the structural features of molecules.

    :param mol: Molecule object for which the atom pair fingerprint will be generated.
    :type mol: Chem.Mol
    :return: Atom pair fingerprint of the input molecule.
    :rtype: RDKit fingerprint object
    """
    gen = FPG.GetAtomPairGenerator()
    return gen.GetFingerprint(mol)


def generate_maccs_keys(mol: Chem.Mol):
    """
    Generates the MACCS keys fingerprint for the given molecular object.

    This function computes the MACCS keys fingerprint for a specified
    molecule using the RDKit library. MACCS keys are a specific type of
    molecular fingerprint that encode structural information of a molecule
    into a binary vector representation. They are often employed in cheminformatics
    for molecular similarity and structure-based analyses.

    :param mol: A molecule object represented as an RDKit `Chem.Mol` instance. The molecule
        for which the MACCS keys fingerprint will be generated.
    :return: An RDKit `ExplicitBitVect` object representing the MACCS keys fingerprint
        of the input molecule.
    """
    return rdMolDescriptors.GetMACCSKeysFingerprint(mol)


def get_stereo_centers(mol: Chem.Mol) -> int:
    """
    Calculate the number of defined atom stereocenters in a molecule.

    This function utilizes RDKit's `CalcNumAtomStereoCenters` method to compute
    the count of stereocenters for the input molecular structure. It returns the
    number of stereocenters as an integer.

    :param mol: A molecule instance of type `Chem.Mol` for which the
        stereocenters are to be computed.
    :return: The number of defined atom stereocenters in the molecule.
    :rtype: int
    """
    # Number of defined atom stereocenters
    return int(rdMolDescriptors.CalcNumAtomStereoCenters(mol))


def get_net_charge(mol: Chem.Mol) -> int:
    """
    Get the net formal charge of a molecule.

    This function calculates the total formal charge of a molecule by
    summing the formal charges of all the atoms in the molecule.

    :param mol: A molecule object (Chem.Mol) from which the net formal charge
        will be calculated.
    :return: The net formal charge of the molecule as an integer.
    """
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def get_logp(mol: Chem.Mol) -> float:
    """
    Calculates and returns the logarithm of the partition coefficient (LogP) for the
    given molecular structure.

    LogP is a measure of the hydrophobicity (lipophilicity or hydrophilicity) of
    a molecule, calculated as the ratio of the concentrations of the compound in
    a mixture of two immiscible solvents, such as octanol and water. This value
    is commonly used in the field of chemistry and drug development.

    :param mol: Molecular structure for which the LogP value will be calculated.
    :type mol: Chem.Mol
    :return: The calculated LogP value as a floating-point number.
    :rtype: float
    """
    return float(Crippen.MolLogP(mol))


def fingerprint_to_string(fp) -> str:
    """
    Converts a fingerprint object to its string representation. The method first tries to use the RDKit
    library's `ToBitString` method if available, which is an efficient implementation for bitstrings.
    If not available, it falls back to converting the fingerprint to a list, mapping its elements to
    string representations, and joining them to form the result. In case of an error during the conversion,
    an empty string is returned.

    :param fp: The fingerprint object to be converted. The object should either have a `ToBitString`
        method for optimized string extraction or be iterable.
    :return: The string representation of the fingerprint. If the conversion fails, an empty string
        is returned.
    :rtype: str
    """
    # Prefer RDKit’s fast bitstring if available
    if hasattr(fp, "ToBitString"):
        return fp.ToBitString()
    try:
        return "".join(map(str, list(fp)))
    except Exception:
        return ""
