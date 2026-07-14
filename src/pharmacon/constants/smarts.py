"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Constants and regex patterns for biomolecular data handling.

This module provides constants and regular expressions useful for parsing and
working with biomolecular data, including amino acids, nucleotides, water
residues, and ion representations. Additionally, it contains patterns for
identifying backbone atom names relevant in structural biology.

Constants defined in this module include amino acid three-letter codes, nucleic
acid residue names, common water residue names, ion-to-element mappings, and
specific backbone and sidechain atom type groupings. These are critical for
analyzing molecular structures or performing computational simulations.

Attributes:
    AA3 (Final[Tuple[str, ...]]): Tuple of three-letter amino acid codes
        including standard and modified residues.
    NA_RES (Final[Tuple[str, ...]]): Tuple of nucleic acid residue names,
        including DNA and RNA nucleotides in various formats.
    WATER_RESIDUES (Final[Tuple[str, ...]]): Tuple of common water residue
        identifiers in molecular structures.
    ION_RES_TO_ELT (Final[Dict[str, str]]): Dictionary mapping common ion
        residue names to their corresponding chemical elements.
    PROTEIN_CARBONS (Final[Tuple[str, ...]]): Tuple of typical side-chain
        carbon atom names found in amino acids.
    SIDECHAIN_OXY (Final[Tuple[str, ...]]): Tuple of sidechain oxygen atom
        names present in amino acids.
    SIDECHAIN_NITRO (Final[Tuple[str, ...]]): Tuple of sidechain nitrogen
        atom names present in amino acids.
    NUC_OXY (Final[Tuple[str, ...]]): Tuple of nucleotide backbone oxygen
        atom names.
    TWO_LETTER (Final[Tuple[str, ...]]): Tuple of two-letter element
        abbreviations for commonly occurring elements in biomolecules.
    BB_RE (Final[re.Pattern[str]]): Regular expression pattern for backbone
        atom name matching, designed to identify backbone atoms (such as N,
        CA, C, O) and various forms of hydrogen atoms.
"""


import re
from typing import Final, Dict, Tuple




__all__ = [
    "AA3",
    "AA3_to_AA1",
    "NA_RES",
    "WATER_RESIDUES",
    "WATER_RESIDUES_MDA_SELECTION_STR",
    "ION_RES_TO_ELT",
    "PROTEIN_CARBONS",
    "SIDECHAIN_OXY",
    "SIDECHAIN_NITRO",
    "NUC_OXY",
    "TWO_LETTER",
    "BB_RE",
    "HYDROPHOBIC_PATTERNS",
    "HYDROPHOBIC_ELEMENTS",
    "HYDROGEN_BOND_DONOR_PATTERNS",
    "HYDROGEN_BOND_ACCEPTOR_PATTERNS",
    "HYDROGEN_BOND_DONOR_ELEMENTS",
    "HYDROGEN_BOND_ACCEPTOR_ELEMENTS",
    "NEGATIVELY_CHARGED_PATTERNS",
    "POSITIVELY_CHARGED_PATTERNS",
    "METAL_ELEMENTS",
    "HALOGEN_PATTERNS",
    "HALOGEN_ELEMENTS",
    "AROMATIC_PATTERNS",
]


AA3: Final[Tuple[str, ...]] = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HSD", "HSE", "HSP",
    "HIP", "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL", "ASH", "GLH", "CYM", "CYX", "MSE", "SEC", "PYL",
)

AA3_to_AA1: Final[Dict[str, str]] = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E",
    "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I",
    "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H",
    "HIE": "H", "HIP": "H", "ASH": "D", "GLH": "E",
    "CYX": "C", "CYM": "C", "MSE": "M",
}

NA_RES: Final[Tuple[str, ...]] = (
    "A", "C", "G", "U", "DA", "DC", "DG", "DT", "ADE", "CYT", "GUA", "URA", "THY",
    "RGU", "RCY", "RAN")

WATER_RESIDUES: Final[Tuple[str, ...]] = (
    "HOH", "WAT", "TIP3", "TIP4", "TIP5", "SPC", "SPE", "T3P", "T4P", "T5P",
)

WATER_RESIDUES_MDA_SELECTION_STR = "resname HOH or resname WAT or resname TIP3 or resname TIP4 or resname TIP5 or resname SPC or resname SPE or resname T3P or resname T4P or resname T5P"

ION_RES_TO_ELT: Final[Dict[str, str]] = {
    "CL": "Cl", "CL-": "Cl", "CLA": "Cl", "CL1": "Cl",
    "BR": "Br", "BR-": "Br", "IOD": "I", "I": "I", "I-": "I",
    "F": "F", "F-": "F", "NA": "Na", "NA+": "Na", "SOD": "Na",
    "K": "K", "K+": "K", "POT": "K", "LI": "Li", "LI+": "Li",
    "MG": "Mg", "MG2": "Mg", "MG2+": "Mg", "CA": "Ca", "CA2": "Ca",
    "CA2+": "Ca", "CAL": "Ca", "SR": "Sr", "BA": "Ba", "RB": "Rb",
    "CS": "Cs", "ZN": "Zn", "FE": "Fe", "FE2": "Fe", "MN": "Mn",
    "CO": "Co", "CU": "Cu", "NI": "Ni", "AG": "Ag", "CD": "Cd",
    "HG": "Hg", "PB": "Pb", "SN": "Sn", "PT": "Pt", "PD": "Pd",
    "IR": "Ir", "RU": "Ru", "RH": "Rh", "OS": "Os", "AU": "Au",
}

PROTEIN_CARBONS: Final[Tuple[str, ...]] = ("CA", "CB", "CG", "CD", "CE", "CZ")
SIDECHAIN_OXY: Final[Tuple[str, ...]] = ("OG", "OG1", "OH", "OE", "OE1", "OE2", "OD", "OD1", "OD2")
SIDECHAIN_NITRO: Final[Tuple[str, ...]] = ("NE", "NE1", "NE2", "ND", "ND1", "ND2", "NZ")
NUC_OXY: Final[Tuple[str, ...]] = ("OP1", "OP2", "O5", "O3", "O4", "O2")

TWO_LETTER: Final[Tuple[str, ...]] = (
    "CL", "BR", "NA", "MG", "CA", "FE", "ZN", "MN", "CO", "CU", "NI", "HG", "CD", "AG", "AU", "PT", "PD", "TI", "CR",
    "SN", "SB", "PB", "SR", "BA", "CS", "RB", "LI", "AL", "SI", "SE", "AS", "GA", "GE", "IR", "OS", "RH", "RU", "XE",
    "KR", "NE", "AR", "MO", "ZR", "NB", "TA", "RE", "HF"
)


# Backbone atom name matcher
BB_RE: Final[re.Pattern[str]] = re.compile(
    r"""^(?:
     N|CA|C|O|
     OXT|OT[12]|O[12]|
     HA(?:[123])?|
     HN|H(?:[123])?|
     HT(?:[123])?
     )$""",
    re.X | re.I,
)

HYDROPHOBIC_PATTERNS: Final[Tuple[str, ...]] = (
    "[c;!$([c]-[O,N,S,P])]",
    "[C;!$([C]-[O,N,S,P])]",
    "[S;X2,X3;!$([S](=O))]",
    "[P;!$([P](=O))]",
    "[Cl,Br,I]",
    "[H][c;!$([c]-[O,N,S,P])]",
    "[H][C;!$([C]-[O,N,S,P])]",
)

HYDROPHOBIC_ELEMENTS: Final[Tuple[str, ...]] = ("C", "H", "S", "BR", "I", "F", "CL")

HYDROGEN_BOND_DONOR_PATTERNS: Final[Tuple[str, ...]] = ("[#1][O;X2]",
                                                        "[#1]S[#6]",
                                                        "[#1][C;X2]#[C;X2]",
                                                        "[#1][NX3]C(=[NX2])[#6]",
                                                        "[#1][#7]"
                                                        )

HYDROGEN_BOND_ACCEPTOR_PATTERNS: Final[Tuple[str, ...]] = ("[N;X1]#[#6]",
                                                           "[N;X1]#CC",
                                                           "[N;X2](=C~[C,c])C",
                                                           "[N;X2](O)=N[a]",
                                                           "[N;X2](=N-O)[a]",
                                                           "[n;X2]1ccccc1",
                                                           "[n;X2]([a])([a])",
                                                           "[N;X2](=C~[C,c])(~[*])",
                                                           "[N;X3](C)(C)[N;X3]C",
                                                           "[N;X2](=C)(~[*])",
                                                           "[N;X2](~[C,c])=[N;X2]",
                                                           "[n;X2]1c[nH]cc1",
                                                           "O=[S;X4](=O)([!#8])([!#8])",
                                                           "[O;X2]C",
                                                           "[O;X2]N",
                                                           "[O;X1]=[C,c]",
                                                           "o","[O;X2](C)C",
                                                           "[O;X2]c1ncccc1",
                                                           "[O;X2]~[a]",
                                                           "O=PO([!#1])",
                                                           "[O;X2]",
                                                           "[S;X2](C)C",
                                                           "[S;X2](=C)N",
                                                           "O=C[O-]"
                                                           )

HYDROGEN_BOND_DONOR_ELEMENTS: Final[Tuple[str, ...]] = ("O", "S", "N")

HYDROGEN_BOND_ACCEPTOR_ELEMENTS: Final[Tuple[str, ...]] = ("O", "N", "S")

NEGATIVELY_CHARGED_PATTERNS: Final[Tuple[str, ...]] = ("[-]",)

POSITIVELY_CHARGED_PATTERNS: Final[Tuple[str, ...]] = ("[+]",)

METAL_ELEMENTS: Final[Tuple[str, ...]] = ("LI", "NA", "K", "MG", "CA", "ZN", "FE", "MN", "CO", "CU", "NI", "GA")

HALOGEN_PATTERNS: Final[Tuple[str, ...]] = ("[#9]", "[#17]", "[#35]")

HALOGEN_ELEMENTS: Final[Tuple[str, ...]] = ("BR", "CL", "I", "F")

AROMATIC_PATTERNS: Final[Tuple[str, ...]] = ("[a;r4]1aaa1",
                                             "[a;r5]1aaaa1",
                                             "[a;r6]1aaaaa1",
                                             "[a;r7]1aaaaaa1",
                                             "[a;r8]1aaaaaaa1")
