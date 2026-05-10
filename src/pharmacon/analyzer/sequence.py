"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Provides functionality to extract protein sequence information from an MDAnalysis
Universe object and convert sequence data into FASTA formatted strings.

This module defines two main functions. The `get_sequence` function is used to
extract amino acid sequences, residue identifiers, and topology information for
protein chains in the input Universe. The `sequence_dict_to_fasta` function is used
to convert the extracted sequence data into a standardized FASTA format for output.

Functions:
- get_sequence: Extracts chain-wise sequence data from an MDAnalysis Universe.
- sequence_dict_to_fasta: Converts sequence data into FASTA format.
"""


import warnings
import MDAnalysis
import numpy as np
from typing import Dict, Iterable, Optional


from pharmacon.constants import AA3_to_AA1




__all__ = [
    "get_sequence",
    "sequence_dict_to_fasta",
]


warnings.filterwarnings("ignore")


def get_sequence(u: MDAnalysis.Universe) -> Dict[str, Dict]:
    """
    Extracts protein sequence information from an MDAnalysis Universe object. This function
    retrieves amino acid sequences, residue identifiers, and topology information based on CA
    atoms for chains found within the Universe. If no chain identifiers exist, the Universe is
    treated as a single chain.

    dict[chainid -> {
            'aa1_seq'         : str,
            'aa3_list'        : list[str],
            'resids_topology' : np.ndarray,
            'resid_seq'       : np.ndarray,
        }]

    :param u: The MDAnalysis Universe object containing the molecular data.
    :type u: MDAnalysis.Universe
    :return: A dictionary with chain IDs as keys. Each key maps to another dictionary, containing
             the fields:
             - "aa1_seq" (str): A single-letter amino acid sequence.
             - "aa3_list" (List[str]): A list of three-letter amino acid codes.
             - "resids_topology" (np.ndarray): Residue identifiers in topology order.
             - "resid_seq" (np.ndarray): Residue sequence numbers starting from 1.
    :rtype: Dict[str, Dict]
    """

    def _extract_chain(chainid: Optional[str]):
        """
        Extracts chain information from a provided chain identifier.

        This function selects alpha-carbon atoms (CA) from a molecular system, applies
        a chain-specific filter if a chain ID is provided, and processes the data to
        construct amino acid sequences and residue identifiers. If the chain ID is not
        specified, all CA atoms are considered.

        :param chainid: A string representing the chain ID to filter the alpha-carbon
                        atoms in the molecular system, or None to select all chains.
        :type chainid: Optional[str]
        :return: A dictionary containing processed chain information, including:
                 - "aa1_seq": A string of the single-letter amino acid sequence.
                 - "aa3_list": A list of three-letter amino acid identifiers.
                 - "resids_topology": A NumPy array of residue IDs sorted by topology.
                 - "resid_seq": A NumPy array of sequential residue indices.
        :rtype: Optional[dict]
        """
        if chainid is None:
            ca = u.select_atoms("protein and name CA")
        else:
            ca = u.select_atoms(f"protein and name CA and chainid {chainid}")

        if ca.n_atoms == 0:
            return None

        # Sort by topology resid
        try:
            order = ca.resids.argsort(kind="stable")
            ca = ca[order]
        except Exception:
            pass

        aa3_list = []
        aa1_list = []
        resids_top = []

        for atom in ca:
            aa3 = (atom.resname or "").strip().upper() or "UNK"
            aa1 = AA3_to_AA1.get(aa3, "X")
            aa3_list.append(aa3)
            aa1_list.append(aa1)
            resids_top.append(int(atom.resid) if atom.resid is not None else -1)

        return {
            "aa1_seq": "".join(aa1_list),
            "aa3_list": aa3_list,
            "resids_topology": np.asarray(resids_top, dtype=int),
            "resid_seq": np.arange(1, len(aa1_list) + 1, dtype=int),
        }

    results: Dict[str, Dict] = {}

    # Collect unique chain IDs (ignore empty / None)
    chainids = sorted(
        {cid for cid in u.atoms.chainIDs if cid and cid.strip()}
    )

    if chainids:
        for cid in chainids:
            data = _extract_chain(cid)
            if data is not None:
                results[cid] = data

        if not results:
            raise RuntimeError("No protein CA atoms found in any chain.")
        return results

    # No chain IDs at all → treat as single chain
    data = _extract_chain(None)
    if data is None:
        raise RuntimeError("No protein CA atoms found in the Universe.")
    return {"NO_CHAINID": data}


def sequence_dict_to_fasta(sequences: Dict[str, Dict], *, wrap: int = 80, header_prefix: str = "chain") -> str:
    """
    Converts a dictionary of sequence data into a FASTA formatted string. Each sequence
    is represented by a header and a wrapped sequence. The header can include an optional
    residue ID range if available in the input data.

    :param sequences: Dictionary where each key is a chain ID (string) and the value is
        another dictionary containing sequence data for the given chain ID.
    :param wrap: Maximum number of characters per line in the FASTA sequence. Defaults
        to 80. If set to 0 or a negative value, the sequence will not be wrapped.
    :type wrap: int
    :param header_prefix: Prefix to include in the FASTA headers to identify chains.
        Defaults to "chain".
    :type header_prefix: str
    :return: A string in FASTA format, containing the sequences and their respective headers.
    :rtype: str
    """

    def _wrap(seq: str, width: int) -> Iterable[str]:
        """
        Wraps a string sequence into a specified width, dividing the string into
        substrings of the given width. If the width is less than or equal to zero,
        the entire sequence is returned in a single iteration.

        :param seq: A string sequence to be wrapped.
        :param width: An integer specifying the width to wrap the sequence into.
        :return: An iterable of strings, each with a maximum length equal to the
            specified width, unless the width is less than or equal to zero.
        """
        if width <= 0:
            yield seq
        else:
            for i in range(0, len(seq), width):
                yield seq[i : i + width]

    fasta_lines: list[str] = []

    for chainid in sorted(sequences.keys()):
        data = sequences[chainid]

        aa1_seq = data.get("aa1_seq")
        if not aa1_seq:
            continue

        resids = data.get("resids_topology")
        resid_range = ""
        if resids is not None and len(resids) > 0:
            resid_range = f"{resids[0]}-{resids[-1]}"

        header = f">{header_prefix}:{chainid}"
        if resid_range:
            header += f" resid={resid_range}"

        fasta_lines.append(header)

        for line in _wrap(aa1_seq, wrap):
            fasta_lines.append(line)

    return "\n".join(fasta_lines)
