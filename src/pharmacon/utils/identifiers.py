"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Utility module to support generation and manipulation of timestamps, UUIDs, hashes,
and secure tokens.
"""

from __future__ import annotations

import json
import hmac
import base64
import hashlib
import numpy as np
from typing import TYPE_CHECKING, Any, Dict, Union
from collections import Counter

if TYPE_CHECKING:
    import MDAnalysis




__all__ = [
    "generate_mda_blueprint",
    "create_mda_artifact_token",
    "validate_mda_artifact_token",
]


def generate_mda_blueprint(*,
                           u: MDAnalysis.Universe,
                           id_bits: int = 128,
                           encoding: str = "base64url",
                           round_nd: int = 4,
                           include_debug: bool = False,
                           relax_resids: bool = False,
                           **kwargs: Any,) -> str | Dict:
    """
    Generate a unique blueprint ID or detailed payload based on the topology and
    other properties of the MDAnalysis Universe. This is achieved by hashing
    normalized input data using a controlled and deterministic methodology.

    :param u: The MDAnalysis Universe containing the molecular topology and
        trajectory data.
    :type u: MDAnalysis.Universe
    :param id_bits: The number of bits to encode in the resulting ID. Must be
        a multiple of 8 within the range [8, 256]. Default is 128.
    :param encoding: The encoding format for the generated ID. Available
        options are 'base64url', 'base32', and 'hex'. Default is 'base64url'.
    :param round_nd: The number of decimal places for rounding attributes
        such as masses and charges. Default is 4.
    :param include_debug: If True, return a detailed payload with additional
        debug information. If False, return only the compact ID. Default is False.
    :param relax_resids: If True, relax the residue ID consideration in the
        topology block to improve tolerance to certain differences in atom groupings.
        Default is False.
    :param kwargs: Additional key-value arguments to control and modify the
        hashing behavior, such as the `__hash_keys__` and `__runtime__` policies.

    :return: A unique ID string or detailed payload dictionary, depending on the
        value of `include_debug`.
    :rtype: str | Dict
    """

    def _safe_len(x) -> int:
        """
        Safely calculates the length of the input object. If the object doesn't
        support the length operation, it will return 0 instead of raising an
        exception.

        :param x: The object whose length is to be determined. Can be any type.
        :type x: Any
        :return: The length of the object if it supports the `len()` operation,
            otherwise returns 0.
        :rtype: int
        """
        try:
            return len(x)
        except Exception:
            return 0

    def _sorted_pairs(pairs):
        """
        Sorts a list of pairs, ensuring each pair is ordered and the list of pairs
        is sorted in ascending order.

        This function takes a list of tuples (pairs), orders elements within each
        pair such that the smaller element comes first, and then sorts the list
        of pairs in ascending order.

        :param pairs: A list of tuples where each tuple contains two comparable
                      elements to be sorted within and as a list.
        :return: A list of tuples where each tuple is sorted in ascending order,
                 and the entire list is also sorted.
        :rtype: list
        """
        pairs = [(min(i, j), max(i, j)) for (i, j) in pairs]
        pairs.sort()
        return pairs

    def _round_list(arr, nd):
        """
        Rounds elements of a list or array to a specified number of decimal places.

        This function accepts a list or array-like input and attempts to round each of
        its elements to the specified number of decimal places. If the input array is
        empty, an empty list is returned. In cases where the input cannot be entirely
        resolved to a numerical format, it handles the input on an element-wise basis,
        rounding each item individually or replacing non-numeric values with None.

        :param arr: The input list or array-like object containing numeric values. May
            include None or non-numeric entries.
        :type arr: list or array-like
        :param nd: The number of decimal places to which elements should be rounded.
        :type nd: int
        :return: A list with the rounded numeric elements, or an empty list if the
            input array is empty, or None if rounding could not be performed.
        :rtype: list or None
        """
        try:
            a = np.asarray(arr, dtype=float)
            if a.size == 0:
                return []
            a = np.round(a, nd)
            return a.tolist()
        except Exception:
            try:
                return [None if v is None else round(float(v), nd) for v in arr]
            except Exception:
                return None

    def _canon_jsonable(v):
        """
        Converts a given value into a JSON-serializable canonical form. This function ensures
        the order of dictionary keys, sorts elements for sets, frozensets, and recursively
        formats nested structures to maintain consistency. If the input contains types such
        as numpy scalars, attempts to convert them using the `item` method. Fails gracefully
        if the conversion of a scalar is not possible.

        :param v: Input value to be converted into a JSON-serializable canonical form. It can
            be of any type including dictionaries, lists, tuples, sets, frozensets, or numpy
            scalar-like objects.
        :return: A JSON-serializable canonical form of the input value with sorted structures
            where applicable and recursively processed nested components.
        """
        if isinstance(v, dict):
            return {str(k): _canon_jsonable(v[k]) for k in sorted(v, key=lambda x: str(x))}
        if isinstance(v, (list, tuple)):
            return [_canon_jsonable(x) for x in v]
        if isinstance(v, (set, frozenset)):
            return [_canon_jsonable(x) for x in sorted(v, key=lambda x: str(x))]
        # numpy scalars
        try:
            if hasattr(v, "item") and callable(v.item):
                return v.item()
        except Exception:
            pass
        return v

    # Dynamic kwargs policy
    # Optional reserved controls:
    #   __hash_keys__ : iterable[str]  -> only these kwargs affect the hash
    #   __runtime__   : dict          -> never affects the hash (debug/UI/etc)
    #
    # If __hash_keys__ is absent, ALL non-reserved kwargs affect the hash.
    hash_keys_raw = kwargs.pop("__hash_keys__", None)
    runtime_kwargs = kwargs.pop("__runtime__", {}) or {}

    if hash_keys_raw is None:
        hash_keys = None  # means "hash everything (except reserved)"
    else:
        try:
            hash_keys = {str(k).strip() for k in hash_keys_raw}
        except Exception:
            raise TypeError("__hash_keys__ must be an iterable of strings")

    # Normalize + canonicalize kwargs that participate in hashing
    hashed_kwargs: Dict[str, Any] = {}
    ignored_kwargs: Dict[str, Any] = {}  # kwargs not in hash_keys (when provided)

    for k, v in kwargs.items():
        kk = str(k).strip()
        if not kk:
            continue

        if hash_keys is None or kk in hash_keys:
            hashed_kwargs[kk] = _canon_jsonable(v)
        else:
            ignored_kwargs[kk] = _canon_jsonable(v)

    # Frames
    try:
        total_frames_int = int(getattr(u.trajectory, "n_frames", 0) or 0)
    except Exception:
        try:
            total_frames_int = int(len(u.trajectory))
        except Exception:
            total_frames_int = 0

    # Topology-derived bits
    n_atoms = _safe_len(u.atoms)
    n_residues = _safe_len(u.residues)
    n_segments = _safe_len(u.segments)
    n_bonds = _safe_len(getattr(u, "bonds", []))

    resnames = [r.resname for r in u.residues] if n_residues else []
    elements = [getattr(a, "element", None) for a in u.atoms] if n_atoms else []

    # Make composition dicts order-stable
    resname_counts = dict(sorted(Counter(resnames).items()))
    element_counts = dict(sorted(Counter([e for e in elements if e]).items()))

    # ordered per-atom tuples (stable with MDAnalysis atom order)
    atom_tuples = []
    for a in u.atoms:
        segid = getattr(a, "segid", "") or ""
        resname = getattr(a, "resname", "") or ""
        resid = getattr(a, "resid", None)
        aname = getattr(a, "name", "") or ""
        elem = getattr(a, "element", None) or ""

        if relax_resids:
            atom_tuples.append((str(segid), str(resname), str(aname), str(elem)))
        else:
            atom_tuples.append(
                (
                    str(segid),
                    str(resname),
                    int(resid) if resid is not None else None,
                    str(aname),
                    str(elem),
                )
            )

    try:
        bond_pairs = _sorted_pairs(
            [(int(b.atoms[0].index), int(b.atoms[1].index)) for b in u.bonds]
        )
    except Exception:
        bond_pairs = []

    try:
        masses = _round_list(u.atoms.masses, round_nd)
    except Exception:
        masses = None

    try:
        charges = _round_list(u.atoms.charges, round_nd)
    except Exception:
        charges = None

    try:
        n_fragments = len(u.atoms.fragments())
    except Exception:
        n_fragments = None

    # Canonical payload (hash-participating only)
    payload = {
        "_schema": "bp_exact_v1",
        "analysis_kwargs": dict(sorted(hashed_kwargs.items(), key=lambda kv: kv[0])),
        "total_frames": total_frames_int,
        "counts": {
            "atoms": n_atoms,
            "residues": n_residues,
            "segments": n_segments,
            "bonds": n_bonds,
        },
        "composition": {
            "resname_counts": resname_counts,
            "element_counts": element_counts,
        },
        "topology_block": {
            "atom_tuples": atom_tuples,
            "bond_pairs": bond_pairs,
            "masses": masses,
            "charges": charges,
            "n_atoms": n_atoms,
            "n_residues": n_residues,
            "n_fragments": n_fragments,
            "n_segments": n_segments,
            "relax_resids": bool(relax_resids),
            "round_nd": int(round_nd),
        },
    }

    # hash -> compact ID
    payload_bytes = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    digest = hashlib.sha256(payload_bytes).digest()

    if id_bits % 8 != 0 or not (8 <= id_bits <= 256):
        raise ValueError("id_bits must be a multiple of 8 within [8, 256].")

    short = digest[: id_bits // 8]

    enc = (encoding or "").lower()
    if enc in ("base64url", "b64url", "b64"):
        bid = base64.urlsafe_b64encode(short).decode("ascii").rstrip("=")
    elif enc in ("hex", "hexlower"):
        bid = short.hex()
    elif enc in ("base32", "b32"):
        bid = base64.b32encode(short).decode("ascii").rstrip("=")
    else:
        raise ValueError("encoding must be one of: 'base64url', 'hex', 'base32'")

    if not include_debug:
        return bid
    else:
        return {
            "id": bid,
            "full_sha256_hex": digest.hex(),
            "counts": payload["counts"],
            "topology_sizes": {
                "atom_tuples": len(atom_tuples),
                "bond_pairs": len(bond_pairs),
            },
            "hash_policy": {
                "hash_all_kwargs": (hash_keys is None),
                "hash_keys": None if hash_keys is None else sorted(hash_keys),
            },
            "analysis_kwargs_hashed": hashed_kwargs,
            "analysis_kwargs_ignored": ignored_kwargs,
            "runtime_kwargs": _canon_jsonable(runtime_kwargs),
        }


def create_mda_artifact_token(*,
                              blueprint: Union[str, Dict],
                              secret: str,
                              namespace: str = "mda_artifact",
                              token_bits: int = 128,
                              encoding: str = "base64url",
                              include_debug: bool = False):
    """
    Generate an artifact token with the provided blueprint and secret using HMAC
    and specified encoding. This utility is designed for creating securely encoded
    tokens that include a reference to a namespace and blueprint identifier.
    Additionally, it supports optional debugging information.

    :param blueprint: A blueprint identifier, which can be a string or dictionary
                      containing the "id" key. The value will be used to generate
                      the token.
    :type blueprint: Union[str, Dict]
    :param secret: The secret key used to calculate the HMAC digest for token
                   generation.
    :param namespace: An optional namespace for the token, defaulting to
                      "mda_artifact".
    :type namespace: str
    :param token_bits: The bit-length for the generated token, must be a multiple
                       of 8 and within the range [8, 256]. Defaults to 128.
    :type token_bits: int
    :param encoding: Specifies the output encoding format for the token. Permitted
                     values are "base64url", "hex", or "base32". Defaults to
                     "base64url".
    :type encoding: str
    :param include_debug: Optional flag to include debugging information such as
                          blueprint ID, namespace, and token settings in the
                          output. Defaults to False.
    :type include_debug: bool
    :return: If `include_debug` is False, returns the encoded artifact token as a
             string. If True, returns a dictionary containing the artifact token
             and related debug details.
    :rtype: Union[str, Dict]
    """

    # extract blueprint_id
    if isinstance(blueprint, dict):
        blueprint_id = blueprint.get("id")
        if not blueprint_id:
            raise ValueError("Blueprint dict must contain 'id'")
    elif isinstance(blueprint, str):
        blueprint_id = blueprint
    else:
        raise TypeError("blueprint must be str or dict")

    # validate bits
    if token_bits % 8 != 0 or not (8 <= token_bits <= 256):
        raise ValueError("token_bits must be multiple of 8 within [8, 256]")

    # HMAC
    message = f"{namespace}:{blueprint_id}".encode()
    digest = hmac.new(secret.encode(), message, hashlib.sha256).digest()
    short = digest[: token_bits // 8]

    # encode
    enc = encoding.lower()
    if enc in ("base64url", "b64", "b64url"):
        token = base64.urlsafe_b64encode(short).decode().rstrip("=")
    elif enc in ("hex",):
        token = short.hex()
    elif enc in ("base32", "b32"):
        token = base64.b32encode(short).decode().rstrip("=")
    else:
        raise ValueError("encoding must be base64url, hex, or base32")

    if not include_debug:
        return token

    return {
        "artifact_token": token,
        "blueprint_id": blueprint_id,
        "namespace": namespace,
        "token_bits": token_bits,
    }



def validate_mda_artifact_token(*,
                                artifact_token: str,
                                blueprint: Union[str, Dict[str, Any]],
                                secret: str,
                                namespace: str = "mda_artifact",
                                token_bits: int = 128,
                                encoding: str = "base64url",
                                return_debug: bool = False):
    """
    Validate an MDA artifact token by comparing it with an expected token that is
    generated using the provided blueprint, secret, and other parameters.

    :param artifact_token: The artifact token to validate.
    :param blueprint: The blueprint, either as a string or a dictionary, used for
        token generation.
    :param secret: The secret key used for HMAC generation of the token.
    :param namespace: Optional namespace to scope the token generation. Defaults
        to "mda_artifact".
    :param token_bits: Optional length of the token in bits. Defaults to 128.
    :param encoding: Optional encoding format used for the token. Defaults to
        "base64url".
    :param return_debug: Optional flag to indicate whether to return debugging
        information. Defaults to False.
    :return: A boolean value indicating whether the validation was successful if
        `return_debug` is False. If `return_debug` is True, a dictionary with
        validation details is returned.
    """
    expected = create_mda_artifact_token(
        blueprint=blueprint,
        secret=secret,
        namespace=namespace,
        token_bits=token_bits,
        encoding=encoding,
        include_debug=False,
    )

    valid = hmac.compare_digest(artifact_token, expected)

    if not return_debug:
        return valid

    if isinstance(blueprint, dict):
        blueprint_id = blueprint.get("id")
    else:
        blueprint_id = blueprint

    return {
        "valid": valid,
        "provided_token": artifact_token,
        "expected_token": expected,
        "blueprint_id": blueprint_id,
        "namespace": namespace,
        "token_bits": token_bits,
        "encoding": encoding,
    }
