"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Module :mod:`pharmacon.utils.fingerprint`.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Final

from pharmacon.constants import __version__




__all__ = [
    "PharmaconFileSignature",
    "create_pharmacon_signature",
]


_PHARMACON_MAGIC: Final[str] = "PHARMACON"
_DEFAULT_SCHEMA_VERSION: Final[str] = __version__
_INTERNAL_NAMESPACE_SALT: Final[str] = "PHARMACON-FILE-SIGNATURE"


@dataclass(frozen=True, slots=True)
class PharmaconFileSignature:
    """
    Represents the file signature details for a Pharmacon format file.

    This class encapsulates the attributes used to define and identify a Pharmacon
    file signature. It includes information about the file's magic identifier,
    format nomenclature, version, and specific commands along with its unique
    signature and fingerprint. The class is immutable and employs slots to
    optimize memory usage.

    :ivar magic: Magic identifier of the Pharmacon file.
    :type magic: str
    :ivar format_name: Name of the file's format.
    :type format_name: str
    :ivar version: Version of the Pharmacon file format.
    :type version: str
    :ivar command: Command category associated with the file.
    :type command: str
    :ivar subcommand: Subcommand specified within the file.
    :type subcommand: str
    :ivar signature: Unique signature of the file.
    :type signature: str
    :ivar fingerprint: Unique fingerprint of the Pharmacon file.
    :type fingerprint: str
    """

    magic: str
    format_name: str
    version: str
    command: str
    subcommand: str
    signature: str
    fingerprint: str


def _normalize_token(value: str, field_name: str) -> str:
    """
    Normalizes and validates a string token by performing several checks and transformations.
    The function trims any leading or trailing whitespace, ensures the string is neither
    empty nor contains invalid control characters, and converts the string to uppercase.

    :param value: The input string to be normalized.
    :type value: str
    :param field_name: A descriptive name identifying the purpose of the input string,
                       used in error messages for better readability.
    :type field_name: str
    :return: The normalized and uppercase version of the input string.
    :rtype: str
    :raises TypeError: If the `value` parameter is not of type `str`.
    :raises ValueError: If the `value` parameter is an empty string or contains invalid
                        control characters like NULL, carriage return, or newline.
    """
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string.")

    normalized: str = value.strip()
    if not normalized:
        raise ValueError(f"{field_name} cannot be empty.")

    if any(char in normalized for char in ("\x00", "\r", "\n")):
        raise ValueError(f"{field_name} contains invalid control characters.")

    return normalized.upper()


def _chunk_text(value: str, chunk_size: int = 4, separator: str = "-") -> str:
    """
    Splits the given text into chunks of a specified size, separated by a
    defined separator. If the chunk size is set to a non-positive number,
    a `ValueError` will be raised.

    :param value:
        The string that needs to be divided into chunks.
    :param chunk_size:
        The length of each chunk. Default is 4. Must be greater than 0.
    :param separator:
        The string used to join the chunks. Default is "-".
    :return:
        A string composed of chunks of the original text separated by
        the defined separator.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")

    return separator.join(
        value[index:index + chunk_size]
        for index in range(0, len(value), chunk_size)
    )


def create_pharmacon_signature(*,
                               format_name: str,
                               command: str,
                               subcommand: str,
                               version: str = _DEFAULT_SCHEMA_VERSION
                               ) -> PharmaconFileSignature:
    """
    Creates a Pharmacon-specific file signature by combining and encoding necessary
    information such as format name, command, subcommand, and version. This function
    normalizes the inputs, generates a payload, calculates hash digests, and constructs
    the signature and fingerprint for the file.

    :param format_name: The name identifying the specific format of the Pharmacon file.
    :type format_name: str
    :param command: Primary command associated with the Pharmacon file.
    :type command: str
    :param subcommand: Subcommand related to the Pharmacon file.
    :type subcommand: str
    :param version: Version of the file schema. Defaults to _DEFAULT_SCHEMA_VERSION
                     if not provided. Cannot be an empty string.
    :type version: str
    :return: A PharmaconFileSignature instance containing all necessary signature
             details for the file.
    :rtype: PharmaconFileSignature
    :raises ValueError: If the provided version is empty.
    """
    normalized_format: str = _normalize_token(format_name, "format_name")
    normalized_command: str = _normalize_token(command, "command")
    normalized_subcommand: str = _normalize_token(subcommand, "subcommand")
    normalized_version: str = version.strip()

    if not normalized_version:
        raise ValueError("version cannot be empty.")

    base_payload: str = "::".join(
        (
            _PHARMACON_MAGIC,
            normalized_format,
            normalized_version,
            normalized_command,
            normalized_subcommand,
            _INTERNAL_NAMESPACE_SALT,
        )
    )

    digest: str = hashlib.sha256(base_payload.encode("utf-8")).hexdigest().upper()
    signature_body: str = _chunk_text(digest[:32], chunk_size=4)
    fingerprint: str = hashlib.blake2b(
        base_payload.encode("utf-8"),
        digest_size=16,
    ).hexdigest().upper()

    signature: str = (
        f"{_PHARMACON_MAGIC}::"
        f"{normalized_command}::"
        f"{normalized_subcommand}::"
        f"{signature_body}"
    )

    return PharmaconFileSignature(
        magic=_PHARMACON_MAGIC,
        format_name=normalized_format,
        version=normalized_version,
        command=normalized_command,
        subcommand=normalized_subcommand,
        signature=signature,
        fingerprint=fingerprint,
    )
