"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Read Pharmacon INI files into nested :class:`~types.SimpleNamespace` objects.

Pharmacon uses :mod:`configobj` as the INI backend because it natively supports
arbitrarily nested sections via bracket depth::

    [section]
    key = value

    [[subsection]]
    key = value

    [[[subsubsection]]]
    key = value

This module converts the resulting :class:`configobj.ConfigObj` tree into a
:class:`types.SimpleNamespace` tree so that callers can reach nested values with
attribute access (``cfg.section.subsection.key``) instead of dict indexing.

Scalars are best-effort coerced into their natural Python types
(``int``, ``float``, ``bool``, ``None``). Comma-separated values and configobj
list syntax are preserved as Python lists with per-element coercion.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any




__all__ = [
    "read_ini",
    "namespace_to_dict",
]


_TRUE_LITERALS: frozenset[str] = frozenset({"true", "yes", "on", "1"})
_FALSE_LITERALS: frozenset[str] = frozenset({"false", "no", "off", "0"})
_NONE_LITERALS: frozenset[str] = frozenset({"none", "null", "~"})


def _coerce_scalar(value: Any) -> Any:
    """
    Coerces a scalar value into its appropriate type, such as None, boolean, integer,
    or float, based on its content. If the value cannot be interpreted as any of
    those specific types, it is returned as a stripped string.

    :param value: The value to coerce. Can be of any type.
    :return: The coerced value, interpreted into the appropriate type (NoneType,
        bool, int, float, or str).
    """
    if not isinstance(value, str):
        return value

    text = value.strip()
    if text == "":
        return ""

    lowered = text.lower()
    if lowered in _NONE_LITERALS:
        return None
    if lowered in _TRUE_LITERALS:
        return True
    if lowered in _FALSE_LITERALS:
        return False

    try:
        return int(text)
    except ValueError:
        pass

    try:
        return float(text)
    except ValueError:
        pass

    return text


def _coerce_value(value: Any) -> Any:
    """
    Processes the input value by ensuring it conforms to a certain standard.
    If the input is a list, each element of the list is individually processed.
    Otherwise, the input is treated as a scalar and processed accordingly.

    :param value: The input to process, can be a scalar or a list of scalars.
    :type value: Any
    :return: The processed value, either a single scalar or a list of processed scalars.
    :rtype: Any
    """
    if isinstance(value, list):
        return [_coerce_scalar(item) for item in value]
    return _coerce_scalar(value)


def _section_to_namespace(section: Any) -> SimpleNamespace:
    """
    Recursively converts a structured section object into a SimpleNamespace.
    Each scalar value and nested section becomes an attribute of the returned
    SimpleNamespace instance.

    :param section: The input section object, containing scalars and nested
                    sections to be converted
    :type section: Any
    :return: A SimpleNamespace representation of the input section, including
             all scalar values and nested sections.
    :rtype: SimpleNamespace
    """
    ns = SimpleNamespace()
    for key in section.scalars:
        setattr(ns, key, _coerce_value(section[key]))
    for key in section.sections:
        setattr(ns, key, _section_to_namespace(section[key]))
    return ns


def read_ini(path: str | Path,
             *,
             encoding: str = "utf-8") -> SimpleNamespace:
    """
    Reads an INI file and converts its contents to a SimpleNamespace object.

    The function uses the `configobj` library to parse and read the INI file. It
    validates the file's existence and ensures proper parsing, raising exceptions
    in case of errors. The parsed INI file's sections and values are ultimately
    converted into a SimpleNamespace structure for easier access.

    :param path: The path to the INI file to be read. This can be a string
                 or a pathlib.Path object.
    :param encoding: The encoding used to read the INI file. Defaults to "utf-8".
    :return: A SimpleNamespace representation of the INI file's content.
    :raises ImportError: If the `configobj` library is not installed.
    :raises FileNotFoundError: If the specified INI file does not exist.
    :raises ValueError: If there are errors while parsing the INI file.
    """
    try:
        from configobj import ConfigObj, ConfigObjError
    except ImportError as exc:
        raise ImportError(
            "Reading Pharmacon INI files requires 'configobj'. "
            "Install with: pip install configobj"
        ) from exc

    ini_path = Path(path).expanduser().resolve()
    if not ini_path.is_file():
        raise FileNotFoundError(f"INI file not found: {ini_path}")

    try:
        config = ConfigObj(
            str(ini_path),
            encoding=encoding,
            file_error=True,
            raise_errors=True,
            list_values=True,
            interpolation=False,
        )
    except ConfigObjError as exc:
        raise ValueError(f"Failed to parse INI file '{ini_path}': {exc}") from exc

    return _section_to_namespace(config)


def namespace_to_dict(ns: SimpleNamespace) -> dict[str, Any]:
    """
    Converts a SimpleNamespace object into a dictionary recursively.

    This function allows representing the internal structure of a
    SimpleNamespace object using a nested dictionary format.
    If the namespace contains nested SimpleNamespace objects,
    they are recursively converted into dictionaries.

    :param ns: SimpleNamespace instance to be converted.
    :type ns: SimpleNamespace
    :return: A dictionary representation of the SimpleNamespace.
    :rtype: dict[str, Any]
    """
    out: dict[str, Any] = {}
    for key, value in vars(ns).items():
        if isinstance(value, SimpleNamespace):
            out[key] = namespace_to_dict(value)
        else:
            out[key] = value
    return out
