"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Contract tests between ``pharmacon.constants.plots`` settings classes and the
``pharmacon.plotter`` modules.

These tests catch the class of bug where a plotter reads ``settings.foo`` but
the dataclass field was renamed to ``bar`` (e.g. ``bar_edgecolor`` vs
``bar_edge_color``, ``fontsize_title`` vs ``font_size_title``). Rather than
writing one assertion per attribute, we:

1. Walk every plotter module with ``ast`` and collect every attribute read
   off a variable literally named ``settings``.
2. Build the union of field names across every concrete
   ``PlotSettingsBase`` subclass.
3. Fail if a plotter references a name that no settings class exposes.

We also sanity-check that every settings subclass:

- is a ``dataclass`` (so ``cls()`` and ``from_dict({})`` work),
- has a non-empty ``alias`` tuple,
- survives ``from_dict({})`` and ``validate()`` without exceptions.
"""

from __future__ import annotations

import ast
import inspect
from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import Iterable

import pytest

from pharmacon import plotter as _plotter_pkg
from pharmacon.constants import plots as _plots_mod
from pharmacon.constants.plots import PlotSettingsBase


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Methods/attributes defined on ``PlotSettingsBase`` that plotters legitimately
# touch through a ``settings`` reference but which aren't dataclass fields.
_BASE_ALLOWED: frozenset[str] = frozenset(
    {
        "alias",
        "validate",
        "from_dict",
        "namespace",
        "_validate_fields",
        "_current_warnings",
        "_warn",
        "_registry",
    }
)

# Attribute names the plotters touch on objects that happen to be bound to
# a local called ``settings`` but are not PlotSettings (e.g. matplotlib
# ``Text`` or ``Legend`` handles). Extend this list only when a false
# positive is unavoidable.
_FALSE_POSITIVES: frozenset[str] = frozenset(set())


def _iter_settings_classes() -> Iterable[type[PlotSettingsBase]]:
    for name in dir(_plots_mod):
        obj = getattr(_plots_mod, name)
        if (
            inspect.isclass(obj)
            and issubclass(obj, PlotSettingsBase)
            and obj is not PlotSettingsBase
        ):
            yield obj


def _all_known_field_names() -> set[str]:
    names: set[str] = set(_BASE_ALLOWED)
    for cls in _iter_settings_classes():
        if not is_dataclass(cls):
            continue
        names.update(f.name for f in fields(cls))
    return names


def _collect_settings_attribute_reads(py_path: Path) -> set[str]:
    """Return every attr name read off a variable called ``settings``."""
    tree = ast.parse(py_path.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id == "settings"
        ):
            names.add(node.attr)
    return names


def _plotter_files() -> list[Path]:
    root = Path(_plotter_pkg.__file__).parent
    return sorted(p for p in root.glob("*.py") if p.name != "__init__.py")


# ---------------------------------------------------------------------------
# 1. Every plotter-file reference must resolve to a real settings field.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("py_file", _plotter_files(), ids=lambda p: p.name)
def test_plotter_settings_references_resolve(py_file: Path) -> None:
    used = _collect_settings_attribute_reads(py_file)
    known = _all_known_field_names()
    missing = used - known - _FALSE_POSITIVES
    assert not missing, (
        f"{py_file.name} references unknown settings attribute(s): "
        f"{sorted(missing)}. Either the plotter has a typo or the field "
        f"was renamed without updating the plotter."
    )


# ---------------------------------------------------------------------------
# 2. Every settings class must be constructible and self-validating.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cls",
    list(_iter_settings_classes()),
    ids=lambda c: c.__name__,
)
def test_settings_class_contract(cls: type[PlotSettingsBase]) -> None:
    assert is_dataclass(cls), f"{cls.__name__} must be a @dataclass"

    assert cls.alias, f"{cls.__name__} has empty alias tuple"
    for a in cls.alias:
        assert isinstance(a, str) and a, (
            f"{cls.__name__} alias must be non-empty strings, got {cls.alias!r}"
        )

    # Default construction + no-op override round-trip.
    instance = cls.from_dict({})
    assert instance is not None

    # validate() is idempotent and leaves the instance in a clean state.
    instance.validate()
    instance.validate()


# ---------------------------------------------------------------------------
# 3. Alias registry uniqueness — fail fast if two classes claim the same key.
# ---------------------------------------------------------------------------


def test_alias_registry_has_no_duplicates() -> None:
    seen: dict[str, str] = {}
    for cls in _iter_settings_classes():
        for a in cls.alias:
            key = a.upper()
            if key in seen:
                pytest.fail(
                    f"Duplicate alias {key!r}: both "
                    f"{seen[key]} and {cls.__name__} register it"
                )
            seen[key] = cls.__name__
    assert seen, "No aliases registered — registry looks empty"
