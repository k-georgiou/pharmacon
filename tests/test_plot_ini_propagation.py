"""For every example INI under examples/plot_ini/, verify that every
section maps to a known settings class, every key in every section is a
real field on that class, and after the real CLI pipeline runs
(`read_ini -> namespace_to_dict -> _flatten_section -> cls.from_dict`)
each key's value lands on the settings instance verbatim.

This is the exhaustive end-to-end check: if any INI line is silently
dropped or mutated, one of these tests will fail and name the line.
"""

from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any

import pytest

from pharmacon.utils.ini import read_ini, namespace_to_dict
from pharmacon.constants.plots import PlotSettingsBase
from pharmacon.command_line.plot.pta import _flatten_section


REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = REPO_ROOT / "examples" / "plot_ini"


def _resolve_class(section_name: str):
    """Look up the settings class registered under any of its aliases.

    PlotSettingsBase stores aliases uppercased at registration, so a
    direct uppercased lookup is sufficient.
    """
    return PlotSettingsBase._registry.get(section_name.upper())


def _collect_section_cases() -> list:
    """Build pytest parametrize cases for every (ini, section) pair that
    resolves to a known settings class.
    """
    cases = []
    for ini_path in sorted(EXAMPLES_DIR.glob("*.ini")):
        ns = read_ini(ini_path)
        d = namespace_to_dict(ns)
        for section_name, section_body in d.items():
            if not isinstance(section_body, dict):
                continue
            cls = _resolve_class(section_name)
            if cls is None:
                # Unresolved sections are flagged by a separate test below.
                continue
            flat = _flatten_section(section_body)
            cases.append(
                pytest.param(
                    ini_path.name, section_name, cls, flat,
                    id=f"{ini_path.stem}::{section_name}",
                )
            )
    return cases


SECTION_CASES = _collect_section_cases()


# ---------------------------------------------------------------------------
# 1) every section resolves to a known settings class
# ---------------------------------------------------------------------------

class TestSectionRecognition:
    def test_at_least_one_ini_with_sections(self):
        assert SECTION_CASES, (
            "no parametrize cases collected — examples/plot_ini/ looks empty"
        )

    def test_no_unrecognized_sections_in_examples(self):
        """Catches typos in section headers. Every dict section in every
        example INI must resolve to a PlotSettingsBase subclass."""
        unrecognized = []
        for ini_path in sorted(EXAMPLES_DIR.glob("*.ini")):
            ns = read_ini(ini_path)
            d = namespace_to_dict(ns)
            for section_name, section_body in d.items():
                if not isinstance(section_body, dict):
                    continue
                if _resolve_class(section_name) is None:
                    unrecognized.append((ini_path.name, section_name))
        assert not unrecognized, (
            "sections with no matching settings class:\n"
            + "\n".join(f"  {p}  [{s}]" for p, s in unrecognized)
        )


# ---------------------------------------------------------------------------
# 2) every key in the INI section is a real field on the settings class
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ini_file,section_name,cls,overrides", SECTION_CASES)
def test_no_unknown_keys_in_section(ini_file, section_name, cls, overrides):
    """Catches typos/dead keys in the example INI itself: every key under
    a section must correspond to a dataclass field on the section's
    settings class. Unknown keys would be silently dropped by `from_dict`,
    which is a feature for forward-compat but a bug in shipped examples."""
    field_names = {f.name for f in fields(cls)}
    unknown = sorted(k for k in overrides if k not in field_names)
    assert not unknown, (
        f"{ini_file} [{section_name}]: {len(unknown)} key(s) not on "
        f"{cls.__name__}: {unknown}"
    )


# ---------------------------------------------------------------------------
# 3) every known key propagates from INI to settings instance verbatim
# ---------------------------------------------------------------------------

def _equal_for_settings(a: Any, b: Any) -> bool:
    """Treat list/tuple as equal element-wise so `[..., ...]` from the
    INI matches an attribute declared as a tuple in the dataclass."""
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return list(a) == list(b)
    return a == b


@pytest.mark.parametrize("ini_file,section_name,cls,overrides", SECTION_CASES)
def test_every_key_propagates(ini_file, section_name, cls, overrides):
    """For each key in the section that maps to a real field, the
    rendered settings instance must carry the same value the INI
    declared (after the INI parser's type coercion)."""
    field_names = {f.name for f in fields(cls)}
    instance = cls.from_dict(overrides)

    mismatches = []
    for key, raw_value in overrides.items():
        if key not in field_names:
            continue
        # Empty-string values in the INI are a documented "use the
        # default" sentinel — validate() restores the dataclass default
        # for that field. We're testing propagation of EXPLICIT values,
        # so skip the blanks.
        if isinstance(raw_value, str) and raw_value == "":
            continue
        actual = getattr(instance, key)
        if not _equal_for_settings(actual, raw_value):
            mismatches.append((key, raw_value, actual))

    assert not mismatches, (
        f"{ini_file} [{section_name}]: {len(mismatches)} key(s) did not "
        f"propagate to {cls.__name__}:\n"
        + "\n".join(
            f"  {k!s:30s}  INI={v!r}  ->  instance={a!r}"
            for k, v, a in mismatches
        )
    )
