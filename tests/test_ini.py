"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.utils.ini — read_ini and namespace_to_dict
"""
import pytest
from pathlib import Path
from types import SimpleNamespace

from pharmacon.utils.ini import read_ini, namespace_to_dict


def _write_ini(path: Path, content: str) -> Path:
    path.write_text(content, encoding="utf-8")
    return path



class TestReadIniFileErrors:
    def test_missing_file_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_ini(tmp_path / "nonexistent.ini")

    def test_returns_simple_namespace(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[section]\nkey = value\n")
        result = read_ini(ini)
        assert isinstance(result, SimpleNamespace)


class TestScalarCoercion:
    def test_integer_coerced(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv = 42\n")
        ns = read_ini(ini)
        assert ns.s.v == 42
        assert isinstance(ns.s.v, int)

    def test_float_coerced(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv = 3.14\n")
        ns = read_ini(ini)
        assert abs(ns.s.v - 3.14) < 1e-9
        assert isinstance(ns.s.v, float)

    @pytest.mark.parametrize("literal", ["true", "yes", "on", "True", "YES"])
    def test_true_literals_coerced(self, tmp_path, literal):
        ini = _write_ini(tmp_path / "t.ini", f"[s]\nv = {literal}\n")
        ns = read_ini(ini)
        assert ns.s.v is True

    @pytest.mark.parametrize("literal", ["false", "no", "off", "False", "NO"])
    def test_false_literals_coerced(self, tmp_path, literal):
        ini = _write_ini(tmp_path / "t.ini", f"[s]\nv = {literal}\n")
        ns = read_ini(ini)
        assert ns.s.v is False

    @pytest.mark.parametrize("literal", ["none", "null", "~", "None"])
    def test_none_literals_coerced(self, tmp_path, literal):
        ini = _write_ini(tmp_path / "t.ini", f"[s]\nv = {literal}\n")
        ns = read_ini(ini)
        assert ns.s.v is None

    def test_plain_string_preserved(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv = hello world\n")
        ns = read_ini(ini)
        assert ns.s.v == "hello world"

    def test_empty_string_preserved(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv =\n")
        ns = read_ini(ini)
        assert ns.s.v == ""


class TestListValues:
    def test_comma_separated_list(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv = 1, 2, 3\n")
        ns = read_ini(ini)
        assert ns.s.v == [1, 2, 3]

    def test_mixed_type_list(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nv = hello, 42, true\n")
        ns = read_ini(ini)
        assert ns.s.v == ["hello", 42, True]


class TestNestedSections:
    def test_top_level_section_accessible(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[outer]\nkey = val\n")
        ns = read_ini(ini)
        assert ns.outer.key == "val"

    def test_nested_section_accessible(self, tmp_path):
        content = "[outer]\nk = 1\n[[inner]]\nk2 = 2\n"
        ini = _write_ini(tmp_path / "t.ini", content)
        ns = read_ini(ini)
        assert ns.outer.inner.k2 == 2

    def test_multiple_sections(self, tmp_path):
        content = "[A]\nv = 1\n[B]\nv = 2\n"
        ini = _write_ini(tmp_path / "t.ini", content)
        ns = read_ini(ini)
        assert ns.A.v == 1
        assert ns.B.v == 2


class TestNamespaceToDict:
    def test_returns_dict(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nk = 1\n")
        ns = read_ini(ini)
        result = namespace_to_dict(ns)
        assert isinstance(result, dict)

    def test_flat_namespace_to_dict(self, tmp_path):
        ini = _write_ini(tmp_path / "t.ini", "[s]\nk = 42\n")
        ns = read_ini(ini)
        d = namespace_to_dict(ns)
        assert d["s"]["k"] == 42

    def test_nested_namespace_to_dict(self, tmp_path):
        content = "[outer]\n[[inner]]\nv = 99\n"
        ini = _write_ini(tmp_path / "t.ini", content)
        ns = read_ini(ini)
        d = namespace_to_dict(ns)
        assert d["outer"]["inner"]["v"] == 99

    def test_empty_namespace(self):
        ns = SimpleNamespace()
        assert namespace_to_dict(ns) == {}

    def test_roundtrip_read_and_convert(self, tmp_path):
        content = "[PTA-UNIFIED]\nfig_dpi = 300\nenable_grid = true\n"
        ini = _write_ini(tmp_path / "t.ini", content)
        ns = read_ini(ini)
        d = namespace_to_dict(ns)
        assert d["PTA-UNIFIED"]["fig_dpi"] == 300
        assert d["PTA-UNIFIED"]["enable_grid"] is True
