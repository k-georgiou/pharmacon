"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Regression tests for RMSF ``shading`` and ``colors_by_label`` parsing.

Both fields use comma-separated, multi-value syntax. The INI loader
(configobj, ``list_values=True``) splits any comma value into a Python list,
so the raw value arriving at the validator is a list — which previously got
stringified (``str(list)``) and failed to parse. ``_listish_to_csv`` now
re-joins the list before parsing, mirroring the ``parse_range_dict`` fix.

(A ``#`` hex colour left unquoted in the INI is stripped as a configobj
comment *before* it reaches the settings layer; that is handled by the
documentation telling users to wrap the whole value in double quotes. These
tests cover the parsing once the value reaches the validator.)
"""

from __future__ import annotations

import pytest

from pharmacon.constants.plots import RMSFPlotSettings as R


class TestListishToCsv:
    def test_list_rejoined(self):
        assert R._listish_to_csv(["a", "b", "c"]) == "a, b, c"

    def test_tuple_rejoined(self):
        assert R._listish_to_csv(("0", "10", "gray")) == "0, 10, gray"

    def test_plain_string_passthrough(self):
        assert R._listish_to_csv("0,10,gray") == "0,10,gray"

    def test_none_becomes_str(self):
        assert R._listish_to_csv(None) == "None"


class TestShadingParsing:
    def test_string_form_two_regions(self):
        s = R.from_dict({"shading": "0,3,gray; 4,6,orange"})
        assert s._current_warnings == 0
        assert [(r["start"], r["end"], r["color"]) for r in s.shading_regions] == \
            [(0.0, 3.0, "gray"), (4.0, 6.0, "orange")]

    def test_configobj_split_list_rejoined(self):
        # configobj turns "0,3,gray; 4,6,orange" into this list:
        s = R.from_dict({"shading": ["0", "3", "gray; 4", "6", "orange"]})
        assert s._current_warnings == 0
        assert len(s.shading_regions) == 2
        assert s.shading_regions[0]["start"] == 0.0
        assert s.shading_regions[1]["end"] == 6.0

    def test_alpha_and_label(self):
        s = R.from_dict({"shading": "0,3,gray,0.4,N-term"})
        reg = s.shading_regions[0]
        assert reg["alpha"] == 0.4 and reg["label"] == "N-term"

    def test_empty_is_no_regions(self):
        s = R.from_dict({"shading": ""})
        assert s.shading_regions == []
        assert s._current_warnings == 0


class TestColorsByLabelParsing:
    def test_string_form(self):
        s = R.from_dict({"colors_by_label": "calpha:red, backbone:blue"})
        assert s._current_warnings == 0
        assert s.colors_by_label_map == {"calpha": "red", "backbone": "blue"}

    def test_configobj_split_list_rejoined(self):
        s = R.from_dict({"colors_by_label": ["calpha:red", "backbone:blue"]})
        assert s._current_warnings == 0
        assert s.colors_by_label_map == {"calpha": "red", "backbone": "blue"}

    def test_quoted_hex_string(self):
        # The working INI form arrives here as a single string.
        s = R.from_dict({"colors_by_label": "calpha:#d62728, backbone:#1f77b4"})
        assert s.colors_by_label_map == {"calpha": "#d62728", "backbone": "#1f77b4"}

    def test_malformed_entry_skipped_with_warning(self):
        s = R.from_dict({"colors_by_label": "calpha-no-colon"})
        assert s.colors_by_label_map == {}
        assert s._current_warnings >= 1

    def test_empty_is_empty_map(self):
        s = R.from_dict({"colors_by_label": ""})
        assert s.colors_by_label_map == {}
        assert s._current_warnings == 0
