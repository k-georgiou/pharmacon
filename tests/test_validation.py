"""
Tests for pharmacon.utils.validation — all 11 public functions.
"""
import pytest
from pathlib import Path

from pharmacon.utils.validation import (
    normalize_path,
    validate_existing_input_file,
    validate_output_file,
    normalize_selection,
    validate_string_list,
    validate_non_negative_int,
    validate_positive_int,
    validate_bool_flag,
    validate_logging_level,
    validate_frame_range,
    validate_reference_frame,
)
from pharmacon.command_line.exceptions import ValidationError


# ---------------------------------------------------------------------------
# normalize_path
# ---------------------------------------------------------------------------

class TestNormalizePath:
    def test_string_returns_path(self):
        result = normalize_path("/tmp/test.txt", "path")
        assert isinstance(result, Path)

    def test_path_object_accepted(self):
        result = normalize_path(Path("/tmp/test.txt"), "path")
        assert isinstance(result, Path)

    def test_strips_whitespace(self):
        result = normalize_path("  /tmp/test.txt  ", "path")
        assert str(result) == str(Path("/tmp/test.txt"))

    def test_expands_tilde(self):
        result = normalize_path("~/some_file.txt", "path")
        assert "~" not in str(result)

    @pytest.mark.parametrize("bad", [123, None, [], True])
    def test_non_path_type_raises(self, bad):
        with pytest.raises(ValidationError, match="path-like"):
            normalize_path(bad, "path")

    def test_empty_string_raises(self):
        with pytest.raises(ValidationError, match="cannot be empty"):
            normalize_path("", "path")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValidationError, match="cannot be empty"):
            normalize_path("   ", "path")

    def test_null_char_raises(self):
        with pytest.raises(ValidationError, match="control characters"):
            normalize_path("/tmp/bad\x00file", "path")

    def test_newline_raises(self):
        with pytest.raises(ValidationError, match="control characters"):
            normalize_path("/tmp/bad\nfile", "path")

    def test_carriage_return_raises(self):
        with pytest.raises(ValidationError, match="control characters"):
            normalize_path("/tmp/bad\rfile", "path")


# ---------------------------------------------------------------------------
# validate_existing_input_file
# ---------------------------------------------------------------------------

class TestValidateExistingInputFile:
    def test_valid_file_returns_resolved_path(self, tmp_path):
        f = tmp_path / "topology.tpr"
        f.write_bytes(b"data")
        result = validate_existing_input_file(f, "topology", (".tpr",))
        assert result == f.resolve()

    def test_format_without_leading_dot(self, tmp_path):
        f = tmp_path / "input.pdb"
        f.write_bytes(b"data")
        result = validate_existing_input_file(f, "topology", ("pdb",))
        assert result == f.resolve()

    def test_case_insensitive_suffix(self, tmp_path):
        f = tmp_path / "input.PDB"
        f.write_bytes(b"data")
        result = validate_existing_input_file(f, "topology", (".pdb",))
        assert result == f.resolve()

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(ValidationError, match="does not exist"):
            validate_existing_input_file(tmp_path / "ghost.tpr", "topology", (".tpr",))

    def test_wrong_suffix_raises(self, tmp_path):
        f = tmp_path / "file.xyz"
        f.write_bytes(b"data")
        with pytest.raises(ValidationError, match="Unsupported"):
            validate_existing_input_file(f, "topology", (".tpr", ".pdb"))

    def test_directory_as_file_raises(self, tmp_path):
        with pytest.raises(ValidationError, match="not a regular file"):
            validate_existing_input_file(tmp_path, "topology", ("",))


# ---------------------------------------------------------------------------
# validate_output_file
# ---------------------------------------------------------------------------

class TestValidateOutputFile:
    def test_new_file_accepted(self, tmp_path):
        out = tmp_path / "result.pta"
        result = validate_output_file(out, "output", overwrite=False)
        assert result == out.resolve()

    def test_existing_file_without_overwrite_raises(self, tmp_path):
        out = tmp_path / "result.pta"
        out.write_bytes(b"data")
        with pytest.raises(ValidationError, match="already exists"):
            validate_output_file(out, "output", overwrite=False)

    def test_existing_file_with_overwrite_accepted(self, tmp_path):
        out = tmp_path / "result.pta"
        out.write_bytes(b"data")
        result = validate_output_file(out, "output", overwrite=True)
        assert result == out.resolve()

    def test_suffix_validation_pass(self, tmp_path):
        out = tmp_path / "result.pta"
        result = validate_output_file(out, "output", overwrite=False, allowed_suffixes=(".pta",))
        assert result == out.resolve()

    def test_suffix_validation_fail_raises(self, tmp_path):
        out = tmp_path / "result.xyz"
        with pytest.raises(ValidationError, match="Unsupported"):
            validate_output_file(out, "output", overwrite=False, allowed_suffixes=(".pta",))

    def test_nonexistent_parent_raises(self, tmp_path):
        out = tmp_path / "nonexistent_dir" / "result.pta"
        with pytest.raises(ValidationError, match="Parent directory"):
            validate_output_file(out, "output", overwrite=False)

    def test_pointing_to_directory_raises(self, tmp_path):
        with pytest.raises(ValidationError):
            validate_output_file(tmp_path, "output", overwrite=False)


# ---------------------------------------------------------------------------
# normalize_selection
# ---------------------------------------------------------------------------

class TestNormalizeSelection:
    def test_valid_string_returned(self):
        result = normalize_selection("protein", "sel", required=True)
        assert result == "protein"

    def test_strips_whitespace(self):
        result = normalize_selection("  resname LIG  ", "sel", required=True)
        assert result == "resname LIG"

    def test_none_with_default_returns_default(self):
        result = normalize_selection(None, "sel", required=False, default="all")
        assert result == "all"

    def test_none_required_no_default_raises(self):
        with pytest.raises(ValidationError):
            normalize_selection(None, "sel", required=True)

    def test_blank_required_no_default_raises(self):
        with pytest.raises(ValidationError):
            normalize_selection("   ", "sel", required=True)

    def test_blank_not_required_returns_default(self):
        result = normalize_selection("   ", "sel", required=False, default="backbone")
        assert result == "backbone"

    def test_non_string_raises(self):
        with pytest.raises(ValidationError):
            normalize_selection(123, "sel", required=True)

    def test_control_char_raises(self):
        with pytest.raises(ValidationError, match="control characters"):
            normalize_selection("protein\x00", "sel", required=True)

    def test_none_not_required_no_default_returns_none(self):
        result = normalize_selection(None, "sel", required=False)
        assert result is None


# ---------------------------------------------------------------------------
# validate_string_list
# ---------------------------------------------------------------------------

class TestValidateStringList:
    def test_valid_list_returned(self):
        result = validate_string_list(["alpha", "beta"], "labels")
        assert result == ["alpha", "beta"]

    def test_strips_whitespace_from_items(self):
        result = validate_string_list(["  A  ", "B  "], "labels")
        assert result == ["A", "B"]

    def test_empty_list_raises(self):
        with pytest.raises(ValidationError):
            validate_string_list([], "labels")

    def test_not_a_list_raises(self):
        with pytest.raises(ValidationError):
            validate_string_list("alpha", "labels")

    def test_item_not_string_raises(self):
        with pytest.raises(ValidationError):
            validate_string_list(["alpha", 42], "labels")

    def test_whitespace_only_item_raises(self):
        with pytest.raises(ValidationError):
            validate_string_list(["alpha", "   "], "labels")


# ---------------------------------------------------------------------------
# validate_non_negative_int
# ---------------------------------------------------------------------------

class TestValidateNonNegativeInt:
    @pytest.mark.parametrize("val", [0, 1, 100])
    def test_valid_values(self, val):
        assert validate_non_negative_int(val, "n") == val

    def test_negative_raises(self):
        with pytest.raises(ValidationError, match=">= 0"):
            validate_non_negative_int(-1, "n")

    def test_bool_raises(self):
        with pytest.raises(ValidationError, match="integer"):
            validate_non_negative_int(True, "n")

    def test_float_raises(self):
        with pytest.raises(ValidationError, match="integer"):
            validate_non_negative_int(1.0, "n")

    def test_string_raises(self):
        with pytest.raises(ValidationError):
            validate_non_negative_int("1", "n")


# ---------------------------------------------------------------------------
# validate_positive_int
# ---------------------------------------------------------------------------

class TestValidatePositiveInt:
    @pytest.mark.parametrize("val", [1, 2, 999])
    def test_valid_values(self, val):
        assert validate_positive_int(val, "n") == val

    def test_zero_raises(self):
        with pytest.raises(ValidationError, match="> 0"):
            validate_positive_int(0, "n")

    def test_negative_raises(self):
        with pytest.raises(ValidationError):
            validate_positive_int(-1, "n")

    def test_bool_raises(self):
        with pytest.raises(ValidationError):
            validate_positive_int(True, "n")

    def test_float_raises(self):
        with pytest.raises(ValidationError):
            validate_positive_int(1.0, "n")


# ---------------------------------------------------------------------------
# validate_bool_flag
# ---------------------------------------------------------------------------

class TestValidateBoolFlag:
    def test_true_accepted(self):
        assert validate_bool_flag(True, "flag") is True

    def test_false_accepted(self):
        assert validate_bool_flag(False, "flag") is False

    @pytest.mark.parametrize("bad", [1, 0, "true", "false", None, 1.0])
    def test_non_bool_raises(self, bad):
        with pytest.raises(ValidationError, match="boolean"):
            validate_bool_flag(bad, "flag")


# ---------------------------------------------------------------------------
# validate_logging_level
# ---------------------------------------------------------------------------

class TestValidateLoggingLevel:
    @pytest.mark.parametrize("level", ["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    def test_valid_levels(self, level):
        assert validate_logging_level(level, "level") == level

    def test_case_insensitive(self):
        assert validate_logging_level("debug", "level") == "DEBUG"

    def test_strips_whitespace(self):
        assert validate_logging_level("  INFO  ", "level") == "INFO"

    def test_invalid_level_raises(self):
        with pytest.raises(ValidationError, match="Invalid logging level"):
            validate_logging_level("VERBOSE", "level")

    def test_non_string_raises(self):
        with pytest.raises(ValidationError):
            validate_logging_level(10, "level")


# ---------------------------------------------------------------------------
# validate_frame_range
# ---------------------------------------------------------------------------

class TestValidateFrameRange:
    def test_valid_defaults(self):
        data = {"begin": 0, "end": None, "step": 1}
        begin, end, step = validate_frame_range(data)
        assert begin == 0
        assert end is None
        assert step == 1

    def test_valid_range(self):
        data = {"begin": 5, "end": 100, "step": 2}
        begin, end, step = validate_frame_range(data)
        assert (begin, end, step) == (5, 100, 2)

    def test_mutates_data_dict(self):
        data = {"begin": 0, "end": 10, "step": 1}
        validate_frame_range(data)
        assert data["begin"] == 0
        assert data["end"] == 10
        assert data["step"] == 1

    def test_end_less_than_begin_raises(self):
        data = {"begin": 10, "end": 5, "step": 1}
        with pytest.raises(ValidationError, match="greater than or equal"):
            validate_frame_range(data)

    def test_end_equals_begin_accepted(self):
        data = {"begin": 5, "end": 5, "step": 1}
        begin, end, step = validate_frame_range(data)
        assert begin == end == 5

    def test_negative_begin_raises(self):
        data = {"begin": -1, "end": 10, "step": 1}
        with pytest.raises(ValidationError):
            validate_frame_range(data)

    def test_zero_step_raises(self):
        data = {"begin": 0, "end": 10, "step": 0}
        with pytest.raises(ValidationError):
            validate_frame_range(data)

    def test_negative_step_raises(self):
        data = {"begin": 0, "end": 10, "step": -1}
        with pytest.raises(ValidationError):
            validate_frame_range(data)

    def test_missing_keys_use_defaults(self):
        data = {}
        begin, end, step = validate_frame_range(data)
        assert begin == 0
        assert end is None
        assert step == 1


# ---------------------------------------------------------------------------
# validate_reference_frame
# ---------------------------------------------------------------------------

class TestValidateReferenceFrame:
    def test_valid_reference_frame(self):
        data = {"reference_frame": 5}
        result = validate_reference_frame(data, begin=0, end=10)
        assert result == 5

    def test_mutates_data_dict(self):
        data = {"reference_frame": 3}
        validate_reference_frame(data, begin=0, end=10)
        assert data["reference_frame"] == 3

    def test_default_is_zero(self):
        data = {}
        result = validate_reference_frame(data, begin=0, end=10)
        assert result == 0

    def test_ref_at_begin_accepted(self):
        data = {"reference_frame": 0}
        result = validate_reference_frame(data, begin=0, end=10)
        assert result == 0

    def test_ref_at_end_accepted(self):
        data = {"reference_frame": 10}
        result = validate_reference_frame(data, begin=0, end=10)
        assert result == 10

    def test_ref_below_begin_raises(self):
        data = {"reference_frame": 1}
        with pytest.raises(ValidationError, match=">= begin"):
            validate_reference_frame(data, begin=5, end=10)

    def test_ref_above_end_raises(self):
        data = {"reference_frame": 15}
        with pytest.raises(ValidationError, match="<= end"):
            validate_reference_frame(data, begin=0, end=10)

    def test_end_none_no_upper_check(self):
        data = {"reference_frame": 1000}
        result = validate_reference_frame(data, begin=0, end=None)
        assert result == 1000

    def test_negative_ref_raises(self):
        data = {"reference_frame": -1}
        with pytest.raises(ValidationError):
            validate_reference_frame(data, begin=0, end=10)

    def test_custom_field_name(self):
        data = {"ref": 5}
        result = validate_reference_frame(data, begin=0, end=10, field_name="ref")
        assert result == 5
