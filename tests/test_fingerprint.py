"""
Tests for pharmacon.utils.fingerprint
"""
import pytest

from pharmacon.utils.fingerprint import (
    PharmaconFileSignature,
    _normalize_token,
    _chunk_text,
    create_pharmacon_signature,
)


# ---------------------------------------------------------------------------
# _normalize_token
# ---------------------------------------------------------------------------

class TestNormalizeToken:
    def test_returns_uppercase(self):
        assert _normalize_token("hello", "f") == "HELLO"

    def test_strips_whitespace(self):
        assert _normalize_token("  world  ", "f") == "WORLD"

    def test_mixed_case_and_whitespace(self):
        assert _normalize_token("  TrAjEcToRy  ", "f") == "TRAJECTORY"

    def test_non_string_raises_type_error(self):
        with pytest.raises(TypeError, match="must be a string"):
            _normalize_token(123, "f")

    def test_none_raises_type_error(self):
        with pytest.raises(TypeError):
            _normalize_token(None, "f")

    def test_empty_string_raises_value_error(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            _normalize_token("", "f")

    def test_whitespace_only_raises_value_error(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            _normalize_token("   ", "f")

    def test_null_char_raises_value_error(self):
        with pytest.raises(ValueError, match="control characters"):
            _normalize_token("he\x00llo", "f")

    def test_newline_raises_value_error(self):
        with pytest.raises(ValueError, match="control characters"):
            _normalize_token("he\nllo", "f")

    def test_carriage_return_raises_value_error(self):
        with pytest.raises(ValueError, match="control characters"):
            _normalize_token("he\rllo", "f")


# ---------------------------------------------------------------------------
# _chunk_text
# ---------------------------------------------------------------------------

class TestChunkText:
    def test_default_chunk_size_4(self):
        assert _chunk_text("ABCDEFGH") == "ABCD-EFGH"

    def test_custom_chunk_size(self):
        assert _chunk_text("ABCDEF", chunk_size=2) == "AB-CD-EF"

    def test_custom_separator(self):
        assert _chunk_text("ABCDEF", chunk_size=3, separator=":") == "ABC:DEF"

    def test_empty_separator(self):
        assert _chunk_text("ABCDE", chunk_size=2, separator="") == "ABCDE"

    def test_empty_string_returns_empty(self):
        assert _chunk_text("") == ""

    def test_string_shorter_than_chunk(self):
        assert _chunk_text("AB", chunk_size=8) == "AB"

    def test_exact_multiple(self):
        assert _chunk_text("ABCDEF", chunk_size=3) == "ABC-DEF"

    def test_non_multiple_length(self):
        assert _chunk_text("ABCDE", chunk_size=4) == "ABCD-E"

    def test_zero_chunk_size_raises(self):
        with pytest.raises(ValueError, match="chunk_size must be > 0"):
            _chunk_text("ABC", chunk_size=0)

    def test_negative_chunk_size_raises(self):
        with pytest.raises(ValueError):
            _chunk_text("ABC", chunk_size=-1)


# ---------------------------------------------------------------------------
# create_pharmacon_signature
# ---------------------------------------------------------------------------

class TestCreatePharmaconSignature:
    def test_returns_dataclass_instance(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        assert isinstance(sig, PharmaconFileSignature)

    def test_magic_is_pharmacon(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        assert sig.magic == "PHARMACON"

    def test_inputs_uppercased(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        assert sig.format_name == "PTA"
        assert sig.command == "TRAJECTORY"
        assert sig.subcommand == "RMSD"

    def test_signature_format_pharmacon_cmd_subcmd_body(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        parts = sig.signature.split("::")
        assert len(parts) == 4
        assert parts[0] == "PHARMACON"
        assert parts[1] == "TRAJECTORY"
        assert parts[2] == "RMSD"

    def test_deterministic_same_inputs(self):
        kwargs = dict(format_name="pta", command="trajectory", subcommand="rmsd")
        sig1 = create_pharmacon_signature(**kwargs)
        sig2 = create_pharmacon_signature(**kwargs)
        assert sig1.signature == sig2.signature
        assert sig1.fingerprint == sig2.fingerprint

    def test_different_commands_produce_different_signatures(self):
        sig1 = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        sig2 = create_pharmacon_signature(
            format_name="pta", command="structure", subcommand="rmsd"
        )
        assert sig1.signature != sig2.signature

    def test_different_subcommands_produce_different_signatures(self):
        sig1 = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        sig2 = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="distances"
        )
        assert sig1.signature != sig2.signature

    def test_fingerprint_is_hex_string(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="trajectory", subcommand="rmsd"
        )
        int(sig.fingerprint, 16)  # raises ValueError if not valid hex

    def test_empty_version_raises(self):
        with pytest.raises(ValueError, match="version cannot be empty"):
            create_pharmacon_signature(
                format_name="pta", command="traj", subcommand="rmsd", version=""
            )

    def test_whitespace_version_raises(self):
        with pytest.raises(ValueError, match="version cannot be empty"):
            create_pharmacon_signature(
                format_name="pta", command="traj", subcommand="rmsd", version="   "
            )

    def test_empty_format_name_raises(self):
        with pytest.raises((TypeError, ValueError)):
            create_pharmacon_signature(
                format_name="", command="traj", subcommand="rmsd"
            )

    def test_empty_command_raises(self):
        with pytest.raises((TypeError, ValueError)):
            create_pharmacon_signature(
                format_name="pta", command="", subcommand="rmsd"
            )

    def test_version_stored_on_dataclass(self):
        sig = create_pharmacon_signature(
            format_name="pta", command="traj", subcommand="rmsd", version="2.0.0"
        )
        assert sig.version == "2.0.0"
