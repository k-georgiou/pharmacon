"""
Tests for pharmacon.utils.identifiers:
  generate_mda_blueprint, create_mda_artifact_token, validate_mda_artifact_token
"""
import pytest
import numpy as np
import MDAnalysis as Mda
from MDAnalysis.coordinates.memory import MemoryReader

from pharmacon.utils.identifiers import (
    generate_mda_blueprint,
    create_mda_artifact_token,
    validate_mda_artifact_token,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_universe(n_atoms: int = 3, n_frames: int = 1) -> Mda.Universe:
    u = Mda.Universe.empty(
        n_atoms=n_atoms, n_residues=1, n_segments=1,
        atom_resindex=[0] * n_atoms, residue_segindex=[0], trajectory=True,
    )
    u.add_TopologyAttr("names", [f"A{i}" for i in range(n_atoms)])
    u.add_TopologyAttr("resnames", ["MOL"])
    u.add_TopologyAttr("resids", [1])
    u.add_TopologyAttr("segids", ["SYS"])
    positions = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    u.load_new(positions, format=MemoryReader)
    return u


# ---------------------------------------------------------------------------
# generate_mda_blueprint — return type and format
# ---------------------------------------------------------------------------

class TestGenerateMdaBlueprint:
    def test_returns_string_by_default(self):
        u = _make_universe()
        bp = generate_mda_blueprint(u=u)
        assert isinstance(bp, str)

    def test_non_empty_string(self):
        u = _make_universe()
        bp = generate_mda_blueprint(u=u)
        assert len(bp) > 0

    def test_include_debug_returns_dict(self):
        u = _make_universe()
        bp = generate_mda_blueprint(u=u, include_debug=True)
        assert isinstance(bp, dict)
        assert "id" in bp

    def test_debug_dict_contains_counts(self):
        u = _make_universe(n_atoms=3)
        bp = generate_mda_blueprint(u=u, include_debug=True)
        assert bp["counts"]["atoms"] == 3

    def test_hex_encoding_produces_hex_string(self):
        u = _make_universe()
        bp = generate_mda_blueprint(u=u, encoding="hex")
        assert isinstance(bp, str)
        int(bp, 16)  # must be valid hex

    def test_base32_encoding(self):
        u = _make_universe()
        bp = generate_mda_blueprint(u=u, encoding="base32")
        assert isinstance(bp, str)
        assert len(bp) > 0

    def test_invalid_encoding_raises(self):
        u = _make_universe()
        with pytest.raises(ValueError, match="encoding"):
            generate_mda_blueprint(u=u, encoding="latin1")

    def test_id_bits_not_multiple_of_8_raises(self):
        u = _make_universe()
        with pytest.raises(ValueError, match="id_bits"):
            generate_mda_blueprint(u=u, id_bits=13)

    def test_id_bits_out_of_range_raises(self):
        u = _make_universe()
        with pytest.raises(ValueError):
            generate_mda_blueprint(u=u, id_bits=512)


# ---------------------------------------------------------------------------
# generate_mda_blueprint — determinism
# ---------------------------------------------------------------------------

class TestBlueprintDeterminism:
    def test_same_universe_same_blueprint(self):
        u1 = _make_universe()
        u2 = _make_universe()
        assert generate_mda_blueprint(u=u1) == generate_mda_blueprint(u=u2)

    def test_different_atom_count_different_blueprint(self):
        u1 = _make_universe(n_atoms=3)
        u2 = _make_universe(n_atoms=5)
        assert generate_mda_blueprint(u=u1) != generate_mda_blueprint(u=u2)

    def test_kwargs_affect_blueprint(self):
        u = _make_universe()
        bp1 = generate_mda_blueprint(u=u, selection="protein")
        bp2 = generate_mda_blueprint(u=u, selection="ligand")
        assert bp1 != bp2

    def test_hash_keys_filters_kwargs(self):
        u = _make_universe()
        bp1 = generate_mda_blueprint(u=u, __hash_keys__=["selection"], selection="A", label="X")
        bp2 = generate_mda_blueprint(u=u, __hash_keys__=["selection"], selection="A", label="Y")
        assert bp1 == bp2

    def test_relax_resids_affects_blueprint(self):
        u = _make_universe()
        bp_strict = generate_mda_blueprint(u=u, relax_resids=False)
        bp_relaxed = generate_mda_blueprint(u=u, relax_resids=True)
        assert bp_strict != bp_relaxed


# ---------------------------------------------------------------------------
# create_mda_artifact_token
# ---------------------------------------------------------------------------

class TestCreateArtifactToken:
    def test_returns_string(self):
        token = create_mda_artifact_token(blueprint="abc123", secret="secret")
        assert isinstance(token, str)

    def test_non_empty(self):
        token = create_mda_artifact_token(blueprint="abc123", secret="secret")
        assert len(token) > 0

    def test_dict_blueprint_with_id_key(self):
        token = create_mda_artifact_token(blueprint={"id": "abc123"}, secret="secret")
        assert isinstance(token, str)

    def test_dict_blueprint_missing_id_raises(self):
        with pytest.raises(ValueError, match="id"):
            create_mda_artifact_token(blueprint={"x": "y"}, secret="secret")

    def test_invalid_blueprint_type_raises(self):
        with pytest.raises(TypeError):
            create_mda_artifact_token(blueprint=42, secret="secret")  # type: ignore

    def test_token_bits_not_multiple_of_8_raises(self):
        with pytest.raises(ValueError, match="token_bits"):
            create_mda_artifact_token(blueprint="abc", secret="s", token_bits=13)

    def test_invalid_encoding_raises(self):
        with pytest.raises(ValueError, match="encoding"):
            create_mda_artifact_token(blueprint="abc", secret="s", encoding="latin1")

    def test_include_debug_returns_dict(self):
        result = create_mda_artifact_token(blueprint="abc123", secret="secret", include_debug=True)
        assert isinstance(result, dict)
        assert "artifact_token" in result

    def test_deterministic_for_same_inputs(self):
        t1 = create_mda_artifact_token(blueprint="bp", secret="sec")
        t2 = create_mda_artifact_token(blueprint="bp", secret="sec")
        assert t1 == t2

    def test_different_secrets_produce_different_tokens(self):
        t1 = create_mda_artifact_token(blueprint="bp", secret="sec1")
        t2 = create_mda_artifact_token(blueprint="bp", secret="sec2")
        assert t1 != t2


# ---------------------------------------------------------------------------
# validate_mda_artifact_token
# ---------------------------------------------------------------------------

class TestValidateArtifactToken:
    def test_correct_token_validates_true(self):
        token = create_mda_artifact_token(blueprint="bp", secret="sec")
        assert validate_mda_artifact_token(artifact_token=token, blueprint="bp", secret="sec")

    def test_wrong_token_validates_false(self):
        assert not validate_mda_artifact_token(
            artifact_token="wrong", blueprint="bp", secret="sec"
        )

    def test_wrong_secret_validates_false(self):
        token = create_mda_artifact_token(blueprint="bp", secret="sec")
        assert not validate_mda_artifact_token(
            artifact_token=token, blueprint="bp", secret="wrong"
        )

    def test_wrong_blueprint_validates_false(self):
        token = create_mda_artifact_token(blueprint="bp1", secret="sec")
        assert not validate_mda_artifact_token(
            artifact_token=token, blueprint="bp2", secret="sec"
        )

    def test_return_debug_returns_dict(self):
        token = create_mda_artifact_token(blueprint="bp", secret="sec")
        result = validate_mda_artifact_token(
            artifact_token=token, blueprint="bp", secret="sec", return_debug=True
        )
        assert isinstance(result, dict)
        assert result["valid"] is True

    def test_dict_blueprint_validates(self):
        bp = {"id": "abc123"}
        token = create_mda_artifact_token(blueprint=bp, secret="sec")
        assert validate_mda_artifact_token(artifact_token=token, blueprint=bp, secret="sec")
