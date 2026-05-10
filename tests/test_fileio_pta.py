"""
Tests for pharmacon.fileio.pta.PharmaconPTAFile
"""
import json
import pytest
import numpy as np

from pharmacon.fileio.pta import PharmaconPTAFile


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open(path, mode="a"):
    return PharmaconPTAFile(path, mode=mode)


# ---------------------------------------------------------------------------
# Schema constants
# ---------------------------------------------------------------------------

class TestSchemaConstants:
    def test_interaction_columns_length(self):
        assert len(PharmaconPTAFile.INTERACTION_COLUMNS) == 21

    def test_interaction_columns_first(self):
        assert PharmaconPTAFile.INTERACTION_COLUMNS[0] == "frame_number"

    def test_interaction_columns_last(self):
        assert PharmaconPTAFile.INTERACTION_COLUMNS[-1] == "details"

    def test_detail_schema_has_nine_interaction_types(self):
        assert len(PharmaconPTAFile.INTERACTION_DETAIL_SCHEMA) == 9

    def test_detail_schema_contains_all_types(self):
        expected = {
            "HYDROPHOBIC", "HYDROGEN-BOND", "IONIC", "HALOGEN-BOND",
            "METAL-CONTACT", "WATER-BRIDGE-1", "WATER-BRIDGE-2",
            "PI-CATION", "PI-STACKING",
        }
        assert set(PharmaconPTAFile.INTERACTION_DETAIL_SCHEMA.keys()) == expected

    def test_hydrophobic_detail_fields(self):
        fields = PharmaconPTAFile.INTERACTION_DETAIL_SCHEMA["HYDROPHOBIC"]
        assert "distance" in fields
        assert "is_hydrogen" in fields

    def test_water_bridge_1_has_water_fields(self):
        fields = PharmaconPTAFile.INTERACTION_DETAIL_SCHEMA["WATER-BRIDGE-1"]
        assert "water_index" in fields
        assert "d_g1_water" in fields

    def test_water_bridge_2_has_two_water_sets(self):
        fields = PharmaconPTAFile.INTERACTION_DETAIL_SCHEMA["WATER-BRIDGE-2"]
        assert "water1_index" in fields
        assert "water2_index" in fields


# ---------------------------------------------------------------------------
# write_frame_interactions / read_frame_interactions
# ---------------------------------------------------------------------------

class TestFrameInteractionsRoundtrip:
    def test_write_and_read_empty_interactions(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[], group_name="pl_interactions")
            result = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
        assert result == ()

    def test_write_and_read_single_record(self, tmp_path):
        path = tmp_path / "test.pta"
        rec = ["HYDROPHOBIC", 1, "CA", 2, "C", "BB", "ALA", 1, "A", "PROA",
               10, "CB", 11, "C", "SC", "VAL", 2, "A", "PROA", 3.5, "false"]
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[rec], group_name="pl_interactions")
            result = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
        assert len(result) == 1
        assert result[0][0] == "HYDROPHOBIC"

    def test_write_multiple_records(self, tmp_path):
        path = tmp_path / "test.pta"
        recs = [["IONIC", i] for i in range(5)]
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=recs, group_name="pl_interactions")
            result = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
        assert len(result) == 5

    def test_multiple_frames_stored_independently(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[["A"]], group_name="pl_interactions")
            f.write_frame_interactions(frame_index=1, interactions=[["B"], ["C"]], group_name="pl_interactions")
            r0 = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
            r1 = f.read_frame_interactions(frame_index=1, group_name="pl_interactions")
        assert len(r0) == 1
        assert len(r1) == 2

    def test_roundtrip_preserves_nested_data(self, tmp_path):
        path = tmp_path / "test.pta"
        rec = {"type": "HYDROPHOBIC", "distance": 3.5, "nested": [1, 2, 3]}
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[rec], group_name="pl_interactions")
            result = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
        assert result[0] == rec


# ---------------------------------------------------------------------------
# write_frame_interactions — error paths
# ---------------------------------------------------------------------------

class TestFrameInteractionsErrors:
    def test_negative_frame_index_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            with pytest.raises(ValueError, match="frame_index must be a non-negative integer"):
                f.write_frame_interactions(frame_index=-1, interactions=[], group_name="pl_interactions")

    def test_float_frame_index_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            with pytest.raises(ValueError):
                f.write_frame_interactions(frame_index=1.5, interactions=[], group_name="pl_interactions")  # type: ignore

    def test_duplicate_frame_without_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[], group_name="pl_interactions")
            with pytest.raises(FileExistsError):
                f.write_frame_interactions(frame_index=0, interactions=[], group_name="pl_interactions", overwrite=False)

    def test_duplicate_frame_with_overwrite_succeeds(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[["OLD"]], group_name="pl_interactions")
            f.write_frame_interactions(frame_index=0, interactions=[["NEW"]], group_name="pl_interactions", overwrite=True)
            result = f.read_frame_interactions(frame_index=0, group_name="pl_interactions")
        assert result[0][0] == "NEW"

    def test_read_nonexistent_frame_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            with pytest.raises(KeyError):
                f.read_frame_interactions(frame_index=99, group_name="pl_interactions")

    def test_read_negative_frame_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with _open(path, mode="w") as f:
            f.create_group("pl_interactions")
            with pytest.raises(ValueError):
                f.read_frame_interactions(frame_index=-1, group_name="pl_interactions")
