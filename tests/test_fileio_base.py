"""
Tests for pharmacon.fileio.base.PharmaconHDF5File
"""
import pytest
import numpy as np
from pathlib import Path

from pharmacon.fileio.base import PharmaconHDF5File, PharmaconHDF5Types


# ---------------------------------------------------------------------------
# Initialization and file creation
# ---------------------------------------------------------------------------

class TestFileInit:
    def test_create_pta_file(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert path.exists()

    def test_create_psa_file(self, tmp_path):
        path = tmp_path / "test.psa"
        with PharmaconHDF5File(path, mode="w") as f:
            assert path.exists()

    def test_unsupported_extension_raises(self, tmp_path):
        with pytest.raises(ValueError, match="Unsupported file extension"):
            PharmaconHDF5File(tmp_path / "test.hdf5", mode="w")

    def test_overwrite_true_with_mode_r_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with pytest.raises(ValueError):
            PharmaconHDF5File(path, mode="r", overwrite=True)

    def test_existing_file_no_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w"):
            pass
        with pytest.raises(FileExistsError):
            PharmaconHDF5File(path, mode="a", overwrite=False)

    def test_overwrite_true_replaces_file(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("old_group")
        with PharmaconHDF5File(path, mode="a", overwrite=True) as f:
            assert not f.group_exists("old_group")

    def test_file_type_pta(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert f.file_type == PharmaconHDF5Types.TRAJECTORY_ANALYSIS

    def test_file_type_psa(self, tmp_path):
        path = tmp_path / "test.psa"
        with PharmaconHDF5File(path, mode="w") as f:
            assert f.file_type == PharmaconHDF5Types.STRUCTURE_ANALYSIS

    def test_path_is_resolved(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert f.path == path.resolve()

    def test_context_manager_closes_file(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert f.file.id.valid
        assert not f.file.id.valid

    def test_file_is_valid_when_open(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert f.file.id.valid  # h5py returns int; file is open and valid

    def test_repr_contains_path(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert "test.pta" in repr(f)


# ---------------------------------------------------------------------------
# File-level metadata
# ---------------------------------------------------------------------------

class TestFileMetadata:
    def test_pharmacon_version_written(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert "pharmacon_version" in f.file.attrs

    def test_signature_written_with_command_and_subcommand(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w", command="trajectory", subcommand="rmsd") as f:
            sig = str(f.file.attrs.get("signature", ""))
        assert sig.startswith("PHARMACON::")

    def test_signature_empty_without_command(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            sig = str(f.file.attrs.get("signature", ""))
        assert sig == ""

    def test_fingerprint_written_with_command(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w", command="trajectory", subcommand="rmsd") as f:
            fp = str(f.file.attrs.get("fingerprint", ""))
        assert fp != ""

    def test_add_file_metadata(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.add_file_metadata({"custom_key": "custom_value"}, overwrite=True)
            assert str(f.file.attrs["custom_key"]) == "custom_value"

    def test_duplicate_metadata_key_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            with pytest.raises(KeyError):
                f.add_file_metadata({"pharmacon_version": "0.0"}, overwrite=False)


# ---------------------------------------------------------------------------
# Groups
# ---------------------------------------------------------------------------

class TestGroups:
    def test_create_group(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("mygroup")
            assert f.group_exists("mygroup")

    def test_duplicate_group_without_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("mygroup")
            with pytest.raises(ValueError, match="already exists"):
                f.create_group("mygroup")

    def test_delete_group(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("mygroup")
            f.delete_group("mygroup")
            assert not f.group_exists("mygroup")

    def test_delete_nonexistent_group_is_noop(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.delete_group("ghost")  # must not raise

    def test_get_groups_returns_list(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("g1")
            f.create_group("g2")
            groups = f.get_groups()
            assert "g1" in groups
            assert "g2" in groups

    def test_group_exists_false_for_missing(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            assert not f.group_exists("nonexistent")

    def test_add_group_metadata(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("g1")
            f.add_group_metadata(group_name="g1", metadata={"key": "val"})
            assert str(f.file["g1"].attrs["key"]) == "val"

    def test_add_group_metadata_missing_group_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            with pytest.raises(KeyError):
                f.add_group_metadata(group_name="missing", metadata={"k": "v"})


# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------

class TestDatasets:
    def test_create_dataset(self, tmp_path):
        path = tmp_path / "test.pta"
        data = np.array([1.0, 2.0, 3.0])
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            f.create_dataset(group_name="grp", dataset_name="vals", data=data)
            assert "vals" in f.file["grp"]

    def test_dataset_roundtrip(self, tmp_path):
        path = tmp_path / "test.pta"
        original = np.array([1.5, 2.5, 3.5])
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            f.create_dataset(group_name="grp", dataset_name="vals", data=original)

        with PharmaconHDF5File(path, mode="r") as f:
            read_back = f.file["grp"]["vals"][:]

        np.testing.assert_array_equal(original, read_back)

    def test_dataset_missing_group_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            with pytest.raises(KeyError):
                f.create_dataset(
                    group_name="missing", dataset_name="d", data=np.array([1])
                )

    def test_duplicate_dataset_without_overwrite_raises(self, tmp_path):
        path = tmp_path / "test.pta"
        data = np.array([1, 2, 3])
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            f.create_dataset(group_name="grp", dataset_name="d", data=data)
            with pytest.raises(ValueError, match="already exists"):
                f.create_dataset(group_name="grp", dataset_name="d", data=data)

    def test_dataset_with_metadata(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            f.create_dataset(
                group_name="grp", dataset_name="d",
                data=np.array([1.0]),
                metadata={"unit": "angstrom"},
            )
            assert str(f.file["grp"]["d"].attrs["unit"]) == "angstrom"

    def test_contains_operator(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            assert "grp" in f

    def test_getitem_operator(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("grp")
            assert f["grp"] is not None

    def test_len_counts_top_level_items(self, tmp_path):
        path = tmp_path / "test.pta"
        with PharmaconHDF5File(path, mode="w") as f:
            f.create_group("g1")
            f.create_group("g2")
            assert len(f) >= 2
