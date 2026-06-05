"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.utils.workspace.PharmaconWorkspace
"""
import os
import pytest
from pathlib import Path

from pharmacon.utils.workspace import PharmaconWorkspace



class TestPharmaconWorkspaceInit:
    def test_no_temp_dir(self, tmp_path):
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        assert ws.temp_directory is None

    def test_working_directory_resolved(self, tmp_path):
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        assert ws.working_directory == tmp_path.resolve()

    def test_working_directory_created_if_missing(self, tmp_path):
        new_dir = tmp_path / "subdir" / "nested"
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(new_dir))
        assert new_dir.is_dir()

    def test_temp_dir_created(self, tmp_path):
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        assert ws.temp_directory is not None
        assert ws.temp_directory.is_dir()
        ws.cleanup()

    def test_temp_dir_name_prefix(self, tmp_path):
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        assert ws.temp_directory.name.startswith("pharmacon-")
        ws.cleanup()

    def test_temp_dir_inside_base_path(self, tmp_path):
        base = tmp_path / "tempbase"
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True,
            working_directory=str(tmp_path),
            temp_directory_base=str(base),
            cleanup_on_exit=False,
        )
        assert ws.temp_directory.parent == base.resolve()
        ws.cleanup()

    def test_origin_cwd_stored(self, tmp_path):
        original = Path.cwd().resolve()
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        assert ws.origin_cwd == original


class TestPharmaconWorkspaceCleanup:
    def test_cleanup_removes_temp_dir(self, tmp_path):
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        temp = ws.temp_directory
        assert temp.is_dir()
        ws.cleanup()
        assert not temp.is_dir()

    def test_cleanup_sets_temp_directory_to_none(self, tmp_path):
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        ws.cleanup()
        assert ws.temp_directory is None

    def test_cleanup_noop_when_no_temp_dir(self, tmp_path):
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        ws.cleanup()  # must not raise

    def test_double_cleanup_is_safe(self, tmp_path):
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        ws.cleanup()
        ws.cleanup()  # second call must not raise


class TestChdirAndRestore:
    def test_chdir_to_temp_changes_cwd(self, tmp_path):
        original = Path.cwd().resolve()
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        ws.chdir_to_temp()
        assert Path.cwd().resolve() == ws.temp_directory
        ws.restore_cwd()
        ws.cleanup()
        assert Path.cwd().resolve() == original

    def test_restore_cwd_returns_to_origin(self, tmp_path):
        original = Path.cwd().resolve()
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=False
        )
        ws.chdir_to_temp()
        ws.restore_cwd()
        ws.cleanup()
        assert Path.cwd().resolve() == original

    def test_chdir_without_temp_dir_raises(self, tmp_path):
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        with pytest.raises(RuntimeError):
            ws.chdir_to_temp()

    def test_chdir_to_tmp_on_init(self, tmp_path):
        original = Path.cwd().resolve()
        ws = PharmaconWorkspace(
            is_tmp_dir_needed=True,
            chdir_to_tmp=True,
            working_directory=str(tmp_path),
            cleanup_on_exit=False,
        )
        assert Path.cwd().resolve() == ws.temp_directory
        ws.restore_cwd()
        ws.cleanup()
        assert Path.cwd().resolve() == original


class TestPharmaconWorkspaceContextManager:
    def test_context_manager_restores_cwd(self, tmp_path):
        original = Path.cwd().resolve()
        with PharmaconWorkspace(
            is_tmp_dir_needed=True,
            chdir_to_tmp=True,
            working_directory=str(tmp_path),
            cleanup_on_exit=True,
        ):
            assert Path.cwd().resolve() != original
        assert Path.cwd().resolve() == original

    def test_context_manager_cleans_temp_dir(self, tmp_path):
        with PharmaconWorkspace(
            is_tmp_dir_needed=True, working_directory=str(tmp_path), cleanup_on_exit=True
        ) as ws:
            temp = ws.temp_directory
        assert temp is None or not temp.is_dir()

    def test_context_manager_returns_self(self, tmp_path):
        ws_outer = PharmaconWorkspace(
            is_tmp_dir_needed=False, working_directory=str(tmp_path)
        )
        with ws_outer as ws_inner:
            assert ws_inner is ws_outer

    def test_repr_contains_class_name(self, tmp_path):
        ws = PharmaconWorkspace(is_tmp_dir_needed=False, working_directory=str(tmp_path))
        assert "PharmaconWorkspace" in repr(ws)
