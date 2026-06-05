"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.command_line.dispatcher.dispatch
"""
import sys
import pytest
from unittest.mock import patch, MagicMock

from pharmacon.command_line.dispatcher import dispatch
from pharmacon.command_line.base import CommandSpec, SubcommandSpec
from pharmacon.command_line.exceptions import ValidationError


def _make_registry(*, validate_raises=None, run_raises=None):
    """Build a minimal mock registry with one command and one subcommand."""
    mock_run = MagicMock()
    mock_validate = MagicMock()

    if validate_raises is not None:
        mock_validate.side_effect = validate_raises
    if run_raises is not None:
        mock_run.side_effect = run_raises

    def _build_parser(subparsers, parents=None):
        return subparsers.add_parser("test-sub")

    sub_spec = SubcommandSpec(
        name="test-sub",
        summary="A test subcommand",
        module_path="pharmacon.command_line.test.stub",
        run_fn=mock_run,
        build_parser_fn=_build_parser,
        validate_fn=mock_validate,
    )
    cmd_spec = CommandSpec(
        name="test-cmd",
        summary="A test command",
        subcommands={"test-sub": sub_spec},
    )
    registry = MagicMock()
    registry.commands = {"test-cmd": cmd_spec}
    return registry, mock_run, mock_validate


class TestTopLevelHelp:
    @patch("pharmacon.command_line.dispatcher.fmt.print_top_level_help")
    def test_empty_argv_shows_top_level_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch([], registry, "pharmacon")
        mock_help.assert_called_once_with("pharmacon", registry)

    @patch("pharmacon.command_line.dispatcher.fmt.print_top_level_help")
    def test_h_flag_shows_top_level_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch(["-h"], registry, "pharmacon")
        mock_help.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_top_level_help")
    def test_help_flag_shows_top_level_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch(["--help"], registry, "pharmacon")
        mock_help.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_top_level_help")
    def test_help_word_shows_top_level_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch(["help"], registry, "pharmacon")
        mock_help.assert_called_once()


class TestCommandHelp:
    @patch("pharmacon.command_line.dispatcher.fmt.print_command_help")
    def test_command_only_shows_command_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch(["test-cmd"], registry, "pharmacon")
        mock_help.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_command_help")
    def test_command_with_h_shows_command_help(self, mock_help):
        registry, _, _ = _make_registry()
        dispatch(["test-cmd", "-h"], registry, "pharmacon")
        mock_help.assert_called_once()


class TestDispatchErrors:
    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_unknown_command_exits_1(self, mock_error):
        registry, _, _ = _make_registry()
        with pytest.raises(SystemExit) as exc_info:
            dispatch(["unknown-command"], registry, "pharmacon")
        assert exc_info.value.code == 1
        mock_error.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_unknown_subcommand_exits_1(self, mock_error):
        registry, _, _ = _make_registry()
        with pytest.raises(SystemExit) as exc_info:
            dispatch(["test-cmd", "ghost-sub"], registry, "pharmacon")
        assert exc_info.value.code == 1
        mock_error.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_validation_error_exits_1(self, mock_error):
        registry, mock_run, _ = _make_registry(
            validate_raises=ValidationError("bad input")
        )
        with pytest.raises(SystemExit) as exc_info:
            dispatch(["test-cmd", "test-sub"], registry, "pharmacon")
        assert exc_info.value.code == 1
        mock_run.assert_not_called()

    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_runtime_error_in_run_exits_1(self, mock_error):
        registry, _, _ = _make_registry(run_raises=RuntimeError("crash"))
        with pytest.raises(SystemExit) as exc_info:
            dispatch(["test-cmd", "test-sub"], registry, "pharmacon")
        assert exc_info.value.code == 1

    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    @patch("pharmacon.command_line.dispatcher.fmt.console")
    def test_keyboard_interrupt_exits_130(self, mock_console, mock_error):
        registry, _, _ = _make_registry(run_raises=KeyboardInterrupt())
        with pytest.raises(SystemExit) as exc_info:
            dispatch(["test-cmd", "test-sub"], registry, "pharmacon")
        assert exc_info.value.code == 130


class TestDispatchSuccess:
    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_valid_dispatch_calls_validate_and_run(self, _):
        registry, mock_run, mock_validate = _make_registry()
        dispatch(["test-cmd", "test-sub"], registry, "pharmacon")
        mock_validate.assert_called_once()
        mock_run.assert_called_once()

    @patch("pharmacon.command_line.dispatcher.fmt.print_error")
    def test_run_not_called_when_validate_fails(self, _):
        registry, mock_run, _ = _make_registry(
            validate_raises=ValidationError("oops")
        )
        with pytest.raises(SystemExit):
            dispatch(["test-cmd", "test-sub"], registry, "pharmacon")
        mock_run.assert_not_called()
