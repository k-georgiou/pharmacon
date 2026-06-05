"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.logger — setup_logger, get_logger, TRACE, PharmaconLogger.
"""
import logging
import pytest

from pharmacon.logger import setup_logger, get_logger, TRACE
from pharmacon.logger.levels import PharmaconLogger, register


class TestTraceLevelConstant:
    def test_trace_is_5(self):
        assert TRACE == 5

    def test_trace_below_debug(self):
        assert TRACE < logging.DEBUG

    def test_trace_name_registered(self):
        assert logging.getLevelName(TRACE) == "TRACE"

    def test_debug_level_name_registered(self):
        assert logging.getLevelName(logging.DEBUG) == "DEBUG"


class TestPharmaconLoggerClass:
    def test_is_subclass_of_logger(self):
        assert issubclass(PharmaconLogger, logging.Logger)

    def test_has_trace_method(self):
        assert hasattr(PharmaconLogger, "trace")
        assert callable(PharmaconLogger.trace)

    def test_register_sets_logger_class(self):
        register()
        assert logging.getLoggerClass() is PharmaconLogger


class TestGetLogger:
    def test_returns_pharmacon_logger(self):
        log = get_logger("pharmacon.test.module")
        assert isinstance(log, PharmaconLogger)

    def test_none_returns_root_logger(self):
        log = get_logger(None)
        assert log.name == "pharmacon"

    def test_unqualified_name_prefixed(self):
        log = get_logger("mymodule")
        assert log.name.startswith("pharmacon.")

    def test_already_qualified_name_unchanged(self):
        log = get_logger("pharmacon.analyzer.rmsd")
        assert log.name == "pharmacon.analyzer.rmsd"

    def test_different_names_return_different_loggers(self):
        log1 = get_logger("pharmacon.a")
        log2 = get_logger("pharmacon.b")
        assert log1 is not log2

    def test_same_name_returns_same_logger(self):
        log1 = get_logger("pharmacon.singleton")
        log2 = get_logger("pharmacon.singleton")
        assert log1 is log2


class TestSetupLoggerTerminal:
    def test_returns_pharmacon_logger(self):
        log = setup_logger(terminal=True, terminal_level="info", replace=True)
        assert isinstance(log, PharmaconLogger)

    def test_terminal_handler_added(self):
        log = setup_logger(terminal=True, terminal_level="info", replace=True)
        assert len(log.handlers) >= 1

    def test_no_terminal_no_handlers(self):
        log = setup_logger(terminal=False, replace=True)
        assert len(log.handlers) == 0

    def test_replace_clears_previous_handlers(self):
        setup_logger(terminal=True, replace=True)
        setup_logger(terminal=True, replace=True)
        log = get_logger(None)
        # replace=True means only one terminal handler exists
        assert len(log.handlers) == 1

    def test_propagate_is_false(self):
        log = setup_logger(terminal=True, replace=True)
        assert log.propagate is False

    @pytest.mark.parametrize("level_str,expected", [
        ("trace",    TRACE),
        ("debug",    logging.DEBUG),
        ("info",     logging.INFO),
        ("warning",  logging.WARNING),
        ("warn",     logging.WARNING),
        ("error",    logging.ERROR),
        ("critical", logging.CRITICAL),
    ])
    def test_level_strings_parsed(self, level_str, expected):
        log = setup_logger(terminal=True, terminal_level=level_str, replace=True)
        assert log.level == expected

    def test_integer_level_accepted(self):
        log = setup_logger(terminal=True, terminal_level=logging.DEBUG, replace=True)
        assert log.level == logging.DEBUG

    def test_unknown_level_string_raises(self):
        with pytest.raises(ValueError, match="Unknown log level"):
            setup_logger(terminal=True, terminal_level="verbose", replace=True)


class TestSetupLoggerFile:
    def test_file_true_without_path_raises(self):
        with pytest.raises(ValueError, match="log_file"):
            setup_logger(terminal=False, file=True, log_file=None, replace=True)

    def test_file_handler_creates_log_file(self, tmp_path):
        log_file = tmp_path / "test.log"
        setup_logger(terminal=False, file=True, log_file=log_file, replace=True)
        root = get_logger(None)
        root.info("hello from test")
        for h in root.handlers:
            h.flush()
        assert log_file.exists()

    def test_log_file_contains_message(self, tmp_path):
        log_file = tmp_path / "test.log"
        setup_logger(
            terminal=False, file=True, file_level="debug",
            log_file=log_file, replace=True,
        )
        root = get_logger(None)
        root.debug("unique_marker_xyz")
        for h in root.handlers:
            h.flush()
        assert "unique_marker_xyz" in log_file.read_text()


class TestTraceMethod:
    def test_trace_called_when_level_trace(self, tmp_path):
        log_file = tmp_path / "trace.log"
        setup_logger(
            terminal=False, file=True, file_level="trace",
            log_file=log_file, replace=True,
        )
        log = get_logger("pharmacon.trace_test")
        log.trace("trace_message_xyz")
        root = get_logger(None)
        for h in root.handlers:
            h.flush()
        assert "trace_message_xyz" in log_file.read_text()

    def test_trace_suppressed_when_level_debug(self, tmp_path):
        log_file = tmp_path / "no_trace.log"
        setup_logger(
            terminal=False, file=True, file_level="debug",
            log_file=log_file, replace=True,
        )
        log = get_logger("pharmacon.no_trace_test")
        log.trace("should_not_appear")
        root = get_logger(None)
        for h in root.handlers:
            h.flush()
        content = log_file.read_text() if log_file.exists() else ""
        assert "should_not_appear" not in content
