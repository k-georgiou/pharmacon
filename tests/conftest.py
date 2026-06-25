"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Shared pytest fixtures.
"""

from __future__ import annotations

import logging

import pytest


@pytest.fixture(autouse=True)
def _restore_pharmacon_logger_state():
    """Isolate global logging state between tests.

    Command ``run()`` paths call ``setup_logger``, which mutates the global
    ``pharmacon`` logger — it sets ``propagate = False`` and attaches handlers.
    Without isolation that configuration leaks into later tests (e.g. breaking
    ``caplog``-based assertions, which rely on propagation). Snapshot and
    restore the logger's propagate flag, handlers and level around every test.
    """
    plog = logging.getLogger("pharmacon")
    saved_propagate = plog.propagate
    saved_handlers = list(plog.handlers)
    saved_level = plog.level
    try:
        yield
    finally:
        plog.propagate = saved_propagate
        plog.handlers[:] = saved_handlers
        plog.setLevel(saved_level)
