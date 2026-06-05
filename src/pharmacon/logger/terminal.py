"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Terminal handler — Rich-powered colourised logging to stderr.

Design
------
* Uses ``rich.logging.RichHandler`` for timestamp / level / path decorations.
* A custom ``TerminalFormatter`` prepends the ``>> `` prefix to every message.
* Level colours are driven by a Rich ``Theme`` built from the project palette.
* The handler is completely independent of the file handler; its level may be
  set separately after construction.

Usage::

    from pharmacon.logger.terminal import make_terminal_handler

    handler = make_terminal_handler(level=logging.INFO)
    logging.getLogger("pharmacon").addHandler(handler)
"""

from __future__ import annotations

import logging

from rich.console import Console
from rich.logging import RichHandler
from rich.theme import Theme

from pharmacon.logger.levels import TRACE

# ── Rich theme — derived from the project palette ──────────────────────────────
# Level label colours follow the palette anchors:
#   #808080  grey50       → trace (dim, least important)
#   #008080  cyan         → debug
#   #00af87  dark_cyan    → info
#   #ff8700  dark_orange  → warning
#   red                   → error / critical (kept standard)



__all__ = [
    "TerminalFormatter",
    "PharmaconRichHandler",
    "make_terminal_handler",
]


_LOG_THEME = Theme(
    {
        f"logging.level.{logging.getLevelName(TRACE).lower()}": "dim #808080",
        "logging.level.debug":    "#008080",
        "logging.level.info":     "#00af87",
        "logging.level.warning":  "bold #ff8700",
        "logging.level.error":    "bold red",
        "logging.level.critical": "bold red underline",
        # chrome
        "log.time":    "#808080",
        "log.message": "default",
        "log.path":    "#808080",
    },
    inherit=True,
)

_PREFIX = ">> "


# ── Formatter ─────────────────────────────────────────────────────────────────


class TerminalFormatter(logging.Formatter):
    """Formatter for RichHandler: prepends ``>> `` to every message.

    RichHandler calls ``self.format(record)`` and uses the returned string as
    the *message* portion of the rendered log line (timestamp, level label,
    and source path are added by Rich independently).
    """

    def format(self, record: logging.LogRecord) -> str:
        record.message = record.getMessage()
        return f"{_PREFIX}{record.message}"


class PharmaconRichHandler(RichHandler):
    """RichHandler subclass that intercepts header records.

    When a log record carries a ``terminal_renderable`` attribute (set by the
    header/subheader helpers), the Rich renderable is printed directly through
    the handler's own console.  This ensures the correct terminal width is
    used for rules while still respecting handler-level filtering.
    """

    def emit(self, record: logging.LogRecord) -> None:
        renderable = getattr(record, "terminal_renderable", None)
        if renderable is not None:
            self.console.print(renderable, highlight=False)
        else:
            super().emit(record)


# ── Handler factory ────────────────────────────────────────────────────────────


def make_terminal_handler(
    level: int = logging.INFO,
    *,
    show_path: bool = True,
    show_time: bool = True,
    time_format: str = "[%H:%M:%S]",
    markup: bool = False,
    rich_tracebacks: bool = True,
    console: Console | None = None,
) -> RichHandler:
    """Return a configured :class:`RichHandler` for terminal output.

    Parameters
    ----------
    level:
        Minimum level for this handler (independent of the file handler).
    show_path:
        Display the source file and line number on the right.
    show_time:
        Display the wall-clock timestamp on the left.
    time_format:
        ``strftime`` format for the timestamp column.
    markup:
        Enable Rich markup interpretation in log messages.
    rich_tracebacks:
        Render exception tracebacks through Rich.
    console:
        Custom Rich :class:`Console` to write to.  Defaults to a new
        ``Console(stderr=True)`` with the project theme applied.
    """
    if console is None:
        console = Console(stderr=True, highlight=False, theme=_LOG_THEME)

    handler = PharmaconRichHandler(
        console=console,
        show_time=show_time,
        show_level=True,
        show_path=show_path,
        rich_tracebacks=rich_tracebacks,
        markup=markup,
        log_time_format=time_format,
        omit_repeated_times=False,
    )
    handler.setLevel(level)
    handler.setFormatter(TerminalFormatter())
    return handler
