"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Compact section-header helpers for terminal and file loggers.

Both ``header()`` and ``subheader()`` emit through the ``pharmacon`` root
logger at ``INFO`` level so that terminal and file handlers apply their own
level filtering consistently.

Terminal rendering uses Rich renderables attached to the log record; the
``PharmaconRichHandler`` prints them directly via its own console so they
inherit the correct terminal width.  File output is plain ASCII so log files
remain grep-friendly.

Two levels of visual weight are provided:

``header(title)``
    Major section separator — full-width rule with the title centred.

``subheader(title)``
    Minor sub-section marker — short leading rule, title, no trailing rule.

Usage::

    from pharmacon.logger.headers import header, subheader

    header("Protein Preparation")
    subheader("Fixing missing atoms")
"""

from __future__ import annotations

import logging

from rich.rule import Rule
from rich.text import Text



__all__ = [
    "header",
    "subheader",
]


_RULE_CHAR = "─"
_FILE_RULE_CHAR = "-"
_FILE_WIDTH = 72


# Internal helpers


def _file_rule(title: str, char: str = _FILE_RULE_CHAR, width: int = _FILE_WIDTH) -> str:
    """
    Generates a formatted string based on the provided title, character, and width.
    The function creates a rule line including the title, padded with the specified
    character, ensuring the total width is met. If no title is provided, a line
    consisting solely of the character is returned.

    Example::

        ---- Protein Preparation --------------------------------------------------

    :param title: The title to include in the rule.
    :type title: str
    :param char: The character to be used for padding the rule. Defaults to '_FILE_RULE_CHAR'.
    :type char: str
    :param width: The total width of the rule. Defaults to '_FILE_WIDTH'.
    :type width: int
    :return: A generated string that conforms to the specified inputs.
    :rtype: str
    """
    if not title:
        return char * width
    inner = f" {title} "
    left = char * 4
    right = char * max(0, width - len(inner) - 4)
    return f"{left}{inner}{right}"


# Public API


def header(title: str,
           *,
           style: str = "bold #ff8700",
           width: int | None = None) -> None:
    """
    Logs a formatted header with a styled rule for terminal output. The method
    creates a logging record with a custom terminal-renderable rule for visually
    distinct headers and processes it through the logging system. The styling,
    alignment, and width of the rule can be customized using the parameters.

    Terminal — full-width Rich ``Rule`` with centred *title*.
    File     — plain ASCII rule line.

    :param title: The text to be displayed within the rule.
    :type title: str
    :param style: An optional style configuration for the rule, specified in a
        string format. Defaults to "bold #ff8700".
    :type style: str, optional
    :param width: Optional integer specifying the width of the rule. If None,
        the default width is determined by the terminal width.
    :type width: int | None, optional
    :return: None
    """
    log = logging.getLogger("pharmacon")

    rule_kwargs: dict = {"style": style, "align": "center"}
    if width is not None:
        rule_kwargs["width"] = width

    record = log.makeRecord(
        name="pharmacon",
        level=logging.INFO,
        fn="",
        lno=0,
        msg=_file_rule(title),
        args=(),
        exc_info=None,
    )
    record.is_pharmacon_header = True  # type: ignore[attr-defined]
    record.terminal_renderable = Rule(title, **rule_kwargs)  # type: ignore[attr-defined]

    log.handle(record)


def subheader(title: str,
              *,
              style: str = "#00af87") -> None:
    """
    Constructs and logs a styled subheader text to a specified logger.

    This function creates a styled label with the provided title and style. It then logs
    the labeled text as a formatted message to the "pharmacon" logger. The logged message
    is primarily used for structured logging purposes and clear terminal output.

    Terminal — short leading rule + title in *style*.
    File     — ``---- title`` (no trailing rule).


    :param title: The text content to be displayed in the subheader.
    :type title: str
    :param style: The style string representing the color and/or font styling to apply to
        the title. By default, the color of the style is "#00af87".
    :type style: str, optional
    :return: Does not return any value. The output is logged to the logger.
    :rtype: None
    """
    log = logging.getLogger("pharmacon")

    label = Text()
    label.append(f"{_RULE_CHAR * 3} ", style="#808080")
    label.append(title, style=style)

    record = log.makeRecord(
        name="pharmacon",
        level=logging.INFO,
        fn="",
        lno=0,
        msg=f"---- {title}",
        args=(),
        exc_info=None,
    )
    record.is_pharmacon_header = True  # type: ignore[attr-defined]
    record.terminal_renderable = label  # type: ignore[attr-defined]

    log.handle(record)
