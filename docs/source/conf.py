# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re
from pathlib import Path


def _convert_gfm_alerts(md: str) -> str:
    """Convert GitHub-style ``> [!NOTE]`` alerts into MyST ``:::{note}``
    admonitions, so a GitHub-flavoured README renders with styled callouts
    when included into the Sphinx docs."""
    kinds = {"note", "tip", "important", "warning", "caution"}
    lines = md.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        m = re.match(r"^>\s*\[!(\w+)\]\s*$", lines[i])
        if m and m.group(1).lower() in kinds:
            out.append(f":::{{{m.group(1).lower()}}}")
            i += 1
            while i < len(lines) and lines[i].lstrip().startswith(">"):
                out.append(re.sub(r"^>\s?", "", lines[i]))
                i += 1
            out.append(":::")
        else:
            out.append(lines[i])
            i += 1
    return "\n".join(out)


# Regenerate the included plot-INI reference from the canonical README so the
# docs site and the GitHub README never drift (single source of truth).
_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
_GEN = _HERE / "_generated"
_GEN.mkdir(exist_ok=True)
(_GEN / "plot_ini_readme.md").write_text(
    _convert_gfm_alerts((_REPO / "examples" / "plot_ini" / "README.md").read_text()),
    encoding="utf-8",
)


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pharmacon'
copyright = '2026, Kyriakos Georgiou'
author = 'Kyriakos Georgiou'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser']

# Markdown (MyST) support so docs pages can ``.. include::`` repo READMEs
# directly — single source of truth, no duplicated content.
myst_enable_extensions = [
    'colon_fence',
    'deflist',
    'html_admonition',
    'html_image',
]
# README files use GitHub-style ``> [!NOTE]`` alerts, which MyST renders as
# blockquotes; convert them to MyST admonitions on the fly.
myst_heading_anchors = 3

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/pharmacon_logo.png'

html_theme_options = {
    'navigation_depth': 4,
    'titles_only': False,
    'logo_only': False,
}
