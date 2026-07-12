# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'ABACUS'
copyright = '2024, ABACUS'
author = 'ABACUS'

# The full version, including alpha/beta/rc tags
# release = '2.3.5'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'myst_parser',
        'deepmodeling_sphinx',
]
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 4

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['README.md']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_logo = 'abacus-logo.svg'

# Theme options for sphinx-book-theme
html_theme_options = {
    "show_toc_level": 2,  # Only show h2 (categories) in right sidebar, not h3 (parameters)
    "toc_title": "On this page",
}


# Changes for compatibility with Read the Docs
import os

# Define the canonical URL if you are using a custom domain on Read the Docs
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

latex_engine = 'xelatex'
mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.min.js'
# deepmodeling_current_site = 'Tutorials'
latex_elements = {
    'extraclassoptions':'openany,oneside'
}


# -- Auto-generate INPUT keyword documentation from ABACUS metadata ----------

import shutil
import subprocess
import sys
from pathlib import Path

try:
    from sphinx.util import logging as sphinx_logging
    logger = sphinx_logging.getLogger(__name__)
except Exception:
    logger = None


def warn_input_docs(message):
    if logger is not None:
        logger.warning(message)
    else:
        print(f"Warning: {message}")


def input_docs_refresh_required(env=None):
    env = os.environ if env is None else env
    return env.get('READTHEDOCS') == 'True'


def handle_input_docs_failure(message, warn, fail_on_error):
    if fail_on_error:
        raise RuntimeError(message)
    warn(message)


def candidate_abacus_binaries(repo_root):
    """Return candidate ABACUS binaries for local and Read the Docs builds."""
    env_binary = os.environ.get('ABACUS_BINARY') or os.environ.get('ABACUS_EXECUTABLE')
    if env_binary:
        return [Path(env_binary)]
    candidates = []
    candidates.extend([
        repo_root / 'build-rtd-docs' / 'source' / 'abacus_pw_ser',
        repo_root / 'build-rtd-docs' / 'abacus_pw_ser',
    ])
    for binary_name in ['abacus', 'abacus_pw_ser']:
        path_binary = shutil.which(binary_name)
        if path_binary:
            candidates.append(Path(path_binary))
    return candidates


def find_abacus_binary(repo_root):
    for candidate in candidate_abacus_binaries(repo_root):
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
    return None


def refresh_input_docs(docs_dir, abacus_binary=None, warn=warn_input_docs, fail_on_error=False):
    """Refresh input-main.md from an ABACUS executable, if one is available."""
    docs_dir = Path(docs_dir)
    repo_root = docs_dir.parent
    binary = Path(abacus_binary) if abacus_binary else find_abacus_binary(repo_root)
    if binary is None or not binary.exists() or not os.access(binary, os.X_OK):
        handle_input_docs_failure(
            "ABACUS executable not found; INPUT parameter documentation may not "
            "be up to date. Set ABACUS_BINARY=/path/to/abacus or build the "
            "reduced documentation binary first.",
            warn,
            fail_on_error,
        )
        return False

    env = os.environ.copy()
    env.setdefault('OMP_NUM_THREADS', '1')
    result = subprocess.run(
        [str(binary), '--generate-parameters-yaml'],
        cwd=repo_root,
        env=env,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        handle_input_docs_failure(
            "ABACUS could not generate INPUT parameter metadata; documentation "
            f"may not be up to date. Command: {binary} --generate-parameters-yaml"
            + (f". Output: {detail}" if detail else ""),
            warn,
            fail_on_error,
        )
        return False

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import yaml
    from generate_input_main import generate_from_data
    generate_from_data(
        data=yaml.safe_load(result.stdout),
        output=docs_dir / 'advanced' / 'input_files' / 'input-main.md',
    )
    return True


def generate_input_docs(app):
    refresh_input_docs(
        Path(__file__).resolve().parent,
        fail_on_error=input_docs_refresh_required(),
    )

def setup(app):
    app.connect('builder-inited', generate_input_docs)
