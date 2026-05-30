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
exclude_patterns = []


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


# -- Auto-generate INPUT keyword documentation from YAML parameter dump ------

from pathlib import Path

def generate_input_docs(app):
    """Auto-generate input-main.md from parameters.yaml before building.

    Workflow:
        abacus --generate-parameters-yaml > docs/parameters.yaml
        # Then Sphinx calls this hook, which runs generate_input_main.py
    """
    docs_dir = Path(__file__).resolve().parent
    yaml_path = docs_dir / 'parameters.yaml'
    if not yaml_path.exists():
        print(f"Warning: {yaml_path} not found. "
              "Run: abacus --generate-parameters-yaml > docs/parameters.yaml")
        return
    import sys
    sys.path.insert(0, str(docs_dir))
    from generate_input_main import generate
    generate(
        yaml_path=yaml_path,
        output=docs_dir / 'advanced' / 'input_files' / 'input-main.md',
    )

def setup(app):
    app.connect('builder-inited', generate_input_docs)
