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
import os
import sys
import textwrap
sys.path.append( "/usr/local/lib/python3.9/site-packages/breathe" )

# -- Project information -----------------------------------------------------

project = 'ABACUS-DeePKS'
copyright = '2021, x'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ 'sphinx.ext.todo', 'breathe', 'exhale', 'sphinx.ext.mathjax', "sphinx_rtd_theme", 'myst_parser'
]

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
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for breathe -------------------------------------------------

breathe_projects = { "ABACUS-DeePKS": "../../../doxygen/xml/" }
breathe_default_project = "ABACUS-DeePKS"

#-- Options for exhale-------------------------------------
# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./DeePKS_API",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "DeePKS",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    textwrap.dedent('''
		PROJECT_NAME = module_deepks
		INPUT = ../../../../module_hamilt_lcao/hamilt_lcaodft
		PROJECT_BRIEF = "DeePKS: Generate descriptors, load a model and calculate energy and force."
		FILE_PATTERNS = LCAO_descriptor*
		EXTRACT_ALL = YES
		USE_MATHJAX = YES
		GENERATE_XML = YES
		ENABLE_PREPROCESSING = NO
	''')
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'
