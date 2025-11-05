# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MOSAIK Documentation'
copyright = '2025, Diane Cruiziat, Anthony Baptista'
author = 'Diane Cruiziat, Anthony Baptista'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",     # API docs from docstrings
    "sphinx.ext.napoleon",    # Support for Google/NumPy style docstrings
    "sphinx.ext.viewcode",    # Add links to highlighted source code
    "sphinx.ext.todo",        # TODO directives in docs
    "myst_parser",            # (optional) allow Markdown files
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_logo = '_static/logo.004.png'
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']
html_title = ''
html_theme_options = {
    "logo_only": True,   # show only the logo in the header
}

html_static_path = ["_static"]
html_css_files = [
    "custom.css",
]

import os
import sys

# Add your package to the Python path
sys.path.insert(0, os.path.abspath('..'))  # Adjust based on your structure


# Napoleon settings for Google-style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Theme
html_theme = 'sphinx_rtd_theme'
