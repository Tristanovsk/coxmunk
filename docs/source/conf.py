# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

# Make the package importable without installation (local builds);
# on Read the Docs the package is also pip-installed (see .readthedocs.yaml).
DOCS_SOURCE = Path(__file__).resolve().parent
REPO_ROOT = DOCS_SOURCE.parents[1]
sys.path.insert(0, str(REPO_ROOT))

import coxmunk  # noqa: E402

on_rtd = os.environ.get('READTHEDOCS') == 'True'

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'coxmunk'
copyright = '2019-2026, Tristan Harmel'
author = 'Tristan Harmel'
release = coxmunk.__version__
version = '.'.join(release.split('.')[:2])
today_fmt = '%Y-%m-%d'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    # converts the README GIF (first frame) for the PDF build; needs ImageMagick
    'sphinx.ext.imgconverter',
    'sphinx_copybutton',
    'sphinxcontrib.mermaid',
    'myst_nb',
    'IPython.sphinxext.ipython_console_highlighting',
]

templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '_readme.md', '**.ipynb_checkpoints', 'Thumbs.db', '.DS_Store']

# -- Autodoc / autosummary ---------------------------------------------------

autosummary_generate = True
autoclass_content = 'class'
autodoc_typehints = 'description'
# 'members' is set in the autosummary templates (_templates/) to avoid
# documenting objects twice
autodoc_default_options = {
    'member-order': 'bysource',
    'show-inheritance': True,
}

# NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

todo_include_todos = True

# the README included in index.rst starts at "##" once its title is skipped
suppress_warnings = ['myst.header']


def _cli_usage_as_literal(app, what, name, obj, options, lines):
    '''Render the docopt usage of coxmunk.visu (not reST) as a literal block.'''
    if what == 'module' and name == 'coxmunk.visu':
        usage = ['   ' + line if line else '' for line in lines]
        lines[:] = ['Command line interface and plotting helpers.', '',
                    '.. code-block:: text', ''] + usage


def setup(app):
    app.connect('autodoc-process-docstring', _cli_usage_as_literal)

# -- Math --------------------------------------------------------------------

# number labelled equations and refer to them as "Eq. (n)" with :eq:
math_eqref_format = 'Eq. ({number})'
math_numfig = True
numfig = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'xarray': ('https://docs.xarray.dev/en/stable', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
pygments_style = 'sphinx'

html_theme_options = {
    'repository_url': 'https://github.com/Tristanovsk/coxmunk',
    'repository_branch': 'master',
    'path_to_docs': 'docs/source',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_edit_page_button': True,
    'use_download_button': True,
    'navigation_with_keys': True,
    'show_toc_level': 2,
    'secondary_sidebar_items': ['page-toc', 'edit-this-page'],
}

html_title = f'coxmunk {release}'

html_show_sourcelink = False
html_last_updated_fmt = today_fmt

htmlhelp_basename = 'coxmunk_doc'

# -- Options for LaTeX/PDF output --------------------------------------------

# xelatex handles the Greek characters of the Kumatage chapter; the DejaVu
# fonts cover them and are available on most Linux systems (incl. Read the Docs)
latex_engine = 'xelatex'
# makeindex is more widely available than xindy (default with xelatex)
latex_use_xindy = False
latex_elements = {
    'fontpkg': r'''
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
''',
}

# -- MyST / notebook rendering -----------------------------------------------

# GitHub-style ```math blocks (used in README.md) rendered as display math
myst_fence_as_directive = ['math']

myst_enable_extensions = [
    'amsmath',
    'colon_fence',
    'deflist',
    'dollarmath',
    'html_admonition',
    'html_image',
    'linkify',
    'replacements',
    'smartquotes',
    'substitution',
]

# Notebooks are executed at build time (they are fast and need no input data);
# outputs are cached in docs/build/.jupyter_cache between local builds.
nb_execution_mode = 'cache'
nb_execution_cache_path = str(DOCS_SOURCE.parent / 'build' / '.jupyter_cache')
nb_execution_timeout = 600
nb_execution_raise_on_error = True
nb_merge_streams = True
