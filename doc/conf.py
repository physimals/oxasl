# -*- coding: utf-8 -*-
#
# oxasl documentation build configuration file
#
# This file is execfile()d with the current directory set to its
# containing dir.

# -- General configuration ------------------------------------------------

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode'
]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'oxasl'
copyright = u'2019, University of Oxford'
author = u'Martin Craig'
build_dir = u"_build"

version = u'0.0.1'
release = u'0.0.1'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

# -- Options for HTML output ----------------------------------------------

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
}

latex_documents = [
    (master_doc, 'oxasl.tex', u'oxasl Documentation',
     u'Martin Craig', 'manual'),
]
