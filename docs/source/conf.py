# -*- coding: utf-8 -*-
#
# HAT  documentation build configuration file, created by

import sys, os

extensions = ['sphinx.ext.intersphinx']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'HAT: HESXRD-Analysis-Toolkit (xray-hat) documentation'
copyright = u'2022, Gary Harlow'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '2.0'
# The full version, including alpha/beta/rc tags.
release = '2.0.0'

exclude_trees = []

pygments_style = 'sphinx'

import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
   'collapse_navigation': False,

}
html_static_path = ['_static']


htmlhelp_basename = 'xray-hat-doc'


latex_documents = [
  ('index', 'xray-hat.tex', u'HAT Documentation',
   u'Gary Harlow', 'manual'),
]



# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'http://docs.python.org/': None}
language = 'en'

