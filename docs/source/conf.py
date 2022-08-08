import os
from datetime import date
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Dockey'
copyright = '{}, Lianming Du'.format(date.today().year)
author = 'Lianming Du'

with open(os.path.join(os.path.abspath('../..'), 'src', 'config.py')) as fh:
	for line in fh:
		if line.startswith('DOCKEY_VERSION'):
			__version__ = line.strip().split('=')[1].strip().strip('"')
			break

version = __version__
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
