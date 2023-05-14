# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os.path
from datetime import datetime
import re

project = 'SingletScalar-DM'
copyright = f'{datetime.now().year}, Mattia Di Mauro'
author = 'Mattia Di Mauro'

package_root_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')
# load version number from setup.cfg
with open(os.path.join(package_root_path, 'setup.cfg'), 'r') as f:
    fstring = str(f.read())
    version = re.search(r'version\s*=\s*(?P<version>\d+\.\d+\.\d+)', fstring).group(1)

release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    # 'sphinx_autodoc_typehints',
    'sphinx_codeautolink',
    'sphinx_gallery.gen_gallery'
]

# Napoleon settings
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True

# Sphinx Gallery
sphinx_gallery_conf = {
    'examples_dirs': os.path.join(package_root_path, 'examples'),
    'gallery_dirs': 'examples_gallery',
    'download_all_examples': False,
    'remove_config_comments': True,
}

exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
