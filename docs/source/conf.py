# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))


# -- Project information -----------------------------------------------------

project = 'pdm_utils'
copyright = '2019, Travis Mavrich'
author = 'Travis Mavrich'
version = '0.0.18'
# The full version, including alpha/beta/rc tags
release = '0.0.18'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.extlinks']

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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
html_static_path = []


# Explicitly define the master file to improve readthedoc builds.
master_doc = "index"


# External link library
extlinks = {
    'phagesdb': ('http://phagesdb.org/%s', ''),
    'seaphages': ('https://seaphages.org/%s', ''),
    'webpham': ('https://phamerator.org/%s', ''),
    'gbfff': ('https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html%s', ''),
    'dnamaster': ('http://cobamide2.bio.pitt.edu/%s', ''),
    'pecaan': ('https://discover.kbrinsgd.org/%s', ''),
    'biopython': ('https://biopython.org/%s', ''),
    'git': ('https://git-scm.com/%s', ''),
    'pdmutils': ('https://github.com/SEA-PHAGES/pdm_utils/%s', ''),
    'cdd': ('https://www.ncbi.nlm.nih.gov/cdd/%s', ''),
    'cdd_le': ('ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz%s',''),
    'blastplus': ('https://blast.ncbi.nlm.nih.gov/'
        'Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download%s', ''),
    'hatfullserver': ('http://phamerator.webfactional.com/databases_Hatfull/%s',
        ''),
    'pymysql': ('https://pymysql.readthedocs.io/en/latest/%s', ''),
    'paramiko': ('http://www.paramiko.org/%s', ''),
    'mysql': ('https://www.mysql.com/%s', ''),
    'mmseqs': ('https://www.ncbi.nlm.nih.gov/pubmed/26743509%s', ''),
    'anaconda': ('https://www.anaconda.com/%s', ''),
    'ncbi': ('https://www.ncbi.nlm.nih.gov/%s', '')
}
