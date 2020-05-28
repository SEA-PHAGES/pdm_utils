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
version = '0.6.3'
# The full version, including alpha/beta/rc tags
release = '0.6.3'


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
    'anaconda': ('https://www.anaconda.com/%s', ''),
    'aragorn': ('http://130.235.244.92/ARAGORN/%s', ''),
    'biopython': ('https://biopython.org/%s', ''),
    'blastplus': ('https://blast.ncbi.nlm.nih.gov/'
        'Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download%s', ''),
    'blastlegacy': ('https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.14%s', ''),
    'cdd': ('https://www.ncbi.nlm.nih.gov/cdd/%s', ''),
    'cdd_le': ('ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz%s',''),
    'conda_docs': ('https://conda.io/projects/conda/en/latest/%s', ''),
    'configparser': ('https://docs.python.org/3/library/configparser.html/%s', ''),
    'dnamaster': ('http://cobamide2.bio.pitt.edu/%s', ''),
    'gbfff': ('https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html%s', ''),
    'git': ('https://git-scm.com/%s', ''),
    'hatfullserver': ('http://phamerator.webfactional.com/databases_Hatfull/%s',
        ''),
    'infernal': ('http://eddylab.org/infernal/%s', ''),
    'mysql': ('https://www.mysql.com/%s', ''),
    'mysql_docs': ('https://dev.mysql.com/doc/refman/5.7/en/%s', ''),
    'mmseqs2': ('https://github.com/soedinglab/MMseqs2%s', ''),
    'networkx': ('https://networkx.github.io%s', ''),
    'ncbi': ('https://www.ncbi.nlm.nih.gov/%s', ''),
    'paramiko': ('http://www.paramiko.org/%s', ''),
    'pecaan': ('https://discover.kbrinsgd.org/%s', ''),
    'pdmutils': ('https://github.com/SEA-PHAGES/pdm_utils/%s', ''),
    'phagesdb': ('http://phagesdb.org/%s', ''),
    'pymysql': ('https://pymysql.readthedocs.io/en/latest/%s', ''),
    'seaphages': ('https://seaphages.org/%s', ''),
    'sqlalchemy': ('https://www.sqlalchemy.org/%s', ''),
    'tabulate': ('https://github.com/astanin/python-tabulate/%s', ''),
    'trna_scan_se': ('http://lowelab.ucsc.edu/tRNAscan-SE/%s', ''),
    'webpham': ('https://phamerator.org/%s', '')
}
