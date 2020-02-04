.. _installation:

Installation
============


The ``pdm_utils`` package is written in Python 3 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies that are only needed for a subset of ``pdm_utils`` tools. Below is a general step-by-step guide to installing all dependencies.

Required dependencies
_____________________

.. toctree::
   :maxdepth: 1

   MySQL <./installation/install_mysql>
   Anaconda <./installation/install_conda>
   pdm_utils Python package <./installation/install_pdm_utils_package>
   MySQL phage genomics database instance <./installation/install_database>



Optional tool-specific dependencies
___________________________________

Several ``pdm_utils`` tools have specific dependencies. Install the following tools/files as needed:

.. toctree::
   :maxdepth: 1

   MMseqs2 (creating gene phamilies) <./installation/install_mmseqs2>
   NCBI Blast+ (finding conserved domains) <./installation/install_ncbi_tools>
   NCBI Conserved Domain Database (finding conserved domains) <./installation/install_ncbi_tools>
   pdm_utils repository (source code and MySQL scripts) <./installation/install_pdm_utils_repo>
