.. _installation:

Installation
============

The ``pdm_utils`` package is written for Python >=3.7 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies that are only needed for a subset of ``pdm_utils`` tools. Below is a general step-by-step installation guide.

Required
________

.. toctree::
  :hidden:

  MySQL <./installation/install_mysql>
  Anaconda <./installation/install_conda>
  pdm_utils Python package <./installation/install_pdm_utils_package>
  MySQL phage genomics database instance <./installation/install_database>

The following software should be installed in the order indicated:

1. :ref:`MySQL <install_mysql>`
2. :ref:`Anaconda <install_conda>`
3. :ref:`pdm_utils Python package <install_pdm_utils_package>`
4. :ref:`MySQL phage genomics database instance <install_database>`



Optional
________

.. toctree::
  :hidden:

  MMseqs2 <./installation/install_mmseqs2>
  NCBI Blast+ <./installation/install_ncbi_tools>
  NCBI Conserved Domain Database <./installation/install_ncbi_tools>
  pdm_utils repository <./installation/install_pdm_utils_repo>

Several ``pdm_utils`` tools have specific dependencies. Install the following software as needed:


- :ref:`MMseqs2 (creating gene phamilies) <install_mmseqs2>`
- :ref:`NCBI Blast+ (finding conserved domains) <install_ncbi_tools>`
- :ref:`NCBI Conserved Domain Database (finding conserved domains) <install_ncbi_tools>`
- :ref:`pdm_utils repository (source code, MySQL scripts, tests) <install_pdm_utils_repo>`
