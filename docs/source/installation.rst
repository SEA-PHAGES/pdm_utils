.. _installation:

Installation
============

The ``pdm_utils`` package is written for Python >=3.7 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies that are only needed for a subset of ``pdm_utils`` tools. Below is a general step-by-step installation guide.

Required
________

.. toctree::
  :hidden:

  ./installation/install_mysql
  ./installation/install_conda
  ./installation/install_pdm_utils_package
  ./installation/install_database

The following software should be installed in the order indicated:

1. :ref:`MySQL <install_mysql>`
2. :ref:`Anaconda <install_conda>`
3. :ref:`pdm_utils Python package <install_pdm_utils_package>`
4. :ref:`MySQL phage genomics database instance <install_database>`



Optional
________

.. toctree::
  :hidden:

  ./installation/install_mmseqs2
  ./installation/install_ncbi_tools
  ./installation/install_mcl
  ./installation/install_pdm_utils_repo
  ./installation/install_trnascan_se
  ./installation/install_aragorn
  ./installation/install_infernal


Several ``pdm_utils`` tools have specific dependencies. Install the following software as needed:


- :ref:`MMseqs2 (creating gene phamilies) <install_mmseqs2>`
- :ref:`NCBI Blast+ (finding conserved domains) <install_ncbi_tools>`
- :ref:`NCBI Conserved Domain Database (finding conserved domains) <install_ncbi_tools>`
- :ref:`Markov Clustering (creating gene phamilies) <install_mcl>`
- :ref:`tRNAScan-SE (importing tRNA and tmRNA data) <install_trnascan_se>`
- :ref:`ARAGORN (importing tRNA and tmRNA data) <install_aragorn>`
- :ref:`Infernal (importing tRNA and tmRNA data) <install_infernal>`
- :ref:`pdm_utils repository (source code, MySQL scripts, tests) <install_pdm_utils_repo>`
