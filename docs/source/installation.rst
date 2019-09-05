Installation
============


The ``pdm_utils`` package is written in Python3 and has several dependencies:

Built-in python packages:

    - os
    - shutil
    - sys
    - getpass
    - subprocess
    - re
    - datetime
    - operator
    - random
    - colorsys



Third-party python packages:

    - `Biopython <https://biopython.org/>`_
    - `pymysql <https://pymysql.readthedocs.io/en/latest/>`_
    - `paramiko <http://www.paramiko.org/>`_


Third-party binaries installed locally:

    - `MySQL 5.7 <https://www.mysql.com/>`_
    - `MMSeqs <https://www.ncbi.nlm.nih.gov/pubmed/26743509>`_
    - `NCBI blast+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_
    - `NCBI Conserved Domain Database <https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml>`_




In order to install and manage all python dependencies, the Conda environment manage can be used, which is available through the `Anaconda <https://www.anaconda.com/>`_.

After installing Conda, create an environment for phamerator:

conda create --name phamerator mysql python pip biopython pymysql paramiko sphinx sphinx_rtd_theme
