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

    - :biopython:`Biopython <>`
    - :pymysql:`pymysql <>`
    - :paramiko:`paramiko <>`


Third-party binaries installed locally:

    - :mysql:`MySQL 5.7 <>`
    - :mmseqs:`MMSeqs <>`
    - :blastplus:`NCBI blast+ <>`
    - :cdd:`NCBI Conserved Domain Database <>`




In order to install and manage all python dependencies, the Conda environment manage can be used, which is available through the :anaconda:`Anaconda <>` package.

After installing Conda, create an environment for phamerator::

    conda create --name phamerator mysql python pip biopython pymysql paramiko sphinx sphinx_rtd_theme
