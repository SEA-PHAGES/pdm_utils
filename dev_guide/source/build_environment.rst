Setting up dependencies
=======================

The Sphinx documentation provides a list of python and binary dependencies to use the package. Additionally, there are several third-party python packages required for development:

    - sphinx
    - sphinx_rtd_theme
    - twine

All python code in the repo can be executed using the following Conda environment 'pdm_utils-dev'. First install Anaconda locally. Then execute the following commands::

    > conda create --name pdm_utils-dev python pip biopython pymysql paramiko sphinx sphinx_rtd_theme twine tabulate curl sqlalchemy networkx
    > conda activate pdm_utils-dev


To create a list of pinned Python and non-Python packages and binaries installed in the conda environment::

    (pdm_utils-dev)> conda env export | grep -v "^prefix: " > environment.yml


To create a list of just the pinned Python packages installed in the conda environment::

    (pdm_utils-dev)> pip freeze > requirements.txt


To add ``pdm_utils`` module to python's syspath in a shell environment on-the-fly, use the following command, edited so that it points to your local repo::

    > export PYTHONPATH="${PYTHONPATH}:/path/to/pdm_utils/src/"
