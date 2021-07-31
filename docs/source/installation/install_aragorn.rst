.. _install_aragorn:


ARAGORN
=======

:aragorn:`ARAGORN <>` (:ref:`Laslett & Canback, 2004 <bibliography>`) is used to evaluate tRNA features when importing data into the database.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install Aragorn::

    (pdm_utils)> conda install -c bioconda aragorn -y

#. Test whether ARAGORN has been successfully installed::

    (pdm_utils)> aragorn -h

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

    -bash: aragorn: command not found

