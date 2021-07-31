.. _install_mcl:

Markov Clustering
=================

Markov Clustering (MCL) is a graph clustering program (:ref:`van Dongen & Abreu-Goodger, 2012 <bibliography>`) that
can be used in combination a pairwise blastp output to infer phamilies of homologous proteins.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install MCL::

    (pdm_utils)> conda install -c bioconda mcl -y

#. Test whether MCL has been successfully installed::

    (pdm_utils)> mcl --help

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

        -bash: mcl: command not found

