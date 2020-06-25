.. _install_mcl:

Markov Clustering
=================

Markov Clustering (MCL) is a graph clustering program (:ref:`van Dongen & Abreu-Goodger, 2012 <bibliography>`) that
can be used in combination a pairwise blastp output to infer phamilies of homologous proteins.  Installation is the
same for both MacOS and Linux operating systems.

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils

#. Markov Clustering can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd

#. Download the latest version of the compiled software::

    (pdm_utils)> curl https://micans.org/mcl/src/mcl-14-137.tar.gz > mcl-14-137.tar.gz

#. Unzip the archived package, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xzvf mcl-14-137.tar.gz
        (pdm_utils)> rm mcl-14-137.tar.gz

#. Enter the mcl-14-137 directory, run the configuration script, and install::

    (pdm_utils)> cd mcl-14-137
    (pdm_utils)> ./configure
    (pdm_utils)> sudo make install

You may be prompted for your computer password after this last command. Enter it, and the Markov Clustering programs
should be compiled and installed.

Test to make sure that mcl was properly installed. Open a new Terminal window and enter the following command::

    > mcl --help


If installation was successful, a detailed description of the software's options should be displayed.

Otherwise, an error message should be displayed, such as::

        -bash: mcl: command not found

