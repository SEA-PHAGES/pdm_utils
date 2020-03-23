.. _install_trnascan_se:


tRNAScan-SE v2.0
================

:trna_scan_se:`tRNAScan-SE v2.0 <>` (:ref:`Chan & Lowe, 2019 <bibliography>`) is used to evaluate tRNA features when importing data into the database.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>


#. tRNAScan-SE can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd

#. Download the latest version of the tarball::

    (pdm_utils)> curl http://trna.ucsc.edu/software/trnascan-se-2.0.5.tar.gz > trnascan-se-2.0.5.tar.gz

#. Uncompress the tarball, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xvzf trnascan-se-2.0.5.tar.gz
        (pdm_utils)> rm trnascan-se-2.0.5.tar.gz

#. Navigate to the new subdirectory and compile the program::

    (pdm_utils)> cd tRNAscan-SE-2.0/
    (pdm_utils)> ./configure
    (pdm_utils)> make

#. Make the program globally executable::

    (pdm_utils)> sudo make install

#. Test whether tRNAScan-SE has been successfully installed:

    - Open a new Terminal window.

    - Run the following command::

        (pdm_utils)> tRNAscan-SE -h

    - If successful, a detailed description of the software's options should be displayed.

    - If unsuccessful, an error message should be displayed, such as::

        -bash: tRNAscan-SE: command not found
