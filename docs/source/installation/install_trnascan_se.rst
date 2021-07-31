.. _install_trnascan_se:


tRNAScan-SE v2.0
================

:trna_scan_se:`tRNAScan-SE v2.0 <>` (:ref:`Chan & Lowe, 2019 <bibliography>`) is used to evaluate tRNA features when importing data into the database.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install tRNAscan-SE::

    (pdm_utils)> conda install -c bioconda trnascan-se -y

#. Test whether tRNAscan-SE has been successfully installed::

    (pdm_utils)> tRNAscan-SE -h

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

        -bash: tRNAscan-SE: command not found

