.. _install_mmseqs2:


MMseqs2
=======

MMseqs2 is a suite of tools for rapid searching and clustering of protein and nucleotide sequences (:ref:`Steinegger and Söding, 2017 <bibliography>`). This software is required only if gene phamilies need to be identified using MMseqs2 in the ``phamerate`` pipeline. Complete, detailed installation instructions are provided by the developers on the project's :mmseqs2:`GitHub page <>`. The following instructions provide an example of how it can be installed on your local MacOS or Ubuntu machine.


MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install MMseqs2::

    (pdm_utils)> conda install -c bioconda -c conda-forge mmseqs2=13.45111 -y

#. Test whether MMseqs2 has been successfully installed::

    (pdm_utils)> mmseqs cluster -h

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

        -bash: mmseqs: command not found

