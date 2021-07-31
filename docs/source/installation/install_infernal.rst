.. _install_infernal:


Infernal
========

:infernal:`Infernal <>` (:ref:`Nawrocki & Eddy, 2013 <bibliography>`) is required to properly run tRNAscan-SE.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install Infernal::

    (pdm_utils)> conda install -c bioconda infernal -y

#. Test whether Infernal has been successfully installed::

    (pdm_utils)> cmsearch -h

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

        -bash: cmsearch: command not found

