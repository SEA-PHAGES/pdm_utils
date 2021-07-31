.. _install_ncbi_tools:


NCBI Tools
==========

NCBI BLAST+
___________


The BLAST+ software package is a stand-alone version of BLAST applications and is maintained by the NCBI (:ref:`Camacho et al., 2009 <bibliography>`). This software is required in the ``find_domains`` pipeline to identify conserved domains within genes. Follow the installation instructions specific to your operating system provided in the :blastplus:`BLAST+ user guide <>`.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. The most straightforward option is to use Conda to install blast::

    (pdm_utils)> conda install -c bioconda blast -y

#. Test whether blast has been successfully installed::

    (pdm_utils)> blastp -help

If successful, a detailed description of the software's options should be displayed.

If unsuccessful, an error message should be displayed, such as::

        -bash: blastp: command not found

NCBI Conserved Domain Database
______________________________


The :cdd:`NCBI Conserved Domain Database <>` is a curated set of protein domain families (:ref:`Lu et al., 2020 <bibliography>`). This database is required in the ``find_domains`` pipeline to identify conserved domains within genes.

1. Download the compressed :cdd_le:`NCBI CDD <>`.

2. Expand the archived file into a directory of CDD files. The entire directory represents the database and no files should be added or removed.
