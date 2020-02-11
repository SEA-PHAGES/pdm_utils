.. _install_legacy_blast:


Legacy Blast
============



The ``phamerate`` tool can implement different clustering software to group genes into phamilies. Although MMseqs2 is the primary strategy currently utilized, blastclust is an alternative that ``phamerate`` is capable of implementing as well. Blastclust is part of the NCBI's legacy blast binaries, and they are no longer supported. General installation directions are provided below:


1. Navigate to :blastlegacy:`Legacy BLAST FTP <>`.

2. Download the appropriate compiled binaries based on your operating system:

    - For MacOS::

        blast-2.2.14-universal-macosx.tar.gz

    - For Ubuntu::

        blast-2.2.14-x64-linux.tar.gz


3. Move the downloaded, zipped file to any directory.
4. Unzip the binaries (example below is for MacOS, so modify the filename for Ubuntu)::

    (pdm_utils)> tar -xzvf blast-2.2.14-universal-macosx.tar.gz
    (pdm_utils)> rm blast-2.2.14-universal-macosx.tar.gz


.. warning::
    Since these binaries are no longer supported by NCBI, and since other ``pdm_utils`` tools require the currently-supported NCBI Blast+ software, it is recommended that the legacy software should not be added to the Terminal configuration file as for Blast+. When using legacy Blast with the ``phamerate`` tool, the path to the binaries can be directly indicated.
