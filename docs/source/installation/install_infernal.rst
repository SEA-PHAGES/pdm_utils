.. _install_infernal:


Infernal
========

:infernal:`Infernal <>` (:ref:`Nawrocki & Eddy, 2013 <bibliography>`) is used to evaluate tRNA features when importing data into the database.

MacOS and Ubuntu installation
*****************************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. Infernal can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd

#. Download the latest version of the tarball::

    (pdm_utils)> curl http://eddylab.org/infernal/infernal-1.1.3.tar.gz > infernal-1.1.3.tar.gz

#. Uncompress the tarball, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xvzf infernal-1.1.3.tar.gz
        (pdm_utils)> rm infernal-1.1.3.tar.gz

#. Navigate to the new subdirectory and compile the program::

    (pdm_utils)> cd infernal-1.1.3
    (pdm_utils)> ./configure
    (pdm_utils)> make

#. Make the program globally executable::

    (pdm_utils)> sudo make install

#. Test whether Infernal has been successfully installed:

    - Open a new Terminal window.

    - Run the following command::

        (pdm_utils)> cmsearch -h

    - If successful, a detailed description of the software's options should be displayed.

    - If unsuccessful, an error message should be displayed, such as::

        -bash: cmsearch: command not found
