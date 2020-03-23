.. _install_aragorn:


ARAGORN
=======

ARAGORN is used to evaluate tRNA features when importing data into the database.



MacOS installation
******************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>


#. ARAGORN can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd


#. Download the latest version of the tarball::

    (pdm_utils)> curl http://130.235.244.92/ARAGORN/Downloads/aragorn1.2.38.tgz > aragorn1.2.38.tgz






#. Uncompress the tarball, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xvzf aragorn1.2.38.tgz
        (pdm_utils)> rm aragorn1.2.38.tgz



#. Navigate to the new subdirectory and compile the program::

        (pdm_utils)> cd aragorn1.2.38
        (pdm_utils)> gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.38.c


#. Make ARAGORN globally executable. This can be accomplished one of two ways:

    - Get the full path to the installed software::


        (pdm_utils)> pwd

    - It should be something like::

        /Users/<user name>/aragorn1.2.38

    - Create a file that contains the full path to ARAGORN.

        A. Enter the following command, provide your login password when prompted::

            (pdm_utils)> sudo sh -c "echo '/Users/<user name>/aragorn1.2.38' > /etc/paths.d/aragorn"

        .. note::
             Be sure to use double and single quoting as indicated above.

        B. Exit Terminal.
        C. Open a new Terminal window so that the changes are implemented.


    - Alternatively, edit Terminal's '.bash_profile' configuration file using the Nano text editor.

        A. Navigate to the home directory::

            (pdm_utils)> cd

        B. Open the configuration file::

            (pdm_utils)> nano .bash_profile

        C. Navigate the cursor to the bottom of the file.
        D. Press <enter> to add a new line.
        E. Add the following text::

            export PATH=$PATH:/Users/<user name>/aragorn1.2.38

        F. Exit the file by typing Ctrl+x.
        G. When prompted to save the file type: 'y'. Then hit <enter>, which will return you to the command line prompt.
        H. Exit Terminal.
        I. Open a new Terminal window so that the changes are implemented.


#. Test whether ARAGORN has been successfully installed:

    - Open a new Terminal window.

    - Run the following command::

        (pdm_utils)> aragorn -h

    - If successful, a detailed description of the software's options should be displayed.

    - If unsuccessful, an error message should be displayed, such as::

        -bash: aragorn: command not found
