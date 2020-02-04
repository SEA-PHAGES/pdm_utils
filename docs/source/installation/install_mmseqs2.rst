.. _install_mmseqs2:


MMseqs2
=======

MMseqs2 is a suite of tools for rapid searching and clustering of protein and nucleotide sequences (:ref:`Steinegger and Söding, 2017 <bibliography>`). This software is required only if gene phamilies need to be identified using MMseqs2 in the ``phamerate`` pipeline. Complete, detailed installation instructions are provided by the developers on the project's :mmseqs2:`GitHub page <>`. The following instructions provide an example of how it can be installed on your local MacOS or Ubuntu machine.

    #. Open a Terminal window and start the Conda environment::

        > conda activate pdm_utils
        (pdm_utils)>


    #. MMseqs2 can be installed anywhere on your computer. For example, navigate to your home directory::

        (pdm_utils)> cd

    #. Download the latest version of the compiled software:

        - For MacOS::

            (pdm_utils)> curl https://mmseqs.com/latest/mmseqs-osx-sse41.tar.gz > mmseqs-osx-sse41.tar.gz

        - For Ubuntu::

            (pdm_utils)> curl https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz > mmseqs-linux-sse41.tar.gz

    #. Unzip the archived package, then delete the original archive, by running the following commands:

        - For MacOS::

            (pdm_utils)> tar -xzvf mmseqs-osx-sse41.tar.gz
            (pdm_utils)> rm mmseqs-osx-sse41.tar.gz

        - For Ubuntu::

            (pdm_utils)> tar -xzvf mmseqs-linux-sse41.tar.gz
            (pdm_utils)> rm mmseqs-linux-sse41.tar.gz

        - An 'mmseqs' subdirectory should now exist within the directory.

    #. Get the complete path to the installed MMseqs2 directory:

        - Enter the following commands::

            (pdm_utils)> cd ./mmseqs/bin
            (pdm_utils)> pwd

        - For MacOS, the full path may be similar to::

                /Users/<user name>/mmseqs/bin

        - For Ubuntu, the full path may be similar to::

                /home/<user name>/mmseqs/bin

    #. Navigate to the home directory::

        (pdm_utils)> cd

    #. Use the nano text editor to edit Terminal's configuration file:

        * For MacOS:

            A. Open the file::

                (pdm_utils)> nano .bash_profile

            B. Navigate the cursor to the line above the 'export PATH' statement.
            C. Press <enter> to add a new line.
            D. Add the following text::

                PATH=$PATH:/Users/<user name>/mmseqs/bin

        * For Ubuntu:

            A. Open the file::

                (pdm_utils)> nano .bashrc

            B. Navigate the cursor to the end of the file.
            C. Add the following text::

                export PATH=$PATH:/home/<user name>/mmseqs/bin


    #. Exit the file by typing: <control> + 'X'.
    #. When prompted to save the file type: 'y'. Then hit <enter>, which will return you to the command line prompt.
    #. Exit Terminal.
    #. Open a new Terminal window so that the changes are implemented.
    #. Test whether MMseqs2 has been successfully installed:

        - Run the following command::

            (pdm_utils)> mmseqs cluster -h

        - If successful, a detailed description of the software's options should be displayed.

        - If unsuccessful, an error message should be displayed, such as::

            -bash: mmseqs: command not found
