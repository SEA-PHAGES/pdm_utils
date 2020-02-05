.. _install_mmseqs2:


MMseqs2
=======

MMseqs2 is a suite of tools for rapid searching and clustering of protein and nucleotide sequences (:ref:`Steinegger and Söding, 2017 <bibliography>`). This software is required only if gene phamilies need to be identified using MMseqs2 in the ``phamerate`` pipeline. Complete, detailed installation instructions are provided by the developers on the project's :mmseqs2:`GitHub page <>`. The following instructions provide an example of how it can be installed on your local MacOS or Ubuntu machine.


MacOS installation
******************

#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>


#. MMseqs2 can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd


#. Download the latest version of the compiled software::

    (pdm_utils)> curl https://mmseqs.com/latest/mmseqs-osx-sse41.tar.gz > mmseqs-osx-sse41.tar.gz


#. Unzip the archived package, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xzvf mmseqs-osx-sse41.tar.gz
        (pdm_utils)> rm mmseqs-osx-sse41.tar.gz


    - An 'mmseqs' subdirectory should now exist within the directory.


#. Get the complete path to the installed MMseqs2 directory:

    - Enter the following commands::

        (pdm_utils)> cd ./mmseqs/bin
        (pdm_utils)> pwd

    - The full path may be similar to::

        /Users/<user name>/mmseqs/bin


#. Add the installed MMseqs2 software to your Terminal environment. This can be accomplished one of two ways:

    - Create a file that contains the full path to MMseqs2.

        A. Enter the following command, provide your login password when prompted::

            (pdm_utils)> sudo sh -c "echo '/Users/<user name>/mmseqs/bin' > /etc/paths.d/mmseqs"

        .. note::
             Be sure to use double and single quoting as indicated above.

        B. Exit Terminal.
        C. Open a new Terminal window so that the changes are implemented.


    - Alternatively, edit Terminal's '.bash_profile' configuration file using the Nano text editor.

        A. Navigate to the home directory::

            (pdm_utils)> cd

        B. Open the configuration file::

            (pdm_utils)> nano .bash_profile

        C. Navigate the cursor to the line above the 'export PATH' statement.
        D. Press <enter> to add a new line.
        E. Add the following text::

            PATH=$PATH:/Users/<user name>/mmseqs/bin

        F. Exit the file by typing Ctrl+x.
        G. When prompted to save the file type: 'y'. Then hit <enter>, which will return you to the command line prompt.
        H. Exit Terminal.
        I. Open a new Terminal window so that the changes are implemented.


#. Test whether MMseqs2 has been successfully installed:

    - Run the following command::

        (pdm_utils)> mmseqs cluster -h

    - If successful, a detailed description of the software's options should be displayed.

    - If unsuccessful, an error message should be displayed, such as::

        -bash: mmseqs: command not found




Ubuntu installation
*******************


#. Open a Terminal window and start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>


#. MMseqs2 can be installed anywhere on your computer. For example, navigate to your home directory::

    (pdm_utils)> cd


#. Download the latest version of the compiled software::

    (pdm_utils)> curl https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz > mmseqs-linux-sse41.tar.gz

#. Unzip the archived package, then delete the original archive:

    - Enter the following commands::

        (pdm_utils)> tar -xzvf mmseqs-linux-sse41.tar.gz
        (pdm_utils)> rm mmseqs-linux-sse41.tar.gz

    - An 'mmseqs' subdirectory should now exist within the directory.


#. Get the complete path to the installed MMseqs2 directory:

    - Enter the following commands::

        (pdm_utils)> cd ./mmseqs/bin
        (pdm_utils)> pwd

    - The full path may be similar to::

        /home/<user name>/mmseqs/bin





#. Add the installed MMseqs2 software to your Terminal environment. This can be accomplished one of two ways:


    - Edit Terminal's '.bashrc' configuration file indirectly:

        A. Navigate to the home directory::

            (pdm_utils)> cd

        B. Enter the following command::

            (pdm_utils)> echo "export PATH=$PATH:/home/<user name>/mmseqs/bin" >> .bash_profile

        C. Exit Terminal.
        D. Open a new Terminal window so that the changes are implemented.


    - Alternatively, edit Terminal's '.bashrc' configuration file directly using the Nano text editor:

        A. Navigate to the home directory::

            (pdm_utils)> cd

        B. Open the configuration file::

            (pdm_utils)> nano .bashrc

        C. Navigate the cursor to the end of the file.
        D. Add the following text::

            export PATH=$PATH:/home/<user name>/mmseqs/bin

        E. Exit the file by typing Ctrl+x.
        F. When prompted to save the file type: 'y'. Then hit <enter>, which will return you to the command line prompt.
        G. Exit Terminal.
        H. Open a new Terminal window so that the changes are implemented.



#. Test whether MMseqs2 has been successfully installed:

    - Run the following command::

        (pdm_utils)> mmseqs cluster -h

    - If successful, a detailed description of the software's options should be displayed.

    - If unsuccessful, an error message should be displayed, such as::

        -bash: mmseqs: command not found
