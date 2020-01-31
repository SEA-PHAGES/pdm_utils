.. _installation:

Installation
============


The ``pdm_utils`` package is written in Python 3 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies that are only needed for a subset of ``pdm_utils`` tools. Below is a general step-by-step guide to installing all dependencies.


1. MySQL
________

:mysql:`MySQL<>` is a database management system, and it is required for practically all ``pdm_utils`` tools. ``pdm_utils`` has been developed and tested using MySQL Community Server 5.7. The package may work with newer versions of MySQL, but it has not been tested. The :mysql_docs:`MySQL user guide <>` provides very thorough, detailed installation instructions, and there are many online tutorials to help overcome common installation challenges. Below is a brief summary of installation for MacOS and Ubuntu.


MacOS installation
******************

Installation
^^^^^^^^^^^^

#. Navigate to the :mysql:`MySQL Community Server download site <downloads/mysql/>`.

#. Verify the 'General Availability (GA) Releases' tab is selected.

#. Verify 'macOS' is selected from the 'Select Operating System' menu.

#. By default, the most recent version of MySQL Community Server (version 8.0) is displayed for download. Select the 'Looking for previous GA versions?' option. Verify that version 5.7 is now displayed.

#. Download the 'macOS DMG Archive'.

#. Follow the graphical installation instructions:

    .. figure:: /images/mac_mysql_install/step1.png


#. Agree with the license agreement:

    .. figure:: /images/mac_mysql_install/step2.png

    .. figure:: /images/mac_mysql_install/step3.png

#. Install at the default location:

    .. figure:: /images/mac_mysql_install/step4.png

#. Allow access to 'System Events.app':

    .. figure:: /images/mac_mysql_install/step5.png

#. Make note of the temporary password provided during the installation process (highlighted in red):

    .. figure:: /images/mac_mysql_install/step6.png

    .. warning::
         Be sure to record the temporary password that is generated!!! Once MySQL is installed, this password can only be used ONCE to login.

#. Installation should now be complete:

    .. figure:: /images/mac_mysql_install/step7.png

#. Open Terminal.
#. Log in to MySQL as the root user::

    > mysql -u root -p

#. Enter the temporary password when prompted::

    Enter password: <temporary password>

#. At the mysql prompt, change the password for the root user::

    mysql> ALTER USER 'root'@'localhost' IDENTIFIED BY '<new password>';
    mysql> exit



Password reset
^^^^^^^^^^^^^^

If the MySQL password is lost, it can be reset.

#. In Mac Preferences Pane, turn off MySQL server.
#. Open a Terminal window.
#. Enter the following command::

    > sudo /usr/local/mysql/bin/mysqld_safe --skip-grant-tables

#. Enter the password for the computer user (not MySQL), at the prompt::

    Enter password: <user password>

#. Open a SECOND Terminal window.
#. Enter the following command::

    > sudo /usr/local/mysql/bin/mysql -u root

#. Enter the password for the computer user (not MySQL), at the prompt::

    Enter password: <user password>

#. You should now be logged in to mysql. Execute the following commands::

        mysql> UPDATE mysql.user SET authentication_string=PASSWORD('<new password>') WHERE User='root';
        mysql> FLUSH PRIVILEGES;
        mysql> exit

#. You should now be returned to the bash command line. Enter the following command::

    > sudo /usr/local/mysql/support-files/mysql.server restart

#. Close the second Terminal window.
#. Close the first Terminal window.



Server control
^^^^^^^^^^^^^^

Access to MySQL, even on your local computer, is controlled through a server-client model. The server needs to be turned ON in order to use MySQL. This can be accomplished manually or it can be set to start automatically every time your Mac is restarted.

#. Click on the Apple icon in the top left corner of your desktop.
#. Select 'System Preferences'.
#. Click on the MySQL icon.
#. If 'MySQL Server Instance is stopped' is displayed, then click on 'Start MySQL Server'.
#. To perform this step every time automatically, select 'Automatically Start MySQL Server on Startup'.

If the automatic option is not selected, anytime your Mac is restarted the server is turned OFF, and you will be unable to use any ``pdm_utils`` tools that require access to MySQL until you manually turn the server ON.


Ubuntu installation
*******************

Installing MySQL on Ubuntu is more straightforward. MySQL 5.7 can be downloaded through either the Ubuntu repositories or the official MySQL repositories. Installing MySQL using the Ubuntu repositories is outlined below:

#. Open a Terminal window.
#. Update all available repositories (provide the computer login password when prompted)::

    > sudo apt update

#. Enter the following command to install the MySQL version 5.7 (answer 'yes' to proceed with installing the new packages, when prompted)::

    > sudo apt install mysql-server=5.7.*

#. MySQL Community Server should now be installed, but the server may not be running.

    - Check the server status:

        A. Enter the following command::

            > systemctl status mysql.service

        B. If the server is running, it should display::

            Active: active (running))

        C. If the server is not running, it should display::

            Active: inactive (dead)

    - If the server is not running, it needs to be started::

        > sudo systemctl start mysql


    - Check status again to confirm it is running::

        > systemctl status mysql.service


#. Although MySQL is installed, no password has yet been set for the 'root' user. Login to MySQL without a username (provide the computer login password if prompted)::

    > sudo mysql
    mysql>

#. Now set a password for the 'root' user::

    mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY '<new password>';
    mysql> FLUSH PRIVILEGES;
    mysql> exit;




Create additional users (optional)
**********************************

After MySQL is installed (on MacOS or Ubuntu), additional user accounts with different types of access privileges can be created, if needed.

#. Login to mysql as 'root' (provide the password when prompted)::

    > mysql -u root -p
    mysql>

#. Create a new user 'new_user', and specify the password::

    mysql> CREATE USER 'new_user'@'localhost' IDENTIFIED BY '<new_password>';

#. Grant different levels of access using one of the following commands:

    - Grant unrestricted access to all databases::

        mysql> GRANT ALL ON *.* TO 'new_user'@'localhost' WITH GRANT OPTION;

    - Grant access with all privileges to a specific database (such as Actinobacteriophage)::

        mysql> GRANT ALL ON Actinobacteriophage.* TO 'new_user'@'localhost';

    - Grant access to all databases, but only with the privilege to retrieve data::

        mysql> GRANT SELECT ON *.* TO 'new_user'@'localhost';

#. Implement the changes::

    mysql> FLUSH PRIVILEGES;
    mysql> exit;



2. Anaconda
___________

There are several third-party Python packages required:

    - :biopython:`Biopython <>`
    - :paramiko:`Paramiko <>`
    - :pymysql:`PyMySQL <>`
    - :sqlalchemy:`SQLAlchemy <>`
    - :tabulate:`Tabulate <>`

Some of them also have Python or binary dependencies. Manual installation of these dependencies can be tricky, but the Conda environment manager is a simple, automated alternative. First install Conda, then use Conda to install Python and other dependencies.

Conda is available as part of Anaconda or Miniconda, and complete installation instructions are available in the Conda :conda_docs:`user guide <>`. The directions below highlight installation of Anaconda, but either of these tools is fine since they both install Conda:

#. Install Conda locally through the :anaconda:`Anaconda <>` package.

#. Navigate to the 'Anaconda Distribution' option.

#. Begin the installation:

    - For MacOS: download the Python 3.7 graphical installer and follow the graphical prompts.

    - For Linux:

        A. Download the Python 3.7 x86 Linux installer (e.g. Anaconda3-2019.10-Linux-x86_64.sh) to the Downloads folder.
        B. Open a Terminal window.
        C. Execute the following command::

            > bash ~/Downloads/Anaconda3-2019.10-Linux-x86_64.sh


#. Follow the manufacturer's installation instructions.

    - Accept the license agreement.
    - Install at the default directory.
    - Enter 'yes' when prompted for the installer to run conda init.

#. Optional: execute the following command to prevent Conda from starting up automatically every time a Terminal window is opened::

    > conda config --set auto_activate_base false

#. Close the Terminal window and open a new window.

#. After installing Conda, create an environment to be able to install and use ``pdm_utils`` (the example below creates a Conda environment named 'pdm_utils', but it can be named anything). Enter 'y' when prompted to install all dependencies::

    > conda create --name pdm_utils python pip biopython pymysql paramiko tabulate curl sqlalchemy

#. After the Conda environment is created, it needs to be activated using the following command. The command line prompt will now include '(pdm_utils)', indicating it is operating within this environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. Optional: enter the following command to exit the Conda environment. The default command line prompt will be displayed, and the name of the Conda environment will no longer be displayed::

    (pdm_utils)> conda deactivate
    >


.. note::

    If Conda is used to manage dependencies, the Conda environment must be activated every time you want to use ``pdm_utils``. Otherwise, an error will be encountered.


The 'pdm_utils' Conda environment now contains the necessary dependencies, and the actual ``pdm_utils`` Python package can be installed (see below).


3. The ``pdm_utils`` package
____________________________

Once MySQL and the Conda environment are installed, ``pdm_utils`` can be easily installed:

#. Open a Terminal window.

#. Activate the Conda environment (see above).

#. Install the ``pdm_utils`` package using pip::

    (pdm_utils)> pip install pdm_utils

#. The package is routinely updated, and the most recent version can be retrieved::

    (pdm_utils)> pip install --upgrade pdm_utils


4. MySQL database instance
_______________________________

Many ``pdm_utils`` modules and pipelines require access to a specifically structured MySQL database.

The primary database instance that reflects the most up-to-date actinobacteriophage genomics data in the SEA-PHAGES program is the 'Actinobacteriophage' database. Typically, different versions, or instances, of the database are created ('frozen') for specific studies/publications. The unique name of the database is normally published in the Materials and Methods.

The ``pdm_utils get_db`` installation management tool can be used to retrieve, install, and update these databases, or any custom MySQL database that is compliant with the database schema, from a local file or from the Hatfull lab server (:ref:`get_db <getdb>`).

Alternatively, databases can be manually downloaded and installed, as described below (using Actinobacteriophage as an example):

Manual installation
*******************

#. Open a Terminal window.
#. Create an empty database (enter your password when prompted)::

    > mysql -u root -p --execute "CREATE DATABASE Actinobacteriophage"

#. Retrieve the current version of the database::

    > curl http://phamerator.webfactional.com/databases_Hatfull/Actinobacteriophage.sql > ./Actinobacteriophage.sql

#. Import the database into MySQL (enter your password when prompted)::

    > mysql -u root -p Actinobacteriophage < Actinobacteriophage.sql


Manual update
*************

#. Log in to MySQL (enter your password when prompted)::

    > mysql -u root -p

#. Execute the following query to get the current version::

    mysql> SELECT Version FROM version;
    mysql> exit

#. Download the current version file from the Hatfull lab server::

    > curl http://phamerator.webfactional.com/databases_Hatfull/Actinobacteriophage.version > ./Actinobacteriophage.version

#. If the current version on the server is different from the version in the local MySQL database, there is a new database available on the server. Repeat steps 3-4 listed above in the 'Manual Installation' section.




5. Tool-specific dependencies
_____________________________

Several ``pdm_utils`` tools have specific dependencies. Install the following tools/files as needed.


MMseqs2
*******

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




NCBI BLAST+
***********

The BLAST+ software package is a stand-alone version of BLAST applications and is maintained by the NCBI (:ref:`Camacho et al., 2009 <bibliography>`). This software is required in the ``cdd`` pipeline to identify conserved domains within genes. Follow the installation instructions specific to your operating system provided in the :blastplus:`BLAST+ user guide <>`.


NCBI Conserved Domain Database
******************************


The :cdd:`NCBI Conserved Domain Database <>` is a curated set of protein domain families (:ref:`Lu et al., 2020 <bibliography>`). This database is required in the ``cdd`` pipeline to identify conserved domains within genes.

1. Download the compressed :cdd_le:`NCBI CDD <>`.

2. Expand the archived file into a directory of CDD files. The entire directory represents the database and no files should be added or removed.



pdm_utils source code repository
********************************

Some ``pdm_utils`` tools may require non-Python data files that are not directly installed with the Python package. Instead, these files are available on the ``pdm_utils`` git repository, which can be accessed through :pdmutils:`GitHub <>` and downloaded as follows:

    #. Open a Terminal window.
    #. Navigate to a directory where the source code subdirectory should be installed (for example, your home directory)::

        > cd

    #. Enter the following command::

        > git clone https://github.com/SEA-PHAGES/pdm_utils.git

There should now be a ``pdm_utils`` directory installed locally.
