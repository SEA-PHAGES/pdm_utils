.. _installation:

Installation
============


The ``pdm_utils`` package is written in Python 3 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies specific to a subset of ``pdm_utils`` tools or pipelines that are only needed for specific goals. Below is a general step-by-step guide to installing all dependencies and the ``pdm_utils`` package.


1. MySQL
________

:mysql:`MySQL Community Server 5.7 <>`

    This is required for practically all ``pdm_utils`` tools. ``pdm_utils`` has been developed and tested using MySQL 5.7. The package may work with newer versions of MySQL, but it has not been tested. Below is a brief summary of the installation steps. Refer to the official MySQL documentation for more details.



MacOS installation
******************

Installing MySQL on MacOS can be tricky.

    1. Download the MySQL Community Server 5.7 dmg file from the MySQL website.

        1. Follow the graphical installation instructions:

        .. figure:: /images/mac_mysql_install/step1.png


        2. Agree with the license agreement:

        .. figure:: /images/mac_mysql_install/step2.png

        .. figure:: /images/mac_mysql_install/step3.png

        3. Install at the default location:

        .. figure:: /images/mac_mysql_install/step4.png

        4. Allow access to 'System Events.app':

        .. figure:: /images/mac_mysql_install/step5.png

        5. Make note of the temporary password provided during the installation process (highlighted in red):

        .. figure:: /images/mac_mysql_install/step6.png

        .. warning::
             Be sure to record the temporary password that is generated!!! Once MySQL is installed, this password can only be used ONCE to login.

        6. Installation is complete:

        .. figure:: /images/mac_mysql_install/step7.png


    2. Open Terminal.
    3. Log in to MySQL as the root user::

        > mysql -u root -p

    4. Enter the temporary password when prompted::

        Enter password: <temporary password>

    5. At the mysql prompt, change the password for the root user::

        mysql> ALTER USER 'root'@'localhost' IDENTIFIED BY '<new password>';
        mysql> exit


If the MySQL password is lost, it can be reset.

    1. In Mac Preferences Pane, turn off MySQL server.
    2. Open a Terminal window.
    3. Enter the following command::

        > sudo /usr/local/mysql/bin/mysqld_safe --skip-grant-tables

    4. Enter the password for the computer user (not MySQL), at the prompt::

        Enter password: <user password>

    5. Open a SECOND Terminal window.
    6. Enter the following command::

        > sudo /usr/local/mysql/bin/mysql -u root

    7. Enter the password for the computer user (not MySQL), at the prompt::

        Enter password: <user password>

    8. You should now be logged in to mysql. Execute the following commands::

            mysql> UPDATE mysql.user SET authentication_string=PASSWORD('<new password>') WHERE User='root';
            mysql> FLUSH PRIVILEGES;
            mysql> exit

    9. You should now be returned to the bash command line. Enter the following command::

        > sudo /usr/local/mysql/support-files/mysql.server restart

    10. Close the second Terminal window.
    11. Close the first Terminal window.


Access to MySQL, even on your local computer, is controlled through a server-client model. The server needs to be turned ON in order to use MySQL. This can be accomplished manually or it can be set to start automatically every time your Mac is restarted.

    1. Click on the Apple icon in the top left corner of your desktop.
    2. Select 'System Preferences'.
    3. Click on the MySQL icon.
    4. If 'MySQL Server Instance is stopped' is displayed, then click on 'Start MySQL Server'.
    5. To perform this step every time automatically, select 'Automatically Start MySQL Server on Startup'.

If the automatic option is not selected, anytime your Mac is restarted the server is turned OFF, and you will be unable to use any ``pdm_utils`` tools that require access to MySQL until you manually turn the server ON.


Ubuntu installation
*******************

Installing MySQL on Ubuntu is more straightforward. MySQL 5.7 can be downloaded through either the Ubuntu repositories or the official MySQL repositories. Installing MySQL using the Ubuntu repositories is outlined below:

    1. Open a Terminal window.
    2. Update all available repositories (provide the computer login password when prompted)::

        > sudo apt update

    3. Enter the following command to install the MySQL version 5.7 (answer 'yes' to proceed with installing the new packages, when prompted)::

        > sudo apt install mysql-server=5.7.*

    4. MySQL Community Server should now be installed, but the server may not be running.

        A. Check the server status::

            > systemctl status mysql.service

        If the server is running, it should display::

            Active: active (running))

        If the server is not running, it should display::

            Active: inactive (dead)


        B. If the server is not running, it needs to be started::

            > sudo systemctl start mysql


        C. Check status again to confirm it is running::

            > systemctl status mysql.service


    7. Although MySQL is installed, no password has yet been set for the 'root' user. Login to MySQL without a username (provide the computer login password if prompted)::

        > sudo mysql
        mysql>

    8. Now set a password for the 'root' user::

        mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY '<new password>';
        mysql> FLUSH PRIVILEGES;
        mysql> exit;




Create additional users (optional)
**********************************

After MySQL is installed (on MacOS or Ubuntu), additional user accounts with different types of access privileges can be created, if needed.

    1. Login to mysql as 'root' (provide the password when prompted)::

        > mysql -u root -p
        mysql>

    2. Create a new user 'new_user', and specify the password::

        mysql> CREATE USER 'new_user'@'localhost' IDENTIFIED BY '<new_password>';

    3. Grant different levels of access using one of the following commands:

        - Grant unrestricted access to all databases::

            mysql> GRANT ALL ON *.* TO 'new_user'@'localhost' WITH GRANT OPTION;

        - Grant access with all privileges to a specific database (such as Actinobacteriophage)::

            mysql> GRANT ALL ON Actinobacteriophage.* TO 'new_user'@'localhost';

        - Grant access to all databases, but only with the privilege to retrieve data::

            mysql> GRANT SELECT ON *.* TO 'new_user'@'localhost';

    4. Implement the changes::

        mysql> FLUSH PRIVILEGES;
        mysql> exit;



2. Anaconda
___________

There are several third-party Python packages required:

    - :biopython:`Biopython <>`
    - :pymysql:`pymysql <>`
    - :paramiko:`paramiko <>`
    - tabulate

Some of them also have Python or binary dependencies. Manual installation of these dependencies can be tricky, but the Conda environment manager is a simple, automated alternative. First install Conda, then use Conda to install all Python dependencies. (Conda is available as part of Anaconda or Miniconda. The directions below highlight installation of Anaconda, but either of these tools is fine since they both install Conda):

    1. Install Conda locally through the :anaconda:`Anaconda <>` package.

    2. Navigate to the 'Anaconda Distribution' option.

    3. Begin the installation:

        - For MacOS: download the Python 3.7 graphical installer and follow the graphical prompts.

        - For Linux:

            1. Download the Python 3.7 x86 Linux installer (e.g. Anaconda3-2019.10-Linux-x86_64.sh) to the Downloads folder.
            2. Open a Terminal window.
            3. Execute the following command::

                > bash ~/Downloads/Anaconda3-2019.10-Linux-x86_64.sh


    4. Follow the manufacturer's installation instructions.

        - Accept the license agreement.
        - Install at the default directory.
        - Enter 'yes' when prompted for the installer to run conda init.

    5. Optional: execute the following command to prevent Conda from starting up automatically every time a Terminal window is opened::

        > conda config --set auto_activate_base false

    6. Close the Terminal window and open a new window.

    7. After installing Conda, create an environment to be able to install and use ``pdm_utils`` (the example below creates a Conda environment named 'pdm_utils', but it can be named anything). Enter 'y' when prompted to install all dependencies::

        > conda create --name pdm_utils python pip biopython pymysql paramiko tabulate curl sqlalchemy

    8. After the Conda environment is created, it needs to be activated using the following command. The command line prompt will now include '(pdm_utils)', indicating it is operating within this environment::

        > conda activate pdm_utils
        (pdm_utils)>

    9. Optional: enter the following command to exit the Conda environment. The default command line prompt will be displayed, and the name of the Conda environment will no longer be displayed::

        (pdm_utils)> conda deactivate
        >


.. note::

    If Conda is used to manage dependencies, the Conda environment must be activated every time you want to use ``pdm_utils``. Otherwise, an error will be encountered.


The 'pdm_utils' Conda environment now contains the necessary dependencies, and the actual ``pdm_utils`` Python package can be installed (see below).


3. The ``pdm_utils`` package
____________________________

Once MySQL and the Conda environment are installed, ``pdm_utils`` can be easily installed:

    1. Open a Terminal window.

    2. Activate the Conda environment (see above).

    3. Install the ``pdm_utils`` package using pip::

        (pdm_utils)> pip install pdm_utils

    4. The package is routinely updated, and the most recent version can be retrieved::

        (pdm_utils)> pip install --upgrade pdm_utils


4. MySQL database instance
_______________________________

Many ``pdm_utils`` modules and pipelines require access to a specifically structured MySQL database.

The primary database instance that reflects the most up-to-date actinobacteriophage genomics data in the SEA-PHAGES program is the 'Actinobacteriophage' database. Typically, different versions, or instances, of the database are created ('frozen') for specific studies/publications. The unique name of the database is normally published in the Materials and Methods.

The ``pdm_utils get_db`` installation management tool can be used to retrieve, install, and update these databases, or any custom MySQL database that is compliant with the database schema, from a local file or from the Hatfull lab server (:ref:`getdb <getdb>`).

Alternatively, databases can be manually downloaded and installed, as described below (using Actinobacteriophage as an example):

Manual installation
*******************

    1. Open a Terminal window.
    2. Create an empty database (enter your password when prompted)::

        > mysql -u root -p --execute "CREATE DATABASE Actinobacteriophage"

    3. Retrieve the current version of the database::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actinobacteriophage.sql > ./Actinobacteriophage.sql

    4. Import the database into MySQL (enter your password when prompted)::

        > mysql -u root -p Actinobacteriophage < Actinobacteriophage.sql


Manual update
*************

    1. Log in to MySQL (enter your password when prompted)::

        > mysql -u root -p

    2. Execute the following query to get the current version::

        mysql> SELECT Version FROM version;
        mysql> exit

    3. Download the current version file from the Hatfull lab server::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actinobacteriophage.version > ./Actinobacteriophage.version

    4. If the current version on the server is different from the version in the local MySQL database, there is a new database available on the server. Repeat steps 3-4 listed above in the 'Manual Installation' section.









5. Tool-specific dependencies
_____________________________

Several ``pdm_utils`` tools have specific dependencies. Install the following tools/files as needed.


MMSeqs
******


:mmseqs:`MMSeqs <>`

    Required only if gene phamilies need to be identified using MMSeqs in the 'phamerate' pipeline.

# TODO add installation instructions.


NCBI Blast+ toolkit
*******************

:blastplus:`NCBI blast+ <>`

    Required only if conserved domains within genes need to be identified from the NCBI Conserved Domain Database in the 'cdd' pipeline.

    1. Follow the installation instructions at the link above according to your operating system.
    2. The cdd tool assumes the binaries are installed at "/usr/bin/rpsblast+".


NCBI Conserved Domain Database
******************************

:cdd:`NCBI Conserved Domain Database <>`

    Required only if conserved domains within genes need to be identified from the NCBI Conserved Domain Database in the 'cdd' pipeline.

    1. Download the compressed :cdd_le:`NCBI CDD <>`.
    2. Expand the archived file into a directory of CDD files.




``pdm_utils`` source code repository
************************************

Some ``pdm_utils`` tools may require non-Python data files that are not directly installed with the Python package. Instead, these files are available on the ``pdm_utils`` git repository, which can be accessed through :pdmutils:`GitHub <>`. The repository can be downloaded two ways:

    1. Using git on the command line::

        > git clone https://github.com/SEA-PHAGES/pdm_utils.git

    2. Manually through GitHub.
