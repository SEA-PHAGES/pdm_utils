.. _installation:

Installation
============


The ``pdm_utils`` package is written in Python3 and can be installed on MacOS and Linux platforms. There are several third-party general dependencies that need to be installed in advance, and there are several dependencies specific to a subset of ``pdm_utils`` tools or pipelines that are only needed for specific goals. Below is a general step-by-step guide to installing all dependencies and the ``pdm_utils`` package.


1. MySQL
________

:mysql:`MySQL Community Server 5.7 <>`

    This is required for practically all ``pdm_utils`` tools. Below is a brief summary of the installation steps. Refer to the official MySQL documentation for more details.

MacOS installation
******************

Installing MySQL on MacOS can be tricky.

    1. Download the MySQL Community Server 5.7 dmg file from the MySQL website.
    2. Follow the standard installation instructions. MySQL can be installed at the default location.
    3. Make note of the temporary password provided during the installation process.

    .. warning::
         Be sure to record the temporary password that is generated!!! Once MySQL is installed, this password can only be used ONCE to login.

    4. Open Terminal.
    5. Log in to MySQL as the root user::

        > mysql -u root -p

    6. Enter the temporary password when prompted::

        Enter password: <temporary password>

    7. At the mysql prompt, change the password for the root user::

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

Ubuntu installation
*******************

# TODO Add description



2. Python dependencies
______________________

There are several third-party python packages required:

    - :biopython:`Biopython <>`
    - :pymysql:`pymysql <>`
    - :paramiko:`paramiko <>`
    - tabulate

Some of them also have python or binary dependencies. Manual installation of these dependencies can be tricky, but the Conda environment manager is a simple, automated alternative. First install Conda, then use Conda to install all python dependencies:

    1. Install Conda locally through the :anaconda:`Anaconda <>` package. Follow the manufacturer's installation guide for MacOS or Ubuntu.

    2. After installing Conda, create an environment to be able to install and use ``pdm_utils`` (the example below creates a Conda environment named 'pdm_utils', but it can be named anything)::

        > conda create --name pdm_utils python pip biopython pymysql paramiko tabulate
        > conda activate pdm_utils
        (pdm_utils)>

The command line prompt will now include '(pdm_utils)', indicating it is operating within this environment.


3. The ``pdm_utils`` package
____________________________

Once MySQL and the Conda environment are installed, ``pdm_utils`` can be easily installed:

    1. Execute the following command::

        (pdm_utils)> pip install pdm_utils

    2. The package is routinely updated, and the most recent version can be retrieved::

        (pdm_utils)> pip install --upgrade pdm_utils


4. Phamerator database instance
_______________________________

Many ``pdm_utils`` modules and pipelines require access to a specifically structured MySQL database that can be used by the Phamerator GUI.

The primary database instance that reflects the most up-to-date actinobacteriophage genomics data is the 'Actino_Draft' database. Typically, different versions, or instances, of the Phamerator database are created ('frozen') for specific studies/publications. The unique name of the database is normally published in the Materials and Methods.

The ``pdm_utils`` 'get_db' installation management tool can be used to retrieve, install, and update these databases, or any custom MySQL database that is compliant with the Phamerator database schema, from a local file or from the Hatfull lab server (:ref:`getdb`).

Alternatively, databases can be manually downloaded and installed, as described below (using Actino_Draft as an example):

Manual installation
*******************

    1. Open a Terminal window.
    2. Create an empty database (enter your password when prompted)::

        > mysql -u root -p --execute "CREATE DATABASE Actino_Draft"

    3. Retrieve the current version of the database::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actino_Draft.sql > ./Actino_Draft.sql

    4. Import the database into MySQL (enter your password when prompted)::

        > mysql -u root -p Actino_Draft < Actino_Draft.sql


Manual update
*************

    1. Log in to MySQL (enter your password when prompted)::

        > mysql -u root -p

    2. Execute the following query to get the current version::

        mysql> SELECT Version FROM version;
        mysql> exit

    3. Download the current version file from the Hatfull lab server::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actino_Draft.version > ./Actino_Draft.version

    4. If the current version on the server is different from the version in the local MySQL database, there is a new database available on the server. Repeat steps 3-4 listed above in the 'installing a new Actino_Draft database' section.









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

Some ``pdm_utils`` tools, such as the 'convert' tool, require non-Python data files that are not directly installed with the Python package. Instead, these files are available on the ``pdm_utils`` git repository, which can be accessed through :pdmutils:`GitHub <>`. The repository can be downloaded two ways:

    1. Using git on the command line::

        > git clone https://github.com/SEA-PHAGES/pdm_utils.git

    2. Manually through GitHub.
