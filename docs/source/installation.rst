Installation
============


The ``pdm_utils`` package is written in Python3 and can be installed on MacOS and Linux platforms.

Dependencies
____________


There are also several third-party dependencies that need to be installed locally for certain ``pdm_utils`` tools.

1. MySQL
********

:mysql:`MySQL Community Server 5.7 <>`

    This is required for practically all ``pdm_utils`` tools. Below is a brief summary of the installation steps. Refer to the official MySQL documentation for more details.

MacOS installation
++++++++++++++++++

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
+++++++++++++++++++

# TODO Add description





2. Phamerator database instance
*******************************

Many ``pdm_utils`` modules and pipelines require access to a specifically structured MySQL database that can be used by the Phamerator GUI.


The Actino_Draft database
+++++++++++++++++++++++++

Installing the primary, actinobacteriophage, database instance (Actino_Draft) for the first time can be performed as follows.

    1. Open a Terminal window.
    2. Create an empty database (enter your password when prompted)::

        > mysql -u root -p --execute "CREATE DATABASE Actino_Draft"

    3. Download the current version of the database from the Hatfull lab server::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actino_Draft.sql > ./Actino_Draft.sql

    4. Import the database into MySQL (enter your password when prompted)::

        > mysql -u root -p Actino_Draft < Actino_Draft.sql

    5. The database is now set up for use.

This database is always being updated. Keeping the local database up-to-date can be performed as follows.

    1. Log in to MySQL (enter your password when prompted)::

        > mysql -u root -p

    2. Execute the following query to get the current version::

        mysql> SELECT Version FROM version;
        mysql> exit

    3. Download the current version file from the Hatfull lab server::

        > curl http://phamerator.webfactional.com/databases_Hatfull/Actino_Draft.version > ./Actino_Draft.version

    4. If the current version on the server is different from the version in the local MySQL database, there is a new database available on the server. Repeat steps 3-4 listed above in the 'installing a new Actino_Draft database' section.


Frozen Phamerator databases
+++++++++++++++++++++++++++

Typically, different versions, or instances, of the Phamerator database are created for specific studies/publications. The unique name of the database is normally published in the Materials and Methods. To download this database, follow the same steps as described above, substituting the frozen database name for Actino_Draft.


3. MMSeqs
*********

:mmseqs:`MMSeqs <>`

    Required only if gene phamilies need to be identified using MMSeqs and the mmseqs_phamerate pipeline.

# TODO add installation instructions.


4. NCBI Blast+ toolkit
**********************

:blastplus:`NCBI blast+ <>`

    Required only if conserved domains within genes need to be identified from the NCBI Conserved Domain Database and the cdd pipeline.

# TODO add installation instructions.




5. NCBI Conserved Domain Database
*********************************

:cdd:`NCBI Conserved Domain Database <>`

    Required only if conserved domains within genes need to be identified using the cdd pipeline.

# TODO add installation instructions.





6. Python dependencies in MacOS or Ubuntu
*****************************************

There are several third-party python packages:

    - :biopython:`Biopython <>`
    - :pymysql:`pymysql <>`
    - :paramiko:`paramiko <>`
    - tabulate

Some of the python dependencies themselves have python or binary dependencies. Although these dependencies can be manually installed, it can be tricky to do so. Instead, the Conda environment manager is a simple, automated alternative.

    1. Install Conda locally through the :anaconda:`Anaconda <>` package.
    2. After installing Conda, create an environment to be able to install and use ``pdm_utils``::

        > conda create --name pdm_utils-user python pip biopython pymysql paramiko tabulate
        > source activate pdm_utils-user



The ``pdm_utils`` package
_________________________


    1. ``pdm_utils`` cannot yet be installed through Conda. With the Conda environment activated, execute the following command::

        (pdm_utils-user)> pip install pdm_utils

    2. Update the version::

        (pdm_utils-user)> pip install --upgrade pdm_utils
