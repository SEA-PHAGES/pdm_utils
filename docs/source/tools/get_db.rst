.. _getdb:

get_db: download and install a database
=======================================


The Hatfull lab server hosts the primary actinobacteriophage database instance, Actino_Draft, which is routinely updated with new genomics data, as well as databases that have been frozen for publication. These databases can be downloaded and installed using the ``pdm_utils get_db`` tool.

To view and interactively select the list of available databases on the Hatfull lab server::

    > python3 -m pdm_utils get_db server

The ``pdm_utils get_db`` interactive tool functions as a faux command line interface that allows you to view the database packages at the specified url and navigate through subdirectories.  Navigation through subdirectories is modelled after common CLI commands::

    root@http://databases.hatfull.org/::$ ls
            Abscessus_prophages             PhageDatabase 
            Actino_Draft_v6                 Actino_Draft
            Actino_prophage                 Published/

    root@http://databases.hatfull.org/::$ cd Published
    root@http://databases.hatfull.org/::/Published$ ls 
            Actinobacteriophage_1060        Actinobacteriophage_1321
            Actinobacteriophage_2422        Actinobacteriophage_554
            Actinobacteriophage_685         Actinobacteriophage_692
            ../

With the ``pdm_utils get_db`` interactive tool you can see a description of a database package that is at available at the specified url::

    root@http://databases.hatfull.org/::$ desc PhageDatabase 
    Name:
            PhageDatabase 
    Date:
            202x-xx-xx
    Description:
            This database contains sequenced phages from the John Doe lab.

A database of choice can be selected using the ``pdm_utils get_db`` interactive tool with the following::
    
    root@http://databases.hatfull.org/::$ select PhageDatabase

To download and install the current version of a database, like the Actino_Draft database, without the interactive tool::

    > python3 -m pdm_utils get_db server -db Actino_Draft 

    > python3 -m pdm_utils get_db server --database Actino_Draft 

The -db argument 'Actino_Draft' indicates the name of the database to download from the server, and it will install it under the same name in your local MySQL. The database will be downloaded, installed, and then the file will be removed.

.. note::
    This will overwrite an existing Actino_Draft database, if it exists.


To download and install a database from a non-standard server, specify the URL::

    > python3 -m pdm_utils get_db server -db PhageDatabaseName -u http://custom/server/website 

    > python3 -m pdm_utils get_db server -db PhageDatabaseName --url http://custom/server/website

The ``pdm_utils get_db`` tool checks your local database version against the specified server database package and will not download if the local database version is equal to or higher than the database package to prevent redundancies and/or loss of data.  To ignore this check::

    > python3 -m pdm_utils get_db server -db Actino_Draft -fp

    > python3 -m pdm_utils get_db server -db Actino_Draft --force_pull

Databases can be downloaded and installed in two steps, which can be used to install a database under a new name:

    1. First download the database sql and version files from the Hatfull server using the '-d' and '-v' flags. This will save the database to your local computer as a SQL file (e.g. Actino_Draft.sql) without installing it in MySQL. Also specify where the file should be downloaded using the '-o' flag (if omitted, the default is the /tmp/ directory)::

        > python3 -m pdm_utils get_db server -db Actino_Draft -o ./ -d -v


    2. Next, indicate the new name of the database to be created (e.g. NewDB), indicate that a local file will be used with 'file' , and indicate the path to the downloaded SQL file::

        > python3 -m pdm_utils get_db NewDB file ./downloaded_db/Actino_Draft.sql


Use of a :ref:`config_file` can automate the pipeline so that user input such as MySQL credentials or server URL is not needed::

    > python3 -m pdm_utils get_db server -db Actino_Draft -o ./ -d -v -c config_file.txt
