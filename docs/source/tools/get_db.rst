.. _getdb:

get_db: download and install a database
=======================================


The Hatfull lab server hosts the primary actinobacteriophage database instance, Actinobacteriophage, which is routinely updated with new genomics data, as well as databases that have been frozen for publication. These databases can be downloaded and installed using the ``pdm_utils get_db`` tool.

To download and install the current version of the Actinobacteriophage::

    > python3 -m pdm_utils get_db Actinobacteriophage ./ -a

The argument 'Actinobacteriophage' indicates the name of the database to download from the server, and it will install it under the same name in your local MySQL. The './' indicates the working directory where the file should be downloaded. The '-a' flag indicates that the database will be downloaded, installed, and then the file will be removed.

.. note::
    This will overwrite an existing Actinobacteriophage database, if it exists.


Databases can be downloaded and installed in two steps, which can be used to install a database under a new name:

    1. First download the database from the Hatfull server using the '-d' flag. This will save the database to your local computer as a SQL file (e.g. Actinobacteriophage.sql) without installing it in MySQL::

        > python3 -m pdm_utils get_db Actinobacteriophage ./ -d

    2. Next, indicate the new name of the database to be created (e.g. Actinobacteriophage_new), use the '-i' flag to install the database with this new name, and use the '-f' flag to indicate a local SQL file to use (e.g. the file that was just downloaded)::

        > python3 -m pdm_utils get_db Actinobacteriophage_new ./ -i -f Actinobacteriophage.sql
