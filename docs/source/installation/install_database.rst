.. _install_database:

MySQL database instance
=======================

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
