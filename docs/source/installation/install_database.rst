.. _install_database:

MySQL database instance
=======================

Many ``pdm_utils`` modules and pipelines require access to a specifically structured MySQL database.

The primary database instance that reflects the most up-to-date actinobacteriophage genomics data in the SEA-PHAGES program is the 'Actino_Draft' database. Typically, different versions, or instances, of the database are created ('frozen') for specific studies/publications. The unique name of the database is normally published in the Materials and Methods.

The ``pdm_utils get_db`` installation management tool can be used to retrieve, install, and update these databases, or any custom MySQL database that is compliant with the database schema, from a local file or from the Hatfull lab server (:ref:`get_db <getdb>`).

Alternatively, databases can be manually downloaded and installed, as described below (using Actino_Draft as an example):

pdm_utils installation
*******************

#. Open a Terminal window.
#. Start the Conda environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. Use the pdm_utils package to retrieve the Actino_Draft database (enter your username and password when prompted)::

    (pdm_utils)> python3 -m pdm_utils get_db server -u http://databases.hatfull.org/ -db Actino_Draft -v

