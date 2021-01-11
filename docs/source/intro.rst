Introduction
============

``pdm_utils`` is a Python package designed in combination with a pre-defined MySQL database schema in order to facilitate the creation, management, and manipulation of phage genomics databases in the :seaphages:`SEA-PHAGES program <>`. The package is directly connected to the structure of the MySQL database, and it provides several types of functionality:

    1. :ref:`Python library <package_overview>` including:

        a. Classes to store/parse phage genomes, interact with a local MySQL genomics database, and manage the process of making database changes.

        b. Functions and methods to manipulate those classes as well as interact with several databases and servers, including PhagesDB, GenBank, PECAAN, and MySQL.

    2. A command line :ref:`toolkit <toolkit>` to process data and maintain a phage genomics database.

``pdm_utils`` is useful for:

    1. The Hatfull lab to maintain MySQL phage genomics databases in the SEA-PHAGES program (:ref:`current pipeline <actinobacteriophagedb>`).

    2. Researchers to evaluate new genome annotations (:ref:`flat file QC <flatfileqc>`).

    3. Researchers to directly access and retrieve phage genomics data from any compatible MySQL database (:ref:`tutorial <library_tutorial>`).

    4. Researchers to create :ref:`custom <custom_dbs>`  MySQL phage genomics databases.

    5. Developers to build downstream data analysis tools (:ref:`tutorial <library_tutorial>`).


``pdm_utils`` source code
_________________________

This project is maintained using git and is available on :pdmutils:`GitHub <>`.
