Introduction
============

``pdm_utils`` is a Python package designed to facilitate the creation, management, and manipulation of MySQL phage genomics databases for :webpham:`Phamerator <>` in the :seaphages:`SEA-PHAGES program <>`.

``pdm_utils`` contains several types of tools:

    1. Python library including:

        a. Classes to store/parse phage genome data from a variety of sources, interact with a MySQL database, and manage instructions to add, remove, and update genomes.

        b. Functions and methods to manipulate those classes as well as interact with several databases and servers, including PhagesDB, GenBank, PECAAN, and MySQL.

    2. Data processing pipelines to maintain a Phamerator database.

``pdm_utils`` is useful for:

    1. Database administrators in the Hatfull lab to maintain the primary Phamerator database instance, Actino_Draft.

    2. Researchers to retrieve, access, and export data from a Phamerator database.

    3. Non-Hatfull lab database administrators to generate custom Phamerator databases.

    4. Developers building data analysis tools that rely on a Phamerator database.



``pdm_utils`` source code
_________________________

This project is maintained using git and is available on :pdmutils:`GitHub <>`.
