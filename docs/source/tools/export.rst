.. _export:

export: export data from a database
===================================

This tool is used to export data from a MySQL database in a variety of formats, including:

    1. A SQL file that represents the entire database.
    2. CSV-formatted tables for selected tables.
    3. Formatted files for selected genomes, such as:

        a. GenBank-formatted flat files
        b. Fasta-formatted files


SQL File
________

Export the entire database as a SQL file::

    > python3 -m pdm_utils export sql Actinobacteriophage

In order to distribute an updated database to end-users, the database needs to be exported from MySQL into a single file that can be easily uploaded to a server (e.g. Actinobacteriophage.sql). The database version is tracked as an integer in the Version field of the *version* table. A version file is also generated (e.g. Actinobacteriophage.version), which is a text file that contains a single integer corresponding to the database version.



CSV File
________

Export specific tables from the database into CSV-formatted files::

    > python3 -m pdm_utils export csv Actinobacteriophage ...







GenBank-formatted Flat File
___________________________

Export genomes into GenBank-formatted flat files::

    > python3 -m pdm_utils export gb Actinobacteriophage ...
