Export data from database
=========================

This tool is used to export data from a Phamerator database in a variety of formats, including:

    1. A SQL file that represents the entire database.
    2. CSV-formatted tables for selected tables.
    3. Formatted files for selected genomes, such as:

        a. GenBank-formatted flat files
        b. Fasta-formatted files


SQL File
________

Export the entire database as a SQL file::

    > python3 -m pdm_utils export sql Actino_Draft ...


CSV File
________

Export specific tables from the database::

    > python3 -m pdm_utils export csv Actino_Draft ...


GenBank-formatted Flat File
___________________________

Export genomes into GenBank-formatted flat files::

    > python3 -m pdm_utils export gb Actino_Draft ...
