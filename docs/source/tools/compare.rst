.. _compare:

compare: compare data between databases
=======================================

In the SEA-PHAGES program, genomics data may be stored in three separate databases (a local MySQL database, PhagesDB, and GenBank). As a result, data inconsistencies may arise without a mechanism to compare and synchronize different data sources. Although the import tool implements many QC checks for this purpose, it only ensures consistency for the specific genomes being imported and it only cross-checks databases as directed by the import ticket. In order to ensure comprehensive consistency, the ``pdm_utils compare`` tool can perform an all-against-all data assessment::

    > python3 -m pdm_utils compare Actino_Draft -p -g -c config_file.txt

The argument 'Actino_Draft' indicates the name of the local MySQL database that serves as the central data source for the comparison. The '-p' and '-g' arguments indicate that data from PhagesDB and GenBank should be compared to the MySQL database as well as to each other (an all-against-all-against-all comparison).

    ..note:
        It does not directly compare PhagesDB and GenBank data unless it also compares MySQL database data.

The '-c config_file.txt' argument indicates a local file containing login information for accessing MySQL and NCBI.

Selected subsets of data can be compared using the filtering argument. For instance, to only compare phages that infect Mycobacterium and that are 'final' status::

    > python3 -m pdm_utils compare Actino_Draft -p -g -c config_file.txt -f "phage.HostGenus=Mycobacterium AND phage.Status=final"

The ``compare`` tool retrieves data stored in the *phage* and *gene* tables. PhagesDB data for all sequenced phages are retrieved from: http://phagesdb.org/api/sequenced_phages/. All GenBank records are retrieved using accession numbers stored in the Accession field of the *phage* table. All data is matched between the local MySQL database and PhagesDB using the PhageID field in the *phage* and the phage name in PhagesDB. All data is matched between the local MySQL database and GenBank using the Accession field in the *phage* table. After retrieving and matching data from all databases, the script compares the genome data (e.g. phage name, host strain, genome sequence, etc.) and gene data (e.g. locus tags, coordinates, descriptions, etc.), and generates several results files. Additionally, the script can output genomes retrieved from all three databases for future reference and analysis if selected by the user.
