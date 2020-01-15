.. _freeze:

freeze: create a static database for long-term reference
========================================================

The primary database instance, Actinobacteriophage, contains all available actinobacteriophage data, and is routinely updated and modified. However, specific projects may require a version of the database that:

    1. contains a subset of genome data (such as all genomes that have NOT been auto-annotated),

    2. is no longer routinely modified/updated, and/or

    3. has an identifier distinct from other databases.


The ``pdm_utils freeze`` tool can create these “frozen” databases::

    > python3 -m pdm_utils freeze Actinobacteriophage ./

The argument 'Actinobacteriophage' indicates the name of the local database from which a frozen database will be generated. The './' argument indicates the working directory for manipulating files and storing output files.

The tool copies the database indicated in the first script argument, deletes all data pertaining to draft genomes (and thus retains all final and unknown genomes), and saves the new database with a unique identifier that indicates both the type of phage database and the number of genomes (e.g. Actinobacteriophage_<0123>, where <0123> is an integer). The tool creates a folder for this new database in the directory indicated by the second script argument, and creates three subfolders (current, backup, and updates_history) analogous to the folders regularly used to maintain the Actinobacteriophage instance. Since the new frozen database no longer contains draft genomes, gene phams need to be recomputed. As a result, the script does not export the new database, and it resets the database version to '0'. After re-computing gene phams, the 'update' tool can increment the version number, the 'export' tool can generate a new SQL file, and the 'push' tool can upload it to the server.

Different types of databases may need to be frozen. For instance, sometimes all actinobacteriophages (other than draft genomes) need to be retained. Other times, only the *Mycobacterium* phages need to be retained. As a result, the script prompts the user to choose an appropriate database name (e.g. Actinobacteriophage, Mycobacteriophage, etc.) before it appends the genome count to the database name. A customized database name can be provided if needed. However, this script does not yet provide the functionality of enabling the database administrator to indicate which types of phages should be retained (other than the annotation status), and this step currently needs to be completed with the mysql command line utility. It is also important to note that although this database is regarded as “frozen”, it is still able to be modified. Therefore, if minor errors are encountered in the database that need to be modified, the database can be adjusted in MySQL, the version number can be incremented, the database can be re-exported, and the updated database file can overwrite the older file on the server.
