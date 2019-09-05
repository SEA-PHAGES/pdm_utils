Data export
===========

After all new data is imported into the database, new CDD data is retrieved, and new genes are grouped into phams, the database needs to be uploaded to a public server so that it is accessible by other tools and applications [Figure A-1]. The export_database.py script accomplishes this, which requires four arguments [Table A-12]. There are several steps in this script, and each can be performed separately.

Step 1: export the database
___________________________

The updated PhameratorDB database is stored within the MySQL directory on the database administrator’s local computer. In order to provide the updated database to other users, the database needs to be exported from MySQL into a single file that can be easily uploaded to a server (e.g. Actino_Draft.sql). In PhameratorDB, the database version is tracked as an integer in the Version field of the Version table. During export, the script can also update the version (in which the integer is incremented by 1) and create a version file (e.g. Actino_Draft.version). This text file contains an integer that corresponds to the database version integer. The database and version files are stored in the current directory indicated in the second script argument. Downstream applications rely on the name of the database to remain unchanged, so filenames in this directory do not change. Since the newly-exported files represent the most up-to-date, current version of the PhameratorDB instance, every round of exports overwrites the two files in this folder. As a result, a backup copy of the new database is also stored in the backup directory indicated in the third script argument, but the filename is modified to indicate the version (e.g. Actino_Draft_v[0123].sql, where [0123] is an integer). These backup copies can be re-imported into MySQL if needed to return to an older version of the database.

Step 2: Query the new database for updated genome and gene data
_______________________________________________________________

Managing PhameratorDB often requires quick reference to data in the Phage and Gene tables. This step queries the new database and outputs commonly-referenced fields from the Phage table (such as Cluster2, Subcluster2, HostStrain, DateLastModified, etc.) into a csv-formatted table with a filename that contains the date of retrieval, the database name, and the database version (e.g. 20180101_Actino_Draft_v[0123]_genomes.csv, where [0123] is an integer). It also outputs commonly-referenced fields from the Gene table (such as GeneID, Gene Description, etc.) into a csv-formatted table with a similarly structured filename (e.g. 20180101_Actino_Draft_v[0123]_genes.csv, where [0123] is an integer). The files are stored in the output directory indicated in the fourth script argument. This step is merely for convenience and does not impact any other process of database administration.

Step 3: Upload database to the server
_____________________________________

Once the database and version files are created, they can be uploaded to the `Hatfull lab’s public server <http://phamerator.webfactional.com/databases_Hatfull>`_. The script uploads the two files associated with the database that is indicated in the first script argument from the current directory that is indicated in the second script argument.

Each step in this script is independent of the others. For instance, data can be queried from a particular database without exporting it or incrementing the version number, and a database file that has been previously exported can be uploaded to the server without having to export it a second time.




Managing the record of tickets
______________________________

The successful import tickets from the import script and the update tickets from any of the update scripts represent a record of the diverse types of changes made to the database during the round of updates. All tickets and associated files used in each round of updates can be stored in an updates_history folder in the same directory as the current and backup database folders. Direct changes to the database can be made in MySQL software without using scripts in the k_phamerate repository, but since this does not generate a record of the changes, it is not recommended as common practice. However, descriptions of these types of changes can also be manually recorded in text files and stored in the updates_history folder as well.
