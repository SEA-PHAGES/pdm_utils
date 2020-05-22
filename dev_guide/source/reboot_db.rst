Rebooting the production database
=================================


As changes and improvements are made a) to the database schema, b) how data is persisted, and c) how data is parsed from external sources, it may be necessary to re-import all genomes into the database, effectively 'rebooting' the database. This ensures all genomes have been evaluated, parsed, and persisted the identical way.


At any given time, genomes in the database may be at one of several different storage stages. The stage is based on several factors, and there is no single field that reflects the stage.


Create a base database to begin the reboot::

    python3 -m pdm_utils get_db Actinobacteriophage server -d -o ./

This provides a reliable starting point in case an error is encountered during the process.


Install the base database::

    python3 -m pdm_utils get_db Actinobacteriophage file ./Actinobacteriophage.sql


Reset import date for all genomes::

    mysql -u root -p Actinobacteriophage --execute "UPDATE phage SET DateLastModified = '1900-01-01';""


Download all data from PhagesDB, PECAAN, and GenBank. The '-f' force option will cause all date data in the MySQL database to be ignored. This resuls in all flat files on PhagesDB and on GenBank to be retrieved and all 'draft' genomes currently in the database to be re-downloaded from PECAAN. The '-gr' option will report which records could not be retrieved from GenBank::

    python3 -m pdm_utils get_data Actinobacteriophage -a -c ./ncbi_creds.txt -o ./ -gr


Note: re-downloading thousands of genomes can take a long time. The above step could be run separately with the '-d', '-f', and '-g' options.





#
