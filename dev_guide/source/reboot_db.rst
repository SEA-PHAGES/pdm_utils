Rebooting the production database
=================================


As changes and improvements are made a) to the database schema, b) how data is persisted, and c) how data is parsed from external sources, it may be necessary to re-import all genomes into the database, effectively 'rebooting' the database. This ensures all genomes have been evaluated, parsed, and persisted the identical way.


At any given time, genomes in the database may be at one of several different storage stages. The stage is based on several factors, and there is no single field that reflects the stage.


Create a 'prep' database to begin the reboot::

    python3 -m pdm_utils get_db Actinobacteriophage server -d -o ./

Install the 'prep' database::

    python3 -m pdm_utils get_db Actinobacteriophage file ./Actinobacteriophage.sql


Perform a round of updates to ensure the database is in a clean state prior to rebooting (i.e. it contains the most up-to-date information for genomes)::

    python3 -m pdm_utils get_data Actinobacteriophage -a -c ./ncbi_creds.txt -o ./

Update and/or import each category of new data::

    python3 -m pdm_utils update Actinobacteriophage ...
    python3 -m pdm_utils import Actinobacteriophage ...


There is no need to run phamerate or find_domains pipelines, since that will be done later on. Update the version and export a copy of this database::

    python3 -m pdm_utils update Actinobacteriophage -v
    python3 -m pdm_utils export Actinobacteriophage sql

This is the 'base' database for all following steps. This provides a reliable starting point in case an error is encountered during the process.


It is easier to assess the reboot success if different categories of phages are processed separately. To do this, the get_data pipeline needs to be run on a derivative 'sub' database containing a sub-selection from the base database, so that a new staged directory of that category of phages is created. For each of the categories below, delete the indicated genomes from the base database, run the indicated get_data pipeline, then re-install the 'base' database.

When running the get_db pipeline, the '-fd' force option will cause all date data in the MySQL database to be ignored. For 'draft' genomes, this results in all 'draft' genomes in the database to be re-auto-annotated with PECAAN (although with any new genomes available on PhagesDB. For 'non-draft' genomes, also results in all flat files in PhagesDB and GenBank to be retrieved regardless of import date. The '-gr' option will report which records could/could not not be retrieved from GenBank.



Category 1: Preparing draft genomes
-----------------------------------

Get auto-annotated genomes from PECAAN for new genomes in PhagesDB and all 'draft' genomes currently in the database::

    python3 -m pdm_utils get_data Actinobacteriophage -fd -d -o ./

Reinstall the 'base' database.




Category 2: Preparing SEA-PHAGES genomes that are auto-updated from GenBank and that have valid Genbank flat files
------------------------------------------------------------------------------------------------------

Delete all genomes that are not SEA-PHAGES and not set to be auto-updated::


    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE AnnotationAuthor != 1;
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE RetrieveRecord != 1;


Now download all data from GenBank::

    python3 -m pdm_utils get_data Actinobacteriophage -fd -g -c ./ncbi_creds.txt -o ./ -gr


Review the genbank_results file to assess which accessions are not active/valid. From this list, construct an update ticket table (aa1_acc.csv) using the subset of phages with active Accessions, and change the DateLastModiifed to an unused value (such as '3000-01-01').




Category 3: Preparing SEA-PHAGES genomes do not have valid Genbank flat files
-----------------------------------------------------------------

Using the list provided above, generate a list of SEA-PHAGES genomes that should be retrieved from PhagesDB::

    python3 -m pdm_utils update Actinobacteriophage -f aa1_acc.csv
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE AnnotationAuthor != 1;
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE DateLastModified = '3000-01-01';


Now run get_db to retrieve as many records from PhagesDB that are available::

    python3 -m pdm_utils get_data Actinobacteriophage -fd -f -o ./



Category 4: Prepare subsets of SEA-PHAGES genomes that are in temporary intermediate stages
-------------------------------------------------------------------------------

Some phages may have valid accessions, but are not set to be auto-updated yet, usually because there is an error in the publicly available flat file. The genome should then be retrieved from either PhagesDB as a new final, or from PECAAN as a new draft. If that list isn't readily available, it will become apparent when the error-prone database is imported and errors are encountered.

For these genomes, manually stage them in the genomes folder, manually create the import ticket table, and set the evaluation mode appropriately.

Also prepare a SQL script to delete this list of genomes from the database (delete_phages.sql).






Category 5: Prepare non-SEA-PHAGES genomes (assuming that none are set to be auto-updated) that do not have valid Genbank flat files (such as prophage genomes manually prepared) or that have been manually edited:
-------------------------------------------------------------------------------

This is a manually process, but it is a pretty static list. Currently, this list of phages is:

Bfk20, E3 = flat files have been manually edited.
ISF9, mu16, phiSAV, Shyg, SPB78, Sros11, StrepC, VWB = genomes that have been manually re-oriented.





Category 6: Preparing non-SEA-PHAGES genomes (assuming that none are set to be auto-updated) that have valid Genbank flat files
-------------------------------------------------------------------------------

Currently, for the non-SEA-PHAGES genomes that have flat files that have been manually edited for import, they have a comment in phage.Notes. Remove all genomes except for non-SEA-PHAGES that have an accession and that have not been manually edited::

    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE AnnotationAuthor != 0;"
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE Accession != '';"
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE Notes is not NULL;"

Alternatively, a pre-defined list could be provided::

    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE AnnotationAuthor != 0;"
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE Accession != '';"
    mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE PhageID in ('Bfk20', 'E3', 'ISF9', 'mu16', 'phiSAV', 'Shyg', 'SPB78', 'Sros11', 'StrepC', 'VWB');"


Now use get_data to retrieve those GenBank flat files::

    python3 -m pdm_utils get_data Actinobacteriophage -g -c ./ncbi_creds.txt -o ./ -gr







######

#

# To reboot the database

.. Remove all 'draft' genomes, since the import pipeline will log an error if a 'draft' genome is replacing a 'draft' genome::
..
..     mysql -u root -p Actinobacteriophage --execute "DELETE FROM phage WHERE Status = 'draft';"


.. Reset the import date for all genomes::
..
..     mysql -u root -p Actinobacteriophage --execute "UPDATE phage SET DateLastModified = '1900-01-01';"
