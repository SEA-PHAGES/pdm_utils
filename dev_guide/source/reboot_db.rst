Rebooting the production database
=================================


As changes and improvements are made a) to the database schema, b) how data is persisted, and c) how data is parsed and evaluated from external sources, it may be necessary to re-import all genomes into the database, effectively 'rebooting' the database. This ensures all genomes have been evaluated, parsed, and persisted the identical way.


Genome categories
*****************

At any given time, genomes in the database may be at one of several different stages. The stage is based on several factors, and there is no single field that reflects the stage. Currently there are 8 different stages, but these are not comprehensive, and in the future there could be more or fewer stages:

SEA-PHAGES genomes

1. draft status with auto-update from GenBank, so auto-annotate from PECAAN (easy to stage, easy to import)
2. draft status  without auto-update from GenBank, so auto-annotate from PECAAN but change RetrieveRecord = 0 in ticket (tricky to stage, easy to import)
3. non-draft status  without auto-update from GenBank, so retrieve data from PhagesDB (tricky to stage, tricky to import)
4. non-draft status with no GenBank accession, so retrieve data from PhagesDB (easy to stage, tricky to import)
5. non-draft status with valid GenBank accessions, so retrieve data from GenBank (easy to stage, easy to import)
6. non-draft status with invalid GenBank accessions (only determined by trying to download from GenBank), so retrieve data from PhagesDB (tricky to stage, tricky to import)

Non-SEA-PHAGES genomes

7. no manual genome modifications, so retrieve data from GenBank (tricky to stage, tricky to import)
8. with manual genome modifications, so need to re-import from local files (tricky to stage, tricky to import)





Preparing the database for staging the genomes
**********************************************

Perform a full round of updates (get data, update, import, phamerate, etc.) to ensure that the database contains the most up-to-date data.

Download the most up-to-date database from the server. This is the 'base' database::

    python3 -m pdm_utils get_db Actino_Draft server -d -o ./
    mysql -u root -p --execute "CREATE DATABASE staging"
    python3 -m pdm_utils get_db staging file ./Actino_Draft.sql

The base database provides a reliable starting point in case an error is encountered during the staging and import steps.



Staging genomes for import
**************************

It is easier to assess the reboot success if different categories of phages are processed separately (when problems are encountered, knowing the phage category helps to interpret the errors). To do this, the get_data pipeline needs to be run on a derivative 'sub' database containing a sub-selection from the base database, which will then generate a staged directory only of that category of phages.

For each genome category:

1. query the base database to know how many genomes are expected
2. delete the indicated unwanted genomes from the base database
3. run the indicated get_data pipeline
4. verify all genomes were retrieved and staged
5. re-install the 'base' database.

When running the get_db pipeline, the '-fd' force option is needed. In general, get_data relies on several data flags to determine whether the genome data should be retrieved or not. The '-fd' option will ignore those flags. For example, when getting 'draft' genomes, '-fd' results in all 'draft' status genomes in the database to be re-auto-annotated with PECAAN (in addition to any new genomes available on PhagesDB).


Category 1 and 2: Preparing draft genomes
-----------------------------------------

Check how many to expect: category 1, auto-annotate from PECAAN::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 1 AND Status = 'draft' AND RetrieveRecord = 1;

Check how many to expect: category 2, auto-annotate from PECAAN but manually set RetrieveRecord = 0 in ticket::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 1 AND Status = 'draft' AND RetrieveRecord = 0;


Get auto-annotated genomes from PECAAN for new genomes in PhagesDB and all 'draft' genomes currently in the database::

    python3 -m pdm_utils get_data staging -fd -d -o ./

Identify the genomes that are frozen at RetrieveRecord = 0, and manually set RetrieveRecord = 0 in the import tickets for these genomes::

    SELECT PhageID FROM phage WHERE AnnotationAuthor = 1 and Status = 'draft' and RetrieveRecord = 0;



Category 3: SEA-PHAGES genomes with manual annotations but not auto-updated
---------------------------------------------------------------------------

Check how many to expect: non-draft phages without auto-update from GenBank (need to import data from PhagesDB or GenBank)::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 1 AND Status != 'draft' AND RetrieveRecord = 0;


Delete all genomes that are not SEA-PHAGES, are 'draft' status, and are set to be auto-updated::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE AnnotationAuthor != 1;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Status = 'draft';"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE RetrieveRecord != 0;"

Some of these genomes may have valid accessions, so first try to download from GenBank::

    python3 -m pdm_utils get_data staging -fd -g -c ./ncbi_creds.txt -o ./ -gr

If there were any genomes not retrieved from GenBank, then delete genomes that were downloaded, and then download the rest of the data from PhagesDB::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE PhageID IN (<list of PhageIDs);"
    python3 -m pdm_utils get_data staging -fd -f -o ./

In the import table, set eval_mode = 'auto'.



Category 4: SEA-PHAGES genomes with manual annotations but not yet in GenBank
-----------------------------------------------------------------------------

Check how many to expect: phages with no GenBank accession (need to import data from PhagesDB)::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 1 AND Status != 'draft' AND RetrieveRecord = 1 AND Accession = '';


Delete all genomes that are not SEA-PHAGES, are 'draft' status, and are set to be auto-updated::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE AnnotationAuthor != 1;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Status = 'draft';"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE RetrieveRecord != 1;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Accession != '';"

Now download all data from PhagesDB::

    python3 -m pdm_utils get_data staging -fd -f -o ./



Category 5 and 6: SEA-PHAGES genomes that are auto-updated from GenBank and have valid or invalid accessions
------------------------------------------------------------------------------------------------------------

Check how many to expect: phages with valid/invalid GenBank accessions (try to download from GenBank)::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 1 AND Status != 'draft' AND RetrieveRecord = 1 AND Accession != '';


Delete all genomes that are not SEA-PHAGES and not set to be auto-updated::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE AnnotationAuthor != 1;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Status = 'draft';"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE RetrieveRecord != 1;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Accession = '';"

Now download all data from GenBank::

    python3 -m pdm_utils get_data staging -fd -g -c ./ncbi_creds.txt -o ./ -gr

Review the genbank_results file to assess which accessions are not active/valid. From this list, construct an update ticket table (valid_accessions.csv) using the subset of phages with active Accessions, and change the DateLastModified to an invalid value (such as '3000-01-01'). Then update the database::

    python3 -m pdm_utils update staging -f valid_accessions.csv

Check to confirm::

    SELECT COUNT(*) FROM phage WHERE DateLastModified = '3000-01-01';

Remove the genomes with active accessions::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE DateLastModified = '3000-01-01';"


Now run get_db to retrieve as many records from PhagesDB that are available::

    python3 -m pdm_utils get_data staging -fd -f -o ./

In the import table, set eval_mode = 'auto'.



Category 7: Non-modified non-SEA-PHAGES genomes with valid Genbank accessions
-----------------------------------------------------------------------------

Currently, for the non-SEA-PHAGES genomes that have flat files that have been manually edited for import, they have a comment in phage.Notes.

Check how many to expect: non-SEA-PHAGES genomes with no manual genome modifications (can retrieve from GenBank)::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 0 AND Accession != '' AND Notes IS NULL;

Remove all genomes except for non-SEA-PHAGES that have an accession and that have not been manually edited::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE AnnotationAuthor != 0;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Accession = '';"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Notes is not NULL;"

Alternatively, a pre-defined list could be provided::

    mysql -u root -p staging --execute "DELETE FROM phage WHERE AnnotationAuthor != 0;"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE Accession != '';"
    mysql -u root -p staging --execute "DELETE FROM phage WHERE PhageID in ('Bfk20', 'E3', 'ISF9', 'mu16', 'phiSAV', 'Shyg', 'SPB78', 'Sros11', 'StrepC', 'VWB');"


Now use get_data to retrieve those GenBank flat files::

    python3 -m pdm_utils get_data staging -fd -g -c ./ncbi_creds.txt -o ./ -gr





Category 8: Manually-modified non-SEA-PHAGES genomes with/without valid Genbank accessions
------------------------------------------------------------------------------------------

Check how many to expect: non-SEA-PHAGES genomes with manual genome modifications (need to re-import from local files)::

    SELECT COUNT(*) FROM phage WHERE AnnotationAuthor = 0 AND (Accession = '' OR Notes IS NOT NULL);

This is a manual process, but it is a pretty static list. To identify the most up-to-date list::

    SELECT PhageID FROM phage WHERE AnnotationAuthor = 0 and (Accession = '' or Notes != '')

Currently, this list of phages is:

Bfk20, E3 = flat files have been manually edited.
ISF9, mu16, phiSAV, Shyg, SPB78, Sros11, StrepC, VWB = genomes that have been manually re-oriented.

Manually prepare these flat files, stage in a genome folder, and create the import table.

Manually create an update table (phage_notes_update_table.csv) that will be used to add a description of how the genome has been manually modified.


Confirm all genomes have been staged
------------------------------------

Sum all staged genomes from each Categories, and compare to the total number of genomes in the database. If the totals don't match, then some genomes fell through the cracks. Review each category and determine which one(s) don't have the total number of expected phages.



Reboot the database
*******************

Import all genomes
------------------

Prepare the 'reboot' database::

    mysql -u root -p --execute "CREATE DATABASE reboot;"
    python3 -m pdm_utils get_db reboot file ./Actino_Draft.sql

Remove all 'draft' genomes, since the import pipeline will log an error if a 'draft' genome is replacing a 'draft' genome::

    mysql -u root -p reboot --execute "DELETE FROM phage WHERE Status = 'draft';"


Reset the import date for all genomes, else errors will be logged::

    mysql -u root -p reboot --execute "UPDATE phage SET DateLastModified = '1900-01-01';"


Confirm the database configuration::

    SELECT COUNT(*) FROM phage WHERE DateLastModified != '1900-01-01';
    SELECT COUNT(*) FROM phage WHERE Status = 'draft';


For each Category, run import non-interactively::

    python3 -m pdm_utils import reboot ./genomes ./import_table.csv -p -o ./


For some Categories, all genomes will be successfull, while for others (such as Categories 4, 5 and 8), some genomes will fail. Many failed genomes can pass if processed interactively and the appropriate warnings ignored::

    python3 -m pdm_utils import reboot ./genomes ./import_table.csv -p -o ./ -i

All genomes that fail the second import attempt will need to be modified.

For Category 8 genomes, update phage.Notes with how genomes have been modified::

    python3 -m pdm_utils update reboot -f phage_notes_update_table.csv



Review the final database configuration
---------------------------------------

Confirm that all genomes have been re-imported::

    SELECT COUNT(*) FROM phage;
    SELECT COUNT(*) FROM phage WHERE DateLastModified != '1900-01-01';


Confirm that all manually-modified non-SEA-PHAGES genomes have Notes::

    SELECT PhageID, Notes from phage where Notes is not NULL;


Confirm that phages with alternative spellings in GenBank files have the correct spelling::


    SELECT PhageID, Name FROM phage WHERE PhageID in ('pZL12', 'LeBron', 'BPBiebs31', 'CapnMurica', 'Fionnbharth');
    SELECT PhageID, Name FROM phage where Status != 'draft' and PhageID != Name;


Polish the database
-------------------

After all checks pass, proceed with polishing the database by running phamerate and find_domains pipelines, increment the version. Rename the rebooted database to 'Actino_Draft' if needed, and push to server. Copy to 'Actino_Draft' database and convert to the current downgrade schema version.
