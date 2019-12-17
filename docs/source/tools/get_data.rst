.. _getdata:

get_data: get new data to import into the database
==================================================


New genomics data routinely becomes available for adding to the Actino_Draft database, including:

    1. Metadata pertaining to individual phages (such as host, cluster, subcluster, and accession)
    2. Newly-sequenced and auto-annotated 'draft' genomes
    3. New manually-annotated 'final' genomes
    4. Updated annotations from GenBank


These data can be automatically retrieved using the ``pdm_utils`` **get_data** tool::

    > python3 -m pdm_utils get_data Actino_Draft ./ -c ncbi_credentials.txt

The argument 'Actino_Draft' indicates the name of the database from which updates are determined. The './' indicates the working directory where the data should be downloaded. Similar to the 'get_gb_records' tool, retrieving updates from GenBank relies upon the NCBI E-utilities (using a Biopython wrapper), and NCBI requests that you provide information about yourself. The '-c' flag points to a simple text file containing your information (:ref:`ncbicreds`). If only certain types of updates are required, there are separate command line flags for each type of update.

Each type of data is retrieved and staged for import. A new folder is created that contains:

    1. a CSV-formatted import table listing each 'ticket' pertaining to the type of update
    2. a 'genomes' subdirectory containing flat files for import (if applicable)


Metadata updates
----------------


PhagesDB is the primary source for Cluster, Subcluster, Host, and Accession data for phages in the Actino_Draft database. **get_data** compares these data for each phage in the selected database to the corresponding data in PhagesDB, and creates a new update ticket for all discrepant data that need to be corrected in the Actino_Draft database. New metadata is retrieved from PhagesDB: :phagesdb:`sequenced phages <api/sequenced_phages>`. For each phage, PhagesDB stores both GenBank and RefSeq accession data, but only GenBank accession data (stored in the *genbank_accession* field) are stored in the Actino_Draft database.




New 'draft' auto-annotations
----------------------------

PhagesDB is the primary source for new genome sequences generated through the SEA-PHAGES program. Automatically-generated 'draft' annotations can be imported into a MySQL database for immediate reference. **get_data** identifies genomes present in PhagesDB that are not present in the selected MySQL database. For each new genome, a request is sent to the phage genomics tool :pecaan:`PECAAN <>` to generate auto-annotations with the URL: 'https://discoverdev.kbrinsgd.org/phameratoroutput/phage/<PhageID>' (where <PhageID> indicates the specific phage name of interest). PECAAN retrieves the new sequence from PhagesDB, automatically annotates genes, and returns a GenBank-formatted flat file. **get_data** stages the new file and a corresponding import ticket in an import table that are ready to be processed with the 'import' tool.


New 'final' annotations
-----------------------

The 'draft' annotations are eventually replaced with manually-generated 'final' annotations. The refined annotations are submitted by senior annotators to PhagesDB in GenBank-formatted flat files, which are are stored in the *qced_genbank_file* field with a timestamp stored in the *qced_genbank_file_date* field. Similar to metadata updates, **get_data** matches phage data in the selected MySQL database to the corresponding data in PhagesDB, and determines whether there is a new version of a 'final' annotation available for import. It reviews the date that the genome data was uploaded to PhagesDB, and if it is more recent than the date of the annotations stored in the selected MySQL database (indicated in the DateLastModified field of the *phage* table), it stages the new flat file from PhagesDB and a corresponding import ticket in an import table that are ready to be processed with the 'import' tool.


Updated GenBank annotations
---------------------------

The 'final' annotations are eventually submitted to GenBank with a unique accession number. The GenBank records are routinely updated with improved annotation data, so these records are subsequently re-imported into the database to replace the prior annotations.

**get_data** matches phage data in the selected database to the corresponding data in GenBank (indicated in the Accession field of the *phage* table) and assesses whether the date of the record is more recent than the date of the annotations stored in the database (indicated in the DateLastModified field of the *phage* table). If the GenBank record is more recent, **get_data** stages the new flat file from GenBank and a corresponding import ticket in an import table that are ready to be processed with the 'import' tool.

Additionally, a CSV-formatted summary table of all PhageIDs, their accession, and the results of data retrieval from GenBank is generated.


New non-SEA-PHAGES annotations
------------------------------

Actinobacteriophage genomes that have been sequenced and annotated outside of the SEA-PHAGES program occasionally become available in GenBank. If the genomes should be imported into a MySQL database, the GenBank-formatted flat files need to be manually retrieved from GenBank and staged in a local directory with a manually-generated import ticket.
