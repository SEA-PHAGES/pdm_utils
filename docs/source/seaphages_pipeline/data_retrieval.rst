Data retrieval
==============


The first step to update PhameratorDB requires gathering all new genome data (new sequences or annotations), creating corresponding tickets, and systematically staging them in a local directory in preparation for evaluation and import into PhameratorDB. This is accomplished using the retrieve_database_updates.py script, which requires two arguments [Table 3].

.. TODO insert table describing how to use script



This script retrieves four types of data to be imported and creates tickets for them. The database administrator can select to retrieve all, or only specific, types of data. For each type of data, the script stores the retrieved data in a structured directory ready for import: a new folder is created that contains a) one csv-formatted import_table listing each ticket, and b) a genome subdirectory containing flat files.


Metadata updates
----------------

Using the PhagesDB API, genome data for all `sequenced phages <http://phagesdb.org/api/sequenced_phages/>`_ can be directly retrieved. For phage metadata, the script iterates through every phage in PhameratorDB, matches the PhageID to the identical phage name in PhagesDB, and compares metadata stored in PhameratorDB to what is stored in PhagesDB, including: Cluster, Subcluster, Host genus, and Accession. PhagesDB is the primary data source for these fields, so if any information is different between the two databases, a new import ticket is created with the current metadata from PhagesDB, so that PhameratorDB is synchronized accordingly (Table 2, Metadata update).



New auto-annotations
--------------------

When new phage genomes are sequenced, the sequence is uploaded to PhagesDB. The draft auto-annotations are imported into the Actino_Draft instance for immediate reference to end-users (with tools such as the Phamerator GUI and Starterator) until the refined, reliable, final annotations are prepared. PhagesDB tracks which genomes have been imported into Actino_Draft, and it provides a `list <http://phagesdb.org/data/unphameratedlist>`_ of newly sequenced "unphamerated" phages that need to be imported. Automated annotations of these new genomes can be generated through the phage genomics tool, `PECAAN <https://discover.kbrinsgd.org>`_. For each new phage genome, a request for auto-annotation is sent to PECAAN with the URL: https://discoverdev.kbrinsgd.org/phameratoroutput/phage/[PhageID] (where [PhageID] indicates the specific phage name of interest). PECAAN retrieves the new sequence from PhagesDB and automatically annotates coding sequence (CDS) genes using Glimmer [ref Delcher 1999] and GeneMark [ref Borodovsky 2003], and tRNA genes using tRNAscan-SE [ref Lowe and Eddy 1997] and ARAGORN [ref Laslett and Canback 2004]. The script retrieves a GenBank-formatted flat file of the auto-annotations, stores it in the local staging directory, and creates a new import ticket in the import_table [Table 2, New auto-ann.]. When auto-annotated draft genomes are imported into Actino_Draft, PhagesDB removes them from the list of “unphamerated” genomes so that they are not re-processed during subsequent rounds of PhameratorDB updates. It is important to note that since the list of “unphamerated” genomes is created by PhagesDB based on data in the Actino_Draft PhameratorDB instance, this data retrieval step is not reliable if it is used to update alternative PhameratorDB instances.


New preliminary final annotations
---------------------------------

The draft gene annotations are eventually replaced in PhameratorDB with manual, final, gene annotations. The refined annotations are submitted by senior annotators as preliminary final annotations to PhagesDB in GenBank-formatted flat files so that they can be evaluated for quality in the PhameratorDB pipeline. When preliminary final annotations are uploaded to PhagesDB, the flat files are stored on Phagesdb in the qced_genbank_file field with a timestamp stored in the qced_genbank_file_date field. Similar to new metadata retrieval, the retrieval script iterates through every PhageID in PhameratorDB, matches it to the phage name in PhagesDB, reviews the date (if any) that a preliminary final annotation file was uploaded, and if it is more recent than the date of the annotations stored in PhameratorDB (indicated in the DateLastModified field of the Phage table), it retrieves the new flat file from the qced_genbank_file field, stages it in the genome folder, and creates a new import ticket [Table 2, New prelim. final ann.]. The quality of the gene annotations is reviewed during import using the import_phage.py script (see below).

Updated GenBank annotations. After preliminary final annotations are evaluated for quality and imported into PhameratorDB, they are eventually submitted to GenBank with a unique accession number. Annotators provide the accession number to PhagesDB, which gets retrieved for import into PhameratorDB using the metadata retrieval step (see above). In the update GenBank data retrieval step, the script iterates through every PhageID in PhameratorDB, retrieves the accession number from the Accession field (if any), and retrieves the associated genome record from the GenBank nucleotide database using the Entrez python package. For each record, the script identifies the date of the record (Figure 2), and if the record date is more recent than the date of the annotations stored in PhameratorDB (indicated in the DateLastModified field of the Phage table), the script retrieves the updated flat file, stages it in the genome folder, and creates a new import ticket [Table 2, Updated Genbank ann.]. The database administrator is required to provide an email address for this step since the NCBI requests contact information for all users downloading data from GenBank through Entrez. Additionally, the script creates a csv-formatted summary_table of all PhageIDs, their accession, and the results of data retrieval from GenBank.

New non-SEA-PHAGES annotations
------------------------------

Actinobacteriophage genomes that have been sequenced and annotated outside of the SEA-PHAGES program occasionally become available in GenBank. SEA-PHAGES annotators review the quality of these genomes and annotations and assess whether or not they should be imported into the Actino_Draft database. If the genomes should be imported, the database administrator manually retrieves the flat files from GenBank, stages them in a local directory, and creates the appropriate import tickets [Table 2, New non-SEA-PHAGES ann.].

Miscellaneous notes about data retrieval. Currently, the default value populating the Gene description field field is PRODUCT, since this is where descriptions generated from the SEA-PHAGES program are expected to be stored. Also, if GenBank accession numbers are retrieved from PhagesDB, they are retrieved from the genbank_accession field, ensuring that GenBank ACCESSION numbers (and not RefSeq ACCESSION numbers) are stored in PhameratorDB (see above). Last, it is important to note that not all fields are applicable for every ticket type, but the import script requires that all fields are populated. Fields that do not contain relevant data for the ticket type should be populated with none [Table 2].
