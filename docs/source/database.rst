.. _dbstructure:

Phamerator database structure
=============================

The current database schema (schema 5) contains 10 tables:

    1.  phage
    2.  gene
    3.  domain
    4.  gene_domain
    5.  pham
    6.  pham_color
    7.  trna
    8.  tmrna
    9.  scores_summary
    10. version

.. TODO probably could insert image showing how tables are connected.
.. TODO should this be generated from MySQL Workbench?
.. .. csv-table::
..     :file: images/database_structure/database.csv
..     :widths: 10, 10


phage
-----
This table contains information that pertains to the entire phage genome, such as the genome sequence, the host strain, the designated cluster, etc.

.. csv-table::
    :file: images/database_structure/phage_table.csv



**PhageID** This field is the primary key of the *phage* table and is the unique identifier for all phages in the database.  There is a direct correspondence between phage names in PhagesDB or phage names in GenBank records to PhageIDs in the Actino_Draft database (although there are a few exceptions, due to naming restrictions in different databases).

**Name** This field also reflects the phage name, but it is not as constrained as the PhageID, and this name is displayed in the Phamerator GUI. For all 'draft' genomes, the Name contains the PhageID with a '_Draft' suffix appended, indicating the annotations have been automatically annotated. For all other genomes, the Name corresponds to the PhageID.

**Accession** This field is populated and updated directly from import tickets and is used for auto-updating genomes from GenBank records. It is important to note that the NCBI generates RefSeq records that are derived from GenBank records. After data is submitted to GenBank, authors retain control of the GenBank record but not the RefSeq record. As a result, the PhameratorDB Accession field should always store the GenBank ACCESSION number (and not the RefSeq ACCESSION number) for SEA-PHAGES genomes. For non-SEA-PHAGES genomes, either accession number may be stored. In either case, the Accession should not contain the sequence version (represented by the integer to the right of the decimal).

**HostStrain** This field indicates the host genus (e.g. *Mycobacterium*, *Streptomyces*, etc.) from which the phage was isolated.

**Sequence** This genome nucleotide sequence of the phage.

**SequenceLength** The length of the phage’s genome sequence, computed by the import script.

**GC** The GC% of the genome sequence, computed by the import script.

**Cluster2** This field indicates the phage’s cluster designation if it has been clustered. If the phage is a singleton, it remains empty (NULL).

**Subcluster2** This field indicates the phage’s subcluster designation if it has been subclustered, otherwise it remains empty (NULL).

**Cluster** This field combines information from Cluster2 and Subcluster2, and is used by certain applications (such as Phamerator GUI). It remains empty (NULL) if the phage is a singleton (and is not clustered), is populated with the cluster designation if the phage is clustered (but not subclustered), or is populated with the subcluster designation if the phage is clustered and subclustered.

**DateLastModified** This field records the date in which a genome and its annotations have been imported. This is valuable to keep track of which annotation version has been imported, and it also facilitates automated updating of the database. It is important to note that the date stored in this field reflects the date the annotation data were imported, and not the date that the annotation data were created. Although the field is a DATETIME data type, only date data is stored, and no time data is retained.

**AnnotationAuthor** This field indicates if the genome sequence and annotations are (1) or are not (0) maintained by the SEA-PHAGES program, and it facilitates automatic updates from GenBank. If a genome has been sequenced and annotated through the SEA-PHAGES program, its GenBank record is actively updated/maintained.  In this case, “Graham Hatfull” is expected to be a listed author in the GenBank record. (All genomes through the SEA-PHAGES program should have "Graham Hatfull" as a listed author, but not all GenBank records listing "Graham Hatfull" as an author are derived from the SEA-PHAGES program.)

**RetrieveRecord** This field facilitates automatic updates from GenBank records. Most SEA-PHAGES genomes are expected to be automatically updated from GenBank once they are assigned a unique GenBank accession. However, some genomes, such as those generated from non-SEA-PHAGES researchers, may not need to be automatically updated. This field is set to 1 for genomes that are to be automatically updated and set to 0 for those genomes that are not to be automatically updated. Initially, this field is indirectly determined by the AnnotationAuthor field. For newly added genomes, if AnnotationAuthor = hatfull in the import ticket, this field is set to 1, otherwise it is set to 0. For genomes being replaced (by automatic updates from GenBank or by the creation of manual tickets), the value in this field is retained. Note: this field will be either 0 or 1.

**Status** This field indicates whether the gene annotations have automatically (draft) or manually (final) annotated, or whether the annotation strategy is unknown (gbk).


gene
----
This table contains information that pertains to individual genes, including coordinates, orientation, the gene product, etc.

.. csv-table::
    :file: images/database_structure/gene_table.csv


**GeneID** Unique identifier for the gene annotation in the database. This can be derived three ways. First, it can simply be synonymous with the LOCUS_TAG of the CDS feature in the flat file. For SEA-PHAGES flat files, this is usually the case. However, for non-SEA-PHAGES flat files, there may not be a LOCUS_TAG for every, or any, CDS feature. As a result, the GeneID can be computed by concatenating the PhageID with the CDS count (which indicates the order that the CDS was parsed from the feature list during import). However, neither of these naming strategies guarantee a unique identifier, and naming conflicts may arise with features already present in the *gene* table. In this case, a _duplicateID[0123] suffix is appended to the GeneID (where [0123] is an integer).

**Name** This field is an identifier for the annotation but does not need to be unique. Most of the time (but not always), it is a number. This field is displayed on Phamerator GUI genome maps. [Add how this is computed in the script]

**PhageID** The name of the phage genome from which the gene is derived, matching one of the phage names in the PhageID of the *phage* table.

**Start, Stop** These fields store the genomic coordinates marking the coordinate boundaries of the gene. Start and Stop reflect the left and right (respectively) boundaries of the gene based on the genome orientation stored in the database. Note: the coordinates are stored in 0-based half-open format (as opposed to the 1-based closed format used in GenBank records). For practical purposes, the start coordinate has been decreased by 1 nucleotide.

**Orientation** This field indicates the strand in which the feature is encoded.

**Length** This field indicates the nucleotide length of the gene, computed by the length of the amino acid sequence. Note: this field needs to be improved to maintain data integrity.

**Translation** This field contains the translated amino acid sequence and is derived directly from the GenBank record. Note: currently, the maximum length of the translation product is 5,000 amino acids.

**LocusTag** This field facilitates automatic updating of GenBank records. Once a genome has been submitted to GenBank, genes are assigned unique locus tags in the LOCUS_TAG field. These identifiers cannot be changed, and annotators are required to use them when requesting to update details about individual genes. This field provides a direct link to the corresponding GenBank feature. Note: this field is only populated for records retrieved from GenBank.

**Notes** This field contains data on the gene function, and is derived from one of several fields of the GenBank feature. [Add more info on how it is parsed here?]

**DomainStatus** Indicates whether conserved domain data has been retrieved for this feature. When new phage genomes are added to PhameratorDB, the DomainStatus field for each new gene is set to 0. The cdd_script.py script retrieves gene products (stored in the Translation field of the *gene* table) for all genes with DomainStatus < 1. The rpsblast+ package is used to identity conserved domains using BLAST with an e-value threshold = 0.001. For each gene, retrieved CDD data is inserted into the *domain* and *gene_domain* tables, and the DomainStatus field in the *gene* table is set to 1 so that this gene is not re-processed during subsequent rounds of updates. Note: this field will be either 0 or 1.




gene_domain
-----------
This table stores the positions of NCBI-defined conserved domains within each CDS feature in the *gene* table.


.. csv-table::
    :file: images/database_structure/gene_domain_table.csv


**ID** Auto-incrementing values. This is the primary key.

**GeneID** Unique gene identifier matching GeneID in the *gene* table.

**HitID** Identifier to match location of conserved domain in this table to conserved domain data, stored in the *domain* table.

**QueryStart** First amino acid position within the conserved domain.

**QueryEnd** Last amino acid position within the conserved domain.

**Expect** E-value reflecting significance of the domain hit.





domain
------
This table stores information about NCBI-defined conserved domains relevant to CDS features within the database.

.. csv-table::
    :file: images/database_structure/domain_table.csv

**ID** Auto-incrementing values. This is the primary key.

**HitID** Identifier to match conserved domain data in this table to location of conserved domain in the gene, stored in the *gene_domain* table.

**Description** Description of the conserved domain.

**DomainID** Conserved domain identifier in CDD.

**Name** Conserved domain name in CDD.




pham
----
This table contains a list of CDS features from the *gene* table with their computed pham.

.. csv-table::
    :file: images/database_structure/pham_table.csv


**GeneID** Corresponds to unique GeneIDs from *gene* table.

**Name** Unique pham numbers.

**OrderAdded** Auto-incrementing values.




pham_color
----------
This table contains a list of color codes for each unique pham.

.. csv-table::
    :file: images/database_structure/pham_color_table.csv


**ID** The primary key of the table. Auto-incrementing values.

**Name** Unique identifier for each hexrgb color code.

**Color** The hexrgb color code reflecting unique phams, which is used to create phamerator maps. The script attempts to maintain consistency of pham designations and colors between rounds of clustering.




version
-------
This table keeps track of the database version and is updated every time the database is changed.

.. csv-table::
    :file: images/database_structure/version_table.csv



**Version** This field reflects the current version of the database. Every time changes are made to the database, this integer is incremented by 1.

**SchemaVersion** This field indicates the current version of the database structure, or schema and enhances version control of scripts that directly communicate with PhameratorDB. As the structure of the database changes, such as by the addition or removal of tables or fields, the database schema number can be incremented to reflect that changes have been made. This does not occur often, and needs to be manually changed.
