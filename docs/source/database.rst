.. _dbstructure:

MySQL database structure
========================

The current schema (schema version 8) of the MySQL database contains the following tables:

    1.  phage
    2.  gene
    3.  domain
    4.  gene_domain
    5.  pham
    6.  version
    7.  trna (in development)
    8.  trna_structures (in development)
    9.  tmrna (in development)


.. _figschema:

.. figure:: /images/database_structure/schema_8_map.jpg

    Map of the MySQL database schema (schema version 8).

.. .. csv-table::
    :file: images/database_structure/database.csv


.. :widths: 10, 10


phage
-----
This table contains information that pertains to the entire phage genome, such as the genome sequence, the host strain, the designated cluster, etc.

.. csv-table::
    :file: images/database_structure/phage_table.csv



**PhageID** This field is the primary key of the *phage* table and is the unique identifier for all phages in the database.  There is a direct correspondence between phage names in PhagesDB or phage names in GenBank records to PhageIDs in the Actinobacteriophage database (although there are a few exceptions, due to naming restrictions in different external databases).

**Name** This field also reflects the phage name, but it is not as constrained as the PhageID, and does not have to be unique. For all 'draft' genomes, the Name contains the PhageID with a '_Draft' suffix appended, indicating the annotations have been automatically annotated. For all other genomes, the Name corresponds to the PhageID. In some downstream applications, such as Phamerator, this serves as the phage's display name.

**Accession** This field is populated and updated directly from import tickets and is used for auto-updating genomes from GenBank records. It is important to note that the NCBI generates RefSeq records that are derived from GenBank records. After data is submitted to GenBank, authors retain control of the GenBank record but not the RefSeq record. As a result, this field should always store the GenBank ACCESSION number (and not the RefSeq ACCESSION number) for SEA-PHAGES genomes. For non-SEA-PHAGES genomes, either accession number may be stored. In either case, the Accession should not contain the sequence version (represented by the integer to the right of the decimal).

**HostGenus** This field indicates the host genus (e.g. *Mycobacterium*, *Streptomyces*, etc.) from which the phage was isolated.

**Sequence** This genome nucleotide sequence of the phage.

**Length** The length of the phage’s genome sequence.

**GC** The GC% of the genome sequence.

**Cluster** This field indicates the phage’s cluster designation if it has been clustered. If the phage is a singleton, it remains empty (NULL).

**Subcluster** This field indicates the phage’s subcluster designation if it has been subclustered, otherwise it remains empty (NULL).

**DateLastModified** This field records the date in which a genome and its annotations have been imported. This keeps track of which annotation version has been imported, and it facilitates automated updating of the database. It is important to note that the date stored in this field reflects the date the annotation data were imported, and not the date that the annotation data were created. Although the field is a DATETIME data type, only date data is stored, and no time data is retained.

**AnnotationAuthor** This field indicates if the genome sequence and annotations are (1) or are not (0) maintained by the SEA-PHAGES program, and it facilitates automatic updates from GenBank. If a genome has been sequenced and annotated through the SEA-PHAGES program, its GenBank record is actively updated/maintained.  In this case, “Graham Hatfull” is expected to be a listed author in the GenBank record. (All genomes through the SEA-PHAGES program should have "Graham Hatfull" as a listed author, but not all GenBank records listing "Graham Hatfull" as an author are derived from the SEA-PHAGES program.)

**RetrieveRecord** This field will be 0 or 1, and it facilitates automatic updates from GenBank records . Most SEA-PHAGES genomes are expected to be automatically updated from GenBank once they are assigned a unique GenBank accession. However, some genomes, such as those generated from non-SEA-PHAGES researchers, may not need to be automatically updated. This field is set to 1 for genomes that are to be automatically updated and set to 0 for those genomes that are not to be automatically updated. Initially, this field is indirectly determined by the AnnotationAuthor field. For newly added genomes, if AnnotationAuthor = 1 in the import ticket, this field is set to 1, otherwise it is set to 0. For genomes being replaced (by automatic updates from GenBank or by the creation of manual tickets), the value in this field is retained.

**Status** This field indicates whether the gene annotations have automatically (draft) or manually (final) annotated, or whether the annotation strategy is unknown (unknown).


gene
----
This table contains information that pertains to individual genes, including coordinates, orientation, the gene product, etc.


.. csv-table::
    :file: images/database_structure/gene_table.csv


**GeneID** Unique identifier for the gene feature in the entire database. It is distinct from other common gene identifiers in a flat file such as LOCUS_TAG or GENE.

**Name** This field is an identifier for the annotation but does not need to be unique, analogous to the distinction between the PhageID and Name fields in the *phage* table. Most of the time (but not always), it is a number. This field is displayed on Phamerator genome maps. This can be derived two ways. If the CDS feature in the flat file contains a GENE qualifier, this is set to the Name. Not all features contain this qualifier though. If this is the case, if the CDS feature contains an integer in the LOCUS_TAG qualifier, this is set to the Name. Not all features contain an integer in the qualifier, or even a qualifier at all, though. If neither of these conditions are met, this field remains empty.

**PhageID** The name of the phage genome from which the gene is derived, matching one of the phage names in the PhageID of the *phage* table.

**Start, Stop** These fields store the genomic coordinates marking the coordinate boundaries of the gene. The coordinates are stored in '0-based half-open' format (as opposed to the '1-based closed' format used in other representations, such as a GenBank-formatted flat file). For practical purposes, the start coordinate has been decreased by 1 nucleotide. Start and Stop reflect the left and right (respectively) boundaries of the feature based on the genome orientation stored in the database. They do not directly reflect the translational start and stop coordinates of the feature, which are dependent on orientation. Since only two coordinates are stored for each feature, compound features spanning more than one contiguous region of the genome (such as features that wrap-around genome termini or features with a translational frameshift) are not completely represented in the database.

**Orientation** This field indicates the strand in which the feature is encoded.

**Parts** This field indicates the number of regions in the genome that define the feature. Only two coordinates are stored for each feature, which is an accurate representation of the majority of features. However, the definition of some features, such as those that extend across the genome termini or those that contain a frameshift, are not completely represented with this strategy. This field is used to discriminate between these types of features.

**Length** This field indicates the length of the transcribed and translated nucleotide sequence and is directly proportional to Translation. If the feature is comprised of a single ORF, this is equal to Stop - Start. If the feature is comprised of two ORFs (such as due to ribosome frameshifting), Length will be larger than Stop - Start.

**Translation** This field contains the translated amino acid sequence and is derived directly from the GenBank record. Note: currently, the maximum length of the translation product is 5,000 amino acids.

**LocusTag** This field facilitates automatic updating of GenBank records. Once a genome has been submitted to GenBank, genes are assigned unique locus tags in the LOCUS_TAG field. These identifiers cannot be changed, and annotators are required to use them when requesting to update details about individual genes. This field provides a direct link to the corresponding GenBank feature. Note: this field is only populated for records retrieved from GenBank.

**Notes** This field contains data on the gene function, and is derived from one of several fields of the GenBank feature.

**DomainStatus** Indicates whether the ``find_domains`` pipeline has searched for conserved domain data for this feature (0 or 1). When new phage genomes are added to the *gene* table, the DomainStatus field for each new gene is set to 0. The ``find_domains`` pipeline updates this to 1 after searching for conserved domains (regardless of whether the feature contains any conserved domains).

**PhamID** Pham designation for the translation, matching one of the PhamIDs in the *pham* table.



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
This table contains a list of color codes for each unique pham.

.. csv-table::
    :file: images/database_structure/pham_table.csv

**PhamID** The primary key of the table. Unique identifier for each pham.

**Color** The hexrgb color code reflecting unique phams, which is used by downstream applications such as Phamerator. The script attempts to maintain consistency of pham designations and colors between rounds of clustering.




version
-------
This table keeps track of the database version and is updated every time the database is changed.


.. csv-table::
    :file: images/database_structure/version_table.csv



**Version** This field reflects the current version of the database. Every time changes are made to the database, this integer is incremented by 1.

**SchemaVersion** This field indicates the current version of the database structure (schema) and enhances version control of downstream tools that utilize the database. As the structure of the database changes, such as by the addition or removal of tables or fields, the database schema number can be incremented to reflect that changes have been made. This does not occur often, and needs to be manually changed.



trna (in development)
---------------------
This table contains information that pertains to individual tRNA features.

trna_structures (in development)
--------------------------------
This table contains information that pertains to tRNA secondary structure.

tmrna (in development)
----------------------
This table contains information that pertains to individual tmRNA features.
