.. _table_phage:

phage
=====

This table contains information that pertains to the entire phage genome, such as the genome sequence, the host strain, the designated cluster, etc.

.. csv-table::
    :file: ../images/database_structure/phage_table.csv



**PhageID** This field is the primary key of the *phage* table and is the unique identifier for all phages in the database.  There is a direct correspondence between phage names in PhagesDB or phage names in GenBank records to PhageIDs in the Actino_Draft database (although there are a few exceptions, due to naming restrictions in different external databases).

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
