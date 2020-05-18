.. _table_gene:

gene
====

This table contains information that pertains to individual genes, including coordinates, orientation, the gene product, etc.


.. csv-table::
    :file: ../images/database_structure/gene_table.csv


**GeneID** Unique identifier for the gene feature in the entire database. It is distinct from other common gene identifiers in a flat file such as LOCUS_TAG or GENE.

**Name** This field is an identifier for the annotation but does not need to be unique, analogous to the distinction between the PhageID and Name fields in the *phage* table. Most of the time (but not always), it is an integer, but is ocassionally a float, an alphanumeric string, or a strictly alphabetic.  This field is displayed on Phamerator genome maps. This can be derived from several locations in the flat file: the LOCUS_TAG, GENE, PRODUCT, NOTE, and FUNCTION qualifiers. If no gene identifier is present in any of these qualifiers, this field remains empty.

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
