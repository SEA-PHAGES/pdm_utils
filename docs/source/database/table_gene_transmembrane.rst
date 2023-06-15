.. _table_gene_transmembrane:

gene_transmembrane
==================

This table stores the positions of transmembrane domains predicted by either DeepTMHMM or Sosui within each  CDS feature in the *gene* table.


.. csv-table::
   :file: ../images/database_structure/gene_transmembrane_table.csv


**ID** Auto-incrementing values. This is the primary key.

**GeneID** Unique gene identifier matching GeneID in the *gene* table.

**QueryStart** First amino acid position within the conserved domain.

**QueryEnd** Last amino acid position within the conserved domain.

**Type** Type of transmembrane domain detected.

**Source** Software used to detect transmembrane domains.
