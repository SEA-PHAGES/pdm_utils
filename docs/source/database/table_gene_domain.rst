.. _table_gene_domain:

gene_domain
===========

This table stores the positions of NCBI-defined conserved domains within each CDS feature in the *gene* table.


.. csv-table::
    :file: ../images/database_structure/gene_domain_table.csv


**ID** Auto-incrementing values. This is the primary key.

**GeneID** Unique gene identifier matching GeneID in the *gene* table.

**HitID** Identifier to match location of conserved domain in this table to conserved domain data, stored in the *domain* table.

**QueryStart** First amino acid position within the conserved domain.

**QueryEnd** Last amino acid position within the conserved domain.

**Expect** E-value reflecting significance of the domain hit.
