.. _table_domain:

domain
======

This table stores information about NCBI-defined conserved domains relevant to CDS features within the database.

.. csv-table::
    :file: ../images/database_structure/domain_table.csv

**ID** Auto-incrementing values. This is the primary key.

**HitID** Identifier to match conserved domain data in this table to location of conserved domain in the gene, stored in the *gene_domain* table.

**Description** Description of the conserved domain.

**DomainID** Conserved domain identifier in CDD.

**Name** Conserved domain name in CDD.
