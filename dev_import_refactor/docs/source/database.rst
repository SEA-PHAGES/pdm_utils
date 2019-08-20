Phamerator database structure
=============================
Below is a description of how the Phamerator database is structured.

**phage** table
---------------
Below is a description of the ``phage`` table.

.. csv-table:: phage table
    :file: phage_table.csv



**PhageID** This field is the primary key of the *phage* table and is the unique identifier for all phages in the database.  There is a direct correspondence between phage names in PhagesDB or phage names in GenBank records to PhageIDs in the Actino_Draft database, although with some exceptions.

**Name** This field also reflects the phage name, but it is not as constrained as the PhageID, and this name is displayed in the Phamerator GUI. For all 'draft' genomes, the Name contains the PhageID with a '_Draft' suffix appended, indicating the annotations have been automatically annotated. For all other genomes, the Name corresponds to the PhageID.
