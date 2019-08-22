Phamerator database structure
=============================
Below is a description of how the Phamerator database is structured.

phage
-----
Below is a description of the ``phage`` table.

.. csv-table:: phage table
    :file: phage_table.csv



**PhageID** This field is the primary key of the *phage* table and is the unique identifier for all phages in the database.  There is a direct correspondence between phage names in PhagesDB or phage names in GenBank records to PhageIDs in the Actino_Draft database, although with some exceptions.

**Name** This field also reflects the phage name, but it is not as constrained as the PhageID, and this name is displayed in the Phamerator GUI. For all 'draft' genomes, the Name contains the PhageID with a '_Draft' suffix appended, indicating the annotations have been automatically annotated. For all other genomes, the Name corresponds to the PhageID.




gene
----
Below is a description of the ``gene`` table.

.. csv-table:: gene table
    :file: gene_table.csv





gene_domain
-----------
Below is a description of the ``gene_domain`` table.

.. csv-table:: gene_domain table
    :file: gene_domain_table.csv






domain
------
Below is a description of the ``domain`` table.

.. csv-table:: domain table
    :file: domain_table.csv





pham
----
Below is a description of the ``pham`` table.

.. csv-table:: pham table
    :file: pham_table.csv






pham_color
----------
Below is a description of the ``pham_color`` table.

.. csv-table:: pham_color table
    :file: pham_color_table.csv






version
-------
Below is a description of the ``version`` table.

.. csv-table:: version table
    :file: version_table.csv
