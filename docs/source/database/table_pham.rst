.. _table_pham:

pham
====

This table contains a list of color codes for each unique pham.

.. csv-table::
    :file: ../images/database_structure/pham_table.csv

**PhamID** The primary key of the table. Unique identifier for each pham.

**Color** The hexrgb color code reflecting unique phams, which is used by downstream applications such as Phamerator. The script attempts to maintain consistency of pham designations and colors between rounds of clustering.
