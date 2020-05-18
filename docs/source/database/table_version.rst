.. _table_version:

version
=======

This table keeps track of the database version and is updated every time the database is changed.


.. csv-table::
    :file: ../images/database_structure/version_table.csv



**Version** This field reflects the current version of the database. Every time changes are made to the database, this integer is incremented by 1.

**SchemaVersion** This field indicates the current version of the database structure (schema) and enhances version control of downstream tools that utilize the database. As the structure of the database changes, such as by the addition or removal of tables or fields, the database schema number can be incremented to reflect that changes have been made. This does not occur often, and needs to be manually changed.
