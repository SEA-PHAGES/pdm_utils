.. _ticketupdate:

Update tickets
==============

Within an update table, an individual row of data populating 5 columns constructs a unique 'ticket'.

    1. **table**: the name of the database table to update.

    2. **field**: the name of the table column to update.

    3. **value**: the new value to be inserted into the column.

    4. **key_name**: the name of the table column by which MySQL can identify which rows will be updated.

    5. **key_value**: the value of the conditional column by which MySQL can identify which rows will be updated.



For example, the update ticket below...


.. csv-table::
    :file: ../../images/update_table_example.csv


...will result in the following MySQL statement...

    UPDATE phage SET Cluster = 'A' WHERE PhageID = 'Trixie';

...and only the Cluster data for Trixie will be updated.
