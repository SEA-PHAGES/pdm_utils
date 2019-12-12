.. _ticketupdate:

Update tickets
==============

Within an update table, an individual row of data populating 5 columns constructs a unique 'ticket'.

    1. **Table name**: the name of the database table to update.

    2. **Target column name**: the name of the table column to update.

    3. **Target value**: the new value to be inserted into the target column.

    4. **Conditional column name**: the name of the table column by which MySQL can identify which rows will be updated.

    5. **Conditional column value**: the value of the conditional column by which MySQL can identify which rows will be updated.

For example, the update ticket below...

    1. phage
    2. Cluster
    3. A
    4. PhageID
    5. Trixie

...will result in the following MySQL statement...

    UPDATE phage SET Cluster = 'A' WHERE PhageID = 'Trixie';

...and only the Cluster data for Trixie will be updated.
