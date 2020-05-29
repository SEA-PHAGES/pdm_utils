.. _update:

update: make updates to specific database fields
================================================

Sometimes it is necessary to modify or update specific fields for specific phages. This can be accomplished using the ``pdm_utils update`` tool that relies on a specifically structured update ticket (:ref:`ticketupdate`)::

    > python3 -m pdm_utils update Actinobacteriophage -f ./update_table.csv -c config_file.txt

The argument 'Actinobacteriophage' indicates the name of the database in which the updates should be implemented. The './update_table.csv' indicates the CSV-formatted table containing the list of update tickets. Each ticket is implemented with very little QC. As with other pipelines, use of the :ref:`config_file` option can automate accessing MySQL.

Additionally, ``update`` can be used to quickly increment the database version, which needs to be changed every time changes are made in the database::

    > python3 -m pdm_utils update Actinobacteriophage -v
