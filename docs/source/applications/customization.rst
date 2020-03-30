Create customized phage genomics databases
==========================================
Below is a description of how the ``pdm_utils`` package can be used to
create customized phage genomics databases using the same tools and pipelines
implemented in the SEA-PHAGES program.


First, create a new, empty database with the most current database schema::

    > python3 -m pdm_utils get_db PhageDatabase new

A new MySQL database, PhageDatabase, has been created, and it has the same structure as the most current Actinobacteriophage database.

.. note::

    You may see two MySQL notices appear during installation:
    Warning: 3090, "Changing sql mode 'NO_AUTO_CREATE_USER' is deprecated. It will be removed in a future release."
    Warning: 1305, "PROCEDURE test_empty.split_subcluster does not exist"
    These warnings can be ignored.

If a different schema version is needed, this can be selected at the time of installation::

    > python3 -m pdm_utils get_db PhageDatabase new -s 7

Alternatively, the ``convert`` tool can be used to subsequently downgrade the schema to a prior version.

Next, use ``import`` to add new genomics data to PhageDatabase...
