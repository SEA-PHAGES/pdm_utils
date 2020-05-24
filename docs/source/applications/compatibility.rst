Compatibility
=============

As the production phage database schema changes and improves, downstream pipelines and tools designed for earlier schema versions may no longer be compatible. These tools may still be used as long as the newest database is first downgraded to the schema version that is compatible with the code.

Conversely, databases created for specific projects or publications using older schema versions may no longer be accessible by newer data analysis pipelines and tools that were designed for newer schema versions. These newer tools may still be used on older databases as long as the older databases are first upgraded to the schema version that is compatible with the analysis code.

To achieve both of these goals, databases can be converted (i.e. either upgraded or downgraded) along an incremental, linear schema history path using the ``pdm_utils convert`` tool. Refer to a description of this :ref:`conversion pipeline <convert>`.
