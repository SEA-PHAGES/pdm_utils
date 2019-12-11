SEA-PHAGES data pipeline
========================
Below is a description of how the ``pdm_utils`` package is used to update,
maintain, and manipulate data within the Actinobacteriophage Phamerator
database instances.



.. _figpipeline:

.. figure:: /images/data_pipeline.jpg

    Overview of the pipeline to maintain and update Actino_Draft

.. _figflatfile:

.. figure:: /images/flatfile_parsing.jpg

    Summary of how a GenBank-formatted flat file is parsed


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Overview <./seaphages_pipeline/overview>
   Ticket tracking systems <./seaphages_pipeline/ticket_tracking>
   Retrieving new genomics data <./seaphages_pipeline/data_retrieval>
   Import new genomes <./seaphages_pipeline/data_import>
   Update specific fields <./seaphages_pipeline/data_updates>
   Identifying conserved domains <./seaphages_pipeline/conserved_domains>
   Grouping gene products into phamilies <./seaphages_pipeline/create_phams>
   Export updated databases <./seaphages_pipeline/data_export>
   Freeze databases <./seaphages_pipeline/freeze_data>
   Synchronize databases <./seaphages_pipeline/data_compare>
