.. _toolkit:


Database management pipelines
=============================

Below is a description of command line tools that the ``pdm_utils`` package
contains to analyze and manipulate data within a MySQL database:

Some tools are bound to the most current schema version ("yes"), as they need to know where to find types of data in the database. Other tools ("no"), or sub-tools ("yes"/"no"), are schema-agnostic, and can be used for any type of MySQL database.

Most tools provide functionality without respect to SEA-PHAGES-specific assumptions or goals. Some tools are more oriented to the SEA-PHAGES program by providing (optional) interactions with PhagesDB or by encoding SEA-PHAGES-specific assumptions about the data. For instance, a series of evaluations are implemented in the import and compare pipelines to ensure data quality, and while some of these evaluations are SEA-PHAGES-specific, others are not.

.. toctree::
   :maxdepth: 1

   compare <./tools/compare>
   convert <./tools/convert>
   export <./tools/export>
   find_domains <./tools/find_domains>
   freeze <./tools/freeze>
   get_data <./tools/get_data>
   get_db <./tools/get_db>
   get_gb_records <./tools/get_gb_records>
   import <./tools/import>
   phamerate <./tools/phamerate>
   push <./tools/push>
   pham_review <./tools/pham_review>
   revise <./tools/revise>
   update <./tools/update>



.. csv-table:: pdm_utils tools
  :file: ./images/pipeline_stages.csv


The ``pdm_utils`` toolkit can be used to manage different database instances. However, some tools may only be relevant specifically to the primary instance, Actino_Draft.


.. toctree::
   :hidden:

   supplementary files <./tools/supp_files>
