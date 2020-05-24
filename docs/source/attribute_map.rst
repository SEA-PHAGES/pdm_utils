.. _attributemap:


The bio-centric ORM
===================

The ``pdm_utils`` bio-centric ORM is designed to process data that has been stored/retrieved in a MySQL database, GenBank-formatted flat files, PhagesDB, etc. Data is parsed from GenBank-formatted flat files using :biopython:`BioPython <>`. Data is parsed from PhagesDB using the PhagesDB API. As a result, some ``pdm_utils`` objects contain attributes that directly map to each of these different data structures.

Refer to the brief introductory :ref:`library tutorial <library_tutorial>` to coding with the ``pdm_utils`` library.

Below is a map of how genome-level data is stored within these data structures:

.. csv-table::
    :file: images/genome_attribute_map.csv


Below is a map of how CDS-level data is stored within these data structures:

.. csv-table::
    :file: images/cds_attribute_map.csv



Below is a map of how tRNA-level data is stored within these data structures:

.. csv-table::
    :file: images/trna_attribute_map.csv



Below is a map of how tmRNA-level data is stored within these data structures:

.. csv-table::
    :file: images/tmrna_attribute_map.csv


GenBank-formatted flat files contain a Source feature. Although this data is not stored within the MySQL database, it is parsed and evaluated for quality when the genome is imported into the database. Below is a map of how Source-level data is stored within these data structures:

.. csv-table::
    :file: images/source_attribute_map.csv
