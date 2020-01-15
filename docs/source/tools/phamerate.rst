.. _phamerate:

phamerate: create gene phamilies
================================

All gene products stored in the database can be grouped into phamilies (“phams”) using the ``pdm_utils phamerate`` tool.

To phamerate gene products::

    > python3 -m pdm_utils phamerate Actino_Draft


The argument 'Actino_Draft' indicates the name of the database in which gene products will be grouped into phamilies.


In general, the script groups genes into phamilies using a kmer-based strategy implemented with MMSeqs [ref ].

< Insert description of phamerate pipeline >

The new phams and their unique colors are inserted into the *pham* table of the database. The script attempts to maintain consistency of pham designations and colors between rounds of clustering, although this is not strictly enforced.
