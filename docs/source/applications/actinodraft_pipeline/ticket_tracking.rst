.. _tickettracking:

Ticket tracking
===============

With new genome sequences and annotations routinely becoming available for PhameratorDB, and with thousands of researchers routinely using Phamerator, many updates and issues relating to the general management of PhameratorDB, or specifically to the Actino_Draft instance, arise. Together, these diverse actions constitute database “tickets”, such as:

    - New genome sequences need to be auto-annotated and imported
    - New preliminary final annotations need to be reviewed and imported
    - Updated SEA-PHAGES genomes in GenBank need to be imported
    - New non-SEA-PHAGES genomes in GenBank need to be imported
    - Changes need to be made to genome information (e.g. host, cluster, etc.) or to gene information (e.g. coordinates, descriptions, etc.)
    - Genes need to be added or removed

Systematically maintaining and addressing tickets provides a framework for efficient database management. As a result, two separate ticket tracking systems have been developed:

    1. :ref:`ticketemail`
    2. :ref:`ticketimport`



Managing the record of tickets
______________________________

The successful import tickets from the import script and the update tickets from any of the update scripts represent a record of the diverse types of changes made to the database during a round of updates. All tickets and associated files used in each round of updates can be stored in an updates_history folder. Direct changes to the database can be made in MySQL software without using scripts in the ``pdm_utils`` repository, but since this does not generate a record of the changes, it is not recommended as common practice. However, descriptions of these types of changes can also be manually recorded in text files and stored in the updates_history folder as well.
