.. _flatfileqc:

Reviewing genome annotations
============================

The Actino_Draft database is routinely updated with new genomics data. When new genome annotations need to to be imported into the database, they are processed using ``pdm_utils`` **import** that reviews the quality of the annotations and how they relate to data already present in the database. The :ref:`import pipeline <import>` processes on GenBank-formatted flat files and checks for a variety of potential errors, including:

    1.	Accidental changes to the genome sequence
    2.	Phage name typos
    3.	Host name typos
    4.	Missing CDS locus tags
    5.	Incorrect protein translation tables used
    6.	Missing protein translations
    7.	Incorrect field used to store the gene descriptions
    8.	Incorrect tRNA genes

After creating the GenBank-formatted flat file, annotators can follow the steps below to review their files using this pipeline to verify that it contains all the necessary information to be successfully imported into the Actino_Draft database:

    1. Ensure that the newest version of ``pdm_utils`` is installed, along with the MySQL and python dependencies (:ref:`installation`).

    2. Ensure you have the most recent version of the Actino_Draft database, using the **get_db** tool (:ref:`getdb`).

    3. Open a Terminal window, create a folder to work in and navigate to it::

        > mkdir validation
        > cd ./validation

    4. Within this new folder, create a csv-formatted import table (such as 'import_table.csv') of import tickets (:ref:`ticketimport`). A template table is provided on the ``pdm_utils`` source code repository on GitHub. Below are tips to structure tickets for routine review your flat files:

        1. Ticket Type should always be “replace”.
        2. Host, Cluster, Subcluster, and Accession: these can be set to “retrieve”.
        3. Annotation Status should always be “final”.
        4. Annotation Author: this should always be “hatfull”.
        5. Gene Description Field should always be “product”.
        6. Run mode should always be “phagesdb”.
        7. Only the Primary PhageID and Secondary PhageID need to be changed for each flat file to be processed.

        Example ticket in ticket table:

        .. csv-table::
            :file: ../images/import_table.csv


    5.	Create a new folder (such as 'genomes') within the validation folder to contain all flat files you would like to check. Ideally with no other types of files should be present::

        > mkdir genomes

    6. Manually move all flat files into that folder.

    7.	Run the **import** pipeline. The pipeline requires you indicate the name of the database, the folder of flat files, and the import table. Below is an example of the command that executes the script, assuming you are still in the ‘validation’ folder::

        > python3 -m pdm_utils import Actino_Draft ./genomes/ ./import_table.csv

    .. note::

        By default, the pipeline runs in 'test' mode so it does not actually make any changes to the database.

    8.	When prompted, provide your MySQL username and password to access your local Actino_Draft database.

    9.	Monitor the output as the file is processed.

    10.	Check log file for warnings and errors.

    11.	Repeat process if needed. After any errors are identified, re-create the flat files with the appropriate corrections, and repeat the import process to ensure the corrected file now passes validation.

    12.	Once everything is correct, upload the flat file to PhagesDB for  official import into the database.
