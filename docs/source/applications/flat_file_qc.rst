.. _flatfileqc:

Reviewing genome annotations
============================

The Actinobacteriophage database is routinely updated with new genomics data. When new genome annotations need to to be imported into the database, they are processed using ``pdm_utils import`` which reviews the quality of the annotations and how they relate to data already present in the database. The :ref:`import pipeline <import>` processes GenBank-formatted flat files and checks for a variety of potential errors, including:

    1.	Accidental changes to the genome sequence
    2.	Phage name typos
    3.	Host name typos
    4.	Missing CDS locus tags
    5.	Incorrect protein translation tables used
    6.	Missing protein translations
    7.	Incorrect field used to store the gene descriptions
    8.	Incorrect tRNA genes

.. note::

    ``pdm_utils import`` checks the quality of flat files for the purpose of maintaining database consistency and integrity. Although it does check the quality of many aspects of the genome, it is not intended for comprehensive evaluation of the quality, validity, or biological accuracy of all data stored in the flat file.

After creating the GenBank-formatted flat file, annotators can follow the steps below to review their files using this pipeline to verify that it contains all the necessary information to be successfully imported into the Actinobacteriophage database:

    1. Ensure that MySQL is installed. If using a Mac, also ensure that the MySQL server is turned ON (:ref:`installation <installation>`).

    2. Open a Terminal window.

    3. If Conda is used to manage dependencies, activate the Conda environment. If Conda is not used, ensure all dependencies are installed (:ref:`installation <installation>`).

    4. Ensure that the newest version of ``pdm_utils`` is installed (:ref:`installation <installation>`).

    5. Ensure you have the most recent version of the Actinobacteriophage database, using :ref:`getdb <getdb>`.

    6. Create a folder (such as 'validation') to work in and navigate to it::

        > mkdir validation
        > cd ./validation

    7. Within this new folder, create a csv-formatted import table (such as 'import_table.csv') of :ref:`import tickets <ticketimport>`. A template table is provided on the ``pdm_utils`` source code repository on GitHub. Below are tips to structure tickets for routine review your flat files:

        1. Ticket Type should be set to “replace”.
        2. Host, Cluster, Subcluster, and Accession should be set to “retrieve”.
        3. Annotation Status should be set to “final”.
        4. Annotation Author this should be set to “hatfull”.
        5. Gene Description Field should be set to “product”.
        6. Run mode should be set to “phagesdb”.
        7. Only the Primary PhageID and Secondary PhageID need to be changed for each flat file.

        Example ticket in ticket table (columns labeled only for illustration):

        .. csv-table::
            :file: ../images/import_table.csv


    8.	Create a new folder (such as 'genomes') within the validation folder to contain all flat files you would like to check::

        > mkdir genomes

    9. Manually move all flat files into that folder. No other files should be present.

    10.	Run ``import``. The pipeline requires you to indicate the name of the database, the folder of flat files, and the import table. Below is an example of the command that executes the script, assuming you are still in the ‘validation’ folder::

        > python3 -m pdm_utils import Actinobacteriophage ./genomes/ ./import_table.csv

    .. note::

        By default, the pipeline runs in 'test' mode so it does not actually make any changes to the database.

    11.	When prompted, provide your MySQL username and password to access your local Actinobacteriophage database.

    12.	Monitor the output as the file is processed.

    13.	After the evaluation is complete, review specific warnings and errors in the log file if needed.

    14.	Repeat process if needed. After any errors are identified, re-create the flat files with the appropriate corrections, and repeat the import process to ensure the corrected file now passes validation.

    15.	Once everything is correct, upload the flat file to PhagesDB for official import into the database.