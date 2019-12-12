Reviewing genome annotations
============================

The Actino_Draft database is routinely updated with new genomics data. When new genome annotations need to to be imported into the database, they are processed using the ``pdm_utils`` 'import' pipeline that reviews the quality of the annotations and how they relate to data already present in the database. The genome import pipeline processes on GenBank-formatted flat files and checks for a variety of potential errors, including:

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

    2. Ensure you have the most recent version of the Actino_Draft database, using the 'get_db' tool (:ref:`getdb`).

    3. Open a Terminal window, create a folder to work in and navigate to it::

        > mkdir validation
        > cd ./validation

    4. Within this new folder, create a csv-formatted import table (such as 'import_table.csv') of import tickets (:ref:`ticketimport`). A template table is provided on the ``pdm_utils`` source code repository on GitHub. Below are tips to structure tickets for routine review your flat files:

        1. The Database Action should always be “replace”.
        2. Host, Cluster, Subcluster, and Accession: these can be set to “retrieve”.
        3. Annotation Status should always be “final”.
        4. Annotation Author: this should always be “hatfull”.
        5. Gene Description Field should always be “product”.
        6. Run mode should always be “phagesdb”.
        7. Only the Primary PhageID and Secondary PhageID need to be changed for each flat file to be processed.

        < insert ticket example >

    5.	Create a new folder (such as 'genomes') within the validation folder to contain all flat files you would like to check. Ideally with no other types of files should be present::

        > mkdir genomes

    6. Manually move all flat files into that folder.

    7.	Run the 'import' pipeline. The pipeline requires you indicate the name of the database, the folder of flat files, and the import table. Below is an example of the command that executes the script, assuming you are still in the ‘validation’ folder::

        > python3 -m pdm_utils import Actino_Draft ./genomes/ ./import_table.csv

    .. note::

        By default, the pipeline runs in 'test' mode so it does not actually make any changes to the database.

    8.	When prompted, provide your MySQL username and password to access your local Actino_Draft database.

    9.	Monitor the output as the file is processed:

        As the script iterates through each flat file, it checks numerous fields in the file for accuracy. The script will pause at two types of issues and request input from the user to proceed.

            a.	“Errors”: there are many fields that, if they do not contain the correct information, automatically throws “error” messages, in which the file will not be successfully processed.
            b.	“Warnings”: alternatively, there are many fields that may throw “warnings”, in which there could possibly be a mistake, but the user can indicate this.

        Each time the script highlights warnings or errors, it will require input from the user to proceed. The script may ask “Is this correct?”. It is designed such that if the current values are what the user wants, then simply type “yes”. If the script has identified a bona fide mistake, typing “no” will throw an error and the file will not pass review.

        A table is printed in the terminal window reporting the information found for each gene (locus tag, descriptions found in the product, function, and note fields, translation table, and the first several amino acids of the translation), as well as what the assigned gene names and gene description will be in the database.

    10.	Check output file for warnings and errors. As the script runs, it creates two folders (“failed files” and “successful files”) within the folder that contains the flat files. Each file, after it is processed, is moved into one of those two folders based on whether no errors (“success”) or errors (“fail”) were encountered. In each folder, a csv file is created that contains a list of the failed or successful import actions, as a quick reference. Additionally, in the success folder, a log file is produced that outputs data from the import process, as a record of what was wrong with each file. Searching for “warnings” or “errors” in the file can quickly highlight the potential problems.

    11.	Repeat process if needed. After any errors are identified, re-create the flat files with the appropriate corrections, and repeat the import process to ensure the corrected file now passes validation.

    12.	Once everything is correct, upload the flat file to PhagesDB for  official import into the database.
