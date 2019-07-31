"""Use this script to import data into Phamerator database using the
UNIX command line interface so that the import process provides
interactive feedback about the import process.
"""



# Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib
from datetime import datetime
from classes import Genome

from functions import phamerator
from functions import flat_files
from functions import tickets
from classes import MySQLConnectionHandler


# Import third-party modules
try:
    pass
except ModuleNotFoundError as err:
    print(err)
	sys.exit(1)


script_description = """
    This is a script to import and update phage genomes in
    the Phamerator database.
    It requires three arguments:
    First argument: name of MySQL database that will be updated
        (e.g. 'Actino_Draft').
    Second argument: directory path to the folder of genome files that
        will be uploaded (genbank-formatted).
    Third argument: directory path to the import table file with the
        following columns (csv-formatted):

        1. Action to implement on the database (add, remove, replace, update)
        2. PhageID to add or update
        3. Host genus of the updated phage
        4. Cluster of the updated phage
        5. Subcluster of the updated phage
        6. Annotation status of the updated phage (draft, final, gbk)
        7. Annotation authorship of the updated phage (hatfull, gbk)
        8. Gene description field of the updated phage (product, note, function)
        9. Accession of the updated phage
        10. Run mode of the updated phage
        11. PhageID that will be removed or replaced
"""



# TODO re-structure using argparse.
# TODO confirm arguments are structured properly.
# Get the command line parameters.
try:
    database = sys.argv[1]
    genome_folder = sys.argv[2]
    import_table = sys.argv[3]
except:
    print(script_description)
	sys.exit(1)


# Confirm that genome folder exists.
genome_folder = basic.expand_path(genome_folder)


# TODO not sure if this is needed anymore.
# #Add '/' at the end if it's not there
# if genome_folder[-1] != "/":
#     genome_folder = genome_folder + "/"


if not basic.verify_path(genome_folder, "dir"):
    print "\n\nInvalid input for genome folder.\n\n"
    sys.exit(1)


# Confirm that import table exists.
import_table = basic.expand_path(import_table)
if not basic.verify_path(import_table, "file"):
    print "\n\nInvalid input for import table file.\n\n"
    sys.exit(1)


# TODO replace with new function?
# Create output directories
date = time.strftime("%Y%m%d")

failed_folder = '%s_failed_upload_files' % date
success_folder = '%s_successful_upload_files' % date

try:
    os.mkdir(os.path.join(genome_folder,failed_folder))
except:
    print "\nUnable to create output folder: %s" % \
        os.path.join(genome_folder,failed_folder)
    sys.exit(1)


try:
    os.mkdir(os.path.join(genome_folder,success_folder))
except:
    print "\nUnable to create output folder: %s" % \
        os.path.join(genome_folder,success_folder)
    sys.exit(1)





# TODO command now should include an argument that specifies id_field
# from which id should be assigned as flat files are parsed.






# TODO create SQL connector object using parsed arguments.
# TODO add sql connection option to script arguments.
if sql_argument:
    sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
    sql_handle.username = "" # Populate from script arguments.
    sql_handle.password = "" # Populate from script arguments.
    sql_handle.database = "" # Populate from script arguments.
else:
    sql_handle = None



# TODO test the sql connection. Verify that the user-provided parameters
# are correct. Otherwise, exit the script.


# Identify valid files in folder for evaluation.
files_in_folder = basic.identify_files(genome_dir)



# TODO match record to ticket in bundle object.


# TODO parsing from import table:
# 1. parse ticket data from table. = prepare_tickets()
# 2. set case for all fields. = prepare_tickets()
# 3. confirm all tickets have a valid type. = check_ticket_structure()
# 4. populate Genome objects as necessary.
# 5. retrieve data if needed.
# 6. check for PhageID conflicts.
# 7. confirm correct fields are populated based on ticket type.


# Retrieve import ticket data.
# lists_of_ticket_data = read.csv(ticket_filename)
lists_of_ticket_data # TODO run a function to read the import table and retrieve the data.





# TODO not sure how many elements (or what types) are returned.
results = import_main.main1(lists_of_ticket_data, files_in_folder, sql_handle)




# TODO after evaluations, if sql argument option is True,
# update the database as needed...




# Now that all flat files and tickets have been evaluated,
# provide summary of results...



###
