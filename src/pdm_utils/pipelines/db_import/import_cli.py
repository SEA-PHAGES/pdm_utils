"""Use this script to import data into Phamerator database using the
UNIX command line interface so that the import process provides
interactive feedback about the import process.
"""
# Built-in libraries
import time
import sys
import os
import getpass
import csv
import re
import shutil
import json
import urllib
import argparse
from datetime import datetime
# Import third-party modules
try:
    pass
except ModuleNotFoundError as err:
    print(err)
    sys.exit(1)

# Package modules
from pdm_utils.functions import phamerator
from pdm_utils.functions import flat_files
from pdm_utils.functions import tickets
from pdm_utils.functions import basic
from pdm_utils.classes import mysqlconnectionhandler
from pdm_utils.classes import genome
from pdm_utils.constants import constants

#
# script_description = """
#     Structure of import ticket table:
#         following columns (csv-formatted):
#         1. Action to implement on the database (add, remove, replace, update)
#         2. PhageID to add or update
#         3. Host genus of the updated phage
#         4. Cluster of the updated phage
#         5. Subcluster of the updated phage
#         6. Annotation status of the updated phage (draft, final, gbk)
#         7. Annotation authorship of the updated phage (hatfull, gbk)
#         8. Gene description field of the updated phage (product, note, function)
#         9. Accession of the updated phage
#         10. Run mode of the updated phage
#         11. PhageID that will be removed or replaced
#     """

parser = argparse.ArgumentParser(
    description="Import genomes into a Phamerator database.")
parser.add_argument("-g", "--genome_folder", type=os.path.abspath,
    required=True,
    help=("Path to the folder containing GenBank-formatted "
          "flat files to be processed."))
parser.add_argument("-t", "--table", type=os.path.abspath, required=True,
    help="Path to the CSV-formatted 12-field table containing instructions "
         "to process each genome. ")
parser.add_argument("-db", "--database", type=str, default=False,
    help="Name of the MySQL database to import the genomes.")
parser.add_argument("-f", "--filename", action="store_true", default=False,
    help="Indicates whether the filename should be used to identify the genome.")
args = parser.parse_args()

# Confirm that genome folder and import table exists.
if not basic.verify_path(args.genome_folder, "dir"):
    print("\n\nInvalid input for genome folder.\n\n")
    sys.exit(1)
if not basic.verify_path(args.table, "file"):
    print("\n\nInvalid input for import table file.\n\n")
    sys.exit(1)

# Create output directories
date = time.strftime("%Y%m%d")
attempts = 3
failed_folder = '%s_failed_upload_files' % date
success_folder = '%s_successful_upload_files' % date
f_valid = False
f_count = 0
while not f_valid and f_count < attempts:
    if f_count > 0:
        failed_folder_mod += "_" + str(f_count)
    failed_dir = os.path.join(args.genome_folder,failed_folder)
    if not basic.verify_path(failed_dir, "dir"):
        f_valid = True
        os.mkdir(failed_dir)
    f_count += 1
if not f_valid:
    print("\nUnable to create failed_folder")
    sys.exit(1)
success_folder = '%s_failed_upload_files' % date
success_folder = '%s_successful_upload_files' % date
success_folder_valid = False
success_folder_count = 0
while not success_folder_valid and success_folder_count < attempts:
    if success_folder_count > 0:
        success_folder += "_" + str(success_folder_count)
    failed_dir = os.path.join(args.genome_folder,success_folder)
    if not basic.verify_path(failed_dir, "dir"):
        success_folder_valid = True
        os.mkdir(failed_dir)
    success_folder_count += 1
if not success_folder_valid:
    print("\nUnable to create success_folder")
    sys.exit(1)


# TODO add sql connection option to script arguments?
if args.database is not False:
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.username = getpass.getpass(prompt="Provide the MySQL username: ")
    sql_handle.password = getpass.getpass(prompt="Provide the MySQL password: ")
    sql_handle.database = args.database

    # TODO test the sql connection. Verify that the user-provided parameters
    # are correct. Otherwise, exit the script.
    sql_handle.validate_credentials()
    sql_handle.validate_database_access()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        print("\nUnable to connect to the database")
        sys.exit(1)
else:
    sql_handle = None





# Identify valid files in folder for evaluation.
files_in_folder = basic.identify_files(args.genome_folder)


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
lists_of_ticket_data = []
with open(args.table,'r') as file:
    file_reader = csv.reader(file)
    for row in file_reader:
        lists_of_ticket_data.append(row)




# TODO not sure how many elements (or what types) are returned.
results = import_main.main1(lists_of_ticket_data, files_in_folder, sql_handle)




# TODO after evaluations, if sql argument option is True,
# update the database as needed...




# Now that all flat files and tickets have been evaluated,
# provide summary of results...



###
