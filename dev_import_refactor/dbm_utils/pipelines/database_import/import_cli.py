""" Use this script to import data into Phamerator database using the
UNIX command line interface so that the import process provides
interactive feedback about the import process.
"""



# Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib
from datetime import datetime
from classes import Genome



# Import third-party modules
try:
    pass
except ModuleNotFoundError as err:
    print(err)
	sys.exit(1)





#TODO replace with new function
#Get the command line parameters
try:
    database = sys.argv[1]
    phageListDir = sys.argv[2]
    updateFile = sys.argv[3]
except:
    print "\n\n\
            This is a python script to import and update phage genomes in the Phamerator database.\n\
            It requires three arguments:\n\
            First argument: name of MySQL database that will be updated (e.g. 'Actino_Draft').\n\
            Second argument: directory path to the folder of genome files that will be uploaded (genbank-formatted).\n\
            Third argument: directory path to the import table file with the following columns (csv-formatted):\n\
                1. Action to implement on the database (add, remove, replace, update)\n\
                2. PhageID to add or update\n\
                3. Host genus of the updated phage\n\
                4. Cluster of the updated phage\n\
                5. Subcluster of the updated phage\n\
                6. Annotation status of the updated phage (draft, final, gbk)\n\
                7. Annotation authorship of the updated phage (hatfull, gbk)\n\
                8. Gene description field of the updated phage (product, note, function)\n\
                9. Accession of the updated phage\n\
                10. Run mode of the updated phage\n\
                11. PhageID that will be removed or replaced\n\n"
	sys.exit(1)











#TODO replace with new functions
#Expand home directory
home_dir = os.path.expanduser('~')


#Verify the genome folder exists

#Add '/' at the end if it's not there
if phageListDir[-1] != "/":
    phageListDir = phageListDir + "/"


#Expand the path if it references the home directory
if phageListDir[0] == "~":
    phageListDir = home_dir + phageListDir[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
phageListDir = os.path.abspath(phageListDir)



if os.path.isdir(phageListDir) == False:
    print "\n\nInvalid input for genome folder.\n\n"
    sys.exit(1)






#Verify the import table path exists
#Expand the path if it references the home directory
if updateFile[0] == "~":
    updateFile = home_dir + updateFile[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
updateFile = os.path.abspath(updateFile)
if os.path.exists(updateFile) == False:
    print "\n\nInvalid input for import table file.\n\n"
    sys.exit(1)











#TODO replace with new function
#Create output directories
date = time.strftime("%Y%m%d")

failed_folder = '%s_failed_upload_files' % date
success_folder = '%s_successful_upload_files' % date

try:
    os.mkdir(os.path.join(phageListDir,failed_folder))
except:
    print "\nUnable to create output folder: %s" % os.path.join(phageListDir,failed_folder)
    sys.exit(1)


try:
    os.mkdir(os.path.join(phageListDir,success_folder))
except:
    print "\nUnable to create output folder: %s" % os.path.join(phageListDir,success_folder)
    sys.exit(1)













### Below - refactored script in progress



import prepare_tickets
from functions import phamerator
from functions import flat_files
from functions import tickets

# TODO command now should include an argument that specifies phage_id_field
# from which phage_id should be assigned as flat files are parsed.


# TODO confirm arguments are structured properly.

# TODO confirm directories and files exist.

# TODO confirm import table file exists.

# TODO create output directories.

# TODO create SQL object.

list_of_ticket_data # TODO run a function to read the import table and retrieve the data.




# Parse flat files and create list of genome objects
#TODO insert real function

# Identify valid files in folder.
files_in_folder = basic.identify_files(genome_dir)


# Iterate through the list of files.
# Parse each file into a Genome object.
# Returns lists of Genome objects, Eval objects, parsed files,
# and failed files.
flat_file_genomes, valid_files, failed_files = \
    flat_files.create_parsed_flat_file_list(files_in_folder)

# TODO check for flat file parsing errors = exit script if there are errors.
if len(all_results) > 0:
    sys.exit(1)



list_of_matched_objects = import_main.main(list_of_ticket_data, flat_file_genomes, sql_obj)



import_main.import_into_database(list_of_matched_objects, sql_obj)












### Unused code below.


#Not sure what to do with this:
failed_actions = []
file_tally = 0
script_warnings = 0
script_errors = 0































###
