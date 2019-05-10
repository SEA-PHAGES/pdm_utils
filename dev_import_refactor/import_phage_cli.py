""" Use this script to import data into Phamerator database using the
UNIX command line interface so that the import process provides
interactive feedback about the import process.
"""



#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib
from datetime import datetime



#Import third-party modules
try:
	from Bio import SeqIO
	from Bio.Alphabet import IUPAC
	from tabulate import tabulate
	import MySQLdb as mdb
except:
	print "\nUnable to import one or more of the following third-party modules: MySQLdb, Biopython, tabulate."
	print "Install modules and try again.\n\n"
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



#Open file to record update information
output_file = open(os.path.join(phageListDir,success_folder,date + "_phage_import_log_" + run_type + "_run.txt"), "w")
write_out(output_file,date + " Phamerator database updates:\n\n\n")
write_out(output_file,"\n\n\n\nBeginning import script...")
write_out(output_file,"\nRun type: " + run_type)






#TODO not sure if I need these counters any more.
# table_errors = 0
# add_total = 0
# remove_total = 0
# replace_total = 0
# update_total = 0
# run_mode_custom_total = 0


# Retrieve import ticket data.
#TODO insert real function
list_of_tickets = prepare_tickets.some_function()



#Tickets will be matched with other genome data.
#Ticket data will be paired with data from PhameratorDB
#and/or a flat file.
#TODO I may want to POP each ticket off this list as I assign to
#matched genome objects.
matched_data_list = []
for ticket in list_of_tickets:
    matched_data_obj = MatchedGenomes()
	matched_data_obj.ticket = ticket
	matched_data_list.append(matched_data_obj)









# Retrieve data from Phamerator.
#TODO insert real function
all_phamerator_data = prepare_phamerator_data.main()






# Now that Phamerator data has been retrieved and
# Phamerator genome objects created, match them to ticket data
for matched_data_obj in matched_data_list:
	phage_id = matched_data_obj.ticket.primary_phage_id

	try:
		matched_genome = phamerator_genome_dict[phage_id]

		#Now add the Phamerator data to the MatchedGenomes object
		matched_data_obj.matched_genomes_dict["phamerator"] = matched_genome

	except:
		#TODO error handling







# Parse flat files and create list of genome objects
#TODO insert real function
all_flat_file_data = prepare_flat_file_data.main()






#TODO need to set strategy variable in advance
# Now that the flat file data has been retrieved and parsed,
# match them to ticket data


flat_file_genome_dict = {}


for flat_file_object in list_of_flat_file_genomes:

	if strategy == "phage_id":
		match_name = flat_file_object.phage_id


	elif strategy == "filename":
		match_name = flat_file_object.filename

	else:
		match_name = ""


	if match_name not in flat_file_genome_dict.keys():

		flat_file_genome_dict[match_name] = flat_file_object

	else:
		pass

		#TODO throw an error - unable to create unique set of objects


for matched_data_obj in matched_data_list:

	match_name = matched_data_obj.ticket.primary_phage_id

	try:
		flat_file_genome = flat_file_genome_dict.pop(match_name)
	except:

		flat_file_genome = None

	#Now add the flat file data to the MatchedGenomes object
	matched_data_obj.matched_genomes_dict["flat_file"] = flat_file_genome







#TODO
# After parsing flat file
# Prepare gene_id and gene_name appropriately















# TODO now that all data is matched, evaluate each ticket


#Evaluate every CDS feature


#Evaluate every tRNA feature



#Evaluate every genome



#Evaluate flat file genome to phagesdb genome




#Evaluate flat file genome to phamerator genome













#Not sure what to do with this:
failed_actions = []
file_tally = 0
script_warnings = 0
script_errors = 0


























#Now that all data has been retrieved, split objects by ticket type
#Create separate lists of ticket based on the indicated action: update, add/replace, remove
#TODO I should pop off each matched_object as I assign to next list.
list_of_update_tickets = []
list_of_remove_tickets = []
list_of_add_replace_tickets = []
list_of_unassigned_tickets = []


for matched_data_obj in matched_data_list:

    ticket_type = matched_data_obj.ticket.type

	if ticket_type == "update":
		list_of_update_tickets.append(matched_data_obj)
	elif ticket_type == "remove":
		list_of_remove_tickets.append(matched_data_obj)
	elif (ticket_type == "add" or ticket_type == "replace"):
		list_of_add_replace_tickets.append(matched_data_obj)

    #TODO error handling
	else:
		write_out(output_file,"\nError: during parsing of actions.")
		table_errors += 1





#TODO Compile all ticket-specific errors and ticket-group errors
# decide how to report errors









###
