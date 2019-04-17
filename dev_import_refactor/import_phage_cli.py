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








#Retrieve import info from indicated import table file and read all lines into a list and verify contents are correctly populated.
#0 = Type of database action to be performed (add, remove, replace, update)
#1 = New PhageID that will be added to database
#2 = Host of new phage
#3 = Cluster of new phage (singletons should be reported as "singleton")
#4 = Subcluster of new phage (no subcluster should be reported as "none")
#5 = Annotation status of new phage
#6 = Annotation author of the new phage
#7 = Feature field containing gene descriptions of new phage
#8 = Accession
#9 = Run mode
#10 = PhageID of genome to be removed from the database



write_out(output_file,"\n\n\n\nRetrieving import info from table in file...")

file_object = open(updateFile,'r')
file_reader = csv.reader(file_object)
import_table_data_list = []
for input_row in file_reader:
    import_table_data_list.append(input_row)
file_object.close()


#TODO not sure if I need these counters any more.
table_errors = 0
add_total = 0
remove_total = 0
replace_total = 0
update_total = 0
run_mode_custom_total = 0




#Convert data from import file into ticket objects
list_of_tickets = parse_import_tickets(import_table_data_list)



#Tickets will be matched with other genome data.
#Ticket data will be paired with data from PhagesDB, PhameratorDB,
#and/or a flat file.
matched_data_list = []
for ticket in list_of_tickets:
    matched_data_obj = MatchedGenomes()
	matched_data_obj.ticket = ticket
	matched_data_list.append(matched_data_obj)


#Verify all data is cased appropriately.
for matched_data_obj in matched_data_list:
	matched_data_obj.ticket.check_case()



# Retrieve all PhagesDB data for each ticket. Then retrieve all
# necessary data from PhagesDB genome data for each ticket.

for matched_data_obj in matched_data_list:

	ticket = matched_data_obj.ticket
	phagesdb_genome = Genome()

	#TODO make sure api_prefix and api_suffix variables are set
	phage_url = api_prefix + ticket.primary_phage_id + api_suffix


	try:
		online_data_json = urllib.urlopen(phage_url)
		online_data_dict = json.loads(online_data_json.read())


        #Returns a genome object
        phagesdb_genome = parse_phagesdb_data(phagesdb_genome,online_data_dict)


	except:
		online_data_json = ""
		online_data_dict = {}

        #TODO handle error better
		write_out(output_file,"\nError: unable to retrieve Host, Cluster, Subcluster, or Accession data for phage %s from phagesdb." %row[1])






	if ticket.host == "retrieve":
		ticket.host = phagesdb_genome.host
	if ticket.cluster == "retrieve":
		ticket.cluster = phagesdb_genome.cluster
	if ticket.subcluster == "retrieve":
		ticket.subcluster = phagesdb_genome.subcluster
	if ticket.accession == "retrieve":
		ticket.accession = phagesdb_genome.accession



	#Now add the parsed PhagesDB data to the MatchedGenomes object
	matched_data_obj.matched_genomes_dict["phagesdb"] = phagesdb_genome



# Now that data from PhagesDB is matched to the ticket,
# validate each ticket by checking each field in the ticket
# that it is populated correctly.
#TODO not sure if I should pass a list of valid types to this function.

for matched_data_obj in matched_data_list:
	ticket = matched_data_obj.ticket

    validate(ticket)




# Now that individual tickets have been validated,
# validate the entire group of tickets.

temp_list_of_tickets = []
for matched_data_obj in matched_data_list:

    temp_list_of_tickets.append(matched_data_obj.ticket)


#TODO this should return information
validate_tickets(temp_list_of_tickets)





#Create separate lists of ticket based on the indicated action: update, add/replace, remove

list_of_update_tickets = []
list_of_remove_tickets = []
list_of_add_replace_tickets = []


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








#Retrieve database version
#Retrieve current data in database
#0 = PhageID
#1 = Name
#2 = HostStrain
#3 = Sequence
#4 = status
#5 = Cluster2
#6 = DateLastModified
#7 = Accession
#8 = Subcluster2
#9 = AnnotationAuthor
#10 = AnnotationQC
#11 = RetrieveRecord
try:
	con = mdb.connect(mysqlhost, username, password, database)
	con.autocommit(False)
	cur = con.cursor()
except:
	print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
	output_file.close()
	sys.exit(1)

try:

	cur.execute("START TRANSACTION")
	cur.execute("SELECT version FROM version")
	db_version = str(cur.fetchone()[0])
	cur.execute("SELECT PhageID,Name,HostStrain,Sequence,status,\
						Cluster2,DateLastModified,Accession,\
						Subcluster2,AnnotationAuthor,\
						AnnotationQC,RetrieveRecord FROM phage")
	current_genome_data_tuples = cur.fetchall()
	cur.execute("COMMIT")
	cur.close()
	con.autocommit(True)
except:
	mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")


con.close()

write_out(output_file,"\nDatabase: " + database)
write_out(output_file,"\nDatabase version: " + db_version)
write_out(output_file,"\nTotal phages in database before changes: " + str(len(current_genome_data_tuples)))

#Create data sets to compare new data with
#Originally, a phageName_set, phageSequence_set, phageGene_set, and phageGene_dict were implemented, but I ended up not using them here.
#The phageGene_set get implemented later in the script.
#The SQL query still returns these values so if I need to re-implement those, I am able to.
phage_id_set = set()
phage_host_set = set()
phage_status_set = set()
phage_cluster_set = set()
phage_accession_set = set()
phage_subcluster_set = set()

modified_genome_data_lists = []
print "Preparing genome data sets from the database..."


#Parse Phamerator genome data
#Now that sets and dictionaries have been made, create a phamerator_data_dict
#Key = phageID
#Value = Genome object of Phamerator data
phamerator_genome_dict = {}

for genome_tuple in current_genome_data_tuples:

    phamerator_genome.Genome()
    phamerator_genome = parse_phamerator_data(phamerator_genome,genome_tuple)
    phamerator_genome_dict[phamerator_genome.phage_id] = phamerator_genome










#Create sets
for genome in phamerator_genome_dict.keys():


	phage_id_set.add(genome.phage_id)
	phage_host_set.add(genome.host)
	phage_status_set.add(genome.status)
	phage_cluster_set.add(genome.cluster)
	phage_subcluster_set.add(genome.subcluster)


    #TODO I will need to modify this, since accessions are now "none" instead of ""
    #Only add to the accession set if there was an accession, and not if it was empty "".
	phage_accession_set.add(genome.accession)







# Now that Phamerator data has been retrieved and
# Phamerator genome objects created, match them to ticket data





















###
