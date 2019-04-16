#!/usr/bin/env python
#Phage Import Script
#Based on original script by Charles Bowman and Damani Brown on 20150216
#University of Pittsburgh
#Revamped by Travis Mavrich on 20160701
#Now it is designed to be a single script to upload genomes, update host/cluster/status data, remove genomes, and check integrity of files.
#WORKS FOR .gb, .gbf, .gbk, or .txt files




#Flow of the import process
#1 Import python modules, set up variables and functions
#2 Retrieve current database information
#3 Retrieve import information from file
#4 Update genome information for genomes already in the database
#5 Remove genomes from database that have no replacements
#6 Add or replace genomes

#The strategy for handling MySQL errors:
#The database is accessed in Steps 2, 4, 5, and 6. Connections are opened and closed at each Step. This ensures that any unforeseen errors encountered elsewhere in the script
#that force the script to exit do not cause ROLLBACK errors or connection errors.



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




















































#########################


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
phageId_set = set()
phageHost_set = set()
phageStatus_set = set()
phageCluster_set = set()
phageAccession_set = set()
phageSubcluster_set = set()

modified_genome_data_lists = []
print "Preparing genome data sets from the database..."

#Modify several fields of the data and create sets
for genome_tuple in current_genome_data_tuples:
	phageId_set.add(genome_tuple[0])
	phageHost_set.add(genome_tuple[2])
	phageStatus_set.add(genome_tuple[4])
	phageCluster_set.add(genome_tuple[5])
	phageSubcluster_set.add(genome_tuple[8])

	#If there is no date in the DateLastModified field, set it to a very early date
	if genome_tuple[6] is None:
		modified_datelastmod = datetime.strptime('1/1/1900','%m/%d/%Y')
	else:
		modified_datelastmod = genome_tuple[6]


	#The Accession field in Phamerator defaults to "". Some accessions have the version number suffix.
	#Process Accession data. Discard the version number.
	if genome_tuple[7] != "":
		modified_accession = genome_tuple[7].split('.')[0]
		phageAccession_set.add(modified_accession) #Only add to the accession set if there was an accession, and not if it was empty "".
	else:
		modified_accession = "none"

	#Add all modified data into new list
	#0 = PhageID
	#1 = Name
	#2 = HostStrain
	#3 = Sequence
	#4 = status
	#5 = Cluster2
	#6 = Modified DateLastModified
	#7 = Modified Accession
	#8 = Subcluster2
	#9 = AnnotationAuthor
	#10 = AnnotationQC
	#11 = RetrieveRecord
	modified_genome_data_lists.append([genome_tuple[0],\
										genome_tuple[1],\
										genome_tuple[2],\
										genome_tuple[3],\
										genome_tuple[4],\
										genome_tuple[5],\
										modified_datelastmod,\
										modified_accession,\
										genome_tuple[8],\
										str(genome_tuple[9]),\
										str(genome_tuple[10]),\
										str(genome_tuple[11])])



#Now that sets and dictionaries have been made, create a phamerator_data_dict
#Key = phageID
#Value = list of modified data retrieved from MySQL query
phamerator_data_dict = {}
for element in modified_genome_data_lists:
	phamerator_data_dict[element[0]] = element
























#Create output file to store successful actions implemented
success_action_file = '%s_successful_import_table_actions.csv' % date
success_action_file_handle = open(os.path.join(phageListDir,success_folder,success_action_file),"w")
success_action_file_writer = csv.writer(success_action_file_handle)



#Update actions implemented
#Prepare to update fields that are not accompanied by adding or replacing a genome.
write_out(output_file,"\n\n\n\nUpdating Host, Cluster, Subcluster, and Status fields...")
updated = 0
update_statements = []
for genome_data in update_data_list:
	write_out(output_file,"\nPreparing: " + str(genome_data))

	if genome_data[7] == "none":
		genome_data[7] = ""

	#HostStrain, status, Accession, Author updates.
	update_statements.append("UPDATE phage SET HostStrain = '" + genome_data[2] + "' WHERE PhageID = '" + genome_data[1] + "';")
	update_statements.append("UPDATE phage SET status = '" + genome_data[4] + "' WHERE PhageID = '" + genome_data[1] + "';")
	update_statements.append("UPDATE phage SET Accession = '" + genome_data[7] + "' WHERE PhageID = '" + genome_data[1] + "';")
	update_statements.append("UPDATE phage SET AnnotationAuthor = '" + genome_data[9] + "' WHERE PhageID = '" + genome_data[1] + "';")

	#Create the statement to update Cluster, Cluster2, and Subcluster2
	update_statements.append(\
			create_cluster2_statement(genome_data[1],genome_data[3]))
	update_statements.append(\
			create_subcluster2_statement(genome_data[1],genome_data[8]))
	update_statements.append(\
			create_cluster_statement(genome_data[1],\
			assign_cluster_field(genome_data[8],genome_data[3])))

	updated += 1

#If it looks like there is a problem with some of the genomes on the list,
#cancel the transaction, otherwise proceed
if updated == update_total:

	if run_type == "production":
		con = mdb.connect(mysqlhost, username, password, database)
		con.autocommit(False)
		cur = con.cursor()

		try:
			cur.execute("START TRANSACTION")
			for statement in update_statements:
				cur.execute(statement)
				write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
			cur.execute("COMMIT")
			write_out(output_file,"\nAll update statements committed.")
			cur.close()
			con.autocommit(True)

		except:
			success_action_file_handle.close()
			mdb_exit("\nError: problem updating genome information.\nNo changes have been made to the database.")

		con.close()

	else:
		write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

else:
	write_out(output_file,"\nError: problem processing data list to update genomes. Check input table format.\nNo changes have been made to the database.")
	write_out(output_file,"\nExiting import script.")
	output_file.close()
	success_action_file_handle.close()
	sys.exit(1)

#Document the update actions
for element in update_data_list:

	if element[7] == "":
		element[7] = "none"

	update_output_list = [element[0],\
							element[1],\
							element[2],\
							element[3],\
							element[8],\
							element[4],\
							author_dictionary[element[9]],\
							element[5],\
							element[7],\
							element[10],\
							element[6]]
	success_action_file_writer.writerow(update_output_list)


write_out(output_file,"\nAll field update actions have been implemented.")
raw_input("\nPress ENTER to proceed to next import stage.")








#Remove actions implemented
#Prepare to remove any genomes that are not accompanied by a new genome
write_out(output_file,"\n\n\n\nRemoving genomes with no replacement...")
removed = 0
removal_statements = []
for genome_data in remove_data_list:
	write_out(output_file,"\nPreparing: " + str(genome_data))
	removal_statements.append("DELETE FROM phage WHERE PhageID = '" + genome_data[6] + "';")
	removed += 1

#If it looks like there is a problem with some of the genomes on the list,
#cancel the transaction, otherwise proceed
if removed == remove_total:

	if run_type == "production":
		con = mdb.connect(mysqlhost, username, password, database)
		con.autocommit(False)
		cur = con.cursor()

		try:
			cur.execute("START TRANSACTION")
			for statement in removal_statements:
				cur.execute(statement)
				write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
			cur.execute("COMMIT")
			write_out(output_file,"\nAll remove statements committed.")
			cur.close()
			con.autocommit(True)

		except:
			success_action_file_handle.close()
			mdb_exit("\nError: problem removing genomes with no replacements.\nNo remove actions have been implemented.")

		con.close()

	else:
		write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

else:
	write_out(output_file,"\nError: problem processing data list to remove genomes. Check input table format.\nNo remove actions have been implemented.")
	output_file.close()
	success_action_file_handle.close()
	sys.exit(1)

#Document the remove actions
for element in remove_data_list:

	remove_output_list = [element[0],\
							element[1],\
							element[2],\
							element[3],\
							element[8],\
							element[4],\
							element[9],\
							element[5],\
							element[7],\
							element[10],\
							element[6]]
	success_action_file_writer.writerow(remove_output_list)

write_out(output_file,"\nAll genome remove actions have been implemented.")
raw_input("\nPress ENTER to proceed to next import stage.")











#Add and replace actions implemented
#Now that the data list does not contain update-only or remove-only data,
#create a dictionary. This serves to verify that all genomes
#to be imported are unique, as well as to
#be able to quickly retrieve information based on the genbank file that is opened.
#Key = PhageID
#Value = Genome data list:
add_replace_data_dict = {}
for genome_data in add_replace_data_list:

	#Verify there are no duplicate PhageIDs. This was checked before in a
	#slightly different way when looking for duplicate actions.
	if genome_data[1] not in add_replace_data_dict:
		add_replace_data_dict[genome_data[1]] = genome_data

	else:
		write_out(output_file,"\nError: problem creating genome data dictionary. PhageID %s already in dictionary. Check input table format." % genome_data[1])
		write_out(output_file,"\nNo add/replace actions have been implemented.")
		output_file.close()
		success_action_file_handle.close()
		sys.exit(1)




#Iterate over each file in the directory
admissible_file_types = set(["gb","gbf","gbk","txt"])
write_out(output_file,"\n\n\n\nAccessing genbank-formatted files for add/replace actions...")
all_files =  [X for X in os.listdir(phageListDir) if os.path.isfile(os.path.join(phageListDir,X))]
genbank_files = []
failed_genome_files = []
failed_actions = []
file_tally = 0
script_warnings = 0
script_errors = 0




#Files in the directory may not have the correct extension, or may contain 0, 1, or >1 parseable records.
#When SeqIO parses files, if there are 0 Genbank-formatted records, it does not throw an error, but simply moves on.
#The code first iterates through each file to check how many actually are valid Genbank files.
#It's advantageous to do this first round of file iteration to identify any files with multiple records present,
#and that the file can be successfully read by Biopython SeqIO (if not, it could crash the script).
#Not sure how often there are multiple records present in non-SEA-PHAGES NCBI records, but it is a good idea to verify there is only one record per file before proceeding.
write_out(output_file,"\nA total of %s file(s) present in the directory." % len(all_files))
for filename in all_files:


	#If the file extension is not admissible, then skip. Otherwise, proceed
	if filename.split('.')[-1] not in admissible_file_types:
		failed_genome_files.append(filename)
		write_out(output_file,"\nError: file %s does not have a valid file extension. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")
		continue

	#This try/except clause prevents the code from crashing if there
	#is a problem with a file that Biopython has trouble parsing.
	try:
		#Keep track of how many records Biopython parses
		parsed_records_tally = 0
		for seq_record in SeqIO.parse(os.path.join(phageListDir,filename), "genbank"):
			parsed_records_tally += 1

	except:
		failed_genome_files.append(filename)
		write_out(output_file,"\nError: Biopython is unable to parse file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")
		continue


	if parsed_records_tally == 0:
		failed_genome_files.append(filename)
		write_out(output_file,"\nError: Biopython was unable to parse any records from file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")

	elif parsed_records_tally > 1:
		failed_genome_files.append(filename)
		write_out(output_file,"\nError: Biopython found two records in file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")

	elif parsed_records_tally == 1:
		genbank_files.append(filename)




#Iterate over each genbank-formatted genome file in the directory and parse all the needed data
write_out(output_file,"\nA total of %s file(s) containing Genbank-formatted records will be parsed." % len(genbank_files))
for filename in genbank_files:

	file_tally += 1
	write_out(output_file,"\n\nProcessing file %s: %s" % (file_tally,filename))

	basename = filename.split('.')[0]

	#The file is parsed to grab the header information
	for seq_record in SeqIO.parse(os.path.join(phageListDir,filename), "genbank"):

		matched_by_basename = ''
		add_replace_statements = []
		record_errors = 0
		record_warnings = 0
		geneID_set = set()
		phage_data_list = []



		#Create a list to hold summary info on the genome record:
		#Record Name
		#Record ID
		#Record Definition
		#Record Source
		#Record Organism
		#Organism qualifier of Source feature
		#List of CDS info:
			#Assigned geneID (how the geneID will look when uploaded to the database)
			#Locus tag
			#Product
			#Note
			#Function
			#Translation table
			#Translation
			#Processed Description (how the descriptions will look when uploaded to the database)
		record_summary_header = [["Header Field","Data"]]


		#Parse the PhageID from the Genbank record
		try:
			record_organism = seq_record.annotations["organism"]
			if record_organism.split(' ')[-1] == "Unclassified.":
				phageName = record_organism.split(' ')[-2]
			else:
				phageName = record_organism.split(' ')[-1]

			#PECAAN auto-annotation adds the Draft suffix to the name
			#in this field. Remove this suffix before proceeding.
			if phageName[-6:].lower() == "_draft":
				phageName = phageName[:-6]


			#Some phage names are spelled differently in phagesdb and Phamerator
			#compared to the Genbank record. Convert the parsed name to the
			#expected name in these databases before proceeding.
			if phageName in phage_name_typo_dict.keys():
				incorrect_phageName = phageName
				phageName = phage_name_typo_dict[phageName]
				record_warnings += 1
				write_out(output_file,"\nWarning: parsed phage name %s converted to %s." \
					% (incorrect_phageName,phageName))


		except:
			write_out(output_file,"\nError: problem retrieving phage name in file %s. This file will not be processed." % filename)
			record_errors += 1
			failed_genome_files.append(filename)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue

		#Sequence, length, and GC
		try:
			phageSeq = seq_record.seq.upper()
			seqLength = len(phageSeq)
			seqGC = 100 * (float(phageSeq.count('G')) + float(phageSeq.count('C'))) / float(seqLength)

		except:
			write_out(output_file,"\nError: problem retrieving DNA sequence information in file %s. This file will not be processed." % filename)
			record_errors += 1
			failed_genome_files.append(filename)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue

		#Check DNA sequence for possible errors
		nucleotide_set = set(phageSeq)
		nucleotide_error_set = nucleotide_set - dna_alphabet_set
		if len(nucleotide_error_set) > 0:
			record_warnings += 1
			write_out(output_file,"\nWarning: phage %s contains unexpected nucleotide(s): %s" % (phageName,str(nucleotide_error_set)))
			for element in nucleotide_error_set:
				print "\nThere are %s unexpected %s nucleotides in %s." % (phageSeq.count(element),element,phageName)
			record_errors += question("\nError: problem with DNA sequence in phage %s." % phageName)




		#File header fields are retrieved to be able to check phageName and HostStrain typos
		#The Accession field, with the appended version number, is stored as the record.id
		#The Locus name at the top of the file is stored as the record.name
		#The base accession number, without the version, is stored in the 'accession' annotation list
		try:
			record_name = str(seq_record.name)
		except:
			record_name = ""
			print "\nRecord does not have record Locus information."
			record_errors += question("\nError: problem with header info of file %s." % filename)


		try:
			record_id = str(seq_record.id)
		except:
			record_id = ""
			print "\nRecord does not have record ID information."
			record_errors += question("\nError: problem with header info of file %s." % filename)


		try:
			record_def = str(seq_record.description)
		except:
			record_def = ""
			print "\nRecord does not have record definition information."
			record_errors += question("\nError: problem with header info of file %s." % filename)


		try:
			record_source = str(seq_record.annotations["source"])
		except:
			record_source = ""
			print "\nRecord does not have record source information."
			record_errors += question("\nError: problem with header info of file %s." % filename)


		#Date of the record
		try:
			seq_record_date = seq_record.annotations["date"]
			seq_record_date = datetime.strptime(seq_record_date,'%d-%b-%Y')

		except:
			write_out(output_file,"\nError: problem retrieving date in file %s. This file will not be processed." % filename)
			record_errors += 1
			failed_genome_files.append(filename)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue


		#PhageNotes
		try:
			phageNotes = str(seq_record.annotations["comment"])
		except:
			phageNotes = ""


		#Author list
		try:
			#The retrieved authors can be stored in multiple Reference elements
			record_references_author_list = []
			for reference in seq_record.annotations["references"]:
				if reference.authors != "":
					record_references_author_list.append(reference.authors)
			if len(record_references_author_list) > 0:
				record_author_string = ";".join(record_references_author_list)
			else:
				record_author_string = ""
		except:
			record_author_string = ""

		#Accession
		#This initiates this variable. Since different things affect this variable
		#depending on whether it is an add or replace action, it is easiest to
		#initiate it in advance to avoid throwing an error.
		accession_to_upload = ""


		try:
			#There may be a list of accessions associated with this file.
			#I think the first accession in the list is the most recent.
			#Discard the version suffix if it is present in the
			#Accession field (it might not be present).
			#If an Accession is present, the script will interpret this to mean
			#the genome record is derived from NCBI. An empty Accession will
			#be interpreted as a manual annotation not retrieved from NCBI.
			parsed_accession = seq_record.annotations["accessions"][0]
			parsed_accession = parsed_accession.split('.')[0]
		except:
			parsed_accession = "none"


		#Old code used to match up to the import ticket data
		# if use_basename == "yes":
		#     matchedData = add_replace_data_dict.pop(basename,"error")
		# else:
		#     matchedData = add_replace_data_dict.pop(phageName,"error")




		#Since the Genbank record and the import ticket hasn't been matched yet
		#it is not clear whether it should be matched by filename or PhageID
		#If the basename is the same as the PhageID, it doesn't matter.
		#Otherwise it first checks by PhageID, then the basename.
		#Once matched, it confirms later that the matching strategy indicated in the
		#import ticket's run mode was the same that was actually used.
		if basename == phageName:
			matchedData = add_replace_data_dict.pop(basename,"error")
			matched_by_basename = 'both'

		elif phageName in add_replace_data_dict.keys():
			matchedData = add_replace_data_dict.pop(phageName,"error")
			matched_by_basename = 'no'

		elif basename in add_replace_data_dict.keys():
			matchedData = add_replace_data_dict.pop(basename,"error")
			matched_by_basename = 'yes'

		else:
			matchedData = "error"
			matched_by_basename = 'neither'


		if matchedData != "error":
			write_out(output_file,"\nPreparing: " + str(matchedData))
			import_action = matchedData[0]
			import_host = matchedData[2]
			import_cluster = matchedData[3]
			import_status = matchedData[4]
			import_cds_qualifier = matchedData[5]
			import_genome_replace = matchedData[6]
			import_accession = matchedData[7]
			import_subcluster = matchedData[8]
			import_author = matchedData[9]
			import_run_mode_dict = run_mode_options_dict[matchedData[10]]



			#Assign run mode parameters according to import ticket.
			use_basename = import_run_mode_dict['use_basename']
			custom_gene_id = import_run_mode_dict['custom_gene_id']
			ignore_gene_id_typo = import_run_mode_dict['ignore_gene_id_typo']
			ignore_description_field_check = import_run_mode_dict['ignore_description_field_check']
			ignore_replace_warning = import_run_mode_dict['ignore_replace_warning']
			ignore_trna_check = import_run_mode_dict['ignore_trna_check']
			ignore_locus_tag_import = import_run_mode_dict['ignore_locus_tag_import']
			ignore_phage_name_typos = import_run_mode_dict['ignore_phage_name_typos']
			ignore_host_typos = import_run_mode_dict['ignore_host_typos']
			ignore_generic_author = import_run_mode_dict['ignore_generic_author']
			ignore_description_check = import_run_mode_dict['ignore_description_check']


		else:
			write_out(output_file,"\nError: problem matching phage %s in file %s to genome data from table. This genome was not added. Check input table format." % (phageName,filename))
			record_errors += 1
			failed_genome_files.append(filename)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue



		#Now that the import ticket is matched, verify that the matching strategy
		#indicated in the tickets is equivalent to the strategy actually used
		if use_basename != matched_by_basename and matched_by_basename != 'both':
			record_errors += 1
			write_out(output_file,\
			"\nError: the genome in file %s was not matched to the import ticket as expected in the indicated run mode." \
			% (filename))



		#Check the database to make sure the genome sequence and name matches correctly.
		con = mdb.connect(mysqlhost, username, password, database)
		con.autocommit(False)
		cur = con.cursor()
		try:
			cur.execute("START TRANSACTION")
			cur.execute("""SELECT PhageID,status FROM phage WHERE Sequence = "%s" """ % phageSeq)
			query_results = cur.fetchall()
			cur.execute("SELECT GeneID,PhageID FROM gene")
			current_gene_data_tuples = cur.fetchall()
			cur.execute("COMMIT")
			cur.close()
			con.autocommit(True)


		except:
			success_action_file_handle.close()
			mdb_exit("\nError: retrieving genome information from database while processing file %s.\nNot all genbank files were processed." % filename)

		con.close()


		#Create a set of GeneIDs. If a genome will be replaced,
		#do not add those GeneIDs to the set.
		all_GeneID_set = set()
		for gene_tuple in current_gene_data_tuples:
			if (import_action == "replace" and gene_tuple[1] == import_genome_replace):
				continue
			all_GeneID_set.add(gene_tuple[0])



		#Cross-check the import action against the current state of
		#the database and create SQL statements

		#If adding a new genome, no genome sequence in database is expected
		#to match the current genome sequence
		if (import_action == "add" and len(query_results) > 0):
			record_errors += 1
			write_out(output_file,\
			"\nError: these genome(s) in the database currently contain the same genome sequence as %s: %s." \
			% (phageName,query_results))


		#If replacing a genome:
		elif import_action == "replace":

			#Retrieve current phamerator data (if it is a 'replace' ticket)
			try:
				matched_phamerator_data = phamerator_data_dict[import_genome_replace]
			except:
				matched_phamerator_data = "error"


			if matched_phamerator_data != "error":
				write_out(output_file,"\nRetrieved phamerator data for phage: %s" % phageName)

				phamerator_host = matched_phamerator_data[2]
				phamerator_cluster = matched_phamerator_data[5]
				phamerator_status = matched_phamerator_data[4]
				phamerator_datelastmod = matched_phamerator_data[6]
				phamerator_accession = matched_phamerator_data[7]
				phamerator_subcluster = matched_phamerator_data[8]
				phamerator_author = matched_phamerator_data[9]
				phamerator_annotation_qc = matched_phamerator_data[10]
				phamerator_retrieve_record = matched_phamerator_data[11]

			else:
				write_out(output_file,"\nError: problem matching phage %s in file %s to phamerator data. This genome was not added. Check input table format." % (phageName,filename))
				record_errors += 1
				failed_genome_files.append(filename)
				script_warnings += record_warnings
				script_errors += record_errors
				raw_input("\nPress ENTER to proceed to next file.")
				continue


			#Import and Phamerator host data check
			if import_host != phamerator_host:

				record_warnings += 1
				write_out(output_file,"\nWarning: There is conflicting host data for genome %s" % phageName)
				print "Phamerator host: %s" % phamerator_host
				print "Import ticket host: %s" % import_host
				print "The new host data will be imported."
				record_errors += question("\nError: incorrect host data for %s." % phageName)


			#Import and Phamerator status data check
			if import_status != phamerator_status:

				#It is not common to change from 'gbk' or 'final' to anything else
				if phamerator_status != 'draft':

					record_warnings += 1
					write_out(output_file,"\nWarning: There is conflicting status data for genome %s" % phageName)
					print "Phamerator status: %s" % phamerator_status
					print "Import ticket status: %s" % import_status
					print "The new status data will be imported."
					record_errors += question("\nError: incorrect status data for %s." % phageName)

				#It is common to change status from 'draft' to 'final', but not anything else
				elif import_status != "final":

					record_warnings += 1
					write_out(output_file,"\nWarning: There is conflicting status data for genome %s" % phageName)
					print "Phamerator status: %s" % phamerator_status
					print "Import ticket status: %s" % import_status
					print "The new status data will be imported."
					record_errors += question("\nError: incorrect status data for %s." % phageName)


			#Verify there are no accession conflicts.
			#Phamerator may or may not have accession data.
			#Import table may or may not have accession data.
			#Genbank-formatted file may or may not have accession data.

			accession_comparison_set = set()
			if phamerator_accession != "none":
				accession_comparison_set.add(phamerator_accession)
			if import_accession != "none":
				accession_comparison_set.add(import_accession)
			if parsed_accession != "none":
				accession_comparison_set.add(parsed_accession)


			if len(accession_comparison_set) == 1:
				accession_to_upload = list(accession_comparison_set)[0]
			elif len(accession_comparison_set) > 1:
				record_warnings += 1
				write_out(output_file,"\nWarning: There is conflicting accession data for genome %s" % phageName)
				print "Phamerator accession: %s" % phamerator_accession
				print "Import ticket accession: %s" % import_accession
				print "Parsed accession from file: %s" % parsed_accession
				print "If the parsed accession is not None, it will be imported."
				print "If the parsed accession is None, but the import ticket accession is not None, it will be imported."
				record_errors += question("\nError: incorrect accession data for %s." % phageName)

				if parsed_accession != "none":
					accession_to_upload = parsed_accession
				elif import_accession != "none":
					accession_to_upload = import_accession
				elif phamerator_accession != "none":
					accession_to_upload = phamerator_accession
				else:
					accession_to_upload = ""
			else:
				accession_to_upload = ""



			#Author check
			#It is not common for authorship to change
			if import_author != phamerator_author:

				record_warnings += 1
				write_out(output_file,"\nWarning: there is conflicting author data for genome %s." % phageName)
				print "Phamerator author: %s" % author_dictionary[phamerator_author]
				print "Import ticket author: %s" % author_dictionary[import_author]
				print "The new author data will be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)



			#Exactly one and only one genome in the database is
			#expected to have the same sequence.
			if len(query_results) > 1:
				record_errors += 1
				write_out(output_file,"\nError: the following genomes in the database currently contain the same genome sequence as %s: %s).\nUnable to perform replace action." % (phageName,query_results))

			elif len(query_results) == 0:

				record_warnings += 1
				write_out(output_file,"\nWarning: %s appears to be a different genome sequence than %s. These genomes do not match." % (phageName,import_genome_replace))
				print "The genome will still be replaced."
				record_errors += question("\nError: %s and %s have different genome sequences." % (phageName,import_genome_replace))

			elif len(query_results) == 1:

				#If the new genome is a Final or Gbk, this code block
				#alerts the user, unless this warning is turned off.
				if query_results[0][1].lower() != "draft" and ignore_replace_warning != 'yes':

					record_warnings += 1
					write_out(output_file,"\nWarning: The genome in the database with matching sequence, %s, is listed as %s status." % (query_results[0][0],query_results[0][1]))
					print "The genome will still be replaced."
					record_errors +=  question("\nError: the genome to be removed, %s, was incorrect status." %import_genome_replace)

				#The genome to be replaced does not match the genome
				#name in the database with the same sequence.
				if query_results[0][0] != import_genome_replace:
					write_out(output_file,"\nError: the genome to be removed, %s, does not match the genome name in the database, %s, that has the matching genome sequence to %s." % (import_genome_replace,query_results[0][0],phageName))
					record_errors += 1
			else:
				pass

			#Check to see if the date in the new record is more recent
			#than when the old record was uploaded into Phamerator
			#(stored in DateLastModified)
			if not seq_record_date > phamerator_datelastmod:
				record_warnings += 1
				write_out(output_file,"\nWarning: The date %s in file %s is not more recent than the Phamerator date %s." %(seq_record_date,filename,phamerator_datelastmod))
				print 'Despite it being an older record, the phage %s will continue to be imported.' % phageName
				record_errors +=  question("\nError: the date %s in file %s is not more recent than the Phamerator date %s." %(seq_record_date,filename,phamerator_datelastmod))

			#Create the DELETE command
			add_replace_statements.append("DELETE FROM phage WHERE PhageID = '" + import_genome_replace + "';")

		else:
			#At this point, if the genome is being added, simply assign the accession data
			#from the import table to the accession_to_upload variable. The assignment
			#logic is more complex if the genome is being replaced.
			if import_accession != "none":
				accession_to_upload = import_accession






		#Author list check
		#For annotation author = hatfull
		#If annotation status is draft, author field can be missing Hatfull
		#If annotation status is final, author field should have Hatfull
		#For annotation author = gbk, author should NOT be in either field
		pattern5 = re.compile("hatfull")
		search_result = pattern5.search(record_author_string.lower())


		#For Hatfull authored final annotations, Hatfull is an expected author
		if import_author == '1' and import_status == 'final':

			if search_result == None:
				record_warnings += 1
				write_out(output_file,"\nWarning: Graham Hatfull is not a listed author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)


		#For Hatfull authored draft annotations, it doesn't matter whether
		#there are authors, since they will be added at the Final stage.
		elif import_author == '1' and import_status == 'draft':
			pass

		#For non-Hatfull authored annotations of any kind (draft, final, gbk),
		#Graham may or may not be an author
		else:
			if search_result != None:
				record_warnings += 1
				write_out(output_file,"\nWarning: Graham Hatfull is a listed author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)





		#Check for generic author that sometimes gets added by DNA Master
		if ignore_generic_author != 'yes':
			pattern6 = re.compile("lastname")
			search_result6 = pattern6.search(record_author_string.lower())

			pattern7 = re.compile("firstname")
			search_result7 = pattern7.search(record_author_string.lower())


			if search_result6 != None or search_result7 != None:
				record_warnings += 1
				write_out(output_file,"\nWarning: the author list appears to contain a generic Lastname, Firstname author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)






		#Determine AnnotationQC and RetrieveRecord settings
		#AnnotationQC and RetrieveRecord settings can be carried over
		#from the previous settings in Phamerator under specific
		#circumstances.
		#This enables manual updates made to the database for these two fields
		#to be carried over as the genome is manually or automatically updated.

		#Since a change in AnnotatioQC or RetrieveRecord from default settings
		#is expected to be rare, these two fields can be updated without
		#requiring two additional fields in the import table.

		#For genome replacements, RetrieveRecord is determined based on
		#previous setting. AnnotationQC is determined from previous setting
		#UNLESS it is a Final genome replacing a Draft, in which the AnnotationQC
		#must be changed from 0 to 1.
		if import_action == "replace":

			ncbi_update_status = phamerator_retrieve_record

			#Only under the specific event of a Final replacing a Draft
			#should the AnnotationQC be set to 1. All other events should simply
			#carry it over from the previous setting.
			if phamerator_status == 'draft' and import_status == 'final':
				annotation_qc = '1'

			else:
				annotation_qc = phamerator_annotation_qc

		#For adding a new genome, AnnotationQC is determined by the status,
		#and the RetrieveRecord is determined by the author.
		else:

			if import_author == '1':
				ncbi_update_status = '1'
			else:
				ncbi_update_status = '0'

			if import_status == 'draft' or import_status == 'gbk':
				annotation_qc = '0'
			else:
				annotation_qc = '1'








		#Create list of phage data, then append it to the SQL statement
		#0 = PhageID or basename
		#1 = accession
		#2 = Phage Name
		#3 = import_host
		#4 = phageSeq
		#5 = seqLength
		#6 = seqGC
		#7 = import_status
		#8 = date
		#9 = ncbi_update_status
		#10 = annotation_qc
		#11 = import_author

		if use_basename == "yes":
			phage_data_list.append(basename) #[0]
		else:
			phage_data_list.append(phageName) #[0]



		phage_data_list.append(accession_to_upload) #[1]

		#For the Phage Name, append the Draft suffix if it is draft status
		if use_basename == "yes":
			phage_data_list.append(basename) #[2]
		elif import_status == 'draft':
			phage_data_list.append(phageName + "_Draft") #[2]
		else:
			phage_data_list.append(phageName) #[2]

		phage_data_list.append(import_host) #[3]
		phage_data_list.append(str(phageSeq)) #[4]
		phage_data_list.append(seqLength) #[5]
		phage_data_list.append(round(seqGC,4)) #[6] Trim down to last 4 digits
		phage_data_list.append(import_status) #[7]
		phage_data_list.append(date) #[8]
		phage_data_list.append(ncbi_update_status) #[9]
		phage_data_list.append(annotation_qc) #[10]
		phage_data_list.append(import_author) #[11]

		add_replace_statements.append("""INSERT INTO phage (PhageID, Accession, Name, HostStrain, Sequence, SequenceLength, GC,status, DateLastModified, RetrieveRecord, AnnotationQC, AnnotationAuthor) VALUES ("%s","%s","%s","%s","%s",%s,%s,"%s","%s","%s","%s","%s")""" \
										% (phage_data_list[0],\
										phage_data_list[1],\
										phage_data_list[2],\
										phage_data_list[3],\
										phage_data_list[4],\
										phage_data_list[5],\
										phage_data_list[6],\
										phage_data_list[7],\
										phage_data_list[8],\
										phage_data_list[9],\
										phage_data_list[10],\
										phage_data_list[11]))


		add_replace_statements.append(\
				create_cluster_statement(phage_data_list[0],\
				assign_cluster_field(import_subcluster,import_cluster)))
		add_replace_statements.append(\
				create_cluster2_statement(phage_data_list[0],import_cluster))
		add_replace_statements.append(\
				create_subcluster2_statement(phage_data_list[0],import_subcluster))


		#Next each CDS feature is parsed from the file
		#The cdsCount will increment for each CDS processed, even if it does not pass the QC filters. This way, the genome still retains the info that another CDS was originally present.
		cdsCount = 0
		addCount = 0
		transl_table_set = set()
		missing_transl_table_tally = 0
		missing_locus_tag_tally = 0
		assigned_description_tally = 0
		feature_note_tally = 0
		feature_product_tally = 0
		feature_function_tally = 0
		feature_source_organism = ""
		feature_source_lab_host = ""
		feature_source_host = ""
		all_features_data_list = []
		all_coordinates_set = set()
		record_summary_cds = [["Locus Tag",\
								"Product",\
								"Function",\
								"Note",\
								"Translation Table",\
								"Translation",\
								"Assigned GeneID",\
								"Assigned Description"]]

		for feature in seq_record.features:


			if feature.type != "CDS":

				#Retrieve the Source Feature info
				if feature.type == "source":

					try:
						feature_source_organism = str(feature.qualifiers["organism"][0])
					except:
						feature_source_organism = ""

					try:
						feature_source_host = str(feature.qualifiers["host"][0])
					except:
						feature_source_host = ""

					try:
						feature_source_lab_host = str(feature.qualifiers["lab_host"][0])
					except:
						feature_source_lab_host = ""

				#Evaluate if tRNA is properly structured
				#Only evaluate if it is a new manually annotated SEA-PHAGES genomes
				#that is replacing an auto-annotated genome.
				#Eventually, this restriction can be loosened to evaluate tRNAs
				#in any type of genome.
				elif feature.type == "tRNA" and ignore_trna_check != "yes":


					#Retrieve tRNA coordinates
					try:

						#Biopython converts coordinates to 0-index
						#Start(left) coordinates are 0-based inclusive (feature starts there)
						#Stop (right) coordinates are 0-based exclusive (feature stops 1bp prior to coordinate)
						tRNA_left = str(feature.location.start)
						tRNA_right = str(feature.location.end)

					except:
						write_out(output_file,"\nError: a tRNA has incorrect coordinates in phage %s."\
								% phageName)
						record_errors += 1
						continue

					#Now that start and stop have been parsed, check if coordinates are fuzzy or not
					if (tRNA_left.isdigit() and tRNA_right.isdigit()):
						tRNA_left = int(tRNA_left)
						tRNA_right = int(tRNA_right)
					else:
						write_out(output_file,"\nError: tRNA starting at %s has fuzzy coordinates in phage %s."\
								% (tRNA_left,phageName))
						record_errors += 1
						continue

					#Retrieve top strand of tRNA feature. It is NOT necessarily
					#in the correct orientation
					tRNA_size = abs(tRNA_right - tRNA_left)
					tRNA_seq = phageSeq[tRNA_left:tRNA_right].upper()
					if len(tRNA_seq) != tRNA_size:
						write_out(output_file,"\nError: unable to retrieve sequence for tRNA starting at %s in phage %s."\
								% (tRNA_left + 1,phageName))
						record_errors += 1
						continue


					#Convert sequence to reverse complement if it is on bottom strand
					if feature.strand == 1:
						pass
					elif feature.strand == -1:
						tRNA_seq = tRNA_seq.reverse_complement()
					else:
						record_errors += 1
						write_out(output_file,"\Error: tRNA starting at %s does not have proper orientation in %s phage." \
								% (tRNA_left + 1,phageName))
						continue

					#Check to see if forward strand terminal nucleotide is correct = A or C
					if tRNA_seq[-1] != 'A' and tRNA_seq[-1] != 'C':
						record_warnings += 1
						write_out(output_file,"\nWarning: tRNA starting at %s does not appear to have correct terminal nucleotide in %s phage." \
								% (tRNA_left + 1,phageName))
						record_errors += question("\nError: tRNA starting at %s has incorrect terminal nucleotide in %s phage." \
								% (tRNA_left + 1,phageName))

					if tRNA_size < 60 or tRNA_size > 100:
						record_warnings += 1
						write_out(output_file,"\nWarning: tRNA starting at %s does not appear to be the correct size in %s phage."  \
								% (tRNA_left + 1,phageName))
						record_errors += question("\nError: tRNA starting at %s is incorrect size in %s phage." \
								% (tRNA_left + 1,phageName))


					#Retrieve and check product
					try:
						tRNA_product = feature.qualifiers['product'][0].lower().strip()

						if len(tRNA_product) > 0:
							if check_tRNA_product(tRNA_product) > 0:
								write_out(output_file,"\nError: tRNA starting at %s has incorrect amino acid or anticodon in %s." \
									% (tRNA_left + 1, phageName))
								record_errors += 1
						else:
							write_out(output_file,"\nError: tRNA starting at %s has incorrect product in %s." \
								% (tRNA_left + 1, phageName))
							record_errors += 1

					except:
						write_out(output_file,"\nError: tRNA starting at %s is missing product field in phage %s." \
							% (tRNA_left + 1,phageName))
						record_errors += 1
						tRNA_product = ''


					#Retrieve note
					#In the future, this field may need to be parsed in a similar
					#manner as the product field. For now, do nothing.
					try:
						tRNA_note = feature.qualifiers['note'][0].lower().strip()
					except:
						tRNA_note = ''



				#If feature is not CDS, Source, or tRNA, skip it
				else:
					pass

				continue

			else:
				cdsCount += 1
				typeID = feature.type

			#This will store all data for this feature that will be imported
			feature_data_list = []


			#GeneID
			#Feature_locus_tag is a record of the locus tag found in the file.
			#GeneID is what will be assigned in the database.
			try:
				feature_locus_tag = feature.qualifiers["locus_tag"][0]
			except:
				feature_locus_tag = ""
				missing_locus_tag_tally += 1

			#If user selected customized GeneIDs, OR if there is no locus tag,
			#then create a concatenated GeneID
			if custom_gene_id == 'yes' or feature_locus_tag == '':
				if use_basename == 'yes':
					geneID = basename.upper() + "_" + str(cdsCount)
				else:
					geneID = phageName.upper() + "_" + str(cdsCount)
			else:
				geneID = feature_locus_tag




			#See if the geneID is already in the database
			duplicate = False
			if geneID in geneID_set:
				duplicate = True
				record_warnings += 1
				write_out(output_file,"\nWarning: there is a duplicate geneID %s in phage %s." % (geneID,phageName))

			elif geneID in all_GeneID_set:
				duplicate = True
				record_warnings += 1
				write_out(output_file,"\nWarning: there is a duplicate geneID %s in the current database." % geneID)
			else:
				geneID_set.add(geneID)

			#If there is a geneID duplication conflict, try to resolve it
			old_ID = geneID
			dupe_value = 1
			while duplicate == True and dupe_value < 99:

				geneID = old_ID + "_duplicateID" + str(dupe_value)
				record_warnings += 1
				#Check to see if the new geneID is found in the set of all geneIDs
				if (geneID not in geneID_set and geneID not in all_GeneID_set):
					duplicate = False
					write_out(output_file,"\nGeneID %s duplication has been automatically resolved by renaming ID to %s." % (old_ID,geneID))
					duplicate_answer = question("\nError: feature %s of %s is a duplicate geneID." % (old_ID,phageName))

					#If user indicates the feature with the new geneID should not be added, add to record_errors. Otherwise, assign new geneID to the geneID_set
					if duplicate_answer == 1:
						record_errors += 1
					else:
						geneID_set.add(geneID)
				dupe_value += 1

			#Once the while loop exits, check if the duplication was resolved
			if duplicate == True:
				record_warnings += 1
				write_out(output_file,"\nWarning: unable to resolve duplicate geneID %s conflict. This CDS will be skipped, but processing of the other genes will continue." % old_ID)
				record_errors += question("\nError: feature %s of phage %s is a duplicate geneID and cannot be renamed to %s." % (old_ID,phageName,geneID))
				continue





			#Name
			if (geneID.split('_')[-1].isdigit()):
				geneName = geneID.split('_')[-1]
			else:
				geneName = cdsCount



			#Orientation
			if feature.strand == 1:
				orientation = "Forward"
			elif feature.strand == -1:
				orientation = "Reverse"
			#ssRNA phages
			elif feature.strand is None:
				orientation = "Forward"
			else:
				record_warnings += 1
				write_out(output_file,"\nWarning: feature %s of %s does not have a common orientation. This CDS will be skipped, but processing of the other genes will continue." % (geneID,phageName))
				record_errors += question("\nError: feature %s of %s does not have correct orientation." % (geneID,phageName))
				continue






			#Gene boundary coordinates
			#Compound features are tricky to parse.
			if str(feature.location)[:4] == "join":


				#Skip this compound feature if it is comprised of more than two features (too tricky to parse).
				if len(feature.location.parts) > 2:

					strStart = ""
					strStop = ""
					record_warnings += 1
					write_out(output_file,"\nWarning: gene %s is a compound feature that is unable to be parsed. This CDS will be skipped, but processing of the other genes will continue." % geneID)
					record_errors += question("\nError: unable to parse gene %s of phage %s." % (geneID,phageName))
					continue

				else:

					#Retrieve compound feature positions based on strand
					if feature.strand == 1:

						strStart = str(feature.location.parts[0].start)
						strStop = str(feature.location.parts[1].end)

					elif feature.strand == -1:

						strStart = str(feature.location.parts[1].start)
						strStop = str(feature.location.parts[0].end)

					#If strand is None...
					else:
						strStart = ""
						strStop = ""

			else:
				strStart = str(feature.location.start)
				strStop = str(feature.location.end)

			#Now that start and stop have been parsed, check if coordinates are fuzzy or not
			if (strStart.isdigit() and strStop.isdigit()):
				startCoord = int(strStart)
				stopCoord = int(strStop)
			else:
				record_warnings += 1
				write_out(output_file,"\nWarning: gene %s start %s and stop %s are non-traditional coordinates. This CDS will be skipped, but processing of the other genes will continue." % (geneID,strStart,strStop))
				record_errors += question("\nError: feature %s of %s does not have correct coordinates." % (geneID,phageName))
				continue

			#Test if there is a gene with the same coordinates already parsed.
			coordinate_tuple = tuple([startCoord,stopCoord,orientation])
			if coordinate_tuple not in all_coordinates_set:
				all_coordinates_set.add(coordinate_tuple)
			else:
				record_warnings += 1
				write_out(output_file,"\nWarning: multiple genes have coordinates %s. This is likely a gene feature duplication." % str(coordinate_tuple))
				record_errors += question("\nError: gene coordinates %s are duplicated in this genome." % str(coordinate_tuple))









			#Translation, Gene Length (via Translation)
			try:
				translation = feature.qualifiers["translation"][0].upper()
				geneLen = (len(translation) * 3) + 3  #Add 3 for the stop codon...
			except:
				translation = ""
				geneLen = 0
				record_warnings += 1
				write_out(output_file,"\nWarning: gene %s has no translation. This CDS will be skipped, but processing of the other genes will continue." % geneID)
				record_errors += question("\nError: problem with %s translation in phage %s." % (geneID,phageName))
				continue

			#Check translation for possible errors
			amino_acid_set = set(translation)
			amino_acid_error_set = amino_acid_set - protein_alphabet_set
			if len(amino_acid_error_set) > 0:
				record_warnings += 1
				write_out(output_file,"\nWarning: feature %s of %s appears to have unexpected amino acid(s)." % (geneID,phageName))
				print "Unexpected amino acids: " + str(amino_acid_error_set)
				record_errors += question("\nError: problem with %s translation in phage %s." % (geneID,phageName))

			#Translation table used
			try:
				feature_transl_table = feature.qualifiers["transl_table"][0]
				transl_table_set.add(feature_transl_table)
			except:
				feature_transl_table = ""
				missing_transl_table_tally += 1






			#Gene Description
			#Use the feature qualifier from the import file.
			#If it is a generic description, leave field empty.
			#For generic 'gp#' descriptions, make sure there is no other info in the product that is valuabe by checking for other words present.


			#First retrieve gene function, note, and product fields.
			try:
				feature_product = retrieve_description(feature,"product")
				if feature_product != "":
					feature_product_tally += 1
			except:
				feature_product = ""

			try:
				feature_function = retrieve_description(feature,"function")
				if feature_function != "":
					feature_function_tally += 1
			except:
				feature_function = ""

			try:
				feature_note = retrieve_description(feature,"note")
				if feature_note != "":
					feature_note_tally += 1
			except:
				feature_note = ""




			#Now assign the appropriate description info to the
			#assigned_description variable, as indicated from the import table.
			try:

				if import_cds_qualifier == "product":
					assigned_description = feature_product

				elif import_cds_qualifier == "function":
					assigned_description = feature_function

				elif import_cds_qualifier == "note":
					assigned_description = feature_note

				#This clause allows the user to specify an uncommon
				#feature qualifier to retrieve the gene description from.
				else:
					assigned_description = retrieve_description(feature,import_cds_qualifier)

			except:
				assigned_description = ""

			if assigned_description != "":
				assigned_description_tally += 1



			#Check to verify misc details about the assigned gene description field
			#This is distinct from the ignore_description_field_check, which
			#is used to determine which field (product, function, note) is parsed.

			if ignore_description_check != 'yes':
				if assigned_description[-1:] == '.':

					record_warnings += 1
					write_out(output_file,"\nWarning: the description for feature %s of %s contains a period." % (geneID,phageName))
					print assigned_description
					record_errors += question("\nError: problem with description for feature %s of %s." % (geneID,phageName))









			#Now that it has acquired all gene feature info, create list of
			#gene data and append to list of all gene feature data
			#0 = geneID
			#1 = PhageID or basename
			#2 = startCoord
			#3 = stopCoord
			#4 = geneLen
			#5 = geneName
			#6 = typeID
			#7 = translation
			#8 = orientation[0]
			#9 = assigned_description
			#10 = feature_product
			#11 = feature_function
			#12 = feature_note
			#13 = feature_locus_tag
			addCount+= 1

			feature_data_list.append(geneID) #0


			if use_basename == "yes":
				feature_data_list.append(basename) #1
			else:
				feature_data_list.append(phageName) #1

			feature_data_list.append(startCoord) #2
			feature_data_list.append(stopCoord) #3
			feature_data_list.append(geneLen) #4
			feature_data_list.append(geneName) #5
			feature_data_list.append(typeID) #6
			feature_data_list.append(translation) #7
			feature_data_list.append(orientation[0]) #8
			feature_data_list.append(assigned_description) #9
			feature_data_list.append(feature_product) #10
			feature_data_list.append(feature_function) #11
			feature_data_list.append(feature_note) #12

			#If there is a parsed accession, the file is interpreted to have
			#been retrieved from NCBI, which means the locus tags are the official
			#locus tags. Otherwise, do not retain the locus tags in the record.
			if ignore_locus_tag_import != 'yes':
				feature_data_list.append(feature_locus_tag) #13
			else:
				feature_data_list.append('') #13

			all_features_data_list.append(feature_data_list)


			#Retrieve summary info to verify quality of file
			feature_product_trunc = feature_product
			if len(feature_product_trunc) > 15:
				feature_product_trunc = feature_product_trunc[:15] + "..."

			feature_note_trunc = feature_note
			if len(feature_note_trunc) > 15:
				feature_note_trunc = feature_note_trunc[:15] + "..."

			feature_function_trunc = feature_function
			if len(feature_function_trunc) > 15:
				feature_function_trunc = feature_function_trunc[:15] + "..."

			translation_trunc = translation
			if len(translation_trunc) > 5:
				translation_trunc = translation_trunc[:5] + "..."

			assigned_description_trunc = assigned_description
			if len(assigned_description_trunc) > 15:
				assigned_description_trunc = assigned_description_trunc[:15] + "..."

			record_summary_cds.append([feature_locus_tag,\
										feature_product_trunc,\
										feature_function_trunc,\
										feature_note_trunc,\
										feature_transl_table,\
										translation_trunc,\
										geneID,\
										assigned_description_trunc])






		#Now that all CDS features are processed, run quality control on the data:


		#Check to see if there are any CDS features processed. If not, then the genbank record does not have any called genes.
		#The record_summary_cds list contains the column headers, so at minimum, it is length == 1
		if len(record_summary_cds) == 1:
			print "\nNo CDS features were found in this record. The genome will still be added to the database."
			record_errors += question("\nError: no CDS features found in %s." % filename)


		#Process the source and organism fields to look for problems

		#Print the summary of the header information
		record_summary_header.append(["Record Name",record_name])
		record_summary_header.append(["Record ID",record_id])
		record_summary_header.append(["Record Defintion",record_def])
		record_summary_header.append(["Record Source",record_source])
		record_summary_header.append(["Record Organism",record_organism])
		record_summary_header.append(["Source Feature Organism",feature_source_organism])
		record_summary_header.append(["Source Feature Host",feature_source_host])
		record_summary_header.append(["Source Feature Lab Host",feature_source_lab_host])
		print "\nSummary of record header information for %s from file %s:\n" % (phageName,filename)
		print tabulate(record_summary_header,headers = "firstrow")
		print "\n\n\n"





		#See if there are any phage name typos in the header block
		if ignore_phage_name_typos != 'yes':


			pattern1 = re.compile('^' + phageName + '$')
			pattern2 = re.compile('^' + phageName)


			#REVIEW This QC check may not be necessary. It really doesn't matter
			#for Draft genomes, since the fields are auto-populated. For NCBI genomes,
			#the field should be an Accession. For Final genomes yet to be submitted
			#to NCBI, even if there is a typo, it should not be retained since this
			#field gets populated with the Accession.
			# if find_name(pattern1,record_name.split(' ')) == 0:
			#
			#     if record_name.split('.')[0] != parsed_accession:
			#         print "\nRecord name does not have the accession number or the identical phage name as found in the record organism field."
			#         record_errors += question("\nError: problem with header info of file %s." % filename)


			#Record definition QC
			#It can contain ", complete genome." or "." at the end,
			#so remove this before doing search.
			if record_def[-1:].lower() == '.':
				if record_def[-18:].lower() == ', complete genome.':
					record_def_trimmed = record_def[:-18]
				else:
					record_def_trimmed = record_def[:-1]
			else:
				record_def_trimmed = record_def


			if find_name(pattern1,record_def_trimmed.split(' ')) == 0:
				print "\nRecord definition does not have identical phage name as found in the record organism field."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			#REVIEW the above Record definition check replaces the code block below
			# if find_name(pattern1,record_def[:-18].split(' ')) == 0:
			#     print record_def[:-18].split(' ')
			#     print "\nRecord definition does not have identical phage name as found in the record organism field."
			#     record_errors += question("\nError: problem with header info of file %s." % filename)
			#
			# else:
			#     if find_name(pattern2,record_def.split(' ')) == 0:
			#         print "\nRecord definition does not have identical phage name as found in the record organism field."
			#         record_errors += question("\nError: problem with header info of file %s." % filename)

			#Record source QC
			if find_name(pattern1,record_source.split(' ')) == 0:
				print "\nRecord source does not have identical phage name as found in the record organism field."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			#Source feature organism QC
			if find_name(pattern1,feature_source_organism.split(' ')) == 0:
				print "\nSource feature organism does not have identical phage name as found in the record organism field."
				record_errors += question("\nError: problem with header info of file %s." % filename)







		#See if there are any host name typos in the header block.
		#Skip this step if it is a Draft genome, because it won't correctly have this information.
		if import_status != 'draft' and ignore_host_typos != 'yes':
			import_host_trim = import_host
			if import_host_trim == "Mycobacterium":
				import_host_trim = import_host_trim[:-3]

			pattern3 = re.compile('^' + import_host_trim)


			if (find_name(pattern3,record_def.split(' ')) == 0 and record_def.split(' ')[0].lower() not in host_ignore):
				print "\nRecord definition does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			if (find_name(pattern3,record_source.split(' ')) == 0 and record_source.split(' ')[0].lower() not in host_ignore):

				print "\nRecord source does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			if (find_name(pattern3,record_organism.split(' ')) == 0 and record_organism.split(' ')[0].lower() not in host_ignore):

				print "\nRecord organism does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			if (find_name(pattern3,feature_source_organism.split(' ')) == 0 and feature_source_organism.split(' ')[0].lower() not in host_ignore):

				print "\nSource feature organism does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			#Host and Lab_Host data may not have been present, so skip if it is blank
			if (feature_source_host != "" and find_name(pattern3,feature_source_host.split(' ')) == 0 and feature_source_host.split(' ')[0].lower() not in host_ignore):

				print "\nSource feature host does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)

			if (feature_source_lab_host != "" and find_name(pattern3,feature_source_lab_host.split(' ')) == 0 and feature_source_lab_host.split(' ')[0].lower() not in host_ignore):

				print "\nSource feature lab host does not appear to have same host data as found in import table."
				record_errors += question("\nError: problem with header info of file %s." % filename)



		#Print record summary for all CDS information for quality control
		print "\nSummary of CDS summary information for %s from file %s:\n" % (phageName,filename)
		print tabulate(record_summary_cds,headers = "firstrow")
		print "\n\n\n"


		#Check locus tag info:
		if missing_locus_tag_tally > 0:
			record_warnings += 1
			write_out(output_file,"\nWarning: phage %s from file %s is missing %s CDS locus tag(s)." % (phageName, filename, missing_locus_tag_tally))
			record_errors += question("\nError: problem with locus tags in file  %s." % filename)


		#Check the phage name spelling in the locus tags.
		pattern4 = re.compile(phageName.lower())
		geneID_typo_tally = 0
		geneID_typo_list = []

		if ignore_gene_id_typo != "yes":
			for geneID in geneID_set:

			   search_result = pattern4.search(geneID.lower())
			   if search_result == None:
					geneID_typo_tally += 1
					geneID_typo_list.append(geneID)

			if geneID_typo_tally > 0:
				record_warnings += 1
				write_out(output_file,"\nWarning: there are %s geneID(s) that do not have the identical phage name included." % geneID_typo_tally)
				print geneID_typo_list
				record_errors += question("\nError: problem with locus tags of file %s." % filename)





		#Check all translation table info:
		if len(transl_table_set) > 1:
			write_out(output_file,"\nError: more than one translation table used in file %s." % filename)
			record_errors += 1

		elif len(transl_table_set) == 1:
			transl_table_list = list(transl_table_set)
			if transl_table_list[0] != '11':
				write_out(output_file,"\nThe translation table used for %s is: %s." % (phageName,transl_table_list[0]))
				record_errors += question("\nError: phage %s does not use correct translation table." % phageName)
		else:
			pass

		if missing_transl_table_tally > 0:
			record_warnings += 1
			write_out(output_file,"\nWarning: there are %s genes with no translation table for phage %s." % (missing_transl_table_tally,phageName))
			record_errors += question("\nError: phage %s has missing translation table information." % phageName)



		#Check to ensure the best gene description field was retained
		#Element indices for feature data:
		#0 = geneID
		#1 = PhageID or basename
		#2 = startCoord
		#3 = stopCoord
		#4 = geneLen
		#5 = geneName
		#6 = typeID
		#7 = translation
		#8 = orientation[0]
		#9 = assigned_description
		#10 = feature_product
		#11 = feature_function
		#12 = feature_note

		if import_cds_qualifier not in description_set:
			write_out(output_file,"\nNumber of gene %s descriptions found for phage %s: %s" % (import_cds_qualifier,phageName, assigned_description_tally))
		write_out(output_file,"\nNumber of gene product descriptions found for phage %s: %s" % (phageName, feature_product_tally))
		write_out(output_file,"\nNumber of gene function descriptions found for phage %s: %s" % (phageName, feature_function_tally))
		write_out(output_file,"\nNumber of gene note descriptions found for phage %s: %s" % (phageName, feature_note_tally))


		#If other CDS fields contain descriptions, they can be chosen to
		#replace the default import_cds_qualifier descriptions.
		#Then provide option to verify changes.
		#This block is skipped if user selects to do so.
		if ignore_description_field_check != 'yes':


			changed = ""


			if (import_cds_qualifier != "product" and feature_product_tally > 0):

			   print "\nThere are %s CDS products found." % feature_product_tally
			   change_descriptions()

			   if question("\nCDS products will be used for phage %s in file %s." % (phageName,filename)) == 1:

					for feature in all_features_data_list:
						feature[9] = feature[10]
					changed = "product"

			if (import_cds_qualifier != "function" and feature_function_tally > 0):

				print "\nThere are %s CDS functions found." % feature_function_tally
				change_descriptions()


				if question("\nCDS functions will be used for phage %s in file %s." % (phageName,filename)) == 1:

					for feature in all_features_data_list:
						feature[9] = feature[11]
					changed = "function"


			if (import_cds_qualifier != "note" and feature_note_tally > 0):

				print "\nThere are %s CDS notes found." % feature_note_tally
				change_descriptions()

				if question("\nCDS notes will be used for phage %s in file %s." % (phageName,filename)) == 1:

					for feature in all_features_data_list:
						feature[9] = feature[12]
					changed = "note"

			if changed != "":
				record_warnings += 1
				write_out(output_file,"\nWarning: CDS descriptions only from the %s field will be retained." % changed)
				record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)


		#Add all updated gene feature data to the add_replace_statements list
		for feature in all_features_data_list:
			add_replace_statements.append("""INSERT INTO gene (GeneID, PhageID, Start, Stop, Length, Name, TypeID, translation, Orientation, Notes, LocusTag) VALUES ("%s","%s",%s,%s,%s,"%s","%s","%s","%s","%s","%s");""" \
									% (feature[0],\
									feature[1],\
									feature[2],\
									feature[3],\
									feature[4],\
									feature[5],\
									feature[6],\
									feature[7],\
									feature[8],\
									feature[9],\
									feature[13]))


		#If errors were encountered with the file parsing, do not add to the genome. Otherwise, proceed.
		if record_errors == 0:
			write_out(output_file,"\nNo errors encountered while parsing: %s." % filename)
			diffCount = cdsCount - addCount
			if record_warnings > 0:
				write_out(output_file,"\nWarning summary: there are %s warning(s) with phage %s and %s CDS feature(s) were not added." % (record_warnings,phageName,diffCount))
		else:
			write_out(output_file,"\n%s errors were encountered with file %s. %s genome was not added to the database." % (record_errors,filename,phageName))
			failed_genome_files.append(filename)
			failed_actions.append(matchedData)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue


		#Execute SQL transactions
		if run_type == "production":
			con = mdb.connect(mysqlhost, username, password, database)
			con.autocommit(False)
			cur = con.cursor()
			try:
				cur.execute("START TRANSACTION")
				for statement in add_replace_statements:
					cur.execute(statement)
					write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
				cur.execute("COMMIT")
				write_out(output_file,"\nAll add/replace statements for %s committed." % matchedData)
				cur.close()
				con.autocommit(True)

			except:
				success_action_file_handle.close()
				mdb_exit("\nError: problem importing the file %s with the following add/replace action: %s.\nNot all genbank files were processed." % (filename,matchedData))

			con.close()

		else:
			write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

		#Now that genome has been successfully uploaded, proceed
		try:
			shutil.move(os.path.join(phageListDir,filename),os.path.join(phageListDir,success_folder,filename))
		except:
			print "Unable to move file %s to success file folder." % filename


		#Add the action data to the success output file, update tally of total script warnings and errors, then proceed
		add_replace_output_list = [matchedData[0],\
								matchedData[1],\
								matchedData[2],\
								matchedData[3],\
								matchedData[8],\
								matchedData[4],\
								author_dictionary[matchedData[9]],\
								matchedData[5],\
								matchedData[7],\
								matchedData[10],\
								matchedData[6]]


		success_action_file_writer.writerow(add_replace_output_list)
		script_warnings += record_warnings
		script_errors += record_errors
		print "Processing of %s is complete." %filename



write_out(output_file,"\nAll files have been iterated through.")
raw_input("\nPress ENTER to proceed to next import stage.")


#Final verifications
write_out(output_file,"\n\n\n\nFinal checks...")
write_out(output_file,"\n\nAll update actions have been implemented.")
write_out(output_file,"\n\nAll remove actions have been implemented.")

#If there were failed actions during the genome upload stage, add it back to the add_replace_data_dictionary.
#The add_replace_data_dictionary now contains a list of failed actions as well as unaddressed actions (actions that did not have matching genome files)
if len(failed_actions) > 0:
	for genome_data in failed_actions:
		add_replace_data_dict[genome_data[1]] = genome_data


#Verify all add/replace actions from the import csv file were addressed.
if len(add_replace_data_dict) > 0:

	failed_action_file = '%s_failed_import_table_actions.csv' % date
	failed_action_file_handle = open(os.path.join(phageListDir,failed_folder,failed_action_file),"w")
	failed_action_file_writer = csv.writer(failed_action_file_handle)
	write_out(output_file,"\n\nThe following add/replace action(s) in the import table were NOT successfully implemented:")

	for key in add_replace_data_dict:

		failed_output_list = [add_replace_data_dict[key][0],\
								add_replace_data_dict[key][1],\
								add_replace_data_dict[key][2],\
								add_replace_data_dict[key][3],\
								add_replace_data_dict[key][8],\
								add_replace_data_dict[key][4],\
								author_dictionary[add_replace_data_dict[key][9]],\
								add_replace_data_dict[key][5],\
								add_replace_data_dict[key][7],\
								add_replace_data_dict[key][10],\
								add_replace_data_dict[key][6]]

		failed_action_file_writer.writerow(failed_output_list)
		write_out(output_file,"\n" + str(failed_output_list))

	failed_action_file_handle.close()

else:
	write_out(output_file,"\nAll add/replace actions have been implemented.")


#Verify that none of the genbank files failed to be uploaded.
if len(failed_genome_files) > 0:

	write_out(output_file,"\n\nThe following genbank files were NOT successfully processed:")
	for filename in failed_genome_files:
		write_out(output_file,"\n" + filename)

		try:
			shutil.move(os.path.join(phageListDir,filename),os.path.join(phageListDir,failed_folder,filename))
		except:
			print "Unable to move file %s to failed file folder." % filename


else:
	write_out(output_file,"\n\nAll genbank files were successfully processed.")



#Output the total number of warnings and errors encountered during the import process
write_out(output_file,"\n\nTotal number of warnings encountered: %s" % script_warnings)
write_out(output_file,"\nTotal number of errors encountered: %s" % script_errors)


#Retrieve the total number of genomes now in the updated database
try:
	con = mdb.connect(mysqlhost, username, password, database)
	con.autocommit(False)
	cur = con.cursor()
except:
	print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password.\nImport script was not completed."
	output_file.close()
	success_action_file_handle.close()
	sys.exit(1)

try:
	cur.execute("START TRANSACTION")
	cur.execute("SELECT count(*) FROM phage")
	final_tally = cur.fetchall()
	cur.execute("COMMIT")
	cur.close()
	con.autocommit(True)

except:
	output_file.close()
	success_action_file_handle.close()
	mdb_exit("\nUnable to access the database to retrieve genome information.\nImport script was not completed.")

write_out(output_file,"\n\nTotal phages in database after changes: " + str(final_tally[0][0]))


#Close script.
success_action_file_handle.close()
write_out(output_file,"\n\n\n\nImport script completed.")
output_file.close()
