



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
genome_data_list = []
table_errors = 0
add_total = 0
remove_total = 0
replace_total = 0
update_total = 0
run_mode_custom_total = 0

for input_row in file_reader:


	#Verify the row of information has the correct number of fields to parse.
	if len(input_row) != 11:
		write_out(output_file,"\nRow in import table is not formatted correctly: " + str(input_row))
		table_errors += 1
		continue


	#Once the row length has been verified, rearrange the elements into another order.
	#This is a straightforward, but not ideal, solution to when columnes are added or removed from the import table.
	#This prevents the entire code for needing to be re-factored to account for the new indices.
	#Internal import row structure:
	#0 = Import action (unchanged)
	#1 = New PhageID (unchanged)
	#2 = Host (unchanged)
	#3 = Cluster (unchanged)
	#4 = Status
	#5 = Feature field (unchanged)
	#6 = PhageID to be removed
	#7 = Accession
	#8 = Subcluster
	#9 = AnnotationAuthor
	#10 = Run Mode
	row = []
	row.append(input_row[0]) #Import action
	row.append(input_row[1]) #New PhageID
	row.append(input_row[2]) #Host
	row.append(input_row[3]) #Cluster
	row.append(input_row[5]) #Status
	row.append(input_row[7]) #Feature field
	row.append(input_row[10]) #PhageID to be removed
	row.append(input_row[8]) #Accession
	row.append(input_row[4]) #Subcluster
	row.append(input_row[6]) #AnnotationAuthor
	row.append(input_row[9]) #Run mode








	#Modify fields if needed

	#Modify Host if needed
	if row[2] == "retrieve":

		#On phagesdb, phages should always have the Genus info of the isolation host.
		try:
			row[2] = online_data_dict['isolation_host']['genus']
		except:
			write_out(output_file,"\nError: unable to retrieve Host data for phage %s from phagesdb." %row[1])
			row[2] = "none"
			table_errors += 1

	if row[2] != "none":
		row[2] = row[2].split(' ')[0] #Keep only the genus in the host data field and discard the rest
		if row[2] not in phageHost_set:
			print "The host strain %s is not currently in the database." % row[2]
			table_errors +=  question("\nError: %s is not the correct host for %s." % (row[2],row[1])) #errors will be incremented if host was not correct


	#Modify Cluster and Subcluster if needed
	if row[3] == "retrieve":
		try:

			#On phagesdb, phages may have a Cluster and no Subcluster info (which is set to None).
			#If the phage has a Subcluster, it should also have a Cluster.
			#If by accident no Cluster or Subcluster info is added at the time the
			#genome is added to phagesdb, the Cluster may automatically be set to
			#"Unclustered". This will be filtered out later in the script due to its character length.

			#Retrieve Cluster data
			if online_data_dict['pcluster'] is None:

				#Sometimes cluster information is not present. In the phagesdb database, it is is recorded as NULL.
				#When phages data is downloaded from phagesdb, NULL cluster data is converted to "Unclustered".
				#In these cases, leaving the cluster as NULL in phamerator won't work,
				#because NULL means Singleton. Therefore, assign the cluster as Unknown.
				row[3] = 'UNK'
			else:
				row[3] = online_data_dict['pcluster']['cluster']

			#Retrieve Subcluster data
			if online_data_dict['psubcluster'] is None:

				#Subcluster could be empty if by error no Cluster/Subcluster data
				#has yet been entered on phagesdb. But it may be empty because
				#there is no subcluster designation yet for members of the Cluster.
				row[8] = "none"

			else:
				row[8] = online_data_dict['psubcluster']['subcluster']

		except:
			write_out(output_file,"\nError: unable to retrieve Cluster and Subcluster data for phage %s from phagesdb." %row[1])
			row[3] = "none"
			row[8] = "none"
			table_errors += 1


	#Check Subcluster data
	if row[8] != "none":
		if row[8] not in phageSubcluster_set:
			print "The Subcluster %s is not currently in the database." % row[8]
			table_errors +=  question("\nError: %s is not the correct Subcluster for %s." % (row[8],row[1]))

		if len(row[8]) > 5:
			write_out(output_file,"\nError: phage %s Subcluster designation %s exceeds character limit." % (row[1],row[8]))
			table_errors += 1


	#Check Cluster data
	if row[3] != "none":
		if row[3].lower() == "singleton":
			row[3] = row[3].lower()

		if (row[3] not in phageCluster_set and row[3] != "singleton"):
			print "The Cluster %s is not currently in the database." % row[3]
			table_errors +=  question("\nError: %s is not the correct Cluster for %s." % (row[3],row[1]))

		if (row[3] != "singleton" and len(row[3]) > 5):
			write_out(output_file,"\nError: phage %s Cluster designation %s exceeds character limit." % (row[1],row[3]))
			table_errors += 1

		#If Singleton of Unknown Cluster, there should be no Subcluster
		if (row[3] == "singleton" or row[3] == "UNK"):
			if row[8] != "none":
				write_out(output_file,"\nError: phage %s Cluster and Subcluster discrepancy." % row[1])
				table_errors += 1

		#If not Singleton or Unknown or none, then Cluster should be part
		#of Subcluster data and the remainder should be a digit
		elif row[8] != "none":

			if (row[8][:len(row[3])] != row[3] or \
				row[8][len(row[3]):].isdigit() == False):

				write_out(output_file,"\nError: phage %s Cluster and Subcluster discrepancy." % row[1])
				table_errors += 1
		else:
			pass



	#Modify Status if needed
	if (row[4] not in phageStatus_set and row[4] != "none"):
			print "The status %s is not currently in the database." % row[4]
			table_errors +=  question("\nError: %s is not the correct status for %s." % (row[4],row[1]))
	if len(row[4]) > 5:
		write_out(output_file,"\nError: the status %s exceeds character limit." % row[4])
		table_errors += 1


	#Modify Description Qualifier if needed
	if (row[5] not in description_set and row[5] != "none"):
		print row[5] + " is an uncommon qualifier."
		table_errors += question("\nError: %s is an incorrect qualifier." % row[5])


	#Modify Accession if needed
	if row[7] == "retrieve":

		#On phagesdb, phages should always have a Genbank Accession field. If no accession, it will be ""
		try:
			row[7] = online_data_dict['genbank_accession']
			if row[7] != "":
				row[7] = row[7].strip() #Sometimes accession data from phagesdb have whitespace characters
				row[7] = row[7].split('.')[0] #Sometimes accession data from phagesdb have version suffix
			else:
				row[7] = "none"

		except:
			write_out(output_file,"\nError: unable to retrieve Accession data for phage %s from phagesdb." %row[1])
			row[7] = "none"
			table_errors += 1

	elif row[7].strip() == "":
		row[7] = "none"


	#Modify AnnotationAuthor
	#Author should only be 'hatfull','gbk', or 'none'.
	if row[9] == "hatfull":
		row[9] = "1"
	elif row[9] == "gbk":
		row[9] = "0"
	elif row[9] == "none":
		row[9] = "none"
	else:
		row[9] = "error"


	#Make sure run mode is permissible
	if row[10] not in set(run_mode_options_dict.keys()):
		write_out(output_file,"\nError: run mode is not valid for phage %s." %row[1])
		table_errors += 1
	elif row[10] == 'custom':
		run_mode_custom_total += 1
	else:
		pass




	#Rules for how each field is populated differs depending on each specific action



	#Update
	if row[0] == "update":

		#FirstPhageID
		if row[1] not in phageId_set:
			write_out(output_file,"\nError: %s is not a valid PhageID in the database." %row[1])
			table_errors += 1

		#Host, Cluster, Status
		if (row[2] == "none" or \
			row[3] == "none" or \
			row[4] == "none"):

			write_out(output_file,"\nError: %s does not have correctly populated HostStrain, Cluster, Subcluster, or Status fields." %row[1])
			table_errors += 1

		#Description
		if row[5] != "none":
			write_out(output_file,"\nError: %s does not have correctly populated Description field." %row[1])
			table_errors += 1

		#SecondPhageID
		if row[6] != "none":
			write_out(output_file,"\nError: %s should not have a genome listed to be removed." %row[1])
			table_errors += 1

		#Accession = it will either be an accession or it will be "none"
		#Subcluster = it will either be a Subcluster or it will be "none"

		#Author
		if row[9] != '1' and row[9] != '0':
			write_out(output_file,"\nError: %s does not have correctly populated Author field." %row[1])
			table_errors += 1

		#Run Mode
		if row[10] != 'none':
			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %row[1])
			table_errors += 1


	#Add
	elif row[0] == "add":

		#FirstPhageID
		if row[1] in phageId_set:
			write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %row[1])
			table_errors += 1

		#FirstPhageID, Host, Cluster, Status, Description
		if (row[1] == "none" or \
			row[2] == "none" or \
			row[3] == "none" or \
			row[4] == "none" or \
			row[5] == "none"):

			write_out(output_file,"\nError: %s does not have correctly populated fields." %row[1])
			table_errors += 1

		#Status
		if row[4] == "final":
			print row[1] + " to be added is listed as Final status, but no Draft (or other) genome is listed to be removed."
			table_errors +=  question("\nError: %s is not the correct status for %s." % (row[4],row[1]))

		#SecondPhageID
		if row[6] != "none":
			write_out(output_file,"\nError: %s to be added should not have a genome indicated for removal." %row[1])
			table_errors += 1

		#Accession = it will either be an accession or it will be "none"
		#Subcluster = it will either be a Subcluster or it will be "none"

		#Author
		if row[9] != '1' and row[9] != '0':
			write_out(output_file,"\nError: %s does not have correctly populated Author field." %row[1])
			table_errors += 1

		#Run Mode
		if row[10] == 'none':
			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %row[1])
			table_errors += 1


	#Remove
	elif row[0] == "remove":

		#FirstPhageID,Host, Cluster, Subcluster, Status, Description, Accession, Author, Run Mode
		if (row[1] != "none" or \
			row[2] != "none" or \
			row[3] != "none" or \
			row[4] != "none" or \
			row[5] != "none" or \
			row[7] != "none" or \
			row[8] != "none" or \
			row[9] != "none" or \
			row[10] != "none"):


			write_out(output_file,"\nError: %s to be removed does not have correctly populated fields." %row[6])
			table_errors += 1

		#SecondPhageID
		if row[6] not in phageId_set:
			write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %row[6])
			table_errors += 1


	#Replace
	elif row[0] == "replace":

		#FirstPhageID
		if row[1] == "none":
			write_out(output_file,"\nError: %s is not a valid PhageID. %s genome cannot be replaced." % (row[1],row[6]))
			table_errors += 1

		#FirstPhageID. If replacing a genome, ensure that if the genome to
		#be removed is not the same, that the new genome added has a unique name
		if (row[1] in phageId_set and row[1] != row[6]):
			write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %row[1])
			table_errors += 1

		#Host,Cluster,Status,Description
		if (row[2] == "none" or \
			row[3] == "none" or \
			row[4] == "none" or \
			row[5] == "none"):

			write_out(output_file,"\nError: %s does not have correctly populated fields." %row[1])
			table_errors += 1

		#SecondPhageID
		if row[6] not in phageId_set:
			write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %row[6])
			table_errors += 1

		if row[1] != row[6]:
			print "%s to replace %s is not spelled the same." %(row[1],row[6])
			table_errors +=  question("\nError: Phage %s is not spelled the same as phage %s." % (row[1],row[6]))

		#Accession = it will either be an accession or it will be "none"
		#Subcluster = it will either be a Subcluster or it will be "none"

		#Author
		if row[9] != '1' and row[9] != '0':
			write_out(output_file,"\nError: %s does not have correctly populated Author field." %row[1])
			table_errors += 1

		#Run Mode
		if row[10] == 'none':
			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %row[1])
			table_errors += 1



	else:
		pass

	genome_data_list.append(row)
	#genome_data_list elements are lists with follow index:
	#0 = Import action
	#1 = New PhageID
	#2 = Host
	#3 = Cluster
	#4 = Status
	#5 = Feature field
	#6 = PhageID to be removed
	#7 = Accession
	#8 = Subcluster
	#9 = Author
	#10 = Run mode

file_object.close()





#Now that all rows have been added to the list, verify there are no duplicate actions
add_set = set()
remove_set = set()
action_add_set = set()
action_remove_set = set()
action_add_remove_set = set()

#Create each set and do initial checks for duplications.
#If the Add name or Remove name is "none", skip that because there are
#probably duplicates of those.
for genome_data in genome_data_list:
	current_add = (genome_data[1],)
	current_remove = (genome_data[6],)
	current_action_add = (genome_data[0],genome_data[1])
	current_action_remove = (genome_data[0],genome_data[6])
	current_action_add_remove = (genome_data[0],genome_data[1],genome_data[6])


	#First check the one-field and two-field combinations
	if current_add[0] != "none":
		if current_add in add_set:
			print genome_data[1] + " appears to be involved in more than one step."
			table_errors += question("\nError: %s is duplicated" % str(current_add))

		else:
			add_set.add(current_add)

		if current_action_add in action_add_set:
			write_out(output_file,"\nError: %s is duplicated" % str(current_action_add))
			table_errors += 1
		else:
			action_add_set.add(current_action_add)

	if current_remove[0] != "none":
		if current_remove in remove_set:
			print genome_data[6] + " appears to be involved in more than one step."
			table_errors += question("\nError: %s is duplicated" % str(current_remove))

		else:
			remove_set.add(current_remove)

		if current_action_remove in action_remove_set:
			write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
			table_errors += 1
		else:
			action_remove_set.add(current_action_remove)

	#Now check the three-field combinations
	if current_action_add_remove in action_add_remove_set:
		write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
		table_errors += 1
	else:
		action_add_remove_set.add(current_action_add_remove)

#Once the sets are created, also check if genomes to be removed are
#found in the Add field and vice versa.
for genome_data in genome_data_list:
	current_add = (genome_data[1],)
	current_remove = (genome_data[6],)

	#If the phage name is not replacing itself, the Add name is not expected
	#to be in the Remove set and vice versa.
	if current_add != current_remove:
		if (current_add in remove_set and current_add != "none"):
			print genome_data[1] + " appears to be involved in more than one step."
			table_errors += question("\nError: %s is duplicated" % str(current_add))


		if (current_remove in add_set and current_remove != "none"):
			print genome_data[6] + " appears to be involved in more than one step."
			table_errors += question("\nError: %s is duplicated" % str(current_remove))




#Verify there are no duplicate accessions in the import table
importAccession_set = set()
for genome_data in genome_data_list:

	if genome_data[7] != "none":
		if genome_data[7] in importAccession_set:
			write_out(output_file,"\nError: Accession %s is duplicated in the import table." %genome_data[7])
			table_errors += 1
		else:
			importAccession_set.add(genome_data[7])






#Create separate lists of genome_data based on the indicated action: update, add/replace, remove
update_data_list = []
remove_data_list = []
add_replace_data_list = []
for genome_data in genome_data_list:
	if genome_data[0] == "update":
		update_data_list.append(genome_data)
	elif genome_data[0] == "remove":
		remove_data_list.append(genome_data)
	elif (genome_data[0] == "add" or genome_data[0] == "replace"):
		add_replace_data_list.append(genome_data)
	else:
		write_out(output_file,"\nError: during parsing of actions.")
		table_errors += 1


#Check to see if there are any inconsistencies with the
#update data compared to current phamerator data
for genome_data in update_data_list:

	#Initialize variable
	matched_phamerator_data = ''

	#Now that the Draft suffix is no longer appended to the import ticket,
	#this is less complications with matching to Phamerator PhageIDs
	try:
		matched_phamerator_data = phamerator_data_dict[genome_data[1]]
	except:
		write_out(output_file,"\nError: unable to retrieve phamerator data for %s." %genome_data[1])
		table_errors += 1
		continue

	#Host data check
	if genome_data[2] != matched_phamerator_data[2]:

		print "\n\nThere is conflicting host data for genome %s" % genome_data[1]
		print "Phamerator host: %s" % matched_phamerator_data[2]
		print "Import ticket host: %s" % genome_data[2]
		print "The new host data will be imported."
		table_errors += question("\nError: incorrect host data for %s." % genome_data[1])


	#Status data check
	if genome_data[4] != matched_phamerator_data[4]:

		#It is not common to change from 'gbk' or 'final' to anything else
		if matched_phamerator_data[4] != 'draft':

			print "\n\nThere is conflicting status data for genome %s" % genome_data[1]
			print "Phamerator status: %s" % matched_phamerator_data[4]
			print "Import ticket status: %s" % genome_data[4]
			print "The new status data will be imported."
			table_errors += question("\nError: incorrect status data for %s." % genome_data[1])

		#It is common to change status from 'draft' to 'final', but not anything else
		elif genome_data[4] != "final":

			print "\n\nThere is conflicting status data for genome %s" % genome_data[1]
			print "Phamerator status: %s" % matched_phamerator_data[4]
			print "Import ticket status: %s" % genome_data[4]
			print "The new status data will be imported."
			table_errors += question("\nError: incorrect status data for %s." % genome_data[1])


	#Accession data check
	if genome_data[7] == "none" and matched_phamerator_data[7] != "none":

		print "\n\nThere is conflicting accession data for genome %s" % genome_data[1]
		print "Phamerator accession: %s" % matched_phamerator_data[7]
		print "Import ticket accession: %s" % genome_data[7]
		print "The new accession data will be imported."
		table_errors += question("\nError: incorrect accession data for %s." % genome_data[1])

	elif genome_data[7] != "none" and matched_phamerator_data[7] != "none" and genome_data[7] != matched_phamerator_data[7]:

		print "\n\nThere is conflicting accession data for genome %s" % genome_data[1]
		print "Phamerator accession: %s" % matched_phamerator_data[7]
		print "Import ticket accession: %s" % genome_data[7]
		print "The new accession data will be imported."
		table_errors += question("\nError: incorrect accession data for %s." % genome_data[1])

	#Cluster, Subcluster check = no need to check this, as this data may be
	#more frequently updated than other fields.

	#Author
	#It is not common for authorship to change
	if genome_data[9] != matched_phamerator_data[9]:

		print "\n\nThere is conflicting author data for genome %s" % genome_data[1]
		print "Phamerator author: %s" % author_dictionary[matched_phamerator_data[9]]
		print "Import ticket author: %s" % author_dictionary[genome_data[9]]
		print "The new author data will be imported."
		table_errors += question("\nError: incorrect author data for %s." % genome_data[1])


#Check to see if genomes to be removed are the correct status
for genome_data in remove_data_list:

	#The next QC check relies on the PhageID being present in the Phamerator.
	#If it is not, this error will have already been identified, and a table_error
	#will already have been added. However, the code will crash without
	#here without adding an exception. So no need to add another table_error
	#in the except clause.
	try:
		matched_phamerator_data = phamerator_data_dict[genome_data[6]]
		if matched_phamerator_data[4] != "draft":
			print "The genome %s to be removed is currently %s status." % \
					(genome_data[6],matched_phamerator_data[4])
			table_errors += question("\nError: %s is %s status and should not be removed." % \
					(genome_data[6],matched_phamerator_data[4]))
	except:
		pass


#If no errors encountered, print list of action items to be
#implemented and continue. Otherwise, exit the script.
if table_errors == 0:
	write_out(output_file,"\nImport table verified with 0 errors.")
	write_out(output_file,"\nList of all actions to be implemented:")
	for genome_data in genome_data_list:
		write_out(output_file,"\n" + str(genome_data))
	raw_input("\nPress ENTER to proceed to next import stage.")

else:
	write_out(output_file,"\n%s error(s) encountered with import file.\nNo changes have been made to the database." % table_errors)
	write_out(output_file,"\nExiting import script.")
	output_file.close()
	sys.exit(1)
