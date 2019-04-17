




#Define several functions

#Print out statements to both the terminal and to the output file
#For SQL statements that may be long (>150 characters), don't print entire statements.
def write_out(filename,statement):
	if (statement[:7] == "\nINSERT" or statement[:7] == "\nUPDATE" or statement[:7] == "\nDELETE"):
		if len(statement) > 150:
			print statement[:150] + "...(statement truncated)"
			filename.write(statement[:150] + "...(statement truncated)")
		else:
			print statement
			filename.write(statement)
	else:
		print statement
		filename.write(statement)


#For questionable data, user is requested to clarify if the data is correct or not
def question(message):
	number = -1
	while number < 0:
		value = raw_input("Is this correct? (yes or no): ")
		if (value.lower() == "yes" or value.lower() == "y"):
			number = 0
		elif (value.lower() == "no" or value.lower() == "n"):
			write_out(output_file,message)
			number = 1
		else:
			print "Invalid response."
	#This value will be added to the current error total. If 0, no error was encountered. If 1, an error was encountered.
	return number

#Exits MySQL
def mdb_exit(message):
	write_out(output_file,"\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
	write_out(output_file,message)
	write_out(output_file,"\nThe import script did not complete.")
	write_out(output_file,"\nExiting MySQL.")
	cur.execute("ROLLBACK")
	cur.execute("SET autocommit = 1")
	cur.close()
	con.close()
	write_out(output_file,"\nExiting import script.")
	output_file.close()
	sys.exit(1)




#This function decides whether Cluster2 or Subcluster2 data
#gets assigned to the Cluster field
def assign_cluster_field(subcluster,cluster):

	cluster_field_data = ""
	if subcluster != "none":
		cluster_field_data = subcluster
	else:
		cluster_field_data = cluster

	return cluster_field_data


#If phage Cluster is Singleton, make sure MySQL statement is created correctly
def create_cluster_statement(phage_name,cluster):
	cluster_statement = ""
	if cluster == "singleton":
		cluster_statement = "UPDATE phage SET Cluster = NULL" + " WHERE PhageID = '" + phage_name + "';"
	else:
		cluster_statement = "UPDATE phage SET Cluster = '" + cluster + "' WHERE PhageID = '" + phage_name + "';"
	return cluster_statement




#If phage Cluster is Singleton, make sure MySQL statement is created correctly
def create_cluster2_statement(phage_name,cluster):
	cluster2_statement = ""
	if cluster == "singleton":
		cluster2_statement = "UPDATE phage SET Cluster2 = NULL" + " WHERE PhageID = '" + phage_name + "';"
	else:
		cluster2_statement = "UPDATE phage SET Cluster2 = '" + cluster + "' WHERE PhageID = '" + phage_name + "';"
	return cluster2_statement


#If phage Subcluster is empty ("none"), make sure MySQL statement is created correctly
def create_subcluster2_statement(phage_name,subcluster):
	subcluster2_statement = ""
	if subcluster == "none":
		subcluster2_statement = "UPDATE phage SET Subcluster2 = NULL" + " WHERE PhageID = '" + phage_name + "';"
	else:
		subcluster2_statement = "UPDATE phage SET Subcluster2 = '" + subcluster + "' WHERE PhageID = '" + phage_name + "';"
	return subcluster2_statement








#Function to split gene description field
def retrieve_description(genbank_feature,description_field):
	description = genbank_feature.qualifiers[description_field][0].lower().strip()
	split_description = description.split(' ')
	if description == "hypothetical protein":
		description = ""

	elif description == "phage protein":
		description = ""

	elif description == "unknown":
		description = ""

	elif description == "conserved hypothetical protein":
		description = ""

	elif description.isdigit():
		description = ""

	elif len(split_description) == 1:

		if (split_description[0][:2] == "gp" and split_description[0][2:].isdigit()):
			description = ""

		elif (split_description[0][:3] == "orf" and split_description[0][3:].isdigit()):
			description = ""

		else:
			description = genbank_feature.qualifiers[description_field][0].strip()

	elif len(split_description) == 2:

		if (split_description[0] == "orf" and split_description[1].isdigit()):
			description = ""

		elif (split_description[0] == "putative" and split_description[1][:7] == "protein"):
			description = ""

		else:
			description = genbank_feature.qualifiers[description_field][0].strip()

	else:
		description = genbank_feature.qualifiers[description_field][0].strip()
	return description



#Function to search through a list of elements using a regular expression
def find_name(expression,list_of_items):
	search_tally = 0
	for element in list_of_items:
		search_result = expression.search(element)
		if search_result:
			search_tally += 1
	return search_tally




def change_descriptions():
   print "These will be ignored, unless this is NOT correct."
   print "If it is NOT correct, no error will be generated."
   print "Instead, only gene descriptions in this field will be retained."




def check_tRNA_product(product_field):

	#This is an initial attempt at checking the tRNA product description
	#Ultimately, a regular expression would be better to use
	#tRNA product example = 'tRNA-Ser (AGC)'

	#The code block below functions, but it does not fully account for
	#tRNA-OTHER descriptions, tRNA-Stop descriptions,
	#and it does not check the accuracy of
	#the amino acid and anticodon pairing.
	#The biggest problem is that the expected product and note descriptions
	#are expected to change after they reach NCBI, so it is not clear
	#how to best address that issue here, since nothing in the import
	#table reflects WHERE the annotated genome came from.



	product_error = 0


	#product starts off as lowercase 'trna-ser (agc)'
	#split1_list = 'trna' and 'ser (agc)'
	tRNA_product_split1_list = product_field.split('-')

	#If product is missing, an error will already have been thrown.
	#The product should have a hypthen, so only parse product if it can be
	#split into two elements.
	if len(tRNA_product_split1_list) == 2:

		tRNA_product_split1_prefix = tRNA_product_split1_list[0].strip() #'trna'

		#split2_list = 'ser' and 'agc)'
		tRNA_product_split2_list = tRNA_product_split1_list[1].split('(')
		tRNA_product_amino_acid_three = tRNA_product_split2_list[0].strip() #'ser'

		if tRNA_product_amino_acid_three != 'other' and \
			tRNA_product_amino_acid_three != 'stop' and \
			len(tRNA_product_amino_acid_three) != 3:
				product_error += 1

		#The code block below checks for the presence of an anticodon.
		#No need to use it currently, since there is so much variability
		#at the tRNA product field, but retain for future use.
		# if len(tRNA_product_split2_list) == 2:
		#
		#     #split3_list = 'agc' and ''
		#     tRNA_product_split3_list = tRNA_product_split2_list[1].split(')')
		#
		#     #Only check the anticodon if the amino acid is NOT 'other'
		#     if tRNA_product_amino_acid_three != 'other' and \
		#         len(tRNA_product_split3_list) == 2:
		#
		#         tRNA_product_anticodon = tRNA_product_split3_list[0].strip() #'agc'
		#         if len(tRNA_product_anticodon) != 3:
		#             product_error += 1
		#     else:
		#         product_error += 1
		#
		# else:
		#     product_error += 1


	else:
		product_error += 1

	return product_error


#Allows user to select specific options
def select_option(message,valid_response_set):

	response_valid = False
	while response_valid == False:
		response = raw_input(message)
		if response.isdigit():
			response = int(response)
		else:
			response = response.lower()

		if response in valid_response_set:
			response_valid = True
			if response == 'y':
				response  = 'yes'
			elif response == 'n':
				response  = 'no'
		else:
			print 'Invalid response.'
	return response











#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"







#Set up run type
run_type_options = [\
	'none',\
	'test',\
	'production']
print '\n\nThe following run types are available:\n'
#print '0: ' + run_type_options[0]
# print '1: ' + run_type_options[1]
# print '2: ' + run_type_options[2]
print '1: ' + run_type_options[1] + ' (checks flat files for accuracy, but the database is not changed.)'
print '2: ' + run_type_options[2] + ' (after testing files, the database is updated.)'
run_type = select_option(\
	"\nWhich run type do you want? ", \
	set([1,2]))
run_type = run_type_options[run_type]

























###Pasted misc functions from compare_databases.py


#Exits MySQL
def mdb_exit(message):
    print "\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`
    print "\nThe import script did not complete."
    print "\nExiting MySQL."
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print "\nExiting import script."
    sys.exit(1)

#Closes all file handles currently open
def close_all_files(file_list):
    for file_handle in file_list:
        file_handle.close()

#Make sure there is no "_Draft" suffix
def remove_draft_suffix(value):
    # Is the word "_Draft" appended to the end of the name?
    value_truncated = value.lower()
    if value_truncated[-6:] == "_draft":
        value_truncated = value_truncated[:-6]
    return value_truncated

def parse_strand(value):
    value = value.lower()
    if value == "f" or value == "forward":
        value = "forward"
    elif value == "r" or value == "reverse":
        value = "reverse"
    else:
        value = "NA"
    return value


#Function to split gene description field
def retrieve_description(description):
    if description is None:
        description = ''
    else:
        description = description.lower().strip()

    split_description = description.split(' ')
    if description == 'hypothetical protein':
        search_description = ''

    elif description == 'phage protein':
        search_description = ''

    elif description == 'unknown':
        search_description = ''

    elif description == '\\n':
        search_description = ''

    elif description.isdigit():
        search_description = ''

    elif len(split_description) == 1:

        if split_description[0][:2] == 'gp' and split_description[0][2:].isdigit():
            search_description = ''

        elif split_description[0][:3] == 'orf' and split_description[0][3:].isdigit():
            search_description = ''

        else:
            search_description = description

    elif len(split_description) == 2:

        if split_description[0] == 'gp' and split_description[1].isdigit():
            search_description = ''

        elif split_description[0] == 'orf' and split_description[1].isdigit():
            search_description = ''

        else:
            search_description = description

    else:
        search_description = description
    return description,search_description



#Function to search through a list of elements using a regular expression
def find_name(expression,list_of_items):
    search_tally = 0
    for element in list_of_items:
        search_result = expression.search(element)
        if search_result:
            search_tally += 1
    return search_tally


#Allows user to select specific options
def select_option(message,valid_response_set):

    response_valid = False
    while response_valid == False:
        response = raw_input(message)
        if response.isdigit():
            response = int(response)
        else:
            response = response.lower()

        if response in valid_response_set:
            response_valid = True
            if response == 'y':
                response  = 'yes'
            elif response == 'n':
                response  = 'no'
        else:
            print 'Invalid response.'
    return response



#Output list to file
def output_to_file(data_list,filename,genome_status_selected,database_string,genome_author_selected):
    filename_fh = open(os.path.join(main_output_path,date + "_" + filename), 'w')
    filename_writer = csv.writer(filename_fh)
    filename_writer.writerow([date + ' Database comparison'])
    filename_writer.writerow([database_string])
    filename_writer.writerow([genome_author_selected])
    filename_writer.writerow([genome_status_selected])
    for element in data_list:
        filename_writer.writerow([element])
    filename_fh.close()





#Ensure the output filename is unique
def create_unique_filename(filename_directory,filename_base,filename_ext):

    file_exists = True
    rename_counter = 0
    while file_exists == True:

        if rename_counter == 0:
            unique_filename = filename_base + filename_ext
        else:
            unique_filename = filename_base + '_' + str(rename_counter) + filename_ext

        unique_filename_path = os.path.join(filename_directory,unique_filename)
        file_exists = os.path.isfile(unique_filename_path)

        if file_exists == True:
            print 'Warning: duplicate output file:'
            print unique_filename_path
            print 'Filename will be modified.'
            raw_input('Press ENTER to proceed')
            rename_counter += 1

    return unique_filename_path










######

#New functions for new import script




def parse_import_tickets(list_of_data):

	list_of_tickets = []
	for input_row in list_of_data:

		#Verify the row of information has the correct number of fields to parse.
		if len(input_row) != 11:

			#TODO error handling
			write_out(output_file,"\nRow in import table is not formatted correctly: " + str(input_row))
			table_errors += 1
			continue


		ticket = ImportTicket()
		ticket.type = input_row[0] #Import action
		ticket.primary_phage_id = input_row[1] #New PhageID
		ticket.host = input_row[2] #Host
		ticket.cluster = input_row[3] #Cluster
		ticket.subcluster = input_row[4] #Subcluster
		ticket.status = input_row[5] #Status
		ticket.description_field = input_row[7] #Feature field
		ticket.accession = input_row[8] #Accession
		ticket.annotation_author = input_row[6] #AnnotationAuthor
		ticket.secondary_phage_id = input_row[10] #PhageID to be removed
		ticket.run_mode = input_row[9] #Run mode

		list_of_tickets.append(ticket)

		return(list_of_tickets)








def parse_phagesdb_data(phagesdb_genome,data_dict):

	#Name and Search Name
	phagesdb_genome.set_phage_name(online_data_dict['phage_name'])


    #Host, Accession
    phagesdb_genome.set_host(online_data_dict['isolation_host']['genus'])
    phagesdb_genome.set_accession(online_data_dict['genbank_accession'])


    #Cluster
	# On phagesdb, phages may have a Cluster and no Subcluster info
    # (which is set to None). If the phage has a Subcluster,
    # it should also have a Cluster. If by accident no Cluster or
    # Subcluster info is added at the time the genome is added to
    # phagesdb, the Cluster may automatically be set to NULL,
	# which gets converted to "Unclustered". This is problematic
	# because in Phamerator NULL means Singleton, and the long
	# form "Unclustered" will be filtered out later in the script
	# due to its character length, so it needs to be abbreviated.
	if online_data_dict['pcluster'] is None:
		phagesdb_genome.cluster = 'UNK'
	elif online_data_dict['pcluster']['cluster'].lower() == "singleton"
		phagesdb_genome.cluster = online_data_dict['pcluster']['cluster'].lower()
	else:
		phagesdb_genome.cluster = online_data_dict['pcluster']['cluster']





    #Subcluster
	# Subcluster could be empty if by error no Cluster or
    # Subcluster data has yet been entered on phagesdb. But
    # it may be empty because there is no subcluster
    # designation yet for members of the Cluster.
	if online_data_dict['psubcluster'] is None:
		phagesdb_genome.subcluster = "none"
	else:
		phagesdb_genome.subcluster = online_data_dict['psubcluster']['subcluster']


    #Check to see if there is a fasta file stored on phagesdb for this phage
    if online_data_dict['fasta_file'] is not None:
        fastafile_url = online_data_dict['fasta_file']

        response = urllib2.urlopen(fastafile_url)

        retrieved_fasta_file = response.read()
        response.close()

		#TODO convert to Biopython object?
        #All sequence rows in the fasta file may not have equal widths, so some processing of the data is required
        #If you split by newline, the header is retained in the first list element
        split_fasta_data = retrieved_fasta_file.split('\n')
        pdb_sequence = ''
        index = 1
        while index < len(split_fasta_data):

			#Strip off potential whitespace before appending, such as '\r'
            pdb_sequence = pdb_sequence + split_fasta_data[index].strip()
            index += 1
        phagesdb_genome.sequence = pdb_sequence

		#TODO do I need to return this object?
		return(phagesdb_genome)







def validate_ticket(ticket):

	ticket.check_retrieve_status()
	ticket.check_type()
	ticket.check_host()
	ticket.check_cluster_subcluster()
	ticket.check_status()
	ticket.check_description_field()
	ticket.check_accession()
	ticket.check_annotation_author()
	ticket.check_run_mode()
	ticket.check_ticket()

	#TODO do I need to return this object?
	return(ticket)






#TODO need to complete
# Verify there are no duplicate actions

def validate_tickets(list_of_tickets):


    # Initialize all calculated attributes:
    add_set = set()
    remove_set = set()
    action_add_set = set()
    action_remove_set = set()
    action_add_remove_set = set()
    import_accession_set = set()



    #Create each set and do initial checks for duplications.
    #If the Add name or Remove name is "none", skip that because there are
    #probably duplicates of those.
    for ticket in list_of_tickets:
    	current_add = (ticket.primary_phage_id,)
    	current_remove = (ticket.secondary_phage_id,)
    	current_action_add = (ticket.type,ticket.primary_phage_id)
    	current_action_remove = (ticket.type,ticket.secondary_phage_id)
    	current_action_add_remove = (ticket.type,ticket.primary_phage_id,ticket.secondary_phage_id)


    	#First check the one-field and two-field combinations
    	if current_add[0] != "none":
    		if current_add in add_set:
    			print ticket.primary_phage_id + " appears to be involved in more than one step."

                #TODO error handling
                table_errors += question("\nError: %s is duplicated" % str(current_add))

    		else:
    			add_set.add(current_add)

    		if current_action_add in action_add_set:

                #TODO error handling
    			write_out(output_file,"\nError: %s is duplicated" % str(current_action_add))
    			table_errors += 1

    		else:
    			action_add_set.add(current_action_add)

    	if current_remove[0] != "none":

    		if current_remove in remove_set:

                #TODO error handling
    			print ticket.secondary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_remove))

    		else:
    			remove_set.add(current_remove)

    		if current_action_remove in action_remove_set:

                #TODO error handling
    			write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
    			table_errors += 1

    		else:
    			action_remove_set.add(current_action_remove)

    	#Now check the three-field combinations
    	if current_action_add_remove in action_add_remove_set:

            #TODO error handling
    		write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
    		table_errors += 1

    	else:
    		action_add_remove_set.add(current_action_add_remove)


    #Once the sets are created, also check if genomes to be removed are
    #found in the Add field and vice versa.
    for ticket in list_of_tickets:
    	current_add = (ticket.primary_phage_id,)
    	current_remove = (ticket.secondary_phage_id,)

    	#If the phage name is not replacing itself, the Add name is not expected
    	#to be in the Remove set and vice versa.
    	if current_add != current_remove:
    		if (current_add in remove_set and current_add != "none"):

                #TODO error handling
    			print ticket.primary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_add))


    		if (current_remove in add_set and current_remove != "none"):

                #TODO error handling
    			print ticket.secondary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_remove))




    #Verify there are no duplicate accessions in the import table
    for ticket in list_of_tickets:

    	if ticket.accession != "none":
    		if ticket.accession in import_accession_set:

                #TODO error handling
    			write_out(output_file,"\nError: Accession %s is duplicated in the import table." %ticket.accession)
    			table_errors += 1
    		else:
    			import_accession_set.add(ticket.accession)




	#TODO return proper information
	return(pass)






def parse_phamerator_data(phamerator_genome,data_tuple):


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

	phamerator_genome.set_phage_id(genome_tuple[0])
	phamerator_genome.set_phage_name(genome_tuple[1])
	phamerator_genome.set_host(genome_tuple[2])
	phamerator_genome.set_sequence(genome_tuple[3])
	phamerator_genome.set_status(genome_tuple[4])

	#TODO singletons handled properly?
	phamerator_genome.set_cluster(genome_tuple[5])

	#TODO non-subclustered handled properly
	phamerator_genome.set_subcluster(genome_tuple[8])
	phamerator_genome.set_date_last_modified_filled(genome_tuple[6])
	phamerator_genome.set_accession_filled(genome_tuple[7])
	phamerator_genome.annotation_author = str(genome_tuple[9])
	phamerator_genome.annotation_qc = str(genome_tuple[10])
	phamerator_genome.retrieve_record = str(genome_tuple[11])


	return(phamerator_genome)

























#Count # of ticket types
#TODO not sure when I should re-implement this
if row[0] in action_set:
	if row[0] == "add":
		add_total += 1
	elif row[0] == "remove":
		remove_total += 1
	elif row[0] == "replace":
		replace_total += 1
	elif row[0] == "update":
		update_total += 1
	else:
		pass










###
