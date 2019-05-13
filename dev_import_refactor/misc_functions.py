




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
    if value[-6:].lower() == "_draft":
        value = value[:-6]
    return(value)

def parse_strand(value):

    value = value.lower()
    if value == "f" or value == "forward":
        value = "forward"
    elif value == "r" or value == "reverse":
        value = "reverse"
    else:
        value = "NA"
    return(value)



def parse_strand_for_import(value):

    #Orientation
    if value == 1:
		value = "forward"
	elif value == -1:
		value = "reverse"
	else:
	    #Some features have no orientation
		value = "NA"

	return(value)











#Function to split gene description field
def retrieve_description(description):




    if description is None:
        description = ''

    description = description.strip()
    processed_description = description.lower()
    split_description = processed_description.split(' ')

    if processed_description == 'hypothetical protein':
        processed_description = ''

    elif processed_description == 'phage protein':
        processed_description = ''

    elif processed_description == 'unknown':
        processed_description = ''

    elif processed_description == 'conserved hypothetical protein':
        processed_description = ''

    elif processed_description == '\\n':
        processed_description = ''

    elif processed_description == '.':
        processed_description = ''

    elif processed_description.isdigit():
        processed_description = ''

    elif len(split_description) == 1:

        if split_description[0][:2] == 'gp' and split_description[0][2:].isdigit():
            processed_description = ''

        elif split_description[0][:3] == 'orf' and split_description[0][3:].isdigit():
            processed_description = ''

        else:
            processed_description = description

    elif len(split_description) == 2:

        if split_description[0] == 'gp' and split_description[1].isdigit():
            processed_description = ''

        elif split_description[0] == 'orf' and split_description[1].isdigit():
            processed_description = ''

        elif (split_description[0] == "putative" and split_description[1][:7] == "protein"):
            processed_description = ""

        else:
            processed_description = description

    else:
        processed_description = description


    return(description,processed_description)

















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










match_phamerator_to_tickets(list_of_matched_objects, all_phamerator_data):

	for matched_data_obj in list_of_matched_objects:
		phage_id = matched_data_obj.ticket.primary_phage_id

		try:
			matched_genome = phamerator_genome_dict[phage_id]

			#Now add the Phamerator data to the MatchedGenomes object
			matched_data_obj.matched_genomes_dict["phamerator"] = matched_genome

		except:
			#TODO error handling

	return list_of_matched_objects






def assign_match_strategy(list_of_matched_objects):

	for matched_object in list_of_matched_objects:
		ticket = matched_object.ticket


		#TODO complete this function
	pass



def match_flat_files_to_tickets(list_of_matched_objects, all_flat_file_data):

	flat_file_genome_dict = {}



	# TODO:
	# Iterate through all tickets and determine what the strategy is to
	# match tickets to flat files.
	strategy = assign_match_strategy(list_of_matched_objects)

	#TODO need to set strategy variable in advance
	# Now that the flat file data has been retrieved and parsed,
	# match them to ticket data

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


	return list_of_matched_objects





def create_matched_object_dict(list_of_matched_objects):

	dictionary = {}
	list_of_update_objects = []
	list_of_remove_objects = []
	list_of_add_replace_objects = []


	for matched_object in list_of_matched_objects:
		type = matched_object.ticket.ticket_type

		if type == "update":
			list_of_update_objects.append(matched_object)

		elif type == "remove":
			list_of_remove_objects.append(matched_object)

		elif (type = "add" or type == "replace"):
			list_of_add_replace_objects.append(matched_object)

		else:
			pass
			# TODO if ticket type is none of the above, thrown an error?

	dictionary["update"] = list_of_update_objects
	dictionary["remove"] = list_of_remove_objects
	dictionary["add_replace"] = list_of_add_replace_objects

	return dictionary










###copied below from compare_databases script



#
# #Count # of ticket types
# #TODO not sure when I should re-implement this
# if row[0] in action_set:
# 	if row[0] == "add":
# 		add_total += 1
# 	elif row[0] == "remove":
# 		remove_total += 1
# 	elif row[0] == "replace":
# 		replace_total += 1
# 	elif row[0] == "update":
# 		update_total += 1
# 	else:
# 		pass
#









###
