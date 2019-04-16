




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
