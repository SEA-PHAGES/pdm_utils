"""Misc. base/simple functions. These should not require import of other
k_phamerate modules to prevent circular imports."""

#Note: used to be 'find_name' function.
def find_expression(expression,list_of_items):
    """Searches through a list of items and counts the number of items
    that match the expression.
    """

    search_tally = 0
    for element in list_of_items:
        search_result = expression.search(element)
        if search_result:
            search_tally += 1
    return search_tally



def remove_draft_suffix(value):
    """Removes the '_Draft' suffix if present."""

    if value[-6:].lower() == "_draft":
        value = value[:-6]
    return(value)



#Old 'parse_strand' function and 'parse_strand_for_import' function combined.
# Note: parse_strand_for_import used to convert numeric to long string format.
# E.g. 1 == forward
def reformat_strand(input_value, format, case = False):
    """Converts common strand orientation formats, including:
    'fr_long' ('forward', 'reverse')
    'fr_short' ('f', 'r')
    'tb_long' ('top', 'bottom')
    'tb_short' ('t', 'b')
    'wc_long' ('watson', 'crick')
    'wc_short' ('w','c')
    'operator' ('+', '-')
    'numeric' (1, -1)."""

    format_dict = {"fr_long":["forward", "reverse"], \
                    "fr_short":["f", "r"], \
                    "tb_long":["top", "bottom"], \
                    "tb_short":["t", "b"], \
                    "wc_long":["watson", "crick"], \
                    "wc_short":["w", "c"], \
                    "operator":["+", "-"], \
                    "numeric":[1, -1]}

    strand1_values = set()
    strand2_values = set()
    for values in format_dict.values():

        strand1_values.add(values[0])
        strand2_values.add(values[1])


    if format in format_dict.keys():

        if isinstance(input_value, str):
            new_value = input_value.lower()
        else:
            new_value = input_value

        if new_value in strand1_values:
            new_value = format_dict[format][0]

        elif new_value in strand2_values:
            new_value = format_dict[format][1]

        else:
            new_value = "NA"

    else:
        new_value = "NA"


    if isinstance(new_value, str):
        if case == True:
            new_value = new_value.capitalize()

    return new_value




#Note: used to be 'retrieve_description' function.
def reformat_description(description):
    """Removes leading and trailing whitespace from text and returns this value.
    Then it assesses whether the description is informative or not, processes
    it if it is too generic (converting it to an empty value), then returning
    this second value.
    """

    if description is None:
        description = ""

    description = description.strip()
    processed_description = description.lower()
    split_description = processed_description.split(" ")

    if processed_description == "hypothetical protein":
        processed_description = ""
    elif processed_description == "phage protein":
        processed_description = ""
    elif processed_description == "unknown":
        processed_description = ""
    elif processed_description == "conserved hypothetical protein":
        processed_description = ""
    elif processed_description == "\\n":
        processed_description = ""
    elif processed_description == ".":
        processed_description = ""
    elif processed_description.isdigit():
        processed_description = ""

    elif len(split_description) == 1:

        if split_description[0][:2] == "gp" and \
            split_description[0][2:].isdigit():
            processed_description = ""
        elif split_description[0][:3] == "orf" and \
            split_description[0][3:].isdigit():
            processed_description = ""
        else:
            processed_description = description

    elif len(split_description) == 2:

        if split_description[0] == "gp" and \
            split_description[1].isdigit():
            processed_description = ""
        elif split_description[0] == "orf" and \
            split_description[1].isdigit():
            processed_description = ""
        elif (split_description[0] == "putative" and \
            split_description[1][:7] == "protein"):
            processed_description = ""
        else:
            processed_description = description

    else:
        processed_description = description

    return(description,processed_description)


def parse_fasta_file(fasta_file):
    """Parses sequence data from a fasta-formatted file.
    """

    # TODO convert to Biopython object?
    # All sequence rows in the fasta file may not have equal widths,
    # so some processing of the data is required. If you split by newline,
    # the header is retained in the first list element.
    split_fasta_data = fasta_file.split('\n')

    sequence = ""
    index = 1
    while index < len(split_fasta_data):

        # Strip off potential whitespace before appending, such as '\r'.
        sequence = sequence + split_fasta_data[index].strip()
        index += 1

    return sequence




def identify_unique_items(complete_list):
    """Iterate over a list of items and generate two lists.
    The first list contains items that are unique and non-duplicated
    in the original list.
    The second list contains items that are not unique in the original
    list, although non-duplicated items are returned."""

    unique_set = set() # Set of all non-duplicated unique items.
    duplicate_set = set() # Set of all duplicated items.

    for item in complete_list:
        if item not in unique_set:
            unique_set.add(item)
        else:
            duplicate_set.add(item)

    # Remove items from the unique list that eventually were duplicated.
    unique_set = unique_set - duplicate_set

    return (list(unique_set), list(duplicate_set))




def trim_characters(string):
    """Remove leading and trailing generic characters from a string.
    """

    generic_characters = set([".", ",", ";", "-"])

    if len(string) > 0:

        # Trim non-specific trailing characters.
        value1 = True
        while value1:
            if string[-1] in generic_characters:
                string = string[:-1]
            else:
                value1 = False

        # Trim non-specific leading characters.
        value2 = True
        while value2:
            if string[0] in generic_characters:
                string = string[1:]
            else:
                value2 = False

    return string



# TODO this can probably be improved.
def parse_names_from_record_field(description):
    """Parse string of text from GenBank-formatted flat file to
    identify the phage name and host name.
    """

    generic_words = set(["complete", "genome", "phage"])

    host_name = ""
    phage_name = ""

    # Remove trailing whitespace and split into a list.
    description = description.strip()
    split_description = description.split()


    # Trim leading and trailing generic characters.
    index = 0
    while index < len(split_description):
        split_description[index] = trim_characters(split_description[index])
        index += 1


    # Iterate through the list of processed words and attempt to
    # identify the host name and phage name.
    index = 0
    while index < len(split_description):

        word = split_description[index]

        word_lower = word.lower()

        # Attempt to identify the host.

        # Sometimes the host name is the word preceding 'phage' or 'virus'.
        # e.g. 'Mycobacterium phage'.
        if index > 0:
            if (word_lower == "phage" or word_lower == "virus"):
                host_name = split_description[index - 1]

        # Sometimes the host name is merged with 'phage'.
        # e.g. 'Mycobacteriophage'.
        elif (len(word) > 5 and word_lower[-5:] == "phage"):
            host_name = word[:-5]

        else:
            pass


        # Attempt to identify the phage.

        # If there is one non-generic word in the string, assign it as the
        # phage name.
        if len(split_description) == 1:

            if (len(word) > 5 and word_lower[-5:] == "phage"):
                host_name = word[:-5]

            elif word_lower not in generic_words:
                phage_name = word

            else:
                pass

        # Sometimes the phage name is the word following 'phage' or 'virus'.
        # e.g. 'Mycobacterium phage Trixie' or 'Mycobacteriophage Trixie'
        if index < (len(split_description) - 1):
            if (word_lower[-5:] == "phage" or \
                word_lower[-5:] == "virus"):

                phage_name = split_description[index + 1]

        index += 1


    return (phage_name, host_name)

































#TODO Unit test below


# TODO finish revamping code for matching features.
# TODO unit test.
def analyze_sets(set1, set2):
    """Compute the intersection and differences between two sets."""

    set_intersection = set1 & set2
    set1_diff = set1 - set2
    set2_diff = set2 - set1

    return (set_intersection, set1_diff, set2_diff)

# TODO finish revamping code for matching features.
# TODO unit test.
def match_items(self, list1, list2):
    """Match values of two lists. Return the matched value list,
    and a list of unmatched values from each original list."""

    # Determine the unique values in each list.
    list1_items_unique, list1_items_duplicate = \
        identify_unique_items(list1_items)

    list2_items_unique, list2_items_duplicate = \
        identify_unique_items(list2_items)

    # Create matched and difference sets.
    items_matched, list1_items_unmatched, list2_items_unmatched = \
        analyze_sets( \
            list1_items_unique, list2_items_unique)


    items_matched = list(items_matched)
    list1_items_unmatched = list(list1_items_unmatched)
    list2_items_unmatched = list(list2_items_unmatched)

    return (items_matched, list1_items_unmatched, list2_items_unmatched)









#TODO set up and unit test
# #Print out statements to both the terminal and to the output file
# #For SQL statements that may be long (>150 characters), don't print
# # entire statements.
# def write_out(filename,statement):
# 	if (statement[:7] == "\nINSERT" or statement[:7] == "\nUPDATE" or statement[:7] == "\nDELETE"):
# 		if len(statement) > 150:
# 			print statement[:150] + "...(statement truncated)"
# 			filename.write(statement[:150] + "...(statement truncated)")
# 		else:
# 			print statement
# 			filename.write(statement)
# 	else:
# 		print statement
# 		filename.write(statement)





#TODO set up and unit test
# #For questionable data, user is requested to clarify if the data is correct or not
# def question(message):
# 	number = -1
# 	while number < 0:
# 		value = raw_input("Is this correct? (yes or no): ")
# 		if (value.lower() == "yes" or value.lower() == "y"):
# 			number = 0
# 		elif (value.lower() == "no" or value.lower() == "n"):
# 			write_out(output_file,message)
# 			number = 1
# 		else:
# 			print "Invalid response."
# 	#This value will be added to the current error total. If 0, no error was encountered. If 1, an error was encountered.
# 	return number



#TODO set up and unit test
# def change_descriptions():
#
#
#     print "These will be ignored, unless this is NOT correct."
#     print "If it is NOT correct, no error will be generated."
#     print "Instead, only gene descriptions in this field will be retained."



#TODO set up and unit test
# #Allows user to select specific options
# def select_option(message,valid_response_set):
#
# 	response_valid = False
# 	while response_valid == False:
# 		response = raw_input(message)
# 		if response.isdigit():
# 			response = int(response)
# 		else:
# 			response = response.lower()
#
# 		if response in valid_response_set:
# 			response_valid = True
# 			if response == 'y':
# 				response  = 'yes'
# 			elif response == 'n':
# 				response  = 'no'
# 		else:
# 			print 'Invalid response.'
# 	return response
#





#TODO unit test - this may no longer be needed.
# #Closes all file handles currently open
# def close_all_files(file_list):
#     for file_handle in file_list:
#         file_handle.close()


#TODO unit test - this may no longer be needed
# def choose_run_type():
#
#     run_type_options = [\
#     	'none',\
#     	'test',\
#     	'production']
#     print '\n\nThe following run types are available:\n'
#     #print '0: ' + run_type_options[0]
#     # print '1: ' + run_type_options[1]
#     # print '2: ' + run_type_options[2]
#     print '1: ' + run_type_options[1] + ' (checks flat files for accuracy, but the database is not changed.)'
#     print '2: ' + run_type_options[2] + ' (after testing files, the database is updated.)'
#     run_type = select_option(\
#     	"\nWhich run type do you want? ", \
#     	set([1,2]))
#     run_type = run_type_options[run_type]



#TODO unit test - this may no longer be needed
# #Output list to file
# def output_to_file(data_list,filename,genome_status_selected,database_string,genome_author_selected):
#     filename_fh = open(os.path.join(main_output_path,date + "_" + filename), 'w')
#     filename_writer = csv.writer(filename_fh)
#     filename_writer.writerow([date + ' Database comparison'])
#     filename_writer.writerow([database_string])
#     filename_writer.writerow([genome_author_selected])
#     filename_writer.writerow([genome_status_selected])
#     for element in data_list:
#         filename_writer.writerow([element])
#     filename_fh.close()




#TODO unit test - this may no longer be needed.
# #Ensure the output filename is unique
# def create_unique_filename(filename_directory,filename_base,filename_ext):
#
#     file_exists = True
#     rename_counter = 0
#     while file_exists == True:
#
#         if rename_counter == 0:
#             unique_filename = filename_base + filename_ext
#         else:
#             unique_filename = filename_base + '_' + str(rename_counter) + filename_ext
#
#         unique_filename_path = os.path.join(filename_directory,unique_filename)
#         file_exists = os.path.isfile(unique_filename_path)
#
#         if file_exists == True:
#             print 'Warning: duplicate output file:'
#             print unique_filename_path
#             print 'Filename will be modified.'
#             raw_input('Press ENTER to proceed')
#             rename_counter += 1
#
#     return unique_filename_path


























#Functions that Christian has probably created

#TODO unit test - this may have already been set up
# #Exits MySQL
# def mdb_exit(message):
# 	write_out(output_file,"\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
# 	write_out(output_file,message)
# 	write_out(output_file,"\nThe import script did not complete.")
# 	write_out(output_file,"\nExiting MySQL.")
# 	cur.execute("ROLLBACK")
# 	cur.execute("SET autocommit = 1")
# 	cur.close()
# 	con.close()
# 	write_out(output_file,"\nExiting import script.")
# 	output_file.close()
# 	sys.exit(1)





#TODO unit test - this may already have been set up
#Set up MySQL parameters
# def setup_mysql():
#
#     mysqlhost = 'localhost'
#     print "\n\n"
#     username = getpass.getpass(prompt='mySQL username:')
#     print "\n\n"
#     password = getpass.getpass(prompt='mySQL password:')
#     print "\n\n"












###
