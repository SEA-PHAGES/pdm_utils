"""Misc. base/simple functions. These should not require import of other
k_phamerate modules to prevent circular imports."""

from datetime import datetime
import os

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




def edit_draft_suffix(value, option):
    """Adds or removes the '_Draft' suffix to a genome name.
    The suffix is not added if the input value already has the suffix."""

    if option.lower() == "add":
        if value[-6:].lower() != "_draft":
            value += "_Draft"

    elif option.lower() == "remove":
        if value[-6:].lower() == "_draft":
            value = value[:-6]

    else:
        pass

    return(value)












#Old 'parse_strand' function and 'parse_strand_for_import' function combined.
# Note: parse_strand_for_import used to convert numeric to long string format.
# E.g. 1 == forward
def reformat_strand(input_value, format, case = False):
    """Converts common strand orientation formats, including:
    'fr_long' ('forward', 'reverse')
    'fr_short' ('f', 'r')
    'fr_abbrev1' ('for', 'rev')
    'fr_abbrev2' ('fwd', 'rev')
    'tb_long' ('top', 'bottom')
    'tb_short' ('t', 'b')
    'wc_long' ('watson', 'crick')
    'wc_short' ('w','c')
    'operator' ('+', '-')
    'numeric' (1, -1)."""

    format_dict = {"fr_long":["forward", "reverse"],
                    "fr_short":["f", "r"],
                    "fr_abbrev1":["for", "rev"],
                    "fr_abbrev2":["fwd", "rev"],
                    "tb_long":["top", "bottom"],
                    "tb_short":["t", "b"],
                    "wc_long":["watson", "crick"],
                    "wc_short":["w", "c"],
                    "operator":["+", "-"],
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




# TODO this function can probably be improved. Maybe it should return
# an eval object?
# TODO also, coordinates in CdsFeature object can be stored as strings,
# so this function needs to take that into account.
def reformat_coordinates(left, right, current, new):
    """Converts common coordinate formats, including:
    '0_half_open' = 0-based half-open intervals that is the common format
                    for BAM files and UCSC Browser database. This format seems
                    to be more efficient when performing genomics
                    computations.
    '1_closed' = 1-based closed intervals that is the common format
                    for the Phamerator Database, UCSC Browser,
                    the Ensembl genomics database,
                    VCF files, GFF files. This format seems to be more
                    intuitive and used for visualization.
    The function assumes coordinates reflect the left and right
    boundaries (where the left coordinates is smaller than the right
    coordinate), instead of gene start and stop coordinates."""


    format_set = set(["0_half_open", "1_closed"])


    if (current in format_set and new in format_set):

        if current == "0_half_open":

            if new == "1_closed":
                new_left = left + 1
                new_right = right
            elif new == "0_half_open":
                new_left = left
                new_right = right
            else:
                new_left = ""
                new_right = ""

        elif current == "1_closed":

            if new == "1_closed":
                new_left = left
                new_right = right
            elif new == "0_half_open":
                new_left = left - 1
                new_right = right
            else:
                new_left = ""
                new_right = ""

        else:
            new_left = ""
            new_right = ""

    else:
        new_left = ""
        new_right = ""


    return new_left, new_right







def check_empty(value, lower = True):
    """Checks if the value represents a null value."""

    empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
    empty_set = set(["",
                    "none",
                    "null",
                    None,
                    "not applicable",
                    "na",
                    "n/a",
                    "0",
                    0,
                    empty_date])

    if (lower == True and isinstance(value, str)):
        value = value.lower()

    if value in empty_set:
        result = True
    else:
        result = False
    return result


def convert_empty(input_value, format, upper = False):
    """Converts common NULL value formats, including:
    'empty_string' = ''
    'none_string' = 'none'
    'null_string' = 'null'
    'none_object' = None
    'na_long' = 'not applicable'
    'na_short' = 'na'
    'n/a' = 'n/a'
    'zero_string' = '0'
    'zero_num' = 0
    'empty_datetime_obj' = datetime object with arbitrary early date, '1/1/0001'
    """

    empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')

    format_dict = {"empty_string": "",
                    "none_string": "none",
                    "null_string": "null",
                    "none_object": None,
                    "na_long": "not applicable",
                    "na_short": "na",
                    "n/a": "n/a",
                    "zero_string": "0",
                    "zero_num": 0,
                    "empty_datetime_obj": empty_date}

    output_values = set(format_dict.values())

    if format in format_dict.keys():

        # Convert input value to ensure it matches the formatting of
        # all possible output values.
        if isinstance(input_value, str):
            new_value = input_value.lower()
        else:
            new_value = input_value

        # Now check to see if the value is present in among all possible
        # output values, and convert.
        if new_value in output_values:
            new_value = format_dict[format]

        else:
            new_value = input_value

    else:
        new_value = input_value

    if isinstance(new_value, str):
        if upper == True:
            new_value = new_value.upper()

    return new_value



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

    return (description, processed_description)







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

    return (unique_set, duplicate_set)




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

    generic_words = set(["complete", "genome", "phage", "unclassified"])

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



def compare_sets(set1, set2):
    """Compute the intersection and differences between two sets."""

    set_intersection = set1 & set2
    set1_diff = set1 - set2
    set2_diff = set2 - set1

    return (set_intersection, set1_diff, set2_diff)


def match_items(list1, list2):
    """Match values of two lists and return several results.
    First, return the set of matched unique values.
    Second, return the set of unmatched unique values from the first list.
    Third, return the set of unmatched unique values from the second list.
    Fourth, return the set of duplicate values from the first list.
    Fifth, return the set of unmatched unique values from the second list.
    """

    # Determine the unique values in each list.
    set1_unique_items, set1_duplicate_items = \
        identify_unique_items(list1)

    set2_unique_items, set2_duplicate_items = \
        identify_unique_items(list2)

    # Create matched and difference sets.
    matched_unique_items, \
    set1_unmatched_unique_items, \
    set2_unmatched_unique_items = \
        compare_sets(set1_unique_items, set2_unique_items)

    return (matched_unique_items,
            set1_unmatched_unique_items,
            set2_unmatched_unique_items,
            set1_duplicate_items,
            set2_duplicate_items)




def split_string(string):
    """Iterates through a string, identifies the first position
    in which the character is a digit, and creates two strings at this
    position. The first string returned contains only alphabetic characters.
    The second string returned contains only numberic characters.
    If there are no numeric characters present,
    the second string will be empty.
    If there are no alphabetic characters present,
    the first string will be empty."""

    left = ""
    right = ""

    if string.isalpha():
        left = string
    elif string.isdigit():
        right = string
    else:
        index = 0
        value = False
        while (value == False and index < len(string)):
            if (string[:index].isalpha() and \
                string[index:].isdigit()):

                left = string[:index]
                right = string[index:]
                value = True

            index += 1

    return (left, right)




















def identify_one_list_duplicates(item_list):
    """Identify duplicate items within a list."""

    duplicate_items = set()
    item_set = set(item_list)
    for item1 in item_set:
        count = 0
        for item2 in item_list:
            if item1 == item2:
                count += 1
        if count > 1:
            duplicate_items.add(item1)

    return duplicate_items

def identify_two_list_duplicates(item1_list, item2_list):
    """Identify duplicate items between two lists.
    It does not identify duplicate items within each list."""

    item1_set = set(item1_list)
    item2_set = set(item2_list)
    item3_set = item1_set & item2_set

    return item3_set










#TODO Unit test below



# TODO unit test.
def identify_files(path_to_folder):
    """Create a list of filenames from an indicated directory."""

    files_in_folder = []
    for item in os.listdir(path_to_folder):
        item_path = os.path.join(path_to_folder, item)

        if os.path.isfile(item_path):
            files_in_folder.append(item)

    return files_in_folder










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
