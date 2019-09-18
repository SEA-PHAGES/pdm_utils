"""Misc. base/simple functions. These should not require import of other
modules in this package to prevent circular imports."""

from pdm_utils.constants import constants
import os
import csv


def find_expression(expression, list_of_items):
    """Counts the number of items with matches to a regular expression.

    :param expression: Regular expression object
    :type expression: re
    :param list_of_items:
        List of items that will be searched with the regular expression.
    :type list_of_items: list
    :returns:
        Number of times the regular expression was identified in the list.
    :rtype: int
    """
    search_tally = 0
    for element in list_of_items:
        search_result = expression.search(element)
        if search_result:
            search_tally += 1
    return search_tally


def edit_suffix(value, option, suffix=constants.NAME_SUFFIX):
    """Adds or removes the indicated suffix to an input value.

    :param value: Value that will be edited.
    :type value: str
    :param option:
        Indicates what to do with the value and suffix ('add', 'remove').
    :type option: str
    :param suffix: The suffix that will be added or removed.
    :type suffix: str
    :returns:
        The edited value. The suffix is not added if the input
        value already has the suffix.
    :rtype: str
    """

    if option.lower() == "add":
        if not value.lower().endswith(suffix.lower()):
            value = value + suffix
    elif option.lower() == "remove":
        if value.lower().endswith(suffix.lower()):
            value = value.strip(suffix)
    else:
        pass
    return value


def reformat_strand(input_value, format, case=False):
    """Converts common strand orientation formats.


    :param input_value: Value that will be edited.
    :type input_value: str, int
    :param format:
        Indicates how the value should be edited.
        Valid format types include:
        'fr_long' ('forward', 'reverse')
        'fr_short' ('f', 'r')
        'fr_abbrev1' ('for', 'rev')
        'fr_abbrev2' ('fwd', 'rev')
        'tb_long' ('top', 'bottom')
        'tb_short' ('t', 'b')
        'wc_long' ('watson', 'crick')
        'wc_short' ('w','c')
        'operator' ('+', '-')
        'numeric' (1, -1).
    :type format: str
    :param case:
        Indicates whether the output value should be capitalized.
    :type case: bool
    :returns: The re-formatted value as indicated by 'format'.
    :rtype: str, int
    """

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


def reformat_coordinates(left, right, current, new):
    """Converts common coordinate formats.

    The type of coordinate formats include:

    '0_half_open':

        0-based half-open intervals that is the common format
        for BAM files and UCSC Browser database. This format seems
        to be more efficient when performing genomics
        computations.

    '1_closed':

        1-based closed intervals that is the common format
        for the Phamerator Database, UCSC Browser,
        the Ensembl genomics database,
        VCF files, GFF files. This format seems to be more
        intuitive and used for visualization.

    The function assumes coordinates reflect the left and right
    boundaries (where the left coordinates is smaller than the right
    coordinate), instead of gene start and stop coordinates.

    :param left: Left coordinate
    :type left: int
    :param right: Left coordinate
    :type right: int
    :param current: Indicates the indexing format of the input coordinates.
    :type current: str
    :param new: Indicates the indexing format of the output coordinates.
    :type new: str
    :returns: The re-formatted left and right coordinates.
    :rtype: int
    """
    format_set = set(["0_half_open", "1_closed"])
    if (current in format_set and new in format_set):
        if current == "0_half_open":
            if new == "1_closed":
                new_left = left + 1
                new_right = right
            else:
                new_left = left
                new_right = right
        else:
            if new == "0_half_open":
                new_left = left - 1
                new_right = right
            else:
                new_left = left
                new_right = right
    else:
        # new_left = left
        # new_right = right
        raise ValueError("Format for CDS coordinate formats must "
                         "either be '0_half_open' or '1_closed'")
    return (new_left, new_right)


def check_empty(value, lower=True):
    """Checks if the value represents a null value.

    :param value: Value to be checked against the empty set.
    :type value: misc.
    :param lower:
        Indicates whether the input value should be lowercased
        prior to checking.
    :type lower: bool
    :returns: Indicates whether the value is present in the empty set.
    :rtype: bool
    """
    if (lower == True and isinstance(value, str)):
        value = value.lower()
    if value in constants.EMPTY_SET:
        result = True
    else:
        result = False
    return result


def convert_empty(input_value, format, upper=False):
    """Converts common null value formats.


    :param input_value: Value to be re-formatted.
    :type input_value: str, int, datetime
    :param format:
        Indicates how the value should be edited.
        Valid format types include:
        'empty_string' = ''
        'none_string' = 'none'
        'null_string' = 'null'
        'none_object' = None
        'na_long' = 'not applicable'
        'na_short' = 'na'
        'n/a' = 'n/a'
        'zero_string' = '0'
        'zero_num' = 0
        'empty_datetime_obj' = datetime object with arbitrary date, '1/1/0001'
    :type format: str
    :param upper:
        Indicates whether the output value should be uppercased.
    :type upper: bool
    :returns: The re-formatted value as indicated by 'format'.
    :rtype: str, int, datetime
    """
    format_dict = {"empty_string": "",
                    "none_string": "none",
                    "null_string": "null",
                    "none_object": None,
                    "na_long": "not applicable",
                    "na_short": "na",
                    "n/a": "n/a",
                    "zero_string": "0",
                    "zero_num": 0,
                    "empty_datetime_obj": constants.EMPTY_DATE}
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
    """Reformat a gene description.

    :param description: Input value to be reformatted.
    :type description: str
    :returns:
        tuple (description, processed_description)
        WHERE
        description(str) is the original value stripped of leading
        and trailing whitespace.
        processed_description(str) is the reformatted value, in which
        non-informative/generic data is removed.
    :rtype: tuple
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
    """Identify unique and non-unique items in a list.

    :param complete_list: List of items that will be evaluated.
    :type complete_list: list
    :returns:
        tuple (unique_set, duplicate_set)
        WHERE
        unique_set(set) is a set of all unique/non-duplicated items.
        duplicate_set(set) is a set of non-unique/duplicated items.
        non-informative/generic data is removed.
    :rtype: tuple
    """
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


# TODO implement? this may no longer be needed.
def identify_nested_items(complete_list):
    """Identify nested and non-nested two-element tuples in a list.

    :param complete_list: List of tuples that will be evaluated.
    :type complete_list: list
    :returns:
        tuple (not_nested_set, nested_set)
        WHERE
        not_nested_set(set) is a set of non-nested tuples.
        nested_set(set) is a set of nested tuples.
    :rtype: tuple
    """
    not_nested_set = set()
    nested_set = set()
    return (not_nested_set, nested_set)




def trim_characters(string):
    """Remove leading and trailing generic characters from a string.

    :param string:
        Value that will be trimmed.
        Characters that will be removed include: '.', ',', ';', '-', '_'.
    :type string: str
    :returns: Edited value.
    :rtype: str
    """
    generic_characters = set([".", ",", ";", "-", "_"])
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
    """Parse string of text to identify the phage name and host genus.

    :param description: Input value to be parsed.
    :type description: str
    :returns:
        tuple (name, host_genus)
        WHERE
        name(str) is the parsed phage name.
        host_genus(str) is the parsed host_genus.
    :rtype: tuple
    """
    generic_words = set(["complete", "genome", "phage", "unclassified"])
    host_genus = ""
    name = ""
    description = description.strip()
    split_description = description.split()

    # Trim leading and trailing generic characters.
    index = 0
    while index < len(split_description):
        split_description[index] = trim_characters(split_description[index])
        index += 1

    # Iterate through the list of processed words and attempt to
    # identify the host_genus name and phage name.
    index = 0
    while index < len(split_description):
        word = split_description[index]
        word_lower = word.lower()

        # Attempt to identify the host_genus.

        # Sometimes the host_genus name is the word preceding 'phage' or 'virus'.
        # e.g. 'Mycobacterium phage'.
        if index > 0:
            if (word_lower == "phage" or word_lower == "virus"):
                host_genus = split_description[index - 1]
        # Sometimes the host_genus name is merged with 'phage'.
        # e.g. 'Mycobacteriophage'.
        elif (len(word) > 5 and word_lower[-5:] == "phage"):
            host_genus = word[:-5]
        else:
            pass

        # Attempt to identify the phage.

        # If there is one non-generic word in the string, assign it as the
        # phage name.
        if len(split_description) == 1:
            if (len(word) > 5 and word_lower[-5:] == "phage"):
                host_genus = word[:-5]
            elif word_lower not in generic_words:
                name = word
            else:
                pass

        # Sometimes the phage name is the word following 'phage' or 'virus'.
        # e.g. 'Mycobacterium phage Trixie' or 'Mycobacteriophage Trixie'
        if index < (len(split_description) - 1):
            if (word_lower[-5:] == "phage" or \
                word_lower[-5:] == "virus"):
                name = split_description[index + 1]
        index += 1
    return (name, host_genus)


def compare_sets(set1, set2):
    """Compute the intersection and differences between two sets.

    :param set1: The first input set.
    :type set1: set
    :param set2: The second input set.
    :type set2: set
    :returns:
        tuple (set_intersection, set1_diff, set2_diff)
        WHERE
        set_intersection(set) is the set of shared values.
        set1_diff(set) is the set of values unique to the first set.
        set2_diff(set) is the set of values unique to the second set.
    :rtype: tuple
    """
    set_intersection = set1 & set2
    set1_diff = set1 - set2
    set2_diff = set2 - set1
    return (set_intersection, set1_diff, set2_diff)


def match_items(list1, list2):
    """Match values of two lists and return several results.

    :param list1: The first input list.
    :type list1: list
    :param list2: The second input list.
    :type list2: list
    :returns:
        tuple (matched_unique_items, set1_unmatched_unique_items,
        set2_unmatched_unique_items, set1_duplicate_items,
        set2_duplicate_items)
        WHERE
        matched_unique_items(set) is the set of matched unique values.
        set1_unmatched_unique_items(set) is the set of
        unmatched unique values from the first list.
        set2_unmatched_unique_items(set) is the set of
        unmatched unique values from the second list.
        set1_duplicate_items(set) is the the set of
        duplicate values from the first list.
        set2_duplicate_items(set) is the set of
        unmatched unique values from the second list.
    :rtype: tuple
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
    """Split a string based on alphanumeric characters.

    Iterates through a string, identifies the first position
    in which the character is a digit, and creates two strings at this
    position.

    :param string: The value to be split.
    :type string: str
    :returns:
        tuple (left, right)
        WHERE
        left(str) is the left portion of the input value prior to the
        first numeric character and only contains alphabetic characters
        (or will be '').
        right(str) is the right portion of the input value after the
        first numeric character and only contains numeric characters
        (or will be '').
    :rtype: tuple
    """
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


def compare_cluster_subcluster(cluster, subcluster):
    """Check if a cluster and subcluster designation are compatible.

    :param cluster:
        The cluster value to be compared.
        'Singleton' and 'UNK' are lowercased.
    :type cluster: str
    :param subcluster: The subcluster value to be compared.
    :type subcluster: str
    :returns:
        The result of the evaluation, indicating whether the two
        values are compatible.
    :rtype: bool
    """
    result = True
    # If Singleton or Unknown Cluster, there should be no Subcluster
    if (cluster.lower() == "singleton" or
        cluster == "UNK" or
        cluster.lower() == "none"):
        if subcluster != "none":
            result = False

    # If not Singleton or Unknown or none, then if
    # there is a Subcluster designation, the Cluster should be part
    # of Subcluster data and the remainder should be a digit
    elif subcluster != "none":
        left, right = split_string(subcluster)
        if (left != cluster or right.isdigit() == False):
            result = False
    else:
        pass
    return result


def identify_one_list_duplicates(item_list):
    """Identify duplicate items within a list.

    :param item_list: The input list to be checked.
    :type item_list: list
    :returns: The set of non-unique/duplicated items.
    :rtype: set
    """
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

    :param item1_list: The first input list to be checked.
    :type item1_list: list
    :param item2_list: The second input list to be checked.
    :type item2_list: list
    :returns:
        The set of non-unique/duplicated items between the two lists
        (but not duplicate items within each list).
    :rtype: set
    """
    item1_set = set(item1_list)
    item2_set = set(item2_list)
    item3_set = item1_set & item2_set
    return item3_set


def check_value_expected_in_set(value, set1, expect=True):
    """Check if a value is present within a set and if it is expected.

    :param value: The value to be checked.
    :type value: misc.
    :param set1: The reference set of values.
    :type set1: set
    :param expect: Indicates if 'value' is expected to be present in 'set1'.
    :type expect: bool
    :returns: The result of the evaluation.
    :rtype: bool
    """
    # Check if the value is present in the set.
    if value in set1:
        present = True
    else:
        present = False
    # Compare the presence/absence with what was expected.
    if (expect and present):
        result = True
    elif (not expect and not present):
        result = True
    else:
        result = False
    return result


def check_value_in_two_sets(value, set1, set2):
    """Check if a value is present within two sets.

    :param value: The value to be checked.
    :type value: misc.
    :param set1: The first reference set of values.
    :type set1: set
    :param set2: The second reference set of values.
    :type set2: set
    :returns:

        The result of the evaluation, indicating whether the
        value is present within:

            1. only the 'first' set
            2. only the 'second' set
            3. 'both' sets
            4. 'neither' set

    :rtype: str
    """
    present1 = False
    present2 = False
    if value in set1:
        present1 = True
    if value in set2:
        present2 = True
    if (present1 and present2):
        result = "both"
    elif (present1 and not present2):
        result = "first"
    elif (not present1 and present2):
        result = "second"
    else:
        result = "neither"
    return result



# TODO this is probably no longer needed.
# def convert_author(input_value):
#     """Converts author string to author integer using the
#     author dictionary.
#     """
#     input_value = input_value.lower()
#     if input_value in constants.AUTHOR_DICTIONARY[1]:
#         new_value = 1
#     else:
#         new_value = 0
#     return new_value


def lower_case(value):
    """Return the value lowercased if it is within a specific set of values.

    :param value: The value to be checked.
    :type value: str
    :returns:
        The lowercased value if it is equivalent to
        'none', 'retrieve', or 'retain'.
    :rtype: str
    """
    lower_set = set(["none", "retrieve", "retain"])
    if isinstance(value, str):
        if value.lower() in lower_set:
            value = value.lower()
    return value


def close_files(list_of_filehandles):
    """Closes all the files in a list of open file handles.

    :param list_of_filehandles: A list of open file handles.
    :type list_of_filehandles: list
    """
    for handle in list_of_filehandles:
        handle.close()


def ask_yes_no(prompt="", response_attempt=1):
    """Function to get the user's yes/no response to a question.

    Accepts variations of yes/y, true/t, no/n, false/f.

    :param prompt: the question to ask the user.
    :type prompt: str
    :param response_attempt:

        The number of the number of attempts allowed before the
        function exits. This prevents the script from getting stuck in a loop.

    :type response_attempt: int
    :returns:

        default is False (e.g. user hits Enter w/o typing
        anything else), but variations of yes or true responses will return
        True instead.

    :rtype: bool, None
    """
    response = None
    response_valid = False
    while response_valid is False and response_attempt > 0:
        response = get_input(prompt)
        response_attempt -= 1
        if response.lower() in set(["yes", "y", "t", "true"]):
            response = True
            response_valid = True
        elif response.lower() in set(["no", "n", "f", "false", ""]):
            response = False
            response_valid = True
        else:
            response = None
            print("Invalid response.")
    return response


# TODO the associated tests can probably be re-factored to not rely on this.
def get_input(prompt=""):
    """Wrapper function to improve testing.

    :param prompt: the question to ask the user.
    :type prompt: str
    :returns: The response from the user.
    :rtype: misc.
    """
    return input(prompt)


def identify_files(path_to_folder):
    """Create a list of filenames from an indicated directory.

    :param path_to_folder: A valid directory path.
    :type path_to_folder: str
    :returns: List of valid files in the directory.
    :rtype: list
    """
    files_in_folder = []
    for item in os.listdir(path_to_folder):
        item_path = os.path.join(path_to_folder, item)
        if os.path.isfile(item_path):
            files_in_folder.append(item)
    return files_in_folder


def expand_path(input_path):
    """Convert a non-absolute path into an absolute path.

    :param input_path: The path to be expanded.
    :type input_path: str
    :returns: The expanded path.
    :rtype: str
    """
    # "~" needs to be expanded first
    home_dir = os.path.expanduser("~")
    if input_path[0] == "~":
        expanded_path = home_dir + input_path[1:]
    else:
        expanded_path = input_path
    # os.path.abspath will resolve ./ or ../ present in the path
    expanded_path = os.path.abspath(expanded_path)
    return expanded_path


def verify_path(filepath, kind=None):
    """Verifies that a given path exists.

    :param filepath: full path to the desired file/directory.
    :type filepath: str
    :param kind:

        ("file", "dir"), corresponding with paths to be
        checked as either files or directories.

    :type kind: str
    :return Boolean: True if path is verified, False otherwise.
    """
    if kind == "file":
        if os.path.isfile(filepath) is True:
            return True
        else:
            return False
    elif kind == "dir":
        if os.path.isdir(filepath) is True:
            return True
        else:
            return False
    elif kind is None:
        if os.path.exists(filepath) is True:
            return True
        else:
            return False
    else:
        # TODO not sure if this should print a statement, or if it
        # should only return False.
        print("{} is not a valid kind for this function.".format(kind),
              "Please try again",
              "use one of (None, dir, file).")
        return False


def make_new_dir(output_dir, new_dir, attempt=1):
    """Make a new directory.

    Checks to verify the new directory name is valid and does not
    already exist. If it already exists, it attempts to extend
    the name with an integer suffix.

    :param output_dir:
        Full path to the directory where the new directory will be created.
    :type output_dir: str
    :param new_dir: Name of the new directory to be created.
    :type new_dir: str
    :param attempt: Number of attempts to create the directory.
    :type attempt: int
    :returns:
        If successful, the name of the created directory.
        If unsuccessful, empty string "".
    :rtype: str
    """
    valid = False
    count = 0
    while (not valid and count < attempt):
        if count > 0:
            new_dir_mod = new_dir + "_" + str(count)
        else:
            new_dir_mod = new_dir
        new_path = os.path.join(output_dir, new_dir_mod)
        if not verify_path(new_path, "dir"):
            valid = True
            os.mkdir(new_path)
        count += 1
    if not valid:
        return ""
    else:
        return new_dir_mod


def parse_flag_file(flag_file):
    """Parse a file to an evaluation flag dictionary.

    :param flag_file:
        A two-column csv-formatted file
        WHERE
        1. evaluation flag
        2. 'True' or 'False'
    :type flag_file: str
    :returns:
        A dictionary
        WHERE
        keys (str) are evaluation flags
        values (bool) indicate the flag setting
        Only flags that contain boolean values are returned.
    :rtype: dict
    """
    eval_flags = {}
    with open(flag_file,'r') as file:
        file_reader = csv.reader(file)
        for row in file_reader:
            if row[1].lower() == "true":
                row[1] = True
                eval_flags[row[0].lower()] = row[1]
            elif row[1].lower() == "false":
                row[1] = False
                eval_flags[row[0].lower()] = row[1]
            else:
                pass
    return eval_flags





#TODO Unit test below






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
#     run_type_options = [
#     	'none',\
#     	'test',\
#     	'production']
#     print '\n\nThe following run types are available:\n'
#     #print '0: ' + run_type_options[0]
#     # print '1: ' + run_type_options[1]
#     # print '2: ' + run_type_options[2]
#     print '1: ' + run_type_options[1] + ' (checks flat files for accuracy, but the database is not changed.)'
#     print '2: ' + run_type_options[2] + ' (after testing files, the database is updated.)'
#     run_type = select_option(
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
