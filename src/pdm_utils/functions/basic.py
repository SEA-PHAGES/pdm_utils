"""Misc. base/simple functions. These should not require import of other
modules in this package to prevent circular imports."""

from pdm_utils.constants import constants
import sys
import os
import csv
import getpass
from pathlib import Path


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
            value = value[:-len(suffix)]
    else:
        pass
    return value


def create_indices(input_list, batch_size):
    """Create list of start and stop indices to split a list into batches.

    :param input_list: List from which to generate batch indices.
    :type input_list: list
    :param batch_size: Size of each batch.
    :type batch_size: int
    :returns: List of 2-element tuples (start index, stop index).
    :rtype: list
    """
    index_list = []
    for start in range(0, len(input_list), batch_size):
        if start + batch_size > len(input_list):
            stop = len(input_list)
        else:
            stop = start + batch_size
        tup = (start, stop)
        index_list.append(tup)
    return index_list


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


def reformat_coordinates(start, stop, current, new):
    """Converts common coordinate formats.

    The type of coordinate formats include:

    '0_half_open':

        0-based half-open intervals that is the common format
        for BAM files and UCSC Browser database. This format seems
        to be more efficient when performing genomics
        computations.

    '1_closed':

        1-based closed intervals that is the common format
        for the MySQL Database, UCSC Browser,
        the Ensembl genomics database,
        VCF files, GFF files. This format seems to be more
        intuitive and used for visualization.

    The function assumes coordinates reflect the start and stop
    boundaries (where the start coordinates is smaller than the stop
    coordinate), instead of transcription start and stop coordinates.

    :param start: Start coordinate
    :type start: int
    :param stop: Stop coordinate
    :type stop: int
    :param current: Indicates the indexing format of the input coordinates.
    :type current: str
    :param new: Indicates the indexing format of the output coordinates.
    :type new: str
    :returns: The re-formatted start and stop coordinates.
    :rtype: int
    """
    format_set = set(["0_half_open", "1_closed"])
    if (current in format_set and new in format_set):
        if current == "0_half_open":
            if new == "1_closed":
                new_start = start + 1
                new_stop = stop
            else:
                new_start = start
                new_stop = stop
        else:
            if new == "0_half_open":
                new_start = start - 1
                new_stop = stop
            else:
                new_start = start
                new_stop = stop
    else:
        # new_start = start
        # new_stop = stop
        raise ValueError("Format for CDS coordinate formats must "
                         "either be '0_half_open' or '1_closed'")
    return (new_start, new_stop)


# TODO this may no longer be needed.
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


def reformat_description(raw_description):
    """Reformat a gene description.

    :param raw_description: Input value to be reformatted.
    :type raw_description: str
    :returns:
        tuple (description, processed_description)
        WHERE
        description(str) is the original value stripped of leading
        and trailing whitespace.
        processed_description(str) is the reformatted value, in which
        non-informative/generic data is removed.
    :rtype: tuple
    """
    if raw_description is None:
        raw_description = ""
    raw_description = raw_description.strip()
    processed_description = raw_description.lower()
    split_description = processed_description.split(" ")

    generic_set = {"hypothetical protein",
                   "phage protein",
                   "unknown",
                   "conserved hypothetical protein",
                   "conserved phage protein",
                   "hypothetical phage protein",
                   "hypothetical phage membrane protein",
                   "conserved phage membrane protein",
                   "hypothetical phageprotein",
                   "phage membrane protein",
                   "membrane protein",
                   "structural protein",
                   "regulatory protein",
                   "integral membrane protein",
                   "\\n",
                   "."
                   }

    if processed_description in generic_set:
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
            processed_description = raw_description
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
            processed_description = raw_description
    else:
        processed_description = raw_description
    return (raw_description, processed_description)


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


# TODO this function needs to be improved.
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
    generic_words = {"complete", "genome", "sequence", "phage", "unclassified"}
    host_genus = ""
    name = ""
    description = description.strip()
    split_description = description.split()

    # Trim leading and trailing generic characters.
    for index in range(len(split_description)):
        split_description[index] = trim_characters(split_description[index])

    # Iterate through the list of processed words and attempt to
    # identify the host_genus name and phage name.

    for index in range(len(split_description)):
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
        else:
            if index < (len(split_description) - 1):
                if (word_lower[-5:] == "phage" or
                        word_lower[-5:] == "virus"):
                    # Sometimes the phage name follows 'phage' or 'virus'.
                    # e.g. 'Mycobacterium phage Trixie' or
                    # 'Mycobacteriophage Trixie'
                    name = split_description[index + 1]
                elif len(split_description) == 2:
                    # Sometimes phage name precedes 'Unclassified'.
                    # e.g. 'Trixie_Draft Unclassified'
                    next_lower = split_description[index + 1].lower()
                    if (word_lower not in generic_words and
                            next_lower == "unclassified"):
                        name = word
                else:
                    pass
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


# TODO this function is specific to genome-level data, so should
# it be move as a Genome class method?
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

    # If Singleton or Unknown Cluster, there should be no Subcluster.
    # Subcluster = '' is considered an error since it is the default
    # attribute value when the Genome class is instantiated.
    if (cluster.capitalize() == "Singleton" or cluster == "UNK"):
        if subcluster != "none":
            result = False

    # If not Singleton or Unknown, then if
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


def ask_yes_no(prompt="", response_attempt=1):
    """Function to get the user's yes/no response to a question.

    Accepts variations of yes/y, true/t, no/n, false/f, exit/quit/q.

    :param prompt: the question to ask the user.
    :type prompt: str
    :param response_attempt:

        The number of the number of attempts allowed before the
        function exits. This prevents the script from getting stuck in a loop.

    :type response_attempt: int
    :returns:

        The default is False (e.g. user hits Enter without typing
        anything else), but variations of yes or true responses will return
        True instead. If the response is 'exit' or 'quit', the loop is exited
        and None is returned.

    :rtype: bool, None
    """
    response_valid = False
    while response_valid is False and response_attempt > 0:
        response = input(prompt)
        response_attempt -= 1
        if response.lower() in set(["yes", "y", "t", "true"]):
            response = True
            response_valid = True
        elif response.lower() in set(["no", "n", "f", "false", ""]):
            response = False
            response_valid = True
        elif response.lower() in set(["exit", "quit", "q"]):
            response = None
            response_valid = True
        else:
            response = None
            print("Invalid response. Enter [yes/no/exit].")
    return response


def identify_contents(path_to_folder, kind=None, ignore_set=set()):
    """Create a list of filenames and/or folders from an indicated directory.

    :param path_to_folder: A valid directory path.
    :type path_to_folder: Path
    :param kind:

        ("file", "dir"), corresponding with paths to be
        checked as either files or directories.
    :type kind: str
    :param ignore_set:
        A set of strings representing file or folder names to ignore.
    :type ignore_set: set
    :returns: List of valid contents in the directory.
    :rtype: list
    """
    contents = []
    for item in path_to_folder.iterdir():
        item_path = Path(path_to_folder, item)
        if kind == "file":
            if (item_path.is_file() and item.name not in ignore_set):
                contents.append(item)
        elif kind == "dir":
            if (item_path.is_dir() and item.name not in ignore_set):
                contents.append(item)
        elif kind == None:
            if item.name not in ignore_set:
                contents.append(item)
        else:
            msg = (f"{kind} is not a valid kind (None, dir, file) "
                  "for this function.")
            print(msg)
            # Returns None so that it is clear it wasn't evaluated.
            contents = None
    return contents


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
        print("{} is not a valid kind for this function.".format(kind),
              "Please try again",
              "use one of (None, dir, file).")
        return False


def verify_path2(path, kind=None, expect=True):
    """Verifies that a given path exists.

    :param path: path
    :type path: Path
    :param kind:

        ("file", "dir"), corresponding with paths to be
        checked as either files or directories.

    :type kind: str
    :param expect: Indicates if the path is expected to the indicated kind.
    :type expect: bool
    :returns:
        tuple (result, message)
        WHERE
        result(bool) indicates if the expectation was satisfied.
        message(str) is a description of the result.
    :rtype: tuple
    """
    if kind == "file":
        exists = path.is_file()
        if (exists and not expect):
            result = False
            msg = f"The file {path} already exists."
        elif (not exists and expect):
            result = False
            msg = f"The file {path} does not exist."
        else:
            result = True
            msg = None
    elif kind == "dir":
        exists = path.is_dir()
        if (exists and not expect):
            result = False
            msg = f"The directory {path} already exists."
        elif (not exists and expect):
            result = False
            msg = f"The directory {path} does not exist."
        else:
            result = True
            msg = None
    elif kind == None:
        exists = path.exists()
        if (exists and not expect):
            result = False
            msg = f"The {path} already exists."
        elif (not exists and expect):
            result = False
            msg = f"The {path} does not exist."
        else:
            result = True
            msg = None
    else:
        result = False
        msg = (f"{kind} is not a valid kind (None, dir, file) "
              "for this function.")

    return result, msg


def set_path(path, kind=None, expect=True):
    """Confirm validity of path argument.

    :param path: path
    :type path: Path
    :param kind:

        ("file", "dir"), corresponding with paths to be
        checked as either files or directories.

    :type kind: str
    :param expect: Indicates if the path is expected to the indicated kind.
    :type expect: bool
    :returns: Absolute path if valid, otherwise sys.exit is called.
    :rtype: Path
    """
    path = path.expanduser()
    path = path.resolve()
    result, msg = verify_path2(path, kind=kind, expect=expect)
    if not result:
        print(msg)
        sys.exit(1)
    else:
        return path


def make_new_dir(output_dir, new_dir, attempt=1):
    """Make a new directory.

    Checks to verify the new directory name is valid and does not
    already exist. If it already exists, it attempts to extend
    the name with an integer suffix.

    :param output_dir:
        Full path to the directory where the new directory will be created.
    :type output_dir: Path
    :param new_dir: Name of the new directory to be created.
    :type new_dir: Path
    :param attempt: Number of attempts to create the directory.
    :type attempt: int
    :returns:
        If successful, the full path of the created directory.
        If unsuccessful, None.
    :rtype: Path, None
    """
    valid = False
    count = 0
    while (not valid and count < attempt):
        if count > 0:
            new_dir_mod = new_dir.stem + "_" + str(count)
            new_dir_mod = Path(new_dir_mod)
        else:
            new_dir_mod = new_dir
        new_path = Path(output_dir, new_dir_mod)
        if not new_path.is_dir():
            valid = True
            new_path.mkdir()
        count += 1
    if not valid:
        return None
    else:
        return new_path

# TODO this may not be needed, if a standard config file format is used.
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


def get_user_pwd(user_prompt="Username: ", pwd_prompt="Password: "):
    """Get username and password.

    :param user_prompt: Displayed description when prompted for username.
    :type user_prompt: str
    :param pwd_prompt: Displayed description when prompted for password.
    :type pwd_prompt: str
    :returns:
        tuple (username, password)
        WHERE
        username(str) is the user-supplied username.
        password(str) is the user-supplied password.
    :rtype: tuple
    """
    username = getpass.getpass(prompt=user_prompt)
    password = getpass.getpass(prompt=pwd_prompt)
    return (username, password)


def choose_from_list(options):
    """Iterate through a list of values and choose a value.

    :param options: List of options to choose from.
    :type options: list
    :returns:
        the user select option of None
    :rtype: option or None
    """
    exit = False
    result = False
    option = None
    x = 0
    print(f"These are the options: {', '.join(options)}.")
    while (exit == False and result is False and x < len(options)):
        option = options[x]
        prompt = f"Would you like to select '{option}'? "
        result = ask_yes_no(prompt=prompt, response_attempt=3)
        if (result is None or result is True):
            exit = True
        else:
            x += 1
    if result is True:
        return option
    else:
        return None


def truncate_value(value, length, suffix):
    """Truncate a string.

    :param value: String that should be truncated.
    :type value: str
    :param length: Final length of truncated string.
    :type length: int
    :param suffix: String that should be appended to truncated string.
    :type suffix: str
    :returns: the truncated string
    :rtype: str
    """
    if len(value) > length:
        if len(suffix) > length:
            suffix = suffix[:length]
        length2 = length - len(suffix)
        value = value[:length2] + suffix
    return value


# TODO this needs to be improved.
# TODO unittest.
def select_option(prompt, valid_response_set):
    """Select an option from a set of options.

    :param prompt: Message to display before displaying option.
    :type prompt: str
    :param valid_response_set: Set of valid options to choose.
    :type valid_response_set: set
    :returns: option
    :rtype: str, int
    """

    response_valid = False
    while response_valid == False:
        response = input(prompt)
        if response.isdigit():
            response = int(response)
        else:
            response = response.lower()

        if response in valid_response_set:
            response_valid = True
            if response == "y":
                response  = "yes"
            elif response == "n":
                response  = "no"
        else:
            print("Invalid response.")
    return response


# TODO could probably improve using configparser module.
# TODO unittest
def parse_config_file(path, delimiter="="):
    """Parse a file to return configuration settings."""

    config_dict = {}
    duplicate_config_keys = set()
    with open(path, "r") as fh:
        config_data = fh.readlines()
    config_data = [i.strip() for i in config_data]
    config_data = [i.split(delimiter) for i in config_data]

    for row in config_data:
        if row[0] not in config_dict.keys():
            config_dict[row[0]] = row[1]
        else:
            duplicate_config_keys.add(row[0])

    if len(duplicate_config_keys) > 0:
        config_dict = None
        print("There are duplicate config keys in the file. "
              "Unable to parse config settings.")
    return config_dict


# TODO unittest.
def get_values_from_dict_list(list_of_dicts):
    """Convert a list of dictionaries to a set of the dictionary values.

    :param list_of_dicts: List of dictionaries.
    :type list_of_dicts: list
    :returns: Set of values from all dictionaries in the list.
    :rtype: set
    """
    output_set = set()
    for dict in list_of_dicts:
        output_set = output_set | set(dict.values())
    return output_set

# TODO unittest.
def get_values_from_tuple_list(list_of_tuples):
    """Convert a list of tuples to a set of the tuple values.

    :param list_of_tuples: List of tuples.
    :type list_of_tuples: list
    :returns: Set of values from all tuples in the list.
    :rtype: set
    """
    output_set = set()
    for tup in list_of_tuples:
        output_set.add(tup[0])
    return output_set


def convert_list_to_dict(data_list, key):
    """Convert list of dictionaries to a dictionary of dictionaries

    :param data_list: List of dictionaries.
    :type data_list: list
    :param key: key in each dictionary to become the returned dictionary key.
    :type key: str
    :returns:
        Dictionary of all dictionaries. Returns an empty dictionary if
        all intended keys are not unique.
    :rtype: dict
    """
    data_dict = {}
    for element_dict in data_list:
        if element_dict[key] not in data_dict.keys():
            data_dict[element_dict[key]] = element_dict
    diff = len(data_list) - len(data_dict)
    if diff != 0:
        print("\nUnable to create dictionary since "
              "some intended keys are duplicated.")
        data_dict = {}
    return data_dict


# TODO unittest
def prepare_filepath(folder_path, file_name, folder_name=None):
    """Prepare path to new file.

    :param folder_path: Path to the directory to contain the file.
    :type folder_path: Path
    :param file_name: Name of the file.
    :type file_name: str
    :param folder_name: Name of sub-directory to create.
    :type folder_name: Path
    :returns: Path to file in directory.
    :rtype: Path
    """
    if folder_name is not None:
        subfolder_path = Path(folder_path, folder_name)
        subfolder_path.mkdir()
    else:
        subfolder_path = folder_path
    file_path = Path(subfolder_path, file_name)
    return file_path


def show_progress(current, end, width=50):
    """Updating progress bar that prints to the command line in-line.

    :param current: current value (0-n)
    :param end: end value (n)
    :param width:
        character width for the progress bar (default 50);
        ideally 100 divided by width should naturally be an integer
    :return:
        progress - the integer percent completion calculated
        as the ratio of current to end multiplied by 100
    """
    progress = int(float(current)/end * 100)
    ratio = int(100/width)
    print("\r[{}{}] {}%".format('#' * int(progress / ratio), ' ' * (width - int(progress / ratio)), progress), end="")
    return progress


def export_data_dict(data_dicts, file_path, headers, include_headers=False):
    """Save a dictionary of data to file using specified column headers.

    Ensures the output file contains a specified number of columns,
    and it ensures the column headers are exported as well.

    :param data_dicts:
        list of elements, where each element is a dictionary.
    :type data_dicts: list
    :param file_path: Path to file to export data.
    :type file_path: Path
    :param headers:
        List of strings to define the column order in the file.
        If include_headers is selected, the first row of the file
        will contain each string.
    :type headers: list
    :param include_headers:
        Indicates whether the file should contain a
        row of column names derived from the headers parameter.
    :type include_headers: bool
    """

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    # with open(file_path, "w") as file_handle:
    with file_path.open("w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        if include_headers:
            file_writer.writerow(headers_dict)
        for data_dict in data_dicts:
            file_writer.writerow(data_dict)


def retrieve_data_dict(filepath):
    """Open file and retrieve a dictionary of data.

    :param filepath:
        Path to file containing data and column names.
    :type filepath: Path
    :returns:
        A list of elements, where each element is a dictionary
        representing one row of data.
        Each key is a column name and each value is the data
        stored in that field.
    :rtype: list
    """
    data_dicts = []
    with filepath.open(mode='r') as file:
        file_reader = csv.DictReader(file)
        for dict in file_reader:
            data_dicts.append(dict)
    return data_dicts


def join_strings(input_list, delimiter=" "):
    """Open file and retrieve a dictionary of data.

    :param input_list: List of values to join.
    :type input_list: list
    :param delimiter: Delimiter used between values.
    :type delimiter: str
    :returns: Concatenated values, excluding all None and '' values.
    :rtype: str
    """
    concat_list = []
    for value in input_list:
        if (value is not None and value != ""):
            concat_list.append(value)
    if len(concat_list) > 0:
        string = delimiter.join(concat_list)
    else:
        string = ""
    return string


def merge_set_dicts(dict1, dict2):
    """Merge two dictionaries of sets.

    :param dict1: First dictionary of sets.
    :type dict1: dict
    :param dict2: Second dictionary of sets.
    :type dict2: dict
    :returns:
        Merged dictionary containing all keys from both dictionaries,
        and for each shared key the value is a set of merged values.
    :rtype: dict
    """
    keys = dict1.keys() | dict2.keys()
    dict3 = {}
    for key in keys:
        if key in dict1.keys():
            set1 = dict1[key]
        else:
            set1 = set()

        if key in dict2.keys():
            set2 = dict2[key]
        else:
            set2 = set()

        dict3[key] = set1 | set2
    return dict3
