"""Misc. functions to manipulate tickets."""

import csv
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import phagesdb
from pdm_utils.functions import run_modes
from pdm_utils.classes import ticket
from pdm_utils.classes import genome
import pathlib

def export_ticket_data(list_of_data_dicts, file_path, headers):
    """Save a dictionary of data to file using specified column headers.

    Ensures the output file contains a specified number of columns,
    and it ensures the column headers are exported as well."""

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    # with open(file_path, "w") as file_handle:
    with file_path.open("w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        file_writer.writerow(headers_dict)
        for data_dict in list_of_data_dicts:
            file_writer.writerow(data_dict)


def retrieve_ticket_data(filepath):
    """Open file and retrieve the ticket data.

    :param filename:
        Path to file containing import ticket data, where
        each row is assumed to contain column names.
    :type filename: Path
    :returns:
        A list of elements, where each element is a dictionary
        representing one row of data.
        Each key is a column name and each value is the data
        stored in that field.
    :rtype: list
    """
    list_of_ticket_data = []
    # with open(filename,'r') as file:
    with filepath.open(mode='r') as file:
        file_reader = csv.DictReader(file)
        for dict in file_reader:
            list_of_ticket_data.append(dict)
    return list_of_ticket_data


# TODO review unittests - are they sufficient?
def modify_import_data(data_dict, required_keys, optional_keys, keywords):
    """Checks and modifies a data dictionary to conform to requirements
    for a ticket object.

    :type data_dict: dict
    :returns:
        If successful, a Ticket object is returned.
        If the keys are not correct, None object is returned.
    :rtype: Ticket or None
    """
    # Some import table fields are required. Others are optional.
    # First check if all required fields are present.
    missing_req_fields = required_keys - data_dict.keys()
    # Second, check if there are any fields that are not expected.
    extra_fields = data_dict.keys() - required_keys - optional_keys

    # If row in import table has more fields than headers, they are
    # assigned None by csv.DictReader.
    # If row in import table contains fields with no data, they are
    # assigned "" by csv.DictReader.
    # "retrieve" = get data from PhagesDB.
    # "retain" = get data from current genome in Phamerator.
    # "none" = there is no applicable data (only for accession and subcluster)
    # "" = source of data should be automatically determined (optional)
    if (len(missing_req_fields) == 0 and len(extra_fields) == 0):

        # Convert None values to "".
        set_empty(data_dict)
        # Confirm that certain keywords are lowercase.
        set_keywords(data_dict, keywords)
        # Determine which fields are not present, and add them to the dict.
        set_missing_keys(data_dict, optional_keys)
        # Certain fields should be lowercase.
        data_dict["type"] = data_dict["type"].lower()
        data_dict["description_field"] = data_dict["description_field"].lower()
        data_dict["run_mode"] = data_dict["run_mode"].lower()
        # For fields that are empty (""), set them with default values.
        set_dict_value(data_dict, "host_genus", "retrieve", "retain")
        set_dict_value(data_dict, "cluster", "retrieve", "retain")
        set_dict_value(data_dict, "subcluster", "retrieve", "retain")
        set_dict_value(data_dict, "annotation_author", "1", "retain")
        set_dict_value(data_dict, "retrieve_record", "1", "retain")
        set_dict_value(data_dict, "annotation_status", "draft", "final")

        # If status = 'draft', accession can remain "".
        if data_dict["annotation_status"] != "draft":
            set_dict_value(data_dict, "accession", "retrieve", "retain")
        return True
    else:
        print(f"The ticket {data_dict} is not formatted correctly.")
        return False


def parse_import_ticket_data(data_dict):
    """Converts import ticket data to a Ticket object.

    :param data_dict:
        A dictionary of data with the following keys:

            0. Import action type
            1. Primary PhageID
            2. Host
            3. Cluster
            4. Subcluster
            5. Status
            6. Annotation Author (int)
            7. Feature field
            8. Accession
            9. Retrieve Record (int)
            10. Run mode

    :type data_dict: dict
    :returns: A pdm_utils Ticket object.
    :rtype: Ticket
    """
    ticket_attributes = constants.IMPORT_TABLE_STRUCTURE["valid_ticket"]
    other_attributes = data_dict.keys() - ticket_attributes

    tkt = ticket.GenomeTicket()
    tkt.data_dict = data_dict
    for attr in ticket_attributes:
        attr_value = data_dict[attr]
        setattr(tkt, attr, attr_value)

    data_retrieve = set()
    data_retain = set()
    data_add = set()

    other_attributes = list(other_attributes)
    x = 0
    while x < len(other_attributes):
        attr = other_attributes[x]
        attr_value = data_dict[attr]
        if attr_value == "retrieve":
            data_retrieve.add(attr)
        elif attr_value == "retain":
            data_retain.add(attr)
        else:
            data_add.add(attr)
        x += 1
    tkt.data_retrieve = data_retrieve
    tkt.data_retain = data_retain
    tkt.data_add = data_add
    return tkt


def set_empty(data_dict):
    """Convert None values to an empty string."""
    for key in data_dict.keys():
        if data_dict[key] is None:
            data_dict[key] = ""

def set_keywords(data_dict, keywords):
    """Convert specific values in a dictionary to lowercase."""
    for key in data_dict.keys():
        value = data_dict[key]
        if isinstance(value, str):
            if value.lower() in keywords:
                data_dict[key] = value.lower()


def set_missing_keys(data_dict, expected_keys):
    """Add a list of keys-values to a dictionary if it doesn't have those keys."""
    missing_keys = expected_keys - data_dict.keys()
    for key in missing_keys:
        data_dict[key] = ""
        # print("added key")
        # print(key)
        # print(data_dict[key])
        # print(type(data_dict[key]))
        # input("check keys")

def set_dict_value(data_dict, key, first, second):
    """Set the value for a specific key based on 'type' key-value.

    It expects that the dictionary contains the indicated key, as well as
    a "type" key."""
    if data_dict[key] == "":
        if data_dict["type"] == "add":
            data_dict[key] = first
        else:
            data_dict[key] = second

# TODO it probably makes more sense to pass the 'run_mode' and
# 'description_field' data into modify_import_data() so that they get added
# to the data_dict. As a result, those values will be exported to the
# revised import table file when the script completes.
def construct_tickets(list_of_data_dict, run_mode_eval_dict, description_field,
                      required_keys, optional_keys, keywords):
    """Construct tickets from parsed data dictionaries."""

    tkt_id = 0
    list_of_tickets = []
    for dict in list_of_data_dict:
        tkt_id += 1
        # Each ticket should contain a distinct eval_flag dictionary object.
        input_run_mode = run_mode_eval_dict["run_mode"]
        input_eval_flag_dict = run_mode_eval_dict["eval_flag_dict"].copy()

        result = modify_import_data(dict, required_keys, optional_keys, keywords)
        if result:
            tkt = parse_import_ticket_data(dict)
            tkt.id = tkt_id

            # Only set description_field and run_mode from parameters
            # if they weren't set within the ticket.
            if tkt.description_field == "":
                tkt.description_field = description_field

            if tkt.run_mode == "":
                tkt.run_mode = input_run_mode
                tkt.eval_flags = input_eval_flag_dict
            else:
                if tkt.run_mode in run_modes.RUN_MODES.keys():
                    tkt.eval_flags = run_modes.get_eval_flag_dict(tkt.run_mode)
            list_of_tickets.append(tkt)
        else:
            print("Unable to create ticket.")
    return list_of_tickets



def identify_duplicates(list_of_tickets, null_set=set()):
    """Compare all tickets to each other to identify ticket conflicts.

    Identifies if the same id, PhageID, and Accession
    is present in multiple tickets.

    :param list_of_tickets:
        A list of pdm_utils Ticket objects.
    :type list_of_tickets: list
    :param null_set:
        A set of values that may be expected to be duplicated, that
        should not throw errors.
    :type null_set: set
    :returns:
        tuple (tkt_id_dupes, phage_id_dupes)
        WHERE
        tkt_id_dupes(set) is a set of duplicate ticket ids.
        phage_id_dupes(set) is a set of duplicate PhageIDs.
    :rtype: tuple
    """
    tkt_id_list = []
    phage_id_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for tkt in list_of_tickets:
        if tkt.id not in null_set:
            tkt_id_list.append(tkt.id)
        if tkt.phage_id not in null_set:
            phage_id_list.append(tkt.phage_id)
    tkt_id_dupe_set = basic.identify_one_list_duplicates(tkt_id_list)
    phage_id_dupe_set = basic.identify_one_list_duplicates(phage_id_list)
    return (tkt_id_dupe_set, phage_id_dupe_set)



# TODO re-evaluate the structure of this function. Replace tickets
# may not instantiate more than one Genome object.
def copy_ticket_to_genome(bndl):
    """Construct genome objects from tickets.

    This function operates on a Bundle object
    instead of a Genome object because some tickets (such as 'replace')
    need to instantiate more than one Genome object.

    :param bndl: A pdm_utils Bundle object.
    :type bndl: Bundle
    """

    tkt = bndl.ticket
    if (tkt.type == "add" or tkt.type == "replace"):
        genome1 = genome.Genome()
        genome1.type = "add"
        genome1.set_id(value=tkt.phage_id)
        genome1.name = tkt.phage_id
        genome1.set_host_genus(tkt.host_genus)
        genome1.set_accession(tkt.accession)
        genome1.annotation_status = tkt.annotation_status
        genome1.set_cluster(tkt.cluster)
        genome1.set_subcluster(tkt.subcluster)
        genome1.set_cluster_subcluster()
        genome1.set_annotation_author(tkt.annotation_author)
        genome1.set_retrieve_record(tkt.retrieve_record)

        bndl.genome_dict[genome1.type] = genome1


    else:
        pass

# TODO this will probably replace the first copy_ticket_to_genome()
def get_genome(tkt, gnm_type=""):
    """Construct a genome object from a ticket.

    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    """
    gnm = genome.Genome()
    gnm.type = gnm_type
    gnm.set_id(value=tkt.phage_id)
    gnm.name = tkt.phage_id
    if "host_genus" in tkt.data_add:
        gnm.set_host_genus(tkt.data_dict["host_genus"])
    if "accession" in tkt.data_add:
        gnm.set_accession(tkt.data_dict["accession"])
    if "annotation_status" in tkt.data_add:
        gnm.annotation_status = tkt.data_dict["annotation_status"]
    if "cluster" in tkt.data_add:
        gnm.set_cluster(tkt.data_dict["cluster"])
    if "subcluster" in tkt.data_add:
        gnm.set_subcluster(tkt.data_dict["subcluster"])
    if "annotation_author" in tkt.data_add:
        gnm.set_annotation_author(tkt.data_dict["annotation_author"])
    if "retrieve_record" in tkt.data_add:
        gnm.set_retrieve_record(tkt.data_dict["retrieve_record"])
    gnm.set_cluster_subcluster()
    return gnm


# TODO re-evaluate if needed, since updates and add/replaces are
# run in different scripts now.
def create_bundle_dict(bndl_list):
    """Create a dictionary of Bundle objects based on their ticket type.

    :param bndl_list: A list of pdm_utils Bundle objects.
    :type bndl_list: list
    :returns:
        A dictionary
        WHERE
        key = ticket type (e.g. 'update', 'add', etc.).
        value = list of Bundle objects.
    :rtype: dict
    """
    type_set = set()
    for bndl in bndl_list:
        type_set.add(bndl.ticket.type)
    ticket_type_dict = {}
    for type in type_set:
        bundle_object_list = []
        index = 0
        while index < len(bndl_list):
            bndl = bndl_list[index]
            if bndl.ticket.type == type:
                bundle_object_list.append(bndl)
            index += 1
        ticket_type_dict[type] = bundle_object_list
    return ticket_type_dict













# TODO implement below. This function may no longer be needed.
# TODO unit test below.
def prepare_tickets(ticket_filename):
    """Parse import table into ticket objects."""

    # Assumes that filename has already been validated.

    # TODO retrieve import data.
    # List of ticket data.


    #Retrieve import info from indicated import table file and
    # read all lines into a list and verify contents are correctly populated.
    #0 = Type of database action to be performed (add, remove, replace, update)
    #1 = New PhageID that will be added to database
    #2 = Host of new phage
    #3 = Cluster of new phage (singletons should be reported as "Singleton")
    #4 = Subcluster of new phage (no subcluster should be reported as "none")
    #5 = Annotation status of new phage
    #6 = Annotation author of the new phage
    #7 = Feature field containing gene descriptions of new phage
    #8 = Accession
    #9 = Run mode
    #10 = PhageID of genome to be removed from the database


    # Parse list of data and construct tickets.

    # Convert data from import file into ticket objects
    ticket_list = parse_import_tickets(import_table_data_list)


    return ticket_list







# TODO this function may no longer be needed. Genome object methods
# can assign the PhageID from either the filename or a flat file record
# field. As a result, all genomes can be matched to tickets using the
# Genome object id.
# def assign_match_strategy(list_of_bundle_objects):
#
#     strategy = ""
#     strategies = set()
#
#     for bundle_obj in list_of_bundle_objects:
#         strategies.add(bundle_obj.ticket.match_strategy)
#
#
#     if len(strategies) > 1:
#         result = "Unable to match genomes in files to tickets, " + \
#                     "since more than one match strategy is indicated."
#         status = "error"
#
#
#     else:
#         strategy = list(strategies)[0]
#         result = "Able to match genomes in files to tickets."
#         status = "correct"
#
#     definition = "Assign matching strategy."
#     evl = eval.Eval("TICKET", definition, result, status)
#
#     return strategy, evl




# TODO this function may no longer be needed. Genome object methods
# can assign the PhageID from either the filename or a flat file record
# field. As a result, all genomes can be matched to tickets using the
# Genome object id.
# TODO fix unit tests.
# def match_genomes_to_tickets(list_of_bundle_objects, all_flat_file_data, key):
#     """Match genome objects parsed from files to tickets."""
#
#     # Determine what the strategy is to match tickets to flat files.
#     strategy, strategy_eval = assign_match_strategy(list_of_bundle_objects)
#
#     # Create list of all ticket identifiers.
#     ticket_id_list = []
#     for bundle_obj in list_of_bundle_objects:
#         ticket_id_list.append(bundle_obj.ticket.phage_id)
#
#
#     # Create list of all genome identifiers based on the strategy.
#     genome_dict = {}
#     genome_id_list = []
#     index1 = 0
#     while index1 < len(all_flat_file_data):
#         genome_object = all_flat_file_data[index1]
#
#         if strategy == "phage_id":
#             genome_id = genome_object.id
#         elif strategy == "filename":
#             genome_id = genome_object.filename
#         else:
#             genome_id = ""
#
#         genome_id_list.append(genome_id)
#         genome_dict[genome_id] = genome_object
#
#         index1 += 1
#
#
#     matched_unique_ids, \
#     ticket_unmatched_unique_ids, \
#     genome_unmatched_unique_ids, \
#     ticket_duplicate_ids, \
#     genome_duplicate_ids = \
#         basic.match_items(ticket_id_list, genome_id_list)
#
#
#     # Match genomes to tickets using the unique identifiers.
#     index2 = 0
#     while index2 < len(list_of_bundle_objects):
#
#         bundle_obj = list_of_bundle_objects[index2]
#         match_id = bundle_obj.ticket.phage_id
#
#         if match_id in matched_unique_ids:
#             genome_object = genome_dict.pop(match_id)
#             bundle_obj.genome_dict[key] = genome_object
#
#         index2 += 1
#
#
#
#     # Check for the various errors that could have been encountered.
#     eval_results = []
#
#     if strategy_eval is not None:
#         eval_results.append(strategy_eval)
#
#     if len(ticket_duplicate_ids) > 0:
#         result = "Unable to match genomes to tickets " + \
#                 "since there are multiple tickets with the same identifier."
#         status = "error"
#         definition = "Match tickets to genome."
#         evl = eval.Eval("TICKET", definition, result, status)
#         eval_results.append(evl)
#
#     if len(genome_duplicate_ids) > 0:
#         result = "Unable to match genomes to tickets " + \
#                 "since there are multiple genomes with the same identifier."
#         status = "error"
#         definition = "Match tickets to genome."
#         evl = eval.Eval("TICKET", definition, result, status)
#         eval_results.append(evl)
#
#     if len(ticket_unmatched_unique_ids) > 0:
#         result = "There is no matching genome for one or more tickets."
#         status = "error"
#         definition = "Match tickets to genome."
#         evl = eval.Eval("TICKET", definition, result, status)
#         eval_results.append(evl)
#
#     if len(genome_unmatched_unique_ids) > 0:
#         result = "There is no matching ticket for one or more genomes."
#         status = "error"
#         definition = "Match tickets to genome."
#         evl = eval.Eval("TICKET", definition, result, status)
#         eval_results.append(evl)
#
#     return eval_results












# TODO this can probably be deleted. It is a more complex compare_tickets
# function that tries to keep track of ticket type.
# # This needs to be changed again. Remove tickets now store
# # the genome to be removed within the phage_id.
# def compare_tickets2(list_of_tickets):
#     """Compare all tickets to each other to determine if
#     there are any ticket conflicts."""
#
#     accession_list = []
#
#     update_primary_id_list = []
#     add_primary_id_list = []
#     replace_primary_id_list = []
#
#     remove_secondary_id_list = []
#     replace_secondary_id_list = []
#
#     result_list = []
#
#     # Create separate lists to check each field for duplications.
#     # Skip "none" values since they are expected to be duplicated.
#     for tkt in list_of_tickets:
#
#         if tkt.phage_id != "none":
#             if tkt.type == "update":
#                 update_primary_id_list.append(tkt.phage_id)
#             elif tkt.type == "add":
#                 add_primary_id_list.append(tkt.phage_id)
#             elif tkt.type == "replace":
#                 replace_primary_id_list.append(tkt.phage_id)
#             else:
#                 pass
#
#         if tkt.secondary_phage_id != "none":
#             if tkt.type == "replace":
#                 replace_secondary_id_list.append(tkt.secondary_phage_id)
#             elif tkt.type == "remove":
#                 remove_secondary_id_list.append(tkt.secondary_phage_id)
#             else:
#                 pass
#
#         if tkt.accession != "none":
#             accession_list.append(tkt.accession)
#
#
#     update_add_primary_id_list = update_primary_id_list + \
#                                 add_primary_id_list
#
#     primary_id_list = update_primary_id_list + \
#                             add_primary_id_list + \
#                             replace_primary_id_list
#
#     secondary_id_list = remove_secondary_id_list + \
#                             replace_secondary_id_list
#
#
#     # Now, iterate through each list to identify duplicate or
#     # conflicting tickets.
#     dupe_set1 = basic.identify_one_list_duplicates(primary_id_list)
#     if len(dupe_set1) > 0:
#         result1 = "Multiple tickets contain the same Primary Phage ID."
#         status1 = "error"
#     else:
#         result1 = "Primary PhageID is fine."
#         status1 = "correct"
#
#     definition1 = "Check if there are multiple tickets with the same Primary PhageID."
#     eval1 = eval.Eval(id = "TICKET", \
#                     definition = definition1, \
#                     result = result1, \
#                     status = status1)
#     result_list.append(eval1)
#
#
#     dupe_set2 = basic.identify_one_list_duplicates(secondary_id_list)
#     if len(dupe_set2) > 0:
#         result2 = "Multiple tickets contain the same Secondary Phage ID."
#         status2 = "error"
#     else:
#         result2 = "Secondary PhageID is fine."
#         status2 = "correct"
#
#     definition2 = "Check if there are multiple tickets with the same Secondary PhageID."
#     eval2 = eval.Eval(id = "TICKET", \
#                     definition = definition2, \
#                     result = result2, \
#                     status = status2)
#     result_list.append(eval2)
#
#
#     dupe_set3 = basic.identify_one_list_duplicates(accession_list)
#     if len(dupe_set3) > 0:
#         result3 = "Multiple tickets contain the same Accession."
#         status3 = "error"
#     else:
#         result3 = "Accession is fine."
#         status3 = "correct"
#
#     definition3 = "Check if there are multiple tickets with the same Accession."
#     eval3 = eval.Eval(id = "TICKET", \
#                     definition = definition3, \
#                     result = result3, \
#                     status = status3)
#     result_list.append(eval3)
#
#
#
#
#
#     # No 'update' or 'add' Primary Phage IDs are expected to be a
#     # 'remove' or 'replace' Secondary Phage IDs.
#     dupe_set4 = basic.identify_two_list_duplicates(
#                         update_add_primary_id_list, \
#                         secondary_id_list)
#     if len(dupe_set4) > 0:
#         result4 = "An update or add Primary Phage ID is also a " + \
#                     "remove or replace Secondary Phage ID."
#         status4 = "error"
#     else:
#         result4 = "There is no PhageID conflict."
#         status4 = "correct"
#
#     definition4 = "Check if an update/add Primary PhageID is used " + \
#                     "as a remove/replace Secondary PhageID."
#     eval4 = eval.Eval(id = "TICKET", \
#                     definition = definition4, \
#                     result = result4, \
#                     status = status4)
#     result_list.append(eval4)
#
#
#     # No 'replace' Primary Phage IDs are expected to be a 'remove'.
#     # Secondary Phage ID.
#     dupe_set5 = basic.identify_two_list_duplicates(replace_primary_id_list, \
#                         remove_secondary_id_list)
#     if len(dupe_set5) > 0:
#         result5 = "A replace Primary Phage ID is also a " + \
#                     "remove Secondary Phage ID."
#         status5 = "error"
#     else:
#         result5 = "There is no PhageID conflict."
#         status5 = "correct"
#
#     definition5 = "Check if a replace Primary PhageID is used " + \
#                     "as a remove Secondary PhageID."
#     eval5 = eval.Eval(id = "TICKET", \
#                     definition = definition5, \
#                     result = result5, \
#                     status = status5)
#     result_list.append(eval5)
#
#
#
#     # Only return results that are Eval objects.
#     # list_of_evals = []
#     # for evl in result_list:
#     #     if evl is not None:
#     #         list_of_evals.append(evl)
#
#     # TODO this used to return only error evals, but now it returns all
#     # evals, so this will need to be corrected in the unit tests.
#     return result_list





###
