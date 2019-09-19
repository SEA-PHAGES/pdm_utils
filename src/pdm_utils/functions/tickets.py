"""Misc. functions to manipulate tickets."""

import csv
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import phagesdb
from pdm_utils.classes import ticket
from pdm_utils.classes import genome


def retrieve_ticket_data(filename):
    """Open file and retrieve the ticket data.

    :param filename:
        Path to file containing import ticket data, where
        each row is assumed to contain column names.
    :type filename: str
    :returns:
        A list of elements, where each element is a dictionary
        representing one row of data.
        Each key is a column name and each value is the data
        stored in that field.
    :rtype: list
    """
    list_of_ticket_data = []
    with open(filename,'r') as file:
        file_reader = csv.DictReader(file)
        for dict in file_reader:
            list_of_ticket_data.append(dict)
    return list_of_ticket_data

def parse_import_ticket_data(tkt=None, data_dict=None,
                              direction="dict_to_ticket"):
    """Converts import ticket data between a list and Ticket object formats.

    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param data_dict:
        A dictionary of data with the following keys:

            0. Import action
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
    :param direction:
        Indicates whether data in the dictionary should populate a
        Ticket object ('dict_to_ticket') or data in a Ticket object should
        populate a dictionary ('ticket_to_dict').
    :type direction: str
    :returns:
        If successful, a Ticket or dict object is returned.
        If copying from a dictionary to Ticket, and the keys are not correct,
        None object is returned.
        If an invalid 'direction' is selected, None object is returned.
    :rtype: Ticket, dict, or None
    """
    # By providing the ability to transfer data from both
    # a dictionary-to-ticket and a ticket-to-dictionary, it helps ensure that
    # conversion between the two formats is consistent.
    if direction == "dict_to_ticket":
        key_check = data_dict.keys() ^ constants.IMPORT_TABLE_DICT.keys()

        # If row in import table has more fields than headers, they are
        # assigned None by csv.DictReader.
        # If row in import table contains fields with no data, they are
        # assigned "" by csv.DictReader.
        value_check = set(["", None]) - set(data_dict.values())
        if (len(key_check) == 0 and len(value_check) == 2):
            tkt = ticket.GenomeTicket()
            tkt.id = data_dict["id"]
            tkt.set_type(data_dict["type"])
            tkt.set_description_field(data_dict["description_field"])
            tkt.set_run_mode(data_dict["run_mode"])
            tkt.set_phage_id(data_dict["phage_id"])
            tkt.set_host(data_dict["host_genus"])
            tkt.set_cluster(data_dict["cluster"])
            tkt.set_subcluster(data_dict["subcluster"])
            tkt.set_annotation_status(data_dict["annotation_status"])
            tkt.set_accession(data_dict["accession"])
            tkt.set_annotation_author(data_dict["annotation_author"])
            tkt.set_retrieve_record(data_dict["retrieve_record"])
            return tkt
        else:
            print("The ticket %s is not formatted correctly." % data_dict)
            return None

    elif direction == "ticket_to_dict":
        data_dict = {}
        data_dict["id"] = tkt.id
        data_dict["type"] = tkt.type
        data_dict["description_field"] = tkt.description_field
        data_dict["run_mode"] = tkt.run_mode
        data_dict["phage_id"] = tkt.phage_id
        data_dict["host_genus"] = tkt.host_genus
        data_dict["cluster"] = tkt.cluster
        data_dict["subcluster"] = tkt.subcluster
        data_dict["annotation_status"] = tkt.annotation_status
        data_dict["annotation_author"] = tkt.annotation_author
        data_dict["accession"] = tkt.accession
        data_dict["retrieve_record"] = tkt.retrieve_record
        return data_dict
    else:
        return None














# TODO this may no longer be needed
# def parse_import_ticket_data_list(tkt, data_list,
#                                   expected_size = constants.IMPORT_TABLE_SIZE,
#                                   id = "", direction="list_to_ticket"):
#     """Converts import ticket data between a list and Ticket object formats.
#
#     :param tkt: A pdm_utils Ticket object.
#     :type tkt: Ticket
#     :param data_list:
#         A list of data with the following structure:
#
#             0. Import action
#             1. Primary PhageID
#             2. Host
#             3. Cluster
#             4. Subcluster
#             5. Status
#             6. Annotation Author
#             7. Feature field
#             8. Accession
#             9. Annotation QC
#             10. Retrieve Record
#             11. Run mode
#
#     :type data_list: list
#     :param expected_size:
#         Indicates how many elements should be in 'data_list'.
#     :type expected_size: int
#     :param direction:
#         Indicates whether data in the list should populate a
#         Ticket object ('list_to_ticket') or data in a Ticket object should
#         populate a list ('ticket_to_list').
#     :type direction: str
#     :param id:
#         A value used to populate the ticket's 'id' attribute, if data is
#         copied to a ticket.
#     :type id: str, int
#     """
#     # By providing the ability to transfer data from both
#     # a list-to-ticket and a ticket-to-list, it helps ensure that
#     # conversion between the two formats is consistent.
#
#     # Verify the row of information has the correct number of fields to parse.
#     if len(data_list) == expected_size:
#
#         if direction == "list_to_ticket":
#             tkt._parsed_fields = len(data_list)
#             tkt.id = id
#
#         if direction == "list_to_ticket":
#             tkt.set_type(data_list[0])
#         elif direction == "ticket_to_list":
#             data_list[0] = tkt.type
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_description_field(data_list[7])
#         elif direction == "ticket_to_list":
#             data_list[7] = tkt.description_field
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_run_mode(data_list[11])
#         elif direction == "ticket_to_list":
#             data_list[11] = tkt.run_mode
#         else:
#             pass
#
#         # This data will eventually populate a Genome object.
#         if direction == "list_to_ticket":
#             tkt.set_phage_id(data_list[1])
#         elif direction == "ticket_to_list":
#             data_list[1] = tkt.phage_id
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_host(data_list[2])
#         elif direction == "ticket_to_list":
#             data_list[2] = tkt.host_genus
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_cluster(data_list[3])
#         elif direction == "ticket_to_list":
#             data_list[3] = tkt.cluster
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_subcluster(data_list[4])
#         elif direction == "ticket_to_list":
#             data_list[4] = tkt.subcluster
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_annotation_status(data_list[5])
#         elif direction == "ticket_to_list":
#             data_list[5] = str(tkt.annotation_status)
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             # Convert to integer if possible.
#             try:
#                 data_list[6] = int(data_list[6])
#             except:
#                 pass
#             tkt.set_annotation_author(data_list[6])
#
#         elif direction == "ticket_to_list":
#             # Convert to string.
#             data_list[6] = str(tkt.annotation_author)
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             tkt.set_accession(data_list[8])
#         elif direction == "ticket_to_list":
#             data_list[8] = tkt.accession
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             # Convert to integer if possible.
#             try:
#                 data_list[9] = int(data_list[9])
#             except:
#                 pass
#             tkt.set_annotation_qc(data_list[9])
#         elif direction == "ticket_to_list":
#             # Convert to string.
#             data_list[9] = str(tkt.annotation_qc)
#         else:
#             pass
#
#         if direction == "list_to_ticket":
#             # Convert to integer if possible.
#             try:
#                 data_list[10] = int(data_list[10])
#             except:
#                 pass
#             tkt.set_retrieve_record(data_list[10])
#         elif direction == "ticket_to_list":
#             # Convert to string.
#             data_list[10] = str(tkt.retrieve_record)
#         else:
#             pass

# TODO this may no longer be needed.
# def parse_import_tickets(list_of_lists):
#     """Parses list of lists of data to a list of Tickets.
#
#     :param list_of_lists:
#         A list of items, where each item is a list of data.
#     :type list_of_lists: list
#     :returns: A list of pdm_utils Ticket objects.
#     :rtype: list
#     """
#
#     counter = 1
#     list_of_tickets = []
#     for list_of_data in list_of_lists:
#         tkt = ticket.GenomeTicket()
#         parse_import_ticket_data_list(tkt, list_of_data, id = counter)
#         list_of_tickets.append(tkt)
#         counter += 1
#     return list_of_tickets


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
        tuple (tkt_id_dupes, phage_id_dupes, accession_dupes)
        WHERE
        tkt_id_dupes(set) is a set of duplicate ticket ids.
        phage_id_dupes(set) is a set of duplicate PhageIDs.
        accession_dupes(set) is a set of duplicate accessions.
    :rtype: tuple
    """
    tkt_id_list = []
    accession_list = []
    phage_id_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for tkt in list_of_tickets:
        if tkt.id not in null_set:
            tkt_id_list.append(tkt.id)
        if tkt.phage_id not in null_set:
            phage_id_list.append(tkt.phage_id)
        if tkt.accession not in null_set:
            accession_list.append(tkt.accession)
    tkt_id_dupe_set = basic.identify_one_list_duplicates(tkt_id_list)
    phage_id_dupe_set = basic.identify_one_list_duplicates(phage_id_list)
    accession_dupe_set = basic.identify_one_list_duplicates(accession_list)
    return (tkt_id_dupe_set, phage_id_dupe_set, accession_dupe_set)



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
        genome1.set_annotation_qc(tkt.annotation_qc)
        genome1.set_retrieve_record(tkt.retrieve_record)

        bndl.genome_dict[genome1.type] = genome1


    else:
        pass




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
    #3 = Cluster of new phage (singletons should be reported as "singleton")
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
