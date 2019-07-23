"""Misc. functions to manipulate tickets."""

from constants import constants
from functions import basic
from functions import phagesdb
from classes import Ticket
from classes import Genome



def parse_import_ticket(
                        ticket,
                        data_list,
                        expected_size = constants.IMPORT_TABLE_SIZE,
                        id = ""):
    """Parses list of data and creates an import ticket.
        Expected data structure:
        0. Import action
        1. Primary PhageID
        2. Host
        3. Cluster
        4. Subcluster
        5. Status
        6. Annotation Author
        7. Feature field
        8. Accession
        9. Annotation QC
        10. Retrieve Record
        11. Run mode
        12. Secondary PhageID
    """

    ticket._parsed_fields = len(data_list)

    # Verify the row of information has the correct number of fields to parse.
    if len(data_list) == expected_size:

        ticket.id = id
        ticket.set_type(data_list[0])
        ticket.set_description_field(data_list[7])
        ticket.set_run_mode(data_list[11])


        # This data will eventually populate a Genome object.
        ticket.set_primary_phage_id(data_list[1])
        ticket.set_host(data_list[2])
        ticket.set_cluster(data_list[3])
        ticket.set_subcluster(data_list[4])
        ticket.set_annotation_status(data_list[5])
        ticket.set_annotation_author(data_list[6])
        ticket.set_accession(data_list[8])
        ticket.set_annotation_qc(data_list[9])
        ticket.set_retrieve_record(data_list[10])

        # This data will eventually populate a second Genome object.
        ticket.set_secondary_phage_id(data_list[12])

    return ticket



def parse_import_tickets(list_of_lists):
    """Parses lists of lists of data from csv file and converts to
    group of parsed tickets. It also returns a lost of errors for tickets
    that failed to be parsed.
    """

    counter = 1
    list_of_tickets = []
    for list_of_data in list_of_lists:
        ticket = Ticket.GenomeTicket()
        parse_import_ticket(ticket, list_of_data, id = counter)
        list_of_tickets.append(ticket)
        counter += 1
    return list_of_tickets





# This function has been substantially simplified if Remove tickets
# store the genome to be removed in the primary_phage_id field, and if
# it does not try to keep track of the types of tickets from which the
# conflict arises.
def compare_tickets(list_of_tickets):
    """Compare all tickets to each other to determine if
    there are any ticket conflicts."""

    accession_list = []
    phage_id_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for ticket in list_of_tickets:

        if ticket.primary_phage_id != "none":
            phage_id_list.append(ticket.primary_phage_id)

        if (ticket.secondary_phage_id != "none" and \
            ticket.primary_phage_id != ticket.secondary_phage_id):

            phage_id_list.append(ticket.secondary_phage_id)

        if ticket.accession != "none":
            accession_list.append(ticket.accession)

    # Identify duplicate values in the group of tickets.
    phage_id_dupe_set = basic.identify_one_list_duplicates(phage_id_list)
    accession_dupe_set = basic.identify_one_list_duplicates(accession_list)

    for ticket in list_of_tickets:

        ticket.check_duplicate_primary_phage_id(phage_id_dupe_set)
        ticket.check_duplicate_secondary_phage_id(phage_id_dupe_set)
        ticket.check_duplicate_accession(accession_dupe_set)



def copy_ticket_to_genome(bundle):
    """Construct genome objects and populate them appropriately using data
    from import ticket. This function operates on a Bundle object
    instead of a Genome object because some tickets (such as 'replace')
    need to instantiate more than one Genome object."""

    ticket = bundle.ticket

    if (ticket.type == "add" or ticket.type == "replace"):
        genome1 = Genome.Genome()
        genome1.type = "add"
        genome1.set_id(ticket.primary_phage_id)
        genome1.phage_name = ticket.primary_phage_id
        genome1.set_host(ticket.host_genus)
        genome1.set_accession(ticket.accession)
        genome1.annotation_status = ticket.annotation_status
        genome1.set_cluster(ticket.cluster)
        genome1.set_subcluster(ticket.subcluster)
        genome1.set_cluster_subcluster()
        genome1.set_annotation_author(ticket.annotation_author)
        genome1.annotation_qc = ticket.annotation_qc
        genome1.retrieve_record = ticket.retrieve_record

        bundle.genome_dict[genome1.type] = genome1

        if ticket.type == "replace":

            genome2 = Genome.Genome()
            genome2.type = "remove"
            genome2.set_id(ticket.secondary_phage_id)

            bundle.genome_dict[genome2.type] = genome2

    # TODO 'update' ticket option will eventually be deleted.
    elif ticket.type == "update":

        # TODO unit test.
        # genome = Genome.Genome()
        # genome.type = "update"
        # genome.set_id(ticket.primary_phage_id)
        # genome.set_host(ticket.host_genus)
        # genome.set_accession(ticket.accession)
        # genome.annotation_status = ticket.annotation_status
        # genome.set_cluster(ticket.cluster)
        # genome.set_subcluster(ticket.subcluster)
        # genome.set_cluster_subcluster()
        # bundle.genome_dict[genome.type] = genome
        pass

    # TODO 'remove' ticket option will eventually be deleted.
    elif ticket.type == "remove":

        # TODO unit test.
        # genome = Genome.Genome()
        # genome.type = "remove"
        # genome.set_id(ticket.primary_phage_id)
        # bundle.genome_dict[genome.type] = genome
        pass

    else:
        pass





def create_bundle_dict(list_of_bundle_objects):
    """Create a dictionary of Bundle objects based on their ticket type.
    Key = ticket type (e.g. update, add, etc.).
    Value = list of Bundle objects."""

    type_set = set()
    for bundle in list_of_bundle_objects:
        type_set.add(bundle.ticket.type)

    ticket_type_dict = {}
    for type in type_set:
        bundle_object_list = []
        index = 0
        while index < len(list_of_bundle_objects):
            bundle = list_of_bundle_objects[index]
            if bundle.ticket.type == type:
                bundle_object_list.append(bundle)
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
#     eval = Eval.Eval("TICKET", definition, result, status)
#
#     return strategy, eval




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
#         ticket_id_list.append(bundle_obj.ticket.primary_phage_id)
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
#         match_id = bundle_obj.ticket.primary_phage_id
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
#         eval = Eval.Eval("TICKET", definition, result, status)
#         eval_results.append(eval)
#
#     if len(genome_duplicate_ids) > 0:
#         result = "Unable to match genomes to tickets " + \
#                 "since there are multiple genomes with the same identifier."
#         status = "error"
#         definition = "Match tickets to genome."
#         eval = Eval.Eval("TICKET", definition, result, status)
#         eval_results.append(eval)
#
#     if len(ticket_unmatched_unique_ids) > 0:
#         result = "There is no matching genome for one or more tickets."
#         status = "error"
#         definition = "Match tickets to genome."
#         eval = Eval.Eval("TICKET", definition, result, status)
#         eval_results.append(eval)
#
#     if len(genome_unmatched_unique_ids) > 0:
#         result = "There is no matching ticket for one or more genomes."
#         status = "error"
#         definition = "Match tickets to genome."
#         eval = Eval.Eval("TICKET", definition, result, status)
#         eval_results.append(eval)
#
#     return eval_results

















# TODO this can probably be deleted. It is a more complex compare_tickets
# function that tries to keep track of ticket type.
# # This needs to be changed again. Remove tickets now store
# # the genome to be removed within the primary_phage_id.
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
#     for ticket in list_of_tickets:
#
#         if ticket.primary_phage_id != "none":
#             if ticket.type == "update":
#                 update_primary_id_list.append(ticket.primary_phage_id)
#             elif ticket.type == "add":
#                 add_primary_id_list.append(ticket.primary_phage_id)
#             elif ticket.type == "replace":
#                 replace_primary_id_list.append(ticket.primary_phage_id)
#             else:
#                 pass
#
#         if ticket.secondary_phage_id != "none":
#             if ticket.type == "replace":
#                 replace_secondary_id_list.append(ticket.secondary_phage_id)
#             elif ticket.type == "remove":
#                 remove_secondary_id_list.append(ticket.secondary_phage_id)
#             else:
#                 pass
#
#         if ticket.accession != "none":
#             accession_list.append(ticket.accession)
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
#     eval1 = Eval.Eval(id = "TICKET", \
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
#     eval2 = Eval.Eval(id = "TICKET", \
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
#     eval3 = Eval.Eval(id = "TICKET", \
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
#     eval4 = Eval.Eval(id = "TICKET", \
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
#     eval5 = Eval.Eval(id = "TICKET", \
#                     definition = definition5, \
#                     result = result5, \
#                     status = status5)
#     result_list.append(eval5)
#
#
#
#     # Only return results that are Eval objects.
#     # list_of_evals = []
#     # for eval in result_list:
#     #     if eval is not None:
#     #         list_of_evals.append(eval)
#
#     # TODO this used to return only error evals, but now it returns all
#     # evals, so this will need to be corrected in the unit tests.
#     return result_list





###
