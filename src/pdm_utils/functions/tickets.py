"""Misc. functions to manipulate tickets."""

from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import phagesdb
from pdm_utils.classes import ticket
from pdm_utils.classes import genome







def parse_import_ticket_data(tkt, data_list,
                             expected_size = constants.IMPORT_TABLE_SIZE,
                             id = "", direction="list_to_ticket"):
    """Converts import ticket data between a list and Ticket object formats.
    'tkt' is a Ticket object.
    'data_list' is a list.
    'expected_size' indicates how many elements should be in the data_list.
    'direction' indicates whether data in the list should populate a
    Ticket object ('list_to_ticket') or data in a Ticket object should
    populate a list ('ticket_to_list').
    The expected data structure of the data list:
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
    """
    # Note: by providing the ability to transfer data from both
    # a list-to-ticket and a ticket-to-list, it helps ensure that
    # conversion between the two formats is consistent.

    # Verify the row of information has the correct number of fields to parse.
    if len(data_list) == expected_size:

        if direction == "list_to_ticket":
            tkt._parsed_fields = len(data_list)
            tkt.id = id


        if direction == "list_to_ticket":
            tkt.set_type(data_list[0])
        elif direction == "ticket_to_list":
            data_list[0] = tkt.type
        else:
            pass


        if direction == "list_to_ticket":
            tkt.set_description_field(data_list[7])
        elif direction == "ticket_to_list":
            data_list[7] = tkt.description_field
        else:
            pass


        if direction == "list_to_ticket":
            tkt.set_run_mode(data_list[11])
        elif direction == "ticket_to_list":
            data_list[11] = tkt.run_mode
        else:
            pass


        # This data will eventually populate a Genome object.

        if direction == "list_to_ticket":
            tkt.set_phage_id(data_list[1])
        elif direction == "ticket_to_list":
            data_list[1] = tkt.phage_id
        else:
            pass

        if direction == "list_to_ticket":
            tkt.set_host(data_list[2])
        elif direction == "ticket_to_list":
            data_list[2] = tkt.host_genus
        else:
            pass

        if direction == "list_to_ticket":
            tkt.set_cluster(data_list[3])
        elif direction == "ticket_to_list":
            data_list[3] = tkt.cluster
        else:
            pass

        if direction == "list_to_ticket":
            tkt.set_subcluster(data_list[4])
        elif direction == "ticket_to_list":
            data_list[4] = tkt.subcluster
        else:
            pass

        if direction == "list_to_ticket":
            tkt.set_annotation_status(data_list[5])
        elif direction == "ticket_to_list":
            data_list[5] = str(tkt.annotation_status)
        else:
            pass

        if direction == "list_to_ticket":
            # Convert to integer if possible.
            try:
                data_list[6] = int(data_list[6])
            except:
                pass
            tkt.set_annotation_author(data_list[6])

        elif direction == "ticket_to_list":
            # Convert to string.
            data_list[6] = str(tkt.annotation_author)
        else:
            pass

        if direction == "list_to_ticket":
            tkt.set_accession(data_list[8])
        elif direction == "ticket_to_list":
            data_list[8] = tkt.accession
        else:
            pass

        if direction == "list_to_ticket":
            # Convert to integer if possible.
            try:
                data_list[9] = int(data_list[9])
            except:
                pass
            tkt.set_annotation_qc(data_list[9])
        elif direction == "ticket_to_list":
            # Convert to string.
            data_list[9] = str(tkt.annotation_qc)
        else:
            pass

        if direction == "list_to_ticket":
            # Convert to integer if possible.
            try:
                data_list[10] = int(data_list[10])
            except:
                pass
            tkt.set_retrieve_record(data_list[10])
        elif direction == "ticket_to_list":
            # Convert to string.
            data_list[10] = str(tkt.retrieve_record)
        else:
            pass








def parse_import_tickets(list_of_lists):
    """Parses lists of lists of data from csv file and converts to
    group of parsed tickets. It also returns a lost of errors for tickets
    that failed to be parsed.
    """

    counter = 1
    list_of_tickets = []
    for list_of_data in list_of_lists:
        tkt = ticket.GenomeTicket()
        parse_import_ticket_data(tkt, list_of_data, id = counter)
        list_of_tickets.append(tkt)
        counter += 1
    return list_of_tickets





# This function has been substantially simplified if Remove tickets
# store the genome to be removed in the phage_id field, and if
# it does not try to keep track of the types of tickets from which the
# conflict arises.
def compare_tickets(list_of_tickets):
    """Compare all tickets to each other to determine if
    there are any ticket conflicts."""

    accession_list = []
    phage_id_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for tkt in list_of_tickets:

        if tkt.phage_id != "none":
            phage_id_list.append(tkt.phage_id)

        if tkt.accession != "none":
            accession_list.append(tkt.accession)

    # Identify duplicate values in the group of tickets.
    phage_id_dupe_set = basic.identify_one_list_duplicates(phage_id_list)
    accession_dupe_set = basic.identify_one_list_duplicates(accession_list)

    for tkt in list_of_tickets:

        tkt.check_duplicate_phage_id(phage_id_dupe_set)
        tkt.check_duplicate_accession(accession_dupe_set)



def copy_ticket_to_genome(bndl):
    """Construct genome objects and populate them appropriately using data
    from import ticket. This function operates on a Bundle object
    instead of a Genome object because some tickets (such as 'replace')
    need to instantiate more than one Genome object."""

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
        genome1.annotation_qc = tkt.annotation_qc
        genome1.retrieve_record = tkt.retrieve_record

        bndl.genome_dict[genome1.type] = genome1


        # TODO probably no need to create a second genome.
        # if tkt.type == "replace":
        #
        #     genome2 = genome.Genome()
        #     genome2.type = "remove"
        #
        #     bndl.genome_dict[genome2.type] = genome2

    # TODO 'update' ticket option will eventually be deleted.
    elif tkt.type == "update":

        # TODO unit test.
        # gnm = genome.Genome()
        # gnm.type = "update"
        # gnm.set_id(tkt.phage_id)
        # gnm.set_host_genus(tkt.host_genus)
        # gnm.set_accession(tkt.accession)
        # gnm.annotation_status = tkt.annotation_status
        # gnm.set_cluster(tkt.cluster)
        # gnm.set_subcluster(tkt.subcluster)
        # gnm.set_cluster_subcluster()
        # bndl.genome_dict[gnm.type] = gnm
        pass

    # TODO 'remove' ticket option will eventually be deleted.
    elif tkt.type == "remove":

        # TODO unit test.
        # gnm = genome.Genome()
        # gnm.type = "remove"
        # gnm.set_id(tkt.phage_id)
        # bndl.genome_dict[gnm.type] = gnm
        pass

    else:
        pass





def create_bundle_dict(list_of_bundle_objects):
    """Create a dictionary of Bundle objects based on their ticket type.
    Key = ticket type (e.g. update, add, etc.).
    Value = list of Bundle objects."""

    type_set = set()
    for bndl in list_of_bundle_objects:
        type_set.add(bndl.ticket.type)

    ticket_type_dict = {}
    for type in type_set:
        bundle_object_list = []
        index = 0
        while index < len(list_of_bundle_objects):
            bndl = list_of_bundle_objects[index]
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













# TODO probably no longer needed now that parse_import_ticket_data is
# reversible.
# def parse_import_ticket(
#                         tkt,
#                         data_list,
#                         expected_size = constants.IMPORT_TABLE_SIZE,
#                         id = ""):
#     """Parses list of data and creates an import ticket.
#         Expected data structure:
#         0. Import action
#         1. Primary PhageID
#         2. Host
#         3. Cluster
#         4. Subcluster
#         5. Status
#         6. Annotation Author
#         7. Feature field
#         8. Accession
#         9. Annotation QC
#         10. Retrieve Record
#         11. Run mode
#         12. Secondary PhageID
#     """
#
#     tkt._parsed_fields = len(data_list)
#
#     # Verify the row of information has the correct number of fields to parse.
#     if len(data_list) == expected_size:
#
#         tkt.id = id
#         tkt.set_type(data_list[0])
#         tkt.set_description_field(data_list[7])
#         tkt.set_run_mode(data_list[11])
#
#
#         # This data will eventually populate a Genome object.
#         tkt.set_phage_id(data_list[1])
#         tkt.set_host(data_list[2])
#         tkt.set_cluster(data_list[3])
#         tkt.set_subcluster(data_list[4])
#         tkt.set_annotation_status(data_list[5])
#         tkt.set_annotation_author(data_list[6])
#         tkt.set_accession(data_list[8])
#         tkt.set_annotation_qc(data_list[9])
#         tkt.set_retrieve_record(data_list[10])
#
#
#     return tkt




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
