"""Misc. functions to manipulate tickets."""


from functions import basic
from functions import phagesdb
from classes import Eval
from classes import Ticket
from classes import Genome




def parse_import_ticket(data_list):
    """Parses list of data and creates an import ticket.
        Expected data structure:
        0. Import action
        1. New PhageID
        2. Host
        3. Cluster
        4. Subcluster
        5. Status
        6. Annotation Author
        7. Feature field
        8. Accession
        9. Run mode
        10. PhageID to be removed
    """


    # Verify the row of information has the correct number of fields to parse.
    if len(data_list) != 11:
        message = "Data is not formatted correctly"
        eval_object = Eval.construct_error(message)
        ticket = None

    else:
        eval_object = None

        ticket = Ticket.GenomeTicket()
        ticket.type = data_list[0]
        ticket.description_field = data_list[7]
        ticket.run_mode = data_list[9]


        # TODO this data should populate a Genome object.
        # genome1 = Genome.Genome()
        ticket.primary_phage_id = data_list[1]
        ticket.host = data_list[2]
        ticket.cluster = data_list[3]
        ticket.subcluster = data_list[4]
        ticket.status = data_list[5]
        ticket.annotation_author = data_list[6]
        ticket.accession = data_list[8]


        # TODO this data should populate a second Genome object.
        # genome2 = Genome.Genome()

        ticket.secondary_phage_id = data_list[10]


        # TODO combine all data into a DataGroup object.
        # matched_data = DataGroup.DataGroup()
        # matched_data.ticket = ticket
        # matched_data.genomes_dict["import_genome1"] = genome1
        # matched_data.genomes_dict["import_genome2"] = genome2
        

    return(ticket, eval_object)







def parse_import_tickets(list_of_lists):
    """Parses lists of lists of data from csv file and converts to
    group of parsed tickets. It also returns a lost of errors for tickets
    that failed to be parsed.
    """

    list_of_evals = []
    list_of_tickets = []

    for list_of_data in list_of_lists:

        ticket, eval = parse_import_ticket(list_of_data)

        if ticket is not None:
            list_of_tickets.append(ticket)
        if eval is not None:
            list_of_evals.append(eval)

    return(list_of_tickets, list_of_evals)




#
# #TODO delete once this function is moved to basic.
# def identify_one_list_duplicates(item_list, message):
#     """Verifies there are no duplicate items in a list."""
#
#     item_set = set(item_list)
#     if len(item_set) != len(item_list):
#         eval_object = Eval.construct_error(message)
#     else:
#         eval_object = None
#     return eval_object


#
# #TODO delete once this function is moved to basic.
# def identify_two_list_duplicates(item1_list, item2_list, message):
#     """Verifies there are no duplicate items between two lists."""
#
#     item1_set = set(item1_list)
#     item2_set = set(item2_list)
#     item3_list = item1_set & item2_set
#
#     if len(item3_list) > 0:
#         eval_object = Eval.construct_error(message)
#     else:
#         eval_object = None
#     return eval_object



# TODO unit test now that I have refactored this function.
def compare_tickets(list_of_tickets):
    """Compare all tickets to each other to determine if
    there are any ticket conflicts."""

    accession_list = []

    update_primary_id_list = []
    add_primary_id_list = []
    replace_primary_id_list = []

    remove_secondary_id_list = []
    replace_secondary_id_list = []

    result_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for ticket in list_of_tickets:

        if ticket.primary_phage_id != "none":
            if ticket.type == "update":
                update_primary_id_list.append(ticket.primary_phage_id)
            elif ticket.type == "add":
                add_primary_id_list.append(ticket.primary_phage_id)
            elif ticket.type == "replace":
                replace_primary_id_list.append(ticket.primary_phage_id)
            else:
                pass

        if ticket.secondary_phage_id != "none":
            if ticket.type == "replace":
                replace_secondary_id_list.append(ticket.secondary_phage_id)
            elif ticket.type == "remove":
                remove_secondary_id_list.append(ticket.secondary_phage_id)
            else:
                pass

        if ticket.accession != "none":
            accession_list.append(ticket.accession)


    update_add_primary_id_list = update_primary_id_list + \
                                add_primary_id_list

    primary_id_list = update_primary_id_list + \
                            add_primary_id_list + \
                            replace_primary_id_list

    secondary_id_list = remove_secondary_id_list + \
                            replace_secondary_id_list


    # Now, iterate through each list to identify duplicate or
    # conflicting tickets.
    dupe_set1 = basic.identify_one_list_duplicates(primary_id_list)
    if len(dupe_set1) > 0:
        result1 = "Multiple tickets contain the same Primary Phage ID."
        status1 = "error"
    else:
        result1 = "Primary PhageID is fine."
        status1 = "correct"

    definition1 = "Check if there are multiple tickets with the same Primary PhageID."
    eval1 = Eval.Eval(id = "TICKET", \
                    definition = definition1, \
                    result = result1, \
                    status = status1)
    result_list.append(eval1)


    dupe_set2 = basic.identify_one_list_duplicates(secondary_id_list)
    if len(dupe_set2) > 0:
        result2 = "Multiple tickets contain the same Secondary Phage ID."
        status2 = "error"
    else:
        result2 = "Secondary PhageID is fine."
        status2 = "correct"

    definition2 = "Check if there are multiple tickets with the same Secondary PhageID."
    eval2 = Eval.Eval(id = "TICKET", \
                    definition = definition2, \
                    result = result2, \
                    status = status2)
    result_list.append(eval2)


    dupe_set3 = basic.identify_one_list_duplicates(accession_list)
    if len(dupe_set3) > 0:
        result3 = "Multiple tickets contain the same Accession."
        status3 = "error"
    else:
        result3 = "Accession is fine."
        status3 = "correct"

    definition3 = "Check if there are multiple tickets with the same Accession."
    eval3 = Eval.Eval(id = "TICKET", \
                    definition = definition3, \
                    result = result3, \
                    status = status3)
    result_list.append(eval3)





    # No 'update' or 'add' Primary Phage IDs are expected to be a
    # 'remove' or 'replace' Secondary Phage IDs.
    dupe_set4 = basic.identify_two_list_duplicates(
                        update_add_primary_id_list, \
                        secondary_id_list)
    if len(dupe_set4) > 0:
        result4 = "An update or add Primary Phage ID is also a " + \
                    "remove or replace Secondary Phage ID."
        status4 = "error"
    else:
        result4 = "There is no PhageID conflict."
        status4 = "correct"

    definition4 = "Check if an update/add Primary PhageID is used " + \
                    "as a remove/replace Secondary PhageID."
    eval4 = Eval.Eval(id = "TICKET", \
                    definition = definition4, \
                    result = result4, \
                    status = status4)
    result_list.append(eval4)


    # No 'replace' Primary Phage IDs are expected to be a 'remove'.
    # Secondary Phage ID.
    dupe_set5 = basic.identify_two_list_duplicates(replace_primary_id_list, \
                        remove_secondary_id_list)
    if len(dupe_set5) > 0:
        result5 = "A replace Primary Phage ID is also a " + \
                    "remove Secondary Phage ID."
        status5 = "error"
    else:
        result5 = "There is no PhageID conflict."
        status5 = "correct"

    definition5 = "Check if a replace Primary PhageID is used " + \
                    "as a remove Secondary PhageID."
    eval5 = Eval.Eval(id = "TICKET", \
                    definition = definition5, \
                    result = result5, \
                    status = status5)
    result_list.append(eval5)



    # Only return results that are Eval objects.
    # list_of_evals = []
    # for eval in result_list:
    #     if eval is not None:
    #         list_of_evals.append(eval)

    # TODO this used to return only error evals, but now it returns all
    # evals, so this will need to be corrected in the unit tests.
    return result_list



def match_genomes_to_tickets1(list_of_group_objects, genome_dict, key):
    """Match genome objects to tickets using phage_id.
    All tickets are expected to be matched."""

    list_of_evals = []

    index = 0
    while index < len(list_of_group_objects):

        group_obj = list_of_group_objects[index]
        phage_id = group_obj.ticket.primary_phage_id

        if phage_id in genome_dict.keys():
            matched_genome = genome_dict[phage_id]
            group_obj.genomes_dict[key] = matched_genome
        else:
            message = "Ticket was not able to be matched with a genome."
            eval_object = Eval.construct_error(message)
            list_of_evals.append(eval_object)

        index += 1

    return list_of_evals

def assign_match_strategy(list_of_group_objects):

    strategy = ""
    strategies = set()
    eval_result = None

    for group_obj in list_of_group_objects:
        strategies.add(group_obj.ticket.match_strategy)


    if len(strategies) > 1:
        message = "Unable to match genomes in files to tickets, " + \
                    "since more than one match strategy is indicated."
        eval_result = Eval.construct_error(message)

    else:
        strategy = list(strategies)[0]

    return strategy, eval_result


# TODO this function may be able to be combined with match_genomes_to_tickets1.
def match_genomes_to_tickets2(list_of_group_objects, all_flat_file_data, key):
    """Match genome objects parsed from files to tickets."""

    # Determine what the strategy is to match tickets to flat files.
    strategy, strategy_eval = assign_match_strategy(list_of_group_objects)

    # Create list of all ticket identifiers.
    ticket_id_list = []
    for group_obj in list_of_group_objects:
        ticket_id_list.append(group_obj.ticket.primary_phage_id)


    # Create list of all genome identifiers based on the strategy.
    genome_dict = {}
    genome_id_list = []
    index1 = 0
    while index1 < len(all_flat_file_data):
        genome_object = all_flat_file_data[index1]

        if strategy == "phage_id":
            genome_id = genome_object.phage_id
        elif strategy == "filename":
            genome_id = genome_object.filename
        else:
            genome_id = ""

        genome_id_list.append(genome_id)
        genome_dict[genome_id] = genome_object

        index1 += 1


    matched_unique_ids, \
    ticket_unmatched_unique_ids, \
    genome_unmatched_unique_ids, \
    ticket_duplicate_ids, \
    genome_duplicate_ids = \
        basic.match_items(ticket_id_list, genome_id_list)


    # Match genomes to tickets using the unique identifiers.
    index2 = 0
    while index2 < len(list_of_group_objects):

        group_obj = list_of_group_objects[index2]
        match_id = group_obj.ticket.primary_phage_id

        if match_id in matched_unique_ids:
            genome_object = genome_dict.pop(match_id)
            group_obj.genomes_dict[key] = genome_object

        index2 += 1



    # Check for the various errors that could have been encountered.
    eval_results = []

    if strategy_eval is not None:
        eval_results.append(strategy_eval)

    if len(ticket_duplicate_ids) > 0:
        message = "Unable to match genomes to tickets " + \
                "since there are multiple tickets with the same identifier."
        eval_result = Eval.construct_error(message)
        eval_results.append(eval_result)

    if len(genome_duplicate_ids) > 0:
        message = "Unable to match genomes to tickets " + \
                "since there are multiple genomes with the same identifier."
        eval_result = Eval.construct_error(message)
        eval_results.append(eval_result)

    if len(ticket_unmatched_unique_ids) > 0:
        message = "There is no matching genome for one or more tickets."
        eval_result = Eval.construct_error(message)
        eval_results.append(eval_result)

    if len(genome_unmatched_unique_ids) > 0:
        message = "There is no matching ticket for one or more genomes."
        eval_result = Eval.construct_error(message)
        eval_results.append(eval_result)

    return eval_results


def create_matched_object_dict(list_of_group_objects):
    """Create a dictionary of DataGroup objects based on their ticket type.
    Key = ticket type (e.g. update, add, etc.).
    Value = list of DataGroup objects."""

    type_set = set()
    for matched_object in list_of_group_objects:
        type_set.add(matched_object.ticket.type)

    ticket_type_dict = {}
    for type in type_set:
        group_object_list = []
        index = 0
        while index < len(list_of_group_objects):
            matched_object = list_of_group_objects[index]
            if matched_object.ticket.type == type:
                group_object_list.append(matched_object)
            index += 1
        ticket_type_dict[type] = group_object_list

    return ticket_type_dict


def complete_ticket(ticket):
    """If the ticket has fields that are set to be auto-completed,
    retrieve the data from PhagesDB to complete the ticket."""

    eval_list1 = []
    if (ticket.host == "retrieve" or \
        ticket.cluster == "retrieve" or \
        ticket.subcluster == "retrieve" or \
        ticket.accession == "retrieve"):

        genome = Genome.Genome()

        phage_url = phagesdb.construct_phage_url(ticket.primary_phage_id)

        data_dict, eval_object1 = phagesdb.retrieve_phagesdb_data(phage_url)

        if eval_object1 is not None:
            eval_list1 += [eval_object1]

        eval_list2 = phagesdb.parse_phagesdb_data(genome, data_dict)

        eval_list1 += eval_list2

        if ticket.host == "retrieve":
            ticket.host = genome.host
        if ticket.cluster == "retrieve":
            ticket.cluster = genome.cluster
        if ticket.subcluster == "retrieve":
            ticket.subcluster = genome.subcluster
        if ticket.accession == "retrieve":
            ticket.accession = genome.accession

    return eval_list1







# TODO implement.
# TODO unit test.


def set_ticket_data(genome_obj, ticket_obj):
    """Populate attributes in a Genome object using data from a ticket."""

    genome_obj.set_phage_id(ticket.primary_phage_id)
    genome_obj.phage_name = ticket.primary_phage_id
    genome_obj.set_host(ticket.host)
    genome_obj.set_accession(ticket.accession)
    genome_obj.author = "" # TODO do I need this in addition to annotation_author?
    genome_obj.status = ticket.status
    genome_obj.set_cluster(ticket.cluster)
    genome_obj.set_subcluster(ticket.subcluster)
    genome_obj.set_cluster_subcluster()
    genome_obj.ncbi_update_flag = ""
    genome_obj.annotation_author = ticket.annotation_author

    # TODO decide how to set this attribute. Check old script.
    if ticket.status == "final":
        genome_obj.annotation_qc = 1
    else:
        genome_obj.annotation_qc = 0

    # TODO decide how to set this attribute. Check old script.
    if ticket.run_mode == "auto_update":
        genome_obj.retrieve_record = 1
    else:
        genome_obj.retrieve_record = 0












# TODO implement below.
# TODO unit test below.




def prepare_tickets(ticket_filename):
    """Parse import table into ticket objects."""

    # Assumes that filename has already been validated.

    # TODO retrieve import data.
    # List of ticket data.


    #Retrieve import info from indicated import table file and read all lines into a list and verify contents are correctly populated.
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



    # TODO rename this function to set_case?
    # Verify all data is cased appropriately.
    for ticket_list in ticket_list:
    	ticket.check_case()


    # Some data may need to be retrieved from PhagesDB.
    for ticket in ticket_list:
        ticket, eval_object = retrieve_online_data(ticket)







    return ticket_list


###















###
