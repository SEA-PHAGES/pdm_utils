"""Misc. functions to manipulate tickets."""


from classes import Eval
from classes import Ticket





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

        ticket = Ticket.ImportTicket()
        ticket.type = data_list[0]
        ticket.primary_phage_id = data_list[1]
        ticket.host = data_list[2]
        ticket.cluster = data_list[3]
        ticket.subcluster = data_list[4]
        ticket.status = data_list[5]
        ticket.annotation_author = data_list[6]
        ticket.description_field = data_list[7]
        ticket.accession = data_list[8]
        ticket.run_mode = data_list[9]
        ticket.secondary_phage_id = data_list[10]

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




#TODO this function should probably be moved to the general functions module.
def identify_duplicates1(item_list, message):
    """Verifies there are no duplicate items in a list."""

    item_set = set(item_list)
    if len(item_set) != len(item_list):
        eval_object = Eval.construct_error(message)
    else:
        eval_object = None
    return eval_object



#TODO this function should probably be moved to the general functions module.
def identify_duplicates2(item1_list, item2_list, message):
    """Verifies there are no duplicate items between two lists."""

    item1_set = set(item1_list)
    item2_set = set(item2_list)
    item3_list = item1_set & item2_set

    if len(item3_list) > 0:
        eval_object = Eval.construct_error(message)
    else:
        eval_object = None
    return eval_object




def validate_tickets(list_of_tickets):
    """Verifies there are no ticket conflicts."""

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
    result1 = identify_duplicates1(primary_id_list,\
        "Multiple tickets contain the same Primary Phage ID.")
    result_list.append(result1)

    result2 = identify_duplicates1(secondary_id_list,\
        "Multiple tickets contain the same Secondary Phage ID.")
    result_list.append(result2)

    result3 = identify_duplicates1(accession_list,\
        "Multiple tickets contain the same Accession.")
    result_list.append(result3)



    # No 'update' or 'add' Primary Phage IDs are expected to be a
    # 'remove' or 'replace' Secondary Phage IDs.
    result4 = identify_duplicates2(update_add_primary_id_list, \
                        secondary_id_list, \
                        "An update or add Primary Phage ID is also a " + \
                        "remove or replace Secondary Phage ID.")
    result_list.append(result4)


    # No 'replace' Primary Phage IDs are expected to be a 'remove'.
    # Secondary Phage ID.
    result5 = identify_duplicates2(replace_primary_id_list, \
                        remove_secondary_id_list, \
                        "A replace Primary Phage ID is also a " + \
                        "remove Secondary Phage ID.")
    result_list.append(result5)


    # Only return results that are Eval objects.
    list_of_evals = []
    for eval in result_list:
        if eval is not None:
            list_of_evals.append(eval)


    return list_of_evals



def match_genomes_to_tickets(list_of_group_objects, genome_dict, key):
    """Match genome objects to tickets using phage_id.
    All tickets are expected to be matched."""

    list_of_evals = []

    for group_obj in list_of_group_objects:
        phage_id = group_obj.ticket.primary_phage_id

        if phage_id in genome_dict.keys():
            matched_genome = genome_dict[phage_id]
            group_obj.genomes_dict[key] = matched_genome

        else:
            message = "Ticket was not able to be matched with a genome."
            eval_object = Eval.construct_error(message)

            list_of_evals.append(eval_object)

    return list_of_evals






#TODO unit test below



#TODO unit test
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




#TODO unit test
def match_files_to_tickets(list_of_datagroup_objects, all_flat_file_data, key):

    eval_results = []
    # Iterate through all tickets and determine what the strategy is to
    # match tickets to flat files.
    strategy, strategy_eval = assign_match_strategy(list_of_datagroup_objects)

    if strategy_eval not None:
        eval_results.append(strategy_eval)


    # TODO use FunctionsSimple.match_items() function to identify whether
    # tickets can be matched to flat file genomes instead of using the text
    # below.

    # Confirm all ticket primary_ids are unique.
    ticket_primary_ids = set()
    for group_obj in list_of_datagroup_objects:

        ticket = group_obj.ticket
        primary_id = ticket.primary_phage_id

        if primary_id not in ticket_primary_ids:
            ticket_primary_ids.add(primary_id)

        else:
            message = "Unable to match genomes to tickets " + \
                    "since there are multiple tickets with the same identifier."
            eval_result = Eval.construct_error(message)
            eval_results.append(eval_result)


    # Confirm that all genomes contain unique matching identifiers.
    flat_file_genome_dict = {}
    for flat_file_object in all_flat_file_data:

        if strategy == "phage_id":
            match_name = flat_file_object.phage_id

        elif strategy == "filename":
            match_name = flat_file_object.filename

        else:
            match_name = ""

        if match_name not in flat_file_genome_dict.keys():
            flat_file_genome_dict[match_name] = flat_file_object
        else:
            message = "Unable to match genome to a ticket " + \
                    "since there is another genome with the same identifier."
            eval_result = Eval.construct_error(message)
            eval_results.append(eval_result)




    # Now match tickets to flat file genomes.
    for group_obj in matched_data_list:
        match_name = group_obj.ticket.primary_phage_id

        try:
            flat_file_genome = flat_file_genome_dict.pop(match_name)
            group_obj.genomes_dict[key] = flat_file_genome


        except:
            message = "There is no matching genome for this ticket."
            eval_result = Eval.construct_error(message)
            eval_results.append(eval_result)


    return eval_results




#TODO unit test
def create_matched_object_dict(list_of_datagroup_objects):

    dictionary = {}
    list_of_update_objects = []
    list_of_remove_objects = []
    list_of_add_replace_objects = []

    for matched_object in list_of_datagroup_objects:
        type = matched_object.ticket.ticket_type

        if type == "update":
            list_of_update_objects.append(matched_object)

        elif type == "remove":
            list_of_remove_objects.append(matched_object)

        elif (type == "add" or type == "replace"):
            list_of_add_replace_objects.append(matched_object)

        else:
            pass
            # TODO if ticket type is none of the above, thrown an error?

    dictionary["update"] = list_of_update_objects
    dictionary["remove"] = list_of_remove_objects
    dictionary["add_replace"] = list_of_add_replace_objects

    return dictionary





#TODO unit test after unit tested nested functions.
# For each ticket, retrieve any data from PhagesDB genome if necessary.
def retrieve_online_data(ticket):

    if (ticket.host == "retrieve" or \
        ticket.cluster == "retrieve" or \
        ticket.subcluster == "retrieve" or \
        ticket.accession == "retrieve"):

        phagesdb_genome = Genome()

        #TODO make sure api_prefix and api_suffix variables are set
        phage_url = api_prefix + ticket.primary_phage_id + api_suffix


        try:
            online_data_json = urllib.urlopen(phage_url)
            online_data_dict = json.loads(online_data_json.read())

            #Returns a genome object
            phagesdb_genome = \
                parse_phagesdb_data(phagesdb_genome, online_data_dict)

            if ticket.host == "retrieve":
                ticket.host = phagesdb_genome.host
            if ticket.cluster == "retrieve":
                ticket.cluster = phagesdb_genome.cluster
            if ticket.subcluster == "retrieve":
                ticket.subcluster = phagesdb_genome.subcluster
            if ticket.accession == "retrieve":
                ticket.accession = phagesdb_genome.accession

        except:
            # online_data_json = ""
            # online_data_dict = {}


            message = "Unable to retrieve data for phage %s from phagesdb." \
                        % ticket.primary_phage_id
            eval_object = Eval.construct_error(message)


    # TODO Does the ticket need to be returned as well?
    # TODO return the online dictionary as well
    return eval_object











###
