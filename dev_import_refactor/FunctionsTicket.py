"""Misc. functions to manipulate tickets."""


import Eval
import Ticket





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

    item_set = set(item_list)
    if len(item_set) != len(item_list):
        eval_object = Eval.construct_error(message)
    else:
        eval_object = None
    return eval_object



#TODO this function should probably be moved to the general functions module.
def identify_duplicates2(item1_list, item2_list, message):

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













#TODO unit test below







#TODO unit test
def match_phamerator_to_tickets(list_of_matched_objects, phamerator_genome_dict):

    for matched_data_obj in list_of_matched_objects:
        phage_id = matched_data_obj.ticket.primary_phage_id

        try:
            matched_genome = phamerator_genome_dict[phage_id]

            #Now add the Phamerator data to the MatchedGenomes object
            matched_data_obj.matched_genomes_dict["phamerator"] = matched_genome

        except:
            pass
            #TODO error handling

    return list_of_matched_objects





#TODO unit test
def assign_match_strategy(list_of_matched_objects):

	for matched_object in list_of_matched_objects:
		ticket = matched_object.ticket


		#TODO complete this function
	pass


#TODO unit test
def match_flat_files_to_tickets(list_of_matched_objects, all_flat_file_data):
    flat_file_genome_dict = {}

    # TODO:
    # Iterate through all tickets and determine what the strategy is to
    # match tickets to flat files.
    strategy = assign_match_strategy(list_of_matched_objects)

    #TODO need to set strategy variable in advance
    # Now that the flat file data has been retrieved and parsed,
    # match them to ticket data


    for flat_file_object in list_of_flat_file_genomes:

        if strategy == "phage_id":
            match_name = flat_file_object.phage_id

        elif strategy == "filename":
            match_name = flat_file_object.filename

        else:
            match_name = ""

        if match_name not in flat_file_genome_dict.keys():
            flat_file_genome_dict[match_name] = flat_file_object
        else:
            pass
            #TODO throw an error - unable to create unique set of objects

    for matched_data_obj in matched_data_list:
        match_name = matched_data_obj.ticket.primary_phage_id

        try:
            flat_file_genome = flat_file_genome_dict.pop(match_name)

        except:
            flat_file_genome = None

        #Now add the flat file data to the MatchedGenomes object
        matched_data_obj.matched_genomes_dict["flat_file"] = flat_file_genome

    return list_of_matched_objects




#TODO unit test
def create_matched_object_dict(list_of_matched_objects):

    dictionary = {}
    list_of_update_objects = []
    list_of_remove_objects = []
    list_of_add_replace_objects = []

    for matched_object in list_of_matched_objects:
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
