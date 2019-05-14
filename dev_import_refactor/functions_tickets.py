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




def identify_duplicate_primary_ids(list_of_tickets):

    id_list = []
    for ticket in list_of_tickets:

        if ticket.primary_phage_id != "none":
            id_list.append(ticket.primary_phage_id)

    id_set = set(id_list)
    if len(id_set) != len(id_list):
        message = "Multiple tickets contain the same Primary Phage ID."
        eval_object = Eval.construct_error(message)
    else:
        eval_object = None
    return eval_object





def identify_duplicate_secondary_ids(list_of_tickets):

    id_list = []
    for ticket in list_of_tickets:

        if ticket.secondary_phage_id != "none":
            id_list.append(ticket.secondary_phage_id)

    id_set = set(id_list)
    if len(id_set) != len(id_list):
        message = "Multiple tickets contain the same Secondary Phage ID."
        eval_object = Eval.construct_error(message)
    else:
        eval_object = None
    return eval_object



#TODO
# create third function to check for duplicated accessions


#TODO create a fourth function to check for Primary-Secondary duplicated IDs 


#TODO the 'identify_duplicate...' functions above can replace segments of
# code below.

def validate_tickets(list_of_tickets):
    """Verifies there are no ticket conflicts."""

    primary_phage_id_list = []
    secondary_phage_id_list = []
    type_primary_id_list = []
    type_secondary_id_list = []
    type_primary_id_secondary_id_lsit = []
    accession_list = []


    list_of_evals = []
    # Create each tuple and do initial checks for duplications.
    # If the Add name or Remove name is "none", skip that because there are
    # probably duplicates of those.
    for ticket in list_of_tickets:

        if ticket.primary_phage_id != "none":
            primary_phage_id_list.append(ticket.primary_phage_id)
            type_primary_id_list.append((ticket.type,ticket.primary_phage_id))


        if ticket.secondary_phage_id != "none":
            secondary_phage_id_list.append(ticket.secondary_phage_id)
            type_secondary_id_list.append((ticket.type,ticket.secondary_phage_id))


        type_primary_id_secondary_id_lsit.append(\
                (ticket.type,ticket.primary_phage_id,ticket.secondary_phage_id))

        if ticket.accession != "none":
            accession_list.append(ticket.accession)

    # Now, iterate through each list to identify duplicate or
    # conflicting tickets.
    primary_id_set = set(primary_phage_id_list)
    if len(primary_id_set) != len(primary_phage_id_list):
        msg1 = "Multiple tickets contain the same Primary Phage ID."
        eval_object = Eval.construct_error(msg1)
        list_of_evals.append(eval_object)

    secondary_id_set = set(secondary_phage_id_list)
    if len(secondary_id_set) != len(secondary_phage_id_list):
        msg2 = "Multiple tickets contain the same Secondary Phage ID."
        eval_object = Eval.construct_error(msg2)
        list_of_evals.append(eval_object)

    accession_set = set(accession_list)
    if len(accession_set) != len(accession_list):
        msg3 = "Multiple tickets contain the same Accession."
        eval_object = Eval.construct_error(msg3)
        list_of_evals.append(eval_object)


    #Once the sets are created, also check if genomes to be removed are
    #found in the Add field and vice versa.

    for ticket in list_of_tickets:

        current_add = (ticket.primary_phage_id,)
        current_remove = (ticket.secondary_phage_id,)


        #If the phage name is not replacing itself, the Add name is not expected
        #to be in the Remove set and vice versa.


        if ticket.primary_phage_id != ticket.secondary_phage_id:


            if (ticket.primary_phage_id in secondary_id_set or \
                ticket.secondary_phage_id in primary_id_set):

                msg4 = "For some tickets, the Primary Phage ID " + \
                        "is also present as a Secondary Phage ID."
                eval_object = Eval.construct_error(msg4)
                list_of_evals.append(eval_object)


    return list_of_evals


#
# ###
#
#     # Initialize all calculated values.
#     add_set = set()
#     remove_set = set()
#     action_add_set = set()
#     action_remove_set = set()
#     action_add_remove_set = set()
#     import_accession_set = set()
#
#     # Create each tuple and do initial checks for duplications.
#     # If the Add name or Remove name is "none", skip that because there are
#     # probably duplicates of those.
#     for ticket in list_of_tickets:
#     	current_add = (ticket.primary_phage_id,)
#     	current_remove = (ticket.secondary_phage_id,)
#     	current_action_add = (ticket.type,ticket.primary_phage_id)
#     	current_action_remove = (ticket.type,ticket.secondary_phage_id)
#     	current_action_add_remove = (ticket.type,ticket.primary_phage_id,ticket.secondary_phage_id)
#
#
#     	#First check the one-field and two-field combinations
#     	if ticket.primary_phage_id != "none":
#     		if current_add in add_set:
#
#                 msg1 = "The Primary Phage ID %s appears to be involved in more than one step."
#                     % ticket.primary_phage_id
#                 eval_object = Eval.construct_warning(msg1,msg1)
#
#     		else:
#     			add_set.add(current_add)
#
#     		if current_action_add in action_add_set:
#
#                 msg2 = "The ticket type %s and Primary Phage ID %s is duplicated."
#                     % (ticket.type, ticket.primary_phage_id)
#                 eval_object = Eval.construct_error(msg2)
#
#     		else:
#     			action_add_set.add(current_action_add)
#
#     	if current_remove[0] != "none":
#
#     		if current_remove in remove_set:
#
#                 #TODO error handling
#     			print ticket.secondary_phage_id + " appears to be involved in more than one step."
#     			table_errors += question("\nError: %s is duplicated" % str(current_remove))
#
#     		else:
#     			remove_set.add(current_remove)
#
#     		if current_action_remove in action_remove_set:
#
#                 #TODO error handling
#     			write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
#     			table_errors += 1
#
#     		else:
#     			action_remove_set.add(current_action_remove)
#
#     	#Now check the three-field combinations
#     	if current_action_add_remove in action_add_remove_set:
#
#             #TODO error handling
#     		write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
#     		table_errors += 1
#
#     	else:
#     		action_add_remove_set.add(current_action_add_remove)
#
#
#     #Once the sets are created, also check if genomes to be removed are
#     #found in the Add field and vice versa.
#     for ticket in list_of_tickets:
#     	current_add = (ticket.primary_phage_id,)
#     	current_remove = (ticket.secondary_phage_id,)
#
#     	#If the phage name is not replacing itself, the Add name is not expected
#     	#to be in the Remove set and vice versa.
#     	if current_add != current_remove:
#     		if (current_add in remove_set and current_add != "none"):
#
#                 #TODO error handling
#     			print ticket.primary_phage_id + " appears to be involved in more than one step."
#     			table_errors += question("\nError: %s is duplicated" % str(current_add))
#
#
#     		if (current_remove in add_set and current_remove != "none"):
#
#                 #TODO error handling
#     			print ticket.secondary_phage_id + " appears to be involved in more than one step."
#     			table_errors += question("\nError: %s is duplicated" % str(current_remove))
#
#
#
#
#     #Verify there are no duplicate accessions in the import table
#     for ticket in list_of_tickets:
#
#     	if ticket.accession != "none":
#     		if ticket.accession in import_accession_set:
#
#                 #TODO error handling
#     			write_out(output_file,"\nError: Accession %s is duplicated in the import table." %ticket.accession)
#     			table_errors += 1
#     		else:
#     			import_accession_set.add(ticket.accession)
#
#
#
#
# 	#TODO return proper information
# 	return pass



#
#
#
# # For each ticket, retrieve any data from PhagesDB genome if necessary.
# def retrieve_online_data(ticket_list):
#
#     for ticket in ticket_list:
#
#         if (ticket.host == "retrieve" or \
#             ticket.cluster == "retrieve" or \
#             ticket.subcluster == "retrieve" or \
#             ticket.accession == "retrieve"):
#
#             phagesdb_genome = Genome()
#
#             #TODO make sure api_prefix and api_suffix variables are set
#             phage_url = api_prefix + ticket.primary_phage_id + api_suffix
#
#
#             try:
#                 online_data_json = urllib.urlopen(phage_url)
#                 online_data_dict = json.loads(online_data_json.read())
#
#                 #Returns a genome object
#                 phagesdb_genome = parse_phagesdb_data(phagesdb_genome,online_data_dict)
#
#             except:
#                 online_data_json = ""
#                 online_data_dict = {}
#
#                 #TODO handle error better
#                 write_out(output_file,"\nError: unable to retrieve Host, Cluster, Subcluster, or Accession data for phage %s from phagesdb." %row[1])
#
#             if ticket.host == "retrieve":
#                 ticket.host = phagesdb_genome.host
#             if ticket.cluster == "retrieve":
#                 ticket.cluster = phagesdb_genome.cluster
#             if ticket.subcluster == "retrieve":
#                 ticket.subcluster = phagesdb_genome.subcluster
#             if ticket.accession == "retrieve":
#                 ticket.accession = phagesdb_genome.accession
