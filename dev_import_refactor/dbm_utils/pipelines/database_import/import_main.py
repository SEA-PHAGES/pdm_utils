"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""











def main(list_of_ticket_data, list_of_flat_file_data, sql_obj = None):
    """The 'list_of_ticket_data' parameter is a list, where each
    element is a list of ticket data.
    The 'list_of_flat_file_data' parameter is a list, where each
    element is data parsed from a flat_file into a Biopython SeqRecord.
    The 'sql_obj' parameter handles MySQL connections."""



    # TODO parsing from import table:
    # 1. parse ticket data from table. = prepare_tickets()
    # 2. set case for all fields. = prepare_tickets()
    # 3. confirm all tickets have a valid type. = check_ticket_structure()
    # 4. populate Genome objects as necessary.
    # 5. retrieve data if needed.
    # 6. check for PhageID conflicts.
    # 7. confirm correct fields are populated based on ticket type.


    # Retrieve import ticket data.
    # Data is returned as a list of validated ticket objects.
    # list_of_ticket_data = read.csv(ticket_filename)
    list_of_tickets = tickets.parse_import_tickets(list_of_ticket_data)



    # Evaluate the tickets to ensure they are structured properly.


    # TODO not sure if I should pass a list of valid types to this function.
    index1 = 0
    while index1 < len(ticket_list):
        evaluate.check_ticket_structure(ticket_list[index1],
                                        constants.TICKET_TYPE_SET,
                                        constants.EMPTY_SET,
                                        constants.RUN_MODE_SET)
        index1 += 1


    # Now that individual tickets have been validated,
    # validate the entire group of tickets.
    tickets.compare_tickets(ticket_list)


    # Tickets will be matched with other genome data.
    # Ticket data will be paired with data from PhagesDB, PhameratorDB,
    # and/or a flat file.
    list_of_matched_objects = []
    index2 = 0
    while index2 < len(list_of_tickets):

        matched_data_obj = MatchedGenomes()
        matched_data_obj.ticket = list_of_tickets[index2]
        list_of_matched_objects.append(matched_data_obj)
        index2 += 1

    # Using each ticket, construct and populate genome objects as needed.
    index3 = 0
    while index3 < len(list_of_matched_objects):
        tickets.copy_ticket_to_genome(list_of_matched_objects[index3])
        index3 += 1


    # Now check to see if there is any missing data for each genome, and
    # retrieve it from phagesdb.
    index4 = 0
    while index4 < len(list_of_matched_objects):

        # If the ticket genome has fields set to 'retrieve', data is
        # retrieved from PhagesDB and populates a new Genome object.
        matched_object = list_of_matched_objects[index4]
        phagesdb.copy_data_from_phagesdb(matched_object, "add")
        index4 += 1




    # Each ticket should now contain all requisite data from PhagesDB.
    # Validate each ticket by checking each field in the ticket
    # that it is populated correctly.



    # TODO check for ticket errors = exit script if not structured correctly.
    # Iterate through tickets and collect all evals.



    if len(list_of_errors) > 0:
        sys.exit(1)






    # TODO at some point annotation_qc and retrieve_record attributes
    # will need to be set. These are dependent on the ticket type.
    # If genomes are being replace, these fields may be carried over from
    # the previous genome, combined with their annotation status.



    # Retrieve data from Phamerator.

    # TODO create SQL connector object using parsed arguments.
    # TODO it may be better to create the SQL object with the
    # prepare_phamerator_data module.


    # Retrieve all data from Phamerator
    # All data is returned as a list of tuples.
    retrieved_phamerator_data = phamerator.retrieve_sql_data(sql_obj)


    # Now iterate through the phamerator genome dictionary and create
    # a second dictionary of parsed Genome objects.

    # Create a dictionary of all data retrieved.
    # Key = PhageID.
    # Value = Genome object with parsed data.
    phamerator_data_dict = \
        phamerator.create_phamerator_dict(retrieved_phamerator_data)



    # Now that Phamerator data has been retrieved and
    # Phamerator genome objects created, match them to ticket data
    index5 = 0
    while index5 < len(list_of_matched_objects):
        matched_obj = list_of_matched_objects[index5]
        if matched_obj.ticket.type == "replace":
            misc.match_genome_by_phage_id(matched_obj,
                            phamerator_genome_dict,
                            "add")
            phamerator.copy_data_from_phamerator(matched_obj, "add")

        index5 += 1












    # TODO when should this be implemented?
    # Create sets of unique values for different data fields.
    phamerator_data_sets = phamerator.create_data_sets(phamerator_genome_dict)








    # TODO check for phamerator data errors = exit script if there are errors
    # TODO is this needed?
    if len(phamerator_errors) > 0:
        sys.exit(1)













    # Match tickets to flat file data





    # TODO check to confirm that genome objects from parsed flat files do not
    # contain any duplicate phage_ids, since that is not gauranteed.


    # Create dictionary of flat file data.
    flat_file_dict = {}
    index100 = 0
    while index100 < len(flat_file_genomes):
        genome = flat_file_genomes[index100]
        flat_file_dict[genome.phage_id] = genome



    # Match flat file genomes.
    index6 = 0
    while index6 < len(list_of_matched_objects):
        matched_obj = list_of_matched_objects[index6]
        misc.match_genome_by_phage_id(matched_obj, flat_file_dict, "add")






    # TODO don't think this is needed anymore, since 'update' and 'remove'
    # tickets aren't valid types for this script.
    # TODO Now that all data has been matched, split matched objects by ticket type.
    # Different types of tickets are evaluated differently.
    # matched_object_dict = create_matched_object_dict(list_of_matched_objects)
    #
    # list_of_update_objects = matched_object_dict["update"]
    # list_of_remove_objects = matched_object_dict["remove"]
    # list_of_add_replace_objects = matched_object_dict["add_replace"]





    # TODO After parsing flat file, prepare gene_id and gene_name appropriately



    # Now that the flat file to be imported is parsed and matched to a ticket,
    # use the ticket to populate specific genome-level fields such as
    # host, cluster, subcluster, etc.
    index7 = 0
    while index7 < len(list_of_add_replace_objects):
        matched_object = list_of_add_replace_objects[index7]
        flat_files.copy_data_to_flat_file(matched_object, "add")
        index7 += 1







    # Perform all evaluations based on the ticket type.






    # TODO after each add_replace ticket is evaluated,
    # should the script re-query the database and re-create the
    # sets of PhageIDs, Sequences, etc?
    index10 = 0
    while index10 < len(list_of_matched_objects):
            evaluate.check_datagroup_for_import(list_of_matched_objects[index10])
            list_of_matched_objects[index10].check_for_errors()












    return list_of_matched_objects








def process_matched_data_object(matched_data_obj, sql_obj):
    """Evaluate data within a single DataGroup object."""
    pass






def import_into_database(list_of_matched_objects, sql_obj):
    """Construct collection of SQL statements for evaluated data and
    execute statements to add the data into the database."""



    # Create all SQL statements for all tickets with no errors.
    index11 = 0
    while index11 < len(list_of_matched_objects):
            matched_object = list_of_matched_objects[index11]
            if matched_object.errors == 0:
                matched_object.create_sql_statements()




    # TODO import all scrubbed add_replace data into Phamerator.



    pass






















###
