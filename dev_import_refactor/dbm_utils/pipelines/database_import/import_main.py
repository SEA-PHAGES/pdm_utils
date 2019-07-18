"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""











def main(lists_of_ticket_data, list_of_flat_file_data, sql_obj = None):
    """The 'lists_of_ticket_data' parameter is a list, where each
    element is a list of ticket data.
    The 'list_of_flat_file_data' parameter is a list, where each
    element is data parsed from a flat_file into a Biopython SeqRecord.
    The 'sql_obj' parameter handles MySQL connections."""






    #return list_of_matched_objects
    pass








def process_matched_data_object(matched_data_obj,
                                sql_obj,
                                host_genera_set = set(),
                                phage_id_set = set(),
                                seq_set = set(),
                                cluster_set = set(),
                                subcluster_set = set()):
    """Evaluate data within a single DataGroup object."""


    # Create sets of unique values for different data fields.
    #phamerator_data_sets = phamerator.create_data_sets(phamerator_genome_dict)
    phage_id_set = phamerator.create_phage_id_set(sql_obj)
    seq_set = phamerator.create_seq_set(sql_obj)


    # Retrieve data from Phamerator.
    if matched_obj.ticket.type == "replace":


        # TODO retrieve data from phamerator.
        misc.match_genome_by_phage_id(matched_obj,
                        phamerator_genome_dict,
                        "add")
        phamerator.copy_data_from_phamerator(matched_obj, "add")





    # TODO not needed to retrieve all data from phamerator at once.
    # # Retrieve all data from Phamerator
    # # All data is returned as a list of tuples.
    # retrieved_phamerator_data = phamerator.retrieve_sql_data(sql_obj)
    #
    # # Now iterate through the phamerator genome dictionary and create
    # # a second dictionary of parsed Genome objects.
    #
    # # Create a dictionary of all data retrieved.
    # # Key = PhageID.
    # # Value = Genome object with parsed data.
    # phamerator_data_dict = \
    #     phamerator.create_phamerator_dict(retrieved_phamerator_data)
    #
    #
    #
    # # Now that Phamerator data has been retrieved and
    # # Phamerator genome objects created, match them to ticket data
    # index5 = 0
    # while index5 < len(list_of_matched_objects):
    #     matched_obj = list_of_matched_objects[index5]
    #     if matched_obj.ticket.type == "replace":
    #         misc.match_genome_by_phage_id(matched_obj,
    #                         phamerator_genome_dict,
    #                         "add")
    #         phamerator.copy_data_from_phamerator(matched_obj, "add")
    #
    #     index5 += 1

    evaluate.check_datagroup_for_import(matched_obj)
    matched_obj.check_for_errors()







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
