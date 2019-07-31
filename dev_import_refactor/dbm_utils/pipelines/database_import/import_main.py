"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""



from functions import tickets
from functions import flat_files
from functions import phagesdb
from functions import phamerator
from classes import Bundle
from pipelines.database_import import evaluate
from constants import constants







def main1(lists_of_ticket_data, files_in_folder, sql_handle = None):
    """The 'lists_of_ticket_data' parameter is a list, where each
    element is a list of ticket data.
    The 'files_in_folder' parameter is a list, where each
    element is a filename to be parsed.
    The 'sql_handle' parameter handles MySQL connections."""


    # Will hold all results data. These can be used by various types of
    # user interfaces to present the summary of import.
    success_ticket_list = []
    failed_ticket_list = []
    success_filename_list = []
    failed_filename_list = []
    evaluation_dict = {}
    query_dict = {}

    # Convert ticket data to Ticket objects.
    # Data is returned as a list of ticket objects.
    ticket_list = tickets.parse_import_tickets(lists_of_ticket_data)
    print("Tickets parsed")


    # Evaluate the tickets to ensure they are structured properly.


    # TODO not sure if I should pass a list of valid types to this function.
    index1 = 0
    while index1 < len(ticket_list):
        evaluate.check_ticket_structure(ticket_list[index1],
                                        constants.TICKET_TYPE_SET,
                                        constants.EMPTY_SET,
                                        constants.RUN_MODE_SET)
        index1 += 1
    print("Ticket structure checked")


    # Now that individual tickets have been validated,
    # validate the entire group of tickets.
    tickets.compare_tickets(ticket_list)
    print("Tickets compared")


    # Create a dictionary of tickets based on the primary_phage_id.
    ticket_dict = {}
    index2 = 0
    while index2 < len(ticket_list):
        ticket_dict[ticket_list[index2].primary_phage_id] = ticket_list[index2]
        index2 += 1
    print("Ticket dictionary created")




    # TODO check for ticket errors = exit script if not structured correctly?






    # If there is at least one file to process, retrieve data from phagesdb to
    # create sets of valid host genera, clusters, and subclusters.
    if len(files_in_folder) > 0:
        phagesdb_host_genera_set = phagesdb.create_host_genus_set()

        phagesdb_cluster_set, \
        phagesdb_subcluster_set = phagesdb.create_cluster_subcluster_sets()

    else:
        phagesdb_host_genera_set = set()
        phagesdb_cluster_set = set()
        phagesdb_subcluster_set = set()

    print("PhagesDB sets retrieved")


    # To minimize memory usage, each flat_file is evaluated one by one.

    for filename in files_in_folder:


        genome = flat_files.create_parsed_flat_file(filename)
        bundle = Bundle.Bundle()
        bundle.genome_dict[genome.type] = genome

        # Match ticket (if available) to flat file.
        matched_ticket = ticket_dict.pop(genome.id, None)
        bundle.ticket = matched_ticket



        # Create sets of unique values for different data fields.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)

        # Perform all evaluations based on the ticket type.
        import_main.main2(bundle = bundle,
                            sql_handle = sql_handle,
                            phage_id_set = phamerator_phage_id_set,
                            seq_set = phamerator_seq_set,
                            host_genera_set = phagesdb_host_genera_set,
                            cluster_set = phagesdb_cluster_set,
                            subcluster_set = phagesdb_subcluster_set)



        # Now that all evaluations have been performed,
        # determine if there are any errors.
        eval_sub_dict = bundle.get_evaluations() # TODO implement this method.
        errors = 0
        for key in eval_sub_dict.keys():
            eval_list = eval_sub_dict[key]
            for eval in eval_list:
                if eval.status == "error"
                    errors += 1



        # TODO construct a basic.get_empty_ticket() function.
        empty_ticket = [None, None, None, None, None, None,
                        None, None, None, None, None, None, ]
        ticket_data =  tickets.parse_import_ticket_data(\
                            bundle.ticket, empty_ticket,
                            direction = "ticket_to_list")))

        if errors == 0:


            # Now import the data into the database if there are no errors and
            # if there is MySQL connection data provided.
            if sql_handle is not None:
                bundle.create_sql_statements()

                # TODO confirm the handler method name and that it can
                # handle a list of queries.
                sql_handle.execute(bundle.sql_queries)
                sql_handle.commit()

                # If successful, keep track of query data.
                query_dict[bundle.ticket.primary_phage_id] = bundle.sql_queries
                success_ticket_list.append(ticket_data)
            else:
                success_ticket_list.append(ticket_data)

        else:
            failed_ticket_list.append(ticket_data)


    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    key_list = set(ticket_dict.keys())
    for key in key_list:
        unmatched_ticket = ticket_dict.pop(key)
        bundle = Bundle.Bundle()
        bundle.ticket = unmatched_ticket

        # TODO run a list of evaluations based on ticket type.





    #return list_of_bundles
    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict, query_dict)




def main2(bundle, sql_handle, host_genera_set = set(), phage_id_set = set(),
            seq_set = set(), cluster_set = set(), subcluster_set = set()):
    """Evaluate data within a single Bundle object."""


    # TODO will need to have an evaluation check to see if a ticket object
    # has been assigned or not.


    # Using each ticket, construct and populate genome objects as needed.
    tickets.copy_ticket_to_genome(bundle)



    # Now that the flat file to be imported is parsed and matched to a ticket,
    # use the ticket to populate specific genome-level fields such as
    # host, cluster, subcluster, etc.
    flat_files.copy_data_to_flat_file(bundle, "add")



    # Now check to see if there is any missing data for each genome, and
    # retrieve it from phagesdb.

    # If the ticket genome has fields set to 'retrieve', data is
    # retrieved from PhagesDB and populates a new Genome object.
    phagesdb.copy_data_from_phagesdb(bundle, "flat_file")


    # Each flat_file genome should now contain all requisite data from PhagesDB.
    # Validate each genome by checking that each field is populated correctly.


    # TODO at some point annotation_qc and retrieve_record attributes
    # will need to be set. These are dependent on the ticket type.
    # If genomes are being replace, these fields may be carried over from
    # the previous genome, combined with their annotation status.



    # TODO After parsing flat file, prepare gene_id and gene_name appropriately


    # If the ticket type is 'replace', retrieve data from phamerator.
    # TODO will need to account for whether the phage_id exists in Phamerator or not.
    if matched_ticket.type == "replace":
        phamerator_genome = \
            phamerator.create_phamerator_genome(sql_handle, genome.id)


        # If any attributes in flat_file are set to 'retain', copy data
        # from the phamerator genome.
        phamerator.copy_data_from_phamerator(bundle, "flat_file")







    evaluate.check_bundle_for_import(bundle)
    bundle.check_for_errors()







def import_into_database(bundle, sql_handle):
    """Construct collection of SQL statements for evaluated data and
    execute statements to add the data into the database."""



    # Create all SQL statements.
    index11 = 0
    while index11 < len(list_of_bundles):
            bundle = list_of_bundles[index11]
            if bundle.errors == 0:
                bundle.create_sql_statements()




    # TODO import all scrubbed add_replace data into Phamerator.



    pass






















###
