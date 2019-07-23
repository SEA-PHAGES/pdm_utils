"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""



from functions import tickets
from functions import flat_files
from functions import phagesdb
from functions import phamerator







def main1(lists_of_ticket_data, list_of_flat_file_data, sql_obj = None):
    """The 'lists_of_ticket_data' parameter is a list, where each
    element is a list of ticket data.
    The 'list_of_flat_file_data' parameter is a list, where each
    element is data parsed from a flat_file into a Biopython SeqRecord.
    The 'sql_obj' parameter handles MySQL connections."""





    # Convert to ticket objects.
    # Data is returned as a list of ticket objects.
    # lists_of_ticket_data = read.csv(ticket_filename)
    ticket_list = tickets.parse_import_tickets(lists_of_ticket_data)



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



    # Create a dictionary of tickets based on the primary_phage_id.
    ticket_dict = {}
    index2 = 0
    while index2 < len(ticket_list):
        ticket_dict[ticket_list[index2].primary_phage_id] = ticket_list[index2]




    # TODO check for ticket errors = exit script if not structured correctly?






    # If there is at least one file to process, retrieve data from phagesdb.
    if len(files_in_folder) > 0:
        phagesdb_host_genera_set = phagesdb.create_host_genus_set()

        phagesdb_cluster_set, \
        phagesdb_subcluster_set = phagesdb.create_cluster_subcluster_sets()

    else:
        phagesdb_host_genera_set = set()
        phagesdb_cluster_set = set()
        phagesdb_subcluster_set = set()



    # To minimize memory usage, each flat_file is evaluated one by one.

    for filename in files_in_folder:


        genome = fasta_files.create_parsed_flat_file(filename)
        bundle = Bundle.Bundle()
        bundle.genome_dict[genome.type] = genome

        # Match ticket (if available) to flat file.
        matched_ticket = ticket_dict.pop(genome.id, None)
        bundle.ticket = matched_ticket



        # Create sets of unique values for different data fields.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_obj)
        phamerator_seq_set = phamerator.create_seq_set(sql_obj)

        # Perform all evaluations based on the ticket type.
        import_main.main2(bundle,
                            sql_obj,
                            phage_id_set = phamerator_phage_id_set,
                            seq_set = phamerator_seq_set,
                            host_genera_set = phagesdb_host_genera_set,
                            cluster_set = phagesdb_cluster_set,
                            subcluster_set = phagesdb_subcluster_set)


        # Now import the data into the database.
        import_main.import_into_database(bundle, sql_obj)



    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    for key in ticket_dict.keys():

        unmatched_ticket = ticket_dict.pop(key)
        bundle = Bundle.Bundle()
        bundle.ticket = unmatched_ticket

        # TODO run a list of evaluations based on ticket type.





    #return list_of_bundles





def main2(bundle, sql_obj, host_genera_set = set(), phage_id_set = set(),
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
            phamerator.create_phamerator_genome(sql_obj, genome.id)


        # If any attributes in flat_file are set to 'retain', copy data
        # from the phamerator genome.
        phamerator.copy_data_from_phamerator(bundle, "flat_file")







    evaluate.check_bundle_for_import(bundle)
    bundle.check_for_errors()







def import_into_database(list_of_bundles, sql_obj):
    """Construct collection of SQL statements for evaluated data and
    execute statements to add the data into the database."""



    # Create all SQL statements for all tickets with no errors.
    index11 = 0
    while index11 < len(list_of_bundles):
            bundle = list_of_bundles[index11]
            if bundle.errors == 0:
                bundle.create_sql_statements()




    # TODO import all scrubbed add_replace data into Phamerator.



    pass






















###
