"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""



from pdm_utils.functions import tickets
from pdm_utils.functions import flat_files
from pdm_utils.functions import phagesdb
from pdm_utils.functions import phamerator
from pdm_utils.classes import bundle
from pdm_utils.pipelines.db_import import evaluate
from pdm_utils.constants import constants





# TODO incorporate this standard phage table query to retrieve all
# genome-level data from PhameratorDB:
# statement1 = "SELECT PhageID, Name, HostStrain, Sequence, status, \
#              Cluster2, DateLastModified, Accession, Subcluster2, \
#              AnnotationAuthor, AnnotationQC, RetrieveRecord \
#              FROM phage"



def main(lists_of_ticket_data, files_in_folder, sql_handle = None):
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
    # print("Tickets parsed")


    # Evaluate the tickets to ensure they are structured properly.
    # At this point, the quality of the ticket data is not evaluated,
    # just that the ticket contains fields populated or empty as expected.
    index1 = 0
    while index1 < len(ticket_list):
        evaluate.check_ticket_structure(ticket_list[index1],
                                        constants.IMPORT_TICKET_TYPE_SET,
                                        constants.EMPTY_DATE,
                                        constants.RUN_MODE_SET),
        index1 += 1
    # print("Ticket structure checked")


    # Now that individual tickets have been validated,
    # validate the entire group of tickets.
    tickets.compare_tickets(ticket_list)
    # print("Tickets compared")


    # Create a dictionary of tickets based on the phage_id.
    ticket_dict = {}
    index2 = 0
    while index2 < len(ticket_list):
        ticket_dict[ticket_list[index2].phage_id] = ticket_list[index2]
        index2 += 1
    # print("Ticket dictionary created")




    # Check for ticket errors.
    ticket_errors = 0
    for key in ticket_dict.keys():
        tkt = ticket_dict[key]
        for evl in tkt.evaluations:
            if evl.status == "error":
                ticket_errors += 1


    # If there are no ticket errors, proceed with evaluating flat files.
    if ticket_errors == 0:

        # If there is at least one file to process, retrieve data from
        # phagesdb to create sets of valid host genera, clusters,
        # and subclusters.
        if (len(files_in_folder) > 0 and ticket_errors == 0):
            phagesdb_host_genera_set = phagesdb.create_host_genus_set()

            phagesdb_cluster_set, \
            phagesdb_subcluster_set = phagesdb.create_cluster_subcluster_sets()

        else:
            phagesdb_host_genera_set = set()
            phagesdb_cluster_set = set()
            phagesdb_subcluster_set = set()

        # print("PhagesDB sets retrieved")


        # To minimize memory usage, each flat_file is evaluated one by one.
        bundle_count = 1
        for filename in files_in_folder:


            gnm = flat_files.create_parsed_flat_file(filename)
            bndl = bundle.Bundle()
            bndl.id = bundle_count
            bndl.genome_dict[gnm.type] = gnm

            # Match ticket (if available) to flat file.
            matched_ticket = ticket_dict.pop(gnm.id, None)
            bndl.ticket = matched_ticket



            # Create sets of unique values for different data fields.
            # Since data from each parsed flat file is imported into the
            # database one file at a time, these sets are not static.
            # So these sets should be recomputed for every flat file evaluated.
            phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
            phamerator_seq_set = phamerator.create_seq_set(sql_handle)


            # TODO implement the main2 function.
            # Perform all evaluations based on the ticket type.
            import_main.evaluate_flat_file(bndl = bndl,
                                    sql_handle = sql_handle,
                                    phage_id_set = phamerator_phage_id_set,
                                    seq_set = phamerator_seq_set,
                                    host_genera_set = phagesdb_host_genera_set,
                                    cluster_set = phagesdb_cluster_set,
                                    subcluster_set = phagesdb_subcluster_set)



            # Now that all evaluations have been performed,
            # determine if there are any errors.




            # TODO construct a basic.get_empty_ticket() function?
            empty_ticket = [None] * 12
            ticket_data =  tickets.parse_import_ticket_data(\
                                bndl.ticket, empty_ticket,\
                                direction = "ticket_to_list")

            bndl.check_for_errors()
            if errors == 0:

                # Now import the data into the database if there are no errors and
                # if there is MySQL connection data provided.
                if sql_handle is not None:
                    bndl.create_sql_statements()

                    # TODO confirm the handler method name and that it can
                    # handle a list of queries.
                    sql_handle.execute(bndl.sql_queries)
                    sql_handle.commit()

                    # If successful, keep track of query data.
                    query_dict[bndl.ticket.phage_id] = bndl.sql_queries
                else:
                    pass
                success_ticket_list.append(ticket_data)
                success_filename_list.append(bndl.genome_dict["add"].filename)

            else:
                failed_ticket_list.append(ticket_data)
                failed_filename_list.append(bndl.genome_dict["add"].filename)

            # TODO implement the get_evaluations() method.
            evaluation_dict[bndl.id] = bndl.get_evaluations()
            bundle_count += 1


    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    key_list = ticket_dict.keys()
    for key in key_list:
        unmatched_ticket = ticket_dict.pop(key)
        bndl = bundle.Bundle()
        bndl.id = bundle_count
        bndl.ticket = unmatched_ticket

        # TODO run a list of evaluations based on ticket type.
        # TODO determine which evaluate function is best.
        evaluate.tickets(unmatched_ticket)
        evaluation_dict[bndl.id] = {"ticket":unmatched_ticket.evaluations}
        bundle_count += 1

        empty_ticket = [None] * 12
        ticket_data =  tickets.parse_import_ticket_data(\
                            unmatched_ticket, empty_ticket,\
                            direction = "ticket_to_list")
        failed_ticket_list.append(ticket_data)



    #return list_of_bundles
    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict, query_dict)




def evaluate_flat_file(bndl, sql_handle, host_genera_set=set(),
                       phage_id_set=set(), seq_set=set(),
                       cluster_set=set(), subcluster_set=set()):
    """Evaluate data within a single Bundle object."""


    # TODO will need to have an evaluation check to see if a ticket object
    # has been assigned or not.


    # Using each ticket, construct and populate genome objects as needed.
    tickets.copy_ticket_to_genome(bndl)



    # Now that the flat file to be imported is parsed and matched to a ticket,
    # use the ticket to populate specific genome-level fields such as
    # host, cluster, subcluster, etc.
    flat_files.copy_data_to(bndl, "add")



    # Now check to see if there is any missing data for each genome, and
    # retrieve it from phagesdb.

    # If the ticket genome has fields set to 'retrieve', data is
    # retrieved from PhagesDB and populates a new Genome object.
    phagesdb.copy_data_from(bndl, "flat_file")





    # Each flat_file genome should now contain all requisite data from PhagesDB.
    # Validate each genome by checking that each field is populated correctly.


    # TODO at some point annotation_qc and retrieve_record attributes
    # will need to be set. These are dependent on the ticket type.
    # If genomes are being replace, these fields may be carried over from
    # the previous genome, combined with their annotation status.



    # TODO After parsing flat file, prepare gene_id and gene_name appropriately


    # If the ticket type is 'replace', retrieve data from phamerator.
    # TODO will need to account for whether the phage_id exists in Phamerator or not.
    if bndl.ticket.type == "replace":
        phamerator_genome = \
            phamerator.create_phamerator_genome(sql_handle, gnm.id)


        # If any attributes in flat_file are set to 'retain', copy data
        # from the phamerator genome.
        phamerator.copy_data_from(bndl, "flat_file")





    # TODO confirm that evaluation checks that fields like
    # annotation_qc are structured properly - int and not str.

    evaluate.check_bundle_for_import(bndl)
    bndl.check_for_errors()







def import_into_database(bndl, sql_handle):
    """Construct collection of SQL statements for evaluated data and
    execute statements to add the data into the database."""



    # Create all SQL statements.
    index11 = 0
    while index11 < len(list_of_bundles):
            bndl = list_of_bundles[index11]
            if bndl.errors == 0:
                bndl.create_sql_statements()




    # TODO import all scrubbed add_replace data into Phamerator.



    pass






















###
