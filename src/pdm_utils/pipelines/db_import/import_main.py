"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""


import time
from datetime import datetime
import csv

from pdm_utils.functions import basic
from pdm_utils.functions import tickets
from pdm_utils.functions import flat_files
from pdm_utils.functions import phagesdb
from pdm_utils.functions import phamerator
from pdm_utils.classes import bundle
from pdm_utils.pipelines.db_import import evaluate
from pdm_utils.constants import constants

def import_io(sql_handle, genome_folder, import_table_file, filename_flag, test_run,
          description_field, run_mode):
    """Set up output directories, log files, etc. for import."""
    # Create output directories
    date = time.strftime("%Y%m%d")
    failed_folder = "%s_failed_upload_files" % date
    success_folder = "%s_successful_upload_files" % date
    new_failed_folder = basic.make_new_dir(genome_folder,failed_folder)
    if new_failed_folder != failed_folder:
        print("\nUnable to create failed_folder")
        sys.exit(1)
    new_success_folder = basic.make_new_dir(genome_folder,success_folder)
    if new_success_folder != success_folder:
        print("\nUnable to create success_folder")
        sys.exit(1)

    # Identify valid files in folder for evaluation.
    files_in_folder = basic.identify_files(genome_folder)

    try:
        eval_flags = constants.RUN_MODES[run_mode.lower()]
    except:
        run_mode = basic.expand_path(run_mode)
        eval_flags = basic.parse_flag_file(run_mode)
        symm_set = eval_flags.keys() ^ constants.RUN_MODE_BASE.keys()
        if len(symm_set) > 0:
            print("Error with the eval_flag_file")
            sys.exit(1)

    ticket_dict = prepare_tickets(import_table_file, eval_flags,
                                  description_field)

    # Proceed if there is at least one file to process.
    if len(files_in_folder) > 0:

        # TODO not sure how many elements (or what types) are returned.
        results = main(ticket_dict, files_in_folder,
                        sql_handle, eval_flags, description_field)








    # Now that all flat files and tickets have been evaluated,
    # provide summary of results...














# TODO incorporate this standard phage table query to retrieve all
# genome-level data from PhameratorDB:
# statement1 = "SELECT PhageID, Name, HostStrain, Sequence, status, \
#              Cluster2, DateLastModified, Accession, Subcluster2, \
#              AnnotationAuthor, AnnotationQC, RetrieveRecord \
#              FROM phage"







###
def prepare_tickets(import_table_file, eval_flags, description_field):
    """Prepare dictionary of pdm_utils Tickets."""

    # TODO parsing from import table:
    # 1. parse ticket data from table. = prepare_tickets()
    # 2. set case for all fields. = prepare_tickets()
    # 3. confirm all tickets have a valid type. = check_ticket_structure()
    # 4. populate Genome objects as necessary.
    # 5. retrieve data if needed.
    # 6. check for PhageID conflicts.
    # 7. confirm correct fields are populated based on ticket type.

    # Retrieve import ticket data.
    list_of_ticket_data = []
    with open(import_table_file,'r') as file:
        file_reader = csv.DictReader(file)
        for dict in file_reader:
            list_of_ticket_data.append(dict)
    ticket_list = []
    for dict in list_of_ticket_data:
        tkt = parse_import_ticket_data(data_dict=dict)
        if tkt is not None:
            ticket_list.append(tkt)
        else:
            # TODO record as error.
            pass




    # TODO unneeded if table is parsed using a dictionary reader.
    # with open(import_table_file,'r') as file:
    #     file_reader = csv.reader(file)
    #     for row in file_reader:
    #         list_of_ticket_data.append(row)
    # Convert ticket data to Ticket objects.
    # Data is returned as a list of ticket objects.
    # ticket_list = tickets.parse_import_tickets(list_of_ticket_data)
    # print("Tickets parsed")

    # Add the eval_flag dictionary and description_field to each ticket.
    for tkt in ticket_list:

        # TODO once pipeline is set up, verify copying is needed.
        # Copy eval_flags dictionary. Each ticket may alter some of the flags
        # based on ticket type, so copy dictionary to instantiate a
        # distinct dictionary.
        tkt.eval_flags = eval_flags.copy()
        tkt.description_field = description_field

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


    # TODO move tickets.compare_tickets() to evalute module since
    # it calls ticket check functions?
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

    # TODO handle ticket errors better.
    if ticket_errors > 0:
        print("Error generating tickets from import table.")
        sys.exit(1)
    else:
        return ticket_dict
###


def main(ticket_dict, files_in_folder, sql_handle=None):
    """The 'ticket_dict' parameter is a dictionary, where each
    element is a Ticket.
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


    # Retrieve data from phagesdb to create sets of
    # valid host genera, clusters, and subclusters.
    phagesdb_host_genera_set = phagesdb.create_host_genus_set()
    phagesdb_cluster_set, \
    phagesdb_subcluster_set = phagesdb.create_cluster_subcluster_sets()

    # TODO add special cluster and subcluster null values to sets,
    # including "UNK".

    # print("PhagesDB sets retrieved")


    # To minimize memory usage, each flat_file is evaluated one by one.
    bundle_count = 1
    for filename in files_in_folder:

        bndl = prepare_bundle(filename, ticket_dict)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)


        # Perform all evaluations based on the ticket type.
        evaluate.check_bundle_for_import(bndl)

        if bndl.tkt.type == "replace":
            genome_pair = bndl.genome_pair_dict["flat_file_phamerator"]
            evaluate.compare_genomes(genome_pair,
                check_replace=bndl.tkt.eval_flags["check_replace"])

        # TODO is there a function to retrieve Phamerator accession set?
        eval_gnm = bndl.genome_dict["flat_file"]
        evaluate.check_genome_for_import(
            eval_gnm,
            bndl.tkt,
            null_set=null_set,
            accession_set=phamerator_accession_set,
            phage_id_set=phamerator_phage_id_set,
            seq_set=phamerator_seq_set,
            host_genera_set=phagesdb_host_genera_set,
            cluster_set=phagesdb_cluster_set,
            subcluster_set=phagesdb_subcluster_set)

        # Check CDS features.
        x = 0
        while x < len(eval_gnm.cds_features):
            check_cds_for_import(eval_gnm.cds_features[x],
                check_locus_tag=tkt.eval_flags["check_locus_tag"],
                check_gene=tkt.eval_flags["check_gene"],
                check_description=tkt.eval_flags["check_description"],
                check_description_field=tkt.eval_flags["check_description_field"])
            x += 1

        # Check tRNA features.
        if tkt.eval_flags["check_trna"]:
            y = 0
            while y < len(gnm.trna_features):
                check_trna_for_import(gnm.trna_features[y])
                y += 1

        # Check Source features.
        z = 0
        while z < len(gnm.source_features):
            check_source_for_import(gnm.source_features[z],
                check_id_typo=tkt.eval_flags["check_id_typo"],
                check_host_typo=tkt.eval_flags["check_host_typo"],)
            z += 1


        # TODO confirm that evaluation checks that fields like
        # annotation_qc are structured properly - int and not str.


        # Now that all evaluations have been performed,
        # determine if there are any errors.
        bndl.check_for_errors()




        # TODO construct a basic.get_empty_ticket() function?
        ticket_data_dict =  tickets.parse_import_ticket_data(
                            tkt=bndl.ticket,
                            direction="ticket_to_dict")




        # TODO after evaluations, if sql argument option is True,
        # update the database as needed...

        bndl.check_for_errors()
        if errors == 0:

            # Now import the data into the database if there are no errors and
            # if there is MySQL connection data provided.
            if sql_handle is not None:
                bndl.create_sql_statements()
                sql_handle.execute_transaction(bndl.sql_queries)

                # If successful, keep track of query data.
                query_dict[bndl.ticket.phage_id] = bndl.sql_queries
            else:
                pass
            success_ticket_list.append(ticket_data_dict)
            success_filename_list.append(bndl.genome_dict["add"].filename)

        else:
            failed_ticket_list.append(ticket_data_dict)
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
        evaluate.check_ticket_structure(unmatched_ticket)
        evaluation_dict[bndl.id] = \
            {unmatched_ticket.id:unmatched_ticket.evaluations}
        bundle_count += 1
        ticket_data_dict =  tickets.parse_import_ticket_data(
                            unmatched_ticket,
                            direction="ticket_to_dict")
        failed_ticket_list.append(ticket_data_dict)



    #return list_of_bundles
    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict, query_dict)



def prepare_bundle(filename, ticket_dict, bundle_count):
    """Gather all genomic data needed to evaluate the flat file."""

    try:
        seqrecords = list(SeqIO.parse(filename, "genbank"))
    except:
        seqrecords = []

    if len(seqrecords) == 1:
        gnm = parse_genome_data(seqrecords[0], filepath=filename)

    else:
        # TODO throw an error of some sort.
        pass

    bndl = bundle.Bundle()
    bndl.id = bundle_count
    bndl.genome_dict[gnm.type] = gnm

    # Match ticket (if available) to flat file.
    matched_ticket = ticket_dict.pop(gnm.id, None)
    bndl.ticket = matched_ticket

    if bndl.ticket is not None:

        # Now that the flat file to be imported is parsed and matched
        # to a ticket, use the ticket to populate specific
        # genome-level fields such as host, cluster, subcluster, etc.
        tickets.copy_ticket_to_genome(bndl)
        flat_files.copy_data_to(bndl, "add")


    # Now check to see if there is any missing data for each genome, and
    # retrieve it from phagesdb.
    # If the ticket genome has fields set to 'retrieve', data is
    # retrieved from PhagesDB and populates a new Genome object.
    phagesdb.copy_data_from(bndl, "flat_file")


    # Each flat_file genome should now contain all requisite
    # data from PhagesDB.
    # Validate each genome by checking that each field is populated correctly.


    # TODO at some point annotation_qc and retrieve_record attributes
    # will need to be set. These are dependent on the ticket type.
    # If genomes are being replaced, these fields may be carried over from
    # the previous genome, combined with their annotation status.



    # If the ticket type is 'replace', retrieve data from phamerator.
    # TODO will need to account for whether the phage_id exists in Phamerator or not.
    if bndl.ticket.type == "replace":
        query = "SELECT * FROM phage"
        phamerator_genomes = \
            phamerator.parse_genome_data(sql_handle, phage_id_list=[gnm.id],
                phage_query=query)

        if len(phamerator_genomes) ==1:
            bndl.genome_dict[phamerator_genomes[0].type] = phamerator_genomes[0]
        else:
            # TODO throw an error of some sort if there is no matching genome in Phamerator.
            pass
        # If any attributes in flat_file are set to 'retain', copy data
        # from the phamerator genome.
        phamerator.copy_data_from(bndl, "flat_file")

    return bndl





























###
