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
        results = main(ticket_dict, files_in_folder, sql_handle, test_run)

    # Now that all flat files and tickets have been evaluated,
    # provide summary of results...




def prepare_tickets(import_table_file, eval_flags, description_field):
    """Prepare dictionary of pdm_utils Tickets."""

    # TODO parsing from import table:
    # 1. parse ticket data from table.
    # 2. set case for all fields.
    # 3. confirm all tickets have a valid type.
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
    tkt_errors = 0
    for dict in list_of_ticket_data:
        tkt = tickets.parse_import_ticket_data(data_dict=dict)
        if tkt is not None:
            ticket_list.append(tkt)
        else:
            tkt_errors += 1

    # Identify duplicated ticket values.
    tkt_id_dupe_set, phage_id_dupe_set, accession_dupe_set = \
        tickets.identify_duplicates(ticket_list, null_set=set(["none"]))

    # For each ticket:
    # 1. Add the eval_flag dictionary and description_field to each ticket.
    #    Copy eval_flags dictionary. Each ticket may alter some of the flags
    #    based on ticket type, so copy dictionary to instantiate a
    #    distinct dictionary.
    # 2. Evaluate the tickets to ensure they are structured properly.
    #    At this point, the quality of the ticket data is not evaluated,
    #    just that the ticket contains fields populated or empty as expected.
    # 3. Create a dictionary of tickets based on the phage_id.
    # 4. Check for ticket errors.
    ticket_dict = {}
    x = 0
    while x < len(ticket_list):
        tkt = ticket_list[x]
        tkt.description_field = description_field
        # TODO once pipeline is set up, verify copying is needed.
        tkt.eval_flags = eval_flags.copy()
        evaluate.check_ticket_structure(
            tkt,
            type_set=constants.IMPORT_TICKET_TYPE_SET,
            description_field_set=constants.DESCRIPTION_FIELD_SET,
            null_set=constants.EMPTY_SET,
            run_mode_set=constants.RUN_MODE_SET,
            id_dupe_set=tkt_id_dupe_set,
            phage_id_dupe_set=phage_id_dupe_set,
            accession_dupe_set=accession_dupe_set)
        for evl in tkt.evaluations:
            if evl.status == "error":
                tkt_errors += 1
        ticket_dict[tkt.phage_id] = tkt
        x += 1

    # TODO handle ticket errors better.
    if tkt_errors > 0:
        print("Error generating tickets from import table.")
        sys.exit(1)
    else:
        return ticket_dict



def main(ticket_dict, files_in_folder, sql_handle=None, test_run=True):
    """Process GenBank-formatted flat files.

    :param ticket_dict:
        A dictionary
        WHERE
        key (str) = The ticket's phage_id
        value (Ticket) = The ticket
    :type ticket_dict: dict
    :param files_in_folder: A list of filenames to be parsed.
    :type files_in_folder: list
    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :param test_run:
        Indicates whether the database should be updated from the
        import tickets.
    :type test_run: bool
    :returns:
    :rtype:
    """
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
    phagesdb_cluster_set.add("UNK")
    # TODO add special cluster and subcluster null values to sets


    # To minimize memory usage, each flat_file is evaluated one by one.
    bundle_count = 1
    for filename in files_in_folder:

        bndl = prepare_bundle(filename, ticket_dict, id=bundle_count)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)

        evaluate.check_bundle_for_import(bndl)
        if (bndl.tkt is not None and bndl.tkt.type == "replace"):
            genome_pair = bndl.genome_pair_dict["flat_file_phamerator"]
            evaluate.compare_genomes(genome_pair,
                check_replace=bndl.tkt.eval_flags["check_replace"])

        if "flat_file" in bndl.genome_dict.keys():
            gnm = bndl.genome_dict["flat_file"]
            evaluate.check_genome_for_import(
                gnm,
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
            while x < len(gnm.cds_features):
                check_cds_for_import(gnm.cds_features[x],
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


        # Now that all evaluations have been performed,
        # determine if there are any errors.
        # TODO confirm that check_for_errors() is implemented.
        bndl.check_for_errors()
        ticket_data_dict =  tickets.parse_import_ticket_data(
                            tkt=bndl.ticket,
                            direction="ticket_to_dict")
        # TODO implement the get_evaluations() method.
        evaluation_dict[bndl.id] = bndl.get_evaluations()
        if errors == 0:
            # Now import the data into the database if there are no errors and
            # if it is not a test run.
            if not test_run:
                bndl.create_sql_statements()
                sql_handle.execute_transaction(bndl.sql_queries)
                query_dict[bndl.ticket.phage_id] = bndl.sql_queries
            else:
                pass
            success_ticket_list.append(ticket_data_dict)
            success_filename_list.append(bndl.genome_dict["flat_file"].filename)
        else:
            failed_ticket_list.append(ticket_data_dict)
            failed_filename_list.append(bndl.genome_dict["flat_file"].filename)
        bundle_count += 1


    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    key_list = list(ticket_dict.keys())
    for key in key_list:
        unmatched_ticket = ticket_dict.pop(key)
        bndl = bundle.Bundle()
        bndl.id = bundle_count
        bndl.ticket = unmatched_ticket
        evaluate.check_bundle_for_import()
        evaluation_dict[bndl.id] = bndl.get_evaluations()
        ticket_data_dict =  tickets.parse_import_ticket_data(
                            unmatched_ticket,
                            direction="ticket_to_dict")
        failed_ticket_list.append(ticket_data_dict)
        bundle_count += 1

    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict, query_dict)



def prepare_bundle(filename, ticket_dict, id=None):
    """Gather all genomic data needed to evaluate the flat file.

    :param filename: Name of a GenBank-formatted flat file.
    :type filename: str
    :param ticket_dict: A dictionary of Tickets.
    :type ticket_dict: dict
    :param id: Identifier to be assigned to the bundled dataset.
    :type id: int
    :returns:
        A pdm_utils Bundle object containing all data required to
        evaluate a flat file.
    :rtype: Bundle
    """

    bndl = bundle.Bundle()
    bndl.id = id
    try:
        seqrecords = list(SeqIO.parse(filename, "genbank"))
    except:
        seqrecords = []
    if len(seqrecords) != 1:
        # TODO throw an error of some sort.
        # Assign some value to the Bundle and break, so that
        # the rest of this function is not executed.
        pass

    gnm = parse_genome_data(seqrecords[0], filepath=filename)
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

    # TODO in progress below. This is partially redundant with
    # phagesdb.copy_data_to()
    bndl.gnm.set_value_flag("retrieve")
    if bndl.gnm._value_flag:
        phagesdb.copy_data_from(bndl, "flat_file")
        evaluate.check_phagesdb_genome()

        # Each flat_file genome should now contain all requisite
        # data from PhagesDB.
        # Validate each genome by checking that each field is populated correctly.
    # TODO in progress above.




    # TODO at some point annotation_qc and retrieve_record attributes
    # will need to be set. These are dependent on the ticket type.
    # If genomes are being replaced, these fields may be carried over from
    # the previous genome, combined with their annotation status.

    # If the ticket type is 'replace', retrieve data from phamerator.
    if bndl.ticket.type == "replace":
        query = "SELECT * FROM phage"
        phamerator_genomes = \
            phamerator.parse_genome_data(sql_handle, phage_id_list=[gnm.id],
                phage_query=query)
        if len(phamerator_genomes) == 1:
            bndl.genome_dict[phamerator_genomes[0].type] = phamerator_genomes[0]
        # If any attributes in flat_file are set to 'retain', copy data
        # from the phamerator genome.
        phamerator.copy_data_from(bndl, "flat_file")
    return bndl

























# TODO implement.
# TODO unit test.
# Cds object now contains a method to reset the primary description based
# on a user-selected choice.
#If other CDS fields contain descriptions, they can be chosen to
#replace the default import_cds_qualifier descriptions.
#Then provide option to verify changes.
#This block is skipped if user selects to do so.
# def check_description_field_choice():
#
#     if ignore_description_field_check != 'yes':
#
#         changed = ""
#         if (import_cds_qualifier != "product" and feature_product_tally > 0):
#            print "\nThere are %s CDS products found." % feature_product_tally
#            change_descriptions()
#
#            if question("\nCDS products will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[10]
#                 changed = "product"
#
#         if (import_cds_qualifier != "function" and feature_function_tally > 0):
#             print "\nThere are %s CDS functions found." % feature_function_tally
#             change_descriptions()
#
#             if question("\nCDS functions will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[11]
#                 changed = "function"
#         if (import_cds_qualifier != "note" and feature_note_tally > 0):
#
#             print "\nThere are %s CDS notes found." % feature_note_tally
#             change_descriptions()
#
#             if question("\nCDS notes will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[12]
#                 changed = "note"
#
#         if changed != "":
#             record_warnings += 1
#             write_out(output_file,"\nWarning: CDS descriptions only from the %s field will be retained." % changed)
#             record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)




###
