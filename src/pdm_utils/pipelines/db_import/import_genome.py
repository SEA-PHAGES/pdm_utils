"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""


import time
from datetime import datetime
import csv
import os
import sys
import argparse
from pdm_utils.functions import basic
from pdm_utils.functions import tickets
from pdm_utils.functions import flat_files
from pdm_utils.functions import phagesdb
from pdm_utils.functions import phamerator
from pdm_utils.classes import bundle
from pdm_utils.constants import constants
from pdm_utils.classes import mysqlconnectionhandler as mch


def run_import(unparsed_args_list):
    """Verify the correct arguments are selected for import new genomes."""

    IMPORT_HELP = ("Pipeline to import new genome data into "
                   "a Phamerator MySQL database.")
    DATABASE_HELP = "Name of the MySQL database to import the genomes."
    INPUT_FOLDER_HELP = ("Path to the folder containing files to be processed.")
    IMPORT_GENOME_FOLDER_HELP = """
        Path to the folder containing
        GenBank-formatted flat files
        to be processed."""
    IMPORT_TABLE_HELP = """
        Path to the CSV-formatted table containing
        instructions to process each genome.
        Structure of import ticket table:
            1. Action to implement on the database (add, remove)
            2. PhageID
            3. Host genus
            4. Cluster
            5. Subcluster
            6. Annotation status (draft, final, unknown)
            7. Annotation authorship (hatfull, gbk)
            8. Gene description field (product, note, function)
            9. Accession
            10. Run mode
        """
    GENOME_ID_FIELD_HELP = """
        Indicates the flat file field that should be used
        as the unique identifier for the genome during import.
        """
    TEST_RUN_HELP = \
        ("Indicates whether the script should make any changes to the database. "
         "If False, the production run will implement all changes in the "
         "indicated database. If True, the test run will not "
         "implement any changes")
    RUN_MODE_HELP = \
        ("Indicates the evaluation configuration "
         "for importing genomes.")
    DESCRIPTION_FIELD_HELP = \
        ("Indicates the field in CDS features that is expected "
         "to store the gene description.")
    parser = argparse.ArgumentParser(description=IMPORT_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("input_folder", type=os.path.abspath,
        help=INPUT_FOLDER_HELP)
    parser.add_argument("import_table", type=os.path.abspath,
        help=IMPORT_TABLE_HELP)
    parser.add_argument("-gf", "--genome_id_field", type=str,
        default="organism_name", choices=["organism_name", "filename"],
        help=GENOME_ID_FIELD_HELP)
    parser.add_argument("-tr", "--test_run", action="store_false", default=True,
        help=TEST_RUN_HELP)
    parser.add_argument("-rm", "--run_mode", type=str.lower,
        choices=list(constants.RUN_MODES.keys()), default="phagesdb",
        help=RUN_MODE_HELP)
    parser.add_argument("-df", "--description_field", default="product",
        choices=list(constants.DESCRIPTION_FIELD_SET),
        help=DESCRIPTION_FIELD_HELP)
    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    # Validate args.
    sql_handle = mch.MySQLConnectionHandler()
    sql_handle.database = args.database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        print("No connection to the selected database.")
        sys.exit(1)

    valid_folder = basic.verify_path(args.input_folder, "dir")
    if not valid_folder:
        print("Invalid input folder.")
        sys.exit(1)

    valid_import_table = basic.verify_path(args.import_table, "file")
    if not valid_import_table:
        print("Invalid import table file.")
        sys.exit(1)

    # TODO run_mode needs to be improved.
    # It should be able to receive a filename that can be parsed
    # into an eval_flag dictionary.

    # If everything checks out, pass args to the main import pipeline:
    input_output(sql_handle=sql_handle,
        genome_folder=args.input_folder, import_table_file=args.import_table,
        genome_id_field=args.genome_id_field, test_run=args.test_run,
        description_field=args.description_field, run_mode=args.run_mode)




# TODO unittest.
def input_output(sql_handle=None, genome_folder="", import_table_file="",
    genome_id_field="", test_run=True, description_field="", run_mode=""):
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
        # TODO check order of parameters.
        results = main(ticket_dict, files_in_folder, sql_handle,
                       test_run, genome_id_field, filename_flag)

    # Now that all flat files and tickets have been evaluated,
    # provide summary of results...



# TODO unittest.
def prepare_tickets(import_table_file="", eval_flags={}, description_field=""):
    """Prepare dictionary of pdm_utils Tickets."""
    # 1. parse ticket data from table.
    # 2. set case for all fields.
    # 3. confirm all tickets have a valid type.
    # 6. check for duplicated values.
    # 7. confirm correct fields are populated based on ticket type.
    ticket_list = []
    tkt_errors = 0
    file_data = tickets.retrieve_ticket_data(import_table_file)
    for dict in file_data:
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
        tkt.eval_flags = eval_flags.copy()
        check_ticket(
            tkt,
            type_set=constants.IMPORT_TICKET_TYPE_SET,
            description_field_set=constants.DESCRIPTION_FIELD_SET,
            null_set=constants.EMPTY_SET,
            run_mode_set=constants.RUN_MODES.keys(),
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
        return None
    else:
        return ticket_dict


# TODO unittest.
def main(ticket_dict, files_in_folder, sql_handle=None, test_run=True,
         genome_id_field=""):
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

        bndl = prepare_bundle(filename, ticket_dict, sql_handle,
                              genome_id_field=genome_id_field, id=bundle_count)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)





        # TODO below create sub-function for the checks? There are a lot of
        # parameters for the genome checks that the others don't need.
        check_bundle(bndl)
        if (bndl.tkt is not None and bndl.tkt.type == "replace"):
            genome_pair = bndl.genome_pair_dict["flat_file_phamerator"]
            compare_genomes(genome_pair,
                check_replace=bndl.tkt.eval_flags["check_replace"])

        if "flat_file" in bndl.genome_dict.keys():
            ff_gnm = bndl.genome_dict["flat_file"]
            check_genome(
                ff_gnm,
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
            while x < len(ff_gnm.cds_features):
                check_cds(ff_gnm.cds_features[x],
                    check_locus_tag=tkt.eval_flags["check_locus_tag"],
                    check_gene=tkt.eval_flags["check_gene"],
                    check_description=tkt.eval_flags["check_description"],
                    check_description_field=tkt.eval_flags["check_description_field"])
                x += 1

            # Check tRNA features.
            if tkt.eval_flags["check_trna"]:
                y = 0
                while y < len(ff_gnm.trna_features):
                    check_trna(ff_gnm.trna_features[y])
                    y += 1

            # Check Source features.
            z = 0
            while z < len(ff_gnm.source_features):
                check_source(ff_gnm.source_features[z],
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
                result = sql_handle.execute_transaction(bndl.sql_queries)
                query_dict[bndl.ticket.phage_id] = bndl.sql_queries
                if result == 1:
                    # Log the error, if there was an issue with
                    # executing the statements.
                    pass
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
        check_bundle()
        evaluation_dict[bndl.id] = bndl.get_evaluations()
        ticket_data_dict =  tickets.parse_import_ticket_data(
                            unmatched_ticket,
                            direction="ticket_to_dict")
        failed_ticket_list.append(ticket_data_dict)
        bundle_count += 1

    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict, query_dict)


def prepare_bundle(filename="", ticket_dict={}, sql_handle=None,
                   genome_id_field="", id=None):
    """Gather all genomic data needed to evaluate the flat file.

    :param filename: Name of a GenBank-formatted flat file.
    :type filename: str
    :param ticket_dict: A dictionary of Tickets.
    :type ticket_dict: dict
    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :param id: Identifier to be assigned to the bundled dataset.
    :type id: int
    :returns:
        A pdm_utils Bundle object containing all data required to
        evaluate a flat file.
    :rtype: Bundle
    """
    bndl = bundle.Bundle()
    bndl.id = id
    seqrecord = flat_files.retrieve_genome_data(filename)
    if seqrecord is None:
        # TODO throw an error of some sort.
        # Assign some value to the Bundle and break, so that
        # the rest of this function is not executed.
        print("No record was retrieved from the file.")
    else:
        ff_gnm = flat_files.parse_genome_data(
                    seqrecord,
                    filepath=filename,
                    genome_id_field=genome_id_field,
                    gnm_type="flat_file")
        bndl.genome_dict[ff_gnm.type] = ff_gnm

        # Match ticket (if available) to flat file.
        bndl.ticket = ticket_dict.pop(ff_gnm.id, None)
        if bndl.ticket is not None:
            # With the flat file parsed and matched
            # to a ticket, use the ticket to populate specific
            # genome-level fields such as host, cluster, subcluster, etc.
            tickets.copy_ticket_to_genome(bndl)
            flat_files.copy_data(bndl, "add", "flat_file", flag="ticket")

            # Check to see if there is any missing data for each genome, and
            # retrieve it from phagesdb.
            # If the ticket genome has fields set to 'retrieve', data is
            # retrieved from PhagesDB and populates a new Genome object.
            ff_gnm.set_value_flag("retrieve")
            if ff_gnm._value_flag:
                phagesdb.copy_data(
                    bndl, "phagesdb", "flat_file", flag="retrieve")

            # If the ticket type is 'replace', retrieve data from phamerator.
            # If any attributes in flat_file are set to 'retain', copy data
            # from the phamerator genome.
            if bndl.ticket.type == "replace":
                if sql_handle is None:
                    print("Ticket %s is a 'replace' ticket but no"
                          " details about how to connect to the"
                          " Phamerator database have been provided."
                          " Unable to retrieve data."
                          % bndl.ticket.id)
                else:
                    query = "SELECT * FROM phage"
                    pmr_genomes =  phamerator.parse_genome_data(
                                       sql_handle=sql_handle,
                                       phage_id_list=[ff_gnm.id],
                                       phage_query=query,
                                       gnm_type="phamerator")
                    if len(pmr_genomes) == 1:
                        pmr_gnm = pmr_genomes[0]
                        bndl.genome_dict[pmr_gnm.type] = pmr_gnm
                        phamerator.copy_data(
                            bndl, "phamerator", "flat_file")
                    else:
                        print("There is no %s genome in the Phamerator"
                              " database. Unable to retrieve data."
                              % ff_gnm.id)
    return bndl




def check_bundle(bndl):
    """Check a Bundle for errors.

    Evaluate whether all genomes have been successfully grouped,
    and whether all genomes have been paired, as expected.
    Based on the ticket type, there are expected to be certain
    types of genomes and pairs of genomes in the bundle.

    :param bndl: A pdm_utils Bundle object.
    :type bndl: Bundle
    """
    bndl.check_ticket(eval_id="BNDL_001")
    if bndl.ticket is not None:
        tkt = bndl.ticket
        bndl.check_genome_dict("add", eval_id="BNDL_002")
        bndl.check_genome_dict("flat_file", eval_id="BNDL_003")
        bndl.check_genome_pair_dict("flat_file_add", eval_id="BNDL_004")

        # There may or may not be data retrieved from PhagesDB.
        tkt.set_value_flag("retrieve")
        if tkt._value_flag:
            bndl.check_genome_dict("phagesdb", eval_id="BNDL_005")
            bndl.check_genome_pair_dict("flat_file_phagesdb",
                                        eval_id="BNDL_006")

        if tkt.type == "replace":
            bndl.check_genome_dict("phamerator", eval_id="BNDL_007")

            # There may or may not be a genome_pair to retain some data.
            tkt.set_value_flag("retain")
            if tkt._value_flag:
                bndl.check_genome_pair_dict("add_phamerator",
                                            eval_id="BNDL_008")

            # There should be a genome_pair between the current phamerator
            # genome and the new flat_file genome.
            bndl.check_genome_pair_dict("flat_file_phamerator",
                                        eval_id="BNDL_009")


def check_ticket(tkt, type_set=set(), description_field_set=set(),
        null_set=set(), run_mode_set=set(), id_dupe_set=set(),
        phage_id_dupe_set=set(), accession_dupe_set=set()):
    """Evaluate a ticket to confirm it is structured appropriately.
    The assumptions for how each field is populated varies depending on
    the type of ticket.

    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param description_field_set: Valid description_field options.
    :type description_field_set: set
    :param null_set: Values that represent an empty field.
    :type null_set: set
    :param run_mode_set: Valid run mode options.
    :type run_mode_set: set
    :param id_dupe_set: Predetermined duplicate ticket ids.
    :type id_dupe_set: set
    :param phage_id_dupe_set: Predetermined duplicate PhageIDs.
    :type phage_id_dupe_set: set
    :param accession_dupe_set: Predetermined duplicate accessions.
    :type accession_dupe_set: set
    """
    # This function simply evaluates whether there is data in the
    # appropriate ticket attributes given the type of ticket.
    # It confirms that ticket attributes 'type', 'run_mode', and
    # 'description_field' are populated with specific values.
    # But it does not evaluate the quality of the data itself for
    # the other fields, since those are genome-specific fields and
    # can be checked within Genome objects.

    # Check for duplicated values.
    tkt.check_duplicate_id(id_dupe_set, eval_id="TKT_001")
    tkt.check_duplicate_phage_id(phage_id_dupe_set, eval_id="TKT_002")
    tkt.check_duplicate_accession(accession_dupe_set, eval_id="TKT_003")

    # Check these fields for specific values.
    tkt.check_type(type_set, True, eval_id="TKT_004")
    tkt.check_description_field(description_field_set, True, eval_id="TKT_005")
    tkt.check_run_mode(run_mode_set, True, eval_id="TKT_006")

    # For these fields, simply check that they are not empty.
    tkt.check_phage_id(null_set, False, eval_id="TKT_007")
    tkt.check_host_genus(null_set, False, eval_id="TKT_008")
    tkt.check_cluster(null_set, False, eval_id="TKT_009")
    tkt.check_annotation_status(null_set, False, eval_id="TKT_010")
    tkt.check_annotation_author(null_set, False, eval_id="TKT_011")
    tkt.check_retrieve_record(null_set, False, eval_id="TKT_012")

    # No need to evaluate the Accession and Subcluster fields
    # since they may or may not be populated.

    # Check if certain combinations of fields make sense.
    tkt.check_compatible_type_and_annotation_status(eval_id="TKT_013")


def check_genome(gnm, tkt, null_set=set(), phage_id_set=set(),
                 seq_set=set(), host_set=set(),
                 cluster_set=set(), subcluster_set=set(),
                 accession_set=set()):
    """Check a Genome object for errors.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param null_set: A set of values representing empty or null data.
    :type null_set: set
    :param phage_id_set: A set PhageIDs.
    :type phage_id_set: set
    :param seq_set: A set of genome sequences.
    :type seq_set: set
    :param host_set: A set of host genera.
    :type host_set: set
    :param cluster_set: A set of clusters.
    :type cluster_set: set
    :param subcluster_set: A set of subclusters.
    :type subcluster_set: set
    :param accession_set: A set of accessions.
    :type accession_set: set
    """

    if tkt.type == "add":
        gnm.check_id(phage_id_set | null_set, False, eval_id="GNM_001")
        gnm.check_name(phage_id_set | null_set, False, eval_id="GNM_002")
        gnm.check_sequence(seq_set | null_set, False, eval_id="GNM_003")

        # It is unusual for 'final' status genomes to be on 'add' tickets.
        gnm.check_annotation_status(check_set=set(["final"]), expect=False,
                                    eval_id="GNM_004")

        # If the genome is being added, and if it has an accession,
        # no other genome is expected to have an identical accession.
        # If the genome is being replaced, and if it has an accession,
        # the prior version of the genome may or may not have had
        # accession data, so no need to check for 'replace' tickets.
        if gnm.accession != "":
            gnm.check_accession(check_set=accession_set, expect=False,
                                eval_id="GNM_005")

    # 'replace' ticket checks.
    else:
        gnm.check_id(phage_id_set, True, eval_id="GNM_006")
        gnm.check_name(phage_id_set, True, eval_id="GNM_007")
        gnm.check_sequence(seq_set, True, eval_id="GNM_008")

        # It is unusual for 'draft' status genomes to be on 'replace' tickets.
        gnm.check_annotation_status(check_set=set(["draft"]), expect=False,
                                    eval_id="GNM_009")
    gnm.check_annotation_status(check_set=constants.ANNOTATION_STATUS_SET,
                                expect=True, eval_id="GNM_010")
    gnm.check_annotation_author(check_set=constants.ANNOTATION_AUTHOR_SET,
                                eval_id="GNM_011")
    gnm.check_retrieve_record(check_set=constants.RETRIEVE_RECORD_SET,
                              eval_id="GNM_012")
    gnm.check_cluster(cluster_set, True, eval_id="GNM_013")
    gnm.check_subcluster(subcluster_set, True, eval_id="GNM_014")
    gnm.check_subcluster_structure(eval_id="GNM_015")
    gnm.check_cluster_structure(eval_id="GNM_016")
    gnm.check_compatible_cluster_and_subcluster(eval_id="GNM_017")
    if tkt.eval_flags["check_seq"]:
        gnm.check_nucleotides(check_set=constants.DNA_ALPHABET,
                              eval_id="GNM_018")
    gnm.check_compatible_status_and_accession(eval_id="GNM_019")
    gnm.check_compatible_status_and_descriptions(eval_id="GNM_020")
    if tkt.eval_flags["check_id_typo"]:
        gnm.check_description_name(eval_id="GNM_021")
        gnm.check_source_name(eval_id="GNM_022")
        gnm.check_organism_name(eval_id="GNM_023")
    if tkt.eval_flags["check_host_typo"]:
        gnm.check_host_genus(host_set, True, eval_id="GNM_024")
        gnm.check_description_host_genus(eval_id="GNM_025")
        gnm.check_source_host_genus(eval_id="GNM_026")
        gnm.check_organism_host_genus(eval_id="GNM_027")
    if tkt.eval_flags["check_author"]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              eval_id="GNM_028")
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                                 expect=False, eval_id="GNM_029")
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET, expect=False,
                              eval_id="GNM_030")

    gnm.check_cds_feature_tally(eval_id="GNM_031")
    gnm.check_feature_ids(cds_ftr=True, trna_ftr=True, tmrna=True,
                          eval_id="GNM_032")

    # TODO confirm that these check_value_flag() are needed here.
    # Currently all "copy_data" functions run the check method
    # to throw an error if not all data was copied.
    gnm.set_value_flag("retrieve")
    gnm.check_value_flag(eval_id="GNM_033")
    gnm.set_value_flag("retain")
    gnm.check_value_flag(eval_id="GNM_034")


def check_source(src_ftr, check_id_typo=True, check_host_typo=True):
    """Check a Source object for errors."""
    if check_id_typo:
        src_ftr.check_organism_name(eval_id="SRC_001")
    if check_host_typo:
        src_ftr.check_organism_host_genus(eval_id="SRC_002")
        src_ftr.check_host_host_genus(eval_id="SRC_003")
        src_ftr.check_lab_host_host_genus(eval_id="SRC_004")


def check_cds(cds_ftr, check_locus_tag=True,
              check_gene=True, check_description=True,
              check_description_field=True):
    """Check a Cds object for errors."""
    cds_ftr.check_amino_acids(eval_id="CDS_001")
    cds_ftr.check_translation(eval_id="CDS_002")
    cds_ftr.check_translation_length(eval_id="CDS_003")
    cds_ftr.check_translation_table(eval_id="CDS_004")
    cds_ftr.check_coordinates(eval_id="CDS_005")
    cds_ftr.check_strand(eval_id="CDS_006")

    # These evaluations vary by genome type, stage of import, etc.
    if check_locus_tag:
        cds_ftr.check_locus_tag_present(eval_id="CDS_007")
        cds_ftr.check_locus_tag_structure(eval_id="CDS_008")
    if check_gene:
        cds_ftr.check_gene_present(eval_id="CDS_009")
        cds_ftr.check_gene_structure(eval_id="CDS_010")
    if check_locus_tag and check_gene:
        cds_ftr.check_compatible_gene_and_locus_tag(eval_id="CDS_011")
    if check_description:
        cds_ftr.check_generic_data(eval_id="CDS_012")
        cds_ftr.check_valid_description(eval_id="CDS_013")
    if check_description_field:
        cds_ftr.check_description_field(eval_id="CDS_014")


def compare_genomes(genome_pair, check_replace=True):
    """Compare two genomes to identify discrepancies."""
    genome_pair.compare_genome_sequence(eval_id="GP_001")
    genome_pair.compare_genome_length(eval_id="GP_002")
    genome_pair.compare_cluster(eval_id="GP_003")
    genome_pair.compare_subcluster(eval_id="GP_004")
    genome_pair.compare_accession(eval_id="GP_005")
    genome_pair.compare_host_genus(eval_id="GP_006")
    genome_pair.compare_author(eval_id="GP_007")
    if check_replace:
        genome_pair.compare_annotation_status("type","phamerator",
            "flat_file","draft","final", eval_id="GP_008")


# TODO implement.
# TODO unit test.
def check_trna(trna_obj):
    """Check a TrnaFeature object for errors."""

    pass




























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
