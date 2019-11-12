"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""


import argparse
import csv
from datetime import datetime, date
import logging
import os
import pathlib
import shutil
import sys
# import time
from tabulate import tabulate
from pdm_utils.functions import basic
from pdm_utils.functions import tickets
from pdm_utils.functions import flat_files
from pdm_utils.functions import phagesdb
from pdm_utils.functions import phamerator
from pdm_utils.classes import bundle
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.functions import run_modes
from pdm_utils.classes import mysqlconnectionhandler as mch

# Add a logger named after this module. Then add a null handler, which
# suppresses any output statements. This allows other modules that call this
# module to define the handler and output formats. If this module is the
# main module being called, the top level main function instantiates
# the root logger and configuration.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())




def main(unparsed_args_list):
    """Runs the complete import pipeline.

    This is the only function of the pipeline that requires user input.
    All other functions can be implemented from other scripts."""

    args = parse_args(unparsed_args_list)

    # Validate folders and files.
    args.input_folder = set_path(args.input_folder, kind="dir", expect=True)
    args.output_folder = set_path(args.output_folder, kind="dir", expect=True)
    args.import_table = set_path(args.import_table, kind="file", expect=True)
    args.log_file = pathlib.Path(args.output_folder, args.log_file)
    args.log_file = set_path(args.log_file, kind="file", expect=False)

    # Set up root logger.
    logging.basicConfig(filename=args.log_file, filemode="w",
                        level=logging.DEBUG,
                        format="%(name)s - %(levelname)s - %(message)s")
    logger.info("Folder and file arguments verified.")

    # Get connection to database.
    sql_handle = setup_sql_handle(args.database)
    logger.info(f"Connected to database: {args.database}.")

    # If everything checks out, pass on args for data input/output:
    data_io(sql_handle=sql_handle,
            genome_folder=args.input_folder,
            import_table_file=args.import_table,
            genome_id_field=args.genome_id_field,
            host_genus_field=args.host_genus_field,
            prod_run=args.prod_run,
            description_field=args.description_field,
            run_mode=args.run_mode,
            output_folder=args.output_folder,
            interactive=args.interactive)
    logger.info("Import complete.")


def setup_sql_handle(database):
    """Connect to a MySQL database."""
    sql_handle = mch.MySQLConnectionHandler()
    sql_handle.database = database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        logger.info(f"No connection to the {database} database. "
                    f"Valid credentials: {sql_handle.credential_status}. "
                    f"Valid database: {sql_handle._database_status}")
        sys.exit(1)
    else:
        return sql_handle


def set_path(path, kind=None, expect=True):
    """Confirm validity of path argument."""
    path = path.expanduser()
    path = path.resolve()
    result, msg = basic.verify_path2(path, kind=kind, expect=expect)
    if not result:
        print(msg)
        sys.exit(1)
    else:
        return path


def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for import new genomes."""

    IMPORT_HELP = ("Pipeline to import new genome data into "
                   "a Phamerator MySQL database.")
    DATABASE_HELP = "Name of the MySQL database to import the genomes."
    INPUT_FOLDER_HELP = ("Path to the folder containing files to be processed.")
    OUTPUT_FOLDER_HELP = ("Path to the folder containing to store results.")
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
    HOST_GENUS_FIELD_HELP = """
        Indicates the flat file field that should be used
        to identify the host genus for the genome during import.
        """
    PROD_RUN_HELP = \
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
    LOG_FILE_HELP = \
        ("Indicates the name of the file to log the import results. "
         "This will be created in the output folder.")
    INTERACTIVE_HELP = \
        ("Indicates whether interactive evaluation of data is permitted.")

    parser = argparse.ArgumentParser(description=IMPORT_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("input_folder", type=pathlib.Path,
        help=INPUT_FOLDER_HELP)
    parser.add_argument("import_table", type=pathlib.Path,
        help=IMPORT_TABLE_HELP)
    parser.add_argument("-g", "--genome_id_field", type=str.lower,
        default="_organism_name",
        choices=["_organism_name", "_description_name",
                 "_source_name", "filename"],
        help=GENOME_ID_FIELD_HELP)
    parser.add_argument("-hg", "--host_genus_field", type=str.lower,
        default="_organism_host_genus",
        choices=["_organism_host_genus", "_description_host_genus",
                 "_source_host_genus"],
        help=GENOME_ID_FIELD_HELP)
    parser.add_argument("-p", "--prod_run", action="store_true",
        default=False, help=PROD_RUN_HELP)
    parser.add_argument("-r", "--run_mode", type=str.lower,
        choices=list(run_modes.RUN_MODES.keys()), default="phagesdb",
        help=RUN_MODE_HELP)
    parser.add_argument("-d", "--description_field", type=str.lower,
        default="product", choices=list(constants.DESCRIPTION_FIELD_SET),
        help=DESCRIPTION_FIELD_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path("/tmp/"), help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-l", "--log_file", type=str, default="import.log",
        help=LOG_FILE_HELP)
    parser.add_argument("-i", "--interactive", action="store_true",
        default=False, help=INTERACTIVE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    return args




def data_io(sql_handle=None, genome_folder=pathlib.Path(),
    import_table_file=pathlib.Path(), genome_id_field="", host_genus_field="",
    prod_run=False, description_field="", run_mode="",
    output_folder=pathlib.Path(), interactive=False):
    """Set up output directories, log files, etc. for import."""
    # Create output directories

    logger.info("Setting up environment.")
    CURRENT_DATE = date.today().strftime("%Y%m%d")
    results_folder = f"{CURRENT_DATE}_results"
    results_folder = pathlib.Path(results_folder)
    results_path = basic.make_new_dir(output_folder, results_folder, attempt=3)
    if results_path is None:
        logger.info("Unable to create output_folder.")
        sys.exit(1)

    # Get the files to process.
    files_to_process = basic.identify_files(genome_folder, set([".DS_Store"]))
    if len(files_to_process) == 0:
        logger.info("There are no flat files to evaluate.")
        sys.exit(1)

    # Get the tickets.
    eval_flags = run_modes.get_eval_flag_dict(run_mode.lower())
    run_mode_eval_dict = {"run_mode": run_mode, "eval_flag_dict": eval_flags}
    ticket_dict = prepare_tickets(import_table_file,
                                  run_mode_eval_dict,
                                  description_field,
                                  constants.IMPORT_TABLE_STRUCTURE)

    if ticket_dict is None:
        logger.info("Invalid import table. Unable to evaluate flat files.")
        sys.exit(1)

    start_count = phamerator.get_phage_table_count(sql_handle)

    # Evaluate files and tickets.
    results_tuple = process_files_and_tickets(
                        ticket_dict, files_to_process,
                        sql_handle=sql_handle,
                        prod_run=prod_run,
                        genome_id_field=genome_id_field,
                        host_genus_field=host_genus_field,
                        interactive=interactive)

    final_count = phamerator.get_phage_table_count(sql_handle)

    success_ticket_list = results_tuple[0]
    failed_ticket_list = results_tuple[1]
    success_filepath_list = results_tuple[2]
    failed_filepath_list = results_tuple[3]
    evaluation_dict = results_tuple[4]

    # Output data.
    logger.info("Logging successful tickets and files.")
    headers = constants.IMPORT_TABLE_STRUCTURE["order"]
    if (len(success_ticket_list) > 0 or len(success_filepath_list) > 0):
        success_path = pathlib.Path(results_path, "success")
        success_path.mkdir()
        if len(success_ticket_list) > 0:
            success_tkt_file = pathlib.Path(success_path, "import_tickets.csv")
            tickets.export_ticket_data(success_ticket_list, success_tkt_file, headers)
        if len(success_filepath_list) > 0:
            success_genomes_path = pathlib.Path(success_path, "genomes")
            success_genomes_path.mkdir()
            for file in success_filepath_list:
                new_file = pathlib.Path(success_genomes_path, file.name)
                shutil.move(str(file), str(new_file))

    logger.info("Logging failed tickets and files.")
    if (len(failed_ticket_list) > 0 or len(failed_filepath_list) > 0):
        failed_path = pathlib.Path(results_path, "fail")
        failed_path.mkdir()
        if len(failed_ticket_list) > 0:
            failed_tkt_file = pathlib.Path(failed_path, "import_tickets.csv")
            tickets.export_ticket_data(failed_ticket_list, failed_tkt_file, headers)
        if len(failed_filepath_list) > 0:
            failed_genomes_path = pathlib.Path(failed_path, "genomes")
            failed_genomes_path.mkdir()
            for file in failed_filepath_list:
                new_file = pathlib.Path(failed_genomes_path, file.name)
                shutil.move(str(file), str(new_file))

    # logger.info("Logging evaluations.")
    # log_evaluations(evaluation_dict)

    logger.info(
        ("Summary of import: "
        f"\n{start_count} genome(s) in the database before import. "
        f"\n{final_count} genome(s) in the database after import. "
        f"\n{len(success_ticket_list)} ticket(s) successfully processed. "
        f"\n{len(failed_ticket_list)} ticket(s) NOT processed. "
        f"\n{len(success_filepath_list)} genome(s) successfully imported. "
        f"\n{len(failed_filepath_list)} genome(s) NOT imported. ")
        )






def log_evaluations(dict_of_dict_of_lists):
    """Export evaluations to log.

    Structure of the evaluation dictionary:
        {1: {"bundle": [eval_object1, ...],
             "ticket": [eval_object1, ...],
             "genome": [eval_object1, ...]},
         2: {...}}
    """
    logger.info("Logging all evaluations.")
    for key1 in dict_of_dict_of_lists:
        msg1 = f"Evaluations for bundle {key1}:"
        logger.info(msg1)
        dict_of_lists = dict_of_dict_of_lists[key1]
        for key2 in dict_of_lists:
            msg2 = f"Evaluations for {key2}:"
            logger.info(msg2)
            evl_list = dict_of_lists[key2]
            for evl in evl_list:
                msg3 = (f"Evaluation: {evl.id}. "
                        f"Status: {evl.status}. "
                        f"Definition: {evl.definition}. "
                        f"Result: {evl.result}.")
                logger.info(msg3)




def prepare_tickets(import_table_file=pathlib.Path(), run_mode_eval_dict=None,
        description_field="", table_structure_dict={}):
    """Prepare dictionary of pdm_utils Tickets."""
    # 1. Parse ticket data from table.
    # 2. Set case for certain fields and set default values for missing fields.
    # 3. Add the eval_flag dictionary and description_field to each ticket
    #    provided from command line arguments if they are not provided in the
    #    import table.
    # 4. Check for duplicated phage_ids and ticket ids.
    # 5. For each ticket, evaluate the tickets to ensure they are
    #    structured properly. At this point, the quality of the
    #    ticket data is not evaluated, just that the ticket contains
    #    fields populated or empty as expected.
    # 6. Create a dictionary of tickets based on the phage_id.
    required_keys = table_structure_dict["required"]
    optional_keys = table_structure_dict["optional"]
    keywords = table_structure_dict["keywords"]
    valid_retain = table_structure_dict["valid_retain"]
    valid_retrieve = table_structure_dict["valid_retrieve"]
    valid_add = table_structure_dict["valid_add"]

    list_of_tkts = []
    tkt_errors = 0
    logger.info("Retrieving ticket data.")
    list_of_data_dicts = tickets.retrieve_ticket_data(import_table_file)
    logger.info("Constructing tickets.")
    list_of_tkts = tickets.construct_tickets(list_of_data_dicts, run_mode_eval_dict,
                    description_field, required_keys, optional_keys,
                    keywords)
    if len(list_of_data_dicts) != len(list_of_tkts):
        tkt_errors += 1

    logger.info("Identifying duplicate tickets.")
    tkt_id_dupes, phage_id_dupes = tickets.identify_duplicates(list_of_tkts)


    ticket_dict = {}
    x = 0
    while x < len(list_of_tkts):
        tkt = list_of_tkts[x]
        tkt_summary = (f"ID: {tkt.id}. "
                       f"Type: {tkt.type}. "
                       f"PhageID: {tkt.phage_id}.")
        logger.info(f"Checking ticket structure for: {tkt_summary}.")
        check_ticket(tkt,
                     type_set=constants.IMPORT_TICKET_TYPE_SET,
                     description_field_set=constants.DESCRIPTION_FIELD_SET,
                     run_mode_set=run_modes.RUN_MODES.keys(),
                     id_dupe_set=tkt_id_dupes,
                     phage_id_dupe_set=phage_id_dupes,
                     retain_set=valid_retain,
                     retrieve_set=valid_retrieve,
                     add_set=valid_add)
        for evl in tkt.evaluations:
            evl_summary = (f"Evaluation: {evl.id}. "
                           f"Status: {evl.status}. "
                           f"Definition: {evl.definition}. "
                           f"Result: {evl.result}.")
            if evl.status == "error":
                tkt_errors += 1
                logger.error(evl_summary)
            else:
                logger.info(evl_summary)
        ticket_dict[tkt.phage_id] = tkt
        x += 1

    if tkt_errors > 0:
        logger.info("Error generating tickets from import table.")
        return None
    else:
        logger.info("Tickets were successfully generated from import table.")
        return ticket_dict

def process_files_and_tickets(ticket_dict, files_in_folder, sql_handle=None,
                              prod_run=False, genome_id_field="",
                              host_genus_field="", interactive=False):
    """Process GenBank-formatted flat files.

    :param ticket_dict:
        A dictionary
        WHERE
        key (str) = The ticket's phage_id
        value (Ticket) = The ticket
    :type ticket_dict: dict
    :param files_in_folder: A list of filepaths to be parsed.
    :type files_in_folder: list
    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :param prod_run:
        Indicates whether the database should be updated from the
        import tickets.
    :type prod_run: bool
    :returns:
    :rtype:
    """
    # Alias for different types of genomes gathered and processed.
    file_ref = "flat_file"
    ticket_ref = "ticket"
    retain_ref = "phamerator"
    retrieve_ref = "phagesdb"

    # Will hold all results data. These can be used by various types of
    # user interfaces to present the summary of import.
    success_ticket_list = []
    failed_ticket_list = []
    success_filepath_list = []
    failed_filepath_list = []
    evaluation_dict = {}


    # TODO should add some sort of flag to set default empty values
    # for phagesdb sets to account for creation of databases that
    # do not rely on phagesdb?
    # Retrieve data from phagesdb to create sets of
    # valid host genera, clusters, and subclusters.
    # Cluster "UNK" may or may not already be present, but it is valid.
    # If there is no subcluster, value may be empty string or "none".
    phagesdb_host_genera_set = phagesdb.create_host_genus_set()
    results_tuple = phagesdb.create_cluster_subcluster_sets()
    phagesdb_cluster_set = results_tuple[0]
    phagesdb_subcluster_set = results_tuple[1]
    phagesdb_cluster_set.add("UNK")


    # To minimize memory usage, each flat_file is evaluated one by one.
    bundle_count = 1
    for filepath in files_in_folder:
        replace_gnm_pair_key = file_ref + "_" + retain_ref
        logger.info(f"Preparing data for file: {filepath.name}.")
        bndl = prepare_bundle(filepath=filepath, ticket_dict=ticket_dict,
                              sql_handle=sql_handle,
                              genome_id_field=genome_id_field,
                              host_genus_field=host_genus_field,
                              id=bundle_count,
                              file_ref=file_ref, ticket_ref=ticket_ref,
                              retrieve_ref=retrieve_ref, retain_ref=retain_ref,
                              interactive=interactive)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)
        logger.info(f"Checking file: {filepath.name}.")
        run_checks(bndl,
                   accession_set=phamerator_accession_set,
                   phage_id_set=phamerator_phage_id_set,
                   seq_set=phamerator_seq_set,
                   host_genus_set=phagesdb_host_genera_set,
                   cluster_set=phagesdb_cluster_set,
                   subcluster_set=phagesdb_subcluster_set,
                   file_ref=file_ref, ticket_ref=ticket_ref,
                   retrieve_ref=retrieve_ref, retain_ref=retain_ref)

        # If interactive is set to True, interate through each eval
        # stored in each object in the bundle.
        # If there is an error status, ask use if it is correct, and
        # downgrade the status as needed.
        if interactive:
            review_evaluations(bndl)
        bndl.check_for_errors()

        logger.info("Logging evaluations.")
        dict_of_eval_lists = bndl.get_evaluations()
        log_evaluations({bndl.id: dict_of_eval_lists})
        evaluation_dict[bndl.id] = dict_of_eval_lists
        # evaluation_dict[bndl.id] = bndl.get_evaluations()


        result = import_into_db(bndl, sql_handle=sql_handle,
                                gnm_key=file_ref, prod_run=prod_run)
        if result:
            success_ticket_list.append(bndl.ticket.data_dict)
            success_filepath_list.append(filepath)
        else:
            if bndl.ticket is not None:
                failed_ticket_list.append(bndl.ticket.data_dict)
            failed_filepath_list.append(filepath)
        bundle_count += 1

    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    if len(ticket_dict.keys()) > 0:
        logger.info("Processing unmatched tickets.")
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)
        key_list = list(ticket_dict.keys())
        for key in key_list:
            bndl = bundle.Bundle()
            bndl.ticket = ticket_dict.pop(key)
            bndl.id = bundle_count
            logger.info("Checking data for ticket: "
                        f"{bndl.ticket.id}, {bndl.ticket.phage_id}.")
            run_checks(bndl,
                      accession_set=phamerator_accession_set,
                      phage_id_set=phamerator_phage_id_set,
                      seq_set=phamerator_seq_set,
                      host_genus_set=phagesdb_host_genera_set,
                      cluster_set=phagesdb_cluster_set,
                      subcluster_set=phagesdb_subcluster_set,
                      file_ref=file_ref, ticket_ref=ticket_ref,
                      retrieve_ref=retrieve_ref, retain_ref=retain_ref)
            evaluation_dict[bndl.id] = bndl.get_evaluations()
            failed_ticket_list.append(bndl.ticket.data_dict)
            bundle_count += 1

    return (success_ticket_list, failed_ticket_list, success_filepath_list,
            failed_filepath_list, evaluation_dict)


def prepare_bundle(filepath=pathlib.Path(), ticket_dict={}, sql_handle=None,
                   genome_id_field="", host_genus_field="", id=None,
                   file_ref="", ticket_ref="", retrieve_ref="", retain_ref="",
                   id_conversion_dict={}, interactive=False):
    """Gather all genomic data needed to evaluate the flat file.

    :param filepath: Name of a GenBank-formatted flat file.
    :type filepath: Path
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
    seqrecord = flat_files.retrieve_genome_data(filepath)
    if seqrecord is None:
        logger.info(f"No record was retrieved from the file: {filepath.name}.")
    else:
        logger.info(f"Parsing record from the file: {filepath.name}.")
        ff_gnm = flat_files.parse_genome_data(
                    seqrecord,
                    filepath=filepath,
                    genome_id_field=genome_id_field,
                    gnm_type=file_ref,
                    host_genus_field=host_genus_field)

        # Some phage names in the flat file are spelled differently
        # than in the database. In order to match the flat file
        # with the ticket (which presumably contains the name spelling
        # in the database), the genome ID needs to be changed.
        # PHAGE_ID_DICT key = incorrect spelling ; value = correct spelling
        if ff_gnm.id in id_conversion_dict.keys():
            ff_gnm.id = id_conversion_dict[ff_gnm.id]

        bndl.genome_dict[ff_gnm.type] = ff_gnm

        # Match ticket (if available) to flat file.
        bndl.ticket = ticket_dict.pop(ff_gnm.id, None)
        if bndl.ticket is None:
            logger.info(f"No matched ticket for file: {filepath.name}.")
        else:
            logger.info(f"Preparing ticket data for file: {filepath.name}.")
            # With the flat file parsed and matched
            # to a ticket, use the ticket to populate specific
            # genome-level fields such as host, cluster, subcluster, etc.

            if len(bndl.ticket.data_add) > 0:
                tkt_gnm = tickets.get_genome(bndl.ticket, gnm_type=ticket_ref)
                bndl.genome_dict[tkt_gnm.type] = tkt_gnm

                # Copy ticket data to flat file. Since the ticket data has been
                # added to a genome object using genome methods, the data
                # can be directly passed from one genome object to another.
                for attr in bndl.ticket.data_add:
                    attr_value = getattr(tkt_gnm, attr)
                    setattr(ff_gnm, attr, attr_value)

            # Check to see if there is any missing data for each genome, and
            # retrieve it from phagesdb.
            # If the ticket genome has fields set to 'retrieve', data is
            # retrieved from PhagesDB and populates a new Genome object.
            if len(bndl.ticket.data_retrieve) > 0:
                pdb_gnm = phagesdb.get_genome(bndl.ticket.phage_id,
                                              gnm_type=retrieve_ref)

                # TODO unit test 'is not None' block.
                # pdb_gnm is None if PhageID not in PhagesDB
                if pdb_gnm is not None:
                    bndl.genome_dict[pdb_gnm.type] = pdb_gnm

                    for attr in bndl.ticket.data_retrieve:
                        attr_value = getattr(pdb_gnm, attr)
                        setattr(ff_gnm, attr, attr_value)

            # If the ticket type is 'replace', retrieve data from phamerator.
            # If any attributes in flat_file are set to 'retain', copy data
            # from the phamerator genome.
            if bndl.ticket.type == "replace":

                if sql_handle is None:
                    logger.info(
                          f"Ticket {bndl.ticket.id} is a 'replace' ticket "
                          "but no details about how to connect to the "
                          "Phamerator database have been provided. "
                          "Unable to retrieve data.")
                else:
                    query = "SELECT * FROM phage"
                    pmr_genomes =  phamerator.parse_genome_data(
                                       sql_handle=sql_handle,
                                       phage_id_list=[ff_gnm.id],
                                       phage_query=query,
                                       gnm_type=retain_ref)
                    if len(pmr_genomes) == 1:
                        pmr_gnm = pmr_genomes[0]
                        bndl.genome_dict[pmr_gnm.type] = pmr_gnm

                        # The ticket may indicate some data should be retained.
                        for attr in bndl.ticket.data_retain:
                            attr_value = getattr(pmr_gnm, attr)
                            setattr(ff_gnm, attr, attr_value)

                        # Pair the genomes for comparison evaluations.
                        gnm_pair = genomepair.GenomePair()
                        bndl.set_genome_pair(gnm_pair, ff_gnm.type, pmr_gnm.type)
                    else:
                        logger.info(f"There is no {ff_gnm.id} genome "
                                    "in the Phamerator database. "
                                    "Unable to retrieve data.")

            ff_gnm.set_cluster_subcluster(value="internal")
            set_cds_descriptions(ff_gnm, bndl.ticket, interactive=interactive)
    return bndl


def run_checks(bndl, accession_set=set(), phage_id_set=set(),
               seq_set=set(), host_genus_set=set(), cluster_set=set(),
               subcluster_set=set(), file_ref="", ticket_ref="",
               retrieve_ref="", retain_ref=""):
    """Run checks on the different elements of a bundle object."""
    logger.info("Checking data.")
    check_bundle(bndl, ticket_ref=ticket_ref, file_ref=file_ref,
                 retrieve_ref=retrieve_ref, retain_ref=retain_ref)
    tkt = bndl.ticket
    if tkt is not None:
        eval_flags = tkt.eval_flags
        gnm_pair_key = file_ref + "_" + retain_ref
        if (tkt.type == "replace" and
                gnm_pair_key in bndl.genome_pair_dict.keys()):
            genome_pair = bndl.genome_pair_dict[gnm_pair_key]
            compare_genomes(genome_pair, eval_flags)

        if file_ref in bndl.genome_dict.keys():
            gnm = bndl.genome_dict[file_ref]
            check_genome(gnm, tkt.type, eval_flags,
                         accession_set=accession_set, phage_id_set=phage_id_set,
                         seq_set=seq_set, host_genus_set=host_genus_set,
                         cluster_set=cluster_set, subcluster_set=subcluster_set)

            # Check CDS features.
            x = 0
            while x < len(gnm.cds_features):
                check_cds(gnm.cds_features[x], eval_flags,
                          description_field=tkt.description_field)
                x += 1

            # Check tRNA features.
            if eval_flags["check_trna"]:
                y = 0
                while y < len(gnm.trna_features):
                    check_trna(gnm.trna_features[y], eval_flags)
                    y += 1

            # Check Source features.
            z = 0
            while z < len(gnm.source_features):
                check_source(gnm.source_features[z], eval_flags,
                             host_genus=gnm.host_genus)
                z += 1


def review_evaluations(bndl):
    """Iterate through all objects stored in the bundle.
    If there are errors, review whether status should be changed."""
    print(f"\n\nReviewing evaluations for bundle: {bndl.id}")
    review_evaluation_list(bndl.evaluations)
    if bndl.ticket is not None:
        print(f"nReviewing evaluations for ticket: {bndl.ticket.id}, "
              f"{bndl.ticket.type}, {bndl.ticket.phage_id}.")
        review_evaluation_list(bndl.ticket.evaluations)
    for key in bndl.genome_dict.keys():
        gnm = bndl.genome_dict[key]
        print(f"Reviewing evaluations for genome: {gnm.id}, {gnm.type}.")
        review_evaluation_list(gnm.evaluations)

        # Capture the exit status for each CDS feature. If user exits
        # the review at any point, skip all of the other CDS features.
        print("Reviewing evaluations for CDS features.")
        print("Answer 'yes'/'no', or 'exit' to exit from review.")
        exit = False
        x = 0
        while (exit is False and x < len(gnm.cds_features)):
            cds_ftr = gnm.cds_features[x]
            exit = review_evaluation_list(cds_ftr.evaluations)
            x += 1

        print(f"Reviewing evaluations for source features.")
        for source_ftr in gnm.source_features:
            review_evaluation_list(source_ftr.evaluations)

        # TODO implement trna and tmrna features

    print(f"Reviewing evaluations for paired genomes.")
    for key in bndl.genome_pair_dict.keys():
        genome_pair = bndl.genome_pair_dict[key]
        review_evaluation_list(genome_pair.evaluations)

def review_evaluation_list(evaluation_list):
    """Iterate through all evaluations and review 'error' results.
    """
    exit = False
    x = 0
    while (exit is False and x < len(evaluation_list)):
        evl = evaluation_list[x]
        if evl.status == "error":
            print("\n\nThe following evaluation is set to 'error':")
            print(f"Evaluation ID: {evl.id}")
            print(f"Status: {evl.status}")
            print(f"Definition: {evl.definition}")
            print(f"Result: {evl.result}")
            prompt = "\nShould this evaluation be downgraded to 'warning'? "
            result = basic.ask_yes_no(prompt=prompt, response_attempt=3)
            if result is None:
                exit = True
            elif result is True:
                evl.status = "warning"
                evl.result = evl.result + " Downgraded from 'error' status."
            else:
                pass
        x += 1
    return exit


def set_cds_descriptions(gnm, tkt, interactive=False):
    """Set the primary CDS descriptions.
    """
    # If interactive is selected, the user can review the CDS descriptions.
    # The ticket indicates where the CDS descriptions are stored.
    if interactive:
        tkt.description_field = review_cds_descriptions(
                                    gnm.cds_features,
                                    tkt.description_field)
    gnm.set_cds_descriptions(tkt.description_field)
    gnm.tally_cds_descriptions()


def review_cds_descriptions(feature_list, description_field):
    """Iterate through all CDS features and review descriptions.
    """
    # Print a table of CDS data for user to review.
    summary = []
    x = 20
    y = "..."
    for cds in feature_list:
        short_product = basic.truncate_value(cds.processed_product, x, y)
        short_function = basic.truncate_value(cds.processed_function, x, y)
        short_note = basic.truncate_value(cds.processed_note, x, y)
        dict = {"CDS ID": cds.id,
                "Start": cds.left,
                "Stop": cds.right,
                "Strand": cds.strand,
                "Product": short_product,
                "Function": short_function,
                "Note": short_note
                }
        summary.append(dict)
    print(tabulate(summary, headers="keys"))
    print(f"Descriptions in the {description_field} field will be imported.")
    prompt = "Is this correct? "
    result = basic.ask_yes_no(prompt=prompt, response_attempt=3)

    # If the data is not correct, let the user select which description
    # field is most appropriate.
    if result == False:
        field_options = constants.DESCRIPTION_FIELD_SET - {description_field}
        print("Select the field in which the descriptions are stored.")
        new_field = basic.choose_from_list(list(field_options))
        if new_field is None:
            new_field = description_field
        else:
            print(f"The '{new_field}' field was selected.")
    else:
        new_field = description_field

    # Return either the original or the new description field.
    return new_field






def check_bundle(bndl, ticket_ref="", file_ref="", retrieve_ref="", retain_ref=""):
    """Check a Bundle for errors.

    Evaluate whether all genomes have been successfully grouped,
    and whether all genomes have been paired, as expected.
    Based on the ticket type, there are expected to be certain
    types of genomes and pairs of genomes in the bundle.

    :param bndl: A pdm_utils Bundle object.
    :type bndl: Bundle
    """
    logger.info(f"Checking bundle: {bndl.id}.")
    bndl.check_ticket(eval_id="BNDL_001") # TODO LOCK
    if bndl.ticket is not None:

        bndl.check_genome_dict(file_ref, expect=True, eval_id="BNDL_002") # TODO LOCK
        if file_ref in bndl.genome_dict.keys():
            bndl.check_compatible_type_and_annotation_status(
                    file_ref, eval_id="BNDL_003")

        tkt = bndl.ticket
        if len(tkt.data_add) > 0:
            bndl.check_genome_dict(ticket_ref, expect=True, eval_id="BNDL_004")# TODO LOCK

        # There may or may not be data retrieved from PhagesDB.
        if len(tkt.data_retrieve) > 0:
            bndl.check_genome_dict(retrieve_ref,
                                   expect=True, eval_id="BNDL_005")# TODO LOCK

        if tkt.type == "replace":
            bndl.check_genome_dict(retain_ref, expect=True, eval_id="BNDL_006")# TODO LOCK

            # There should be a genome_pair between the current phamerator
            # genome and the new flat_file genome.
            pair_key = f"{file_ref}_{retain_ref}"
            bndl.check_genome_pair_dict(pair_key, eval_id="BNDL_007")# TODO LOCK


def check_ticket(tkt, type_set=set(), description_field_set=set(),
        run_mode_set=set(), id_dupe_set=set(), phage_id_dupe_set=set(),
        retain_set=set(), retrieve_set=set(), add_set=set()):
    """Evaluate a ticket to confirm it is structured appropriately.
    The assumptions for how each field is populated varies depending on
    the type of ticket.

    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param description_field_set: Valid description_field options.
    :type description_field_set: set
    :param run_mode_set: Valid run mode options.
    :type run_mode_set: set
    :param id_dupe_set: Predetermined duplicate ticket ids.
    :type id_dupe_set: set
    :param phage_id_dupe_set: Predetermined duplicate PhageIDs.
    :type phage_id_dupe_set: set
    """
    # This function simply evaluates whether there is data in the
    # appropriate ticket attributes given the type of ticket.
    # It confirms that ticket attributes 'type', 'run_mode', and
    # 'description_field' are populated with specific values.
    # But it does not evaluate the quality of the data itself for
    # the other fields, since those are genome-specific fields and
    # can be checked within Genome objects.
    logger.info(f"Checking ticket: {tkt.id}, {tkt.type}, {tkt.phage_id}.")

    # Check for duplicated values.
    tkt.check_attribute("id", id_dupe_set,
                        expect=False, eval_id="TKT_001")
    tkt.check_attribute("phage_id", phage_id_dupe_set,
                        expect=False, eval_id="TKT_002")

    # Check these fields for specific values.
    tkt.check_attribute("type", type_set,
                        expect=True, eval_id="TKT_003")
    tkt.check_attribute("description_field", description_field_set,
                        expect=True, eval_id="TKT_004")
    tkt.check_attribute("run_mode", run_mode_set,
                        expect=True, eval_id="TKT_005")

    # This method may be refactored so that it accepts a list of
    # valid flag dict keys. But this has already been verified earlier
    # in the script, so it could be redundant.
    tkt.check_eval_flags(expect=True, eval_id="TKT_006")

    # For these fields, simply check that they are not empty.
    tkt.check_attribute("phage_id", {""},
                        expect=False, eval_id="TKT_007")

    # Check how genome attributes will be determined.
    tkt.check_compatible_type_and_data_retain(eval_id="TKT_009")
    tkt.check_valid_data_source("data_add", add_set,
                                eval_id="TKT_010")
    tkt.check_valid_data_source("data_retain", retain_set,
                                eval_id="TKT_011")
    tkt.check_valid_data_source("data_retrieve", retrieve_set,
                                eval_id="TKT_012")


def check_genome(gnm, tkt_type, eval_flags, phage_id_set=set(),
                 seq_set=set(), host_genus_set=set(),
                 cluster_set=set(), subcluster_set=set(),
                 accession_set=set()):
    """Check a Genome object for errors.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
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
    logger.info(f"Checking genome: {gnm.id}, {gnm.type}.")

    if tkt_type == "add":
        gnm.check_attribute("id", phage_id_set | {""},
                            expect=False, eval_id="GNM_001") # TODO LOCK
        gnm.check_attribute("name", phage_id_set | {""},
                            expect=False, eval_id="GNM_002") # TODO LOCK
        gnm.check_attribute("seq", seq_set | {constants.EMPTY_GENOME_SEQ},
                            expect=False, eval_id="GNM_003") # TODO LOCK


        # If the genome is being added, and if it has an accession,
        # no other genome is expected to have an identical accession.
        # If the genome is being replaced, and if it has an accession,
        # the prior version of the genome may or may not have had
        # accession data, so no need to check for 'replace' tickets.
        if gnm.accession != "":
            gnm.check_attribute("accession", accession_set,
                                expect=False, eval_id="GNM_004") # TODO LOCK

    # 'replace' ticket checks.
    else:
        gnm.check_attribute("id", phage_id_set,
                            expect=True, eval_id="GNM_005") # TODO LOCK
        gnm.check_attribute("seq", seq_set,
                            expect=True, eval_id="GNM_006") # TODO LOCK

    # Depending on the annotation_status of the genome,
    # CDS features are expected to contain or not contain descriptions.
    # Draft genomes should not have any descriptions.
    # Final genomes should not have any descriptions.
    # There are no expectations for other types of genomes.

    # Also, Draft annotations should not have accession data.

    check_name = basic.edit_suffix(gnm.name, "add",
                                   suffix=constants.NAME_SUFFIX)

    if gnm.annotation_status == "draft":
        gnm.check_attribute("name", {check_name}, expect=True, eval_id="GNM_007")
        gnm.check_magnitude("_cds_processed_descriptions_tally", "=", 0,
                            eval_id="GNM_008")
        gnm.check_attribute("accession", {""}, expect=True, eval_id="GNM_009")

    elif gnm.annotation_status == "final":
        gnm.check_attribute("name", {check_name}, expect=False, eval_id="GNM_010")
        gnm.check_magnitude("_cds_processed_descriptions_tally", ">", 0,
                            eval_id="GNM_011")

    else:
        pass

    check_id = basic.edit_suffix(gnm.name, "add", suffix=constants.NAME_SUFFIX)
    gnm.check_attribute("id", {check_id}, expect=False, eval_id="GNM_012")
    gnm.check_attribute("annotation_status", constants.ANNOTATION_STATUS_SET,
                        expect=True, eval_id="GNM_013") # TODO LOCK
    gnm.check_attribute("annotation_author", constants.ANNOTATION_AUTHOR_SET,
                        expect=True, eval_id="GNM_014") # TODO LOCK
    gnm.check_attribute("retrieve_record", constants.RETRIEVE_RECORD_SET,
                        expect=True, eval_id="GNM_015") # TODO LOCK
    gnm.check_attribute("cluster", cluster_set,
                        expect=True, eval_id="GNM_016")
    gnm.check_attribute("subcluster", subcluster_set | {"none"},
                        expect=True, eval_id="GNM_017")
    gnm.check_attribute("cluster_subcluster", cluster_set | subcluster_set,
                        expect=True, eval_id="GNM_018")
    gnm.check_attribute("translation_table", {11},
                        expect=True, eval_id="GNM_019")
    gnm.check_attribute("host_genus", host_genus_set,
                        expect=True, eval_id="GNM_020")
    gnm.check_cluster_structure(eval_id="GNM_021")
    gnm.check_subcluster_structure(eval_id="GNM_022")
    gnm.check_compatible_cluster_and_subcluster(eval_id="GNM_023")
    gnm.check_magnitude("date", ">", constants.EMPTY_DATE, eval_id="GNM_024")
    gnm.check_magnitude("gc", ">", -0.0001, eval_id="GNM_025")
    gnm.check_magnitude("gc", "<", 100.0001, eval_id="GNM_026")
    gnm.check_magnitude("length", ">", 0, eval_id="GNM_027")
    gnm.check_magnitude("_cds_features_tally", ">", 0, eval_id="GNM_028")

    # TODO set trna=True and tmrna=True after they are implemented.
    gnm.check_feature_coordinates(cds_ftr=True, trna_ftr=False, tmrna=False,
                                  strand=False, eval_id="GNM_029") # TODO LOCK

    if eval_flags["check_seq"]:
        gnm.check_nucleotides(check_set=constants.DNA_ALPHABET,
                              eval_id="GNM_030")

    if eval_flags["check_id_typo"]:
        gnm.compare_two_attributes("id", "_description_name",
                                   expect_same=True, eval_id="GNM_031")
        gnm.compare_two_attributes("id", "_source_name",
                                   expect_same=True, eval_id="GNM_032")
        gnm.compare_two_attributes("id", "_organism_name",
                                   expect_same=True, eval_id="GNM_033")

    if eval_flags["check_host_typo"]:

        # TODO change these checks in parallel with source feature checks.
        # Get a set of host genus synonyms and use check_attribute()

        gnm.compare_two_attributes("host_genus", "_description_host_genus",
                                   expect_same=True, eval_id="GNM_034")
        gnm.compare_two_attributes("host_genus", "_source_host_genus",
                                   expect_same=True, eval_id="GNM_035")
        gnm.compare_two_attributes("host_genus", "_organism_host_genus",
                                   expect_same=True, eval_id="GNM_036")

    if eval_flags["check_author"]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              expect=True, eval_id="GNM_037")
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                              expect=False, eval_id="GNM_038")
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              expect=False, eval_id="GNM_039")




def check_source(src_ftr, eval_flags, host_genus=""):
    """Check a Source object for errors."""
    logger.info(f"Checking source feature: {src_ftr.id}.")

    if eval_flags["check_id_typo"]:
        src_ftr.check_attribute("_organism_name", {src_ftr.genome_id},
                                expect=True, eval_id="SRC_001")

    if eval_flags["check_host_typo"]:
        host_genus_synonyms = basic.get_synonyms(
                                host_genus, constants.HOST_GENUS_SYNONYMS)

        # Only evaluate attributes if the field is not empty, since they
        # are not required to be present.
        if src_ftr.organism != "":
            src_ftr.check_attribute("_organism_host_genus", host_genus_synonyms,
                                    expect=True, eval_id="SRC_002")
        if src_ftr.host != "":
            src_ftr.check_attribute("_host_host_genus", host_genus,
                                    expect=True, eval_id="SRC_003")
        if src_ftr.lab_host != "":
            src_ftr.check_attribute("_lab_host_host_genus", host_genus,
                                    expect=True, eval_id="SRC_004")





def check_cds(cds_ftr, eval_flags, description_field="product"):
    """Check a Cds object for errors."""
    logger.info(f"Checking CDS feature: {cds_ftr.id}.")

    cds_ftr.check_amino_acids(check_set=constants.PROTEIN_ALPHABET,
                              eval_id="CDS_001") # TODO LOCK
    cds_ftr.check_translation(eval_id="CDS_002") # TODO LOCK
    cds_ftr.check_translation_present(eval_id="CDS_003") # TODO LOCK
    cds_ftr.check_translation_table(check_table=11, eval_id="CDS_004")
    cds_ftr.check_coordinates(eval_id="CDS_005") # TODO LOCK
    cds_ftr.check_strand(format="fr_short", case=True, eval_id="CDS_006") # TODO LOCK
    if eval_flags["check_locus_tag"]:
        cds_ftr.check_locus_tag_present(expect=True, eval_id="CDS_007")

        # TODO this check could be improved to take into account the prefix.
        cds_ftr.check_locus_tag_structure(check_value=None, only_typo=True,
            case=True, eval_id="CDS_008")
    if eval_flags["check_gene"]:
        cds_ftr.check_gene_present(expect=True, eval_id="CDS_009")
        cds_ftr.check_gene_structure(eval_id="CDS_010")
    if (eval_flags["check_locus_tag"] and eval_flags["check_gene"]):
        cds_ftr.check_compatible_gene_and_locus_tag(eval_id="CDS_011")

    # if eval_flags["check_description"]:
        # TODO the "check_generic_data" method should be implemented at the genome level.
        # cds_ftr.check_generic_data(eval_id="CDS_012")
        # TODO not implemented yet: cds_ftr.check_valid_description(eval_id="CDS_013")
    if eval_flags["check_description_field"]:
        cds_ftr.check_description_field(attribute=description_field,
                                        eval_id="CDS_014")

def compare_genomes(genome_pair, eval_flags):
    """Compare two genomes to identify discrepancies."""
    logger.info("Comparing data between two genomes: "
                f"Genome 1. ID: {genome_pair.genome1.id}. "
                f"Type: {genome_pair.genome1.type}. "
                f"Genome 1. ID: {genome_pair.genome2.id}. "
                f"Type: {genome_pair.genome2.type}."
                )

    genome_pair.compare_attribute("id",
        expect_same=True, eval_id="GP_001")
    genome_pair.compare_attribute("seq",
        expect_same=True, eval_id="GP_002")
    genome_pair.compare_attribute("length",
        expect_same=True, eval_id="GP_003")
    genome_pair.compare_attribute("cluster",
        expect_same=True, eval_id="GP_004")
    genome_pair.compare_attribute("subcluster",
        expect_same=True, eval_id="GP_005")
    genome_pair.compare_attribute("host_genus",
        expect_same=True, eval_id="GP_007")
    genome_pair.compare_attribute("annotation_author",
        expect_same=True, eval_id="GP_008")
    genome_pair.compare_attribute("translation_table",
        expect_same=True, eval_id="GP_009")
    genome_pair.compare_attribute("retrieve_record",
        expect_same=True, eval_id="GP_010")

    if eval_flags["check_replace"]:
        # The following checks assume that:
        # 'genome1' slot = new genome to be evaluated
        # 'genome2' slot = the current genome in Phamerator

        # The new genome to be evaluated is expected to be
        # newer than the current genome annotations in Phamerator.
        genome_pair.compare_date("newer", eval_id="GP_015")

        if genome_pair.genome2.annotation_status == "draft":
            # It is expected that the replacing genome name no longer
            # retains the "_Draft" suffix, so the name should change.
            genome_pair.compare_attribute("name",
                expect_same=False, eval_id="GP_011")

            # The status should change from 'draft'.
            genome_pair.compare_attribute("annotation_status",
                expect_same=False, eval_id="GP_012")
        else:
            genome_pair.compare_attribute("name",
                expect_same=True, eval_id="GP_013")
            genome_pair.compare_attribute("annotation_status",
                expect_same=True, eval_id="GP_014")
            genome_pair.compare_attribute("accession",
                expect_same=True, eval_id="GP_006")




# TODO implement.
# TODO unit test.
def check_trna(trna_obj, eval_flags):
    """Check a TrnaFeature object for errors."""

    pass


def import_into_db(bndl, sql_handle=None, gnm_key="", prod_run=False):
    """Import data into the MySQL database."""
    IMPORT_DATE = datetime.today().replace(hour=0,
                                           minute=0,
                                           second=0,
                                           microsecond=0)

    if bndl._errors == 0:
        import_gnm = bndl.genome_dict[gnm_key]
        logger.info("Importing data into the database for "
                    f"genome: {import_gnm.id}.")

        # If locus_tag data should not be imported, clear all locus_tag data.
        if bndl.ticket.eval_flags["import_locus_tag"] is False:
            logger.info("Clearing locus_tag data for "
                        f"genome: {import_gnm.id}.")
            import_gnm.clear_locus_tags()

        # Update the date field to reflect the day of import.
        import_gnm.date = IMPORT_DATE
        bndl.sql_statements = phamerator.create_genome_statements(
                                import_gnm, bndl.ticket.type)
        if prod_run:
            result = sql_handle.execute_transaction(bndl.sql_statements)
            if result == 1:
                logger.info("Error executing statements to import data.")
                result = False
            else:
                result = True
                logger.info("Data successfully imported.")
                logger.info("The following SQL statements were executed:")
                for statement in bndl.sql_statements:
                    logger.info(statement[:150] + "...")
        else:
            result = True
            logger.info("Data can be imported if set to production run.")
    else:
        result = False
        logger.info("Data contains errors, so it will not be imported.")
    return result


###
