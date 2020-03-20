"""Primary pipeline to process and evaluate data to be imported
into the MySQL database."""


import argparse
import csv
from datetime import datetime, date
import logging
import pathlib
import shutil
import sys
from tabulate import tabulate
import pdm_utils # to get version number.
from pdm_utils.functions import basic
from pdm_utils.functions import tickets
from pdm_utils.functions import flat_files
from pdm_utils.functions import phagesdb
from pdm_utils.functions import mysqldb
from pdm_utils.classes import bundle
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants, eval_descriptions
from pdm_utils.functions import eval_modes

# Add a logger named after this module. Then add a null handler, which
# suppresses any output statements. This allows other modules that call this
# module to define the handler and output formats. If this module is the
# main module being called, the top level main function configures
# the root logger and primary file handle.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


DEFAULT_OUTPUT_FOLDER = "/tmp/"
IMPORT_DATE = datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_import"
SUCCESS_FOLDER = "success"
FAIL_FOLDER = "fail"
GENOMES_FOLDER = "genomes"
LOGS_FOLDER = "logs"
VERSION = pdm_utils.__version__
EDD = eval_descriptions.EVAL_DESCRIPTIONS

def main(unparsed_args_list):
    """Runs the complete import pipeline.

    This is the only function of the pipeline that requires user input.
    All other functions can be implemented from other scripts.

    :param unparsed_args_list:
        List of strings representing command line arguments.
    :type unparsed_args_list: list
    """

    args = parse_args(unparsed_args_list)

    # Validate folders and files.
    args.input_folder = basic.set_path(args.input_folder, kind="dir",
                                       expect=True)
    args.output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    args.import_table = basic.set_path(args.import_table, kind="file",
                                       expect=True)

    results_folder = pathlib.Path(RESULTS_FOLDER)
    results_path = basic.make_new_dir(args.output_folder,
                                      results_folder, attempt=50)

    if results_path is None:
        print("Unable to create results folder.")
        sys.exit(1)

    args.log_file = pathlib.Path(results_path, args.log_file)

    # Set up root logger.
    logging.basicConfig(filename=args.log_file, filemode="w",
                        level=logging.DEBUG,
                        format="pdm_utils import: %(levelname)s: %(message)s")
                        # format="%(name)s - %(levelname)s - %(message)s")

    logger.info(f"pdm_utils version: {VERSION}")
    logger.info(f"Import run date: {CURRENT_DATE}")
    logger.info(f"Command line arguments: {' '.join(unparsed_args_list)}")
    logger.info(f"Results directory: {results_path}")
    logger.info("Command line arguments verified.")

    # Verify database connection and schema compatibility.
    engine = mysqldb.connect_to_db(args.database)
    logger.info(f"Connected to database: {args.database}.")
    mysqldb.check_schema_compatibility(engine, "the import pipeline")
    logger.info(f"Schema version is compatible.")

    # If everything checks out, pass on args for data input/output.
    data_io(engine=engine,
            genome_folder=args.input_folder,
            import_table_file=args.import_table,
            genome_id_field=args.genome_id_field,
            host_genus_field=args.host_genus_field,
            prod_run=args.prod_run,
            description_field=args.description_field,
            eval_mode=args.eval_mode,
            output_folder=results_path,
            interactive=args.interactive)

    logger.info("Import complete.")


def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for import new genomes.

    :param unparsed_args_list:
        List of strings representing command line arguments.
    :type unparsed_args_list: list
    :returns: ArgumentParser Namespace object containing the parsed args.
    :rtype: Namespace
    """

    IMPORT_HELP = ("Pipeline to import new genome data into "
                   "a MySQL database.")
    DATABASE_HELP = "Name of the MySQL database to import the genomes."
    INPUT_FOLDER_HELP = ("Path to the folder containing GenBank-formatted "
                         "flat files to be processed.")
    OUTPUT_FOLDER_HELP = ("Path to the folder to store results.")
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
            10. Eval mode
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
         "If True, the production run will implement all changes in the "
         "indicated database. If False, the test run will not "
         "implement any changes.")
    EVAL_MODE_HELP = \
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
    parser.add_argument("-e", "--eval_mode", type=str.lower,
        choices=list(eval_modes.EVAL_MODES.keys()), default="final",
        help=EVAL_MODE_HELP)
    parser.add_argument("-d", "--description_field", type=str.lower,
        default="product", choices=list(constants.DESCRIPTION_FIELD_SET),
        help=DESCRIPTION_FIELD_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER), help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-l", "--log_file", type=str, default="import.log",
        help=LOG_FILE_HELP)
    parser.add_argument("-i", "--interactive", action="store_true",
        default=False, help=INTERACTIVE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    return args




def data_io(engine=None, genome_folder=pathlib.Path(),
    import_table_file=pathlib.Path(), genome_id_field="", host_genus_field="",
    prod_run=False, description_field="", eval_mode="",
    output_folder=pathlib.Path(), interactive=False):
    """Set up output directories, log files, etc. for import.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param genome_folder: Path to the folder of flat files.
    :type genome_folder: Path
    :param import_table_file: Path to the import table file.
    :type import_table_file: Path
    :param genome_id_field:
        The SeqRecord attribute that stores the genome identifier/name.
    :type genome_id_field: str
    :param host_genus_field:
        The SeqRecord attribute that stores the host genus identifier/name.
    :type host_genus_field: str
    :param prod_run: Indicates whether MySQL statements will be executed.
    :type prod_run: bool
    :param description_field:
        The SeqFeature attribute that stores the feature's description.
    :type description_field: str
    :param eval_mode:
        Name of the evaluation mode to evaluation genomes.
    :type eval_mode: str
    :param output_folder: Path to the folder to store results.
    :type output_folder: Path
    :param interactive:
        Indicates whether user is able to interact with genome evaluations
        at run time.
    :type interactive: bool
    """

    logger.info("Setting up environment.")
    # Get the files to process.
    files_to_process = basic.identify_contents(genome_folder, kind="file",
                                               ignore_set=set([".DS_Store"]))
    file_msg = f"There are {len(files_to_process)} flat files to evaluate."
    if len(files_to_process) == 0:
        logger.error(file_msg)
        print(file_msg)
        sys.exit(1)
    else:
        log_and_print(file_msg, True)

    # Get the tickets.
    eval_flags = eval_modes.get_eval_flag_dict(eval_mode.lower())
    eval_data_dict = {"eval_mode": eval_mode, "eval_flag_dict": eval_flags}
    ticket_dict = prepare_tickets(import_table_file,
                                  eval_data_dict,
                                  description_field,
                                  constants.IMPORT_TABLE_STRUCTURE)
    if ticket_dict is None:
        logger.error("Invalid import table. Unable to evaluate flat files.")
        sys.exit(1)
    else:
        progress = f"There are {len(ticket_dict.keys())} import ticket(s)."
        log_and_print(progress, True)

    # Set up output folders. Creating all folders in advance, since
    # logging the flatfile log files need a valid path.
    success_path = pathlib.Path(output_folder, SUCCESS_FOLDER)
    success_genomes_path = pathlib.Path(success_path, GENOMES_FOLDER)
    success_logs_path = pathlib.Path(success_path, LOGS_FOLDER)

    failed_path = pathlib.Path(output_folder, FAIL_FOLDER)
    failed_genomes_path = pathlib.Path(failed_path, GENOMES_FOLDER)
    failed_logs_path = pathlib.Path(failed_path, LOGS_FOLDER)

    output_folders = [success_path, success_genomes_path, success_logs_path,
                      failed_path, failed_genomes_path, failed_logs_path]
    for folder in output_folders:
        folder.mkdir()
    log_folder_paths_dict = {SUCCESS_FOLDER: success_logs_path,
                             FAIL_FOLDER: failed_logs_path}
    start_count = mysqldb.get_phage_table_count(engine)

    # Evaluate files and tickets.
    results_tuple = process_files_and_tickets(
                        ticket_dict, files_to_process,
                        engine=engine,
                        prod_run=prod_run,
                        genome_id_field=genome_id_field,
                        host_genus_field=host_genus_field,
                        interactive=interactive,
                        log_folder_paths_dict=log_folder_paths_dict)
    success_ticket_list = results_tuple[0]
    failed_ticket_list = results_tuple[1]
    success_filepath_list = results_tuple[2]
    failed_filepath_list = results_tuple[3]
    evaluation_dict = results_tuple[4]

    final_count = mysqldb.get_phage_table_count(engine)

    # Output data.
    logger.info("Logging successful tickets and files.")
    summary = [
        "Summary of import: ",
        f"{start_count} genome(s) in the database before import.",
        f"{final_count} genome(s) in the database after import.",
        f"{len(success_ticket_list)} ticket(s) successfully processed.",
        f"{len(success_filepath_list)} genome(s) successfully processed."
        ]

    headers = constants.IMPORT_TABLE_STRUCTURE["order"]
    if (len(success_ticket_list) > 0 or len(success_filepath_list) > 0):
        if len(success_ticket_list) > 0:
            success_tkt_file = pathlib.Path(success_path, "import_table.csv")
            basic.export_data_dict(success_ticket_list, success_tkt_file,
                                       headers, include_headers=True)
        if len(success_filepath_list) > 0:
            for file in success_filepath_list:
                new_file = pathlib.Path(success_genomes_path, file.name)
                shutil.copy(str(file), str(new_file))

    logger.info("Logging failed tickets and files.")
    if (len(failed_ticket_list) > 0 or len(failed_filepath_list) > 0):
        if len(failed_ticket_list) > 0:
            failed_tkt_file = pathlib.Path(failed_path, "import_table.csv")
            basic.export_data_dict(failed_ticket_list, failed_tkt_file,
                                       headers, include_headers=True)
            summary.append(f"{len(failed_ticket_list)} ticket(s) NOT implemented:")
            for tkt in failed_ticket_list:
                tkt_summary = (f"Ticket Type: {tkt['type']}. "
                               f"PhageID: {tkt['phage_id']}.")
                summary.append(f"{tkt_summary}")

        if len(failed_filepath_list) > 0:
            for file in failed_filepath_list:
                new_file = pathlib.Path(failed_genomes_path, file.name)
                shutil.copy(str(file), str(new_file))
            summary.append(f"{len(failed_filepath_list)} genome(s) NOT imported:")
            for file in failed_filepath_list:
                summary.append(f"{file.name}")

    # Now check which folders are empty and remove them.
    # List needs to be reversed. For instance, the main 'success' folder
    # can't be removed until the 'genomes' and 'logs' folders are removed.
    for folder in reversed(output_folders):
        if len(basic.identify_contents(folder, kind=None)) == 0:
            folder.rmdir()

    print("\n\n\n")
    for line in summary:
        print(line)
        logger.info(line)

    # Close all connections in the connection pool.
    engine.dispose()


def log_evaluations(dict_of_dict_of_lists, logfile_path=None):
    """Export evaluations to log.

    :param dict_of_dict_of_lists:
        Dictionary of evaluation dictionaries.
        Key1 = Bundle ID.
        Value1 = dictionary for each object in the Bundle.
        Key2 = object type ('bundle', 'ticket', etc.)
        Value2 = List of evaluation objects.
        Example structure:
            {1: {"bundle": [eval_object1, ...],
                 "ticket": [eval_object1, ...],
                 "genome": [eval_object1, ...]},
             2: {...}}
    :type dict_of_dict_of_lists: dict
    :param logfile_path: Path to the log file.
    :type logfile_path: Path
    """

    # Create a new logging filehandle to only capture errors
    # specific to the particular flatfile.
    if logfile_path is not None:
        flatfile_logger = logging.FileHandler(logfile_path, mode="w")
        flatfile_logger.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        flatfile_logger.setFormatter(formatter)
        logger.addHandler(flatfile_logger)

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
                msg3 = str(evl)
                if evl.status == "warning":
                    logger.warning(msg3)
                elif evl.status == "error":
                    logger.error(msg3)
                else:
                    logger.info(msg3)

    if logfile_path is not None:
        flatfile_logger.close()
        logger.removeHandler(flatfile_logger)



def prepare_tickets(import_table_file=pathlib.Path(), eval_data_dict=None,
        description_field="", table_structure_dict={}):
    """Prepare dictionary of pdm_utils ImportTickets.

    :param import_table_file: same as for data_io().
    :param description_field: same as for data_io().
    :param eval_data_dict:
        Evaluation data dictionary
        Key1 = "eval_mode"
        Value1 = Name of the eval_mode
        Key2 = "eval_flag_dict"
        Value2 = Dictionary of evaluation flags.
    :type eval_data_dict: dict
    :param table_structure_dict:
        Dictionary describing structure of the import table.
    :type table_structure_dict: dict
    :returns:
        Dictionary of pdm_utils ImportTicket objects.
        If a problem was encountered parsing the import table, None is returned.
    :rtype: dict
    """
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
    valid_parse = table_structure_dict["valid_parse"]

    list_of_tkts = []
    tkt_errors = 0
    logger.info("Retrieving ticket data.")
    list_of_data_dicts = basic.retrieve_data_dict(import_table_file)
    logger.info("Constructing tickets.")
    list_of_tkts = tickets.construct_tickets(list_of_data_dicts, eval_data_dict,
                    description_field, required_keys, optional_keys,
                    keywords)
    if len(list_of_data_dicts) != len(list_of_tkts):
        tkt_errors += 1

    logger.info("Identifying duplicate tickets.")
    tkt_id_dupes, phage_id_dupes = tickets.identify_duplicates(list_of_tkts)

    ticket_dict = {}
    for x in range(len(list_of_tkts)):
        tkt = list_of_tkts[x]
        tkt_summary = (f"ID: {tkt.id}. "
                       f"Type: {tkt.type}. "
                       f"PhageID: {tkt.phage_id}.")
        logger.info(f"Checking ticket structure for: {tkt_summary}")
        string_list = []
        for key in tkt.eval_flags:
            string_list.append(f"{key}: {tkt.eval_flags[key]}")
        logger.info(f"Eval flags are {', '.join(string_list)}")
        check_ticket(tkt,
                     type_set=constants.IMPORT_TICKET_TYPE_SET,
                     description_field_set=constants.DESCRIPTION_FIELD_SET,
                     eval_mode_set=eval_modes.EVAL_MODES.keys(),
                     id_dupe_set=tkt_id_dupes,
                     phage_id_dupe_set=phage_id_dupes,
                     retain_set=valid_retain,
                     retrieve_set=valid_retrieve,
                     add_set=valid_add,
                     parse_set=valid_parse)
        for evl in tkt.evaluations:
            evl_summary = (f"Evaluation: {evl.id}. "
                           f"Status: {evl.status}. "
                           f"Definition: {evl.definition} "
                           f"Result: {evl.result}")
            if evl.status == "error":
                tkt_errors += 1
                logger.error(evl_summary)
            else:
                logger.info(evl_summary)
        ticket_dict[tkt.phage_id] = tkt

    if tkt_errors > 0:
        logger.error("Error generating tickets from import table.")
        return None
    else:
        logger.info("Tickets were successfully generated from import table.")

        return ticket_dict

def process_files_and_tickets(ticket_dict, files_in_folder, engine=None,
                              prod_run=False, genome_id_field="",
                              host_genus_field="", interactive=False,
                              log_folder_paths_dict=None):
    """Process GenBank-formatted flat files and import tickets.

    :param ticket_dict:
        A dictionary
        WHERE
        key (str) = The ticket's phage_id
        value (Ticket) = The ticket
    :type ticket_dict: dict
    :param files_in_folder: A list of filepaths to be parsed.
    :type files_in_folder: list
    :param engine: same as for data_io().
    :param prod_run: same as for data_io().
    :param genome_id_field: same as for data_io().
    :param host_genus_field: same as for data_io().
    :param interactive: same as for data_io().
    :param log_folder_paths_dict:
        Dictionary indicating paths to success and fail folders.
    :type log_folder_paths_dict: dict
    :returns:
        tuple (success_ticket_list, failed_ticket_list, success_filepath_list,
                failed_filepath_list, evaluation_dict)
        WHERE
        success_ticket_list (list) is a list of successful ImportTickets.
        failed_ticket_list (list) is a list of failed ImportTickets.
        success_filepath_list (list) is a list of successfully parsed flat files.
        failed_filepath_list (list) is a list of unsuccessfully parsed flat files.
        evaluation_dict (dict) is a dictionary from each Bundle,
            containing dictionaries for each bundled object,
            containing lists of evaluation objects.
    :rtype: tuple
    """
    # Alias for different types of genomes gathered and processed.
    file_ref = "flat_file"
    ticket_ref = "ticket"
    retain_ref = "mysql"
    retrieve_ref = "phagesdb"

    # Will hold all results data. These can be used by various types of
    # user interfaces to present the summary of import.
    success_ticket_list = []
    failed_ticket_list = []
    success_filepath_list = []
    failed_filepath_list = []
    evaluation_dict = {}

    # TODO should add a flag to set default empty values
    # for external data if the import doesn't need to rely on PhagesDB.
    # This would be for building non-Actinobacteriophage databases
    # if retrieve_phagesdb == True:
    #    external_ref_data = get_phagesdb_reference_sets()
    # else:
    #   external_ref_data = <empty data dictionary>

    # Retrieve valid cluster, subcluster, host data from PhagesDB.
    external_ref_data = get_phagesdb_reference_sets()

    # To minimize memory usage, each flat_file is evaluated one by one.
    bundle_count = 1
    file_count = 1
    for filepath in files_in_folder:
        replace_gnm_pair_key = file_ref + "_" + retain_ref
        progress = f"Processing data for file #{file_count}: {filepath.name}."
        print("\n\n" + progress)
        logger.info(progress)

        bndl = prepare_bundle(filepath=filepath, ticket_dict=ticket_dict,
                              engine=engine,
                              genome_id_field=genome_id_field,
                              host_genus_field=host_genus_field,
                              id=bundle_count,
                              file_ref=file_ref, ticket_ref=ticket_ref,
                              retrieve_ref=retrieve_ref, retain_ref=retain_ref,
                              interactive=interactive,
                              id_conversion_dict=constants.PHAGE_ID_DICT)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        # Retrieve valid data from MySQL and merge with the valid external data.
        mysql_ref_data = get_mysql_reference_sets(engine)
        ref_data = basic.merge_set_dicts(external_ref_data, mysql_ref_data)
        logger.info(f"Checking file: {filepath.name}.")
        run_checks(bndl,
                   accession_set=ref_data["accession_set"],
                   phage_id_set=ref_data["phage_id_set"],
                   seq_set=ref_data["seq_set"],
                   host_genus_set=ref_data["host_genera_set"],
                   cluster_set=ref_data["cluster_set"],
                   subcluster_set=ref_data["subcluster_set"],
                   file_ref=file_ref, ticket_ref=ticket_ref,
                   retrieve_ref=retrieve_ref, retain_ref=retain_ref)

        review_bundled_objects(bndl, interactive=interactive)

        # TODO this section below could probably be improved.
        # import_into_db may not need to return result, since it now
        # records in an Eval object whether import succeeded or not.
        # It also checks_for_errors twice (before and after attempting to add
        # data to MySQL), which can probably be simplified.
        bndl.check_for_errors()
        result = import_into_db(bndl, engine=engine,
                                gnm_key=file_ref, prod_run=prod_run)
        bndl.check_for_errors()
        dict_of_eval_lists = bndl.get_evaluations()
        logfile_path = get_logfile_path(bndl, paths_dict=log_folder_paths_dict,
                                        filepath=filepath, file_ref=file_ref)
        log_evaluations({bndl.id: dict_of_eval_lists}, logfile_path=logfile_path)
        evaluation_dict[bndl.id] = dict_of_eval_lists

        if result:
            success_ticket_list.append(bndl.ticket.data_dict)
            success_filepath_list.append(filepath)
        else:
            if bndl.ticket is not None:
                failed_ticket_list.append(bndl.ticket.data_dict)
            failed_filepath_list.append(filepath)
        bundle_count += 1
        file_count += 1

    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    if len(ticket_dict.keys()) > 0:
        logger.info("Processing unmatched tickets.")
        key_list = list(ticket_dict.keys())
        for key in key_list:
            bndl = bundle.Bundle()
            bndl.ticket = ticket_dict.pop(key)
            bndl.id = bundle_count
            logger.info("Checking data for ticket: "
                        f"{bndl.ticket.id}, {bndl.ticket.phage_id}.")
            # Since the ticket has no matching genome data,
            # no data sets from MySQL and PhagesDB are needed.
            run_checks(bndl, file_ref=file_ref, ticket_ref=ticket_ref,
                       retrieve_ref=retrieve_ref, retain_ref=retain_ref)
            bndl.check_for_errors()
            dict_of_eval_lists = bndl.get_evaluations()
            logfile_path = get_logfile_path(bndl,
                                            paths_dict=log_folder_paths_dict,
                                            filepath=None, file_ref=None)
            log_evaluations({bndl.id: dict_of_eval_lists}, logfile_path=logfile_path)
            evaluation_dict[bndl.id] = dict_of_eval_lists
            failed_ticket_list.append(bndl.ticket.data_dict)
            bundle_count += 1

    return (success_ticket_list, failed_ticket_list, success_filepath_list,
            failed_filepath_list, evaluation_dict)


def get_phagesdb_reference_sets():
    """Get multiple sets of data from PhagesDB for reference.

    :returns:
        Dictionary of unique clusters, subclusters, and host genera
        stored on PhagesDB.
    :rtype: dict
    """

    # Retrieve data from PhagesDB to create sets of
    # valid host genera, clusters, and subclusters.
    # If there is no subcluster, value may be empty string or "none".
    host_genera = phagesdb.create_host_genus_set()
    results_tuple = phagesdb.create_cluster_subcluster_sets()
    clusters = results_tuple[0]
    subclusters = results_tuple[1]
    dict = {"host_genera_set": host_genera,
            "cluster_set": clusters,
            "subcluster_set": subclusters}
    return dict


def get_mysql_reference_sets(engine):
    """Get multiple sets of data from the MySQL database for reference.

    :param engine: same as for data_io().
    :returns:
        Dictionary of unique PhageIDs, clusters, subclusters,
        host genera, accessions, and sequences stored in the MySQL database.
    :rtype: dict
    """
    phage_ids = mysqldb.get_distinct_data(engine, "phage", "PhageID")
    accessions = mysqldb.get_distinct_data(engine, "phage", "Accession")
    clusters = mysqldb.get_distinct_data(engine, "phage", "Cluster",
                                         null="Singleton")

    # Cluster "UNK" may or may not already be present, but it is valid.
    clusters.add("UNK")
    subclusters = mysqldb.get_distinct_data(engine, "phage", "Subcluster",
                                            null="none")
    host_genera = mysqldb.get_distinct_data(engine, "phage", "HostGenus")
    seqs = mysqldb.create_seq_set(engine)
    dict = {"phage_id_set": phage_ids,
            "accession_set": accessions,
            "seq_set": seqs,
            "host_genera_set": host_genera,
            "cluster_set": clusters,
            "subcluster_set": subclusters}
    return dict


def get_logfile_path(bndl, paths_dict=None, filepath=None, file_ref=None):
    """Choose the path to output the file-specific log.

    :param bndl: same as for run_checks().
    :param paths_dict:
        Dictionary indicating paths to success and fail folders.
    :type paths_dict: dict
    :param filepath: Path to flat file.
    :type filepath: Path
    :param file_ref: same as for prepare_bundle().
    :returns:
        Path to log file to store flat-file-specific evaluations.
        If paths_dict is set to None, then None is returned instead of a path.
    :rtype: Path
    """
    # Filename refers to the file that was attempted to be parsed.
    # The genome object contains the filename data only if the file was
    # successfully parsed. If it is not successfully parsed, there would
    # be no genome object to retrieve the filename.
    if paths_dict is not None:
        if bndl._errors > 0:
            log_folder_path = paths_dict[FAIL_FOLDER]
        else:
            log_folder_path = paths_dict[SUCCESS_FOLDER]

        if filepath is None:
            # If there is no filepath, then it is an unmatched ticket.
            id = bndl.ticket.phage_id
            filename = "no_file"
        else:
            # If there is a filepath, get the parsed genome id.
            filename = filepath.stem
            if file_ref in bndl.genome_dict.keys():
                gnm = bndl.genome_dict[file_ref]
                id = gnm.id
            else:
                id = "no_id"
        new_file_name = id + "__" + filename + ".log"
        logfile_path = pathlib.Path(log_folder_path, new_file_name)
    else:
        logfile_path = None
    return logfile_path


def prepare_bundle(filepath=pathlib.Path(), ticket_dict={}, engine=None,
                   genome_id_field="", host_genus_field="", id=None,
                   file_ref="", ticket_ref="", retrieve_ref="", retain_ref="",
                   id_conversion_dict={}, interactive=False):
    """Gather all genomic data needed to evaluate the flat file.

    :param filepath: Name of a GenBank-formatted flat file.
    :type filepath: Path
    :param ticket_dict: A dictionary of pdm_utils ImportTicket objects.
    :type ticket_dict: dict
    :param engine: same as for data_io().
    :param genome_id_field: same as for data_io().
    :param host_genus_field: same as for data_io().
    :param id: Identifier to be assigned to the Bundle object.
    :type id: int
    :param file_ref: Identifier for Genome objects derived from flat files.
    :type file_ref: str
    :param ticket_ref: Identifier for Genome objects derived from ImportTickets.
    :type ticket_ref: str
    :param retrieve_ref: Identifier for Genome objects derived from PhagesDB.
    :type retrieve_ref: str
    :param retain_ref: Identifier for Genome objects derived from MySQL.
    :type retain_ref: str
    :param id_conversion_dict: Dictionary of PhageID conversions.
    :type id_conversion_dict: dict
    :param interactive: same as for data_io().
    :returns:
        A pdm_utils Bundle object containing all data required to
        evaluate a flat file.
    :rtype: Bundle
    """
    bndl = bundle.Bundle()
    bndl.id = id
    seqrecord = flat_files.retrieve_genome_data(filepath)
    if seqrecord is None:
        logger.error(f"No record was retrieved from the file: {filepath.name}.")
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
        # with the ticket (which must contain the spelling used
        # in the database), the genome ID needs to be changed.
        # id_conversion_dict key = incorrect spelling; value = correct spelling
        if ff_gnm.id in id_conversion_dict.keys():
            ff_gnm.id = id_conversion_dict[ff_gnm.id]

            # Need to recompute the feature ids using the new genome id.
            ff_gnm.set_feature_ids(use_type=True, use_cds=True)
            ff_gnm.set_feature_ids(use_type=True, use_source=True)
            # TODO set tRNA and tmRNA feature ids.

            ff_gnm.set_feature_genome_ids("cds")
            ff_gnm.set_feature_genome_ids("source")
            # TODO set tRNA and tmRNA feature genome_ids.

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
            # Genome attributes from the ticket table that should not be
            # populated from the ticket, from PhagesDB, or from
            # the MySQL database, are stored in data_parse.
            # There is no need to evaluate what is stored in data_parse.

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
            # retrieve it from PhagesDB.
            # If the ticket genome has fields set to 'retrieve', data is
            # retrieved from PhagesDB and populates a new Genome object.
            if len(bndl.ticket.data_retrieve) > 0:
                pdb_gnm = phagesdb.get_genome(bndl.ticket.phage_id,
                                              gnm_type=retrieve_ref)
                if pdb_gnm is not None:
                    bndl.genome_dict[pdb_gnm.type] = pdb_gnm

                    for attr in bndl.ticket.data_retrieve:
                        attr_value = getattr(pdb_gnm, attr)
                        setattr(ff_gnm, attr, attr_value)

            # If the ticket type is 'replace', retrieve data from the MySQL database.
            # If any attributes in flat_file are set to 'retain', copy data
            # from the MySQL genome.
            if bndl.ticket.type == "replace":

                if engine is None:
                    logger.info(
                          f"Ticket {bndl.ticket.id} is a 'replace' ticket "
                          "but no details about how to connect to the "
                          "MySQL database have been provided. "
                          "Unable to retrieve data.")
                else:
                    query = "SELECT * FROM phage"
                    pmr_genomes =  mysqldb.parse_genome_data(
                                       engine=engine,
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
                                    "in the MySQL database. "
                                    "Unable to retrieve data.")

            set_cds_descriptions(ff_gnm, bndl.ticket, interactive=interactive)
    return bndl


def run_checks(bndl, accession_set=set(), phage_id_set=set(),
               seq_set=set(), host_genus_set=set(), cluster_set=set(),
               subcluster_set=set(), file_ref="", ticket_ref="",
               retrieve_ref="", retain_ref=""):
    """Run checks on the different types of data in a Bundle object.

    :param bndl: A pdm_utils Bundle object containing bundled data.
    :type bndl: Bundle
    :param accession_set: Set of accessions to check against.
    :type accession_set: set
    :param phage_id_set: Set of PhageIDs to check against.
    :type phage_id_set: set
    :param seq_set: Set of nucleotide sequences to check against.
    :type seq_set: set
    :param host_genus_set: Set of host genera to check against.
    :type host_genus_set: set
    :param cluster_set: Set of Clusters to check against.
    :type cluster_set: set
    :param subcluster_set: Set of Subclusters to check against.
    :type subcluster_set: set
    :param file_ref: same as for prepare_bundle().
    :param ticket_ref: same as for prepare_bundle().
    :param retrieve_ref: same as for prepare_bundle().
    :param retain_ref: same as for prepare_bundle().
    """
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

            # Check each type of feature.
            for x in range(len(gnm.cds_features)):
                check_cds(gnm.cds_features[x], eval_flags,
                          description_field=tkt.description_field)

            for x in range(len(gnm.trna_features)):
                # TODO make sure this tRNA is implemented correctly.
                check_trna(gnm.trna_features[x], eval_flags)

            for x in range(len(gnm.source_features)):
                check_source(gnm.source_features[x], eval_flags,
                             host_genus=gnm.host_genus)

            # TODO implement tmrna
            # for x in range(len(gnm.tmrna_features)):
            #     check_tmrna(gnm.tmrna_features[x], eval_flags)

        if retain_ref in bndl.genome_dict.keys():
            gnm2 = bndl.genome_dict[retain_ref]
            check_retain_genome(gnm2, tkt.type, eval_flags)


def review_bundled_objects(bndl, interactive=False):
    """Review all evaluations of all bundled objects.

    Iterate through all objects stored in the bundle.
    If there are warnings, review whether status should be changed.

    :param bndl: same as for run_checks().
    :param interactive: same as for data_io().
    """

    # Bundle-level evaluations.
    review_object_list([bndl], "Bundle", ["id"], interactive=interactive)

    # Ticket-level evaluations.
    review_object_list([bndl.ticket], "Ticket", ["id", "type", "phage_id"],
                       interactive=interactive)

    # Genome check.
    if len(bndl.genome_dict.keys()) > 0:
        for key in bndl.genome_dict.keys():
            gnm = bndl.genome_dict[key]

            # Genome-level evaluations.
            # Instead of passing all genomes from the genome_dict at one time,
            # it makes more sense to evaluate only one genome.evaluations list
            # at the same time as its cds_features.evaluations,
            # source.evaluations, etc.
            review_object_list([gnm], "Genome", ["id", "type"],
                               interactive=interactive)

            # CDS feature check.
            review_object_list(gnm.cds_features, "CDS feature",
                          ["id", "start", "stop", "orientation"],
                          interactive=interactive)

            # Source feature check.
            review_object_list(gnm.source_features, "Source feature",
                          ["id"], interactive=interactive)

            # # TODO implement trna features
            # review_object_list(gnm.trna_features, "tRNA feature",
            #               ["id", "start", "stop", "orientation"],
            #               interactive=interactive)
            # # TODO implement tmrna features.
            # review_object_list(gnm.tmrna_features, "tmRNA feature",
            #               ["id", "start", "stop", "orientation"],
            #               interactive=interactive)

    else:
        log_and_print("No genomes to review.", False)

    # Genome-pair check.
    genome_pair_list = list(bndl.genome_pair_dict.values())
    review_object_list(genome_pair_list, "Genome pair",
                       ["genome1","genome2"], interactive=interactive)


def review_object_list(object_list, object_type, attr_list, interactive=False):
    """Determine if evaluations are present and record results.

    :param object_list: List of pdm_utils objects containing evaluations.
    :type object_list: list
    :param object_type: Name of the pdm_utils object.
    :type object_type: str
    :param attr_list:
        List of attributes used to log data about the object instance.
    :type attr_list: list
    :param interactive: same as for data_io().
    """
    # Test for None, since tkt data can be missing and be None.
    if (len(object_list) > 0 and object_list[0] is not None):
        # Capture the exit status for each CDS feature. If user exits
        # the review at any point, skip all of the other CDS features.
        exit = False
        for x in range(len(object_list)):
            object = object_list[x]
            # Compile description of the object being reviewed.
            string = get_result_string(object, attr_list)
            partial_msg = f"evaluations for {object_type}: {string}."
            if len(object.evaluations) > 0:
                log_and_print("Reviewing " + partial_msg, interactive)
                if exit == False:
                    exit = review_evaluation_list(object.evaluations,
                                interactive=interactive)
                else:
                    # If exit=True, all 'warning' evaluations are automatically
                    # changed to 'error'. The exit response is not captured.
                    review_evaluation_list(object.evaluations, interactive=False)
            else:
                log_and_print("No " + partial_msg, False)
        if exit:
            msg = ("Not all evaluations were reviewed, "
                   "since the review was exited.")
            log_and_print(msg, False)
    else:
        log_and_print(f"No evaluations for {object_type}(s)", False)


def review_evaluation_list(evaluation_list, interactive=False):
    """Iterate through all evaluations and review 'warning' results.

    :param evaluation_list: List of pdm_utils Eval objects.
    :type evaluation_list: list
    :param interactive: same as for data_io().
    :returns: Indicates whether user selected to exit the review process.
    :rtype: bool
    """
    exit = False
    y = 0
    for x in range(len(evaluation_list)):
        evl = evaluation_list[x]
        if exit == False:
            exit, correct = review_evaluation(evl, interactive=interactive)
        else:
            # If exit=True, then all 'warning' evaluations are automatically
            # changed to 'error'. The exit2 response is unused and thrown away.
            exit2, correct = review_evaluation(evl, interactive=False)
        if not correct:
            y += 1

    if y == 0:
        log_and_print(f"No warnings or errors encountered.", interactive)
    return exit


def review_evaluation(evl, interactive=False):
    """Review an evaluation object.

    :param evl: A pdm_utils Eval object.
    :type evl: Eval
    :param interactive: same as for data_io().
    :returns:
        tuple (exit, message)
        WHERE
        exit (bool) indicates whether user selected to exit the review process.
        correct(bool) indicates whether the evalution status is accurate.
    :rtype: tuple
    """
    exit = False
    msg = " Status was changed from 'warning' to 'error' {}."
    summary = (f"Evaluation ID: {evl.id}."
               f"\nStatus: {evl.status}."
               f"\nDefinition: {evl.definition}"
               f"\nResult: {evl.result}")
    if evl.status == "warning":
        correct = False
        if interactive == True:
            # If interactive is set to True, ask user if 'warning'
            # is correct, and change the status as needed.
            print("\n\nThe following evaluation is set to 'warning':")
            print(summary)
            prompt = ("\nThis evaluation will remain as a 'warning' "
                      "(instead of it being changed to an 'error'). "
                      "Is this correct? (yes/no/exit) ")
            result = basic.ask_yes_no(prompt=prompt, response_attempt=3)
            # result can be True, False, or None
            if result != True:
                evl.status = "error"
                if result == False:
                    evl.result = evl.result + msg.format("manually")
                else:
                    exit = True
                    evl.result = evl.result + msg.format(
                                    "automatically due to review exit")
        else:
            # If interactive is set to False,
            # change all 'warnings' to 'errors'.
            evl.status = "error"
            evl.result = evl.result + msg.format(
                            "automatically due to no interactivity")
    elif evl.status == "error":
        correct = False
        if interactive == True:
            # Report that an error was encountered.
            print("\n\nThe following evaluation is set to 'error':")
            print(summary)
            input("\nPress ENTER to continue.")
    else:
        correct = True
    return exit, correct


def log_and_print(msg, terminal=False):
    """Print message to terminal in addition to logger if needed.

    :param msg: Message to print.
    :type msg: str
    :param terminal: Indicates if message should be printed to terminal.
    :type terminal: bool
    """
    logger.info(msg)
    if terminal:
        print(f"{msg}")


def get_result_string(object, attr_list):
    """Construct string of values from several object attributes.

    :param object: A object from which to retrieve values.
    :type object: misc
    :param attr_list:
        List of strings indicating attributes to retrieve from the object.
    :type attr_list: list
    :returns: A concatenated string representing values from all attributes.
    :rtype: str
    """
    strings = []
    for attr in attr_list:
        attr_value = getattr(object, attr)
        strings.append(attr + ": " + str(attr_value))
    string = ", ".join(strings)
    return string


def set_cds_descriptions(gnm, tkt, interactive=False):
    """Set the primary CDS descriptions.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt: A pdm_utils ImportTicket object.
    :type tkt: ImportTicket
    :param interactive: same as for data_io().
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

    :param feature_list: A list of pdm_utils Cds objects.
    :type feature_list: list
    :param description_field: same as for data_io().
    :returns: Name of the primary description_field after review.
    :rtype: str
    """
    logger.info("Reviewing CDS description fields.")

    # Print a table of CDS data for user to review.
    summary = []
    x = 20
    y = "..."
    tally = 0
    for cds in feature_list:
        # Get truncated descriptions to report to user.
        short_product = basic.truncate_value(cds.product, x, y)
        short_function = basic.truncate_value(cds.function, x, y)
        short_note = basic.truncate_value(cds.note, x, y)
        dict = {"CDS ID": cds.id,
                "Start": cds.start,
                "Stop": cds.stop,
                "Orientation": cds.orientation,
                "Product": short_product,
                "Function": short_function,
                "Note": short_note
                }
        summary.append(dict)

        # Count the number of descriptions that are not
        # in the selected description field.
        if (description_field != "product" and len(short_product) > 0):
            tally += 1
        if (description_field != "function" and len(short_function) > 0):
            tally += 1
        if (description_field != "note" and len(short_note) > 0):
            tally += 1


    # If there are > 0 descriptions present in unselected description fields,
    # print results and confirm with user.
    result = True
    if tally > 0:
        print(tabulate(summary, headers="keys"))
        print(f"Descriptions in the {description_field} qualifier "
              f"will be imported, but there are {tally} description(s) "
              "detected in other qualifiers.")
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
            msg = ("No new CDS description field was selected. "
                  f"The '{new_field}' CDS description field will be used.")
        else:
            msg = (f"The '{new_field}' CDS description field was selected.")

    else:
        new_field = description_field
        msg = (f"The '{new_field}' CDS description field is correct, "
               "so it will be used.")
    print(msg)
    logger.info(msg)

    # Return either the original or the new description field.
    return new_field



def check_bundle(bndl, ticket_ref="", file_ref="", retrieve_ref="", retain_ref=""):
    """Check a Bundle for errors.

    Evaluate whether all genomes have been successfully grouped,
    and whether all genomes have been paired, as expected.
    Based on the ticket type, there are expected to be certain
    types of genomes and pairs of genomes in the bundle.

    :param bndl: same as for run_checks().
    :param ticket_ref: same as for prepare_bundle().
    :param file_ref: same as for prepare_bundle().
    :param retrieve_ref: same as for prepare_bundle().
    :param retain_ref: same as for prepare_bundle().
    """
    logger.info(f"Checking bundle: {bndl.id}.")

    # Log file of genome, if applicable.
    if file_ref in bndl.genome_dict.keys():
        gnm = bndl.genome_dict[file_ref]
        logger.info(f"Genome from file: {gnm.filename}.")

    bndl.check_ticket(eval_id="BNDL_001", eval_def=EDD["BNDL_001"])
    if bndl.ticket is not None:
        logger.info(f"Ticket: {bndl.ticket.type}, {bndl.ticket.phage_id}.")
        bndl.check_genome_dict(file_ref, expect=True, eval_id="BNDL_002",
                               eval_def=EDD["BNDL_002"])
        tkt = bndl.ticket
        if len(tkt.data_add) > 0:
            bndl.check_genome_dict(ticket_ref, expect=True, eval_id="BNDL_003",
                                   eval_def=EDD["BNDL_003"])
        # There may or may not be data retrieved from PhagesDB.
        if len(tkt.data_retrieve) > 0:
            bndl.check_genome_dict(retrieve_ref,
                                   expect=True, eval_id="BNDL_004",
                                   eval_def=EDD["BNDL_004"])
        if tkt.type == "replace":
            bndl.check_genome_dict(retain_ref, expect=True, eval_id="BNDL_005",
                                   eval_def=EDD["BNDL_005"])
            # There should be a genome_pair between the current MySQL
            # genome and the new flat_file genome.
            pair_key = f"{file_ref}_{retain_ref}"
            bndl.check_genome_pair_dict(pair_key, eval_id="BNDL_006",
                                        eval_def=EDD["BNDL_006"])


def check_ticket(tkt, type_set=set(), description_field_set=set(),
        eval_mode_set=set(), id_dupe_set=set(), phage_id_dupe_set=set(),
        retain_set=set(), retrieve_set=set(), add_set=set(), parse_set=set()):
    """Evaluate a ticket to confirm it is structured appropriately.

    The assumptions for how each field is populated varies depending on
    the type of ticket.

    :param tkt: same as for set_cds_descriptions().
    :param type_set: Set of ImportTicket types to check against.
    :type type_set: set
    :param description_field_set: Set of description fields to check against.
    :type description_field_set: set
    :param eval_mode_set: Set of evaluation modes to check against.
    :type eval_mode_set: set
    :param id_dupe_set: Set of duplicated ImportTicket ids to check against.
    :type id_dupe_set: set
    :param phage_id_dupe_set: Set of duplicated ImportTicket PhageIDs to check against.
    :type phage_id_dupe_set: set
    :param retain_set: Set of retain values to check against.
    :type retain_set: set
    :param retrieve_set: Set of retrieve values to check against.
    :type retrieve_set: set
    :param add_set: Set of add values to check against.
    :type add_set: set
    :param parse_set: Set of parse values to check against.
    :type parse_set: set
    """
    # This function simply evaluates whether there is data in the
    # appropriate ticket attributes given the type of ticket.
    # It confirms that ticket attributes 'type', 'eval_mode', and
    # 'description_field' are populated with specific values.
    # But it does not evaluate the quality of the data itself for
    # the other fields, since those are genome-specific fields and
    # can be checked within Genome objects.
    # logger.info(f"Checking ticket: {tkt.id}, {tkt.type}, {tkt.phage_id}.")

    # Check for duplicated values.
    tkt.check_attribute("id", id_dupe_set,
                        expect=False, eval_id="TKT_001", eval_def=EDD["TKT_001"])
    tkt.check_attribute("phage_id", phage_id_dupe_set,
                        expect=False, eval_id="TKT_002", eval_def=EDD["TKT_002"])

    # Check these fields for specific values.
    tkt.check_attribute("type", type_set,
                        expect=True, eval_id="TKT_003", eval_def=EDD["TKT_003"])
    tkt.check_attribute("description_field", description_field_set,
                        expect=True, eval_id="TKT_004", eval_def=EDD["TKT_004"])
    tkt.check_attribute("eval_mode", eval_mode_set,
                        expect=True, eval_id="TKT_005", eval_def=EDD["TKT_005"])

    # This method may be refactored so that it accepts a list of
    # valid flag dict keys. But this has already been verified earlier
    # in the script, so it could be redundant.
    tkt.check_eval_flags(expect=True, eval_id="TKT_006", eval_def=EDD["TKT_006"])

    # For these fields, simply check that they are not empty.
    tkt.check_attribute("phage_id", {""},
                        expect=False, eval_id="TKT_007", eval_def=EDD["TKT_007"])

    # Check how genome attributes will be determined.
    tkt.check_compatible_type_and_data_retain(eval_id="TKT_008",
                                              eval_def=EDD["TKT_008"])
    tkt.check_valid_data_source("data_add", add_set,
                                eval_id="TKT_009", eval_def=EDD["TKT_009"])
    tkt.check_valid_data_source("data_retain", retain_set,
                                eval_id="TKT_010", eval_def=EDD["TKT_010"])
    tkt.check_valid_data_source("data_retrieve", retrieve_set,
                                eval_id="TKT_011", eval_def=EDD["TKT_011"])
    tkt.check_valid_data_source("data_parse", parse_set,
                                eval_id="TKT_012", eval_def=EDD["TKT_012"])


def check_genome(gnm, tkt_type, eval_flags, phage_id_set=set(),
                 seq_set=set(), host_genus_set=set(),
                 cluster_set=set(), subcluster_set=set(),
                 accession_set=set()):
    """Check a Genome object parsed from file for errors.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param eval_flags: Dictionary of boolean evaluation flags.
    :type eval_flags: dicts
    :param phage_id_set: Set of PhageIDs to check against.
    :type phage_id_set: set
    :param seq_set: Set of genome sequences to check against.
    :type seq_set: set
    :param host_genus_set: Set of host genera to check against.
    :type host_genus_set: set
    :param cluster_set: Set of clusters to check against.
    :type cluster_set: set
    :param subcluster_set: Set of subclusters to check against.
    :type subcluster_set: set
    :param accession_set: Set of accessions to check against.
    :type accession_set: set
    """
    logger.info(f"Checking genome: {gnm.id}, {gnm.type}.")

    if tkt_type == "add":
        gnm.check_attribute("id", phage_id_set | {""}, expect=False,
                            eval_id="GNM_001", eval_def=EDD["GNM_001"])
        gnm.check_attribute("name", phage_id_set | {""}, expect=False,
                            eval_id="GNM_002", eval_def=EDD["GNM_002"])
        gnm.check_attribute("seq", seq_set | {constants.EMPTY_GENOME_SEQ},
                            expect=False, eval_id="GNM_003",
                            eval_def=EDD["GNM_003"])
        gnm.check_attribute("annotation_status", {"final"}, expect=False,
                            eval_id="GNM_004", fail="warning",
                            eval_def=EDD["GNM_004"])



        # If the genome is being added, and if it has an accession,
        # no other genome is expected to have an identical accession.
        # If the genome is being replaced, and if it has an accession,
        # the prior version of the genome may or may not have had
        # accession data, so no need to check for 'replace' tickets.
        # Direct comparison of the accession from the old and new genomes
        # occurs in compare_genomes().
        if gnm.accession != "":
            gnm.check_attribute("accession", accession_set, expect=False,
                                eval_id="GNM_005", eval_def=EDD["GNM_005"])

    # 'replace' ticket checks.
    else:
        gnm.check_attribute("id", phage_id_set, expect=True, eval_id="GNM_006",
                            eval_def=EDD["GNM_006"])
        gnm.check_attribute("seq", seq_set, expect=True,
                            eval_id="GNM_007", eval_def=EDD["GNM_007"])
        gnm.check_attribute("annotation_status", {"draft"}, expect=False,
                            eval_id="GNM_008", fail="warning",
                            eval_def=EDD["GNM_008"])


    # Depending on the annotation_status of the genome,
    # CDS features are expected to contain or not contain descriptions.
    # Draft genomes should not have any descriptions.
    # Final genomes should have any descriptions.
    # There are no expectations for other types of genomes.

    # Also, Draft annotations should not have accession data.

    check_name = basic.edit_suffix(gnm.name, "add",
                                   suffix=constants.NAME_SUFFIX)

    if gnm.annotation_status == "draft":
        gnm.check_attribute("name", {check_name}, expect=True,
                            eval_id="GNM_009", fail="warning",
                            eval_def=EDD["GNM_009"])
        gnm.check_magnitude("_cds_descriptions_tally", "=", 0,
                            eval_id="GNM_010", fail="warning",
                            eval_def=EDD["GNM_010"])
        gnm.check_attribute("accession", {""}, expect=True,
                            eval_id="GNM_011", fail="warning",
                            eval_def=EDD["GNM_011"])

    else:
        if gnm.annotation_status == "final":
            gnm.check_attribute("name", {check_name}, expect=False,
                                eval_id="GNM_012", fail="warning",
                                eval_def=EDD["GNM_012"])

            if eval_flags["check_description_tally"]:
                gnm.check_magnitude("_cds_descriptions_tally", ">", 0,
                                    eval_id="GNM_013", fail="warning",
                                    eval_def=EDD["GNM_013"])

    if eval_flags["check_coords"]:
        # TODO set trna=True and tmrna=True after they are implemented.
        gnm.check_feature_coordinates(cds_ftr=True, trna_ftr=False, tmrna_ftr=False,
                                      strand=False, eval_id="GNM_030",
                                      eval_def=EDD["GNM_030"], fail="warning")


    check_id = basic.edit_suffix(gnm.name, "add", suffix=constants.NAME_SUFFIX)
    gnm.check_attribute("id", {check_id}, expect=False,
                        eval_id="GNM_014", fail="warning", eval_def=EDD["GNM_014"])
    gnm.check_attribute("annotation_status", constants.ANNOTATION_STATUS_SET,
                        expect=True, eval_id="GNM_015", eval_def=EDD["GNM_015"])
    gnm.check_attribute("annotation_author", constants.ANNOTATION_AUTHOR_SET,
                        expect=True, eval_id="GNM_016", eval_def=EDD["GNM_016"])
    gnm.check_attribute("retrieve_record", constants.RETRIEVE_RECORD_SET,
                        expect=True, eval_id="GNM_017", eval_def=EDD["GNM_017"])
    gnm.check_attribute("cluster", cluster_set, expect=True,
                        eval_id="GNM_018", fail="warning", eval_def=EDD["GNM_018"])
    gnm.check_attribute("subcluster", subcluster_set | {"none"}, expect=True,
                        eval_id="GNM_019", fail="warning", eval_def=EDD["GNM_019"])
    gnm.check_attribute("translation_table", {11}, expect=True,
                        eval_id="GNM_020", fail="warning", eval_def=EDD["GNM_020"])
    gnm.check_attribute("host_genus", host_genus_set, expect=True,
                        eval_id="GNM_021", fail="warning", eval_def=EDD["GNM_021"])

    gnm.check_cluster_structure(eval_id="GNM_022", fail="warning",
                                eval_def=EDD["GNM_022"])
    gnm.check_subcluster_structure(eval_id="GNM_023", fail="warning",
                                   eval_def=EDD["GNM_023"])
    gnm.check_compatible_cluster_and_subcluster(eval_id="GNM_024",
                                                eval_def=EDD["GNM_024"])
    gnm.check_magnitude("date", ">", constants.EMPTY_DATE, eval_id="GNM_025",
                        eval_def=EDD["GNM_025"])
    gnm.check_magnitude("gc", ">", -0.0001, eval_id="GNM_026",
                        eval_def=EDD["GNM_026"])
    gnm.check_magnitude("gc", "<", 100.0001, eval_id="GNM_027",
                        eval_def=EDD["GNM_027"])
    gnm.check_magnitude("length", ">", 0, eval_id="GNM_028",
                        eval_def=EDD["GNM_028"])
    gnm.check_magnitude("_cds_features_tally", ">", 0, eval_id="GNM_029",
                        fail="warning", eval_def=EDD["GNM_029"])

    if eval_flags["check_seq"]:
        gnm.check_nucleotides(check_set=constants.DNA_ALPHABET,
                              eval_id="GNM_031", fail="warning",
                              eval_def=EDD["GNM_031"])

    if eval_flags["check_id_typo"]:
        gnm.compare_two_attributes("id", "_description_name", expect_same=True,
                                   eval_id="GNM_032", fail="warning",
                                   eval_def=EDD["GNM_032"])
        gnm.compare_two_attributes("id", "_source_name", expect_same=True,
                                   eval_id="GNM_033", fail="warning",
                                   eval_def=EDD["GNM_033"])
        gnm.compare_two_attributes("id", "_organism_name", expect_same=True,
                                   eval_id="GNM_034", fail="warning",
                                   eval_def=EDD["GNM_034"])

    if eval_flags["check_host_typo"]:

        # TODO change these checks in parallel with source feature checks?
        # Could get a set of host genus synonyms and use check_attribute()
        # Actually, it may be better to make source feature checks parallel to
        # these genome checks. May no longer need host genus synonyms.
        # Host genus can be defined by ticket, which is set to gnm.host_genus.
        gnm.compare_two_attributes("host_genus", "_description_host_genus",
                                   expect_same=True, eval_id="GNM_035",
                                   fail="warning", eval_def=EDD["GNM_035"])
        gnm.compare_two_attributes("host_genus", "_source_host_genus",
                                   expect_same=True, eval_id="GNM_036",
                                   fail="warning", eval_def=EDD["GNM_036"])
        gnm.compare_two_attributes("host_genus", "_organism_host_genus",
                                   expect_same=True, eval_id="GNM_037",
                                   fail="warning", eval_def=EDD["GNM_037"])

    if eval_flags["check_author"]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              expect=True, eval_id="GNM_038", fail="warning",
                              eval_def=EDD["GNM_038"])
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                              expect=False, eval_id="GNM_039", fail="warning",
                              eval_def=EDD["GNM_039"])
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              expect=False, eval_id="GNM_040", fail="warning",
                              eval_def=EDD["GNM_040"])


def check_retain_genome(gnm, tkt_type, eval_flags):
    """Check a Genome object currently in database for errors.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt_type: ImportTicket type
    :type tkt_type: str
    :param eval_flags: Dictionary of boolean evaluation flags.
    :type eval_flags: dicts
    """
    logger.info(f"Checking genome: {gnm.id}, {gnm.type}.")
    if eval_flags["check_replace"]:
        gnm.check_attribute("annotation_status", {"draft"},
                            expect=True, eval_id="GNM2_001", fail="warning",
                            eval_def=EDD["GNM2_001"])


def check_source(src_ftr, eval_flags, host_genus=""):
    """Check a Source object for errors.

    :param src_ftr: A pdm_utils Source object.
    :type src_ftr: Source
    :param eval_flags: Dictionary of boolean evaluation flags.
    :type eval_flags: dicts
    :param host_genus: Host genus to check against.
    :type host_genus: str
    """
    logger.info(f"Checking source feature: {src_ftr.id}.")

    if eval_flags["check_id_typo"]:
        src_ftr.check_attribute("_organism_name", {src_ftr.genome_id},
                                expect=True, eval_id="SRC_001", fail="warning",
                                eval_def=EDD["SRC_001"])

    if eval_flags["check_host_typo"]:

        # TODO not sure if host genus synonyms is needed anymore.
        # Strategy to utilize host genus synonyms, such as 'Mycolicibacterium'.
        # host_genus_synonyms = basic.get_synonyms(host_genus,
        #                                 constants.HOST_GENUS_SYNONYMS)
        # Alternatively, a host_genus_set could be passed into the function,
        # instead of a string, and the decision to pass in a set of synonyms
        # could be pushed up a level.

        # Only evaluate attributes if the field is not empty, since they
        # are not required to be present.
        if src_ftr.organism != "":
            src_ftr.check_attribute("_organism_host_genus", {host_genus},
                                    expect=True, eval_id="SRC_002",
                                    fail="warning", eval_def=EDD["SRC_002"])
        if src_ftr.host != "":
            src_ftr.check_attribute("_host_host_genus", {host_genus},
                                    expect=True, eval_id="SRC_003",
                                    fail="warning", eval_def=EDD["SRC_003"])
        if src_ftr.lab_host != "":
            src_ftr.check_attribute("_lab_host_host_genus", {host_genus},
                                    expect=True, eval_id="SRC_004",
                                    fail="warning", eval_def=EDD["SRC_004"])


def check_cds(cds_ftr, eval_flags, description_field="product"):
    """Check a Cds object for errors.

    :param cds_ftr: A pdm_utils Cds object.
    :type cds_ftr: Cds
    :param eval_flags: Dictionary of boolean evaluation flags.
    :type eval_flags: dicts
    :param description_field: Description field to check against.
    :type description_field: str
    """
    logger.info(f"Checking CDS feature: {cds_ftr.id}.")

    cds_ftr.check_attribute("translation", {constants.EMPTY_PROTEIN_SEQ},
                            expect=False, eval_id="CDS_003",
                            eval_def=EDD["CDS_003"])
    cds_ftr.check_amino_acids(check_set=constants.PROTEIN_ALPHABET,
                              fail="warning", eval_id="CDS_001",
                              eval_def=EDD["CDS_001"])
    cds_ftr.check_translation(eval_id="CDS_002", eval_def=EDD["CDS_002"])
    cds_ftr.check_attribute("translation_table", {11},
                            expect=True, eval_id="CDS_004", fail="warning",
                            eval_def=EDD["CDS_004"])
    # Start can be 0 since coordinate_format = 0_half_open
    cds_ftr.check_magnitude("start", ">", -1, eval_id="CDS_005",
                            eval_def=EDD["CDS_005"])
    cds_ftr.check_magnitude("stop", ">", 0, eval_id="CDS_013",
                            eval_def=EDD["CDS_013"])
    cds_ftr.check_magnitude("parts", ">", 0, eval_id="CDS_014",
                            eval_def=EDD["CDS_014"])
    cds_ftr.check_orientation(format="fr_short", case=True, eval_id="CDS_006",
                              eval_def=EDD["CDS_006"])

    if eval_flags["check_locus_tag"]:
        cds_ftr.check_attribute("locus_tag", {""}, expect=False, eval_id="CDS_007",
                                fail="warning", eval_def=EDD["CDS_007"])

        # TODO this check could be improved to take into account the prefix.
        cds_ftr.check_locus_tag_structure(check_value=None, only_typo=True,
            case=True, eval_id="CDS_008", fail="warning", eval_def=EDD["CDS_008"])
    if eval_flags["check_gene"]:
        cds_ftr.check_attribute("gene", {""}, expect=False, eval_id="CDS_009",
                                fail="warning", eval_def=EDD["CDS_009"])
        cds_ftr.check_gene_structure(eval_id="CDS_010", fail="warning",
                                     eval_def=EDD["CDS_010"])
    if (eval_flags["check_locus_tag"] and eval_flags["check_gene"]):
        cds_ftr.check_compatible_gene_and_locus_tag(eval_id="CDS_011",
                                                    fail="warning",
                                                    eval_def=EDD["CDS_011"])

    # if eval_flags["check_description"]:
        # TODO the "check_generic_data" method should be implemented at the genome level.
        # cds_ftr.check_generic_data(eval_id="CDS_012")
        # TODO not implemented yet: cds_ftr.check_valid_description(eval_id="CDS_013")
    if eval_flags["check_description_field"]:
        cds_ftr.check_description_field(attribute=description_field,
                                        eval_id="CDS_012", fail="warning",
                                        eval_def=EDD["CDS_012"])

def compare_genomes(genome_pair, eval_flags):
    """Compare two genomes to identify discrepancies.

    :param genome_pair: A pdm_utils GenomePair object.
    :type genome_pair: GenomePair
    :param eval_flags: Dictionary of boolean evaluation flags.
    :type eval_flags: dicts
    """
    logger.info("Comparing data between two genomes: "
                f"Genome 1. ID: {genome_pair.genome1.id}. "
                f"Type: {genome_pair.genome1.type}. "
                f"Genome 2. ID: {genome_pair.genome2.id}. "
                f"Type: {genome_pair.genome2.type}."
                )

    genome_pair.compare_attribute("id", expect_same=True, eval_id="GP_001",
                                  fail="warning", eval_def=EDD["GP_001"])
    genome_pair.compare_attribute("seq", expect_same=True, eval_id="GP_002",
                                  fail="warning", eval_def=EDD["GP_002"])
    genome_pair.compare_attribute("length", expect_same=True, eval_id="GP_003",
                                  fail="warning", eval_def=EDD["GP_003"])
    genome_pair.compare_attribute("cluster", expect_same=True, eval_id="GP_004",
                                  fail="warning", eval_def=EDD["GP_004"])
    genome_pair.compare_attribute("subcluster", expect_same=True,
                                  eval_id="GP_005", fail="warning",
                                  eval_def=EDD["GP_005"])
    genome_pair.compare_attribute("host_genus", expect_same=True,
                                  eval_id="GP_006", fail="warning",
                                  eval_def=EDD["GP_006"])
    genome_pair.compare_attribute("annotation_author", expect_same=True,
                                  eval_id="GP_007", fail="warning",
                                  eval_def=EDD["GP_007"])
    genome_pair.compare_attribute("translation_table", expect_same=True,
                                  eval_id="GP_008", fail="warning",
                                  eval_def=EDD["GP_008"])
    genome_pair.compare_attribute("retrieve_record", expect_same=True,
                                  eval_id="GP_009", fail="warning",
                                  eval_def=EDD["GP_009"])

    if eval_flags["check_replace"]:
        # The following checks assume that:
        # 'genome1' slot = new genome to be evaluated
        # 'genome2' slot = the current MySQL genome

        # The new genome to be evaluated is expected to be
        # newer than the current MySQL genome annotations.
        genome_pair.compare_date("newer", eval_id="GP_010", fail="warning",
                                 eval_def=EDD["GP_010"])

        if genome_pair.genome2.annotation_status == "draft":
            # It is expected that the replacing genome name no longer
            # retains the "_Draft" suffix, so the name should change.
            genome_pair.compare_attribute("name", expect_same=False,
                                          eval_id="GP_011", fail="warning",
                                          eval_def=EDD["GP_011"])

            # The status should change from 'draft'.
            genome_pair.compare_attribute("annotation_status", expect_same=False,
                                          eval_id="GP_012", fail="warning",
                                          eval_def=EDD["GP_012"])
        else:
            genome_pair.compare_attribute("name", expect_same=True,
                                          eval_id="GP_013", fail="warning",
                                          eval_def=EDD["GP_013"])

            # This is tricky. Yes, if replacing, you only expect a
            # final -> final, or unknown -> unknown. However, if the
            # check_replace = True, such as when a new final is available from
            # PhagesDB, then you only expect to go from draft -> final,
            # and you don't expect the current MySQL genome to be final.
            genome_pair.compare_attribute("annotation_status", expect_same=True,
                                          eval_id="GP_014", fail="warning",
                                          eval_def=EDD["GP_014"])
            genome_pair.compare_attribute("accession", expect_same=True,
                                          eval_id="GP_015", fail="warning",
                                          eval_def=EDD["GP_015"])




# TODO implement.
# TODO unit test.
def check_trna(trna_obj, eval_flags):
    """Check a TrnaFeature object for errors."""

    if eval_flags["check_trna"]:
        pass
    else:
        pass


def import_into_db(bndl, engine=None, gnm_key="", prod_run=False):
    """Import data into the MySQL database.

    :param bndl: same as for run_checks().
    :param engine: same as for data_io().
    :param gnm_key:
        Identifier for the Genome object in the Bundle's genome dictionary.
    :type gnm_key: str
    :param prod_run: same as for data_io().
    """
    if bndl._errors == 0:
        import_gnm = bndl.genome_dict[gnm_key]

        # If locus_tag data should not be imported, clear all locus_tag data.
        if bndl.ticket.eval_flags["import_locus_tag"] is False:
            logger.info("Clearing locus_tag data for "
                        f"genome: {import_gnm.id}.")
            import_gnm.clear_locus_tags()

        # Update the date field to reflect the day of import.
        import_gnm.date = IMPORT_DATE
        bndl.sql_statements = mysqldb.create_genome_statements(
                                import_gnm, bndl.ticket.type)
        if prod_run:
            logger.info("Importing data into the database for "
                        f"genome: {import_gnm.id}.")
            execute_result, msg = mysqldb.execute_transaction(engine,
                                        bndl.sql_statements)
            if execute_result == 1:
                result = False
                logger.error("Error importing data. " + msg)
            else:
                result = True
                logger.info("Data successfully imported. " + msg)
                logger.info("The following MySQL statements were executed:")
                for statement in bndl.sql_statements:
                    statement = basic.truncate_value(statement, 150, "...")
                    logger.info(statement)

            # Result of statement execution is stored in an eval object
            # so that it can be recorded in the log file.
            # It is stored in the bundle object instead of the genome object
            # since the collection of MySQL statements are constructed
            # based on the ticket data and the genome data combined,
            # and not just the genome data.
            bndl.check_statements(execute_result, msg, eval_id="BNDL_007",
                                  eval_def=EDD["BNDL_007"])
        else:
            result = True
            logger.info("Data can be imported if set to production run.")
    else:
        result = False
        logger.info("Data contains errors, so it will not be imported.")

    return result
