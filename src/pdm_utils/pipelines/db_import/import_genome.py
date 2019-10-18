"""Primary pipeline to process and evaluate data to be imported
into PhameratorDB."""


import time
from datetime import datetime
import csv
import os
import sys
import argparse
import pathlib
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
    parser = argparse.ArgumentParser(description=IMPORT_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("input_folder", type=os.path.abspath,
        help=INPUT_FOLDER_HELP)
    parser.add_argument("import_table", type=os.path.abspath,
        help=IMPORT_TABLE_HELP)
    parser.add_argument("-gf", "--genome_id_field", type=str,
        default="organism_name", choices=["organism_name", "filename"],
        help=GENOME_ID_FIELD_HELP)
    parser.add_argument("-p", "--prod_run", action="store_true", default=False,
        help=PROD_RUN_HELP)
    parser.add_argument("-rm", "--run_mode", type=str.lower,
        choices=list(run_modes.RUN_MODES.keys()), default="phagesdb",
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

    args.input_folder = pathlib.Path(args.input_folder)
    args.input_folder.expanduser()
    args.input_folder.resolve()
    if not args.input_folder.is_dir():
        print("Invalid input folder.")
        sys.exit(1)

    args.import_table = pathlib.Path(args.import_table)
    args.import_table.expanduser()
    args.import_table.resolve()
    if not args.import_table.is_file():
        print("Invalid import table file.")
        sys.exit(1)

    # If everything checks out, pass args to the main import pipeline:
    # TODO add output_folder parameter and args.output_folder.
    setup(sql_handle=sql_handle,
        genome_folder=args.input_folder, import_table_file=args.import_table,
        genome_id_field=args.genome_id_field, prod_run=args.prod_run,
        description_field=args.description_field, run_mode=args.run_mode,
        output_folder="")




# TODO unittest.
def setup(sql_handle=None, genome_folder=pathlib.Path(), import_table_file="",
    genome_id_field="", prod_run=False, description_field="", run_mode="",
    output_folder=pathlib.Path()):
    """Set up output directories, log files, etc. for import."""
    # Create output directories
    date = time.strftime("%Y%m%d")
    results_folder = "%s_import_results" % date
    results_path = basic.make_new_dir(output_folder, results_folder)
    if results_path is None:
        print("\nUnable to create output_folder.")
        sys.exit(1)

    # Get the files to process.
    files_to_process = basic.identify_files(genome_folder, set([".DS_Store"]))
    if len(files_to_process) == 0:
        print("There are no flat files to evaluate.")
        sys.exit(1)

    # Get the tickets.
    eval_flags = run_modes.get_eval_flag_dict(run_mode.lower())
    run_mode_eval_dict = {"run_mode": run_mode, "eval_flag_dict": eval_flags}
    ticket_dict = prepare_tickets(
                    import_table_file, run_mode_eval_dict, description_field,
                    required_keys=constants.IMPORT_TABLE_REQ_DICT.keys(),
                    optional_keys=constants.IMPORT_TABLE_OPT_DICT.keys(),
                    keywords=set(["retrieve", "retain", "none"]))

    if ticket_dict is None:
        print("Invalid import table. Unable to evaluate flat files.")
        sys.exit(1)

    # Evaluate files and tickets.
    results_tuple = process_files_and_tickets(ticket_dict, files_to_process,
                        sql_handle, prod_run, genome_id_field)

    success_ticket_list = results_tuple[0]
    failed_ticket_list = results_tuple[1]
    success_filename_list = results_tuple[2]
    failed_filename_list = results_tuple[3]
    evaluation_dict = results_tuple[4]

    # Prepare output
    failed_path = pathlib.Path(results_path, "fail")
    success_path = pathlib.Path(results_path, "success")
    failed_path.mkdir()
    success_path.mkdir()

    order = list(required_keys) + list(optional_keys)
    if len(success_ticket_list) > 0:
        # TODO this should be improved to retain a specific column order.
        success_tkt_file = pathlib.Path(success_path, "import_tickets.csv")
        with open(success_tkt_file, "w") as file_handle:
            for tkt_dict in success_ticket_list:
                file_writer = csv.DictWriter(file_handle, order)

    if len(failed_ticket_list) > 0:
        failed_tkt_file = pathlib.Path(failed_path, "import_tickets.csv")
        with open(failed_tkt_file, "w") as file_handle:
            for tkt_dict in failed_ticket_list:
                file_writer = csv.DictWriter(file_handle, order)

    if len(success_filename_list) > 0:
        success_genomes_path = pathlib.Path(success_path, "genomes")
        success_genomes_path.mkdir()
        for file in success_filename_list:
            new_file = Path(success_genomes_path, file.name)
            shutil.move(str(file), str(new_file))

    if len(failed_filename_list) > 0:
        failed_genomes_path = pathlib.Path(failed_path, "genomes")
        failed_genomes_path.mkdir()
        for file in failed_filename_list:
            new_file = Path(failed_genomes_path, file.name)
            shutil.move(str(file), str(new_file))

    # TODO write out evaluations to some sort of log file.
    for bundle_id in evaluation_dict:
        print("Evaluations for bundle {}:".format(bundle_id))
        bundle_evals = evaluation_dict[bundle_id]
        for key in bundle_evals:
            print("Evaluations for {}".format(key))
            evl_list = bundle_evals[key]
            for evl in evl_list:
                print("Evaluation: {}. Definition: {}. Status: {}. Result: {}.".format( \
                    evl.id, evl.definition, evl.status, evl.result))


    print("Import complete.")


def prepare_tickets(import_table_file="", run_mode_eval_dict=None,
        description_field="", required_keys=set(), optional_keys=set(),
        keywords=set()):
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

    list_of_tkts = []
    tkt_errors = 0
    list_of_data_dicts = tickets.retrieve_ticket_data(import_table_file)
    list_of_tkts = tickets.construct_tickets(list_of_data_dicts, run_mode_eval_dict,
                    description_field, required_keys, optional_keys,
                    keywords)
    if len(list_of_data_dicts) != len(list_of_tkts):
        tkt_errors += 1
    tkt_id_dupes, phage_id_dupes = tickets.identify_duplicates(list_of_tkts)
    ticket_dict = {}

    x = 0
    while x < len(list_of_tkts):
        tkt = list_of_tkts[x]
        check_ticket(
            tkt,
            type_set=constants.IMPORT_TICKET_TYPE_SET,
            description_field_set=constants.DESCRIPTION_FIELD_SET,
            null_set=constants.EMPTY_SET,
            run_mode_set=run_modes.RUN_MODES.keys(),
            id_dupe_set=tkt_id_dupes,
            phage_id_dupe_set=phage_id_dupes)
        for evl in tkt.evaluations:
            if evl.status == "error":
                tkt_errors += 1
        ticket_dict[tkt.phage_id] = tkt
        x += 1

    if tkt_errors > 0:
        print("Error generating tickets from import table.")
        return None
    else:
        return ticket_dict

def process_files_and_tickets(ticket_dict, files_in_folder, sql_handle=None,
                              prod_run=False, genome_id_field=""):
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
    :param prod_run:
        Indicates whether the database should be updated from the
        import tickets.
    :type prod_run: bool
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

    # Retrieve data from phagesdb to create sets of
    # valid host genera, clusters, and subclusters.
    # Cluster "UNK" may or may not already be present, but it is valid.
    phagesdb_host_genera_set = phagesdb.create_host_genus_set()
    results_tuple = phagesdb.create_cluster_subcluster_sets()
    phagesdb_cluster_set = results_tuple[0]
    phagesdb_subcluster_set = results_tuple[1]
    phagesdb_cluster_set.add("UNK")

    # To minimize memory usage, each flat_file is evaluated one by one.
    bundle_count = 1
    for filename in files_in_folder:
        gnm_type = "flat_file"
        replace_gnm_pair_key = gnm_type + "_phamerator"
        bndl = prepare_bundle(filename=str(filename), ticket_dict=ticket_dict,
                              sql_handle=sql_handle,
                              genome_id_field=genome_id_field, id=bundle_count)

        # Create sets of unique values for different data fields.
        # Since data from each parsed flat file is imported into the
        # database one file at a time, these sets are not static.
        # So these sets should be recomputed for every flat file evaluated.
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)
        run_checks(bndl,
                  null_set=set(),
                  accession_set=phamerator_accession_set,
                  phage_id_set=phamerator_phage_id_set,
                  seq_set=phamerator_seq_set,
                  host_genus_set=phagesdb_host_genera_set,
                  cluster_set=phagesdb_cluster_set,
                  subcluster_set=phagesdb_subcluster_set,
                  gnm_key=gnm_type,
                  gnm_pair_key=replace_gnm_pair_key)
        evaluation_dict[bndl.id] = bndl.get_evaluations()
        result = import_into_db(bndl, sql_handle=sql_handle,
                                gnm_key=gnm_type, prod_run=prod_run)
        if result:
            success_ticket_list.append(bndl.ticket.data_dict)
            success_filename_list.append(filename)
        else:
            if bndl.ticket is not None:
                failed_ticket_list.append(bndl.ticket.data_dict)
            failed_filename_list.append(filename)
        bundle_count += 1

    # Tickets were popped off the ticket dictionary as they were matched
    # to flat files. If there are any tickets left, errors need to be counted.
    if len(ticket_dict.keys()) > 0:
        phamerator_phage_id_set = phamerator.create_phage_id_set(sql_handle)
        phamerator_seq_set = phamerator.create_seq_set(sql_handle)
        phamerator_accession_set = phamerator.create_accession_set(sql_handle)
        key_list = list(ticket_dict.keys())
        for key in key_list:
            bndl = bundle.Bundle()
            bndl.ticket = ticket_dict.pop(key)
            bndl.id = bundle_count
            run_checks(bndl,
                      null_set=set(),
                      accession_set=phamerator_accession_set,
                      phage_id_set=phamerator_phage_id_set,
                      seq_set=phamerator_seq_set,
                      host_genus_set=phagesdb_host_genera_set,
                      cluster_set=phagesdb_cluster_set,
                      subcluster_set=phagesdb_subcluster_set,
                      gnm_key=None, gnm_pair_key=None)
            evaluation_dict[bndl.id] = bndl.get_evaluations()
            failed_ticket_list.append(bndl.ticket.data_dict)
            bundle_count += 1

    return (success_ticket_list, failed_ticket_list, success_filename_list,
            failed_filename_list, evaluation_dict)


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
            if len(bndl.ticket.data_ticket) > 0:
                tkt_gnm = tickets.get_genome(bndl.ticket, gnm_type="ticket")
                bndl.genome_dict[tkt_gnm.type] = tkt_gnm

                # Copy ticket data to flat file. Since the ticket data has been
                # added to a genome object using genome methods, the data
                # can be directly passed from one genome object to another.
                for attr in bndl.ticket.data_ticket:
                    attr_value = getattr(tkt_gnm, attr)
                    setattr(ff_gnm, attr, attr_value)

            # Check to see if there is any missing data for each genome, and
            # retrieve it from phagesdb.
            # If the ticket genome has fields set to 'retrieve', data is
            # retrieved from PhagesDB and populates a new Genome object.
            if len(bndl.ticket.data_retrieve) > 0:
                pdb_gnm = phagesdb.get_genome(bndl.ticket.phage_id,
                                              gnm_type="phagesdb")
                bndl.genome_dict[pdb_gnm.type] = pdb_gnm

                for attr in bndl.ticket.data_retrieve:
                    attr_value = getattr(pdb_gnm, attr)
                    setattr(ff_gnm, attr, attr_value)

            # If the ticket type is 'replace', retrieve data from phamerator.
            # If any attributes in flat_file are set to 'retain', copy data
            # from the phamerator genome.
            if bndl.ticket.type == "replace":

                if sql_handle is None:
                    print("Ticket %s is a 'replace' ticket but no "
                          "details about how to connect to the "
                          "Phamerator database have been provided. "
                          "Unable to retrieve data."
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

                        # The ticket may indicate some data should be retained.
                        for attr in bndl.ticket.data_retain:
                            attr_value = getattr(pmr_gnm, attr)
                            setattr(ff_gnm, attr, attr_value)

                        # Pair the genomes for comparison evaluations.
                        gnm_pair = genomepair.GenomePair()
                        bndl.set_genome_pair(gnm_pair, ff_gnm.type, pmr_gnm.type)
                    else:
                        print("There is no %s genome in the Phamerator"
                              " database. Unable to retrieve data."
                              % ff_gnm.id)
    return bndl


def run_checks(bndl, null_set=set(), accession_set=set(), phage_id_set=set(),
               seq_set=set(), host_genus_set=set(), cluster_set=set(),
               subcluster_set=set(), gnm_key=None, gnm_pair_key=None):
    """Run checks on the different elements of a bundle object."""
    check_bundle(bndl)
    tkt = bndl.ticket
    if tkt is not None:
        eval_flags = tkt.eval_flags
        if (tkt.type == "replace" and \
                gnm_pair_key in bndl.genome_pair_dict.keys()):
            genome_pair = bndl.genome_pair_dict[gnm_pair_key]
            compare_genomes(genome_pair, eval_flags)

        if gnm_key in bndl.genome_dict.keys():
            gnm = bndl.genome_dict[gnm_key]
            check_genome(gnm, tkt.type, eval_flags, null_set=null_set,
                         accession_set=accession_set, phage_id_set=phage_id_set,
                         seq_set=seq_set, host_genus_set=host_genus_set,
                         cluster_set=cluster_set, subcluster_set=subcluster_set)

            # Check CDS features.
            x = 0
            while x < len(gnm.cds_features):
                check_cds(gnm.cds_features[x], eval_flags)
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
                check_source(gnm.source_features[z], eval_flags)
                z += 1
    bndl.check_for_errors()


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
        if len(tkt.data_ticket) > 0:
            bndl.check_genome_dict("ticket", eval_id="BNDL_002")

        bndl.check_genome_dict("flat_file", eval_id="BNDL_003")

        # There may or may not be data retrieved from PhagesDB.
        if len(tkt.data_retrieve) > 0:
            bndl.check_genome_dict("phagesdb", eval_id="BNDL_004")

        if tkt.type == "replace":
            bndl.check_genome_dict("phamerator", eval_id="BNDL_005")

            # There should be a genome_pair between the current phamerator
            # genome and the new flat_file genome.
            bndl.check_genome_pair_dict("flat_file_phamerator",
                                        eval_id="BNDL_006")


def check_ticket(tkt, type_set=set(), description_field_set=set(),
        null_set=set(), run_mode_set=set(), id_dupe_set=set(),
        phage_id_dupe_set=set()):
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

    # Check these fields for specific values.
    tkt.check_type(type_set, True, eval_id="TKT_003")
    tkt.check_description_field(description_field_set, True, eval_id="TKT_004")
    tkt.check_run_mode(run_mode_set, True, eval_id="TKT_005")
    tkt.check_eval_flags(True, eval_id="TKT_006")

    # For these fields, simply check that they are not empty.
    tkt.check_phage_id(null_set, False, eval_id="TKT_007")

    # Check if certain combinations of fields make sense.
    tkt.check_compatible_type_and_annotation_status(eval_id="TKT_008")
    tkt.check_compatible_type_and_data_retain(eval_id="TKT_009")


def check_genome(gnm, tkt_type, eval_flags, null_set=set(), phage_id_set=set(),
                 seq_set=set(), host_genus_set=set(),
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

    if tkt_type == "add":
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
    if eval_flags["check_seq"]:
        gnm.check_nucleotides(check_set=constants.DNA_ALPHABET,
                              eval_id="GNM_018")
    gnm.check_compatible_status_and_accession(eval_id="GNM_019")
    gnm.check_compatible_status_and_descriptions(eval_id="GNM_020")
    if eval_flags["check_id_typo"]:
        gnm.check_description_name(eval_id="GNM_021")
        gnm.check_source_name(eval_id="GNM_022")
        gnm.check_organism_name(eval_id="GNM_023")
    if eval_flags["check_host_typo"]:
        gnm.check_host_genus(host_genus_set, True, eval_id="GNM_024")
        gnm.check_description_host_genus(eval_id="GNM_025")
        gnm.check_source_host_genus(eval_id="GNM_026")
        gnm.check_organism_host_genus(eval_id="GNM_027")
    if eval_flags["check_author"]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              eval_id="GNM_028")
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                                 expect=False, eval_id="GNM_029")
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET, expect=False,
                              eval_id="GNM_030")

    gnm.check_cds_feature_tally(eval_id="GNM_031")

    # TODO set trna and tmrna to True after they are implemented.
    gnm.check_feature_ids(cds_ftr=True, trna_ftr=False, tmrna=False,
                          eval_id="GNM_032")

    # TODO the following checks may no longer be needed, since the
    # copying functions have been changed.
    # TODO confirm that these check_value_flag() are needed here.
    # Currently all "copy_data" functions run the check method
    # to throw an error if not all data was copied.
    gnm.set_value_flag("retrieve")
    gnm.check_value_flag(eval_id="GNM_033")
    gnm.set_value_flag("retain")
    gnm.check_value_flag(eval_id="GNM_034")


def check_source(src_ftr, eval_flags):
    """Check a Source object for errors."""
    if eval_flags["check_id_typo"]:
        src_ftr.check_organism_name(eval_id="SRC_001")
    if eval_flags["check_host_typo"]:
        src_ftr.check_organism_host_genus(eval_id="SRC_002")
        src_ftr.check_host_host_genus(eval_id="SRC_003")
        src_ftr.check_lab_host_host_genus(eval_id="SRC_004")


def check_cds(cds_ftr, eval_flags):
    """Check a Cds object for errors."""
    cds_ftr.check_amino_acids(eval_id="CDS_001")
    cds_ftr.check_translation(eval_id="CDS_002")
    cds_ftr.check_translation_length(eval_id="CDS_003")
    cds_ftr.check_translation_table(eval_id="CDS_004")
    cds_ftr.check_coordinates(eval_id="CDS_005")
    cds_ftr.check_strand(eval_id="CDS_006")
    if eval_flags["check_locus_tag"]:
        cds_ftr.check_locus_tag_present(eval_id="CDS_007")
        cds_ftr.check_locus_tag_structure(eval_id="CDS_008")
    if eval_flags["check_gene"]:
        cds_ftr.check_gene_present(eval_id="CDS_009")
        cds_ftr.check_gene_structure(eval_id="CDS_010")
    if (eval_flags["check_locus_tag"] and eval_flags["check_gene"]):
        cds_ftr.check_compatible_gene_and_locus_tag(eval_id="CDS_011")
    if eval_flags["check_description"]:
        cds_ftr.check_generic_data(eval_id="CDS_012")
        cds_ftr.check_valid_description(eval_id="CDS_013")
    if eval_flags["check_description_field"]:
        cds_ftr.check_description_field(eval_id="CDS_014")


def compare_genomes(genome_pair, eval_flags):
    """Compare two genomes to identify discrepancies."""
    genome_pair.compare_genome_sequence(eval_id="GP_001")
    genome_pair.compare_genome_length(eval_id="GP_002")
    genome_pair.compare_cluster(eval_id="GP_003")
    genome_pair.compare_subcluster(eval_id="GP_004")
    genome_pair.compare_accession(eval_id="GP_005")
    genome_pair.compare_host_genus(eval_id="GP_006")
    genome_pair.compare_author(eval_id="GP_007")
    if eval_flags["check_replace"]:
        genome_pair.compare_annotation_status("type","phamerator",
            "flat_file","draft","final", eval_id="GP_008")


# TODO implement.
# TODO unit test.
def check_trna(trna_obj, eval_flags):
    """Check a TrnaFeature object for errors."""

    pass


def import_into_db(bndl, sql_handle=None, gnm_key="", prod_run=False):
    """Import data into the MySQL database."""
    if bndl._errors == 0:
        phamerator.create_genome_statements(bndl.genome_dict[gnm_key],
                                            bndl.ticket.type)
        if prod_run:
            result = sql_handle.execute_transaction(bndl.sql_queries)
            if result == 1:
                # TODO log any errors.
                print("Error executing statements to import data.")
                result = False
            else:
                result = True
                print("Data successfully imported.")
        else:
            result = True
            print("Data can be imported if set to production run.")
    else:
        result = False
        print("Data contains errors, so it will not be imported.")
    return result

























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
