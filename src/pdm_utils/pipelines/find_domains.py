import argparse
import logging
import os
import pathlib
import platform
import shlex
import sys
from subprocess import Popen, PIPE # import warnings

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbirpsblastCommandline
import sqlalchemy

import pdm_utils
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions.basic import expand_path
from pdm_utils.functions.parallelize import *

# SQL QUERIES
GET_GENES_FOR_CDD = (
    "SELECT GeneID, CONVERT(Translation USING utf8) as Translation "
    "FROM gene WHERE DomainStatus = 0")
GET_UNIQUE_HIT_IDS = "SELECT HitID FROM domain"

# SQL COMMANDS
INSERT_INTO_DOMAIN = (
    """INSERT INTO domain (HitID, DomainID, Name, Description) """
    """VALUES ("{}", "{}", "{}", "{}")""")

INSERT_INTO_GENE_DOMAIN = (
    """INSERT INTO gene_domain """
    """(GeneID, HitID, Expect, QueryStart, QueryEnd) """
    """VALUES ("{}", "{}", {}, {}, {})""")

UPDATE_GENE = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"

CLEAR_GENE_DOMAIN = "TRUNCATE gene_domain"
CLEAR_DOMAIN = "DELETE FROM domain"
CLEAR_GENE_DOMAINSTATUS = "UPDATE gene SET DomainStatus = 0"



# MISC
VERSION = pdm_utils.__version__

# TODO tmp dir stores rpsblast output, while output folder stores logged info
# about the find_domains pipeline. These two directories could possibly be
# combined, so that all tmp dir automatically gets created within output folder
# beside the results folder, or something like that.
DEFAULT_TMP_DIR = "/tmp/find_domains"
DEFAULT_OUTPUT_FOLDER = os.getcwd()
RESULTS_FOLDER = f"{constants.CURRENT_DATE}_find_domains"
MAIN_LOG_FILE = "find_domains.log"
DEFAULT_CDD = "~/Databases/Cdd_LE"

# LOGGING
# Add a logger named after this module. Then add a null handler, which
# suppresses any output statements. This allows other modules that call this
# module to define the handler and output formats. If this module is the
# main module being called, the top level main function configures
# the root logger and primary file handle.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def setup_argparser():
    """
    Builds argparse.ArgumentParser for this script
    :return:
    """
    # Pipeline description
    description = (
        "Uses rpsblast to search the NCBI conserved domain database "
        "for significant domain hits in all new proteins of a "
        "MySQL database.")
    output_folder_help = (
        "Directory where log data can be generated.")
    reset_help = (
        "Clear all domain data in the database before finding domains.")
    cdd_help = (
        "Path to local NCBI Conserved Domain Database. "
        f"Default is {DEFAULT_CDD}.")
    config_file_help = (
        "Path to the file containing user-specific login data.")

    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("database", type=str,
                        help="name of database to find conserved domains for")
    parser.add_argument("-d", "--cdd", type=str, default=DEFAULT_CDD,
                        help=cdd_help)
    parser.add_argument("-t", "--threads", default=mp.cpu_count(), type=int,
                        help="number of concurrent CDD searches to run")
    parser.add_argument("-e", "--evalue", default=0.001, type=float,
                        help="evalue cutoff for rpsblast(+) hits")
    parser.add_argument("-td", "--tmp_dir", default=DEFAULT_TMP_DIR, type=str,
                        help="path to temporary directory for file I/O")
    parser.add_argument("-r", "--rpsblast", default="", type=str,
                        help="path to rpsblast(+) binary")
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
                        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                        help=output_folder_help)
    parser.add_argument("-b", "--batch_size", default=10000, type=int,
                        help="number of translations to search at a time")
    parser.add_argument("-x", "--reset", action="store_true",
                        default=False, help=reset_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)

    return parser


def make_tempdir(tmp_dir):
    """
    Uses pdm_utils.functions.basic.expand_path to expand TMP_DIR; then
    checks whether tmpdir exists - if it doesn't, uses os.makedirs to
    make it recursively.
    :param tmp_dir: location where I/O should take place
    :return:
    """
    try:
        path = expand_path(tmp_dir)
        # If the path doesn't exist yet
        if not os.path.exists(path):
            # Make it recursively
            os.makedirs(path)
    except OSError as err:
        print(f"Error {err.args[0]}: {err.args[1]}")


def search_and_process(rpsblast, cdd_name, tmp_dir, evalue,
                       translation_id, translation):
    """
    Uses rpsblast to search indicated gene against the indicated CDD
    :param rpsblast: path to rpsblast binary
    :param cdd_name: CDD database path/name
    :param tmp_dir: path to directory where I/O will take place
    :param evalue: evalue cutoff for rpsblast
    :param translation_id: unique identifier for the translation sequence
    :param translation: protein sequence for gene to query
    :return: results
    """
    # Currently translation_id is used only for file formatting.
    # Setup I/O variables
    i = "{}/{}.txt".format(tmp_dir, translation_id)
    o = "{}/{}.xml".format(tmp_dir, translation_id)

    # Write the input file
    with open(i, "w") as fh:
        fh.write(">{}\n{}".format(translation_id, translation))

    # Setup, run the rpsblast command, and process results
    rps_command = NcbirpsblastCommandline(cmd=rpsblast, db=cdd_name,
                                          query=i, out=o, outfmt=5,
                                          evalue=evalue)
    rps_command()
    data = process_rps_output(o, evalue)

    # Currently need to return as a list due to a
    # filter in parallelize.start_processes()
    results = [{"Translation": translation, "Data": data}]
    return results


def process_rps_output(filepath, evalue):
    """Process rpsblast output and return list of dictionaries."""
    results = []
    with open(filepath, "r") as fh:
        for record in NCBIXML.parse(fh):
            for align in record.alignments:
                des, d_id, name = process_align(align)
                for hsp in align.hsps:
                    if hsp.expect <= evalue:
                        dict = {"HitID": align.hit_id,
                                "DomainID": d_id,
                                "Name": name,
                                "Description": des,
                                "Expect": float(hsp.expect),
                                "QueryStart": int(hsp.query_start),
                                "QueryEnd": int(hsp.query_end)}
                        results.append(dict)
    return results


def process_align(align):
    """Process alignment data.

    Returns description, domain_id, and name.
    """
    align.hit_def = align.hit_def.replace("\"", "\'")
    des_list = align.hit_def.split(",")
    if len(des_list) == 1:
        description = des_list[0].strip()
        domain_id = None
        name = None
    elif len(des_list) == 2:
        domain_id = des_list[0].strip()
        description = des_list[1].strip()
        name = None
    else:
        domain_id = des_list[0].strip()
        name = des_list[1].strip()
        # Name is occassionally longer than permitted
        # in the database. Truncating avoids a
        # MySQL error.
        # TODO perhaps the database schema should be
        # changed to account for this.
        name = basic.truncate_value(name, 25, "...")
        description = ",".join(des_list[2:]).strip()

    return description, domain_id, name


def learn_cdd_name(cdd_dir):
    cdd_files = os.listdir(cdd_dir)
    cdd_files = [os.path.join(cdd_dir, x.split(".")[0]) for x in cdd_files]
    file_set = set(cdd_files)
    if len(file_set) == 1:
        cdd_name = cdd_files[0]
    else:
        cdd_name = ""
    return cdd_name


def main(argument_list):
    """
    :param argument_list:
    :return:
    """
    # Setup argument parser
    cdd_parser = setup_argparser()

    # Use argument parser to parse argument_list
    args = cdd_parser.parse_args(argument_list[2:])

    # Store arguments in more easily accessible variables
    database = args.database
    cdd_dir = expand_path(args.cdd)
    cdd_name = learn_cdd_name(cdd_dir)
    threads = args.threads
    evalue = args.evalue
    rpsblast = args.rpsblast
    tmp_dir = args.tmp_dir
    output_folder = args.output_folder
    reset = args.reset
    batch_size = args.batch_size

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]

    # Set up directory.
    output_folder = basic.set_path(output_folder, kind="dir", expect=True)
    results_folder = pathlib.Path(RESULTS_FOLDER)
    results_path = basic.make_new_dir(output_folder, results_folder,
                                      attempt=50)
    if results_path is None:
        print("Unable to create output_folder.")
        sys.exit(1)

    log_file = pathlib.Path(results_path, MAIN_LOG_FILE)

    # Set up root logger.
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG,
                        format="pdm_utils find_domains: %(levelname)s: %(message)s")
    logger.info(f"pdm_utils version: {VERSION}")
    logger.info(f"CDD run date: {constants.CURRENT_DATE}")
    logger.info(f"Command line arguments: {' '.join(argument_list)}")
    logger.info(f"Results directory: {results_path}")

    # Early exit if either 1) cdd_name == "" or 2) no rpsblast given and we are
    # unable to find one
    if cdd_name == "":
        msg = (f"Unable to learn CDD database name. Make sure the files in "
              f"{cdd_dir} all have the same basename.")
        logger.error(msg)
        print(msg)
        return

    # Get the rpsblast command and path.
    if rpsblast == "":
        command = get_rpsblast_command()
        rpsblast = get_rpsblast_path(command)

    # Verify database connection and schema compatibility.
    alchemist = AlchemyHandler(database=database,
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    logger.info(f"Connected to database: {database}.")
    mysqldb.check_schema_compatibility(engine, "the find_domains pipeline")
    logger.info(f"Schema version is compatible.")
    logger.info("Command line arguments verified.")

    if reset:
        logger.info("Clearing all domain data currently in the database.")
        clear_domain_data(engine)

    # Get gene data that needs to be processed
    # in dict format where key = column name, value = stored value.
    cdd_genes = mysqldb_basic.query_dict_list(engine, GET_GENES_FOR_CDD)
    msg = f"{len(cdd_genes)} genes to search for conserved domains..."
    logger.info(msg)
    print(msg)

    # Only run the pipeline if there are genes returned that need it
    if len(cdd_genes) > 0:

        log_gene_ids(cdd_genes)
        make_tempdir(tmp_dir)

        # Identify unique translations to process mapped to GeneIDs.
        cds_trans_dict = create_cds_translation_dict(cdd_genes)

        unique_trans = list(cds_trans_dict.keys())
        msg = (f"{len(unique_trans)} unique translations "
               "to search for conserved domains...")
        logger.info(msg)
        print(msg)

        # Process translations in batches. Otherwise, searching could take
        # so long that MySQL connection closes resulting in 1 or more
        # transaction errors.
        batch_indices = basic.create_indices(unique_trans, batch_size)
        total_rolled_back = 0
        for indices in batch_indices:
            start = indices[0]
            stop = indices[1]
            msg = f"Processing translations {start + 1} to {stop}..."
            logger.info(msg)
            print(msg)
            sublist = unique_trans[start:stop]
            batch_rolled_back = search_translations(
                                    rpsblast, cdd_name, tmp_dir, evalue,
                                    threads, engine, sublist, cds_trans_dict)
            total_rolled_back += batch_rolled_back

        search_summary(total_rolled_back)
        engine.dispose()

    return


def search_summary(rolled_back):
    """Print search results."""
    if rolled_back > 0:
        msg = (f"Error executing {rolled_back} transaction(s). "
               "Some genes may still contain unidentified domains.")
        logger.error(msg)
    else:
        msg = "All genes successfully searched for conserved domains."
        logger.info(msg)
    print("\n\n\n" + msg)


def search_translations(rpsblast, cdd_name, tmp_dir, evalue, threads,
                        engine, unique_trans, cds_trans_dict):
    """Search for conserved domains in a list of unique translations.
    """
    # Build jobs list
    jobs = []
    for id, translation in enumerate(unique_trans):
        jobs.append((rpsblast, cdd_name, tmp_dir, evalue, id, translation))
    results_temp = parallelize(jobs, threads, search_and_process)
    # Get rid of unneeded list type which was only needed due to
    # filter restriction in parallelize.start_processes()
    results = [i[0] for i in results_temp]

    # List of dictionaries. Each dictionary:
    # keys: "Translation": translation, "Data": list of results
    # In each list of results, each element is a dictionary:
    # data_dict = {
    #     "HitID": align.hit_id,
    #     "DomainID": domain_id,
    #     "Name": name,
    #     "Description": description,
    #     "Expect": float(hsp.expect),
    #     "QueryStart": int(hsp.query_start),
    #     "QueryEnd": int(hsp.query_end)
    #     }

    results_dict = create_results_dict(results)
    # Returns a dictionary, where:
    # key = unique translation,
    # value = list of dictionaries, each dictionary a unique rpsblast result

    transactions = construct_sql_txns(cds_trans_dict, results_dict)
    rolled_back = insert_domain_data(engine, transactions)
    return rolled_back


def create_cds_translation_dict(cdd_genes):
    """Create a dictionary of genes and translations.

    Returns a dictionary, where:
    key = unique translation
    value = set of GeneIDs with that translation."""
    trans_dict = {}
    for cds in cdd_genes:
        trans = cds["Translation"]
        gene_id = cds["GeneID"]
        if trans in trans_dict.keys():
            trans_dict[trans].add(gene_id)
        else:
            trans_dict[trans] = {gene_id}
    return trans_dict


def construct_sql_txns(cds_trans_dict, rpsblast_results):
    """Construct the list of SQL transactions."""
    # Since identical domain data may be used for > 1 genes, it is simply
    # duplicated for each gene. This could probably be simplified so that it
    # only tries to insert domain data once, and if that fails, then
    # all gene_domain and gene statements that depend on that domain data
    # will also fail.
    transactions = []
    for translation in rpsblast_results.keys():
        # Each translation has a list of rps_data
        rps_data_list = rpsblast_results[translation]
        # Each translation has a set of GeneIDs
        gene_ids = cds_trans_dict[translation]
        for gene_id in gene_ids:
            txn = construct_sql_txn(gene_id, rps_data_list)
            transactions.append(txn)
    return transactions


def construct_sql_txn(gene_id, rps_data_list):
    """Map domain data back to gene_id and create SQL statements for one transaction.

    rps_data_list is a list of dictionaries, where each dictionary reflects
    a significat rpsblast domain hit.
    """
    txn = []
    for rps_hit in rps_data_list:
        domain_stmt = construct_domain_stmt(rps_hit)
        gene_domain_stmt = construct_gene_domain_stmt(rps_hit, gene_id)
        txn.append(domain_stmt)
        txn.append(gene_domain_stmt)
    update_stmt = construct_gene_update_stmt(gene_id)
    txn.append(update_stmt)
    return txn


def construct_domain_stmt(data_dict):
    """Construct the SQL statement to insert data into the domain table."""
    stmt = INSERT_INTO_DOMAIN.format(data_dict["HitID"],
                                    data_dict["DomainID"],
                                    data_dict["Name"],
                                    data_dict["Description"])
    return stmt


def construct_gene_domain_stmt(data_dict, gene_id):
    """Construct the SQL statement to insert data into the gene_domain table."""
    stmt = INSERT_INTO_GENE_DOMAIN.format(gene_id,
                                          data_dict["HitID"],
                                          data_dict["Expect"],
                                          data_dict["QueryStart"],
                                          data_dict["QueryEnd"])
    return stmt


def construct_gene_update_stmt(gene_id):
    """Construct the SQL statement to update data in the gene table."""
    stmt = UPDATE_GENE.format(gene_id)
    return stmt


def create_results_dict(search_results):
    """Create a dictionary of search results

    Input is a list of dictionaries, one dict per translation, where:
    keys = "Translation" and "Data", where
    key = "Translation" has value = translation,
    key = "Data"" has value = list of rpsblast results, where
    Each result element is a dictionary containing domain and gene_domain data.

    Returns a dictionary, where:
    key = unique translation,
    value = list of dictionaries, each dictionary a unique rpsblast result
    """
    dict = {}
    for result in search_results:
        dict[result["Translation"]] = result["Data"]
    return dict


def get_rpsblast_command():
    """Determine rpsblast+ command based on operating system."""
    # See if we're running on a Mac
    if platform.system() == "Darwin":
        msg = "Detected MacOS operating system..."
        logger.info(msg)
        print(msg)
        command = shlex.split("which rpsblast")
    # Otherwise see if we're on a Linux machine
    elif platform.system() == "Linux":
        msg = "Detected Linux operating system..."
        logger.info(msg)
        print(msg)
        command = shlex.split("which rpsblast+")
    # Windows or others - unsupported, leave early
    else:
        msg = (f"Unsupported system '{platform.system()}'; cannot run "
              f"find_domains pipeline.")
        logger.error(msg)
        print(msg)
        sys.exit(1)
    return command


def get_rpsblast_path(command):
    """Determine rpsblast+ binary path."""

    # If we didn't exit, we have a command.
    # Run it, and PIPE stdout into rpsblast_path
    with Popen(args=command, stdout=PIPE) as proc:
        rpsblast = proc.stdout.read().decode("utf-8").rstrip("\n")

    # If empty string, rpsblast not found in globally
    # available executables, otherwise proceed with value.
    if rpsblast == "":
        msg = ("No rpsblast binary found. "
              "If you have rpsblast on your machine, please try "
              "again with the '--rpsblast' flag and provide the "
              "full path to your rpsblast binary.")
        logger.error(msg)
        print(msg)
        sys.exit(1)
    else:
        msg = f"Found rpsblast binary at '{rpsblast}'..."
        logger.info(msg)
        print(msg)
        return rpsblast


def log_gene_ids(cdd_genes):
    """Record names of the genes processed for reference."""
    batch_indices = basic.create_indices(cdd_genes, 20)
    for indices in batch_indices:
        genes_subset = cdd_genes[indices[0]: indices[1]]
        gene_ids = []
        for gene in genes_subset:
            gene_ids.append(gene["GeneID"])
        logger.info("; ".join(gene_ids))


def insert_domain_data(engine, results):
    """Attempt to insert domain data into the database."""
    msg = "Inserting data..."
    logger.info(msg)
    print(msg)

    rolled_back = 0
    connection = engine.connect()
    for result in results:
        exe_result = execute_transaction(connection, result)
        rolled_back += exe_result

    if rolled_back > 0:
        msg = (f"Error executing {rolled_back} transaction(s).")
        logger.error(msg)
        print(msg)
    return rolled_back


def execute_transaction(connection, statement_list=[]):
    trans = connection.begin()
    failed = 0
    index = 0
    # Even though the execution of individual statements are handled within
    # a try/exept block, this try/except block encapsulates all code
    # that is executed while uncommitted changes have been made to the database.
    # Without this encapsulation, even if all MySQL execution errors are caught,
    # other types of errors generated while processing the results would
    # crash the code and possibly result in changes made to the database
    # that shouldn't be persisted.
    try:
        # Try to execute statements as long as none of them failed.
        # Once one statement fails, don't try to execute any other statements.
        while (failed == 0 and index < len(statement_list)):
            statement = statement_list[index]
            result_tup = execute_statement(connection, statement)
            stmt_result = result_tup[0]
            type_error = result_tup[1]
            value_error = result_tup[2]
            msg = result_tup[3]

            msg = msg + "Statement: " + statement
            if stmt_result == 0:
                logger.info(msg)
            else:
                if (type_error == False and value_error == False):
                    logger.error(msg)
                else:
                    # If the insertion failed due to a TypeError, it could be
                    # due to the fact that the string contains '%'.
                    # Some CDD descriptions contain '%', which throws an error
                    # when SQLAlchemy's engine.Connection.execute() is used.
                    # The string gets passed to a pymysql.cursor function which
                    # interprets the % as a string formatting operator,
                    # and since there is no value to insert and format,
                    # the MySQL statement fails and the
                    # entire transactions is rolled back.
                    # For these edge cases, one way around this is to
                    # attempt to replace all '%' with '%%'.
                    # SQLAlchemy provides several different ways to
                    # implement changes to the database, and another strategy
                    # is likely to get around this problem.
                    if "%" not in statement:
                        logger.error(msg)
                    else:
                        logger.warning(msg)
                        logger.info("Attempting to resolve '%' error(s).")
                        statement = statement.replace("%", "%%")
                        result_tup = execute_statement(connection, statement)
                        stmt_result = result_tup[0]
                        type_error = result_tup[1]
                        value_error = result_tup[2]
                        msg = result_tup[3]

                        if stmt_result == 0:
                            logger.info(("The '%' replacement resolved the "
                                         "error(s)."))
                            logger.info(msg + statement)
                        else:
                            logger.info(("The '%' replacement failed to "
                                         "resolve the error(s)."))
                            logger.error(msg + statement)
            failed += stmt_result
            index += 1

        msg = (f"There are {failed} statements that failed execution.")
        if failed == 0:
            logger.info(msg)
            logger.info("Committing all changes.")
            trans.commit()
            txn_result = 0
        else:
            logger.error(msg)
            logger.info("Rolling back transaction.")
            trans.rollback()
            txn_result = 1
    except:
        print("Error executing MySQL statements.")
        print("Rolling back transaction...")
        logger.error("Unable to execute MySQL statements.")
        logger.info("Rolling back transaction.")
        trans.rollback()
        txn_result = 1

    return txn_result


def execute_statement(connection, statement):
    type_error = False
    value_error = False
    try:
        connection.execute(statement)
    except sqlalchemy.exc.DBAPIError as err:
        err_stmt = err.statement
        sqla_err_type = str(type(err))
        pymysql_err_type = str(type(err.orig))
        pymysql_err_code = err.orig.args[0]
        pymysql_err_msg = err.orig.args[1]
        msg = (f"SQLAlchemy Error type: {sqla_err_type}. "
               f"PyMySQL Error type: {pymysql_err_type}. "
               f"PyMYSQL Error code: {pymysql_err_code}. "
               f"PyMySQL Error message: {pymysql_err_msg}.")
        if pymysql_err_code == 1062:
            msg = "Duplicate entry error ignored. " + msg
            result = 0
        else:
            msg = "Unable to execute MySQL statement. " + msg
            result = 1
    except TypeError as err:
        type_error = True
        msg = "Unable to execute statement due to a TypeError. "
        result = 1
    except ValueError as err:
        value_error = True
        msg = "Unable to execute statement due to a ValueError. "
        result = 1
    except:
        # TODO not sure how to test this block. Would need to construct a
        # statement that causes a Python built-in exception other
        # than TypeError or ValueError.
        msg = "Unable to execute statement. "
        result = 1
    else:
        msg = "Successful statement execution. "
        result = 0

    return result, type_error, value_error, msg


def clear_domain_data(engine):
    """Delete all domain data stored in the database."""
    connection = engine.connect()
    stmts = [CLEAR_GENE_DOMAIN, CLEAR_DOMAIN, CLEAR_GENE_DOMAINSTATUS]
    exe_result = execute_transaction(connection, stmts)
    if exe_result == 1:
        logger.error("Unable to clear all domain data.")
    else:
        logger.info("All domain data cleared.")
