"""Pipeline that uses rpsblast to search for NCBI conserved domain database
hits in all new proteins in the phage database."""

import argparse
import logging
import pathlib
import shutil
import sys
import tempfile

import sqlalchemy

import pdm_utils
from pdm_utils.classes.alchemyhandler import AlchemyHandler, Engine
from pdm_utils.constants import constants
from pdm_utils.functions import basic, configfile, mysqldb, mysqldb_basic
from pdm_utils.functions import multiprocess as mp, rpsblast


# Miscellaneous
VERSION = pdm_utils.__version__
EVALUE = 0.001
BATCH = 10000
OUTDIR = pathlib.Path().cwd()
RESULTS_FOLDER = f"{constants.CURRENT_DATE}_find_domains"
MAIN_LOG_FILE = "find_domains.log"
DEFAULT_CDD = pathlib.Path().home().joinpath("Databases/Cdd_LE")

# Add a logger named after this module, with a null handler that will suppress
# any output statements. This allows other modules to define the handler and
# output formats. If this module is the main module being called, the top
# level main function configures the root logger and primary file handle.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# MySQL queries and commands for everything else
GET_NEW_GENES = "SELECT GeneID, Translation FROM gene WHERE DomainStatus = 0"
UPDATE_STATUS = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"
INSERT_DOMAIN = """INSERT INTO domain (HitID, DomainID, Name, Description) """ \
                """VALUES ("{}", "{}", "{}", "{}")"""
INSERT_GENE_DOMAIN = """INSERT INTO gene_domain """ \
                     """(GeneID, HitID, Expect, QueryStart, QueryEnd) """ \
                     """VALUES ("{}", "{}", "{}", "{}", "{}")"""


def parse_args(args=None):
    """Parse commandline arguments."""
    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("database", type=str,
                        help="name of database to find conserved domains for")

    parser.add_argument("--cdd-dir", type=str, default=DEFAULT_CDD,
                        help=f"path to local NCBI Cdd [default: {DEFAULT_CDD}]")
    parser.add_argument("--cpus", default=mp.CPUS, type=int,
                        help=f"number of CPU cores to use [default: {mp.CPUS}]")
    parser.add_argument("-e", "--evalue", default=EVALUE, type=float,
                        help=f"significance threshold to retain hits "
                             f"[default: {EVALUE}]")
    parser.add_argument("-o", "--output-dir", type=pathlib.Path, default=OUTDIR,
                        help=f"directory where log data can be output "
                             f"[default: {OUTDIR}]")
    parser.add_argument("-b", "--batch-size", default=BATCH, type=int,
                        help=f"number of translations to search at a time "
                             f"[default: {BATCH}]")
    parser.add_argument("-r", "--rpsblast", default=None, type=str,
                        help="path to rpsblast(+) binary")
    parser.add_argument("-x", "--reset", action="store_true",
                        help="reset NCBI Cdd hits for all genes")
    parser.add_argument("-c", "--config", type=pathlib.Path,
                        help="path to a pdm_utils configuration file "
                             "containing MySQL login data")

    if args:
        return parser.parse_args(args)

    return parser.parse_args()


def detect_rpsblast_name():
    """Detect whether rpsblast or rpsblast+ is installed locally.

    :return: name
    """
    for program in ["rpsblast", "rpsblast+"]:
        version = rpsblast.get_version(program)
        if version:
            return program


def learn_cdd_name(cdd_dir):
    """Return the NCBI Cdd 'name' learned by analyzing the filename of
    every file in the Cdd directory.

    :param cdd_dir: directory where Cdd database files live
    :type cdd_dir: pathlib.Path
    :return: cdd_name
    """
    filenames = list()
    for filepath in cdd_dir.iterdir():
        if filepath.is_file() and filepath.name != ".DS_Store":
            filenames.append(filepath.stem)

    name_set = set(filenames)
    if len(name_set) == 1:
        return filenames[0]


def find_domains(program, sequences, cdd, evalue, cpus):
    """Find and return NCBI conserved domain database hits for the
    given sequences.

    :param program: name (or path) to the program to run (e.g. 'rpsblast')
    :type program: str or pathlib.Path
    :param sequences: the sequences to search for Cdd hits in
    :type sequences: list[tuple[str, str]]
    :param cdd: path to a local copy of the Cdd database
    :type cdd: str or pathlib.Path
    :param evalue: significance threshold to retain hits
    :type evalue: float
    :param cpus: number of cpu cores to use
    :type cpus: int
    :return: domains
    """
    # Calling code needs to clean up temp dir created this way - see
    # try/finally block below
    tmp_dir = tempfile.mkdtemp(prefix="pdm_utils__find_domains")

    try:
        outfiles = rpsblast.blastall(program=program, sequences=sequences,
                                     db=cdd, tmp_dir=pathlib.Path(tmp_dir),
                                     evalue=evalue, cpus=cpus)
        domains = dict()
        for outfile in outfiles:
            temp_domains = rpsblast.parse_xml(outfile, sequences, evalue)
            for key, value in temp_domains.items():
                domains[key] = value
    finally:
        shutil.rmtree(tmp_dir)

    return domains


def main(argument_list):
    """Commandline interface for this pipeline.

    :param argument_list: the list of CLI args given for this pipeline
    :type argument_list: list[str]
    """
    # Parse commandline arguments
    args = parse_args(argument_list[2:])

    # Store arguments in more easily accessible variables
    database = args.database
    cdd_dir = args.cdd_dir
    cdd_name = learn_cdd_name(cdd_dir)
    cpus = args.cpus
    evalue = args.evalue
    program = args.rpsblast
    output_dir = args.output_dir
    reset = args.reset
    batch_size = args.batch_size

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config)
    mysql_creds = config["mysql"]

    # Set up directory.
    output_folder = basic.set_path(output_dir, kind="dir", expect=True)
    results_folder = pathlib.Path(RESULTS_FOLDER)
    results_path = basic.make_new_dir(output_folder, results_folder,
                                      attempt=50)
    if results_path is None:
        print("Unable to create output_folder.")
        sys.exit(1)

    log_file = pathlib.Path(results_path, MAIN_LOG_FILE)

    # Set up root logger.
    logging.basicConfig(
        filename=log_file, filemode="w", level=logging.DEBUG,
        format="pdm_utils find_domains: %(levelname)s: %(message)s")
    logger.info(f"pdm_utils version: {VERSION}")
    logger.info(f"CDD run date: {constants.CURRENT_DATE}")
    logger.info(f"Command line arguments: {' '.join(argument_list)}")
    logger.info(f"Results directory: {results_path}")

    # Early exit if not cdd_name or no rpsblast given and we can't find one
    if not cdd_name:
        msg = (f"Unable to learn CDD database name. Make sure the files in "
               f"{cdd_dir} all have the same basename.")
        logger.error(msg)
        print(msg)
        return

    cdd = cdd_dir.joinpath(cdd_name)

    # If no rpsblast given, try to discover a globally executable version
    if not program:
        program = detect_rpsblast_name()

    if not program:
        msg = "Unable to auto-detect an rpsblast or rpsblast+ executable on " \
              "$PATH. Please install the NCBI-BLAST+ toolkit. If already " \
              "installed, re-run this pipeline with the '--rpsblast' " \
              "argument and provide the path to the rpsblast(+) executable."
        logger.error(msg)
        print(msg)
        return

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
        msg = "Resetting all domain data currently in the database."
        logger.info(msg)
        print(msg)
        reset_domain_data(engine)

    # Get gene data that needs to be processed
    # in dict format where key = column name, value = stored value.
    cdd_genes = mysqldb_basic.query_dict_list(engine, GET_NEW_GENES)
    msg = f"{len(cdd_genes)} genes to search for conserved domains..."
    logger.info(msg)
    print(msg)

    # Only run the pipeline if there are genes returned that need it
    if not cdd_genes:
        engine.dispose()
        return

    # Log the GeneIDs we're going to search
    log_gene_ids(cdd_genes)

    # Create dictionary keyed on translation
    nr_dict = dict()
    for gene_dict in cdd_genes:
        geneid = gene_dict["GeneID"]
        translation = gene_dict["Translation"].decode("utf-8")
        if translation in nr_dict:
            nr_dict[translation].append(geneid)
        else:
            nr_dict[translation] = [geneid]
    unique_trans = list(nr_dict.keys())

    msg = (f"{len(unique_trans)} unique translations to search for "
           f"conserved domains...")
    logger.info(msg)
    print(msg)

    # Process translations in batches. Otherwise, searching could take
    # so long that MySQL connection closes resulting in 1 or more
    # transaction errors.
    batch_indices = basic.create_indices(unique_trans, batch_size)
    total_rolled_back = 0
    for start, stop in batch_indices:
        msg = f"Processing translations {start + 1} to {stop}..."
        logger.info(msg)
        print(msg)

        batch_sequences = list()
        for translation in unique_trans[start:stop]:
            batch_sequences.append((nr_dict[translation][0], translation))

        batch_domains = find_domains(program, batch_sequences, cdd,
                                     evalue, cpus)

        batch_transactions = list()
        for translation, domains in batch_domains.items():
            for geneid in nr_dict[translation]:
                batch_transactions.append(create_transaction(geneid, domains))

        batch_rolled_back = insert_domain_data(engine, batch_transactions)
        total_rolled_back += batch_rolled_back

    search_summary(total_rolled_back)
    engine.dispose()


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


def create_transaction(geneid, domains):
    """Create a series of MySQL INSERT/UPDATE statements to associate
    a geneid with its newly found conserved domain data.

    :param geneid: the gene to add conserved domain hits for
    :type geneid: str
    :param domains: the conserved domain hit data
    :type domains: list[dict]
    :return: transaction
    """
    transaction = list()

    for hit in domains:
        transaction.append(insert_domain_statement(hit))
        transaction.append(insert_gene_domain_statement(geneid, hit))

    transaction.append(update_gene_status_statement(geneid))

    return transaction


def log_gene_ids(cdd_genes):
    """Record names of the genes processed for reference."""
    batch_indices = basic.create_indices(cdd_genes, 20)
    for indices in batch_indices:
        genes_subset = cdd_genes[indices[0]: indices[1]]
        gene_ids = []
        for gene in genes_subset:
            gene_ids.append(gene["GeneID"])
        logger.info("; ".join(gene_ids))


def insert_domain_statement(domain):
    """Format a SQL statement to insert into the domain table.

    :param domain: a conserved domain hit
    :type domain: dict
    :return: statement
    """
    hit_id, domain_id = domain["HitID"], domain["DomainID"]
    name, description = domain["Name"], domain["Description"]

    if len(name) > 25:
        name = basic.truncate_value(name, 25, "...")

    return INSERT_DOMAIN.format(hit_id, domain_id, name, description)


def insert_gene_domain_statement(geneid, domain):
    """Format a SQL statement to insert into the gene_domain table.

    :param geneid: the gene in which a domain hit was found
    :type geneid: str
    :param domain: a conserved domain hit
    :type domain: dict
    :return: statement
    """
    hit_id, expect = domain["HitID"], domain["Expect"]
    start, end = domain["QueryStart"], domain["QueryEnd"]

    return INSERT_GENE_DOMAIN.format(geneid, hit_id, expect, start, end)


def update_gene_status_statement(geneid):
    """Format a SQL statement to update a gene's domain-search status
    in the gene table.

    :param geneid: a gene for which domain hits were searched
    :type geneid: str
    :return: statement
    """
    """Format the SQL statement to update data in the gene table."""
    return UPDATE_STATUS.format(geneid)


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
        msg = f"Error executing {rolled_back} transaction(s)."
        logger.error(msg)
        print(msg)

    return rolled_back


def reset_domain_data(engine):
    """Reset all domain data stored in the database.

    :param engine: an initialized SQLAlchemy engine
    :type engine: Engine
    :return: failed
    """
    connection = engine.connect()

    transaction = ["TRUNCATE gene_domain", "DELETE FROM domain",
                   "UPDATE gene SET DomainStatus = 0"]

    failed = execute_transaction(connection, transaction)
    if failed:
        logger.error("Unable to clear all domain data.")
    else:
        logger.info("All domain data cleared.")

    return failed


def execute_transaction(connection, statements):
    """Attempt to execute a series of MySQL statements. Failure to
    execute any statement results in rolling back the full transaction.

    :param connection:
    :type connection:
    :param statements:
    :type statements:
    :return: status
    """
    transaction = connection.begin()

    failed = False
    for statement in statements:
        try:
            connection.execute(statement)
        except sqlalchemy.exc.DBAPIError as e:
            pymysql_error_code = e.orig.args[0]
            pymysql_error_message = e.orig.args[1]
            msg = (f"SQLAlchemy Error type: {str(type(e))}. "
                   f"PyMySQL Error type: {str(type(e.orig))}. "
                   f"PyMYSQL Error code: {pymysql_error_code}. "
                   f"PyMySQL Error message: {pymysql_error_message}.")
            if pymysql_error_code == 1062:
                logger.info("Duplicate entry error ignored. " + msg + statement)
            else:
                logger.error("Failed statement execution. " + msg + statement)
                failed = True
                break
        except Exception as e:
            # One cause is that Cdd description contains '%'. Try replacing
            # '%' with '%%' and re-running the statement.
            if "%" in statement:
                logger.warning(f"Failed statement execution ({type(e)}).")
                logger.info("Attempting to resolve '%' error(s).")
                statement = statement.replace("%", "%%")
                try:
                    connection.execute(statement)
                    logger.info("Replacing '%' with '%%' resolved the error.")
                    logger.info("Successful statement execution " + statement)
                except sqlalchemy.exc.DBAPIError as e:
                    pymysql_error_code = e.orig.args[0]
                    if pymysql_error_code == 1062:
                        logger.info("Duplicate entry error ignored. " +
                                    statement)
                except Exception as e:
                    logger.info("Replacing '%' with '%%' did not resolve the "
                                "error.")
                    logger.error(f"Failed statement execution ({type(e)}). " +
                                 statement)
                    failed = True
                    break
            else:
                logger.error(f"Failed statement execution ({type(e)}). " +
                             statement)
                failed = True
                break

    if failed:
        logger.info("Rolling back transaction.")
        transaction.rollback()
    else:
        logger.info("Committing all changes.")
        transaction.commit()

    transaction.close()

    return failed
