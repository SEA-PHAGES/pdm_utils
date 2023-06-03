"""Pipeline that uses DeepTMHMM to identify signal peptide and transmembrane
regions of all new proteins in the phage database."""

# TODO: implement some kind of retry task for remote scanning, defined by
#  --max-retry (number of attempts before giving up) and
#  --retry-interval (time to wait between attempts).
#  The retry interval would be given in minutes, but should probably default
#  to something pretty large, like 360 (6 hours), to account for the fact
#  that the DTU rate limiting seems pretty inconsistent.

# TODO: implement local scanning that depends on availability of a GPU for
#  handling more than for example 100 genes. CPU scanning is VERY slow and
#  likely only works on Linux servers.

import argparse
import pathlib
import shutil
import sys
import time

import sqlalchemy
from sqlalchemy.exc import DBAPIError
from biolib.biolib_errors import BioLibError

import pdm_utils
from pdm_utils.constants import constants
from pdm_utils.functions import (basic, configfile, mysqldb, pipelines_basic)
from pdm_utils.functions.fileio import write_fasta
from pdm_utils.functions.deeptmhmm import run_deeptmhmm

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_membrane"

# Miscellaneous
VERSION = pdm_utils.__version__
BATCH = 100
RESULTS_FOLDER = f"{constants.CURRENT_DATE}_deeptmhmm"

# MySQL queries and commands for everything else
GET_NEW_GENES = "SELECT GeneID, Translation FROM gene WHERE MembraneStatus = 0"
UPDATE_STATUS = "UPDATE gene SET MembraneStatus = 1 WHERE GeneID = '{}'"
INSERT_DOMAIN = """INSERT INTO gene_transmembrane """\
                """(GeneID, QueryStart, QueryEnd, Type, Source) """ \
                """VALUES ("{}", "{}", "{}", "{}", "{}")"""


def parse_args(args=None):
    """Parse commandline arguments."""
    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("database", type=str,
                        help="name of database for which to find secretion "
                             "signals and transmembrane domains ")

    parser.add_argument("-c", "--config", type=pathlib.Path,
                        help="path to a pdm_utils configuration file "
                             "containing MySQL login data")
    parser.add_argument("-m", "--folder_name", type=str,
                        help="set the name of the folder to be exported")
    parser.add_argument("-o", "--folder_path", type=pathlib.Path,
                        help="set the path of the directory where the exported"
                             "files are stored.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="enables progress print statements")

    parser.add_argument("-if", "--import_file",
                        type=pipelines_basic.convert_file_path,
                        help="selection import option that imports values "
                             "from the first column of a csv file.",
                        dest="input", default=[])
    parser.add_argument("-in", "--import_names", nargs="*",
                        help="selection input option that imports values "
                             "from command line input.", dest="input")
    parser.add_argument("-w", "--where", nargs="?", default="",
                        help="Data filtering option that filters data "
                             "by the inputted expressions.", dest="filters")

    parser.add_argument("-b", "--batch_size", default=BATCH, type=int,
                        help=f"number of translations to search at a time "
                             f"[default: {BATCH}]")
    parser.add_argument("-Mb", "--maximum-batches", default=None,
                        type=int,
                        help="number of maximum batches to run in total "
                             "[default: infinite]")
    parser.add_argument("-rm", "--run_machine", default="local",
                        choices=["local", "remote"], type=str,
                        help="runmode option to select whether deeptmhmm "
                             "will be run locally or by using remote servers.")
    parser.add_argument("-x", "--reset", action="store_true",
                        help="reset signal peptide and transmembrane domain "
                             "hits for all genes")

    if args:
        return parser.parse_args(args)

    return parser.parse_args()


def cleanup_and_exit(tmpdir, code):
    """Invoke this function for a clean exit (delete tmpdir and exit
    with the given code).

    :param tmpdir: the path where temporary files are located
    :type tmpdir: pathlib.Path
    :param code: the exit code that other programs could poll for
    :type code: int
    """
    if tmpdir is not None and tmpdir.is_dir():
        shutil.rmtree(tmpdir)

    sys.exit(code)


def main(argument_list):
    """Commandline interface for this pipeline.

    :param argument_list: the list of CLI args given for this pipeline
    :type argument_list: list[str]
    """
    # Parse commandline arguments
    args = parse_args(argument_list[2:])

    verbose = args.verbose

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)
    if not values:
        values = None

    mysqldb.check_schema_compatibility(alchemist.engine,
                                       "the deeptmhmm pipeline", 11)

    execute_deeptmhmm_scan(alchemist, folder_path=args.folder_path,
                           folder_name=DEFAULT_FOLDER_NAME, values=values,
                           verbose=verbose, filters=args.filters,
                           reset=args.reset, run_machine=args.run_machine,
                           batch_size=args.batch_size,
                           max_batches=args.maximum_batches)


def execute_deeptmhmm_scan(alchemist, folder_path=None,
                           folder_name=DEFAULT_FOLDER_NAME, values=None,
                           verbose=False, filters="", reset=False,
                           run_machine="local",
                           batch_size=BATCH, max_batches=None):
    export_path = pipelines_basic.create_working_path(folder_path,
                                                      folder_name)
    pipelines_basic.create_working_dir(export_path)

    if reset:
        if verbose:
            print("Resetting all membrane data currently in the database.")
        reset_domain_data(alchemist.engine)

    db_filter = pipelines_basic.build_filter(alchemist, "gene", filters,
                                             values=values,
                                             verbose=verbose)
    db_filter.update()
    if db_filter.hits() == 0:
        print("No database entries received from gene "
              "for '{export_path}'.")

    values = db_filter.values
    db_filter.reset()
    db_filter.values = values
    db_filter.add("gene.MembraneStatus = 0")
    db_filter.update()
    if db_filter.hits() == 0:
        print("All selected genes have pre-existing search results.")
        return

    cdd_genes = db_filter.select(["gene.GeneID", "gene.Translation"])
    # Get gene data that needs to be processed
    # in dict format where key = column name, value = stored value.
    if verbose:
        print(f"{len(cdd_genes)} genes to search for tranemembrane domains...")

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

    if verbose:
        print(f"{len(unique_trans)} unique translations to search for "
              f"transmembrane domains...")

    # Process translations in batches. Otherwise, searching could take
    # so long that MySQL connection closes resulting in 1 or more
    # transaction errors.
    batch_indices = basic.create_indices(unique_trans, batch_size)
    total_rolled_back = 0
    for i, (start, stop) in enumerate(batch_indices):
        if max_batches:
            if i >= max_batches:
                break

        if verbose:
            print(f"Processing translations {start+1} to {stop}...")

        gs_to_ts = dict()
        for translation in unique_trans[start:stop]:
            gs_to_ts[nr_dict[translation][0]] = translation

        batch_file = export_path.joinpath(f"batch_{i+1}.fasta")
        write_fasta(gs_to_ts, batch_file)

        if verbose:
            print("\tfinding transmembrane domains...")
        try:
            batch_domains = run_deeptmhmm(batch_file, machine=run_machine)
        except BioLibError as err:
            print(f"Unable to continue - {err.message}...")
            return
        except TypeError as err:
            print(f"Unable to continue - {err.message}...")
            return

        if verbose:
            print("\tupdating database...")
        batch_transactions = list()
        for translation, domains in batch_domains.items():
            for geneid in nr_dict[translation]:
                batch_transactions.append(create_transaction(geneid, domains))

        batch_rolled_back = insert_domain_data(alchemist.engine,
                                               batch_transactions)
        total_rolled_back += batch_rolled_back

    search_summary(total_rolled_back)
    alchemist.engine.dispose()


def search_summary(rolled_back):
    """Print search results."""
    if rolled_back > 0:
        msg = (f"Error executing {rolled_back} transaction(s). "
               "Some genes may still contain unidentified domains.")
    else:
        msg = "All genes successfully searched for transmembrane domains."
    print("\n\n\n" + msg)


def create_transaction(geneid, domains, source="deeptmhmm"):
    """Create a series of MySQL INSERT/UPDATE statements to associate
    a geneid with its newly found transmembrane domain data.

    :param geneid: the gene to add transmembrane domain hits for
    :type geneid: str
    :param domains: the transmembrane domain hit data
    :type domains: list[tuple[str, tuple[int, int]]]
    :param source: which tool did the data come from?
    :type source: str
    :return: transaction
    """
    if source not in ("deeptmhmm", "sosui"):
        raise ValueError(f"invalid data source '{source}': valid sources are "
                         f"'deeptmhmm' and 'sosui'")

    transaction = list()

    for domain_type, (query_start, query_end) in domains:
        statement = INSERT_DOMAIN.format(geneid, query_start, query_end,
                                         domain_type, source)
        transaction.append(statement)

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


def update_gene_status_statement(geneid):
    """Format a SQL statement to update a gene's domain-search status
    in the gene table.

    :param geneid: a gene for which domain hits were searched
    :type geneid: str
    :return: statement
    """
    """Format the SQL statement to update data in the gene table."""
    return UPDATE_STATUS.format(geneid)


def insert_domain_data(engine, results, verbose=False):
    """Attempt to insert domain data into the database."""
    if verbose:
        print("\tInserting data...")

    rolled_back = 0
    connection = engine.connect()
    for result in results:
        exe_result = execute_transaction(connection, result)
        rolled_back += exe_result

    if rolled_back > 0:
        print(f"Error executing {rolled_back} transaction(s).")

    return rolled_back


def reset_domain_data(engine):
    """Reset all domain data stored in the database.

    :param engine: an initialized SQLAlchemy engine
    :type engine: Engine
    :return: failed
    """
    connection = engine.connect()

    transaction = ["TRUNCATE gene_transmembrane",
                   "UPDATE gene SET MembraneStatus = 0"]

    failed = execute_transaction(connection, transaction)
    if failed:
        print("Unable to clear all domain data.")
    else:
        print("All domain data cleared.")

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
                pass
            else:
                print(msg)
                failed = True
                break
        except Exception:
            # One cause is that Cdd description contains '%'. Try replacing
            # '%' with '%%' and re-running the statement.
            if "%" in statement:
                statement = statement.replace("%", "%%")
                try:
                    connection.execute(statement)
                except DBAPIError as e:
                    pymysql_error_code = e.orig.args[0]
                    if pymysql_error_code == 1062:
                        pass
                    else:
                        failed = True
                        break
            else:
                failed = True
                break

    if failed:
        transaction.rollback()
    else:
        transaction.commit()

    transaction.close()

    return failed
