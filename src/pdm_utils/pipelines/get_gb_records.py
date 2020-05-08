"""Pipeline to retrieve GenBank records using accessions stored in the MySQL database."""

import argparse
from datetime import date
import os
import pathlib
import sys

from Bio import SeqIO

from pdm_utils.classes.filter import Filter
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import parsing
from pdm_utils.functions import querying


DEFAULT_OUTPUT_FOLDER = os.getcwd()
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_gb_records"
TARGET_TABLE = "phage"

# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting GenBank records."""
    GET_GB_RECORDS_HELP = ("Pipeline to retrieve GenBank records using "
                           "accessions stored in a MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the directory where records will be stored.")
    NCBI_CRED_FILE_HELP = ("Path to the file containing NCBI credentials.")
    FILTERS_HELP = (
        "Indicates which genomes to retain in the new database."
        "Follow selection argument with formatted filter request: "
        "Table.Field=Value"
        )

    parser = argparse.ArgumentParser(description=GET_GB_RECORDS_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
                        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                        help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-c", "--ncbi_credentials_file", type=pathlib.Path,
                        help=NCBI_CRED_FILE_HELP)
    parser.add_argument("-f", "--filters", nargs="?",
                        type=parsing.parse_cmd_string, help=FILTERS_HELP,
                        default=[])

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args



# TODO unittest.
def main(unparsed_args_list):
    """Run main get_gb_records pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)

    # Filters input: phage.Status=draft AND phage.HostGenus=Mycobacterium
    # Args structure: [['phage.Status=draft'], ['phage.HostGenus=Mycobacterium']]
    filters = args.filters
    ncbi_cred_dict = ncbi.get_ncbi_creds(args.ncbi_credentials_file)
    output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    working_dir = pathlib.Path(RESULTS_FOLDER)
    working_path = basic.make_new_dir(output_folder, working_dir,
                                      attempt=50)
    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    alchemist = AlchemyHandler(database=args.database)
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the get_gb_records pipeline")

    # Get SQLAlchemy metadata Table object
    # table_obj.primary_key.columns is a
    # SQLAlchemy ColumnCollection iterable object
    # Set primary key = 'phage.PhageID'
    alchemist.build_metadata()
    table = querying.get_table(alchemist.metadata, TARGET_TABLE)
    for column in table.primary_key.columns:
        primary_key = column

    # Create filter object and then add command line filter strings
    db_filter = Filter(alchemist=alchemist, key=primary_key)
    db_filter.values = []

    # Attempt to add filters and exit if needed.
    add_filters(db_filter, filters)

    # Performs the query
    db_filter.update()

    # db_filter.values now contains list of PhageIDs that pass the filters.
    # Get the accessions associated with these PhageIDs.
    keep_set = set(db_filter.values)


    # Create data sets
    print("Retrieving accessions from the database...")
    query = construct_accession_query(keep_set)
    list_of_dicts = mysqldb_basic.query_dict_list(engine, query)
    id_acc_dict = get_id_acc_dict(list_of_dicts)
    acc_id_dict = get_acc_id_dict(id_acc_dict)
    engine.dispose()
    if len(acc_id_dict.keys()) > 0:
        get_data(working_path, acc_id_dict, ncbi_cred_dict)
    else:
        print("There are no records to retrieve.")



# TODO test.
def get_id_acc_dict(list_of_dicts):
    """Convert list of dictionaries from MySQL query into a single dictionary.

    Receives a list of dictionaries, where Keys are PhageID and Accession.
    Returns a dictionary where Key = PhageID and Value = Accession
    """
    # [{'PhageID':'Trixie', 'Accession':'ABC123'}]
    new_dict = {}
    for dict in list_of_dicts:
        phage_id = dict["PhageID"]
        accession = dict["Accession"]
        new_dict[phage_id] = accession
    return new_dict


# TODO move to basic.
# TODO test.
def get_acc_id_dict(id_acc_dict):
    """Invert a dictionary.

    Receives a dictionary where Key = PhageID and Value = Accession.
    Returns a dictionary where Key = Accession and Value = List[PhageIDs]
    Accessions are not guaranteed to be unique.
    """
    new_dict = {}
    for key in id_acc_dict.keys():
        accession = id_acc_dict[key]
        if (accession is not None and accession != ""):
            if accession not in new_dict.keys():
                new_dict[accession] = [key]
            else:
                new_dict[accession] = new_dict[accession].append(accession)
    return new_dict


# TODO move to ncbi.
# TODO test.
def get_data(output_folder, acc_id_dict, ncbi_cred_dict={}, batch_size=200):
    """Retrieve genomes from GenBank.

    output_folder = Path to where files will be saved.
    acc_id_dict = Dictionary where key = Accession and value = List[PhageIDs]
    """

    # More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    ncbi.set_entrez_credentials(
        tool=ncbi_cred_dict["ncbi_tool"],
        email=ncbi_cred_dict["ncbi_email"],
        api_key=ncbi_cred_dict["ncbi_api_key"])


    # Use esearch to verify the accessions are valid and efetch to retrieve
    # the record
    # Create batches of accessions
    unique_accession_list = list(acc_id_dict.keys())

    # Add [ACCN] field to each accession number
    appended_accessions = \
        [accession + "[ACCN]" for accession in unique_accession_list]


    # When retrieving in batch sizes, first create the list of values
    # indicating which indices of the unique_accession_list should be used
    # to create each batch.
    # For instace, if there are five accessions, batch size of two produces
    # indices = 0,2,4
    batch_indices = basic.create_indices(unique_accession_list, batch_size)
    print(f"There are {len(unique_accession_list)} GenBank accessions to check.")
    for indices in batch_indices:
        start = indices[0]
        stop = indices[1]
        print(f"Checking accessions {start + 1} to {stop}...")

        delimiter = " | "
        esearch_term = delimiter.join(appended_accessions[start:stop])

        # Use esearch for each accession
        search_record = ncbi.run_esearch(db="nucleotide", term=esearch_term,
                            usehistory="y")
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]
        summary_records = ncbi.get_summaries(db="nucleotide",
                                             query_key=search_query_key,
                                             webenv=search_webenv)

        accessions_to_retrieve = ncbi.get_accessions_to_retrieve(summary_records)
        if len(accessions_to_retrieve) > 0:
            records = ncbi.get_records(accessions_to_retrieve,
                                       db="nucleotide",
                                       rettype="gb",
                                       retmode="text")
            for record in records:
                output_data(record, acc_id_dict, output_folder)


# TODO test.
def output_data(seqrecord, acc_id_dict, output_folder):
    """Save retrieved record to file."""
    accession = seqrecord.annotations["accessions"][0]
    phage_id_list = acc_id_dict[accession]
    for phage_id in phage_id_list:
        filename = (f"{phage_id}__{accession}.gb")
        filepath = pathlib.Path(output_folder, filename)
        SeqIO.write(seqrecord, str(filepath), "genbank")

# TODO move to basic or mysqldb.
# TODO test.
def construct_set_string(phage_id_set):
    """Convert set of phage_ids to string formatted for MySQL.

    e.g. set: {'Trixie', 'L5', 'D29'}
    returns: "('Trixie', 'L5', 'D29')""
    """
    string = "('" + "', '".join(phage_id_set) + "')"
    return string

# TODO move to mysqldb.
# TODO test.
def construct_accession_query(phage_id_set):
    """Construct SQL query to retrieve accessions."""
    phage_id_string = construct_set_string(phage_id_set)
    query = (f"SELECT PhageID, Accession FROM phage "
             f"WHERE PhageID IN {phage_id_string}")
    return query

# TODO test.
def add_filters(filter_obj, filters):
    """Add filters from command line to filter object."""
    errors = 0
    for or_filters in filters:
        for filter in or_filters:
            # Catch the error if it is an invalid table.column
            try:
                filter_obj.add(filter)
            except:
                print(f"Invalid filter: filter")
                errors += 1
    if errors > 0:
        print("Unable to create new database.")
        sys.exit(1)
