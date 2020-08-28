"""Pipeline to gather new data to be imported into a MySQL database."""

import argparse
import csv
from datetime import datetime, date
import os
import pathlib
import sys

from Bio import SeqIO

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes import genomepair
from pdm_utils.classes import ticket
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio
from pdm_utils.functions import ncbi
from pdm_utils.functions import mysqldb
from pdm_utils.functions import phagesdb
from pdm_utils.functions import tickets

# Names of folders and files created.
DEFAULT_OUTPUT_FOLDER = os.getcwd()
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_data"
PECAAN_FOLDER = "pecaan"
GENBANK_FOLDER = "genbank"
PHAGESDB_FOLDER = "phagesdb"
UPDATES_FOLDER = "updates"
UPDATE_TABLE = "update_table.csv"
IMPORT_TABLE = "import_table.csv"
GENBANK_RESULTS_TABLE = "genbank_results.csv"
GENOME_FOLDER = "genomes"

# Columns for import ticket table.
IMPORT_COLUMNS = constants.IMPORT_TABLE_STRUCTURE["order"]

# Columns for update ticket table.
UPDATE_COLUMNS = ["table",
                  "field",
                  "value",
                  "key_name",
                  "key_value"]

# Column headers for GenBank results file.
NCBI_RESULTS_COLUMNS = ["phage_id",
                       "phage_name",
                       "accession",
                       "annotation_status",
                       "mysql_date",
                       "genbank_date",
                       "result"]

# GenBank retrieval status options
NOT_AUTO = "not_auto_updated"
NO_ACC = "no_accession"
DUPE_ACC = "duplicate_accession"
OLD_DOCSUM = "docsum_not_new"
OLD_RECORD = "record_not_new"
RET_FAIL = "retrieval_failure"
RETRIEVED = "retrieved_for_import"
GENBANK_STATUS = {NOT_AUTO, NO_ACC, DUPE_ACC, OLD_DOCSUM,
                  OLD_RECORD, RET_FAIL, RETRIEVED}

# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""

    retrieve_help = "Pipeline to get new data to import into a database."
    database_help = "Name of the MySQL database."
    output_folder_help = "Path to the directory where updates will be stored."
    updates_help = "Retrieve field updates from PhagesDB."
    draft_help = "Retrieve auto-annotated 'draft' genomes from PECAAN."
    final_help = "Retrieve manually-annotated 'final' genomes from PhagesDB."
    genbank_help = "Retrieve annotated genomes from GenBank."
    all_help = "Retrieve all types of new data."
    genbank_results_help = "Store results of Genbank record retrieval."
    force_download_help = "Retrieve genomes regardless of date in database."
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=retrieve_help)
    parser.add_argument("database", type=str, help=database_help)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER), help=output_folder_help)
    parser.add_argument("-u", "--updates", action="store_true",
        default=False, help=updates_help)
    parser.add_argument("-d", "--draft", action="store_true",
        default=False, help=draft_help)
    parser.add_argument("-f", "--final", action="store_true",
        default=False, help=final_help)
    parser.add_argument("-g", "--genbank", action="store_true",
        default=False, help=genbank_help)
    parser.add_argument("-a", "--all_data", action="store_true",
        default=False, help=all_help)
    parser.add_argument("-gr", "--genbank_results", action="store_true",
        default=False, help=genbank_results_help)
    parser.add_argument("-fd", "--force_download", action="store_true",
        default=False, help=force_download_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if (args.updates == False and
            args.draft == False and
            args.final == False and
            args.genbank == False):
        args.all_data = True

    if args.all_data == True:
        args.updates = True
        args.draft = True
        args.final = True
        args.genbank = True

    if args.genbank == False:
        args.genbank_results = False

    return args

# TODO unittest.
def main(unparsed_args_list):
    """Run main retrieve_updates pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    force = args.force_download
    args.output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    working_dir = pathlib.Path(RESULTS_FOLDER)
    working_path = basic.make_new_dir(args.output_folder, working_dir,
                                      attempt=50)

    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]
    ncbi_creds = config["ncbi"]

    # Verify database connection and schema compatibility.
    print("Preparing genome data sets from the MySQL database...")
    alchemist = AlchemyHandler(database=args.database,
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the get_data pipeline")

    # Get existing data from MySQL to determine what needs to be updated.
    query = ("SELECT PhageID, Name, HostGenus, Status, Cluster, "
             "DateLastModified, Accession, RetrieveRecord, Subcluster, "
             "AnnotationAuthor FROM phage")

    mysqldb_genome_list =  mysqldb.parse_genome_data(engine=engine,
                                                     phage_query=query,
                                                     gnm_type="mysqldb")
    engine.dispose()
    mysqldb_genome_dict = {}
    for gnm in mysqldb_genome_list:
        # With default date, the date of all records retrieved will be newer.
        if force:
            gnm.date = constants.EMPTY_DATE
        mysqldb_genome_dict[gnm.id] = gnm

    # Get data from PhagesDB
    if (args.updates or args.final or args.draft) is True:
        print("Retrieving data from PhagesDB...")
        phagesdb_phages = phagesdb.get_phagesdb_data(constants.API_SEQUENCED)
        phagesdb_phages_dict = basic.convert_list_to_dict(phagesdb_phages,
                                                          "phage_name")
        phagesdb_genome_dict = phagesdb.parse_genomes_dict(
                                    phagesdb_phages_dict,
                                    gnm_type="phagesdb",
                                    seq=False)

        # Exit if all phage data wasn't retrieved.
        if len(phagesdb_genome_dict) == 0:
            sys.exit(1)

        # Returns a list of tuples.
        tup = match_genomes(mysqldb_genome_dict, phagesdb_genome_dict)
        matched_genomes = tup[0]
        unmatched_phagesdb_ids = tup[1]

    if args.updates is True:
        get_update_data(working_path, matched_genomes)
    if args.final is True:
        get_final_data(working_path, matched_genomes)
    if args.genbank is True:
        get_genbank_data(working_path, mysqldb_genome_dict,
                         ncbi_creds, args.genbank_results, force=force)
    if args.draft is True:
        if force:
            # Add all draft genomes currently in database to the list of
            # draft genomes to be downloaded.
            drafts = get_matched_drafts(matched_genomes)
            unmatched_phagesdb_ids |= drafts
        get_draft_data(working_path, unmatched_phagesdb_ids)


# TODO unittest.
def get_matched_drafts(matched_genomes):
    """Generate a list of matched 'draft' genomes."""
    # matched data = list of gnm_pair objects
    # unmatched data = set of PhageIDs
    drafts = set()
    for gnm_pair in matched_genomes:
        mysqldb_gnm = gnm_pair.genome1
        phagesdb_gnm = gnm_pair.genome2
        if mysqldb_gnm.annotation_status == "draft":
            # Since draft phage ids are used by PECAAN to retrieve data from
            # PhagesDB, use the id from the PhagesDB genome object.
            drafts.add(phagesdb_gnm.id)
    return drafts

# TODO unittest.
def match_genomes(dict1, dict2):
    """Match MySQL database genome data to PhagesDB genome data.

    Both dictionaries:
    Key = PhageID
    Value = pdm_utils genome object"""

    # Generate phage_id sets and match sets.
    d2_keys = dict2.keys()
    d1_keys = dict1.keys()
    matched_keys = d1_keys & d2_keys
    d1_unmatched_keys = d1_keys - d2_keys
    d2_unmatched_keys = d2_keys - d1_keys

    matched_genomes = []
    for key in matched_keys:
        gnm_pair = genomepair.GenomePair()
        gnm_pair.genome1 = dict1[key]
        gnm_pair.genome2 = dict2[key]
        matched_genomes.append(gnm_pair)

    # Only unmatched with AnnotationAuthor = 1
    unmatched_d1_genomes = {}
    for key in d1_unmatched_keys:
        gnm = dict1[key]
        if gnm.annotation_author == 1:
            unmatched_d1_genomes[key] = gnm

    results = {"match_count": len(matched_keys),
               "d1_unmatch_count": len(d1_unmatched_keys),
               "d2_unmatch_count": len(d2_unmatched_keys),
               "d1_unmatch_aa1": unmatched_d1_genomes.keys()
               }
    print_match_results(results)

    return (matched_genomes, d2_unmatched_keys)


# TODO unittest
def print_match_results(dict):
    """Print results of genome matching."""
    print("\nSummary of genome matching:")
    print(f"{dict['match_count']:>6}: genome(s) matched.")
    print(f"{dict['d1_unmatch_count']:>6}: MySQL genome(s) not matched.")
    print(f"{dict['d2_unmatch_count']:>6}: PhagesDB genome(s) not matched.")
    count = len(dict["d1_unmatch_aa1"])
    if count > 0:
        print(f"{count} MySQL genome(s) that are unexpectedly unmatched:")
        for name in dict["d1_unmatch_aa1"]:
            print(name)


# TODO unittest
def get_update_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve field updates from PhagesDB."""
    updates_folder = pathlib.Path(output_folder, UPDATES_FOLDER)
    updates_folder.mkdir()
    update_tickets = []
    for gnm_pair in matched_genomes:
        tkt_list = compare_data(gnm_pair)
        update_tickets.extend(tkt_list)

    # Field updates
    if len(update_tickets) > 0:
        print(f"\n\n{len(update_tickets)} field updates are available.")
        filepath = pathlib.Path(updates_folder, UPDATE_TABLE)
        fileio.export_data_dict(update_tickets, filepath, UPDATE_COLUMNS,
                               include_headers=True)
    else:
        print("\n\nNo field updates.")

    # Now remove empty folders.
    if len(basic.identify_contents(updates_folder, kind=None)) == 0:
        updates_folder.rmdir()



# TODO test.
def compare_data(gnm_pair):
    """Compare data and create update tickets."""
    tkt_list = []
    mysqldb_gnm = gnm_pair.genome1
    phagesdb_gnm = gnm_pair.genome2

    # Compare Cluster
    if mysqldb_gnm.cluster != phagesdb_gnm.cluster:
        # "Singleton" is not a valid cluster in MySQL (should be NULL)
        if phagesdb_gnm.cluster == "Singleton":
            # Convert to "NULL" so pipelines.update_field.update_field()
            # can convert to None
            phagesdb_gnm.cluster = "NULL"
        d1 = create_update_ticket("Cluster", phagesdb_gnm.cluster,
                                  mysqldb_gnm.id)
        tkt_list.append(d1)

    # Compare Subcluster
    if mysqldb_gnm.subcluster != phagesdb_gnm.subcluster:
        d2 = create_update_ticket("Subcluster", phagesdb_gnm.subcluster,
                                  mysqldb_gnm.id)
        tkt_list.append(d2)

    # Compare Host genus
    if mysqldb_gnm.host_genus != phagesdb_gnm.host_genus:
        d3 = create_update_ticket("HostGenus", phagesdb_gnm.host_genus,
                                  mysqldb_gnm.id)
        tkt_list.append(d3)

    # Compare Accession
    # If the genome author is not 1 ("hatfull"), then don't worry about
    # updating the accession. This used to be determined with
    # the status field, but now it is determined with the
    # AnnotationAuthor field.
    if (mysqldb_gnm.accession != phagesdb_gnm.accession and \
            mysqldb_gnm.annotation_author == 1):
        d4 = create_update_ticket("Accession", phagesdb_gnm.accession,
                                  mysqldb_gnm.id)
        tkt_list.append(d4)

    return tkt_list

# TODO test.
def create_update_ticket(field, value, key_value):
    """Create update ticket."""
    dict = {"table":"phage",
            "field":field,
            "value":value,
            "key_name":"PhageID",
            "key_value":key_value}
    return dict


# TODO unittest
def convert_tickets_to_dict(list_of_tickets):
    """Convert list of tickets to list of dictionaries."""
    l = []
    for tkt in list_of_tickets:
        dict = {
            "type": tkt.type,
            "phage_id": tkt.phage_id,
            "host_genus": tkt.data_dict["host_genus"],
            "cluster": tkt.data_dict["cluster"],
            "subcluster": tkt.data_dict["subcluster"],
            "annotation_status": tkt.data_dict["annotation_status"],
            "annotation_author": tkt.data_dict["annotation_author"],
            "description_field": tkt.description_field,
            "accession": tkt.data_dict["accession"],
            "eval_mode": tkt.eval_mode,
            "retrieve_record": tkt.data_dict["retrieve_record"]
            }
        l.append(dict)
    return l


# TODO unittest
def get_final_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve 'final' genomes from PhagesDB."""

    print(f"\n\nDownloading genome(s) from PhagesDB.")
    phagesdb_folder = pathlib.Path(output_folder, PHAGESDB_FOLDER)
    phagesdb_folder.mkdir()
    genome_folder = pathlib.Path(phagesdb_folder, GENOME_FOLDER)
    genome_folder.mkdir()
    import_tickets = []
    failed_list = []

    # Iterate through each phage in the MySQL database
    for gnm_pair in matched_genomes:
        mysqldb_gnm = gnm_pair.genome1
        phagesdb_gnm = gnm_pair.genome2

        # Not all phages have associated Genbank-formatted files
        # available on PhagesDB. Check to see if there is a flatfile for
        # this phage. Download the flatfile only if there is a date tag,
        # and only if that date is more recent than the date stored in
        # the MySQL database for that genome. The tagged date only reflects when
        # the file was uploaded into PhagesDB. The date the actual
        # Genbank record was created is stored within the file,
        # and this too could be less recent than the current version in
        # the MySQL database; however, this part gets checked during the import
        # stage.
        set_phagesdb_gnm_date(phagesdb_gnm)
        set_phagesdb_gnm_file(phagesdb_gnm)
        if (phagesdb_gnm.filename != "" and phagesdb_gnm.date > mysqldb_gnm.date):
            # Save the file on the hard drive with the same name as
            # stored on PhagesDB
            flatfile_data = phagesdb.retrieve_url_data(phagesdb_gnm.filename)
            if flatfile_data == "":
                failed_list.append(mysqldb_gnm.id)
            else:
                save_phagesdb_file(flatfile_data, phagesdb_gnm, genome_folder)
                tkt = create_phagesdb_ticket(mysqldb_gnm.id)
                import_tickets.append(tkt)

    if len(import_tickets) > 0:
        print(f"\n\n{len(import_tickets)} genome(s) "
              "were retrieved from PhagesDB.")
        create_ticket_table(import_tickets, phagesdb_folder)

    if len(failed_list) > 0:
        print(f"{len(failed_list)} genome(s) failed to be retrieved:")
        for element in failed_list:
            print(element)
        input("\n\nPress ENTER to continue.")

    # Now remove empty folders.
    if len(basic.identify_contents(genome_folder, kind=None)) == 0:
        genome_folder.rmdir()
    if len(basic.identify_contents(phagesdb_folder, kind=None)) == 0:
        phagesdb_folder.rmdir()


# TODO unittest.
def save_phagesdb_file(data, gnm, output_folder):
    """Save file retrieved from PhagesDB."""
    filename = gnm.filename.split("/")[-1]
    filepath = pathlib.Path(output_folder, filename)
    with filepath.open("w") as fh:
        fh.write(data)

# TODO test.
def save_pecaan_file(response, name, output_folder):
    """Save data retrieved from PECAAN."""
    filename = f"{name}.txt"
    filepath = pathlib.Path(output_folder, filename)
    with filepath.open("w") as fh:
        fh.write(response)


# TODO test.
def save_genbank_file(seqrecord, accession, name, output_folder):
    """Save retrieved record to file."""
    filename = f"{name}__{accession}.txt"
    filepath = pathlib.Path(output_folder, filename)
    SeqIO.write(seqrecord, str(filepath), "genbank")


# TODO unittest.
def set_phagesdb_gnm_date(gnm):
    """Set the date of a PhagesDB genome object."""

    # Since there may be multiple 'dates' associated with a genome in PhagesDB,
    # this date attribute can't be directly set when parsing the data.
    # Some phages may have a file but no associated
    # date tagged with that file (since date tagging has not
    # always been implemented). If there is a date, it is
    # formatted as: '2017-02-15T10:37:21Z'. If there is no date,
    # it is Null, but change this to a standarized 'empty' value of 1/1/0001.
    date = gnm.misc["qced_genbank_file_date"]
    if date is None:
        # The date used to be set to default (1/1/0001), but the force flag
        # now relies on this early date, so here set to a date
        # more recent than default.
        # gnm.date = constants.EMPTY_DATE
        gnm.date = datetime.strptime("1/1/1000 00:00:00", "%m/%d/%Y %H:%M:%S")
    else:
        #e.g. 2018-09-01T22:19:33-04:00 (string)
        date = date.split("T")[0]
        #e.g. 2018-09-01 (string)
        gnm.date = datetime.strptime(date, "%Y-%m-%d")
        #e.g. 2018-09-01 00:00:00 (datetime.datetime)

# TODO unittest.
def set_phagesdb_gnm_file(gnm):
    """Set the filename of a PhagesDB genome object."""
    # Since there may be multiple 'files'' associated with a genome in PhagesDB,
    # the set filename is not necessarily for the flat file.
    file_url = gnm.misc["qced_genbank_file"]
    if file_url is None:
        gnm.filename = ""
    else:
        gnm.filename = file_url


# TODO unittest.
def get_genbank_data(output_folder, genome_dict, ncbi_cred_dict={},
                     genbank_results=False, force=False):
    """Run sub-pipeline to retrieve genomes from GenBank."""
    # Flow of the NCBI record retrieval process:
    # 1 Create list of phages to check for updates at NCBI (completed above)
    # 2 Using esearch, verify which accessions are valid
    # 3 Using esummary, get update date for each valid accession
    # 4 Using efetch, retrieve flat files for NCBI records newer than
    # the MySQL database date
    # 5 Save new records in a folder and create an import table for them

    print(f"\n\nDownloading genome(s) from GenBank.")
    # Create output folder
    ncbi_folder = pathlib.Path(output_folder, GENBANK_FOLDER)
    ncbi_folder.mkdir()
    ncbi_results_list = []

    # Iterate through each phage in the MySQL database
    tup1 = sort_by_accession(genome_dict, force=force)
    ncbi_results_list.extend(tup1[0])
    accession_dict = tup1[1]

    # More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    ncbi.set_entrez_credentials(tool=ncbi_cred_dict["tool"],
                                email=ncbi_cred_dict["email"],
                                api_key=ncbi_cred_dict["api_key"])

    results = retrieve_records(accession_dict, ncbi_folder, batch_size=200)
    ncbi_results_list.extend(results)

    # Record retrieval results for all phages.
    if genbank_results == True:
        output_genbank_summary(ncbi_folder, ncbi_results_list)

    # Print summary of script
    tallies = compute_genbank_tallies(ncbi_results_list)
    print_genbank_tallies(tallies)

    # Now remove empty folders.
    if len(basic.identify_contents(ncbi_folder, kind=None)) == 0:
        ncbi_folder.rmdir()


# TODO unittest.
def output_genbank_summary(output_folder, results):
    """Save summary of GenBank retrieval results to file."""
    filepath = pathlib.Path(output_folder, GENBANK_RESULTS_TABLE)
    fileio.export_data_dict(results, filepath, NCBI_RESULTS_COLUMNS,
                           include_headers=True)


# TODO unittest.
def compute_genbank_tallies(results):
    """Tally results from GenBank retrieval."""
    dict = {}
    for s in GENBANK_STATUS:
        dict[s] = 0
    for d in results:
        r = d["result"]
        dict[r] += 1
    dict["total"] = len(results)
    dict["auto_updated"] = dict["total"] - dict["not_auto_updated"]
    dict["accession"] = (dict["auto_updated"] -
                         dict["no_accession"] -
                         dict["duplicate_accession"])
    return dict


# TODO unittest.
def print_genbank_tallies(tallies):
    """Print results of GenBank retrieval."""

    print("\n\n\nSummary of GenBank data retrieval:")
    print("Of the genomes in the MySQL database:")
    print(f"{tallies['total']:>6}: total")
    print(f"{tallies['not_auto_updated']:>6}: not auto-updated")
    print(f"{tallies['auto_updated']:>6}: auto-updated")

    print("\nOf the auto-updated genomes:")
    print(f"{tallies['no_accession']:>6}: no accession")
    print(f"{tallies['duplicate_accession']:>6}: duplicated accession")
    print(f"{tallies['accession']:>6}: unique accession")

    print("\nOf the auto-updated genomes with unique accessions:")
    print(f"{tallies['retrieval_failure']:>6}: could not be retrieved")
    print(f"{tallies['docsum_not_new']:>6}: retrieved but docsum not new")
    print(f"{tallies['record_not_new']:>6}: retrieved but record not new")
    print(f"{tallies['retrieved_for_import']:>6}: retrieved for import")


# TODO unittest.
def sort_by_accession(genome_dict, force=False):
    """Sort genome objects based on their accession status.

    Only retain data if genome is set to be automatically updated,
    there is a valid accession, and the accession is unique.
    """
    results = []
    accession_dict = {}
    check_set = {None, "", "none"}
    # Get a list of duplicated accessions.
    unique, duplicate = create_accession_sets(genome_dict)

    # Iterate through each phage in the MySQL database
    for key in genome_dict.keys():
        status = None
        gnm = genome_dict[key]
        if gnm.retrieve_record != 1 and force == False:
            status = NOT_AUTO
        elif gnm.accession in check_set:
            status = NO_ACC
        elif gnm.accession in duplicate:
            status = DUPE_ACC
        else:
            pass

        if status is None:
            # Dictionary of phage data based on unique accessions
            # Key = accession
            # Value = genome object
            accession_dict[gnm.accession] = gnm
        else:
            result = create_results_dict(gnm, "NA", status)
            results.append(result)

    return (results, accession_dict)


# TODO unittest.
def process_failed_retrieval(accession_list, accession_dict):
    """Create list of dictionaries for records that could not be retrieved."""
    l = []
    for acc in accession_list:
        gnm = accession_dict[acc]
        result = create_results_dict(gnm, "NA", RET_FAIL)
        l.append(result)
    return l


# TODO unittest.
def retrieve_records(accession_dict, ncbi_folder, batch_size=200):
    """Retrieve GenBank records."""
    print("\n\nRetrieving records from NCBI")
    genome_folder = pathlib.Path(ncbi_folder, GENOME_FOLDER)
    genome_folder.mkdir()
    retrieval_errors = []
    results = []
    tickets_list = []
    accessions = list(accession_dict.keys())
    mod_accessions = [accession + "[ACCN]" for accession in accessions]

    # When retrieving in batch sizes, first create the list of values
    # indicating which indices of the accessions should be used
    # to create each batch.
    # For instace, if there are five accessions, batch size of two produces
    # indices = 0,2,4

    batch_indices = basic.create_indices(mod_accessions, batch_size)
    print(f"There are {len(mod_accessions)} GenBank accession(s) to check.")
    for indices in batch_indices:
        start = indices[0]
        stop = indices[1]
        print(f"Checking accessions {start + 1} to {stop}...")
        esearch_term = " | ".join(mod_accessions[start:stop])

        # Use esearch for each accession
        # First use esearch to verify the accessions are valid.
        search_record = ncbi.run_esearch(db="nucleotide", term=esearch_term,
                                         usehistory="y")
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]

        # Keep track of the accessions that failed to be located in NCBI
        # Each accession in the error list is formatted "accession[ACCN]"
        current_batch_size = stop - start
        if search_count < current_batch_size:
            search_failure = search_record["ErrorList"]["PhraseNotFound"]
            for accession in search_failure:
                retrieval_errors.append(accession[:-6])

        # Now get summaries for these records using esummary
        summary_records = ncbi.get_summaries(db="nucleotide",
                                             query_key=search_query_key,
                                             webenv=search_webenv)
        results_tuple = get_accessions_to_retrieve(summary_records,
                                                   accession_dict)
        accessions_to_retrieve = results_tuple[0]
        results.extend(results_tuple[1])

        if len(accessions_to_retrieve) > 0:
            # Use efetch to retrieve the record.
            output_list = ncbi.get_records(accessions_to_retrieve,
                                           db="nucleotide",
                                           rettype="gb",
                                           retmode="text")

            # TODO check_record_date may be redundant. It checks date within the
            # record. Earlier in the pipeline, the docsum date has already been
            # checked though. So if docsum date is identical to date in the
            # record, this is redundant.
            tup = check_record_date(output_list, accession_dict)
            new_record_list = tup[0]
            # list of results dictionaries
            results.extend(tup[1])

            if len(new_record_list) > 0:
                tickets = save_and_tickets(new_record_list, accession_dict,
                                           genome_folder)
                tickets_list.extend(tickets)

    if len(tickets_list) > 0:
        create_ticket_table(tickets_list, ncbi_folder)

    # Remove genome folder if empty.
    if len(basic.identify_contents(genome_folder, kind=None)) == 0:
        genome_folder.rmdir()

    # Report the genomes that could not be retrieved.
    failed = process_failed_retrieval(retrieval_errors, accession_dict)
    results.extend(failed)

    return results


# TODO unittest.
def get_accessions_to_retrieve(summary_records, accession_dict):
    """Review GenBank summary to determine which records are new."""
    accessions = []
    results = []
    for doc_sum in summary_records:
        doc_sum_accession = doc_sum["Caption"]
        doc_sum_date = datetime.strptime(doc_sum["UpdateDate"], "%Y/%m/%d")
        gnm = accession_dict[doc_sum_accession]

        # If Document Summary date is newer than the MySQL database date,
        # or if the genome being evaluated is currently draft status but
        # we were able to retrieve a Genbank record (i.e. it's not a
        # draft anymore), append the accession to the list of accessions
        # from this batch to retrieve.
        # Otherwise, record that the GenBank record isn't new.
        if doc_sum_date > gnm.date or gnm.annotation_status == "draft":
            accessions.append(doc_sum_accession)
        else:
            result = create_results_dict(gnm, doc_sum_date, OLD_DOCSUM)
            results.append(result)
    return (accessions, results)


# TODO unittest.
def check_record_date(record_list, accession_dict):
    """Check whether the GenBank record is new."""
    results = []
    new_record_list = []
    # For any records that were retrieved, mark their data in the NCBI
    # results file, and save their records as Genbank files, and create
    # import table entries for them.
    for index in range(len(record_list)):
        record = record_list[index]
        date = record.annotations["date"]
        date = datetime.strptime(date, '%d-%b-%Y')
        accession = record.annotations["accessions"][0]
        # accession = record.name
        # accession = accession.split('.')[0]
        gnm = accession_dict[accession]

        # Save new records in a folder and create an import table row
        # for them. If the genome is currently a draft annotation, create
        # an import ticket for replacement regardless of the date in
        # the Genbank record. This ensures that if a user fails to
        # upload a manual annotation to PhagesDB, once the Genbank
        # accession becomes active MySQL database will get the new version.
        # This should always happen, since we've only retrieved new records.
        if (date > gnm.date or gnm.annotation_status == "draft"):
            data_dict = create_results_dict(gnm, date, RETRIEVED)
            results.append(data_dict)
            new_record_list.append(record)
            # Determine what the status of the genome will be.
            # If the genome in the MySQL database was already 'final' or 'unknown',
            # then keep the status unchanged. If the genome in the MySQL database
            # was 'draft', the status should be changed to 'final'.
            if gnm.annotation_status == "draft":
                gnm.annotation_status = "final"
        else:
            data_dict = create_results_dict(gnm, date, OLD_RECORD)
            results.append(data_dict)

    return (new_record_list, results)


# TODO unittest.
def save_and_tickets(record_list, accession_dict, output_folder):
    """Save flat files retrieved from GenBank and create import tickets."""
    # Both tasks accomplished at the same time to
    # iterate through record list once.
    l = []
    for record in record_list:
        # accession = seqrecord.name
        # accession = accession.split('.')[0]
        accession = record.annotations["accessions"][0]
        gnm = accession_dict[accession]
        save_genbank_file(record, accession, gnm.name, output_folder)
        tkt = create_genbank_ticket(gnm)
        l.append(tkt)
    return l

# TODO test.
def create_ticket_table(tickets, output_folder):
    """Save tickets associated with retrieved from GenBank files."""
    filepath = pathlib.Path(output_folder, IMPORT_TABLE)
    tickets = convert_tickets_to_dict(tickets)
    fileio.export_data_dict(tickets, filepath, IMPORT_COLUMNS,
                           include_headers=True)


# TODO unittest.
def create_accession_sets(genome_dict):
    """Generate set of unique and non-unique accessions.

    Input is a dictionary of pdm_utils genome objects."""
    l = []
    for id in genome_dict.keys():
        gnm = genome_dict[id]
        # In the MySQL databse, empty accession data is stored as ''.
        # There are no NULL accessions.
        if gnm.accession != "":
            l.append(gnm.accession)
    unique, duplicated = basic.identify_unique_items(l)

    if len(duplicated) > 0:
        print("There are duplicated accessions. Some data will not be "
             "retrieved from GenBank:")
        for duplicate in duplicated:
            print(duplicate)
        input("\n\nPress ENTER to continue.")

    return unique, duplicated


# TODO unittest.
def create_results_dict(gnm, genbank_date, result):
    """Create a dictionary of data summarizing NCBI retrieval status."""
    d = {
        "phage_id": gnm.id,
        "phage_name": gnm.name,
        "accession": gnm.accession,
        "annotation_status": gnm.annotation_status,
        "mysql_date": gnm.date,
        "genbank_date": genbank_date,
        "result": result
        }
    return d


# TODO unittest.
def get_draft_data(output_path, phage_id_set):
    """Run sub-pipeline to retrieve auto-annotated 'draft' genomes."""

    if len(phage_id_set) > 0:
        print(f"\n\nRetrieving {len(phage_id_set)} genome(s) from PECAAN")
        phage_id_list = list(phage_id_set)
        pecaan_folder = pathlib.Path(output_path, PECAAN_FOLDER)
        pecaan_folder.mkdir()
        retrieve_drafts(pecaan_folder, phage_id_list)
    else:
        print("No new genomes to retrieve from PECAAN.")


# TODO unittest.
def retrieve_drafts(output_folder, phage_list):
    """Retrieve auto-annotated 'draft' genomes from PECAAN."""

    genome_folder = pathlib.Path(output_folder, GENOME_FOLDER)
    genome_folder.mkdir()
    failed = []
    tickets = []

    # Iterate through each row in the file
    for new_phage in phage_list:
        pecaan_link = constants.PECAAN_PREFIX + new_phage
        response = phagesdb.retrieve_url_data(pecaan_link)
        if response == "":
            print(f"Error: unable to retrieve {new_phage} draft genome.")
            print(pecaan_link)
            failed.append(new_phage)
        else:
            save_pecaan_file(response, new_phage, genome_folder)
            tkt = create_draft_ticket(new_phage)
            tickets.append(tkt)
            print(f"{new_phage} retrieved from PECAAN.")

    if len(tickets) > 0:
        create_ticket_table(tickets, output_folder)
        print(f"{len(tickets)} phage(s) were successfully retrieved")

    if len(failed) > 0:
        print(f"{len(failed)} phage(s) failed to be retrieved:")
        for item in failed:
            print(item)
        input("\n\nPress ENTER to continue.")


# TODO test.
def create_draft_ticket(name):
    """Create ImportTicket for draft genome."""
    tkt = ticket.ImportTicket()
    tkt.type = "add"
    tkt.phage_id = name
    tkt.description_field = "product"
    tkt.eval_mode = "draft"
    tkt.data_dict = {
        "host_genus": "retrieve",
        "cluster": "retrieve",
        "subcluster": "retrieve",
        "annotation_status": "draft",
        "annotation_author": 1,
        "accession": "none",
        "retrieve_record": 1
        }
    return tkt


# TODO test.
def create_phagesdb_ticket(phage_id):
    """Create ImportTicket for PhagesDB genome."""

    # Since the PhagesDB phage has been matched to
    # the MySQL database phage, the AnnotationAuthor field
    # could be assigned from the current mysqldb author
    # variable. However, since this genbank-formatted
    # file is acquired through PhagesDB, both the
    # Annotation status is expected to be 'final' and
    # the Annotation author is expected to be 'hatfull'.
    tkt = ticket.ImportTicket()
    tkt.type = "replace"
    tkt.phage_id = phage_id
    tkt.description_field = "product"
    tkt.eval_mode = "final"
    tkt.data_dict = {
        "host_genus": "retain", # formerly "retrieve",
        "cluster": "retain", # formerly "retrieve",
        "subcluster": "retain", # formerly "retrieve",
        "annotation_status": "final",
        "annotation_author": 1,
        "accession": "retain", # formerly "retrieve",
        "retrieve_record": 1
        }
    return tkt


# TODO test.
def create_genbank_ticket(gnm):
    """Create ImportTicket for GenBank record."""
    # Accession is set to 'parse' to ensure that during import,
    # the file's accession is directly compared to the database
    # record's accession.
    tkt = ticket.ImportTicket()
    tkt.type = "replace"
    tkt.phage_id = gnm.id
    tkt.description_field = "product"
    tkt.eval_mode = "auto"
    tkt.data_dict = {
        "host_genus": "retain", # formerly gnm.host_genus,
        "cluster": "retain", # formerly gnm.cluster,
        "subcluster": "retain", # formerly gnm.subcluster,
        "annotation_status": "retain", # formerly gnm.annotation_status,
        "annotation_author": "retain", # formerly gnm.annotation_author,
        "accession": "parse",
        "retrieve_record": "retain", # formerly 1
        }
    return tkt
