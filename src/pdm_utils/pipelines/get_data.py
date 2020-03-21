"""Pipeline to gather new data to be imported into a MySQL database."""

import argparse
import csv
from datetime import datetime
import pathlib
import sys
import time
from Bio import SeqIO
from pdm_utils.classes import genomepair
from pdm_utils.classes import ticket
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import mysqldb
from pdm_utils.functions import phagesdb
from pdm_utils.functions import tickets

# TODO toggles whether both ticket table files are generated.
BOTH = True

# TODO update column headers for import table
# Old format
IMPORT_COLUMNS1 = ["type",
                   "phage_id",
                   "host_genus",
                   "cluster",
                   "subcluster",
                   "annotation_status",
                   "annotation_author",
                   "description_field",
                   "accession",
                   "eval_mode",
                   "secondary_phage_id"]

# New format
IMPORT_COLUMNS2 = constants.IMPORT_TABLE_STRUCTURE["order"]

# Columns for new format for update tickets to match RandomFieldUpdateHandler
UPDATE_COLUMNS = ["table",
                  "field",
                  "value",
                  "key_name",
                  "key_value"]

# Column headers for NCBI results output file.
NCBI_RESULTS_COLUMNS = ["phage_id",
                       "phage_name",
                       "accession",
                       "annotation_status",
                       "mysql_date",
                       "genbank_date",
                       "result"]

GENOMES_DIR = "genomes"

# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""

    RETRIEVE_HELP = ("Pipeline to retrieve new data to import into a "
                     "MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = "Path to the directory where updates will be stored."
    UPDATES_HELP = ("Retrieve updates to HostGenus, Cluster, "
                    "Subcluster, and Accession field data from PhagesDB.")
    DRAFT_HELP = "Retrieve auto-annotated 'draft' genomes from PECAAN."
    FINAL_HELP = ("Retrieve new manually-annotated 'final' "
                  "genomes from PhagesDB.")
    GENBANK_HELP = "Retrieve revised annotated genomes from GenBank."
    ALL_HELP = "Retrieve all types of new data."
    NCBI_CRED_FILE_HELP = "Path to the file containing NCBI credentials."
    GENBANK_RESULTS_HELP = "Store results of Genbank record retrieval."



    parser = argparse.ArgumentParser(description=RETRIEVE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("output_folder", type=pathlib.Path,
        help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-u", "--updates", action="store_true",
        default=False, help=UPDATES_HELP)
    parser.add_argument("-d", "--draft", action="store_true",
        default=False, help=DRAFT_HELP)
    parser.add_argument("-f", "--final", action="store_true",
        default=False, help=FINAL_HELP)
    parser.add_argument("-g", "--genbank", action="store_true",
        default=False, help=GENBANK_HELP)
    parser.add_argument("-a", "--all_data", action="store_true",
        default=False, help=ALL_HELP)
    parser.add_argument("-c", "--ncbi_credentials_file", type=pathlib.Path,
        help=NCBI_CRED_FILE_HELP)
    parser.add_argument("-gr", "--genbank_results", action="store_true",
        default=False, help=GENBANK_RESULTS_HELP)


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
        args.ncbi_credentials_file = None
        args.genbank_results = False

    return args


# TODO unittest.
def main(unparsed_args_list):
    """Run main retrieve_updates pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    date = time.strftime("%Y%m%d")

    args.output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)

    working_dir = pathlib.Path(f"{date}_get_data")
    working_path = basic.make_new_dir(args.output_folder, working_dir,
                                      attempt=10)

    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    ncbi_cred_dict = ncbi.get_ncbi_creds(args.ncbi_credentials_file)

    # Verify database connection and schema compatibility.
    print("Preparing genome data sets from the MySQL database...")
    engine = mysqldb.connect_to_db(args.database)
    mysqldb.check_schema_compatibility(engine, "the get_data pipeline")

    # Get existing data from MySQL to determine what needs to be updated.
    query = ("SELECT PhageID, Name, HostGenus, Status, Cluster, "
             "DateLastModified, Accession, RetrieveRecord, Subcluster, "
             "AnnotationAuthor FROM phage")

    mysqldb_genome_list =  mysqldb.parse_genome_data(
                       engine=engine,
                       phage_query=query,
                       gnm_type="mysqldb")
    engine.dispose()
    mysqldb_genome_dict = {}
    for gnm in mysqldb_genome_list:
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
        match_output = match_genomes(mysqldb_genome_dict, phagesdb_genome_dict)
        matched_genomes = match_output[0]
        unmatched_phagesdb_ids = match_output[1]

    if args.updates is True:
        get_update_data(working_path, matched_genomes)
    if args.final is True:
        get_final_data(working_path, matched_genomes)
    if args.genbank is True:
        get_genbank_data(working_path, mysqldb_genome_dict,
                         ncbi_cred_dict, args.genbank_results)
    if args.draft is True:
        get_draft_data(working_path, unmatched_phagesdb_ids)
    print("\n\n\nRetrieve updates script completed.")



# TODO unittest.
def match_genomes(mysqldb_dict, phagesdb_dict):
    """Match MySQL database genome data to PhagesDB genome data.

    Both dictionaries:
    Key = PhageID
    Value = pdm_utils genome object"""

    # Generate phage_id sets and match sets.
    phagesdb_ids = phagesdb_dict.keys()
    mysqldb_ids = mysqldb_dict.keys()
    matched_ids = mysqldb_ids & phagesdb_ids
    unmatched_mysqldb_ids = mysqldb_ids - phagesdb_ids
    unmatched_phagesdb_ids = phagesdb_ids - mysqldb_ids

    matched_genomes = []
    for id in matched_ids:
        gnm_pair = genomepair.GenomePair()
        gnm_pair.genome1 = mysqldb_dict[id]
        gnm_pair.genome2 = phagesdb_dict[id]
        matched_genomes.append(gnm_pair)

    unmatched_mysqldb_authored_genomes = {}
    for id in unmatched_mysqldb_ids:
        gnm = mysqldb_dict[id]
        if gnm.annotation_author == 1:
            unmatched_mysqldb_authored_genomes[id] = gnm

    print("\nSummary of genome matching:")
    print(f"{len(matched_ids):>6}: genome(s) matched.")
    print(f"{len(unmatched_mysqldb_ids):>6}: MySQL genome(s) not matched.")
    print(f"{len(unmatched_phagesdb_ids):>6}: PhagesDB genome(s) not matched.")

    count = len(unmatched_mysqldb_authored_genomes.keys())
    if count > 0:
        print(f"{count} Hatfull-authored unmatched MySQL genome(s):")
        for key in unmatched_mysqldb_authored_genomes.keys():
            print(key)

    return (matched_genomes, unmatched_phagesdb_ids)


# TODO unittest
def get_update_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve field updates from PhagesDB."""
    update_tickets = []
    for gnm_pair in matched_genomes:
        mysqldb_gnm = gnm_pair.genome1
        phagesdb_gnm = gnm_pair.genome2


        # Compare Cluster
        if mysqldb_gnm.cluster != phagesdb_gnm.cluster:
            result1 = {
                      "table":"phage",
                      "field":"Cluster",
                      "value":phagesdb_gnm.cluster,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result1)


        # Compare Subcluster
        if mysqldb_gnm.subcluster != phagesdb_gnm.subcluster:
            result3 = {
                      "table":"phage",
                      "field":"Subcluster",
                      "value":phagesdb_gnm.subcluster,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result3)


        # Compare Host genus
        if mysqldb_gnm.host_genus != phagesdb_gnm.host_genus:
            result5 = {
                      "table":"phage",
                      "field":"HostGenus",
                      "value":phagesdb_gnm.host_genus,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result5)

        # Compare Accession
        # If the genome author is not "hatfull", then don't worry about
        # updating the accession. This used to be determined with
        # the status field, but now it is determined with the
        # AnnotationAuthor field.
        if (mysqldb_gnm.accession != phagesdb_gnm.accession and \
                mysqldb_gnm.annotation_author == 1):
            result6 = {
                      "table":"phage",
                      "field":"Accession",
                      "value":phagesdb_gnm.accession,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result6)


    # Field updates
    if len(update_tickets) > 0:
        print("\n\nNew field updates are available.")
        filepath = basic.prepare_filepath(output_folder, "update_table.csv",
                                           folder_name="updates")
        basic.export_data_dict(update_tickets, filepath,
                                   UPDATE_COLUMNS, include_headers=True)
    else:
        print("\n\nNo field updates found.")


# TODO unittest
def convert_tickets_to_dict(list_of_tickets, old_format=False):
    """Convert list of tickets to list of dictionaries."""
    list_of_data = []
    for tkt in list_of_tickets:
        tkt_data = {}
        tkt_data["type"] = tkt.type
        tkt_data["phage_id"] = tkt.phage_id
        tkt_data["host_genus"] = tkt.data_dict["host_genus"]
        tkt_data["cluster"] = tkt.data_dict["cluster"]
        tkt_data["subcluster"] = tkt.data_dict["subcluster"]
        tkt_data["annotation_status"] = tkt.data_dict["annotation_status"]

        if old_format == True:
            if tkt.data_dict["annotation_author"] == 1:
                tkt_data["annotation_author"] = "hatfull"
            else:
                tkt_data["annotation_author"] = "gbk"
        else:
            tkt_data["annotation_author"] = tkt.data_dict["annotation_author"]

        tkt_data["description_field"] = tkt.description_field
        tkt_data["accession"] = tkt.data_dict["accession"]

        if old_format == True:
            if tkt.eval_mode == "auto":
                tkt_data["eval_mode"] = "ncbi_auto"
            elif tkt.eval_mode == "draft":
                tkt_data["eval_mode"] = "pecaan"
            elif tkt.eval_mode == "final":
                tkt_data["eval_mode"] = "phagesdb"
            else:
                tkt_data["eval_mode"] = tkt.eval_mode
        else:
            tkt_data["eval_mode"] = tkt.eval_mode

        if old_format == True:
            tkt_data["secondary_phage_id"] = tkt.data_dict["secondary_phage_id"]
        else:
            tkt_data["retrieve_record"] = tkt.data_dict["retrieve_record"]
        list_of_data.append(tkt_data)
    return list_of_data


# TODO unittest
def get_final_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve 'final' genomes from PhagesDB."""

    phagesdb_folder = pathlib.Path(output_folder, "phagesdb")
    phagesdb_folder.mkdir()
    genome_folder = pathlib.Path(phagesdb_folder, GENOMES_DIR)
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
                flatfile_filename = phagesdb_gnm.filename.split("/")[-1]
                flatfile_path = pathlib.Path(genome_folder,
                                             flatfile_filename)
                with flatfile_path.open("w") as fh:
                    fh.write(flatfile_data)
                # Create the new import ticket
                # Since the PhagesDB phage has been matched to
                # the MySQL database phage, the AnnotationAuthor field
                # could be assigned from the current mysqldb author
                # variable. However, since this genbank-formatted
                # file is acquired through PhagesDB, both the
                # Annotation status is expected to be 'final' and
                # the Annotation author is expected to be 'hatfull'.
                tkt = ticket.ImportTicket()
                tkt.type = "replace"
                tkt.phage_id = mysqldb_gnm.id
                tkt.data_dict["host_genus"] = "retrieve"
                tkt.data_dict["cluster"] = "retrieve"
                tkt.data_dict["subcluster"] = "retrieve"
                tkt.data_dict["annotation_status"] = "final"
                tkt.data_dict["annotation_author"] = 1
                tkt.description_field = "product"
                tkt.data_dict["accession"] = "retrieve"
                tkt.eval_mode = "final"
                # TODO secondary_phage_id data is for old ticket format.
                tkt.data_dict["secondary_phage_id"] = mysqldb_gnm.id
                tkt.data_dict["retrieve_record"] = 1
                import_tickets.append(tkt)


    count1 = len(import_tickets)
    if count1 > 0:
        print(f"\n\n{count1} phage(s) were retrieved from PhagesDB.")
        filepath = basic.prepare_filepath(phagesdb_folder, "legacy_import_table.csv")
        import_tickets1 = convert_tickets_to_dict(import_tickets, old_format=True)
        basic.export_data_dict(import_tickets1, filepath, IMPORT_COLUMNS1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        if BOTH:
            filepath2 = basic.prepare_filepath(phagesdb_folder, "import_table.csv")
            import_tickets2 = convert_tickets_to_dict(import_tickets)
            basic.export_data_dict(import_tickets2, filepath2,
                                   IMPORT_COLUMNS2, include_headers=True)

    if len(failed_list) > 0:
        print(f"{len(failed_list)} phage(s) failed to be retrieved:")
        for element in failed_list:
            print(element)
        input("\n\nPress ENTER to continue.")

    # Now remove empty folders.
    if len(basic.identify_contents(genome_folder, kind=None)) == 0:
        genome_folder.rmdir()
    if len(basic.identify_contents(phagesdb_folder, kind=None)) == 0:
        phagesdb_folder.rmdir()


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
        gnm.date = constants.EMPTY_DATE
    else:
        date = date.split("T")[0]
        gnm.date = datetime.strptime(date, "%Y-%m-%d")



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
                     genbank_results=False):
    """Run sub-pipeline to retrieve genomes from GenBank."""
    # Flow of the NCBI record retrieval process:
    # 1 Create list of phages to check for updates at NCBI (completed above)
    # 2 Using esearch, verify which accessions are valid
    # 3 Using esummary, get update date for each valid accession
    # 4 Using efetch, retrieve flat files for NCBI records newer than
    # the MySQL database date
    # 5 Save new records in a folder and create an import table for them

    # Create output folder
    ncbi_folder = pathlib.Path(output_folder, f"genbank")
    ncbi_folder.mkdir()

    ncbi_results_list = []
    tallies = {}
    tallies["total"] = len(genome_dict.keys())

    # Iterate through each phage in the MySQL database
    result_tuple1 = sort_by_accession(genome_dict)
    tallies["not_auto_updated"] = result_tuple1[0]
    tallies["no_accession"] = result_tuple1[1]
    tallies["duplicate_accession"] = result_tuple1[2]
    ncbi_results_list.extend(result_tuple1[3])
    unique_accession_dict = result_tuple1[4]

    # More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    ncbi.set_entrez_credentials(
        tool=ncbi_cred_dict["ncbi_tool"],
        email=ncbi_cred_dict["ncbi_email"],
        api_key=ncbi_cred_dict["ncbi_api_key"])

    results_tuple2 = retrieve_records(unique_accession_dict, batch_size=200)
    tallies["docsum_not_new"] = results_tuple2[0]
    retrieved_record_list = results_tuple2[1]
    retrieval_error_list = results_tuple2[2]
    ncbi_results_list.extend(results_tuple2[3])

    # Report the genomes that could not be retrieved.
    results3 = process_failed_retrieval(retrieval_error_list, unique_accession_dict)
    ncbi_results_list.extend(results3)
    tallies["retrieval_failure"] = len(retrieval_error_list)

    results_tuple4 = check_record_date(retrieved_record_list, unique_accession_dict)
    new_record_list = results_tuple4[0]
    ncbi_results_list.extend(results_tuple4[1])

    tallies["retrieved_for_import"] = len(new_record_list)
    tallies["record_not_new"] = (len(retrieved_record_list)
                                 - len(new_record_list))

    if len(new_record_list) > 0:
        save_files_and_tkts(new_record_list, unique_accession_dict, ncbi_folder)

    # Record retrieval results for all phages.
    if genbank_results == True:
        filepath3 = basic.prepare_filepath(ncbi_folder, "genbank_results.csv")
        basic.export_data_dict(ncbi_results_list, filepath3,
                                   NCBI_RESULTS_COLUMNS, include_headers=True)

    # Print summary of script
    tallies["auto_updated"] = tallies["total"] - tallies["not_auto_updated"]
    tallies["accession"] = (tallies["auto_updated"]
                            - tallies["no_accession"]
                            - tallies["duplicate_accession"])

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

    # Now remove empty folders.
    if len(basic.identify_contents(ncbi_folder, kind=None)) == 0:
        ncbi_folder.rmdir()




# TODO unittest.
def sort_by_accession(genome_dict):
    """Sort genome objects based on their accession status."""

    tally_not_auto_updated = 0
    tally_no_accession = 0
    tally_duplicate_accession = 0
    ncbi_results_list = []
    accession_dict = {}

    # Get a list of duplicated accessions.
    unique_accessions, duplicate_accessions = create_accession_sets(genome_dict)

    # Iterate through each phage in the MySQL database
    for key in genome_dict.keys():
        gnm = genome_dict[key]
        # For NCBI retrieval, add to dictionary if
        # 1) the genome is set to be automatically updated and
        # 2) if there is an accession number
        if gnm.retrieve_record != 1:
            tally_not_auto_updated += 1
            ncbi_results = create_results_dict(gnm, "NA", "no automatic update")
            ncbi_results_list.append(ncbi_results)

        elif (gnm.accession is None or
                gnm.accession == "none" or
                gnm.accession == ""):
            tally_no_accession += 1
            ncbi_results = create_results_dict(gnm, "NA", "no accession")
            ncbi_results_list.append(ncbi_results)

        elif gnm.accession in duplicate_accessions:
            tally_duplicate_accession += 1
            ncbi_results = create_results_dict(gnm, "NA", "duplicate accession")
            ncbi_results_list.append(ncbi_results)

        else:
            # Dictionary of phage data based on unique accessions
            # Key = accession
            # Value = genome object
            accession_dict[gnm.accession] = gnm

    return (tally_not_auto_updated, tally_no_accession,
            tally_duplicate_accession, ncbi_results_list,
            accession_dict)



# TODO unittest.
def process_failed_retrieval(accession_list, accession_dict):
    """Create list of data dictionaries for records that could not be retrieved."""
    results = []
    for accession in accession_list:
        gnm = accession_dict[accession]
        result = create_results_dict(gnm, "NA", "retrieval failure")
        results.append(result)
    return results



# TODO unittest.
def retrieve_records(accession_dict, batch_size=200):
    """Retrieve GenBank records."""
    # First use esearch to verify the accessions are valid.
    # Seoncd use efetch to retrieve the record.
    print("\n\nRetrieving records from NCBI")
    retrieved_records = [] # GenBank records that have been retrieved.
    retrieval_errors = []
    tally_not_new = 0 # Keeps track if docsum date is new or not.
    results = [] # Summary of retrieval results.
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
        delimiter = " | "
        esearch_term = delimiter.join(mod_accessions[start:stop])

        # Use esearch for each accession
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
        tally_not_new += len(summary_records) - len(accessions_to_retrieve)

        if len(accessions_to_retrieve) > 0:
            output_list = ncbi.get_records(accessions_to_retrieve,
                                           db="nucleotide",
                                           rettype="gb",
                                           retmode="text")
            retrieved_records.extend(output_list)

    return (tally_not_new, retrieved_records, retrieval_errors, results)



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
            result = create_results_dict(gnm, doc_sum_date, "record not new")
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
        accession = record.name
        accession = accession.split('.')[0]
        gnm = accession_dict[accession]

        # Save new records in a folder and create an import table row
        # for them. If the genome is currently a draft annotation, create
        # an import ticket for replacement regardless of the date in
        # the Genbank record. This ensures that if a user fails to
        # upload a manual annotation to PhagesDB, once the Genbank
        # accession becomes active MySQL database will get the new version.
        # This should always happen, since we've only retrieved new records.
        if (date > gnm.date or gnm.annotation_status == "draft"):
            data_dict = create_results_dict(gnm, date, "retrieved for import")
            results.append(data_dict)
            new_record_list.append(record)
            # Determine what the status of the genome will be.
            # If the genome in the MySQL database was already 'final' or 'unknown',
            # then keep the status unchanged. If the genome in the MySQL database
            # was 'draft', the status should be changed to 'final'.
            if gnm.annotation_status == "draft":
                gnm.annotation_status = "final"
        else:
            data_dict = create_results_dict(gnm, date, "record date not new")
            results.append(data_dict)

    return (new_record_list, results)


# TODO unittest.
def save_files_and_tkts(record_list, accession_dict, output_folder):
    """Save flat files retrieved from GenBank and create import tickets."""
    import_tickets = []
    genome_folder = pathlib.Path(output_folder, GENOMES_DIR)
    genome_folder.mkdir()
    for record in record_list:
        accession = record.name
        accession = accession.split('.')[0]
        gnm = accession_dict[accession]
        ncbi_filename = f"{gnm.name.lower()}__{accession}.gb"
        flatfile_path = pathlib.Path(genome_folder, ncbi_filename)
        SeqIO.write(record, str(flatfile_path), "genbank")

        tkt = ticket.ImportTicket()
        tkt.type = "replace"
        tkt.phage_id = gnm.id
        tkt.data_dict["host_genus"] = gnm.host_genus
        tkt.data_dict["cluster"] = gnm.cluster
        tkt.data_dict["subcluster"] = gnm.subcluster
        tkt.data_dict["annotation_status"] = gnm.annotation_status
        tkt.data_dict["annotation_author"] = gnm.annotation_author
        tkt.description_field = "product"
        # Accession is set to 'parse' to ensure that during import,
        # the file's accession is directly compared to the database
        # record's accession.
        # tkt.data_dict["accession"] = gnm.accession
        tkt.data_dict["accession"] = "parse"
        tkt.eval_mode = "auto"
        # TODO secondary_phage_id data is for old ticket format.
        tkt.data_dict["secondary_phage_id"] = gnm.id
        tkt.data_dict["retrieve_record"] = 1
        import_tickets.append(tkt)

    # Now make the import table.
    if len(import_tickets) > 0:
        filepath = basic.prepare_filepath(output_folder, "legacy_import_table.csv")
        import_tickets1 = convert_tickets_to_dict(import_tickets, old_format=True)
        basic.export_data_dict(import_tickets1, filepath, IMPORT_COLUMNS1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        if BOTH:
            filepath2 = basic.prepare_filepath(output_folder, "import_table.csv")
            import_tickets2 = convert_tickets_to_dict(import_tickets)
            basic.export_data_dict(import_tickets2, filepath2,
                                   IMPORT_COLUMNS2, include_headers=True)


# TODO unittest.
def create_accession_sets(genome_dict):
    """Generate set of unique and non-unique accessions.

    Input is a dictionary of pdm_utils genome objects."""
    accessions = []
    for id in genome_dict.keys():
        gnm = genome_dict[id]
        # In the MySQL databse, empty accession data is stored as ''.
        # There are no NULL accessions.
        if gnm.accession != "":
            accessions.append(gnm.accession)
    unique, duplicated = basic.identify_unique_items(accessions)
    return unique, duplicated




# TODO unittest.
def create_results_dict(gnm, genbank_date, retrieve_result):
    """Create a dictionary of data summarizing NCBI retrieval status."""
    data_dict = {}
    data_dict["phage_id"] = gnm.id
    data_dict["phage_name"] = gnm.name
    data_dict["accession"] = gnm.accession
    data_dict["annotation_status"] = gnm.annotation_status
    data_dict["mysql_date"] = gnm.date
    data_dict["genbank_date"] = genbank_date
    data_dict["result"] = retrieve_result
    return data_dict


# TODO unittest.
def get_draft_data(output_path, phage_id_set):
    """Run sub-pipeline to retrieve auto-annotated 'draft' genomes."""

    if len(phage_id_set) > 0:
        phage_id_list = list(phage_id_set)
        pecaan_folder = pathlib.Path(output_path, "pecaan")
        pecaan_folder.mkdir()
        retrieve_drafts(pecaan_folder, phage_id_list)
    else:
        print("No new 'draft' genomes available.")



# TODO unittest.
def retrieve_drafts(output_folder, phage_list):
    """Retrieve auto-annotated 'draft' genomes from PECAAN."""

    print(f"\n\nRetrieving {len(phage_list)} new phages from PECAAN")
    genome_folder = pathlib.Path(output_folder, GENOMES_DIR)
    genome_folder.mkdir()

    # Keep track of how many genomes were retrieved from PECAAN
    retrieved_tally = 0
    failed_list = []
    import_tickets = []

    # Iterate through each row in the file
    for new_phage in phage_list:
        pecaan_link = constants.PECAAN_PREFIX + new_phage
        response = phagesdb.retrieve_url_data(pecaan_link)
        if response == "":
            print(f"Error: unable to retrieve {new_phage} draft genome.")
            print(pecaan_link)
            failed_list.append(new_phage)
        else:
            pecaan_filename = f"{new_phage}.txt"
            pecaan_filepath = pathlib.Path(genome_folder, pecaan_filename)
            with pecaan_filepath.open("w") as fh:
                fh.write(response)

            tkt = ticket.ImportTicket()
            tkt.type = "add"
            tkt.phage_id = new_phage
            tkt.data_dict["host_genus"] = "retrieve"
            tkt.data_dict["cluster"] = "retrieve"
            tkt.data_dict["subcluster"] = "retrieve"
            tkt.data_dict["annotation_status"] = "draft"
            tkt.data_dict["annotation_author"] = 1
            tkt.description_field = "product"
            tkt.data_dict["accession"] = "none"
            tkt.eval_mode = "draft"
            # TODO secondary_phage_id data is for old ticket format.
            tkt.data_dict["secondary_phage_id"] = "none"
            tkt.data_dict["retrieve_record"] = 1
            import_tickets.append(tkt)

            print(f"{new_phage} retrieved from PECAAN.")
            retrieved_tally += 1

    # Now make the import table.
    if len(import_tickets) > 0:
        filepath = basic.prepare_filepath(output_folder, "legacy_import_table.csv")
        import_tickets1 = convert_tickets_to_dict(import_tickets, old_format=True)
        basic.export_data_dict(import_tickets1, filepath, IMPORT_COLUMNS1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        if BOTH:
            filepath2 = basic.prepare_filepath(output_folder, "import_table.csv")
            import_tickets2 = convert_tickets_to_dict(import_tickets)
            basic.export_data_dict(import_tickets2, filepath2,
                                   IMPORT_COLUMNS2, include_headers=True)

    # Report results
    if retrieved_tally > 0:
        print(f"{retrieved_tally} phage(s) were successfully retrieved")

    if len(failed_list) > 0:
        print(f"{len(failed_list)} phage(s) failed to be retrieved:")
        for element in failed_list:
            print(element)
        input("\n\nPress ENTER to continue.")
