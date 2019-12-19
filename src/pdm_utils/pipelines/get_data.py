"""Pipeline to gather new data to be imported into a MySQL database."""

import argparse
import csv
from datetime import datetime
import json
import pathlib
import sys
import time
from urllib import request, error
from Bio import SeqIO
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils.classes import genomepair
from pdm_utils.classes import ticket
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import mysqldb
from pdm_utils.functions import phagesdb
from pdm_utils.functions import tickets

# TODO update column headers for import table
# Old format
import_table_columns1 = ["type",
                        "phage_id",
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "annotation_status",
                        "annotation_author",
                        "description_field",
                        "accession",
                        "run_mode",
                        "secondary_phage_id"]

# New format
import_table_columns2 = constants.IMPORT_TABLE_STRUCTURE["order"]

# Columns for new format for update tickets to match RandomFieldUpdateHandler
update_columns = ["table",
                  "field",
                  "value",
                  "key_name",
                  "key_value"]

# Column headers for NCBI results output file.
ncbi_results_header = ["phage_id",
                       "phage_name",
                       "accession",
                       "annotation_status",
                       "mysql_date",
                       "genbank_date",
                       "result"]


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""

    RETRIEVE_HELP = ("Pipeline to retrieve new data to import into a "
                            "MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the directory where updates will be stored.")
    UPDATES_HELP = ("Retrieve updates to HostStrain, Cluster, "
                           "Subcluster, and Accession field data from PhagesDB.")
    DRAFT_HELP = ("Retrieve auto-annotated 'draft' genomes from PECAAN.")
    FINAL_HELP = ("Retrieve new manually-annotated 'final' "
                         "genomes from PhagesDB.")
    GENBANK_HELP = ("Retrieve revised annotated genomes from GenBank.")
    ALL_HELP = ("Retrieve all types of new data.")
    NCBI_CRED_FILE_HELP = ("Path to the file containing NCBI credentials.")

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


    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if args.all_data == True:
        args.updates = True
        args.draft = True
        args.final = True
        args.genbank = True

    if args.genbank == False:
        args.ncbi_credentials_file = None

    return args


# TODO unittest.
def main(unparsed_args_list):
    """Run main retrieve_updates pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    date = time.strftime("%Y%m%d")

    args.output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)

    working_dir = pathlib.Path(f"{date}_new_data")
    working_path = basic.make_new_dir(args.output_folder, working_dir,
                                      attempt=10)

    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    ncbi_cred_dict = ncbi.get_ncbi_creds(args.ncbi_credentials_file)

    # Create data sets
    print("Preparing genome data sets from the MySQL database...")
    sql_handle = mysqldb.connect_to_db(args.database)

    # Parse existing MySQL database genome data to assess what needs to be updated.
    query = ("SELECT PhageID, Name, HostStrain, Status, Cluster2, "
             "DateLastModified, Accession, RetrieveRecord, Subcluster2, "
             "AnnotationAuthor FROM phage")

    mysqldb_genome_list =  mysqldb.parse_genome_data(
                       sql_handle=sql_handle,
                       phage_query=query,
                       gnm_type="mysqldb")
    mysqldb_genome_dict = {}
    for gnm in mysqldb_genome_list:
        mysqldb_genome_dict[gnm.id] = gnm


    # Get data from PhagesDB
    if (args.updates or args.final or args.draft) is True:
        print("Retrieving data from PhagesDB...")
        sequenced_phages_list = phagesdb.get_phagesdb_data(
                                    constants.API_SEQUENCED)
        phagesdb_data_dict = basic.convert_list_to_dict(
                                sequenced_phages_list, "phage_name")
        phagesdb_genome_dict = phagesdb.parse_genomes_dict(
                                    phagesdb_data_dict,
                                    gnm_type="phagesdb",
                                    seq=False)


        # TODO exit if all phage data wasn't retrieved.
        if len(phagesdb_genome_dict) == 0:
            sys.exit(1)

        # Returns a list of tuples.
        match_output = match_genomes(mysqldb_genome_dict, phagesdb_genome_dict)
        matched_genomes = match_output[0]
        unmatched_phagesdb_ids = match_output[1]


    # Option 1: Determine if any fields need to be updated
    if args.updates is True:
        get_update_data(working_path, matched_genomes)

    # Option 2: Determine if any new manually annotated genomes are available.
    if args.final is True:
        get_final_data(working_path, matched_genomes)

    # Option 3: Retrieve updated records from NCBI
    if args.genbank is True:
        get_genbank_data(working_path, mysqldb_genome_dict, ncbi_cred_dict)

    # Option 4: Retrieve auto-annotated genomes from PECAAN
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

    print("\nSummary of genomes matched between the MySQL database and PhagesDB.")
    print(f"{len(matched_ids)} genomes matched.")
    print(f"{len(unmatched_mysqldb_ids)} MySQL genomes not matched.")
    print(f"{len(unmatched_phagesdb_ids)} PhagesDB genomes not matched.")

    unmatched_hatfull_count = len(unmatched_mysqldb_authored_genomes.keys())
    if unmatched_hatfull_count > 0:
        print(f"{unmatched_hatfull_count} Hatfull-authored "
              "unmatched MySQL genomes:")
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


        # Compare Cluster2
        if mysqldb_gnm.cluster != phagesdb_gnm.cluster:
            result1 = {
                      "table":"phage",
                      "field":"Cluster2",
                      "value":phagesdb_gnm.cluster,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result1)
            result2 = {
                      "table":"phage",
                      "field":"Cluster",
                      "value":phagesdb_cluster,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result2)


        # Compare Subcluster2
        if mysqldb_gnm.subcluster != phagesdb_gnm.subcluster:
            result3 = {
                      "table":"phage",
                      "field":"Subcluster2",
                      "value":phagesdb_gnm.subcluster,
                      "key_name":"PhageID",
                      "key_value":mysqldb_gnm.id}
            update_tickets.append(result3)


            if phagesdb_gnm.subcluster != "none":
                result4 = {
                          "table":"phage",
                          "field":"Cluster",
                          "value":phagesdb_gnm.subcluster,
                          "key_name":"PhageID",
                          "key_value":mysqldb_gnm.id}
                update_tickets.append(result4)

        # Compare Host genus
        if mysqldb_gnm.host_genus != phagesdb_gnm.host_genus:
            result5 = {
                      "table":"phage",
                      "field":"HostStrain",
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
        filepath = prepare_output_filepath(output_folder, "update_table.csv",
                                           folder_name="updates")
        tickets.export_ticket_data(update_tickets, filepath, update_columns, include_headers=True)
    else:
        print("\n\nNo field updates found.")


# TODO unittest
def prepare_output_filepath(output_folder, filename, folder_name=None):
    """Prepare output folder."""
    if folder_name is not None:
        new_folder_path = pathlib.Path(output_folder, folder_name)
        new_folder_path.mkdir()
    else:
        new_folder_path = output_folder
    filepath = pathlib.Path(new_folder_path, filename)
    return filepath


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
            if tkt.run_mode == "sea_auto":
                tkt_data["run_mode"] = "ncbi_auto"
            else:
                tkt_data["run_mode"] = tkt.run_mode
        else:
            tkt_data["run_mode"] = tkt.run_mode

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
    phagesdb_genome_folder = pathlib.Path(phagesdb_folder, "genomes")
    phagesdb_genome_folder.mkdir()


    # Initialize PhagesDB retrieval variables - Final updates
    phagesdb_ticket_list = []
    phagesdb_retrieved_tally = 0
    phagesdb_failed_tally = 0
    phagesdb_retrieved_list = []
    phagesdb_failed_list = []

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

        if (phagesdb_gnm.filename == "" or phagesdb_gnm.date < mysqldb_gnm.date):
            phagesdb_failed_tally += 1
            phagesdb_failed_list.append(mysqldb_gnm.id)
        else:
            # Save the file on the hard drive with the same name as
            # stored on PhagesDB
            flatfile_data = phagesdb.retrieve_fasta_data(phagesdb_gnm.filename)
            if flatfile_data == "":
                phagesdb_failed_tally += 1
                phagesdb_failed_list.append(mysqldb_gnm.id)
            else:
                flatfile_filename = phagesdb_gnm.filename.split("/")[-1]
                flatfile_path = pathlib.Path(phagesdb_genome_folder,
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
                tkt = ticket.GenomeTicket()
                tkt.type = "replace"
                tkt.phage_id = mysqldb_gnm.id
                tkt.data_dict["host_genus"] = "retrieve"
                tkt.data_dict["cluster"] = "retrieve"
                tkt.data_dict["subcluster"] = "retrieve"
                tkt.data_dict["annotation_status"] = "final"
                tkt.data_dict["annotation_author"] = 1
                tkt.description_field = "product"
                tkt.data_dict["accession"] = "retrieve"
                tkt.run_mode = "phagesdb"
                # TODO secondary_phage_id data is for old ticket format.
                tkt.data_dict["secondary_phage_id"] = mysqldb_gnm.id
                tkt.data_dict["retrieve_record"] = 1
                phagesdb_ticket_list.append(tkt)
                phagesdb_retrieved_tally += 1
                phagesdb_retrieved_list.append(mysqldb_gnm.id)


    count1 = len(phagesdb_ticket_list)
    if count1 > 0:
        print(f"\n\n{count1} phage(s) were retrieved from PhagesDB.")
        filepath = prepare_output_filepath(phagesdb_folder, "import_table.csv")
        phagesdb_ticket_list = convert_tickets_to_dict(phagesdb_ticket_list, old_format=True)
        tickets.export_ticket_data(phagesdb_ticket_list, filepath, import_table_columns1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        # filepath2 = prepare_output_filepath(phagesdb_folder, "import_table2.csv")
        # phagesdb_ticket_list3 = convert_tickets_to_dict(phagesdb_ticket_list)
        # tickets.export_ticket_data(phagesdb_ticket_list3, filepath2, import_table_columns2, include_headers=True)

    else:
        print("\n\nNo new phages were retrieved from PhagesDB.")






# TODO unittest.
def set_phagesdb_gnm_date(gnm):
    """Set the date of a PhagesDB genome object."""

    # Since there may be multiple 'dates' associated with a genome in PhagesDB,
    # this date attribute can't be directly set when parsing the data.
    # Some phages may have a file but no associated
    # date tagged with that file (since date tagging has only
    # recently been implemented). If there is a date, it is
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
def get_genbank_data(output_folder, mysqldb_genome_dict, ncbi_cred_dict={}):
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
    genome_folder = pathlib.Path(ncbi_folder, "genomes")
    genome_folder.mkdir()

    # Results file
    ncbi_results_list = []

    tally_total = len(mysqldb_genome_dict.keys())
    tally_not_auto_updated = 0
    tally_no_accession = 0
    tally_retrieved_not_new = 0
    tally_retrieved_for_update = 0
    tally_duplicate_accession = 0


    import_ticket_lists = []
    unique_accession_dict = {}

    # Determine if any accessions are duplicated.
    mysqldb_accession_set, mysqldb_duplicate_accessions = get_accessions(mysqldb_genome_dict)

    # Iterate through each phage in the MySQL database
    for key in mysqldb_genome_dict.keys():
        gnm = mysqldb_genome_dict[key]
        # For NCBI retrieval, add to dictionary if
        # 1) the genome is set to be automatically updated and
        # 2) if there is an accession number
        if gnm.retrieve_record != 1:
            tally_not_auto_updated += 1
            ncbi_results = get_ncbi_results_dict(gnm, "NA", "no automatic update")
            ncbi_results_list.append(ncbi_results)

        elif (gnm.accession is None or
                gnm.accession == "none" or
                gnm.accession == ""):
            tally_no_accession += 1
            ncbi_results = get_ncbi_results_dict(gnm, "NA", "no accession")
            ncbi_results_list.append(ncbi_results)

        elif gnm.accession in mysqldb_duplicate_accessions:
            tally_duplicate_accession += 1
            ncbi_results = get_ncbi_results_dict(gnm, "NA", "duplicate accession")
            ncbi_results_list.append(ncbi_results)

        else:
            # Dictionary of phage data based on unique accessions
            # Key = accession
            # Value = genome object
            unique_accession_dict[gnm.accession] = gnm


    batch_size = 200

    # More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    ncbi.set_entrez_credentials(
        tool=ncbi_cred_dict["ncbi_tool"],
        email=ncbi_cred_dict["ncbi_email"],
        api_key=ncbi_cred_dict["ncbi_api_key"])

    print("\n\nRetrieving updated records from NCBI")

    # Use esearch to verify the accessions are valid and efetch to retrieve
    # the record
    # Create batches of accessions
    unique_accession_list = list(unique_accession_dict.keys())

    # Add [ACCN] field to each accession number
    appended_accessions = \
        [accession + "[ACCN]" for accession in unique_accession_list]

    # Keep track of specific records
    retrieved_record_list = []
    retrieval_error_list = []

    # When retrieving in batch sizes, first create the list of values
    # indicating which indices of the unique_accession_list should be used
    # to create each batch.
    # For instace, if there are five accessions, batch size of two produces
    # indices = 0,2,4
    batch_indices = basic.create_indices(unique_accession_list, batch_size)
    print(f"There are {len(unique_accession_list)} GenBank accessions to check.")
    for indices in batch_indices:
        batch_index_start = indices[0]
        batch_index_stop = indices[1]
        print("Checking accessions "
              f"{batch_index_start + 1} to {batch_index_stop}...")
        current_batch_size = batch_index_stop - batch_index_start
        delimiter = " | "
        esearch_term = delimiter.join(appended_accessions[
                                      batch_index_start:batch_index_stop])

        # Use esearch for each accession
        search_record = ncbi.run_esearch(db="nucleotide", term=esearch_term,
                            usehistory="y")
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]


        # Keep track of the accessions that failed to be located in NCBI
        if search_count < current_batch_size:
            search_accession_failure = search_record["ErrorList"][
                "PhraseNotFound"]

            # Each element in this list is formatted "accession[ACCN]"
            for element in search_accession_failure:
                # print(element, element[:-6])
                retrieval_error_list.append(element[:-6])

        # Now get summaries for these records using esummary
        accessions_to_retrieve = []
        summary_records = ncbi.get_summaries(db="nucleotide",
                                             query_key=search_query_key,
                                             webenv=search_webenv)
        for doc_sum in summary_records:
            doc_sum_name = doc_sum["Title"]
            doc_sum_accession = doc_sum["Caption"]
            doc_sum_date = datetime.strptime(doc_sum["UpdateDate"],
                                                      "%Y/%m/%d")
            gnm = unique_accession_dict[doc_sum_accession]

            # If Document Summary date is newer than the MySQL database date,
            # or if the genome being evaluated is currently draft status but
            # we were able to retrieve a Genbank record (i.e. it's not a
            # draft anymore), append the accession to the list of accessions
            # from this batch to retrieve.
            if doc_sum_date > gnm.date or gnm.annotation_status == "draft":
                accessions_to_retrieve.append(doc_sum_accession)
            # Otherwise, if the MySQL database date is newer or the MySQL database
            # status isn't draft, mark down in the ncbi_results file that
            # the record in Genbank isn't new.
            else:
                # We need more information about the MySQL database data for
                # this genome
                tally_retrieved_not_new += 1
                ncbi_results = get_ncbi_results_dict(gnm, doc_sum_date, "record not new")
                ncbi_results_list.append(ncbi_results)

        if len(accessions_to_retrieve) > 0:
            output_list = ncbi.get_records(accessions_to_retrieve,
                                           db="nucleotide",
                                           rettype="gb",
                                           retmode="text")
            retrieved_record_list.extend(output_list)



    # Report the genomes that could not be retrieved.
    tally_retrieval_failure = len(retrieval_error_list)
    for retrieval_error_accession in retrieval_error_list:
        gnm = unique_accession_dict[retrieval_error_accession]
        ncbi_results = get_ncbi_results_dict(gnm, "NA", "retrieval failure")
        ncbi_results_list.append(ncbi_results)

    # For any records that were retrieved, mark their data in the NCBI
    # results file, and save their records as Genbank files, and create
    # import table entries for them.
    for retrieved_record in retrieved_record_list:

        # Pull out the accession to get the matched MySQL database data
        retrieved_record_accession = retrieved_record.name
        retrieved_record_accession = retrieved_record_accession.split('.')[0]

        # Convert date to datetime object
        retrieved_record_date = retrieved_record.annotations["date"]
        retrieved_record_date = datetime.strptime(
            retrieved_record_date, '%d-%b-%Y')

        # MySQL outputs the DateLastModified as a datetime object
        gnm = unique_accession_dict[retrieved_record_accession]

        # Save new records in a folder and create an import table row
        # for them. If the genome is currently a draft annotation, create
        # an import ticket for replacement regardless of the date in
        # the Genbank record. This ensures that if a user fails to
        # upload a manual annotation to phagesdb, once the Genbank
        # accession becomes active MySQL database will get the new version.
        # This should always happen, since we've only retrieved new records.
        if (retrieved_record_date > gnm.date or
                gnm.annotation_status == "draft"):
            tally_retrieved_for_update += 1
            ncbi_results = get_ncbi_results_dict(
                                gnm,
                                retrieved_record_date,
                                "record retrieved for import")
            ncbi_results_list.append(ncbi_results)

            # Determine what the status of the genome will be.
            # If the genome in the MySQL database was already 'final' or 'unknown',
            # then keep the status unchanged. If the genome in the MySQL database
            # was 'draft', the status should be changed to 'final'.
            if gnm.annotation_status == "draft":
                gnm.annotation_status = "final"

            # Now output the file and create the import ticket.
            ncbi_filename = (f"{gnm.name.lower()}__"
                            f"{retrieved_record_accession}.gb")

            flatfile_path = pathlib.Path(genome_folder, ncbi_filename)
            SeqIO.write(retrieved_record, str(flatfile_path), "genbank")

            tkt = ticket.GenomeTicket()
            tkt.type = "replace"
            tkt.phage_id = gnm.id
            tkt.data_dict["host_genus"] = gnm.host_genus
            tkt.data_dict["cluster"] = gnm.cluster
            tkt.data_dict["subcluster"] = gnm.subcluster
            tkt.data_dict["annotation_status"] = gnm.annotation_status
            tkt.data_dict["annotation_author"] = gnm.annotation_author
            tkt.description_field = "product"
            tkt.data_dict["accession"] = gnm.accession
            tkt.run_mode = "sea_auto"
            # TODO secondary_phage_id data is for old ticket format.
            tkt.data_dict["secondary_phage_id"] = gnm.id
            tkt.data_dict["retrieve_record"] = 1
            import_ticket_lists.append(tkt)

        else:
            print("A Genbank record was retrieved "
                  f"for {gnm.id} even though it wasn't new")


    # Now make the import table.
    if len(import_ticket_lists) > 0:
        filepath = prepare_output_filepath(ncbi_folder, "import_table.csv")
        import_ticket_lists = convert_tickets_to_dict(import_ticket_lists, old_format=True)
        tickets.export_ticket_data(import_ticket_lists, filepath, import_table_columns1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        # filepath2 = prepare_output_filepath(ncbi_folder, "import_table2.csv")
        # import_ticket_lists = convert_tickets_to_dict(import_ticket_lists)
        # tickets.export_ticket_data(import_ticket_lists, filepath2, import_table_columns2, include_headers=True)



    # Record all results.
    filepath3 = prepare_output_filepath(ncbi_folder, "ncbi_results.csv")
    tickets.export_ticket_data(ncbi_results_list, filepath3,
                               ncbi_results_header, include_headers=True)




    # Print summary of script
    print(f"Number of genomes in the MySQL database: {tally_total}")
    print("Number of genomes that are NOT set to be updated: "
          f"{tally_not_auto_updated}")
    print("Number of auto-updated genomes with no accession: "
          f"{tally_no_accession}")
    print("Number of auto-updated genomes with a duplicated accession: "
          f"{tally_duplicate_accession}")
    print("Number of records that failed to be retrieved: "
          f"{tally_retrieval_failure}")
    print("Number of records retrieved that are NOT more recent than "
          f"the MySQL record: {tally_retrieved_not_new}")
    print("Number of records retrieved that should be updated in "
          f"the MySQL database: {tally_retrieved_for_update}")
    # input("\n\nPress ENTER to continue.")




# TODO unittest.
def get_accessions(genome_dict):
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
def get_ncbi_results_dict(gnm, genbank_date, retrieve_result):
    """Create a dictionary of data summarizing NCBI retrieval status."""
    result_dict = {}
    result_dict["phage_id"] = gnm.id
    result_dict["phage_name"] = gnm.name
    result_dict["accession"] = gnm.accession
    result_dict["annotation_status"] = gnm.annotation_status
    result_dict["mysql_date"] = gnm.date
    result_dict["genbank_date"] = genbank_date
    result_dict["result"] = retrieve_result
    return result_dict


# TODO unittest.
def get_draft_data(output_path, unmatched_phagesdb_ids):
    """Run sub-pipeline to retrieve auto-annotated 'draft' genomes."""

    if len(unmatched_phagesdb_ids) > 0:
        phagesdb_new_phages_list = list(unmatched_phagesdb_ids)
        pecaan_folder = pathlib.Path(output_path, f"pecaan")
        pecaan_folder.mkdir()
        retrieve_drafts(pecaan_folder, phagesdb_new_phages_list)
    else:
        print("No new 'draft' genomes available.")



# TODO unittest.
def retrieve_drafts(output_folder, phage_list):
    """Retrieve auto-annotated 'draft' genomes from PECAAN."""

    print("\n\nRetrieving new phages from PECAAN")
    pecaan_genome_folder = pathlib.Path(output_folder, "genomes")
    pecaan_genome_folder.mkdir()

    # Keep track of how many genomes were retrieved from PECAAN
    pecaan_retrieved_tally = 0
    pecaan_failed_tally = 0
    pecaan_retrieved_list = []
    pecaan_failed_list = []
    import_ticket_lists = []

    # Iterate through each row in the file
    for new_phage in phage_list:
        pecaan_link = constants.PECAAN_PREFIX + new_phage
        try:
            # gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1) ==> required for
            # urllib2.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED]
            # certificate verify failed (_ssl.c:590)> ==> creating new TLS
            # context tells urllib2 to ignore certificate chain
            # NOTE this is BAD SECURITY, prone to man-in-the-middle attacks
            response = request.urlopen(pecaan_link)
            response_str = str(response.read().decode("utf-8"))
            response.close()

        except:
            print(f"Error: unable to retrieve {new_phage} draft genome.")
            print(pecaan_link)
            pecaan_failed_tally += 1
            pecaan_failed_list.append(new_phage)
            response_str = None

        if response_str is not None:
            pecaan_filename = f"{new_phage}.txt"
            pecaan_filepath = pathlib.Path(pecaan_genome_folder, pecaan_filename)
            with pecaan_filepath.open("w") as fh:
                fh.write(response_str)

            tkt = ticket.GenomeTicket()
            tkt.type = "add"
            tkt.phage_id = new_phage
            tkt.data_dict["host_genus"] = "retrieve"
            tkt.data_dict["cluster"] = "retrieve"
            tkt.data_dict["subcluster"] = "retrieve"
            tkt.data_dict["annotation_status"] = "draft"
            tkt.data_dict["annotation_author"] = 1
            tkt.description_field = "product"
            tkt.data_dict["accession"] = "none"
            tkt.run_mode = "pecaan"
            # TODO secondary_phage_id data is for old ticket format.
            tkt.data_dict["secondary_phage_id"] = "none"
            tkt.data_dict["retrieve_record"] = 1
            import_ticket_lists.append(tkt)

            print(f"{new_phage} retrieved from PECAAN.")
            pecaan_retrieved_tally += 1
            pecaan_retrieved_list.append(new_phage)

    # Now make the import table.
    if len(import_ticket_lists) > 0:
        filepath = prepare_output_filepath(output_folder, "import_table.csv")
        import_ticket_lists = convert_tickets_to_dict(import_ticket_lists, old_format=True)
        tickets.export_ticket_data(import_ticket_lists, filepath, import_table_columns1)

        # TODO new dictwriter. Use this block instead of above once the
        # new import script is functioning.
        # filepath2 = prepare_output_filepath(output_folder, "import_table2.csv")
        # import_ticket_lists = convert_tickets_to_dict(import_ticket_lists)
        # tickets.export_ticket_data(import_ticket_lists, filepath2, import_table_columns2, include_headers=True)

    # Report results
    if pecaan_retrieved_tally > 0:
        print(f"{pecaan_retrieved_tally} phage(s) were successfully retrieved")

    if pecaan_failed_tally > 0:
        print(f"{pecaan_failed_tally} phage(s) failed to be retrieved:")
        for element in pecaan_failed_list:
            print(element)
        input("\n\nPress ENTER to continue.")
