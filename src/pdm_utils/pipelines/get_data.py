"""Pipeline to gather new data to be imported into PhameratorDB."""

import argparse
import csv
import datetime
import json
import pathlib
import sys
import time
from urllib import request, error
from Bio import SeqIO
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import phamerator
from pdm_utils.functions import tickets
from pdm_utils.classes import mysqlconnectionhandler as mch


# TODO Column headers for import table - new version:
# import_table_columns2 = constants.IMPORT_TABLE_STRUCTURE["order"]
 # New order in constants table structure:
 # ["type",
 # "phage_id",
 # "description_field",
 # "run_mode",
 # "host_genus",
 # "cluster",
 # "subcluster",
 # "accession",
 # "annotation_author",
 # "retrieve_record",
 # "annotation_status"
 # ]

# The only diff = old version had secondary_phage_id,
# new version has retrieve_record
import_table_columns2 = ["type",
                        "phage_id",
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "annotation_status",
                        "annotation_author",
                        "description_field",
                        "accession",
                        "run_mode",
                        "retrieve_record"]

# Columns for new format for update tickets to match RandomFieldUpdateHandler
update_columns2 = ["table",
                  "field",
                  "value",
                  "key_name",
                  "key_value"]


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""

    RETRIEVE_HELP = ("Pipeline to retrieve new data to import into a "
                            "MySQL Phamerator database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the directory where updates will be stored.")
    UPDATES_HELP = ("Retrieve updates to HostStrain, Cluster, "
                           "Subcluster, and Accession field data from PhagesDB.")
    DRAFT_HELP = ("Retrieve auto-annotated 'draft' genomes from PECAAN.")
    FINAL_HELP = ("Retrieve new manually-annotated 'final' "
                         "genomes from PhagesDB.")
    GENBANK_HELP = ("Retrieve revised annotated genomes from GenBank.")
    ALL_HELP = ("Retrieve all types of new data.")

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

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if args.all_data == True:
        args.updates = True
        args.draft = True
        args.final = True
        args.genbank = True

    return args


# TODO not tested, but nearly identical function in import_genome.py tested.
def connect_to_db(database):
    """Connect to a MySQL database."""
    sql_handle, msg = phamerator.setup_sql_handle(database)
    if sql_handle is None:
        print(msg)
        sys.exit(1)
    else:
        return sql_handle


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


    # Create data sets
    print("Preparing genome data sets from the phamerator database...")
    sql_handle = connect_to_db(args.database)

    # Parse existing Phamerator genome data to assess what needs to be updated.
    query = (
        "SELECT PhageID, Name, HostStrain, Status, Cluster2, "
        "DateLastModified, Accession, RetrieveRecord, Subcluster2, "
        "AnnotationAuthor FROM phage")

    # Returns a list of elements.
    # Each element is a dictionary from each row in the phage table.
    # Each key is the column name.
    current_genome_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    modified_genome_data_list = modify_phamerator_data(current_genome_list)

    # Get data from PhagesDB
    if (args.updates or args.final or args.draft) is True:

        phagesdb_data_dict = get_phagesdb_data(constants.API_SEQUENCED)
        # Returns a list of tuples.
        match_output = match_genomes(modified_genome_data_list, phagesdb_data_dict)
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
        get_genbank_data(working_path, modified_genome_data_list)

    # Option 4: Retrieve auto-annotated genomes from PECAAN
    if args.draft is True:
        get_draft_data(working_path, unmatched_phagesdb_ids)

    print("\n\n\nRetrieve updates script completed.")



def modify_phamerator_data(input_list):
    """Modify certain types of data retrieved from Phamerator."""
    mod_list = []
    for genome_dict in input_list:
        phamerator_id = genome_dict["PhageID"]
        phamerator_name = genome_dict["Name"]
        phamerator_host = genome_dict["HostStrain"]
        phamerator_status = genome_dict["Status"]
        phamerator_cluster = genome_dict["Cluster2"]
        phamerator_date = genome_dict["DateLastModified"]
        phamerator_accession = genome_dict["Accession"]
        phamerator_retrieve = genome_dict["RetrieveRecord"]
        phamerator_subcluster = genome_dict["Subcluster2"]
        phamerator_author = genome_dict["AnnotationAuthor"]

        # In Phamerator, Singleton Clusters are recorded as '\N', but in
        # phagesdb they are recorded as "Singleton".
        if phamerator_cluster is None:
            phamerator_cluster = 'Singleton'

        # In Phamerator, if Subcluster has not been assigned,
        # Subcluster2 is recorded as '\N'.
        if phamerator_subcluster is None:
            phamerator_subcluster = "none"

        # Accession data may have version number (e.g. XY12345.1)
        if phamerator_accession is None:
            phamerator_accession = "none"

        elif phamerator_accession.strip() == "":
            phamerator_accession = "none"

        if phamerator_accession != "none":
            phamerator_accession = phamerator_accession.split(".")[0]

        # Make sure there is a date in the DateLastModified field
        if phamerator_date is None:
            phamerator_date = datetime.datetime.strptime(
                "1/1/1900", "%m/%d/%Y")

        # Annotation authorship is stored as 1 (Hatfull) or 0 (Genbank/Other)
        if phamerator_author == 1:
            phamerator_author = "hatfull"
        else:
            phamerator_author = "gbk"


        # Output modified genome data
        # 0 = PhageID
        # 1 = PhageName
        # 2 = Host
        # 3 = Status
        # 4 = Cluster2
        # 5 = DateLastModified
        # 6 = Accession
        # 7 = RetrieveRecord
        # 8 = Subcluster2
        # 9 = AnnotationAuthor
        mod_list.append([phamerator_id,
                          phamerator_name,
                          phamerator_host,
                          phamerator_status,
                          phamerator_cluster,
                          phamerator_date,
                          phamerator_accession,
                          phamerator_retrieve,
                          phamerator_subcluster,
                          phamerator_author])
    return mod_list


# TODO unittest.
def match_genomes(modified_genome_data_list, phagesdb_data_dict):
    """Match Phamerator genome data to PhagesDB genome data."""
    unmatched_hatfull_count = 0
    unmatched_phage_id_list = []
    unmatched_hatfull_phage_id_list = []

    # Generate phage_id sets and match sets.
    phagesdb_ids = phagesdb_data_dict.keys()
    phamerator_ids = set()
    for genome_data in modified_genome_data_list:
        phamerator_ids.add(genome_data[0])
    matched_ids = phamerator_ids & phagesdb_ids
    unmatched_phamerator_ids = phamerator_ids - phagesdb_ids
    unmatched_phagesdb_ids = phagesdb_ids - phamerator_ids


    # Match Phamerator data to PhagesDB data.
    # Iterate through each phage in Phamerator
    matched_genomes = []
    for genome_data in modified_genome_data_list:
        phamerator_id = genome_data[0]
        phamerator_author = genome_data[9]
        if phamerator_id in matched_ids:
            matched_phagesdb_data = phagesdb_data_dict[phamerator_id]
            matched_genomes.append((genome_data, matched_phagesdb_data))
        else:
            if phamerator_author == 'hatfull':
                unmatched_hatfull_count += 1
                unmatched_hatfull_phage_id_list.append(phamerator_id)

    print("\nSummary of genomes matched between PhameratorDB and PhagesDB.")
    print(f"{len(matched_ids)} genomes matched.")
    print(f"{len(unmatched_phamerator_ids)} Phamerator genomes not matched.")
    print(f"{len(unmatched_phagesdb_ids)} PhagesDB genomes not matched.")
    if unmatched_hatfull_count > 0:
        print(f"{unmatched_hatfull_count} Hatfull-authored "
              "unmatched PhameratorDB genomes:")
        for element in unmatched_hatfull_phage_id_list:
            print(element)

    return (matched_genomes, unmatched_phagesdb_ids)


# TODO unittest.
def get_phagesdb_data(url):
    """Retrieve all sequence genome data from PhagesDB."""
    print("Retrieving data from phagesdb...")
    sequenced_phages_json = request.urlopen(url)
    # Response is a bytes object that json.loads can't read without first
    # being decoded to a UTF-8 string.
    sequenced_phages_dict = json.loads(sequenced_phages_json.read().decode(
        'utf-8'))
    sequenced_phages_json.close()

    # Data for each phage is stored in a dictionary per phage, and all
    # dictionaries are stored in a list under "results"
    phagesdb_data_dict = {}
    for element_dict in sequenced_phages_dict["results"]:
        phagesdb_data_dict[element_dict["phage_name"]] = element_dict

    # Make sure all PhagesDB data was retrieved.
    exp_count = len(sequenced_phages_dict["results"])
    diff1 = exp_count - sequenced_phages_dict["count"]
    diff2 = exp_count - len(phagesdb_data_dict)
    if (diff1 != 0 or diff2 != 0):
        print("\nUnable to retrieve all phage data from PhagesDB due to "
              "default parameters.")
        sys.exit(1)

    # Standardize genome data for matching.
    for key in phagesdb_data_dict.keys():
        genome_dict = phagesdb_data_dict[key]

        # Accession
        # Sometimes accession data from phagesdb have
        # whitespace characters and version suffix.
        phagesdb_accession = genome_dict["genbank_accession"]
        if phagesdb_accession.strip() != "":
            phagesdb_accession = phagesdb_accession.strip()
            phagesdb_accession = phagesdb_accession.split('.')[0]
        else:
            phagesdb_accession = "none"
        genome_dict["accession_mod"] = phagesdb_accession

        # Cluster
        # Sometimes cluster information is not present. In the
        # phagesdb database, it is is recorded as NULL. When phages
        # data is downloaded from phagesdb, NULL cluster data is
        # converted to "Unclustered". In these cases, leaving the
        # cluster as NULL in phamerator won't work, because NULL
        # means Singleton. Therefore, the phamerator cluster is
        # listed as 'UNK' (Unknown).
        if genome_dict["pcluster"] is None:
            phagesdb_cluster = "UNK"
        else:
            phagesdb_cluster = genome_dict["pcluster"]["cluster"]
        genome_dict["cluster_mod"] = phagesdb_cluster

        # Subcluster
        # If a phage has a cluster, but not a subcluster,
        # set subcluster to 'none'
        if genome_dict["psubcluster"] is None:
            phagesdb_subcluster = "none"
        else:
            phagesdb_subcluster = genome_dict["psubcluster"]["subcluster"]
        genome_dict["subcluster_mod"] = phagesdb_subcluster
        phagesdb_data_dict[key] = genome_dict
    else:
        return phagesdb_data_dict


# TODO unittest
def get_update_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve field updates from PhagesDB."""

    field_update_list = []

    # Iterate through each pair of matched genomes.
    for matched_genome_tuple in matched_genomes:
        genome_data = matched_genome_tuple[0]
        matched_phagesdb_data = matched_genome_tuple[1]

        phamerator_id = genome_data[0]
        phamerator_name = genome_data[1]
        phamerator_host = genome_data[2]
        phamerator_status = genome_data[3]
        phamerator_cluster = genome_data[4]
        phamerator_date = genome_data[5]
        phamerator_accession = genome_data[6]
        phamerator_retrieve = genome_data[7]
        phamerator_subcluster = genome_data[8]
        phamerator_author = genome_data[9]

        # Matched data
        phagesdb_name = matched_phagesdb_data['phage_name']
        phagesdb_host = matched_phagesdb_data['isolation_host']['genus']
        phagesdb_accession = matched_phagesdb_data["accession_mod"]
        phagesdb_cluster = matched_phagesdb_data["cluster_mod"]
        phagesdb_subcluster = matched_phagesdb_data["subcluster_mod"]

        # Compare Cluster2
        if phamerator_cluster != phagesdb_cluster:
            result1 = ["phage",
                       "Cluster2",
                       phagesdb_cluster,
                       "PhageID",
                       phamerator_id]
            field_update_list.append(result1)
            result2 = ["phage",
                       "Cluster",
                       phagesdb_cluster,
                       "PhageID",
                       phamerator_id]
            field_update_list.append(result2)


        # Compare Subcluster2
        if phamerator_subcluster != phagesdb_subcluster:
            result3 = ["phage",
                       "Subcluster2",
                       phagesdb_subcluster,
                       "PhageID",
                       phamerator_id]
            field_update_list.append(result3)


            if phagesdb_subcluster != "none":
                result4 = ["phage",
                           "Cluster",
                           phagesdb_subcluster,
                           "PhageID",
                           phamerator_id]
                field_update_list.append(result4)

        # Compare Host genus
        if phamerator_host != phagesdb_host:
            result5 = ["phage",
                       "HostStrain",
                       phagesdb_host,
                       "PhageID",
                       phamerator_id]
            field_update_list.append(result5)

        # Compare Accession
        # If the genome author is "gbk", then don't worry about
        # updating the accession. This used to be determined with
        # the status field, but now it is determined with the
        # AnnotationAuthor field.
        if phamerator_accession != phagesdb_accession and \
                phamerator_author == "hatfull":
            result6 = ["phage",
                       "Accession",
                       phagesdb_accession,
                       "PhageID",
                       phamerator_id]
            field_update_list.append(result6)


    # Field updates
    if len(field_update_list) > 0:
        print("\n\nNew field updates are available.")
        field_folder = pathlib.Path(output_folder, "updates")
        field_folder.mkdir()
        filename2 = "import_table.csv"
        filepath2 = pathlib.Path(field_folder, filename2)
        with filepath2.open("w") as fh:
            writer = csv.writer(fh)
            writer.writerow(update_columns2)
            writer.writerows(field_update_list)
    else:
        print("\n\nNo field updates found.")



# TODO unittest
def get_final_data(output_folder, matched_genomes):
    """Run sub-pipeline to retrieve 'final' genomes from PhagesDB."""

    phagesdb_folder = pathlib.Path(output_folder, "phagesdb")
    phagesdb_folder.mkdir()
    phagesdb_genome_folder = pathlib.Path(phagesdb_folder, "genomes")
    phagesdb_genome_folder.mkdir()


    # Initialize phagesdb retrieval variables - Final updates
    phagesdb_ticket_list = []
    phagesdb_ticket_list2 = [] # TODO for updated import pipeline
    phagesdb_retrieved_tally = 0
    phagesdb_failed_tally = 0
    phagesdb_retrieved_list = []
    phagesdb_failed_list = []

    # Iterate through each phage in Phamerator
    for matched_genome_tuple in matched_genomes:
        genome_data = matched_genome_tuple[0]
        matched_phagesdb_data = matched_genome_tuple[1]

        phamerator_id = genome_data[0]
        phamerator_name = genome_data[1]
        phamerator_host = genome_data[2]
        phamerator_status = genome_data[3]
        phamerator_cluster = genome_data[4]
        phamerator_date = genome_data[5]
        phamerator_accession = genome_data[6]
        phamerator_retrieve = genome_data[7]
        phamerator_subcluster = genome_data[8]
        phamerator_author = genome_data[9]

        # Matched data
        phagesdb_name = matched_phagesdb_data['phage_name']
        phagesdb_host = matched_phagesdb_data['isolation_host']['genus']
        phagesdb_accession = matched_phagesdb_data["accession_mod"]
        phagesdb_cluster = matched_phagesdb_data["cluster_mod"]
        phagesdb_subcluster = matched_phagesdb_data["subcluster_mod"]

        # Retrieve the qced_genbank_file_date data and properly
        # format it. Some phages may have a file but no associated
        # date tagged with that file (since date tagging has only
        # recently been implemented). If there is a date, it is
        # formatted as: '2017-02-15T10:37:21Z'. If there is no date,
        # it is Null, but change this to 1/1/1900.
        phagesdb_flatfile_date = \
            matched_phagesdb_data["qced_genbank_file_date"]

        if phagesdb_flatfile_date is None:
            phagesdb_flatfile_date = datetime.datetime.strptime(
                "1/1/1900", "%m/%d/%Y")
        else:
            phagesdb_flatfile_date = phagesdb_flatfile_date.split("T")[0]
            phagesdb_flatfile_date = datetime.datetime.strptime(
                phagesdb_flatfile_date, "%Y-%m-%d")

        # Not all phages have associated Genbank-formatted files
        # available on phagesdb. Check to see if there is a flatfile for
        # this phage. Download the flatfile only if there is a date tag,
        # and only if that date is more recent than the date stored in
        # Phamerator for that genome. The tagged date only reflects when
        # the file was uploaded into phagesdb. The date the actual
        # Genbank record was created is stored within the file,
        # and this too could be less recent than the current version in
        # Phamerator; however, this part gets checked during the import
        # stage.
        if (matched_phagesdb_data['qced_genbank_file'] is None or \
                phagesdb_flatfile_date < phamerator_date):
            phagesdb_failed_tally += 1
            phagesdb_failed_list.append(phamerator_id)
        else:
            # Save the file on the hard drive with the same name as
            # stored on phagesdb
            phagesdb_flatfile_url = \
                matched_phagesdb_data["qced_genbank_file"]

            try:
                response5 = request.urlopen(
                    phagesdb_flatfile_url)
                # response comes back as a byte string that won't be
                # processed correctly without being decoded to UTF-8
                response5_str = response5.read().decode("utf-8")
                response5.close()
            except error:  # urllib has 2 main errors, HTTP or URL error
                print("Error: unable to retrieve or read flatfile for "
                      f"phageID {phamerator_id}.")
                phagesdb_failed_tally += 1
                phagesdb_failed_list.append(phamerator_id)
                response5_str = None

            if response5_str is not None:
                phagesdb_filename = phagesdb_flatfile_url.split("/")[-1]
                flatfile_path = pathlib.Path(phagesdb_genome_folder,
                                             phagesdb_filename)
                with flatfile_path.open("w") as fh:
                    fh.write(response5_str)
                # Create the new import ticket
                # Since the phagesdb phage has been matched to
                # the phamerator phage, the AnnotationAuthor field
                # could be assigned from the current phamerator_author
                # variable. However, since this genbank-formatted
                # file is acquired through phagesdb, both the
                # Annotation status is expected to be 'final' and
                # the Annotation author is expected to be 'hatfull'.
                result7 = ["replace",
                           phamerator_id,
                           "retrieve",
                           "retrieve",
                           "retrieve",
                           "final",
                           "hatfull",
                           "product",
                           "retrieve",
                           "phagesdb",
                           phamerator_id]
                phagesdb_ticket_list.append(result7)

                result8 = ["replace",
                           phamerator_id,
                           "retrieve",
                           "retrieve",
                           "retrieve",
                           "final",
                           "1",
                           "product",
                           "retrieve",
                           "phagesdb",
                           "1"]
                phagesdb_ticket_list2.append(result8)
                phagesdb_retrieved_tally += 1
                phagesdb_retrieved_list.append(phamerator_id)


    count1 = len(phagesdb_ticket_list)
    if count1 > 0:
        print(f"\n\n{count1} phage(s) were retrieved from PhagesDB.")
        filename3 = "import_table.csv"
        filepath3 = pathlib.Path(phagesdb_folder, filename3)
        with filepath3.open("w") as fh:
            writer = csv.writer(fh)
            writer.writerows(phagesdb_ticket_list)

    # TODO new dictwriter. Use this block instead of above once the
    # new import script is functioning.
    count2 = len(phagesdb_ticket_list2)
    if count2 > 0:
        print(f"\n\n{count2} phage(s) were retrieved from PhagesDB.")
        filename4 = "dev_import_table.csv"
        filepath4 = pathlib.Path(phagesdb_folder, filename4)
        # with filepath4.open("w") as fh:
        #     writer = csv.writer(fh)
        #     writer.writerow(import_table_columns2)
        #     writer.writerows(phagesdb_ticket_list2)

    else:
        print("\n\nNo new phages were retrieved from PhagesDB.")

    # input("\n\nPress ENTER to continue.")


# TODO unittest.
def get_genbank_data(output_folder, list_of_genomes):
    """Run sub-pipeline to retrieve genomes from GenBank."""

    # Flow of the NCBI record retrieval process:
    # 1 Create list of phages to check for updates at NCBI (completed above)
    # 2 Using esearch, verify which accessions are valid
    # 3 Using esummary, get update date for each valid accession
    # 4 Using efetch, retrieve flat files for NCBI records newer than
    # phamerator date
    # 5 Save new records in a folder and create an import table for them

    # Create output folder
    ncbi_folder = pathlib.Path(output_folder, f"genbank")
    ncbi_folder.mkdir()
    genome_folder = pathlib.Path(ncbi_folder, "genomes")
    genome_folder.mkdir()

    # Results file
    ncbi_results_list = []

    tally_total = len(list_of_genomes)
    tally_not_auto_updated = 0
    tally_no_accession = 0
    tally_retrieved_not_new = 0
    tally_retrieved_for_update = 0
    tally_duplicate_accession = 0


    import_ticket_lists = []
    import_ticket_lists2 = [] # New ticket format

    phamerator_accession_set = set()
    phamerator_duplicate_accessions = []
    unique_accession_dict = {}

    # Determine if any accessions are duplicated.
    for genome_data in list_of_genomes:
        phamerator_accession = genome_data[6]

        # Check for accession duplicates
        if phamerator_accession in phamerator_accession_set:
            phamerator_duplicate_accessions.append(phamerator_accession)
        else:
            phamerator_accession_set.add(phamerator_accession)

    # Iterate through each phage in Phamerator
    for genome_data in list_of_genomes:
        phamerator_id = genome_data[0]
        phamerator_name = genome_data[1]
        phamerator_host = genome_data[2]
        phamerator_status = genome_data[3]
        phamerator_cluster = genome_data[4]
        phamerator_date = genome_data[5]
        phamerator_accession = genome_data[6]
        phamerator_retrieve = genome_data[7]
        phamerator_subcluster = genome_data[8]
        phamerator_author = genome_data[9]

        # For NCBI retrieval, add to dictionary if
        # 1) the genome is set to be automatically updated and
        # 2) if there is an accession number
        if phamerator_retrieve != 1:
            tally_not_auto_updated += 1
            results = [phamerator_id,
                       phamerator_name,
                       phamerator_accession,
                       phamerator_status,
                       phamerator_date,
                       "NA",
                       "no automatic update"]
            ncbi_results_list.append(results)

        elif (phamerator_accession == "none" or phamerator_accession is None):
            tally_no_accession += 1
            results = [phamerator_id,
                       phamerator_name,
                       phamerator_accession,
                       phamerator_status,
                       phamerator_date,
                       "NA",
                       "no accession"]
            ncbi_results_list.append(results)

        elif phamerator_accession in phamerator_duplicate_accessions:
            tally_duplicate_accession += 1
            results = [phamerator_id,
                       phamerator_name,
                       phamerator_accession,
                       phamerator_status,
                       phamerator_date,
                       "NA",
                       "duplicate accession"]
            ncbi_results_list.append(results)

        else:
            # Dictionary of phage data based on unique accessions
            # Key = accession
            # Value = phage data list
            # 0 = PhageID
            # 1 = PhageName
            # 2 = Host
            # 3 = Status
            # 4 = Cluster2
            # 5 = DateLastModified
            # 6 = Accession
            # 7 = RetrieveRecord
            # 8 = Subcluster2
            # 9 = AnnotationAuthor
            unique_accession_dict[phamerator_accession] = [
                phamerator_id, phamerator_name, phamerator_host,
                phamerator_status, phamerator_cluster, phamerator_date,
                phamerator_accession, phamerator_retrieve,
                phamerator_subcluster, phamerator_author]

    batch_size = 200


    # More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    email = input("\nPlease provide email address for NCBI: ")
    ncbi.set_entrez_credentials(
        tool="NCBIRecordRetrievalScript",
        email=email,
        api_key="3b6b113d973599ce1b30c2f94a38508c5908")




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
            doc_sum_date = datetime.datetime.strptime(doc_sum["UpdateDate"],
                                                      "%Y/%m/%d")
            genome_data = unique_accession_dict[doc_sum_accession]
            phamerator_date = genome_data[5]
            phamerator_status = genome_data[3]

            # If Document Summary date is newer than the phamerator date,
            # or if the genome being evaluated is currently draft status but
            # we were able to retrieve a Genbank record (i.e. it's not a
            # draft anymore), append the accession to the list of accessions
            # from this batch to retrieve.
            if doc_sum_date > phamerator_date or phamerator_status == "draft":
                accessions_to_retrieve.append(doc_sum_accession)
            # Otherwise, if the phamerator date is newer or the phamerator
            # status isn't draft, mark down in the ncbi_results file that
            # the record in Genbank isn't new.
            else:
                # We need more information about the phamerator data for
                # this genome

                phamerator_id = genome_data[0]
                phamerator_name = genome_data[1]
                phamerator_host = genome_data[2]
                phamerator_cluster = genome_data[4]
                phamerator_accession = genome_data[6]
                phamerator_retrieve = genome_data[7]
                phamerator_subcluster = genome_data[8]
                phamerator_author = genome_data[9]

                tally_retrieved_not_new += 1
                results = [phamerator_id,
                           phamerator_name,
                           phamerator_accession,
                           phamerator_status,
                           phamerator_date,
                           doc_sum_date,
                           "record not new"]
                ncbi_results_list.append(results)

        if len(accessions_to_retrieve) > 0:
            output_list = ncbi.get_records(accessions_to_retrieve,
                                           db="nucleotide",
                                           rettype="gb",
                                           retmode="text")
            retrieved_record_list.extend(output_list)



    # Report the genomes that could not be retrieved.
    tally_retrieval_failure = len(retrieval_error_list)
    for retrieval_error_accession in retrieval_error_list:
        genome_data = unique_accession_dict[retrieval_error_accession]
        phamerator_id = genome_data[0]
        phamerator_name = genome_data[1]
        phamerator_host = genome_data[2]
        phamerator_status = genome_data[3]
        phamerator_cluster = genome_data[4]
        phamerator_date = genome_data[5]
        phamerator_accession = genome_data[6]
        phamerator_retrieve = genome_data[7]
        phamerator_subcluster = genome_data[8]
        phamerator_author = genome_data[9]
        results = [phamerator_id,
                   phamerator_name,
                   phamerator_accession,
                   phamerator_status,
                   phamerator_date,
                   "NA",
                   "retrieval failure"]
        ncbi_results_list.append(results)

    # For any records that were retrieved, mark their data in the NCBI
    # results file, and save their records as Genbank files, and create
    # import table entries for them.
    for retrieved_record in retrieved_record_list:

        # Pull out the accession to get the matched Phamerator data
        retrieved_record_accession = retrieved_record.name
        retrieved_record_accession = retrieved_record_accession.split('.')[0]

        # Convert date to datetime object
        retrieved_record_date = retrieved_record.annotations["date"]
        retrieved_record_date = datetime.datetime.strptime(
            retrieved_record_date, '%d-%b-%Y')

        # MySQL outputs the DateLastModified as a datetime object
        genome_data = unique_accession_dict[retrieved_record_accession]
        phamerator_id = genome_data[0]
        phamerator_name = genome_data[1]
        phamerator_host = genome_data[2]
        phamerator_status = genome_data[3]
        phamerator_cluster = genome_data[4]
        phamerator_date = genome_data[5]
        phamerator_accession = genome_data[6]
        phamerator_retrieve = genome_data[7]
        phamerator_subcluster = genome_data[8]
        phamerator_author = genome_data[9]

        # Save new records in a folder and create an import table row
        # for them. If the genome is currently a draft annotation, create
        # an import ticket for replacement regardless of the date in
        # the Genbank record. This ensures that if a user fails to
        # upload a manual annotation to phagesdb, once the Genbank
        # accession becomes active Phamerator will get the new version.
        # This should always happen, since we've only retrieved new records.
        if retrieved_record_date > phamerator_date or \
                phamerator_status == 'draft':
            tally_retrieved_for_update += 1
            results = [phamerator_id,
                       phamerator_name,
                       phamerator_accession,
                       phamerator_status,
                       phamerator_date,
                       retrieved_record_date,
                       "record retrieved for import"]
            ncbi_results_list.append(results)

            # Remove the "_Draft" suffix if it is present.
            if phamerator_id[-6:].lower() == '_draft':
                import_table_name = phamerator_id[:-6]
            else:
                import_table_name = phamerator_id

            # Determine what the status of the genome will be.
            # If the genome in Phamerator was already 'final' or 'unknown',
            # then keep the status unchanged. If the genome in Phamerator
            # was 'draft', the status should be changed to 'final'.
            if phamerator_status == 'draft':
                phamerator_status = 'final'

            # Now output the file and create the import ticket.
            ncbi_filename = (f"{phamerator_name.lower()}__"
                            f"{retrieved_record_accession}.gb")

            flatfile_path = pathlib.Path(genome_folder, ncbi_filename)
            SeqIO.write(retrieved_record, str(flatfile_path), "genbank")

            # TODO modify output
            import_ticket_data = ["replace",
                                   import_table_name,
                                   phamerator_host,
                                   phamerator_cluster,
                                   phamerator_subcluster,
                                   phamerator_status,
                                   phamerator_author,
                                   'product',
                                   phamerator_accession,
                                   'ncbi_auto',
                                   phamerator_id]
            import_ticket_lists.append(import_ticket_data)

            # TODO the following code block should replace above once
            # new import pipeline is functioning.
            # Annotation authorship is stored as 1 (Hatfull) or 0 (Genbank/Other)
            if phamerator_author == "hatfull":
                phamerator_author = "1"
            else:
                phamerator_author = "0"
            import_ticket_data2 = ["replace",
                                   import_table_name,
                                   phamerator_host,
                                   phamerator_cluster,
                                   phamerator_subcluster,
                                   phamerator_status,
                                   phamerator_author,
                                   "product",
                                   phamerator_accession,
                                   "sea_auto",
                                   "1"]
            import_ticket_lists2.append(import_ticket_data2)

        else:
            print("For some reason a Genbank record was retrieved "
                  f"for {phamerator_id} even though it wasn't new")


    # Now make the import table.
    if len(import_ticket_lists) > 0:
        filename1 = f"import_table.csv"
        filepath1 = pathlib.Path(ncbi_folder, filename1)
        with filepath1.open("w") as fh:
            writer = csv.writer(fh)
            writer.writerows(import_ticket_lists)


    # TODO new dictwriter. Use this block instead of above once the
    # new import script is functioning.
    if len(import_ticket_lists2) > 0:
        filename2 = "dev_ncbi_updates_import_table.csv"
        filepath2 = pathlib.Path(ncbi_folder, filename2)
        # with filepath2.open("w") as fh:
        #     writer = csv.writer(fh)
        #     writer.writerow(import_table_columns2)
        #     writer.writerows(import_ticket_lists2)



    # Record all results.

    filename3 = "ncbi_results.csv"
    filepath3 = pathlib.Path(ncbi_folder, filename3)
    with filepath3.open("w") as fh:
        writer = csv.writer(fh)
        ncbi_results_header = ["PhageID", "PhageName", "Accession", "Status",
                                "PhameratorDate", "GenBankDate", "Result"]
        writer.writerow(ncbi_results_header)
        writer.writerows(ncbi_results_list)




    # Print summary of script
    print(f"Number of genomes in Phamerator: {tally_total}")
    print("Number of genomes that are NOT set to be updated: "
          f"{tally_not_auto_updated}")
    print("Number of auto-updated genomes with no accession: "
          f"{tally_no_accession}")
    print("Number of auto-updated genomes with a duplicated accession: "
          f"{tally_duplicate_accession}")
    print("Number of records that failed to be retrieved: "
          f"{tally_retrieval_failure}")
    print("Number of records retrieved that are NOT more recent than "
          f"Phamerator record: {tally_retrieved_not_new}")
    print("Number of records retrieved that should be updated in "
          f"Phamerator: {tally_retrieved_for_update}")
    # input("\n\nPress ENTER to continue.")



# TODO unittest.
def get_draft_data(output_path, unmatched_phagesdb_ids):
    """Run sub-pipeline to retrieve auto-annotated 'draft' genomes."""
    # Note: the 'unphamerated_phage_list' generated by PhagesDB only
    # reflects new genomes relative to the most current version of the
    # Actino_Draft database.
    # phagesdb_new_phages_list = \
    #     get_unphamerated_phage_list(constants.UNPHAMERATED_PHAGE_LIST)

    if len(unmatched_phagesdb_ids) > 0:
        phagesdb_new_phages_list = list(unmatched_phagesdb_ids)
        pecaan_folder = pathlib.Path(output_path, f"pecaan")
        pecaan_folder.mkdir()
        retrieve_drafts(pecaan_folder, phagesdb_new_phages_list)
    else:
        print("No new 'draft' genomes available.")

# TODO unittest.
def get_unphamerated_phage_list(url):
    """Retreive list of unphamerated phages from PhagesDB.

    Retrieved file is a tab-delimited text file.
    Each row is a newly-sequenced phage.
    """
    response = request.urlopen(url)
    processed_list = []
    for new_phage in response:
        new_phage = new_phage.strip()  # Remove \t at end of each row
        new_phage = new_phage.decode("utf-8")  # convert bytes object to str
        processed_list.append(new_phage)
    return processed_list




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
    import_ticket_lists2 = [] # New ticket format

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

            # Create the new import ticket
            import_ticket_data = ["add",
                                  new_phage,
                                  "retrieve",
                                  "retrieve",
                                  "retrieve",
                                  "draft",
                                  "hatfull",
                                  "product",
                                  "none",
                                  "pecaan",
                                  "none"]
            import_ticket_lists.append(import_ticket_data)

            # TODO modified output for new import pipeline.
            import_ticket_data2 = ["add",
                                   new_phage,
                                   "retrieve",
                                   "retrieve",
                                   "retrieve",
                                   "draft",
                                   "1",
                                   "product",
                                   "none",
                                   "pecaan",
                                   "1"]
            import_ticket_lists2.append(import_ticket_data2)
            print(f"{new_phage} retrieved from PECAAN.")
            pecaan_retrieved_tally += 1
            pecaan_retrieved_list.append(new_phage)

    # Now make the import table.
    if len(import_ticket_lists) > 0:
        filename1 = "import_table.csv"
        filepath1 = pathlib.Path(output_folder, filename1)
        with filepath1.open("w") as fh:
            writer = csv.writer(fh)
            writer.writerows(import_ticket_lists)

    # TODO new dictwriter. Use this block instead of above once the
    # new import script is functioning.
    if len(import_ticket_lists2) > 0:
        filename2 = "dev_import_table.csv"
        filepath2 = pathlib.Path(output_folder, filename2)
        # with filepath2.open("w") as fh:
        #     writer = csv.writer(fh)
        #     writer.writerow(import_table_columns2)
        #     writer.writerows(import_ticket_lists2)


    # Report results
    if pecaan_retrieved_tally > 0:
        print(f"{pecaan_retrieved_tally} phage(s) were successfully retrieved")

    if pecaan_failed_tally > 0:
        print(f"{pecaan_failed_tally} phage(s) failed to be retrieved:")
        for element in pecaan_failed_list:
            print(element)
        input("\n\nPress ENTER to continue.")
