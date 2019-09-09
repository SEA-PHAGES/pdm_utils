import argparse
import csv
import datetime
import getpass
import json
import os
import shutil
import sys
import time
from urllib import request, error

# Import third-party modules
try:
    import pymysql as pms
    import Bio
    from Bio import Entrez, SeqIO
except ModuleNotFoundError as err:
    print(err)
    sys.exit(1)

# from misc_functions import ask_yes_no, close_files
from pdm_utils.functions.basic import ask_yes_no, close_files

# set up argparse to interact with users at the command line interface.
script_description = "Retrieve Phamerator database updates from phagesdb, " \
                     "PECAAN, and Genbank."

parser = argparse.ArgumentParser(description=script_description)

parser.add_argument("database", metavar="db", type=str, nargs=1,
                    help="name of the Phamerator database to find updates for")
parser.add_argument("working_dir", metavar="wd", type=str, nargs=1,
                    help="path to the directory where updates will be stored")

# Parse command line arguments
args = parser.parse_args()
database = args.database[0]
working_dir = os.path.abspath(args.working_dir[0])
print(database, working_dir)

# Check if the working directory is real. If it isn't, kill the script.
if not os.path.exists(working_dir):
    print("Invalid working directory '{}'".format(working_dir))
    sys.exit(1)

# Get username and password for MySQL
username = getpass.getpass(prompt="MySQL username: ")
password = getpass.getpass(prompt="MySQL password: ")

# Connect to MySQL database
try:
    con = pms.connect("localhost", username, password, database)
except pms.err.Error as err:
    print("Error connecting to MySQL database")
    print("Error {}: {}".format(err.args[0], err.args[1]))
    sys.exit(1)

# Set up other variables to be used later in the script
date = time.strftime("%Y%m%d")
phagesdb_unphamerated_url = "https://phagesdb.org/data/unphameratedlist"
pecaan_prefix = "https://discoverdev.kbrinsgd.org/phameratoroutput/phage/"
open_file_handles_list = []

# API at phagesdb has you specify how many results to return. 1 page at
# length 100000 will return everything.
sequenced_phages_url = "https://phagesdb.org/api/sequenced_phages/?page=1" \
                       "&page_size=100000"

# Get user input with respect to which updates to look for.
retrieve_field_updates = ask_yes_no("\nDo you want to retrieve Host, "
                                    "Cluster, Subcluster, and Accession "
                                    "updates? (yes or no) ")
retrieve_phagesdb_updates = ask_yes_no("\nDo you want to retrieve manually "
                                       "annotated genomes from phagesdb? ("
                                       "yes or no) ")
retrieve_pecaan_updates = ask_yes_no("\nDo you want to retrieve "
                                     "auto-annotated genomes from PECAAN? ("
                                     "yes or no) ")
retrieve_ncbi_updates = ask_yes_no("\nDo you want to retrieve updated NCBI "
                                   "records? (yes or no) ")

# More setup variables if NCBI updates are desired.  NCBI Bookshelf resource
# "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
# Sayers, recommends that a single request not contain more than about 200
# UIDS so we will use that as our batch size, and all Entrez requests must
# include the user's email address and tool name.
if retrieve_ncbi_updates is True:
    batch_size = 200
    Entrez.tool = "NCBIRecordRetrievalScript"
    Entrez.email = input("\nPlease provide email address for NCBI: ")
    Entrez.api_key = "3b6b113d973599ce1b30c2f94a38508c5908"

# Create appropriate output directories.  Each update type gets its own
# directory in working_dir, as well as a sub-directory for genomes and a
# file for the import table.  NCBI updates get an additional file to log the
# status of each genome.
if retrieve_field_updates is True:
    # Create output folder
    field_folder = os.path.join(working_dir, "{}_field_updates".format(date))
    try:
        exists = os.path.exists(field_folder)
        if exists is False:
            os.mkdir(field_folder)
        else:
            print("'{}' already exists.".format(field_folder))
            overwrite = ask_yes_no("Do you want to overwrite the contents of "
                                   "this folder? ")
            if overwrite is True:
                try:
                    shutil.rmtree(field_folder)
                    os.mkdir(field_folder)
                except OSError:
                    print("Failed to overwrite '{}'".format(field_folder))
                    print("Exiting script")
                    sys.exit(1)
    except OSError:
        print("Couldn't create '{}'".format(field_folder))
        print("Exiting script")
        sys.exit(1)
    # Genomes sub-folder
    os.mkdir(os.path.join(field_folder, "genomes"))

    # Now make the import table, and add to list of open file handles
    field_import_table = os.path.join(field_folder, date +
                                      "_field_updates_import_table.csv")
    field_import_table_handle = open(field_import_table, "w")
    open_file_handles_list.append(field_import_table_handle)
    field_import_table_writer = csv.writer(field_import_table_handle)

if retrieve_phagesdb_updates is True:
    # Create output folder
    phagesdb_folder = os.path.join(working_dir,
                                   "{}_phagesdb_updates".format(date))
    try:
        exists = os.path.exists(phagesdb_folder)
        if exists is False:
            os.mkdir(phagesdb_folder)
        else:
            print("'{}' already exists.".format(phagesdb_folder))
            overwrite = ask_yes_no("Do you want to overwrite the contents of "
                                   "this folder? ")
            if overwrite is True:
                try:
                    shutil.rmtree(phagesdb_folder)
                    os.mkdir(phagesdb_folder)
                except OSError:
                    print("Failed to overwrite '{}'".format(phagesdb_folder))
                    print("Exiting script")
                    sys.exit(1)
    except OSError:
        print("Couldn't create '{}'".format(phagesdb_folder))
        print("Exiting script")
        sys.exit(1)
    # Genomes sub-folder
    os.mkdir(os.path.join(phagesdb_folder, "genomes"))

    # Now make the import table, and add to list of open file handles
    phagesdb_import_table = os.path.join(phagesdb_folder, date +
                                         "_phagesdb_updates_import_table.csv")
    phagesdb_import_table_handle = open(phagesdb_import_table, "w")
    open_file_handles_list.append(phagesdb_import_table_handle)
    phagesdb_import_table_writer = csv.writer(phagesdb_import_table_handle)

if retrieve_pecaan_updates is True:
    # Create output folder
    pecaan_folder = os.path.join(working_dir, "{}_pecaan_updates".format(date))
    try:
        exists = os.path.exists(pecaan_folder)
        if exists is False:
            os.mkdir(pecaan_folder)
        else:
            print("'{}' already exists.".format(pecaan_folder))
            overwrite = ask_yes_no("Do you want to overwrite the contents of "
                                   "this folder? ")
            if overwrite is True:
                try:
                    shutil.rmtree(pecaan_folder)
                    os.mkdir(pecaan_folder)
                except OSError:
                    print("Failed to overwrite '{}'".format(pecaan_folder))
                    print("Exiting script")
                    sys.exit(1)
    except OSError:
        print("Couldn't create '{}'".format(pecaan_folder))
        print("Exiting script")
        sys.exit(1)
    # Genomes sub-folder
    os.mkdir(os.path.join(pecaan_folder, "genomes"))

    # Now make the import table, and add to list of open file handles
    pecaan_import_table = os.path.join(pecaan_folder, date +
                                       "_pecaan_updates_import_table.csv")
    pecaan_import_table_handle = open(pecaan_import_table, "w")
    open_file_handles_list.append(pecaan_import_table_handle)
    pecaan_import_table_writer = csv.writer(pecaan_import_table_handle)

if retrieve_ncbi_updates is True:
    # Create output folder
    ncbi_folder = os.path.join(working_dir, "{}_ncbi_updates".format(date))
    try:
        exists = os.path.exists(ncbi_folder)
        if exists is False:
            os.mkdir(ncbi_folder)
        else:
            print("'{}' already exists.".format(ncbi_folder))
            overwrite = ask_yes_no("Do you want to overwrite the contents of "
                                   "this folder? ")
            if overwrite is True:
                try:
                    shutil.rmtree(ncbi_folder)
                    os.mkdir(ncbi_folder)
                except OSError:
                    print("Failed to overwrite '{}'".format(ncbi_folder))
                    print("Exiting script")
                    sys.exit(1)
    except OSError:
        print("Couldn't create '{}'".format(ncbi_folder))
        print("Exiting script")
        sys.exit(1)
    # Genomes sub-folder
    os.mkdir(os.path.join(ncbi_folder, "genomes"))

    # Now make the import table, and add to list of open file handles
    ncbi_import_table = os.path.join(ncbi_folder, date +
                                     "_ncbi_updates_import_table.csv")
    ncbi_import_table_handle = open(ncbi_import_table, "w")
    open_file_handles_list.append(ncbi_import_table_handle)
    ncbi_import_table_writer = csv.writer(ncbi_import_table_handle)

    # Results file
    ncbi_results_table = os.path.join(ncbi_folder, date + "_ncbi_results.csv")
    ncbi_results_handle = open(ncbi_results_table, "w")
    open_file_handles_list.append(ncbi_results_handle)
    ncbi_results_writer = csv.writer(ncbi_results_handle)
    ncbi_results_header = ["PhageID", "PhageName", "Accession", "Status",
                           "PhameratorDate", "GenbankDate", "Result"]
    ncbi_results_writer.writerow(ncbi_results_header)

# Parse existing Phamerator genome data so we can compare them to phagesdb
# or NCBI records.
if (retrieve_field_updates or retrieve_phagesdb_updates or
    retrieve_ncbi_updates) is True:
    try:
        cur = con.cursor()
        cur.execute("SELECT PhageID, Name, HostStrain, status, Cluster2, "
                    "DateLastModified, Accession, RetrieveRecord, Subcluster2, "
                    "AnnotationAuthor FROM phage")
        current_genome_tuples = cur.fetchall()
        cur.close()
    except pms.err.Error as err:
        print("Error {}: {}".format(err.args[0], err.args[1]))
        con.close()
        sys.exit(1)

    # Initialize phagesdb field updates variables
    field_update_tally = 0

    # Initialize phagesdb retrieval variables
    phagesdb_retrieved_tally = 0
    phagesdb_failed_tally = 0
    phagesdb_retrieved_list = []
    phagesdb_failed_list = []

    # Initialize tally variables
    tally_total = len(current_genome_tuples)
    tally_not_auto_updated = 0
    tally_no_accession = 0
    tally_retrieved_not_new = 0
    tally_retrieved_for_update = 0

    # Initialize variables to match Phamerator and phagesdb data
    matched_count = 0
    unmatched_count = 0
    unmatched_hatfull_count = 0
    unmatched_phage_id_list = []
    unmatched_hatfull_phage_id_list = []

    # Initialize set variables
    phamerator_id_set = set()
    phamerator_name_set = set()
    phamerator_host_set = set()
    phamerator_status_set = set()
    phamerator_cluster_set = set()
    phamerator_accession_set = set()
    phamerator_subcluster_set = set()

    # Initialize data processing variables
    modified_genome_data_list = []
    phamerator_duplicate_accessions = []
    phamerator_duplicate_phage_names = []
    unique_accession_dict = {}

    # Create data sets
    print("Preparing genome data sets from the phamerator database...")
    for genome_tuple in current_genome_tuples:

        phamerator_id = genome_tuple[0]
        phamerator_name = genome_tuple[1]
        phamerator_host = genome_tuple[2]
        phamerator_status = genome_tuple[3]
        phamerator_cluster = genome_tuple[4]
        phamerator_date = genome_tuple[5]
        phamerator_accession = genome_tuple[6]
        phamerator_retrieve = genome_tuple[7]
        phamerator_subcluster = genome_tuple[8]
        phamerator_author = genome_tuple[9]

        # In Phamerator, Singleton Clusters are recorded as '\N', but in
        # phagesdb they are recorded as "Singleton".
        if phamerator_cluster is None:
            phamerator_cluster = 'Singleton'

        # In Phamerator, if Subcluster has not been assigned,
        # Subcluster2 is recorded as '\N'.
        if phamerator_subcluster is None:
            phamerator_subcluster = 'none'

        # Accession data may have version number (e.g. XY12345.1)
        if phamerator_accession is None:
            phamerator_accession = 'none'

        elif phamerator_accession.strip() == "":
            phamerator_accession = 'none'

        if phamerator_accession != 'none':
            phamerator_accession = phamerator_accession.split('.')[0]

            # Check for accession duplicates
            if phamerator_accession in phamerator_accession_set:
                phamerator_duplicate_accessions.append(phamerator_accession)
            else:
                phamerator_accession_set.add(phamerator_accession)

        # Check for phage name duplicates
        if phamerator_name in phamerator_name_set:
            phamerator_duplicate_phage_names.append(phamerator_name)
        else:
            phamerator_name_set.add(phamerator_name)

        # Make sure there is a date in the DateLastModified field
        if phamerator_date is None:
            phamerator_date = datetime.datetime.strptime(
                '1/1/1900', '%m/%d/%Y')

        # Annotation authorship is stored as 1 (Hatfull) or 0 (Genbank/Other)
        if phamerator_author == 1:
            phamerator_author = 'hatfull'
        else:
            phamerator_author = 'gbk'

        phamerator_id_set.add(phamerator_id)
        phamerator_host_set.add(phamerator_host)
        phamerator_cluster_set.add(phamerator_cluster)
        phamerator_subcluster_set.add(phamerator_subcluster)

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
        modified_genome_data_list.append([phamerator_id,
                                          phamerator_name,
                                          phamerator_host,
                                          phamerator_status,
                                          phamerator_cluster,
                                          phamerator_date,
                                          phamerator_accession,
                                          phamerator_retrieve,
                                          phamerator_subcluster,
                                          phamerator_author])

    # phagesdb relies on the phageName, and not the phageID.
    # But Phamerator does not require phageName values to be unique.
    # Check if there are any phageName duplications. If there are,
    # they will not be able to be compared to phagesdb data.
    if len(phamerator_duplicate_phage_names) > 0:
        print("Error: Data is not able to be matched to phagesdb because of "
              "the following non-unique phage Names in phamerator:")
        for element in phamerator_duplicate_phage_names:
            print(element)
        retrieve_field_updates = "no"
        retrieve_phagesdb_genomes = "no"
        input("\n\nPress ENTER to proceed")

    if len(phamerator_duplicate_accessions) > 0:
        print("Error: There are duplicate accessions in Phamerator. Unable "
              "to proceed with NCBI record retrieval.")
        for accession in phamerator_duplicate_accessions:
            print(accession)
        retrieve_ncbi_genomes = "no"
        input("\n\nPress ENTER to proceed")

    # Retrieve a list of all sequenced phages listed on phagesdb
    print("Retrieving data from phagesdb...")
    sequenced_phages_json = request.urlopen(sequenced_phages_url)
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

    if len(sequenced_phages_dict["results"]) != sequenced_phages_dict[
        "count"]:
        print("\nUnable to retrieve all phage data from phagesdb due to "
              "default parameters.")
        print("Unable to retrieve field updates or manually-annotated "
              "flatfiles from phagesdb.")
        print("Update parameters in script to enable these functions.")
        retrieve_field_updates = "no"
        retrieve_phagesdb_genomes = "no"
        input("\n\nPress ENTER to proceed")

    elif len(sequenced_phages_dict["results"]) != len(phagesdb_data_dict):
        print("\nUnable to retrieve all phage data from phagesdb due to "
              "default parameters.")
        print("Unable to retrieve field updates or manually-annotated "
              "flatfiles from phagesdb.")
        print("Update parameters in script to enable these functions.")
        retrieve_field_updates = "no"
        retrieve_phagesdb_genomes = "no"
        input("\n\nPress ENTER to proceed")

    # Iterate through each phage in Phamerator
    for genome_data in modified_genome_data_list:
        field_corrections_needed = 0

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
        if retrieve_ncbi_updates is True:

            if phamerator_retrieve != 1:
                # print "PhageID %s is not set to be automatically updated
                # by NCBI record." %phamerator_id
                tally_not_auto_updated += 1
                ncbi_results_writer.writerow([phamerator_id,
                                              phamerator_name,
                                              phamerator_accession,
                                              phamerator_status,
                                              phamerator_date,
                                              'NA',
                                              'no automatic update'])

            elif phamerator_accession == "none" or phamerator_accession is \
                    None:
                # print "PhageID %s is set to be automatically updated,
                # but it does not have an accession number." %phamerator_id
                tally_no_accession += 1
                ncbi_results_writer.writerow([phamerator_id,
                                              phamerator_name,
                                              phamerator_accession,
                                              phamerator_status,
                                              phamerator_date,
                                              'NA',
                                              'no accession'])

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

        # The next code block is only applicable if all phage data was
        # successfully retrieved from phagesdb. If incomplete data was
        # retrieved from phagesdb, the retrieve_field_updates and
        # retrieve_phagesdb_genomes flags should have been set to "no"
        if retrieve_field_updates is True or retrieve_phagesdb_updates is True:

            # Ensure the phageID does not have Draft appended
            if phamerator_id[-6:].lower() == "_draft":
                phage_id_search_name = phamerator_id[:-6]
            else:
                phage_id_search_name = phamerator_id

            # Ensure the phage name does not have Draft appended
            if phamerator_name[-6:].lower() == "_draft":
                phage_name_search_name = phamerator_name[:-6]
            else:
                phage_name_search_name = phamerator_name

            # First try to match up the phageID, and if that doesn't work,
            # try to match up the phageName
            if phage_id_search_name in phagesdb_data_dict.keys():
                matched_phagesdb_data = phagesdb_data_dict[
                    phage_id_search_name]
                matched_count += 1

            elif phage_name_search_name in phagesdb_data_dict.keys():
                matched_phagesdb_data = phagesdb_data_dict[
                    phage_name_search_name]
                matched_count += 1

            else:
                matched_phagesdb_data = ""
                unmatched_count += 1
                unmatched_phage_id_list.append(phamerator_id)

                # Only add Hatfull-author unmatched phages to list
                if phamerator_author == 'hatfull':
                    unmatched_hatfull_count += 1
                    unmatched_hatfull_phage_id_list.append(phamerator_id)

                continue

            # Matched name and host
            phagesdb_name = matched_phagesdb_data['phage_name']
            phagesdb_host = matched_phagesdb_data['isolation_host']['genus']

            # Matched accession
            phagesdb_accession = matched_phagesdb_data['genbank_accession']
            if phagesdb_accession.strip() != "":
                # Sometimes accession data from phagesdb have whitespace characters
                phagesdb_accession = phagesdb_accession.strip()
                # Sometimes accession data from phagesdb have version suffix
                phagesdb_accession = phagesdb_accession.split('.')[0]
            else:
                phagesdb_accession = "none"

            # Matched cluster
            if matched_phagesdb_data['pcluster'] is None:
                # Sometimes cluster information is not present. In the
                # phagesdb database, it is is recorded as NULL. When phages
                # data is downloaded from phagesdb, NULL cluster data is
                # converted to "Unclustered". In these cases, leaving the
                # cluster as NULL in phamerator won't work, because NULL
                # means Singleton. Therefore, the phamerator cluster is
                # listed as 'UNK' (Unknown).
                phagesdb_cluster = 'UNK'

            else:
                phagesdb_cluster = matched_phagesdb_data['pcluster']['cluster']

            # Matched subcluster
            if matched_phagesdb_data['psubcluster'] is None:
                # If a phage has a cluster, but not a subcluster,
                # set subcluster to 'none'
                phagesdb_subcluster = 'none'

            else:
                phagesdb_subcluster = matched_phagesdb_data['psubcluster'][
                    'subcluster']

        # Determine if any fields need updated
        if retrieve_field_updates is True:

            # Compare Cluster2
            if phamerator_cluster != phagesdb_cluster:
                field_corrections_needed += 1

            # Compare Subcluster2
            if phamerator_subcluster != phagesdb_subcluster:
                field_corrections_needed += 1

            # Compare Host genus
            if phamerator_host != phagesdb_host:
                field_corrections_needed += 1

            # Compare Accession
            # If the genome author is "gbk", then don't worry about updating
            # the accession. This used to be determined with the status
            # field, but now it is determined with the AnnotationAuthor field.
            if phamerator_accession != phagesdb_accession and \
                    phamerator_author == "hatfull":
                field_corrections_needed += 1

            # If errors in the Host, Cluster, or Subcluster information were
            # identified, create an import ticket for the import script to
            # implement.
            if field_corrections_needed > 0:
                field_update_tally += 1

                field_import_table_writer.writerow(["update",
                                                    phamerator_id,
                                                    phagesdb_host,
                                                    phagesdb_cluster,
                                                    phagesdb_subcluster,
                                                    phamerator_status,
                                                    phamerator_author,
                                                    "none",
                                                    phagesdb_accession,
                                                    "none",
                                                    "none"])

        # Determine if any new Genbank-formatted files are available
        if retrieve_phagesdb_updates is True:

            # Retrieve the qced_genbank_file_date data and properly format it.
            # Some phages may have a file but no associated date tagged with
            # that file (since date tagging has only recently been
            # implemented). If there is a date, it is formatted as:
            # '2017-02-15T10:37:21Z'. If there is no date, it is Null,
            # but change this to 1/1/1900.
            phagesdb_flatfile_date = matched_phagesdb_data[
                'qced_genbank_file_date']

            if phagesdb_flatfile_date is None:

                phagesdb_flatfile_date = datetime.datetime.strptime(
                    '1/1/1900', '%m/%d/%Y')

            else:

                phagesdb_flatfile_date = phagesdb_flatfile_date.split('T')[0]
                phagesdb_flatfile_date = datetime.datetime.strptime(
                    phagesdb_flatfile_date, '%Y-%m-%d')

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
            if matched_phagesdb_data['qced_genbank_file'] is None or \
                    phagesdb_flatfile_date < phamerator_date:

                # print "No flatfile is available that is more recent than
                # current phamerator version for phageID %s." % phamerator_id
                phagesdb_failed_tally += 1
                phagesdb_failed_list.append(phamerator_id)

            else:

                # Save the file on the hard drive with the same name as
                # stored on phagesdb
                phagesdb_flatfile_url = matched_phagesdb_data[
                    'qced_genbank_file']
                phagesdb_filename = phagesdb_flatfile_url.split('/')[-1]

                try:
                    phagesdb_flatfile_response = request.urlopen(
                        phagesdb_flatfile_url)
                    phagesdb_file_handle = open(os.path.join(phagesdb_folder,
                                                             "genomes",
                                                             phagesdb_filename),
                                                'w')
                    # response comes back as a byte string that won't be
                    # processed correctly without being decoded to UTF-8
                    phagesdb_file_handle.write(
                        phagesdb_flatfile_response.read().decode('utf-8'))
                    phagesdb_flatfile_response.close()
                    phagesdb_file_handle.close()

                    # Create the new import ticket
                    # Since the phagesdb phage has been matched to the phamerator
                    # phage, the AnnotationAuthor field could be assigned from
                    # the current phamerator_author variable. However, since
                    # this genbank-formatted file is acquired through phagesdb,
                    # both the Annotation status is expected to be 'final'
                    # and the Annotation author is expected to be 'hatfull'.
                    phagesdb_import_table_writer.writerow(["replace",
                                                           phage_id_search_name,
                                                           "retrieve",
                                                           "retrieve",
                                                           "retrieve",
                                                           "final",
                                                           "hatfull",
                                                           "product",
                                                           "retrieve",
                                                           "phagesdb",
                                                           phamerator_id])

                    phagesdb_retrieved_tally += 1
                    phagesdb_retrieved_list.append(phamerator_id)

                except error:  # urllib has 2 main errors, HTTP or URL error
                    print("Error: unable to retrieve or read flatfile for "
                          "phageID {}.".format(phamerator_id))
                    phagesdb_failed_tally += 1
                    phagesdb_failed_list.append(phamerator_id)

    # At this point all genomes in Phamerator have been iterated through and
    # matched to phagesdb data
    # All field updates and manually-annotated flatfiles have been retrieved
    # from phagesdb.
    # Report retrieval results.

    print("\nPhamerator-phagesdb matched phage tally: {}.".format(
        matched_count))
    print("\nPhamerator-phagesdb unmatched phage tally: {}.".format(
        unmatched_count))
    # Only print out Hatfull-author unmatched phages.
    if unmatched_hatfull_count > 0:
        print("\nPhamerator-phagesdb Hatfull-author unmatched phages:")
        for element in unmatched_hatfull_phage_id_list:
            print(element)

    # Field updates
    if field_update_tally > 0:
        print("\n\nNew field updates are available.")
    else:
        print("\n\nNo field updates found.")

    # New flatfiles
    if retrieve_phagesdb_updates is True:

        if phagesdb_retrieved_tally > 0:
            print("\n\n{} phage(s) were retrieved from phagesdb.".format(
                phagesdb_retrieved_tally))
        # print "\n\n%s phage(s) failed to be retrieved from phagesdb:"
        # %phagesdb_failed_tally
        # for element in phagesdb_failed_list:
        #     print element
        else:
            print("\n\nNo new phages were retrieved from phagesdb.")

    input("\n\nPress ENTER to continue.")

# Option 3: Retrieve updated records from NCBI
if retrieve_ncbi_updates is True:
    # Flow of the NCBI record retrieval process:
    # 1 Create list of phages to check for updates at NCBI (completed above)
    # 2 Using esearch, verify which accessions are valid
    # 3 Using esummary, get update date for each valid accession
    # 4 Using efetch, retrieve flat files for NCBI records newer than
    # phamerator date
    # 5 Save new records in a folder and create an import table for them

    print("\n\nRetrieving updated records from NCBI")

    # Use esearch to verify the accessions are valid and efetch to retrieve
    # the record
    # Create batches of accessions
    unique_accession_list = list(unique_accession_dict.keys())

    # Add [ACCN] field to each accession number
    appended_accessions = [accession + "[ACCN]" for accession in
                           unique_accession_list]

    # Keep track of specific records
    retrieved_record_list = []
    retrieval_error_list = []

    # When retrieving in batch sizes, first create the list of values
    # indicating which indices of the unique_accession_list should be used
    # to create each batch.
    # For instace, if there are five accessions, batch size of two produces
    # indices = 0,2,4
    for batch_index_start in range(0, len(unique_accession_list), batch_size):

        if batch_index_start + batch_size > len(unique_accession_list):
            batch_index_stop = len(unique_accession_list)
        else:
            batch_index_stop = batch_index_start + batch_size

        current_batch_size = batch_index_stop - batch_index_start
        delimiter = " | "
        esearch_term = delimiter.join(appended_accessions[
                                      batch_index_start:batch_index_stop])

        # Use esearch for each accession
        search_handle = Entrez.esearch(db="nucleotide", term=esearch_term,
                                       usehistory="y")
        search_record = Entrez.read(search_handle)
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]

        search_handle.close()

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
        summary_handle = Entrez.esummary(db="nucleotide",
                                         query_key=search_query_key,
                                         webenv=search_webenv)
        summary_records = Entrez.read(summary_handle)
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
                ncbi_results_writer.writerow([phamerator_id, phamerator_name,
                                              phamerator_accession,
                                              phamerator_status,
                                              phamerator_date,
                                              doc_sum_date,
                                              'record not new'])

        summary_handle.close()

        if len(accessions_to_retrieve) > 0:
            fetch_query = ",".join(accessions_to_retrieve)
            fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_query,
                                         rettype="gb", retmode="text")
            fetch_records = SeqIO.parse(fetch_handle, "genbank")
            for record in fetch_records:
                retrieved_record_list.append(record)
            fetch_handle.close()

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

        # ncbi_results_headers = 'PhageID','PhageName','Accession','Status',
        # 'PhameratorDate','RetrievedRecordDate','Result'
        ncbi_results_writer.writerow([phamerator_id, phamerator_name,
                                      phamerator_accession, phamerator_status,
                                      phamerator_date, 'NA', 'retrieval '
                                                             'failure'])

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

        # Save new records in a folder and create an import table row for them.
        # If the genome is currently a draft annotation, create an import
        # ticket for replacement regardless of the date in the Genbank record.
        # This ensures that if a user fails to upload a manual annotation to
        # phagesdb, once the Genbank accession becomes active Phamerator will
        # get the new version.
        # This should always happen, since we've only retrieved new records.
        if retrieved_record_date > phamerator_date or phamerator_status == 'draft':
            tally_retrieved_for_update += 1
            ncbi_results_writer.writerow([phamerator_id,
                                          phamerator_name,
                                          phamerator_accession,
                                          phamerator_status,
                                          phamerator_date,
                                          retrieved_record_date,
                                          'record retrieved for import'])

            # Remove the "_Draft" suffix if it is present.
            if phamerator_id[-6:].lower() == '_draft':
                import_table_name = phamerator_id[:-6]
            else:
                import_table_name = phamerator_id

            # Determine what the status of the genome will be.
            # If the genome in Phamerator was already 'final' or 'gbk',
            # then keep the status unchanged. If the genome in Phamerator
            # was 'draft', the status should be changed to 'final'.
            if phamerator_status == 'draft':
                phamerator_status = 'final'

            # Now output the file and create the import ticket.
            ncbi_filename = phamerator_name.lower() + "__" + retrieved_record_accession + ".gb"
            SeqIO.write(retrieved_record, os.path.join(
                ncbi_folder, "genomes", ncbi_filename), "genbank")

            ncbi_import_table_writer.writerow(['replace',
                                               import_table_name,
                                               phamerator_host,
                                               phamerator_cluster,
                                               phamerator_subcluster,
                                               phamerator_status,
                                               phamerator_author,
                                               'product',
                                               phamerator_accession,
                                               'ncbi_auto',
                                               phamerator_id])

        else:
            print("For some reason I retrieved a Genbank record for {} even "
                  "though it wasn't new...".format(phamerator_id))

    # Print summary of script
    print("Number of genomes in Phamerator: {}".format(tally_total))
    print("Number of genomes that are NOT set to be updated: {}".format(
        tally_not_auto_updated))
    print("Number of auto-updated genomes with no accession: {}".format(
        tally_no_accession))
    print("Number of records that failed to be retrieved: {}".format(
        tally_retrieval_failure))
    print("Number of records retrieved that are NOT more recent than "
          "Phamerator record: {}".format(tally_retrieved_not_new))
    print("Number of records retrieved that should be updated in " \
          "Phamerator: {}".format(tally_retrieved_for_update))
    input("\n\nPress ENTER to continue.")

# Option 4: Retrieve auto-annotated genomes from PECAAN
if retrieve_pecaan_updates is True:

    print("\n\nRetrieving new phages from PECAAN")

    # Keep track of how many genomes were retrieved from PECAAN
    pecaan_retrieved_tally = 0
    pecaan_failed_tally = 0
    pecaan_retrieved_list = []
    pecaan_failed_list = []

    # Retrieve list of unphamerated genomes
    # Retrieved file should be tab-delimited text file, each row is a newly
    # sequenced phage
    phagesdb_new_phages_response = request.urlopen(phagesdb_unphamerated_url)

    # Iterate through each row in the file
    for new_phage in phagesdb_new_phages_response:

        # PECAAN should be able to generate any phage that is listed on
        # phagesdb
        new_phage = new_phage.strip()  # Remove \t at end of each row
        new_phage = new_phage.decode("utf-8")  # convert bytes object to str
        pecaan_link = pecaan_prefix + new_phage

        pecaan_filename = new_phage + ".txt"

        try:
            # gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1) ==> required for
            # urllib2.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED]
            # certificate verify failed (_ssl.c:590)> ==> creating new TLS
            # context tells urllib2 to ignore certificate chain
            # NOTE this is BAD SECURITY, prone to man-in-the-middle attacks

            pecaan_response = request.urlopen(pecaan_link)

            pecaan_file_handle = open(os.path.join(
                pecaan_folder, "genomes", pecaan_filename), 'w')
            pecaan_file_handle.write(str(pecaan_response.read().decode(
                "utf-8")))

            pecaan_response.close()
            pecaan_file_handle.close()

            # Create the new import ticket
            pecaan_import_table_writer.writerow(["add",
                                                 new_phage,
                                                 "retrieve",
                                                 "retrieve",
                                                 "retrieve",
                                                 "draft",
                                                 "hatfull",
                                                 "product",
                                                 "none",
                                                 "pecaan",
                                                 "none"])

            print("Retrieved {} from PECAAN.".format(new_phage))
            pecaan_retrieved_tally += 1
            pecaan_retrieved_list.append(new_phage)

        except error:
            print("Error: unable to retrieve {} draft genome.".format(
                new_phage))
            pecaan_failed_tally += 1
            pecaan_failed_list.append(new_phage)

    phagesdb_new_phages_response.close()

    # Report results
    if pecaan_retrieved_tally > 0:
        print("{} phage(s) were successfully retrieved".format(
            pecaan_retrieved_tally))

    if pecaan_failed_tally > 0:
        print("The following {} phage(s) failed to be retrieved:".format(
            pecaan_failed_tally))
        for element in pecaan_failed_list:
            print(element)

    input("\n\nPress ENTER to continue.")

# Close script.
print("Closing {} open file handle(s).".format(len(open_file_handles_list)))
close_files(open_file_handles_list)
print("\n\n\nRetrieve updates script completed.")
