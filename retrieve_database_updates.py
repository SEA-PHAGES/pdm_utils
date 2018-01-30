#!/usr/bin/env python
#Database update script
#University of Pittsburgh
#Travis Mavrich
#20170212
#The purpose of this script is to provide an interactive environment to
#retrieve several types of Phamerator updates
#that will be imported into Phamerator using the import_phage.py script.
#This script has combined the following independent scripts:
#1. compare_databases.py
#2. retrieve_draft_genomes.py
#3. retrieve_phagesdb_flatfiles.py
#4. retrieve_ncbi_phage.py



#Flow of the import process:
#1 Import python modules, set up variables and functions
#2 Determine which types of retrieval that need to be performed:
    #Option 1: compare Phamerator and phagesdb databases
    #Option 2: retrieve auto-annotated files from PECAAN
    #Option 3: retrieve updated files from NCBI
    #Option 4: retrieve manually-annotated files from phagesdb

#3 Parse Phamerator and phagesdb data to assess what needs to be retrieved
#4 Retrieve field updates
#5 Retrieve new annotations from phagesdb
#6 Retrieve new annotations from NCBI
#7 Retrieve new auto-annotations from PECAAN



#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib, urllib2
from datetime import datetime


#Third-party libraries
try:
    import MySQLdb as mdb
    from Bio import SeqIO, Entrez
except:
    print "\nUnable to import one or more of the following third-party modules: MySQLdb."
    print "Install modules and try again.\n\n"
    sys.exit(1)






#Get the command line parameters
try:
    database = sys.argv[1] #What Phamerator database should be compared to phagesdb?
    working_dir = sys.argv[2] #What is the directory into which the report should go

except:
    print "\n\n\
            This is a python script to retrieve several types of Phamerator database updates.\n\
                1. It retrieves new Host, Cluster, and Subcluster data from phagesdb.\n\
                2. It retrieves manually annotated Genbank-formatted flatfiles from phagesdb.\n\
                3. It retrieves updated Genbank-formatted flatfiles from NCBI.\n\
                4. It retrieves auto-annotated Genbank-formatted flatfiles from PECAAN.\n\n\n\
            It requires two arguments:\n\
            First argument: name of MySQL database that will be checked (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where all reports, update import tables, and retrieved files will be generated.\n\
            All retrieval options create a genomes folder and an import table, if updates are available, with the following format (csv-formatted):\n\
                1. Action to implement on the database (add, remove, replace, update)\n\
                2. PhageID to add or update\n\
                3. Host genus of the updated phage\n\
                4. Cluster of the updated phage\n\
                5. Subcluster of the updated phage\n\
                6. Annotation status of the updated phage (draft, final, gbk)\n\
                7. Annotation authorship of the updated phage (hatfull, gbk)\n\
                8. Gene description field of the update phage (product, note, function)\n\
                9. Accession of the updated phage\n\
                10. PhageID that will be removed or replaced\n\n"

    sys.exit(1)




#Expand home directory
home_dir = os.path.expanduser('~')


#Verify the folder for the consistency report exists

#Add '/' at the end if it's not there
if working_dir[-1] != "/":
    working_dir = working_dir + "/"


#Expand the path if it references the home directory
if working_dir[0] == "~":
    working_dir = home_dir + working_dir[1:]

#Expand the path, to make sure it is a complete directory path
#(in case user inputted path with './path/to/folder')
working_dir = os.path.abspath(working_dir)


if os.path.isdir(working_dir) == False:
    print "\n\nInvalid input for output folder.\n\n"
    sys.exit(1)









#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"




#Set up other variables
date = time.strftime("%Y%m%d")
genomes_folder = "genomes"
new_phage_list_url = 'http://phagesdb.org/data/unphameratedlist'
pecaan_prefix = 'https://discoverdev.kbrinsgd.org/phameratoroutput/phage/'
open_file_handles_list = []




#You have to specify how many results to return at once.
#If you set it to 1 page long and 100,000 genomes/page,
#then this will return everything.
sequenced_phages_url = "http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000"








#Define several functions

#Exits MySQL
def mdb_exit(message):
    print "\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`
    print "\nThe import script did not complete."
    print "\nExiting MySQL."
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print "\nExiting import script."
    print "Closing %s open file handle(s)." %len(open_file_handles_list)
    close_all_files(open_file_handles_list)
    sys.exit(1)


#Allows user to select which retrieval/update options to perform
def select_option(message):
    response = "no"
    response_valid = False
    while response_valid == False:
        response = raw_input(message)
        if (response.lower() == "yes" or response.lower() == "y"):
            response  = "yes"
            response_valid = True
        elif (response.lower() == "no" or response.lower() == "n"):
            response = "no"
            response_valid = True
        else:
            print "Invalid response."
    return response



#Closes all file handles currently open
def close_all_files(file_list):
    for file_handle in file_list:
        file_handle.close()












#Determine which type of updates will be performed.
retrieve_field_updates = select_option("\nDo you want to retrieve Host, Cluster, Subcluster, and Accession updates? (yes or no) ")
retrieve_phagesdb_genomes = select_option("\nDo you want to retrieve manually-annotated genomes from phagesdb? (yes or no) ")
retrieve_pecaan_genomes = select_option("\nDo you want to retrieve auto-annotated genomes from PECAAN? (yes or no) ")
retrieve_ncbi_genomes = select_option("\nDo you want to retrieve updated NCBI records? (yes or no) ")





if retrieve_ncbi_genomes == "yes":

    #Get email infor for NCBI
    contact_email = raw_input("\n\nProvide email for NCBI: ")



    #According to NCBI Entrez documentation, if large numbers of records are requested,
    #it is recommended to split the list into smaller batches of 'around 500 records'.
    #Instead of providing the user the option of setting the batch size, it should be
    #set to a default size of something smaller than 500.
    batch_size = 300




#Create all appropriate output directories
#Different folder for each type of update
#Within each folder, a genomes folder is created to stored only new genome files
#Create all appropriate import table files
if retrieve_field_updates == "yes":


    #Output folder
    field_updates_folder = '%s_field_updates' % date
    field_updates_path = os.path.join(working_dir,field_updates_folder)


    try:
        os.mkdir(field_updates_path)
    except:
        print "\nUnable to create output folder: %s" % field_updates_path
        close_all_files(open_file_handles_list)
        sys.exit(1)

    #This is a dummy folder, since field updates don't have any associated genome files.
    #However, creating the folder makes it easier to run the import_script on
    #the corrections_import_table, since this script relies on
    #the presence of a genomes folder.
    os.mkdir(os.path.join(field_updates_path,genomes_folder))


    #Import table
    field_import_table_handle = open(os.path.join(field_updates_path,date + "_field_update_import_table.csv"), "w")
    open_file_handles_list.append(field_import_table_handle)
    field_import_table_writer = csv.writer(field_import_table_handle)


if retrieve_phagesdb_genomes == "yes":

    #Output folder
    phagesdb_output_folder = '%s_phagesdb_flatfiles' % date
    phagesdb_output_path = os.path.join(working_dir,phagesdb_output_folder)

    try:
        os.mkdir(phagesdb_output_path)
    except:
        print "\nUnable to create output folder: %s" %phagesdb_output_path
        close_all_files(open_file_handles_list)
        sys.exit(1)

    #Genomes folder
    os.mkdir(os.path.join(phagesdb_output_path,genomes_folder))

    #Import table
    phagesdb_import_table_handle = open(os.path.join(phagesdb_output_path,date + "_phagesdb_import_table.csv"), "w")
    open_file_handles_list.append(phagesdb_import_table_handle)
    phagesdb_import_table_writer = csv.writer(phagesdb_import_table_handle)


if retrieve_pecaan_genomes == "yes":

    #Output folder
    pecaan_output_folder = '%s_pecaan_flatfiles' % date
    pecaan_output_path = os.path.join(working_dir,pecaan_output_folder)

    try:
        os.mkdir(pecaan_output_path)
    except:
        print "\nUnable to create output folder: %s" %pecaan_output_path
        close_all_files(open_file_handles_list)
        sys.exit(1)

    #Genomes folder
    os.mkdir(os.path.join(pecaan_output_path,genomes_folder))

    #Import table
    pecaan_import_table_handle = open(os.path.join(pecaan_output_path,date + "_pecaan_import_table.csv"), "w")
    open_file_handles_list.append(pecaan_import_table_handle)
    pecaan_import_table_writer = csv.writer(pecaan_import_table_handle)


if retrieve_ncbi_genomes == "yes":

    #Output folder
    ncbi_output_folder = '%s_ncbi_flatfiles' % date
    ncbi_output_path = os.path.join(working_dir,ncbi_output_folder)

    try:
        os.mkdir(ncbi_output_path)
    except:
        print "\nUnable to create output folder: %s" % ncbi_output_path
        close_all_files(open_file_handles_list)
        sys.exit(1)

    #Genomes folder
    os.mkdir(os.path.join(ncbi_output_path,genomes_folder))

    #Import table
    ncbi_import_table_handle = open(os.path.join(ncbi_output_path,date + "_ncbi_import_table.csv"), "w")
    open_file_handles_list.append(ncbi_import_table_handle)
    ncbi_import_table_writer = csv.writer(ncbi_import_table_handle)


    #Results file
    ncbi_results_handle = open(os.path.join(ncbi_output_path,date + "_ncbi_results.csv"),"w")
    open_file_handles_list.append(ncbi_results_handle)
    ncbi_results_writer = csv.writer(ncbi_results_handle)
    ncbi_results_headers = ['PhageID','PhageName','Accession','Status','PhameratorDate','RetrievedRecordDate','Result']
    ncbi_results_writer.writerow(ncbi_results_headers)






#The following code block parses Phamerator and phagesdb data that
#are required for several types of retrievals.
if (retrieve_field_updates == "yes" or retrieve_phagesdb_genomes == "yes" or retrieve_ncbi_genomes == "yes"):

    print "Retrieving Phamerator data from MySQL database..."
    #Retrieve current data in database
    #0 = PhageID
    #1 = Name
    #2 = HostStrain
    #3 = status
    #4 = Cluster2
    #5 = DateLastModified
    #6 = Accession
    #7 = RetrieveRecord
    #8 = Subcluster2
    #9 = AnnotationAuthor
    try:
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()
    except:
        print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
        close_all_files(open_file_handles_list)
        sys.exit(1)

    try:

        cur.execute("START TRANSACTION")
        cur.execute("SELECT version FROM version")
        db_version = str(cur.fetchone()[0])
        cur.execute("SELECT PhageID,Name,HostStrain,status,Cluster2,DateLastModified,\
                            Accession,RetrieveRecord,Subcluster2,AnnotationAuthor FROM phage")
        current_genome_data_tuples = cur.fetchall()
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)

    except:
        mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")

    con.close()



    #Initialize phagesdb field updates variables
    field_update_tally = 0


    #Initialize phagesdb retrieval variables
    phagesdb_retrieved_tally = 0
    phagesdb_failed_tally = 0
    phagesdb_retrieved_list = []
    phagesdb_failed_list = []


    #Initialize tally variables
    tally_total = len(current_genome_data_tuples)
    tally_not_auto_updated = 0
    tally_no_accession = 0
    tally_retrieved_not_new = 0
    tally_retrieved_for_update = 0


    #Initialize variables to match Phamerator and phagesdb data
    matched_count = 0
    unmatched_count = 0
    unmatched_hatfull_count = 0
    unmatched_phage_id_list = []
    unmatched_hatfull_phage_id_list = []


    #Initialize set variables
    phamerator_id_set = set()
    phamerator_name_set = set()
    phamerator_host_set = set()
    phamerator_status_set = set()
    phamerator_cluster_set = set()
    phamerator_accession_set = set()
    phamerator_subcluster_set = set()

    #Initialize data processing variables
    modified_genome_data_list = []
    phamerator_duplicate_accessions = []
    phamerator_duplicate_phage_names = []
    unique_accession_dict = {}


    #Create data sets
    print "Preparing genome data sets from the phamerator database..."
    for genome_tuple in current_genome_data_tuples:

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



        #In Phamerator, Singleton Clusters are recorded as '\N',
        #but in phagesdb they are recorded as "Singleton".
        if phamerator_cluster is None:
            phamerator_cluster = 'Singleton'

        #In Phamerator, if Subcluster has not been assigned,
        #Subcluster2 is recorded as '\N'.
        if phamerator_subcluster is None:
            phamerator_subcluster = 'none'

        #Accession data may have version number (e.g. XY12345.1)
        if phamerator_accession is None:
            phamerator_accession = 'none'

        elif phamerator_accession.strip() == "":
            phamerator_accession = 'none'

        if phamerator_accession != 'none':
            phamerator_accession = phamerator_accession.split('.')[0]

            #Check for accession duplicates
            if phamerator_accession in phamerator_accession_set:
                phamerator_duplicate_accessions.append(phamerator_accession)
            else:
                phamerator_accession_set.add(phamerator_accession)

        #Check for phage name duplicates
        if phamerator_name in phamerator_name_set:
            phamerator_duplicate_phage_names.append(phamerator_name)
        else:
            phamerator_name_set.add(phamerator_name)

        #Make sure there is a date in the DateLastModified field
        if phamerator_date is None:
            phamerator_date = datetime.strptime('1/1/1900','%m/%d/%Y')


        #Annotation authorship is stored as 1 (Hatfull) or 0 (Genbank/Other)
        if phamerator_author == 1:
            phamerator_author = 'hatfull'
        else:
            phamerator_author = 'gbk'


        phamerator_id_set.add(phamerator_id)
        phamerator_host_set.add(phamerator_host)
        phamerator_cluster_set.add(phamerator_cluster)
        phamerator_subcluster_set.add(phamerator_subcluster)


        #Output modified genome data
        #0 = PhageID
        #1 = PhageName
        #2 = Host
        #3 = Status
        #4 = Cluster2
        #5 = DateLastModified
        #6 = Accession
        #7 = RetrieveRecord
        #8 = Subcluster2
        #9 = AnnotationAuthor
        modified_genome_data_list.append([phamerator_id,\
                                            phamerator_name,\
                                            phamerator_host,\
                                            phamerator_status,\
                                            phamerator_cluster,\
                                            phamerator_date,\
                                            phamerator_accession,\
                                            phamerator_retrieve,\
                                            phamerator_subcluster,\
                                            phamerator_author])



    #phagesdb relies on the phageName, and not the phageID.
    #But Phamerator does not require phageName values to be unique.
    #Check if there are any phageName duplications. If there are,
    #they will not be able to be compared to phagesdb data.
    if len(phamerator_duplicate_phage_names) > 0:
        print "Error: Data is not able to be matched to phagesdb because of the following non-unique phage Names in phamerator:"
        for element in phamerator_duplicate_phage_names:
            print element
        retrieve_field_updates = "no"
        retrieve_phagesdb_genomes = "no"
        raw_input("\n\nPress ENTER to proceed")


    if len(phamerator_duplicate_accessions) > 0:
        print "Error: There are duplicate accessions in Phamerator. Unable to proceed with NCBI record retrieval."
        for accession in phamerator_duplicate_accessions:
            print accession
        retrieve_ncbi_genomes = "no"

        raw_input("\n\nPress ENTER to proceed")


    #Retrieve a list of all sequenced phages listed on phagesdb
    print "Retrieving data from phagesdb..."
    sequenced_phages_json = urllib.urlopen(sequenced_phages_url)
    sequenced_phages_dict = json.loads(sequenced_phages_json.read())
    sequenced_phages_json.close()

    #Data for each phage is stored in a dictionary per phage, and all
    #dictionaries are stored in a list under "results"
    phagesdb_data_dict = {}
    for element_dict in sequenced_phages_dict["results"]:
        phagesdb_data_dict[element_dict["phage_name"]] = element_dict


    if (len(sequenced_phages_dict["results"]) != sequenced_phages_dict["count"] or len(sequenced_phages_dict["results"]) != len(phagesdb_data_dict)):
        print "\nUnable to retrieve all phage data from phagesdb due to default parameters."
        print "Unable to retrieve field updates or manually-annotated flatfiles from phagesdb."
        print "Update parameters in script to enable these functions."
        retrieve_field_updates = "no"
        retrieve_phagesdb_genomes = "no"
        raw_input("\n\nPress ENTER to proceed")




    #Iterate through each phage in Phamerator
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




        #For NCBI retrieval, add to dictionary if
        #1) the genome is set to be automatically updated and
        #2) if there is an accession number
        if retrieve_ncbi_genomes == "yes":

            if phamerator_retrieve != 1:
                #print "PhageID %s is not set to be automatically updated by NCBI record." %phamerator_id
                tally_not_auto_updated += 1
                ncbi_results_writer.writerow([phamerator_id,phamerator_name,phamerator_accession,phamerator_status,phamerator_date,'NA','no automatic update'])

            elif phamerator_accession == "none" or phamerator_accession is None:
                #print "PhageID %s is set to be automatically updated, but it does not have an accession number." %phamerator_id
                tally_no_accession += 1
                ncbi_results_writer.writerow([phamerator_id,phamerator_name,phamerator_accession,phamerator_status,phamerator_date,'NA','no accession'])

            else:

                #Dictionary of phage data based on unique accessions
                #Key = accession
                #Value = phage data list
                #0 = PhageID
                #1 = PhageName
                #2 = Host
                #3 = Status
                #4 = Cluster2
                #5 = DateLastModified
                #6 = Accession
                #7 = RetrieveRecord
                #8 = Subcluster2
                #9 = AnnotationAuthor
                unique_accession_dict[phamerator_accession] = [phamerator_id,\
                                                                phamerator_name,\
                                                                phamerator_host,\
                                                                phamerator_status,\
                                                                phamerator_cluster,\
                                                                phamerator_date,\
                                                                phamerator_accession,\
                                                                phamerator_retrieve,\
                                                                phamerator_subcluster,\
                                                                phamerator_author]


        #The next code block is only applicable if all phage data was successfully retrieved from phagesdb
        #If incomplete data was retrieved from phagesdb, the retrieve_field_updates and
        #retrieve_phagesdb_genomes flags should have been set to "no"
        if (retrieve_field_updates == "yes" or retrieve_phagesdb_genomes == "yes"):

            #Ensure the phageID does not have Draft appended
            if phamerator_id[-6:].lower() == "_draft":
                phage_id_search_name = phamerator_id[:-6]
            else:
                phage_id_search_name = phamerator_id

            #Ensure the phage name does not have Draft appended
            if phamerator_name[-6:].lower() == "_draft":
                phage_name_search_name = phamerator_name[:-6]
            else:
                phage_name_search_name = phamerator_name

            #First try to match up the phageID, and if that doesn't work, try to match up the phageName
            if phage_id_search_name in phagesdb_data_dict.keys():
                matched_phagesdb_data = phagesdb_data_dict[phage_id_search_name]
                matched_count += 1

            elif phage_name_search_name in phagesdb_data_dict.keys():
                matched_phagesdb_data = phagesdb_data_dict[phage_name_search_name]
                matched_count += 1

            else:
                matched_phagesdb_data = ""
                unmatched_count += 1
                unmatched_phage_id_list.append(phamerator_id)

                #Only add Hatfull-author unmatched phages to list
                if phamerator_author == 'hatfull':
                    unmatched_hatfull_count += 1
                    unmatched_hatfull_phage_id_list.append(phamerator_id)

                continue


            #Matched name and host
            phagesdb_name = matched_phagesdb_data['phage_name']
            phagesdb_host = matched_phagesdb_data['isolation_host']['genus']


            #Matched accession
            phagesdb_accession = matched_phagesdb_data['genbank_accession']
            if phagesdb_accession.strip() != "":
                phagesdb_accession = phagesdb_accession.strip() #Sometimes accession data from phagesdb have whitespace characters
                phagesdb_accession = phagesdb_accession.split('.')[0] #Sometimes accession data from phagesdb have version suffix
            else:
                phagesdb_accession = "none"




            #Matched cluster
            if matched_phagesdb_data['pcluster'] is None:
                #Sometimes cluster information is not present. In the phagesdb database, it is is recorded as NULL.
                #When phages data is downloaded from phagesdb, NULL cluster data is converted to "Unclustered".
                #In these cases, leaving the cluster as NULL in phamerator won't work, because NULL means Singleton. Therefore, the phamerator cluster is listed as 'UNK' (Unknown).
                phagesdb_cluster = 'UNK'

            else:
                phagesdb_cluster = matched_phagesdb_data['pcluster']['cluster']

            #Matched subcluster
            if matched_phagesdb_data['psubcluster'] is None:
                #If a phage has a cluster, but not a subcluster, set subcluster to 'none'
                phagesdb_subcluster = 'none'

            else:
                phagesdb_subcluster = matched_phagesdb_data['psubcluster']['subcluster']


        #Determine if any fields need updated
        if retrieve_field_updates == "yes":


            #Compare Cluster2
            if phamerator_cluster != phagesdb_cluster:
                field_corrections_needed += 1

            #Compare Subcluster2
            if phamerator_subcluster != phagesdb_subcluster:
                field_corrections_needed += 1


            #Compare Host genus
            if phamerator_host != phagesdb_host:
                field_corrections_needed += 1


            #Compare Accession
            #If the genome author is "gbk", then don't worry about updating the accession
            #This used to be determined with the status field, but now it is
            #determined with the AnnotationAuthor field.
            if phamerator_accession != phagesdb_accession and phamerator_author == "hatfull":
                field_corrections_needed += 1

            #If errors in the Host, Cluster, or Subcluster information were
            #identified, create an import ticket for the import script to implement.
            if field_corrections_needed > 0:
                field_update_tally += 1

                field_import_table_writer.writerow(["update",\
                                                    phamerator_id,\
                                                    phagesdb_host,\
                                                    phagesdb_cluster,\
                                                    phagesdb_subcluster,\
                                                    phamerator_status,\
                                                    phamerator_author,\
                                                    "none",\
                                                    phagesdb_accession,\
                                                    "none"])





        #Determine if any new Genbank-formatted files are available
        if retrieve_phagesdb_genomes == "yes":


            #Retrieve the qced_genbank_file_date data and properly format it.
            #Some phages may have a file but no associated date tagged with that file (since date tagging has only recently been implemented).
            #If there is a date, it is formatted as: '2017-02-15T10:37:21Z'
            #If there is no date, it is Null, but change this to 1/1/1900.
            phagesdb_flatfile_date = matched_phagesdb_data['qced_genbank_file_date']

            if phagesdb_flatfile_date is None:

                phagesdb_flatfile_date = datetime.strptime('1/1/1900','%m/%d/%Y')

            else:

                phagesdb_flatfile_date = phagesdb_flatfile_date.split('T')[0]
                phagesdb_flatfile_date = datetime.strptime(phagesdb_flatfile_date,'%Y-%m-%d')


            #Not all phages have associated Genbank-formatted files available on phagesdb.
            #Check to see if there is a flatfile for this phage.
            #Download the flatfile only if there is a date tag, and only if that date is more recent than the date stored in Phamerator for that genome.
            #The tagged date only reflects when the file was uploaded into phagesdb. The date the actual Genbank record was created is stored within the file,
            #and this too could be less recent than the current version in Phamerator; however, this part gets checked during the import stage.
            if (matched_phagesdb_data['qced_genbank_file'] is None or not phagesdb_flatfile_date > phamerator_date):

                #print "No flatfile is available that is more recent than current phamerator version for phageID %s." % phamerator_id
                phagesdb_failed_tally += 1
                phagesdb_failed_list.append(phamerator_id)

            else:

                #Save the file on the hard drive with the same name as stored on phagesdb
                phagesdb_flatfile_url = matched_phagesdb_data['qced_genbank_file']
                phagesdb_filename = phagesdb_flatfile_url.split('/')[-1]

                try:
                    phagesdb_flatfile_response = urllib2.urlopen(phagesdb_flatfile_url)
                    phagesdb_file_handle = open(os.path.join(phagesdb_output_path,genomes_folder,phagesdb_filename),'w')
                    phagesdb_file_handle.write(phagesdb_flatfile_response.read())
                    phagesdb_flatfile_response.close()
                    phagesdb_file_handle.close()

                    #Create the new import ticket
                    #Since the phagesdb phage has been matched to the phamerator
                    #phage, the AnnotationAuthor field could be assigned from
                    #the current phamerator_author variable. However, since
                    #this genbank-formatted file is acquired through phagesdb,
                    #both the Annotation status is expected to be 'final'
                    #and the Annotation author is expected to be 'hatfull'.
                    phagesdb_import_table_writer.writerow(["replace",\
                                                            phage_id_search_name,\
                                                            "retrieve",\
                                                            "retrieve",\
                                                            "retrieve",\
                                                            "final",\
                                                            "hatfull",\
                                                            "product",\
                                                            "retrieve",\
                                                            phamerator_id])
                    phagesdb_retrieved_tally += 1
                    phagesdb_retrieved_list.append(phamerator_id)

                except:
                    print "Error: unable to retrieve or read flatfile for phageID %s." %phamerator_id
                    phagesdb_failed_tally += 1
                    phagesdb_failed_list.append(phamerator_id)






    #At this point all genomes in Phamerator have been iterated through and matched to phagesdb data
    #All field updates and manually-annotated flatfiles have been retrieved from phagesdb.
    #Report retrieval results.

    print "\nPhamerator-phagesdb matched phage tally: %s." %matched_count
    print "\nPhamerator-phagesdb unmatched phage tally: %s." %unmatched_count

    #Only print out Hatfull-author unmatched phages.
    if unmatched_hatfull_count > 0:
        print "\nPhamerator-phagesdb Hatfull-author unmatched phages:"
        for element in unmatched_hatfull_phage_id_list:
            print element

    #Field updates
    if field_update_tally > 0:
        print "\n\nNew field updates are available."
    else:
        print "\n\nNo field updates found."

    #New flatfiles
    if retrieve_phagesdb_genomes == "yes":

        if phagesdb_retrieved_tally > 0:
            print "\n\n%s phage(s) were retrieved from phagesdb." %phagesdb_retrieved_tally
            #print "\n\n%s phage(s) failed to be retrieved from phagesdb:" %phagesdb_failed_tally
            #for element in phagesdb_failed_list:
                #print element
        else:
            print "\n\nNo new phages were retrieved from phagesdb."

    raw_input("\n\nPress ENTER to continue.")





#Option 3: Retrieve updated records from NCBI
if retrieve_ncbi_genomes == "yes":


    #Flow of the NCBI record retrieval process:
    #1 Create list of phages to check for updates at NCBI (completed in section above)
    #2 Using esearch, verify the accessions are valid
    #3 Retrieve valid records in batches
    #4 Check which records are newer than the upload date of the current version in phamerator
    #5 Save new records in a folder and create an import table for them

    print "\n\nRetrieving updated records from NCBI"



    #Use esearch to verify the accessions are valid and efetch to retrieve the record
    Entrez.email = contact_email
    Entrez.tool = "NCBIRecordRetrievalScript"



    #Create batches of accessions
    unique_accession_list = unique_accession_dict.keys()

    #Add [ACCN] field to each accession number
    index = 0
    while index < len(unique_accession_list):
        unique_accession_list[index] = unique_accession_list[index] + "[ACCN]"
        index += 1

    #Keep track of specific records
    retrieved_record_list = []
    retrieval_error_list = []



    #When retrieving in batch sizes, first create the list of values indicating which indices of the unique_accession_list should be used to create each batch
    #For instace, if there are five accessions, batch size of two produces indices = 0,2,4
    for batch_index_start in range(0,len(unique_accession_list),batch_size):


        if batch_index_start + batch_size > len(unique_accession_list):
            batch_index_stop = len(unique_accession_list)
        else:
            batch_index_stop = batch_index_start + batch_size

        current_batch_size = batch_index_stop - batch_index_start
        delimiter = " | "
        esearch_term = delimiter.join(unique_accession_list[batch_index_start:batch_index_stop])


        #Use esearch for each accession
        search_handle = Entrez.esearch(db = "nucleotide", term = esearch_term,usehistory="y")
        search_record = Entrez.read(search_handle)
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]



        #Keep track of the accessions that failed to be located in NCBI
        if search_count < current_batch_size:
            search_accession_failure = search_record["ErrorList"]["PhraseNotFound"]

            #Each element in this list is formatted "accession[ACCN]"
            for element in search_accession_failure:
                retrieval_error_list.append(element[:-6])



        #Now retrieve all records using efetch
        fetch_handle = Entrez.efetch(db = "nucleotide", rettype = "gb", retmode = "text", retstart = 0,retmax = search_count, webenv = search_webenv,query_key = search_query_key)
        fetch_records = SeqIO.parse(fetch_handle,"genbank")

        for record in fetch_records:
            retrieved_record_list.append(record)

        search_handle.close()
        fetch_handle.close()


    #Report the genomes that could not be retrieved.
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


        #ncbi_results_headers = 'PhageID','PhageName','Accession','Status','PhameratorDate','RetrievedRecordDate','Result'
        ncbi_results_writer.writerow([phamerator_id,\
                                    phamerator_name,\
                                    phamerator_accession,\
                                    phamerator_status,\
                                    phamerator_date,\
                                    'NA',\
                                    'retrieval failure'])


    #Now that all records have been retrieved, check which records are newer
    #than the upload date of the current version in phamerator.
    #Create the genbank-formatted file only if it is a newer genome.
    #Also create an import table.
    for retrieved_record in retrieved_record_list:

        #Pull out the accession to get the matched Phamerator data
        retrieved_record_accession = retrieved_record.name
        retrieved_record_accession = retrieved_record_accession.split('.')[0]

        #Convert date to datetime object
        retrieved_record_date = retrieved_record.annotations["date"]
        retrieved_record_date = datetime.strptime(retrieved_record_date,'%d-%b-%Y')


        #MySQL outputs the DateLastModified as a datetime object
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


        #Save new records in a folder and create an import table row for them.
        #REVIEW If the genome is currently a draft annotation, create an import
        #ticket for replacement regardless of the date in the Genbank record.
        #This ensures that if a user fails to upload a manual annotation to
        #phagesdb, once the Genbank accession becomes active Phamerator will
        #get the new version.
        if retrieved_record_date > phamerator_date or phamerator_status == 'draft':

            tally_retrieved_for_update += 1
            ncbi_results_writer.writerow([phamerator_id,\
                                        phamerator_name,\
                                        phamerator_accession,\
                                        phamerator_status,\
                                        phamerator_date,\
                                        retrieved_record_date,\
                                        'record retrieved for import'])



            #Remove the "_Draft" suffix if it is present.
            if phamerator_id[-6:].lower() == '_draft':
                import_table_name = phamerator_id[:-6]
            else:
                import_table_name = phamerator_id

            #Determine what the status of the genome will be.
            #If the genome in Phamerator was already 'final' or 'gbk', then keep the status unchanged.
            #If the genome in Phamerator was 'draft', the status should be changed to 'final'.
            if phamerator_status == 'draft':
                phamerator_status = 'final'


            #Now output the file and create the import ticket.
            ncbi_filename = phamerator_name.lower() + "__" + retrieved_record_accession + ".gb"
            SeqIO.write(retrieved_record,os.path.join(ncbi_output_path,genomes_folder,ncbi_filename),"genbank")

            ncbi_import_table_writer.writerow(['replace',\
                                                import_table_name,\
                                                phamerator_host,\
                                                phamerator_cluster,\
                                                phamerator_subcluster,\
                                                phamerator_status,\
                                                phamerator_author,\
                                                'product',\
                                                phamerator_accession,\
                                                phamerator_id])


        else:
            tally_retrieved_not_new += 1
            ncbi_results_writer.writerow([phamerator_id,\
                                            phamerator_name,\
                                            phamerator_accession,\
                                            phamerator_status,\
                                            phamerator_date,\
                                            retrieved_record_date,\
                                            'record not new'])



    #Print summary of script
    print "Number of genomes in Phamerator: %s" %tally_total
    print "Number of genomes that are NOT set to be updated: %s" %tally_not_auto_updated
    print "Number of auto-updated genomes with no accession: %s" %tally_no_accession
    print "Number of records that failed to be retrieved: %s" %tally_retrieval_failure
    print "Number of records retrieved that are NOT more recent than Phamerator record: %s" %tally_retrieved_not_new
    print "Number of records retrieved that should be updated in Phamerator: %s" %tally_retrieved_for_update
    raw_input("\n\nPress ENTER to continue.")






#Option 4: Retrieve auto-annotated genomes from PECAAN
if retrieve_pecaan_genomes == "yes":

    print "\n\nRetrieving new phages from PECAAN"


    #Keep track of how many genomes were retrieved from PECAAN
    pecaan_retrieved_tally = 0
    pecaan_failed_tally = 0
    pecaan_retrieved_list = []
    pecaan_failed_list = []


    #Retrieve list of unphamerated genomes
    #Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
    phagesdb_new_phages_response = urllib2.urlopen(new_phage_list_url)




    #Iterate through each row in the file
    for new_phage in phagesdb_new_phages_response:


        #PECAAN should be able to generate any phage that is listed on phagesdb
        new_phage = new_phage.strip() #Remove \t character at the end of each row
        pecaan_link = pecaan_prefix + new_phage
        pecaan_filename = new_phage + "_Draft.txt"


        try:
            pecaan_response = urllib2.urlopen(pecaan_link)
            pecaan_file_handle = open(os.path.join(pecaan_output_path,genomes_folder,pecaan_filename),'w')
            pecaan_file_handle.write(pecaan_response.read())
            pecaan_response.close()
            pecaan_file_handle.close()


            #Create the new import ticket
            pecaan_import_table_writer.writerow(["add",\
                                                new_phage,\
                                                "retrieve",\
                                                "retrieve",\
                                                "retrieve",\
                                                "draft",\
                                                "hatfull",\
                                                "product",\
                                                "none",\
                                                "none"])
            print "Retrieved %s from PECAAN." %new_phage
            pecaan_retrieved_tally += 1
            pecaan_retrieved_list.append(new_phage)

        except:
            print "Error: unable to retrieve %s draft genome." %new_phage
            pecaan_failed_tally += 1
            pecaan_failed_list.append(new_phage)


    phagesdb_new_phages_response.close()


    #Report results
    if pecaan_retrieved_tally > 0:
        print "%s phage(s) were successfully retrieved" %pecaan_retrieved_tally

    if pecaan_failed_tally > 0:
        print "The following %s phage(s) failed to be retrieved:" %pecaan_failed_tally
        for element in pecaan_failed_list:
            print element

    raw_input("\n\nPress ENTER to continue.")






#Close script.

print "Closing %s open file handle(s)." % len(open_file_handles_list)
close_all_files(open_file_handles_list)
print "\n\n\nRetrieve updates script completed."
