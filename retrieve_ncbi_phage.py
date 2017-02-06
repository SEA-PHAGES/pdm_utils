#!/usr/bin/env python
#Script to retrieve updated phage genomes from NCBI
#University of Pittsburgh
#Travis Mavrich
#20161023




#Flow of the import process
#1 Import python modules, set up variables and functions
#2 Retrieve current database information and create list of phages to check for updates at NCBI
#3 Using esearch, verify the accessions are valid
#4 Retrieve valid records in batches
#5 Check which records are newer than the upload date of the current version in phamerator
#6 Save new records in a folder and create an import table for them



#Third-party libraries
from Bio import SeqIO, Entrez
import MySQLdb as mdb

#Built-in libraries
import time, sys, os, getpass, csv
from datetime import datetime


#Get the command line parameters
try:
    database = sys.argv[1]
    output_path = sys.argv[2]

except:
    print "\n\n\
            This is a python script to determine which phage genomes have been updated in NCBI.\n\
            It requires two argument(s):\n\
            First argument: name of MySQL database that will be used (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where output should be created.\n"

    sys.exit(1)
    




#Expand home directory
home_dir = os.path.expanduser('~')



#Verify the output path exists

#Add '/' at the end if it's not there
if output_path[-1] != "/":
    output_path = output_path + "/"

#Expand the path if it references the home directory
if output_path[0] == "~":
    output_path = home_dir + output_path[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
output_path = os.path.abspath(output_path)


if os.path.isdir(output_path) == False:
    print "\n\nInvalid input for output folder.\n\n"
    sys.exit(1)



#Create output directory and processing file
date = time.strftime("%Y%m%d")
output_folder = '%s_retrieved_files' % date
new_dir = os.path.join(output_path,output_folder)

try:
    os.mkdir(new_dir)
except:
    print "\nUnable to create output folder: %s" % new_dir
    sys.exit(1)

os.chdir(new_dir)












#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"




#Get email infor for NCBI
contact_email = raw_input("Provide email for NCBI: ")
print "\n\n"


batch_size = ""
batch_size_valid = False
while batch_size_valid == False:
    batch_size = raw_input("Record retrieval batch size (must be greater than 0 and recommended is 100-200): ")
    print "\n\n"
    if batch_size.isdigit():
        batch_size = int(batch_size)
        if batch_size > 0:
            batch_size_valid = True
        else:
            print "Invalid choice."
            print "\n\n"

    else:
        print "Invalid choice."
        print "\n\n"




#Define several functions

#Exits MySQL
def mdb_exit(message):
    print "\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`
    print "\nThe record retrieveal script did not complete."
    print "\nExiting MySQL."
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print "\nExiting record retrieval script."
    sys.exit(1)









#2 Retrieve current database information and create list of phages to check for updates at NCBI


#Retrieve database version
#Retrieve current data in database
#0 = PhageID
#1 = Name
#2 = HostStrain
#3 = Cluster
#4 = status
#5 = accession
#6 = date last modified

try:
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
    sys.exit(1)

try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT version FROM version")
    db_version = str(cur.fetchone()[0])
    cur.execute("SELECT PhageID,Name,HostStrain,Cluster,status,Accession,DateLastModified FROM phage")
    current_genome_data_tuples = cur.fetchall()
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)

except:
    mdb_exit("\nUnable to access the database to retrieve genome information.")

con.close()




#Initialize tally variables
tally_total = 0
tally_not_final = 0
tally_no_accession = 0
tally_duplicate_accession = 0
tally_retrieval_failure = 0
tally_retrieved_not_new = 0
tally_retrieved_for_update = 0



tally_total = len(current_genome_data_tuples)



processing_results_file = '%s_processing_results.csv' % date
processing_results_file_handle = open(processing_results_file,"w")
processing_results_file_writer = csv.writer(processing_results_file_handle)
file_headers = ['PhageID','PhageName','Accession','Status','PhameratorDate','RetrievedRecordDate','Note']
processing_results_file_writer.writerow(file_headers)










#Create dictionary of phage data based on unique accessions
#Key = accession
#Value = phage data list
#Create list of phage data with duplicate accession info
unique_accession_dict = {}
duplicate_accession_list = []

#Add to dictionary if status is 'final', and if there is an accession number
for phage_tuple in current_genome_data_tuples:

    phage_list = list(phage_tuple)
    

    #Edit some of the phage data fields    
    #When querying NCBI with Accession numbers, efetch retrieves the most updated version. So you can drop the version number after the decimal (e.g. 'XY99999.1')
    if phage_list[5] != "":
        print phage_list[5]
        phage_list[5] = phage_list[5].split('.')[0]
        print phage_list[5]


    #Singleton Cluster values should be converted from None to 'Singleton'
    if phage_list[3] is None:
        phage_list[3] = "Singleton"
        print "PhageID %s Cluster converted to Singleton." %phage_list[0]


 
    #Make sure there is a date in the DateLastModified field
    print phage_list
    if phage_list[6] is None:
        print phage_list[6]
        phage_list[6] = datetime.strptime('1/1/1900','%m/%d/%Y')
        print phage_list[6]

   
   
   
    #Now determine what to do with the data    
    if phage_list[4] != 'final':
        print "PhageID %s is not 'final' status." %phage_list[0]
        tally_not_final += 1
        processing_results_file_writer.writerow([phage_list[0],phage_list[1],phage_list[5],phage_list[4],phage_list[6],'NA','not final status'])

    elif phage_list[5] == "" or phage_list[5] is None:
        print "PhageID %s does not have accession number." %phage_list[0]
        tally_no_accession += 1
        processing_results_file_writer.writerow([phage_list[0],phage_list[1],phage_list[5],phage_list[4],phage_list[6],'NA','no accession'])


    
    elif phage_list[5] in unique_accession_dict.keys():
        print "PhageID %s accession %s is duplicated in the Phamerator database." %(phage_list[0],phage_list[5])
        duplicate_accession_list.append(phage_list)

    else:
        unique_accession_dict[phage_list[5]] = phage_list



#For values that were not unique, remove all accession numbers from the dictionary
temp_list = []
for element in duplicate_accession_list:
    print element
    if element[5] in unique_accession_dict.keys():
        temp_list.append(unique_accession_dict.pop(element[5]))




#Now add these elements from dictionary to the duplicate data list
for element in temp_list:
    duplicate_accession_list.append(element)
    

#Output the duplicate data
tally_duplicate_accession = len(duplicate_accession_list)    
for data_list in duplicate_accession_list:
    processing_results_file_writer.writerow([data_list[0],data_list[1],data_list[5],data_list[4],data_list[6],'NA','duplicate accession'])

    





#3 & #4 In batches of 100, use esearch to verify the accessions are valid and efetch to retrieve the record


Entrez.email = contact_email
Entrez.tool = "GenbankRecordRetrievalScript"



#Create batches of accessions
unique_accession_list = unique_accession_dict.keys()

print "List of accessions to be retrieved:"
print unique_accession_list

#Add [ACCN] field to each accession number
index = 0
while index < len(unique_accession_list):
    unique_accession_list[index] = unique_accession_list[index] + "[ACCN]"
    index += 1

print unique_accession_list


retrieved_record_list = []
retrieval_error_list = []


print len(unique_accession_list)

print range(0,len(unique_accession_list),batch_size)


#When retrieving in batch sizes, first create the list of values indicating which indices of the unique_accession_list should be used to create each batch
#For instace, if there are five accessions, batch size of two produces indices = 0,2,4
for batch_index_start in range(0,len(unique_accession_list),batch_size):

    
    if batch_index_start + batch_size > len(unique_accession_list):
        batch_index_stop = len(unique_accession_list)
    else:
        batch_index_stop = batch_index_start + batch_size
    
    current_batch_size = batch_index_stop - batch_index_start
    print batch_index_start
    print batch_index_stop
    print current_batch_size
    
    delimiter = " | "
    esearch_term = delimiter.join(unique_accession_list[batch_index_start:batch_index_stop])


    print esearch_term
    
    print "Ready to retrieve"
    

    #Use esearch for each accession
    search_handle = Entrez.esearch(db = "nucleotide", term = esearch_term,usehistory="y")
    #time.sleep(5)
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



#5 Now that all records have been retrieved, check which records are newer than the upload date of the current version in phamerator.
# Create the genbank-formatted file only if it is a newer genome
# Also create an import table
#0 = Action = replace
#1 = Name = current phamerator PhageID
#2 = HostStrain = current phamerator hoststrain
#3 = Cluster = current phamerator cluster
#4 = status = final
#5 = Gene Description Field = product
#6 = Genome to replace = current phamerator PhageID
import_table_file = '%s_retrieved_records_import_table.csv' % date
import_table_file_handle = open(import_table_file,"w")
import_table_file_writer = csv.writer(import_table_file_handle)






tally_retrieval_failure = len(retrieval_error_list)
for retrieval_error_accession in retrieval_error_list:

    phamerator_data = unique_accession_dict[retrieval_error_accession]
    processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],'NA','retrieval failure'])




for retrieved_record in retrieved_record_list:
    retrieved_record_accession = retrieved_record.name

    #Convert date date to datetime object
    retrieved_record_date = retrieved_record.annotations["date"]
    retrieved_record_date = datetime.strptime(retrieved_record_date,'%d-%b-%Y')


    #phamerator_date_obj = datetime.strptime(phamerator_date,'%m/%d/%Y')
    #
    #MySQL outputs the DateLastModified as a datetime object
    phamerator_data = unique_accession_dict[retrieved_record_accession]

    #6 Save new records in a folder and create an import table row for them
    if retrieved_record_date > phamerator_data[6]:

        print 'Retrieved record date %s is more recent than phamerator date %s.' %(retrieved_record_date,phamerator_data[6])
        tally_retrieved_for_update += 1
        processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],retrieved_record_date,'record to be updated'])


        #Now output genbank-formatted file to be uploaded to Phamerator and create the import table action
        SeqIO.write(retrieved_record, phamerator_data[1].lower() + "__" + retrieved_record_accession + ".gb","genbank")
        import_table_data_list = ['replace',phamerator_data[0],phamerator_data[2],phamerator_data[3],'final','product',phamerator_data[0]]
        import_table_file_writer.writerow(import_table_data_list)


    else:
        print 'Phamerator date %s is more recent than retrieved record date %s.' %(phamerator_data[6],retrieved_record_date)
        tally_retrieved_not_new += 1
        processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],retrieved_record_date,'record not new'])        
        
import_table_file_handle.close()
processing_results_file_handle.close()











#Print summary of script
print "Number of genomes in Phamerator: %s" %tally_total
print "Number of genomes that are NOT final: %s" %tally_not_final
print "Number of final genomes with no accession: %s" %tally_no_accession
print "Number of duplicate accessions: %s" %tally_duplicate_accession
print "Number of records that failed to be retrieved: %s" %tally_retrieval_failure
print "Number of records retrieved that are not more recent than Phamerator record: %s" %tally_retrieved_not_new
print "Number of records retrieved that should be updated in Phamerator: %s" %tally_retrieved_for_update


processing_check = tally_total - tally_not_final - tally_no_accession - tally_duplicate_accession - tally_retrieval_failure - tally_retrieved_not_new - tally_retrieved_for_update
if processing_check != 0:
    print "Processing check: %s" %processing_check
    print "Error: the processing of phages was not tracked correctly."
    print "\n\n\n"


#Close script.
print "\n\n\n\nRecord retrieval script completed."




