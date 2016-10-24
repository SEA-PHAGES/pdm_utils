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
            First argument: name of MySQL database that will be updated (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where output should be created.\n\

    sys.exit(1)
    





#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"




#Get email infor for NCBI
contact_email = input("Provide email for NCBI: ")
batch_size = int(input("Record retrieval batch size: "))






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




#Create output directories
date = time.strftime("%Y%m%d")

output_folder = '%s_retrieved_files' % date
new_dir = os.path.join(output_path,output_folder)

try:
    os.mkdir(new_dir)
except:
    print "\nUnable to create output folder: %s" % new_dir)
    sys.exit(1)


os.chdir(new_dir)



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
    output_file.close()
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





#Create dictionary of phage data based on unique accessions
#Key = accession
#Value = phage data tuple
#Create list of phage data with duplicate accession info
unique_accession_dict = {}
duplicate_accession_list = []

#Add to dictionary if status is 'final', and if there is an accession number
for phage_tuple in current_genome_data_tuples:

    
    if phage_tuple[4] != 'final':
        print "PhageID %s is not 'final' status." %phage_tuple[0]

    elif phage_tuple[5] == "":
        print "PhageID does not have accession number." %phage_tuple[0]
    
    elif phage_tuple[5] in unique_accession_dict.keys():
        print "PhageID %s accession %s is duplicated in the Phamerator database." %(phage_tuple[0],phage_tuple[5])
        duplicate_accession_list.append(phage_tuple)

    else:
        unique_accession_dict[phage_tuple[5]] = phage_tuple

#For values that were not unique, remove all accession numbers from the dictionary
temp_list = []
for element in duplicate_accession_list:
    if element in unique_accession_dict.keys():
        temp_list.append(unique_accession_dict.pop(element))

#Now add these elements from dictionary to the duplicate data list
for element in temp_list:
    duplicate_accession_list.append(element)

#Output the duplicate data
if len(duplicate_accession_list) > 0:

    duplicate_data_file = '%s_duplicate_accession_data.csv' % date
    duplicate_data_file_handle = open(duplicate_data_file,"w")
    duplicate_data_file_writer = csv.writer(duplicate_data_file_handle)
    file_headers = ['PhageID','PhageName','Accession']
    duplicate_data_file_writer.write(file_headers)

    for data_tuple in duplicate_accession_list:
        output_data_list = [data_tuple[0],data_tuple[1],data_tuple[5]]
        duplicate_data_file_writer.write(output_data_list)

    duplicate_data_file_handle.close()








#3 & #4 In batches of 100, use esearch to verify the accessions are valid and efetch to retrieve the record


Entrez.email = contact_email
Entrez.tool = "GenbankRecordRetrievalScript"



#Create batches of accessions
unique_accession_list = unique_accession_dict.keys()


#Add [ACCN] field to each accession number
index = 0
while index < len(unique_accession_list):
    unique_accession_list[index] = unique_accession_list[index] + "[ACCN]"
    index += 1




retrieved_record_list = []
accession_error_list = []
for batch_index_start in range(0,len(unique_accession_list),batch_size):

    
    if batch_index_start + batch_size > len(unique_accession_list):
        batch_index_stop = len(unique_accession_list)
    else:
        batch_index_stop = batch_index_start + batch_size
    
    
    
    delimiter = " | "
    esearch_term = delimiter.join(unique_accession_list[batch_index_start:batch_index_stop]




    #Use esearch for each accession
    search_handle = Entrez.esearch(db = "nucleotide", term = esearch_term,usehistory="y")
    search_record = Entrez.read(search_handle)
    search_count = int(search_record["Count"])
    search_webenv = search_record["WebEnv"]
    search_query_key = search_record["QueryKey"]


    
    #Keep track of the accessions that failed to be located in NCBI

    if search_count < batch_size:
        search_accession_failure = search_record["ErrorList"]["PhraseNotFound"]
        #Each element in this list is formatted "accession[ACCN]"
        for element in search_accession_failure:
            accession_error_list.append(element[:-6])
    
    
    
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



for retrieved_record in retrieved_record_list:
    retrieved_record_accession = retrieved_record.name

    #Convert date date to datetime object
    retrieved_record_date = retrieved_record.annotations["date"]
    retrieved_record_date_obj = datetime.strptime(retrieved_record_date,'%d-%b-%Y')

    phamerator_data = unique_accession_dict[retrieved_record_accession]
    phamerator_date = phamerator_data[6].split(' ')[0]
    phamerator_date_obj = datetime.strptime(phamerator_date,'%m/%d/%Y')


    #6 Save new records in a folder and create an import table row for them
    if retrieved_record_date_obj > phamerator_date_obj:
        SeqIO.write(retrieved_record, phamerator_data[1].lower() + "__" + retrieved_record_accession + ".gb","genbank")

        import_table_data_list = ['replace',phamerator_data[0],phamerator_data[2],phamerator_data[3],'final','product',phamerator_data[0]]
        import_table_file_writer.writerow(import_table_data_list)
        
import_table_file_handle.close()


#Close script.
print "\n\n\n\nRecord retrieval script completed.")  




