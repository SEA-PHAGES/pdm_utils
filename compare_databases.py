#!/usr/bin/env python
#Database comparison script
#University of Pittsburgh
#Travis Mavrich
#20170203
#The purpose of this script is to compare the Phamerator and phagesdb databases for inconsistencies and report what needs to be updated.

#Third-party libraries
import MySQLdb as mdb


#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib




#Get the command line parameters
try:
    database = sys.argv[1] #What Phamerator database should be compared to phagesdb?
    updateFileDir = sys.argv[2] #What is the directory into which the report should go
except:
    print "\n\n\
            This is a python script to compare the Phamerator and phagesdb databases for inconsistencies.\n\
            It requires two arguments:\n\
            First argument: name of MySQL database that will be checked (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where the consistency report should be made (csv-formatted).\n\
                    1. Action to implement on the database (add, remove, replace, update)\n\
                    2. PhageID to add or update\n\
                    3. Host genus of the updated phage\n\
                    4. Cluster of the updated phage\n\
                    5. Field that contains the gene description information (product, note, function)\n\
                    6. PhageID that will be removed or replaced\n\n"
    sys.exit(1)




#Expand home directory
home_dir = os.path.expanduser('~')


#Verify the folder for the consistency report exists

#Add '/' at the end if it's not there
if updateFileDir[-1] != "/":
    updateFileDir = updateFileDir + "/"


#Expand the path if it references the home directory
if updateFileDir[0] == "~":
    updateFileDir = home_dir + updateFileDir[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
updateFileDir = os.path.abspath(updateFileDir)


if os.path.isdir(updateFileDir) == False:
    print "\n\nInvalid input for output folder.\n\n"
    sys.exit(1)







#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"







#Define several functions

#Print out statements to both the terminal and to the output file
#For SQL statements that may be long (>150 characters), don't print entire statements.
def write_out(filename,statement):
    if (statement[:7] == "\nINSERT" or statement[:7] == "\nUPDATE" or statement[:7] == "\nDELETE"): 
        if len(statement) > 150:
            print statement[:150] + "...(statement truncated)"
            filename.write(statement[:150] + "...(statement truncated)")
        else:
            print statement
            filename.write(statement)
    else:
        print statement
        filename.write(statement)


#For questionable data, user is requested to clarify if the data is correct or not
def question(message):
    number = -1
    while number < 0:
        value = raw_input("Is this correct? (yes or no): ")
        if (value.lower() == "yes" or value.lower() == "y"):
            number = 0
        elif (value.lower() == "no" or value.lower() == "n"):                         
            write_out(report_file,message)
            number = 1
        else:
            print "Invalid response."
    #This value will be added to the current error total. If 0, no error was encountered. If 1, an error was encountered.
    return number
            
#Exits MySQL
def mdb_exit(message):
    write_out(report_file,"\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
    write_out(report_file,message)
    write_out(report_file,"\nThe import script did not complete.")
    write_out(report_file,"\nExiting MySQL.")
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    write_out(report_file,"\nExiting import script.")
    report_file.close()
    sys.exit(1)
    
    
    
    
    
    
    
#Create output directories
date = time.strftime("%Y%m%d")

output_folder = '%s_database_comparison' % date


try:
    os.mkdir(os.path.join(updateFileDir,output_folder))
except:
    print "\nUnable to create output folder: %s" % os.path.join(updateFileDir,output_folder)
    sys.exit(1)





#Open file to record update information
report_file = open(os.path.join(updateFileDir,output_folder,date + "_database_comparison.txt"), "w")
write_out(report_file,date + " Database comparison:\n\n\n")


#Open file to create import table with changes that need to be implemented
import_table_file = open(os.path.join(updateFileDir,output_folder,date + "_corrections_import_table.csv"), "w")
import_table_writer = csv.writer(import_table_file)


#Retrieve database version
#Retrieve current data in database
#0 = PhageID
#1 = Name
#2 = HostStrain
#3 = Sequence
#4 = status
#5 = Cluster
try:
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
    report_file.close()
    sys.exit(1)

try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT version FROM version")
    db_version = str(cur.fetchone()[0])
    cur.execute("SELECT PhageID,Name,HostStrain,Sequence,status,Cluster FROM phage")
    current_genome_data_tuples = cur.fetchall()
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)

except:
    mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")

con.close()

write_out(report_file,"\nPhamerator database: " + database)
write_out(report_file,"\nPhamerator database version: " + db_version)
    
    
#Variable to track number of warnings and total_errors encountered
warnings = 0
total_errors = 0

#Create data sets
phageId_set = set()
phageName_set = set()
phageHost_set = set()
phageStatus_set = set()
phageCluster_set = set()
print "Preparing genome data sets from the phamerator database..."
for genome_tuple in current_genome_data_tuples:
    phageId_set.add(genome_tuple[0])
    phageName_set.add(genome_tuple[1])
    phageHost_set.add(genome_tuple[2])
    phageStatus_set.add(genome_tuple[4])
    phageCluster_set.add(genome_tuple[5])


#phagesdb relies on the phageName, and not the phageID. But Phamerator does not require phageName values to be unique.
#Check if there are any phageName duplications. If there are, they will not be able to be compared to phagesdb data.
if len(phageId_set) != len(phageName_set):
    print "There appears to be duplicate phageNames in Phamerator. Data is not able to be matched to phagesdb."
    total_errors += question("\nError: phageNames are not unique")
    
    
    
#Phagesdb API to retrieve genome information
api_prefix = "http://phagesdb.org/api/phages/"
api_suffix = "/?format=json"




#Retrieve a list of all sequenced phages listed on phagesdb
#You have to specify how many results to return at once. If you set it to 1 page long and 100,000 genomes/page, then this will return everything
online_phage_list_json = urllib.urlopen("http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000")
online_phage_list_dict = json.loads(online_phage_list_json.read())

#Data for each phage is stored in a dictionary per phage, and all dictionaries are stored in the "results" list
phagesdb_data_dict = {}
for element_dict in online_phage_list_dict["results"]:
    print element_dict["phage_name"]
    phagesdb_data_dict[element_dict["phage_name"]] = element_dict
    
    
if (len(online_phage_list_dict["results"]) != online_phage_list_dict["count"] or len(online_phage_list_dict["results"]) != len(phagesdb_data_dict)):
    write_out(report_file,"\nError: not all phage data retrieved from phagesdb. Change default parameters in script to proceed.")
    total_errors += 1
    sys.exit()
    
    
    



#Now that all phagesdb data retrieved, match up to Phamerator data

matched_count = 0
unmatched_count = 0
unmatched_phageId_list = []

#Iterate through each phage in Phamerator
for genome_tuple in current_genome_data_tuples:

    corrections_needed = 0

    phameratorId = genome_tuple[0]
    phameratorName = genome_tuple[1]
    phameratorHost = genome_tuple[2]
    phameratorSequence = genome_tuple[3]
    phameratorStatus = genome_tuple[4]
    phameratorCluster= genome_tuple[5]

    #In Phamerator, Singleton Clusters are recorded as '\N', but in phagesdb they are recorded as "Singleton"
    if phameratorCluster is None:
        phameratorCluster = 'Singleton'
        

    matched_phagesdb_data = ""
    #print "PhageID: %s" %genome_tuple[0]
    #print "PhageName: %s" %genome_tuple[1]
    
    #Ensure the phageID does not have Draft appended    
    if phameratorId[-6:].lower() == "_draft":
        phageId_search_name = phameratorId[:-6]
    else:
        phageId_search_name = phameratorId

    #Ensure the phage name does not have Draft appended    
    if phameratorName[-6:].lower() == "_draft":
        phageName_search_name = phameratorName[:-6]
    else:
        phageName_search_name = phameratorName   


    #First try to match up the phageID, and if that doesn't work, try to match up the phageName
    if phageId_search_name in phagesdb_data_dict.keys():
        matched_phagesdb_data = phagesdb_data_dict[phageId_search_name]
        matched_count += 1
  
    elif phageName_search_name in phagesdb_data_dict.keys():
        matched_phagesdb_data = phagesdb_data_dict[phageName_search_name]
        matched_count += 1

    else:
        write_out(report_file,"\nError: unable to find phageID %s or phageName %s from phagesdb." %(genome_tuple[0],genome_tuple[1]))
        matched_phagesdb_data = ""
        unmatched_count += 1
        unmatched_phageId_list.append(phameratorId)
        total_errors += 1
        continue


    #Retrieve specific matched data
    phagesdbName = matched_phagesdb_data['phage_name']
    phagesdbHost = matched_phagesdb_data['isolation_host']['genus']
    #phagesdbSequence = 


    if matched_phagesdb_data['pcluster'] is None:
        #Sometimes cluster information is not present. In the phagesdb database, it is is recorded as NULL.
        #When phages data is downloaded from phagesdb, NULL cluster data is converted to "Unclustered".
        #In these cases, leaving the cluster as NULL in phamerator won't work, because NULL means Singleton. Therefore, the phamerator cluster is listed as 'UKN' (Unknown). 
        phagesdbCluster= 'UKN'

    else: 
        phagesdbCluster= matched_phagesdb_data['pcluster']['cluster']



    if matched_phagesdb_data['psubcluster'] is None:
        #If a phage has a cluster, but not a subcluster, set subcluster to Unspecified
        phagesdbSubcluster = 'Unspecified'
    
    else:
        phagesdbSubcluster = matched_phagesdb_data['psubcluster']['subcluster']


    #If the Host and/or cluster data needs updated in Phamerator, decide what the value will be to update the Cluster data.
    if phagesdbSubcluster == 'Unspecified':
        phagesdbClusterUpdate = phagesdbCluster
    else:
        phagesdbClusterUpdate = phagesdbSubcluster


    #Compare Host genus
    if phameratorHost != phagesdbHost:
        write_out(report_file,"\nError: Phamerator host %s and phagesdb host %s do not match for phageID %s." %(phameratorHost,phagesdbHost,phameratorId))
        total_errors += 1
        corrections_needed += 1
     


    #Compare Cluster and Subcluster
    
    if phagesdbSubcluster == 'Unspecified':
    
        if phameratorCluster != phagesdbCluster:
            print "Phamerator Cluster: " + phameratorCluster
            print "Phagesdb Cluster: " + phagesdbCluster
            write_out(report_file,"\nError: Phamerator Cluster %s does not match with phagesdb Cluster %s for phageID %s." %(phameratorCluster,phagesdbCluster,phameratorId))
            total_errors += 1
            corrections_needed += 1
        
    elif phameratorCluster != phagesdbSubcluster:
            print "Phamerator Cluster: " + phameratorCluster
            print "Phagesdb Subcluster: " + phagesdbSubcluster
            write_out(report_file,"\nError: Phamerator Cluster %s does not match with phagesdb Subcluster %s for phageID %s." %(phameratorCluster,phagesdbSubcluster,phameratorId))
            total_errors += 1
            corrections_needed += 1


    #If errors in the Host or Cluster information were identified, create an import ticket to for the import script to implement.
    if corrections_needed > 0:
        import_table_writer.writerow(["update",phameratorId,phagesdbHost,phagesdbClusterUpdate,phameratorStatus,"none","none"])
    
          



write_out(report_file,"\nMatched phage tally: %s." %matched_count)
write_out(report_file,"\nUnmatched phage tally: %s." %unmatched_count)
write_out(report_file,"\nUnmatched phages:")
for element in unmatched_phageId_list:
    write_out(report_file,"\n%s" %element)

    
    
    
#Close script.
write_out(report_file,"\n\n\n\nImport script completed.")  
report_file.close()
import_table_file.close()    
    
    
    
