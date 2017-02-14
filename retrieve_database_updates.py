#!/usr/bin/env python
#Database comparison script
#University of Pittsburgh
#Travis Mavrich
#20170212
#The purpose of this script is to provide an interactive environment to retrieve several types of Phamerator updates\
#that will be implemented through import_phage.py script.
#This script has combined the following independent scripts:
#1. compare_databases.py
#2. retrieve_draft_genomes.py
#3. retrieve_phagesdb_flatfiles.py
#4. retrieve_ncbi_phage.py



#Flow of the import process:
#1 Import python modules, set up variables and functions
#2 Option 1: compare Phamerator and phagesdb databases
#3 Option 2: retrieve auto-annotated files from PECAAN
#4 Option 3: retrieve manually-annotated files from phagesdb
#5 Option 4: retrieve updated files from NCBI





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
    updateFileDir = sys.argv[2] #What is the directory into which the report should go
    phage_file = sys.argv[3] #What is the name of the file that contains the list of new available phage annotations?
    
except:
    print "\n\n\
            This is a python script to retrieve several types of Phamerator database updates.\n\
                1. It compares Phamerator and phagesdb databases to identify inconsistencies.\n\
                2. It retrieves auto-annotated Genbank-formatted flatfiles from PECAAN.\n\
                3. It retrieves manually annotated Genbank-formatted flatfiles from phagesdb.\n\
                4. It retrieves updated Genbank-formatted flatfiles from NCBI.\n\n\n\
            It requires three arguments:\n\
            First argument: name of MySQL database that will be checked (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where all reports, update import tables, and retrieved files will be generated.\n\
            Third argument: list of phages that have manually-annotated files available on phagesdb (csv-formatted):\n\
                    1. phagesdb phage name\n\
                    (Indicate 'none' if this option will not be requested)\n\n\n\
            All retrieval options create a genomes folder and an import table if updates are available (csv-formatted):\n\
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




#Check to see if the user has indicated a phage file.
#If so, verify the path exists for the phage list file
if phage_file.lower() != "none":

    #Expand the path if it references the home directory
    if phage_file[0] == "~":
        phage_file = home_dir + phage_file[1:]

    #Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
    phage_file = os.path.abspath(phage_file)

    if os.path.exists(phage_file) == False:
        print "\n\nInvalid input for phage list file.\n\n"
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





#Define several functions

#Print out statements to both the terminal and to the output file
def write_out(filename,statement):
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






#Option 1: Compare Phamerator and phagesdb databases
compare_databases = select_option("\nDo you want to compare Phamerator and phagesdb databases for updates and changes? ")

if compare_databases == "yes":


    #Decide whether to compare nucleotide sequences or not.
    #Since this may take a long time, and generate thousands of requests to phagesdb, it may not need to be performed very often.
    compare_sequences = select_option("\nDo you want to compare all genome sequences? (yes or no) ")
        
    
    #Decide whether to retrieve all Genbank-formatted files from phagesdb or not.
    retrieve_phagesdb_genomes = select_option("\nDo you want to retrieve manually-annotated genomes from phagesdb? (yes or no) ")

    if retrieve_phagesdb_genomes == "yes":

        #Create output directories
        phagesdb_output_folder = '%s_retrieved_phagesdb_flatfiles' % date
        phagesdb_output_path = os.path.join(updateFileDir,phagesdb_output_folder)

        try:
            os.mkdir(phagesdb_output_path)
        except:
            print "\nUnable to create output folder: %s" %phagesdb_output_path
            sys.exit(1)
        os.chdir(phagesdb_output_path)


        #Open file to create import table with changes that need to be implemented
        phagesdb_import_table_file = open(os.path.join(phagesdb_output_path,date + "_phagesdb_import_table.csv"), "w")
        phagesdb_import_table_writer = csv.writer(phagesdb_import_table_file)

        #Create the output folder to hold the genome files
        os.mkdir(os.path.join(phagesdb_output_folder,genomes_folder)

 
        #Initialize phagesdb retrieval variables 
        phagesdb_retrieved_tally = 0
        phagesdb_failed_tally = 0
        phagesdb_retrieved_list = []
        phagesdb_failed_list = []

 
 
 
        
    #Create output directories

    comparison_output_folder = '%s_database_comparison' % date
    comparison_output_path = os.path.join(updateFileDir,comparison_output_folder)


    try:
        os.mkdir(os.path.join(updateFileDir,comparison_output_folder))
    except:
        print "\nUnable to create output folder: %s" % os.path.join(updateFileDir,comparison_output_folder)
        sys.exit(1)
        
    os.chdir(comparison_output_path)



    #Create a folder named "genomes"
    #This is a dummy folder, since the comparison script does not retrieve any files.
    #However, creating the folder makes it easier to run the import_script on the corrections_import_table, since this script relies on the presence of a genomes folder.
    os.mkdir(genomes_folder)


    #Open file to record update information
    report_file = open(os.path.join(comparison_output_path,date + "_database_comparison.txt"), "w")
    write_out(report_file,date + " Database comparison:\n\n\n")


    #Open file to create import table with changes that need to be implemented
    import_table_file = open(os.path.join(comparison_output_path,date + "_corrections_import_table.csv"), "w")
    import_table_writer = csv.writer(import_table_file)


    #Retrieve database version
    #Retrieve current data in database
    #0 = PhageID
    #1 = Name
    #2 = HostStrain
    #3 = Sequence
    #4 = status
    #5 = Cluster
    #6 = SequenceLength
    try:
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()
    except:
        print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
        import_table_file.close()
        report_file.close()
        sys.exit(1)

    try:
        cur.execute("START TRANSACTION")
        cur.execute("SELECT version FROM version")
        db_version = str(cur.fetchone()[0])
        cur.execute("SELECT PhageID,Name,HostStrain,Sequence,status,Cluster,SequenceLength,DateLastModified FROM phage")
        current_genome_data_tuples = cur.fetchall()
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)

    except:
        import_table_file.close()
        report_file.close()
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
        phameratorSize = int(genome_tuple[6])
        phameratorDate = genome_tuple[7]

        #In Phamerator, Singleton Clusters are recorded as '\N', but in phagesdb they are recorded as "Singleton"
        if phameratorCluster is None:
            phameratorCluster = 'Singleton'
            
        matched_phagesdb_data = ""
        
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
            write_out(report_file,"\nError: unable to find phageID %s or phageName %s from phagesdb." %(phameratorId,phameratorName))
            matched_phagesdb_data = ""
            unmatched_count += 1
            unmatched_phageId_list.append(phameratorId)
            total_errors += 1
            continue


        #Retrieve specific matched data
        phagesdbName = matched_phagesdb_data['phage_name']
        phagesdbHost = matched_phagesdb_data['isolation_host']['genus']
        phagesdbSize = int(matched_phagesdb_data['genome_length'])


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
        
              

        #Compare recorded genome sizes
        if phameratorSize != phagesdbSize:
            write_out(report_file,"\nError: Phamerator genome size %s does not match phagesdb genome size %s for phageID %s." %(phameratorSize,phagesdbSize,phameratorId))
            total_errors += 1








        ###IN PROGRESS
        ###IF NEW GB FILES SHOULD BE RETRIEVED FROM PHAGESDB, ADD THAT CODE HERE:
        
        #Retrieve all Genbank-formatted files that have been uploaded to phagesdb more recently than the genome was uploaded into Phamerator (if user selected this option)
        if retrieve_phagesdb_genomes == "yes":
        
        
        


            #os.chdir(genomes_folder)



            ###Check if gb file is available
            phagesdbName
            phameratorName

            #Not all phages have associated Genbank-formatted files available on phagesdb
            #Check to see if there is a flatfile for this phage
            if online_data_dict["qced_genbank_file"] is None:

                print "No flatfile found for phageID %s and phageName %s." %(phameratorId,phameratorName)
                phagesdb_failed_tally += 1
                phagesdb_failed_list.append(phameratorId)

            else:

                #Some phages may have a file but no associated date tagged with that file (since date tagging has only recently been implemented).
                #Download the flatfile only if there is a date tag, and only if that date is more recent than the date stored in Phamerator for that genome.
                #The tagged date only reflects when the file was uploaded into phagesdb. The date the actual Genbank record was created is stored within the file,
                #and this too could be less recent than the current version in Phamerator; however, this part gets checked during the import stage.
                if online_data_dict[NAME_OF_NEW_DATE_FIELD] is None:

                    print "Available flatfile does not have a date tag for phageID %s and phageName %s." %(phameratorId,phameratorName)
                    phagesdb_failed_tally += 1
                    phagesdb_failed_list.append(phameratorId)

                elif not online_data_dict[NAME_OF_NEW_DATE_FIELD] > phameratorDate:

                    print "The date of the available flatfile is not more recent than the current version in Phamerator for phageID %s and phageName %s." %(phameratorId,phameratorName)
                    phagesdb_failed_tally += 1
                    phagesdb_failed_list.append(phameratorId)

                else:

                    #Save the file on the hard drive with the same name as stored on phagesdb
                    flatfile_url = online_data_dict["qced_genbank_file"]
                    phagesdb_file = flatfile_url.split('/')[-1]

                
                
                ###CURRENT STOPPING POINT
                
                
                    try:
                        response = urllib2.urlopen(flatfile_url)
                        phagesdb_file_handle = open(phagesdb_file,'w')
                        phagesdb_file_handle.write(response.read())
                        response.close()
                        phagesdb_file_handle.close()
                        
                        
                        #Create the new import ticket
                        #0 = Import action
                        #1 = New phageID
                        #2 = HostStrain
                        #3 = Cluster
                        #4 = Status
                        #5 = Gene description field
                        #6 = Phage to replace
                        import_table_writer.writerow(["replace",phameratorId,"retrieve","retrieve","final","product","unspecified"])

                        retrieved_tally += 1
                        retrieved_list.append(phameratorId)

                    except:
                        print "Error: unable to retrieve %s genome." %phameratorId
                        failed_tally += 1
                        failed_list.append(phameratorId)




            #Report retrieval results
            if retrieved_tally > 0:
                print "\n\nThe following %s phage(s) were successfully retrieved:" %retrieved_tally
                for element in retrieved_list:
                    print element
            else:
                print "No new flatfiles available."


            if failed_tally > 0:
                print "\n\nThe following %s phage(s) failed to be retrieved:" %failed_tally
                for element in failed_list:
                    print element
            else:
                print "No phages failed to be retrieved."
            


            print "\nDone retrieving manually-annotated genomes from phagesdb.\n\n\n"    
            phage_file_handle.close()
            import_table_file.close()        
            os.chdir('..')        




        
        
        
        
        ###IN PROGRESS



        #Compare genome sequence (if this option was selected)
        #phagesdb stores the 'official' nucleotide sequence, so make sure the sequence in Phamerator matches the sequence in phagesdb
        if compare_sequences == "yes":
            
            #First retrieve the fasta file containing the sequence
            
            #Check to see if there is a fasta file stored on phagesdb for this phage
            if matched_phagesdb_data["fasta_file"] is None:
                write_out(report_file,"\nError: no fasta file found on phagesdb for phageID %s." %phameratorID)
                total_errors += 1            
            else:
            
                fastafile_url = matched_phagesdb_data["fasta_file"]                        
                response = urllib2.urlopen(fastafile_url)
                retrieved_fasta_file = response.read()
                response.close()
                
                #All sequence rows in the fasta file may not have equal widths, so some processing of the data is required
                #If you split by newline, the header is retained in the first list element
                split_fasta_data = retrieved_fasta_file.split('\n')
                phagesdbSequence = ""
                index = 1
                while index < len(split_fasta_data):
                    phagesdbSequence = phagesdbSequence + split_fasta_data[index].strip() #Remove any whitespace before appending, such as '\r'
                    index += 1

                phagesdbSequence_size = len(phagesdbSequence)
                
                #Compare phagesdb recorded genome size with genome size based on the stored fasta sequence
                if phagesdbSize != phagesdbSequence_size:
                    write_out(report_file,"\nError: phagesdb genome size %s does not match fasta sequence genome size %s for phagesdb phageName %s." %(phagesdbSize,phagesdbSequence_size,phagesdbName))
                    total_errors += 1
                
                #Compare phamerator recorded genome size with genome size based on the stored fasta sequence
                if phameratorSize != phagesdbSequence_size:
                    write_out(report_file,"\nError: Phamerator genome size %s does not match fasta sequence genome size %s for phageID %s." %(phameratorSize,phagesdbSequence_size,phameratorId))
                    total_errors += 1

                #Compare genome sequences stored in Phamerator and phagesdb fasta file
                if phameratorSequence.lower() != phagesdbSequence.lower():
                    write_out(report_file,"\nError: Genome sequences stored in Phamerator and phagesdb do not match for phageID %s." %phameratorId)
                    total_errors += 1

    write_out(report_file,"\nMatched phage tally: %s." %matched_count)
    write_out(report_file,"\nUnmatched phage tally: %s." %unmatched_count)
    write_out(report_file,"\nUnmatched phages:")
    for element in unmatched_phageId_list:
        write_out(report_file,"\n%s" %element)

    
    print "\nDone comparing databases.\n\n\n"
    report_file.close()
    import_table_file.close()    
    os.chdir('..')











#Option 2: Retrieve auto-annotated genomes from PECAAN
retrieve_pecaan_genomes = select_option("\nDo you want to retrieve auto-annotated genomes from PECAAN? (yes or no) ")


if retrieve_pecaan_genomes == "yes":


    #Create output directories
    pecaan_output_folder = '%s_retrieved_pecaan_genomes' % date
    pecaan_output_path = os.path.join(updateFileDir,pecaan_output_folder)


    try:
        os.mkdir(pecaan_output_path)
    except:
        print "\nUnable to create output folder: %s" %pecaan_output_path
        sys.exit(1)

    os.chdir(pecaan_output_path)

    #Retrieve list of unphamerated genomes
    #Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
    phage_list_url = 'http://phagesdb.org/data/unphameratedlist'
    phagesdb_response = urllib2.urlopen(phage_list_url)



        
    #Open file to create import table with changes that need to be implemented
    import_table_file = open(os.path.join(pecaan_output_path,date + "_pecaan_import_table.csv"), "w")
    import_table_writer = csv.writer(import_table_file)




    #Create the output folder to hold the genome files
    os.mkdir(genomes_folder)
    os.chdir(genomes_folder)




    #Retrieve auto-annotated genomes from PECAAN
    pecaan_prefix = 'https://discoverdev.kbrinsgd.org/phameratoroutput/phage/'
    retrieved_tally = 0
    failed_tally = 0
    retrieved_list = []
    failed_list = []


    #Iterate through each row in the file
    for new_phage in phagesdb_response:


        #PECAAN should be able to generate any phage that is listed on phagesdb
        new_phage = new_phage.strip() #Remove \t character at the end of each row
        pecaan_link = pecaan_prefix + new_phage
        pecaan_file = new_phage + "_Draft.txt"
        #print pecaan_link
        try:
            pecaan_response = urllib2.urlopen(pecaan_link)
            pecaan_file_handle = open(pecaan_file,'w')
            pecaan_file_handle.write(pecaan_response.read())
            pecaan_response.close()
            pecaan_file_handle.close()
            
            
            #Create the new import ticket
            #0 = Import action
            #1 = New phageID
            #2 = HostStrain
            #3 = Cluster
            #4 = Status
            #5 = Gene description field
            #6 = Phage to replace
            import_table_writer.writerow(["add",new_phage,"retrieve","retrieve","draft","product","none"])
            print "Retrieved %s from PECAAN." %new_phage
            retrieved_tally += 1
            retrieved_list.append(new_phage)

        except:
            print "Error: unable to retrieve %s draft genome." %new_phage
            failed_tally += 1
            failed_list.append(new_phage)


    phagesdb_response.close()


    #Report results
    if retrieved_tally > 0:
        print "The following %s phage(s) were successfully retrieved:" %retrieved_tally
        for element in retrieved_list:
            print element
    else:
        print "No new draft genomes available."


    if failed_tally > 0:
        print "The following %s phage(s) failed to be retrieved:" %failed_tally
        for element in failed_list:
            print element
    else:
        print "No phages failed to be retrieved."


    print "\nDone retrieving auto-annotated genomes from PECAAN.\n\n\n"
    import_table_file.close()
    os.chdir('..')







#Option 3: Retrieve manually-annotated genomes from phagesdb
if phage_file.lower() != "none":
    retrieve_phagesdb_genomes = select_option("\nDo you want to retrieve manually-annotated genomes from phagesdb? (yes or no) ")
else:
    retrieve_phagesdb_genomes = "no"


if retrieve_phagesdb_genomes == "yes":


    #Create output directories
    phagesdb_output_folder = '%s_retrieved_phagesdb_flatfiles' % date
    phagesdb_output_path = os.path.join(updateFileDir,phagesdb_output_folder)


    try:
        os.mkdir(phagesdb_output_path)
    except:
        print "\nUnable to create output folder: %s" %phagesdb_output_path
        sys.exit(1)

    os.chdir(phagesdb_output_path)


    #Retrieve list of genomes with new annotations available
    #Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
    phage_file_handle = open(phage_file,'r')
    phage_file_reader = csv.reader(phage_file_handle)




    #Open file to create import table with changes that need to be implemented
    import_table_file = open(os.path.join(phagesdb_output_path,date + "_phagesdb_import_table.csv"), "w")
    import_table_writer = csv.writer(import_table_file)



    #Create the output folder to hold the genome files
    os.mkdir(genomes_folder)
    os.chdir(genomes_folder)



    #Phagesdb API to retrieve genome information
    api_prefix = "http://phagesdb.org/api/phages/"
    api_suffix = "/?format=json"

    retrieved_tally = 0
    failed_tally = 0
    retrieved_list = []
    failed_list = []



    for new_phage in phage_file_reader:

        #Retrieve phage-specific data from phagesdb
        #Note: urlopen is not case sensitive, so while this script tests for correct spelling, it does not test for correct capitalization.
        #However, the import script tests for that.    
        phage_url = api_prefix + new_phage[0] + api_suffix
        online_data_json = urllib.urlopen(phage_url)
        online_data_dict = json.loads(online_data_json.read())


        #If the phage is not found in phagesdb, then the phage url will return a dictionary with a single key ("detail") and value ("not found.")
        #Not sure if "detail" key is present in other phage urls, but they probably do not have the same value
        try:
            if online_data_dict["detail"] == "Not found.":
                print "Unable to locate URL for phage %s." %new_phage[0]
                failed_tally += 1
                failed_list.append(new_phage[0])
                continue
            else:
                print "URL found for phage %s." %new_phage[0]
        except:
            print "URL found for phage %s." %new_phage[0]

        #Check to see if there is a flatfile stored on phagesdb for this phage
        if online_data_dict["qced_genbank_file"] is None:

            print "Error: no flatfile found for phage %s." %new_phage[0]
            failed_tally += 1
            failed_list.append(new_phage[0])

        else:

            flatfile_url = online_data_dict["qced_genbank_file"]

            #Save the file on the hard drive with the same name as stored on phagesdb
            phagesdb_file = flatfile_url.split('/')[-1]

        
            try:
                response = urllib2.urlopen(flatfile_url)
                phagesdb_file_handle = open(phagesdb_file,'w')
                phagesdb_file_handle.write(response.read())
                response.close()
                phagesdb_file_handle.close()
                
                
                #Create the new import ticket
                #0 = Import action
                #1 = New phageID
                #2 = HostStrain
                #3 = Cluster
                #4 = Status
                #5 = Gene description field
                #6 = Phage to replace
                import_table_writer.writerow(["replace",new_phage[0],"retrieve","retrieve","final","product","unspecified"])

                retrieved_tally += 1
                retrieved_list.append(new_phage[0])

            except:
                print "Error: unable to retrieve %s genome." %new_phage[0]
                failed_tally += 1
                failed_list.append(new_phage[0])




    #Report retrieval results
    if retrieved_tally > 0:
        print "\n\nThe following %s phage(s) were successfully retrieved:" %retrieved_tally
        for element in retrieved_list:
            print element
    else:
        print "No new flatfiles available."


    if failed_tally > 0:
        print "\n\nThe following %s phage(s) failed to be retrieved:" %failed_tally
        for element in failed_list:
            print element
    else:
        print "No phages failed to be retrieved."
    


    print "\nDone retrieving manually-annotated genomes from phagesdb.\n\n\n"    
    phage_file_handle.close()
    import_table_file.close()        
    os.chdir('..')
    
    
    
    


#Option 4: Retrieve updated records from NCBI
retrieve_ncbi_genomes = select_option("\nDo you want to retrieve updated NCBI records? (yes or no) ")

if retrieve_ncbi_genomes == "yes":


    #Flow of the NCBI record retrieval process:
    #1 Retrieve current database information and create list of phages to check for updates at NCBI
    #2 Using esearch, verify the accessions are valid
    #3 Retrieve valid records in batches
    #4 Check which records are newer than the upload date of the current version in phamerator
    #5 Save new records in a folder and create an import table for them




    #Create output directory and processing file
    ncbi_output_folder = '%s_retrieved_ncbi_flatfiles' % date
    ncbi_output_path = os.path.join(updateFileDir,ncbi_output_folder)

    try:
        os.mkdir(ncbi_output_path)
    except:
        print "\nUnable to create output folder: %s" % ncbi_output_path
        sys.exit(1)

    os.chdir(ncbi_output_path)





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



    #1 Retrieve current database information and create list of phages to check for updates at NCBI

    #Retrieve database version
    #Retrieve current data in database
    #0 = PhageID
    #1 = Name
    #2 = HostStrain
    #3 = Cluster
    #4 = status
    #5 = accession
    #6 = date last modified
    #7 = retrieve record

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
        cur.execute("SELECT PhageID,Name,HostStrain,Cluster,status,Accession,DateLastModified,RetrieveRecord FROM phage")
        current_genome_data_tuples = cur.fetchall()
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)

    except:
        mdb_exit("\nUnable to access the database to retrieve genome information.")

    con.close()




    #Initialize tally variables
    tally_total = 0
    tally_not_updated = 0
    tally_no_accession = 0
    tally_duplicate_accession = 0
    tally_retrieval_failure = 0
    tally_retrieved_not_new = 0
    tally_retrieved_for_update = 0
    tally_total = len(current_genome_data_tuples)


    #Create file to write results to
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

    #Add to dictionary if 1) the genome is set to be automatically updated and 2) if there is an accession number
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
        if phage_list[7] != 1:
            print "PhageID %s is not set to be automatically updated by NCBI record." %phage_list[0]
            tally_not_updated += 1
            processing_results_file_writer.writerow([phage_list[0],phage_list[1],phage_list[5],phage_list[4],phage_list[6],'NA','no automatic update'])

        elif phage_list[5] == "" or phage_list[5] is None:
            print "PhageID %s is set to be automatically update, but it does not have accession number." %phage_list[0]
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
        if element[5] in unique_accession_dict.keys():
            temp_list.append(unique_accession_dict.pop(element[5]))


    #Now add these elements from dictionary to the duplicate data list
    for element in temp_list:
        duplicate_accession_list.append(element)
        

    #Output the duplicate data
    tally_duplicate_accession = len(duplicate_accession_list)    
    for data_list in duplicate_accession_list:
        processing_results_file_writer.writerow([data_list[0],data_list[1],data_list[5],data_list[4],data_list[6],'NA','duplicate accession'])

        



    #2 & #3 Use esearch to verify the accessions are valid and efetch to retrieve the record
    Entrez.email = contact_email
    Entrez.tool = "NCBIRecordRetrievalScript"



    #Create batches of accessions
    unique_accession_list = unique_accession_dict.keys()

    #Add [ACCN] field to each accession number
    index = 0
    while index < len(unique_accession_list):
        unique_accession_list[index] = unique_accession_list[index] + "[ACCN]"
        index += 1


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



    #4 Now that all records have been retrieved, check which records are newer than the upload date of the current version in phamerator.
    # Create the genbank-formatted file only if it is a newer genome
    # Also create an import table
    #0 = Action = replace
    #1 = Name = current phamerator PhageID
    #2 = HostStrain = current phamerator hoststrain
    #3 = Cluster = current phamerator cluster
    #4 = status = final
    #5 = Gene Description Field = product
    #6 = Genome to replace = current phamerator PhageID
    import_table_file = open(os.path.join(ncbi_output_path,date + "_ncbi_import_table.csv"), "w")
    import_table_file_writer = csv.writer(import_table_file)


    #Create the output folder to hold the genome files
    os.mkdir(genomes_folder)
    os.chdir(genomes_folder)



    tally_retrieval_failure = len(retrieval_error_list)
    for retrieval_error_accession in retrieval_error_list:

        phamerator_data = unique_accession_dict[retrieval_error_accession]
        processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],'NA','retrieval failure'])




    for retrieved_record in retrieved_record_list:
        retrieved_record_accession = retrieved_record.name

        #Convert date date to datetime object
        retrieved_record_date = retrieved_record.annotations["date"]
        retrieved_record_date = datetime.strptime(retrieved_record_date,'%d-%b-%Y')


        #MySQL outputs the DateLastModified as a datetime object
        phamerator_data = unique_accession_dict[retrieved_record_accession]

        #5 Save new records in a folder and create an import table row for them
        if retrieved_record_date > phamerator_data[6]:

            print 'Retrieved record date %s is more recent than phamerator date %s.' %(retrieved_record_date,phamerator_data[6])
            tally_retrieved_for_update += 1
            processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],retrieved_record_date,'update record'])


            #Now output genbank-formatted file to be uploaded to Phamerator and create the import table action
            #First remove the "_Draft" suffix if it is present
            if phamerator_data[0][-6:].lower() == '_draft':
                import_table_name = phamerator_data[0][:-6]
            else:
                import_table_name = phamerator_data[0]
            
            SeqIO.write(retrieved_record, phamerator_data[1].lower() + "__" + retrieved_record_accession + ".gb","genbank")
            import_table_data_list = ['replace',import_table_name,phamerator_data[2],phamerator_data[3],'final','product',phamerator_data[0]]
            import_table_file_writer.writerow(import_table_data_list)


        else:
            print 'Phamerator date %s is more recent than retrieved record date %s.' %(phamerator_data[6],retrieved_record_date)
            tally_retrieved_not_new += 1
            processing_results_file_writer.writerow([phamerator_data[0],phamerator_data[1],phamerator_data[5],phamerator_data[4],phamerator_data[6],retrieved_record_date,'record not new'])        
            


    #Print summary of script
    print "Number of genomes in Phamerator: %s" %tally_total
    print "Number of genomes that are NOT set to be updated: %s" %tally_not_updated
    print "Number of auto-updated genomes with no accession: %s" %tally_no_accession
    print "Number of duplicate accessions: %s" %tally_duplicate_accession
    print "Number of records that failed to be retrieved: %s" %tally_retrieval_failure
    print "Number of records retrieved that are NOT more recent than Phamerator record: %s" %tally_retrieved_not_new
    print "Number of records retrieved that should be updated in Phamerator: %s" %tally_retrieved_for_update


    processing_check = tally_total - tally_not_updated - tally_no_accession - tally_duplicate_accession - tally_retrieval_failure - tally_retrieved_not_new - tally_retrieved_for_update
    if processing_check != 0:
        print "Processing check: %s" %processing_check
        print "Error: the processing of phages was not tracked correctly."
        print "\n\n\n"



    print "\nDone retrieving updated records from NCBI.\n\n\n"
    import_table_file.close()
    processing_results_file_handle.close()
    os.chdir('..')









    
    
####Close script.
print "\n\n\nRetrieve updates script completed."  

















