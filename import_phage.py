#!/usr/bin/env python
#Phage Import Script
#Based on original script by Charles Bowman and Damani Brown on 20150216
#University of Pittsburgh
#Revamped by Travis Mavrich on 20160701
#Now it is designed to be a single script to upload genomes, update host/cluster/status data, remove genomes, and check integrity of files.
#WORKS FOR .gb, .gbf, .gbk, or .txt files




#Flow of the import process
#1 Import python modules, set up variables and functions
#2 Retrieve current database information
#3 Retrieve import information from file
#4 Update genome information for genomes already in the database
#5 Remove genomes from database that have no replacements
#6 Add or replace genomes

#The strategy for handling MySQL errors:
#The database is accessed in Steps 2, 4, 5, and 6. Connections are opened and closed at each Step. This ensures that any unforeseen errors encountered elsewhere in the script
#that force the script to exit do not cause ROLLBACK errors or connection errors.


import time, sys, os, getpass, csv, re
from Bio import SeqIO
import MySQLdb as mdb
from tabulate import tabulate

#Get the command line parameters
try:
    database = sys.argv[1]
    phageListDir = sys.argv[2]
    updateFile = sys.argv[3]
except:
    #print "Incorrect Parameters: ./import_phage.py DATABASE GENOME_DIRECTORY IMPORT_FILE"
    print "\n\n\
            This is a python script to import and update phage genomes in the Phamerator database.\n\
            It requires three arguments:\n\
            First argument: name of MySQL database that will be updated (e.g. 'Actino_Draft').\n\
            Second argument: directory path to the folder of genome files that will uploaded (genbank-formatted).\n\
            Third argument: directory path to the import table file with the following columns (csv-formatted):\n\
                    1. Action to implement on the database (add, remove, replace, update)\n\
                    2. PhageID to add or update\n\
                    3. Host genus of the updated phage\n\
                    4. Cluster of the updated phage\n\
                    5. Field that contains the gene description information (product, note, function)\n\
                    6. PhageID that will be removed or replaced\n\n"
    sys.exit(1)

#Set up MySQL parameters
mysqlhost = 'localhost'
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')



#Set up run type
run_type = ""
while (run_type != "test" and run_type != "production"):
    run_type = raw_input("Indicate run type (test or production): ")
    run_type = run_type.lower()
    









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
            write_out(output_file,message)
            number = 1
        else:
            print "Invalid response."
    #This value will be added to the current error total. If 0, no error was encountered. If 1, an error was encountered.
    return number
            
#Exits MySQL
def mdb_exit(message):
    write_out(output_file,"\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
    write_out(output_file,message)
    write_out(output_file,"\nThe import script did not complete.")
    write_out(output_file,"\nExiting MySQL.")
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    write_out(output_file,"\nExiting import script.")
    output_file.close()
    sys.exit(1)

#If phage Cluster is Singleton, make sure MySQL statement is created correctly
def create_cluster_statement(phage_name,cluster):
    cluster_statement = ""
    if cluster == "singleton":
        cluster_statement = "UPDATE phage SET Cluster = NULL" + " WHERE PhageID = '" + phage_name + "';"	
    else:
        cluster_statement = "UPDATE phage SET Cluster = '" + cluster + "' WHERE PhageID = '" + phage_name + "';"
    return cluster_statement



#Function to split gene description field
def retrieve_description(genbank_feature,description_field):
    description = genbank_feature.qualifiers[description_field][0].lower().strip()
    split_description = description.split(' ')                
    if description == "hypothetical protein":
        description = ""
    elif (split_description[0][:2] == "gp" and len(split_description) == 1):
        description = ""
    else:
        description = genbank_feature.qualifiers[description_field][0].strip()    
    return description













#Open file to record update information
date = time.strftime("%Y%m%d")
output_file = open("/tmp/" + date + "_phage_import_log_" + run_type + "_run.txt", "w")
write_out(output_file,date + " Phamerator database updates:\n\n\n")
write_out(output_file,"\n\n\n\nBeginning import script...")
write_out(output_file,"\nRun type: " + run_type)


#Retrieve database version
#Retrieve current data in database
#0 = PhageID
#1 = Name
#2 = HostStrain
#3 = Sequence
#4 = status
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
    cur.execute("SELECT PhageID,Name,HostStrain,Sequence,status,Cluster FROM phage")
    current_genome_data_tuples = cur.fetchall()
    cur.execute("SELECT GeneID,PhageID FROM gene")
    current_gene_data_tuples = cur.fetchall()
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)

except:
    mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")

con.close()

write_out(output_file,"\nDatabase: " + database)
write_out(output_file,"\nDatabase version: " + db_version)
write_out(output_file,"\nTotal phages in database before changes: " + str(len(current_genome_data_tuples)))

#Create data sets to compare new data with
#Originally, a phageName_set and phageSequence_set was implemented, but I ended up not using it. The SQL query still returns these values
#so if I need to re-implement those, I am able to.
phageId_set = set()
phageHost_set = set()
phageStatus_set = set()
phageCluster_set = set()
phageGene_set = set()
phageGene_dict = {}
print "Preparing genome data sets from the database..."
for genome_tuple in current_genome_data_tuples:
    phageId_set.add(genome_tuple[0])
    phageHost_set.add(genome_tuple[2])
    phageStatus_set.add(genome_tuple[4])
    phageCluster_set.add(genome_tuple[5])

    #Create a dictionary of GeneIDs
    #Key = PhageID
    #Value = set of GeneIDs
    tempGeneID_set = set()
    for gene_tuple in current_gene_data_tuples:
        if gene_tuple[1] == genome_tuple[0]:
            tempGeneID_set.add(gene_tuple[0])
    phageGene_dict[genome_tuple[0]] = tempGeneID_set


#Create a set of GeneIDs, regardless of corresponding PhageID
for gene_tuple in current_gene_data_tuples:
    phageGene_set.add(gene_tuple[0])

  

    


#Create set of all types of actions allowed using this script
#Add = add a new genome without removing another.
#Remove = delete a genome without adding another.
#Replace = delete a genome and replace it with another. Genome names can be different, but the DNA sequence cannot be different.
#Update = make changes to HostStrain, Cluster, or status fields of phages already in the database.
action_set = set(["add","remove","replace","update"]) 

#Create set of most common gene description genbank qualifiers
description_set = set(["product","note","function"])


#Retrieve import info from indicated import table file and read all lines into a list and verify contents are correctly populated.
#0 = Type of database action to be performed (add, remove, replace, update)
#1 = New phage name that will be added to database
#2 = Host of new phage
#3 = Cluster of new phage (singletons should be reported as "singleton")
#4 = Status of new phage
#5 = Feature field containing gene descriptions of new phage
#6 = PhageID of genome to be removed from the database
column_action_headers = ["All","Add/Replace/Update","Add/Replace/Update","Add/Replace/Update","Add/Replace/Update","Add/Replace","Replace/Remove"]
column_headers = ["Action","PhageID","HostStrain","Cluster","Status","DescriptionField","PhageID"]

write_out(output_file,"\n\n\n\nRetrieving import info from table in file...")

file_object = open(sys.argv[3],'r')
file_reader = csv.reader(file_object)
genome_data_list = []
table_errors = 0
add_total = 0
remove_total = 0
replace_total = 0
update_total = 0

for row in file_reader:

    #Make sure "none" indications are lowercase, as well as "action", "status", and "feature" fields are lowercase
    row[0] = row[0].lower()
    if row[1].lower() == "none":
        row[1] = row[1].lower()
    if row[2].lower() == "none":
        row[2] = row[2].lower()        
    if row[3].lower() == "none":
        row[3] = row[3].lower()        
    row[4] = row[4].lower()
    row[5] = row[5].lower()        
    if row[6].lower() == "none":
        row[6] = row[6].lower()
        
        
    #Make sure the requested action is permissible
    #This script currently only allows four actions: add, remove, replace, and update
    if row[0] in action_set:
        if row[0] == "add":
            add_total += 1
        elif row[0] == "remove":
            remove_total += 1
        elif row[0] == "replace":
            replace_total += 1
        elif row[0] == "update":
            update_total += 1
        else:
            pass
    else:
        write_out(output_file,"\nError: %s is not a permissible action." %row[0])
        table_errors += 1


        
    #Modify Host, Cluster, Status, and PhageName fields if needed
    if row[2] != "none":
        row[2] = row[2].split(' ')[0] #Keep only the genus in the host data field and discard the rest
        if row[2] not in phageHost_set:
            print "The host strain " + row[2] + " is not currently in the database."
            table_errors +=  question("\nError: %s is not the correct host for %s." % (row[2],row[1])) #errors will be incremented if host was not correct

    if row[3] != "none":
        if row[3].lower() == "singleton":
            row[3] = row[3].lower()
        if (row[3] not in phageCluster_set and row[3] != "singleton"):
            print "The Cluster " + row[3] + " is not currently in the database."
            table_errors +=  question("\nError: %s is not the correct Cluster for %s." % (row[3],row[1])) #errors will be incremented if host was not correct   

    if (row[4] not in phageStatus_set and row[4] != "none"):
            print "The status " + row[4] + " is not currently in the database."
            table_errors +=  question("\nError: %s is not the correct status for %s." % (row[4],row[1])) #errors will be incremented if host was not correct   
    if (row[4] == "draft" and row[1][-6:].lower() != "_draft"):
        row[1] = row[1] + "_Draft"

    if (row[5] not in description_set and row[5] != "none"):     
        print row[5] + " is an uncommon qualifier."
        table_errors += question("\nError: %s is an incorrect qualifier." % row[5])


    



    #Rules for how each field is populated differs depending on each specific action
    
    #Update
    if row[0] == "update":
        #FirstPhageID
        if row[1] not in phageId_set:
            write_out(output_file,"\nError: %s is not a valid PhageID in the database." %row[1])
            table_errors += 1
        #Host, Cluster, Status
        if (row[2] == "none" or row[3] == "none" or row[4] == "none"):
            write_out(output_file,"\nError: %s does not have correctly populated HostStrain, Cluster, or Status fields." %row[1])
            table_errors += 1
        
        #Description
        if row[5] != "none":
            write_out(output_file,"\nError: %s does not have correctly populated Description field." %row[1])
            table_errors += 1

        #SecondPhageID
        if row[6] != "none":
            write_out(output_file,"\nError: %s should not have a genome listed to be removed." %row[1])
            table_errors += 1

    
    #Add
    elif row[0] == "add":
        #FirstPhageID
        if row[1] in phageId_set:
            write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %row[1])
            table_errors += 1
        #FirstPhageID, Host, Cluster, Status, Description
        if (row[1] == "none" or row[2] == "none" or row[3] == "none" or row[4] == "none" or row[5] == "none"):
            write_out(output_file,"\nError: %s does not have correctly populated fields." %row[1])
            table_errors += 1
        #Status
        if row[4] == "final":
            print row[1] + " to be added is listed as Final status, but no Draft (or other) genome is listed to be removed."
            table_errors +=  question("\nError: %s is not the correct status for %s." % (row[4],row[1]))
            
        #SecondPhageID
        if row[6] != "none":
            write_out(output_file,"\nError: %s to be added should not have a genome indicated for removal." %row[1])
            table_errors += 1
    
    #Remove
    elif row[0] == "remove":
        #FirstPhageID,Host, Cluster, Status, Description
        if (row[1] != "none" or row[2] != "none" or row[3] != "none" or row[4] != "none" or row[5] != "none"):
            write_out(output_file,"\nError: %s to be removed does not have correctly populated fields." %row[6])
            table_errors += 1        
        #SecondPhageID
        if row[6] not in phageId_set:
            write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %row[6])
            table_errors += 1
        
    
    
    #Replace
    elif row[0] == "replace":

        #FirstPhageID
        if row[1] == "none":
            write_out(output_file,"\nError: %s is not a valid PhageID. %s genome cannot be replaced." % (row[1],row[6]))
            table_errors += 1
            
        #FirstPhageID. If replacing a genome, ensure that if the genome to be removed is not the same, that the new genome added has a unique name
        if (row[1] in phageId_set and row[1] != row[6]):
            write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %row[1])
            table_errors += 1
        #Host,Cluster,Status,Description
        if (row[2] == "none" or row[3] == "none" or row[4] == "none" or row[5] == "none"):
            write_out(output_file,"\nError: %s does not have correctly populated fields." %row[1])
            table_errors += 1
        #SecondPhageID
        if row[6] not in phageId_set:
            write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %row[6])
            table_errors += 1
            
    else:
        pass

        
    genome_data_list.append(row)
file_object.close









#Now that all rows have been added to the list, verify there are no duplicate actions
add_set = set()
remove_set = set()
action_add_set = set()
action_remove_set = set()
action_add_remove_set = set()

#Create each set and do initial checks for duplications. If the Add name or Remove name is "none", skip that because there are probably duplicates of those.
for genome_data in genome_data_list:
    current_add = (genome_data[1],)
    current_remove = (genome_data[6],)
    current_action_add = (genome_data[0],genome_data[1])
    current_action_remove = (genome_data[0],genome_data[6])
    current_action_add_remove = (genome_data[0],genome_data[1],genome_data[6])
    
    
    #First check the one-field and two-field combinations
    if current_add[0] != "none":
        if current_add in add_set:
            print genome_data[1] + " appears to be involved in more than one step."
            table_errors += question("\nError: %s is duplicated" % str(current_add)) #errors will be incremented if name is used more than once
        else:
            add_set.add(current_add)

        if current_action_add in action_add_set:
            write_out(output_file,"\nError: %s is duplicated" % str(current_action_add))
            table_errors += 1
        else:
            action_add_set.add(current_action_add)
        
    if current_remove[0] != "none":            
        if current_remove in remove_set:
            print genome_data[6] + " appears to be involved in more than one step."
            table_errors += question("\nError: %s is duplicated" % str(current_remove)) #errors will be incremented if name is used more than once
        else:
            remove_set.add(current_remove)
        
        if current_action_remove in action_remove_set:
            write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
            table_errors += 1
        else:
            action_remove_set.add(current_action_remove)

    #Now check the three-field combinations        
    if current_action_add_remove in action_add_remove_set:
        write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
        table_errors += 1
    else:
        action_add_remove_set.add(current_action_add_remove)
        
#Once the sets are created, also check if genomes to be removed are found in the Add field and vice versa.
for genome_data in genome_data_list:
    current_add = (genome_data[1],)
    current_remove = (genome_data[6],)
    
    #If the phage name is not replacing itself, the Add name is not expected to be in the Remove set and vice versa.
    if current_add != current_remove:
        if (current_add in remove_set and current_add != "none"):
            print genome_data[1] + " appears to be involved in more than one step."
            table_errors += question("\nError: %s is duplicated" % str(current_add)) #errors will be incremented if name is used more than once


        if (current_remove in add_set and current_remove != "none"):
            print genome_data[6] + " appears to be involved in more than one step."
            table_errors += question("\nError: %s is duplicated" % str(current_remove)) #errors will be incremented if name is used more than once

             







#Create separate lists of genome_data based on the indicated action: update, add/replace, remove
update_data_list = []
remove_data_list = []
add_replace_data_list = []
for genome_data in genome_data_list:
    if genome_data[0] == "update":
        update_data_list.append(genome_data)
    elif genome_data[0] == "remove":
        remove_data_list.append(genome_data)
    elif (genome_data[0] == "add" or genome_data[0] == "replace"):
        add_replace_data_list.append(genome_data)
    else:
        write_out(output_file,"\nError: during parsing of actions.")
        table_errors += 1




#Check to see if any genomes to be removed are the correct status
for genome_data in remove_data_list:
    for genome_tuple in current_genome_data_tuples:
        if (genome_data[6] == genome_tuple[0] and genome_tuple[4] != "draft"):
            print "The genome %s to be removed is currently %s status." % (genome_data[6],genome_tuple[4])
            table_errors += question("\nError: %s is %s status and should not be removed." % (genome_data[6],genome_tuple[4])) #errors will be incremented if name is used more than once



#If no erros encountered, print list of action items to be implemented and continue. Otherwise, exit the program.
if table_errors == 0:
    write_out(output_file,"\nImport table verified with 0 errors.")
    write_out(output_file,"\nList of all actions to be implemented:")
    write_out(output_file,"\nApplicable actions per field: " + str(column_action_headers))
    write_out(output_file,"\nField headers: " + str(column_headers))      
    for genome_data in genome_data_list:
        write_out(output_file,"\n" + str(genome_data))
    raw_input("Press ENTER to proceed to next import stage.")
else:
    write_out(output_file,"\n%s error(s) encountered with import file.\nNo changes have been made to the database." % table_errors)
    write_out(output_file,"\nExiting import script.")
    sys.exit(1)








#Update actions implemented
#Prepare to update fields that are not accompanied by adding or replacing a genome.
write_out(output_file,"\n\n\n\nUpdating Host, Cluster, and Status fields...")
updated = 0
update_statements = []
for genome_data in update_data_list:
    write_out(output_file,"\nPreparing: " + str(genome_data))

    #HostStrain and status updates.		
    update_statements.append("UPDATE phage SET HostStrain = '" + genome_data[2] + "' WHERE PhageID = '" + genome_data[1] + "';")
    update_statements.append("UPDATE phage SET status = '" + genome_data[4] + "' WHERE PhageID = '" + genome_data[1] + "';")

    #Create the statement to update Cluster.
    update_statements.append(create_cluster_statement(genome_data[1],genome_data[3]))

    updated += 1

#If it looks like there is a problem with some of the genomes on the list, cancel the transaction, otherwise proceed
if updated == update_total:
    
    if run_type == "production":
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()

        try:
            cur.execute("START TRANSACTION")
            for statement in update_statements:
                cur.execute(statement)
                write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
            cur.execute("COMMIT")
            write_out(output_file,"\nAll update statements committed.")
            cur.close()
            con.autocommit(True)
            
        except:
            mdb_exit("\nError: problem updating genome information.\nNo changes have been made to the database.")
            
        con.close()
    
    else:
        write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

else:
    write_out(output_file,"\nError: problem processing data list to update genomes. Check input table format.\nNo changes have been made to the database.")
    write_out(output_file,"\nExiting import script.")
    sys.exit(1)
    
write_out(output_file,"\nAll field update actions have been implemented.")
raw_input("Press ENTER to proceed to next import stage.")
    







#Remove actions implemented
#Prepare to remove any genomes that are not accompanied by a new genome
write_out(output_file,"\n\n\n\nRemoving genomes with no replacement...")
removed = 0
removal_statements = []
for genome_data in remove_data_list:
    write_out(output_file,"\nPreparing: " + str(genome_data))
    removal_statements.append("DELETE FROM phage WHERE PhageID = '" + genome_data[6] + "';")
    
    #Remove the GeneIDs associated with this genome from the reference set of all GeneIDs that is used later in the script. (It is NOT removing GeneIDs from the database).
    phageGene_set = phageGene_set - phageGene_dict.pop(genome_data[6])
    
    removed += 1

#If it looks like there is a problem with some of the genomes on the list, cancel the transaction, otherwise proceed	    
if removed == remove_total:

    if run_type == "production":
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()

        try:
            cur.execute("START TRANSACTION")
            for statement in removal_statements:
                cur.execute(statement)
                write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
            cur.execute("COMMIT")
            write_out(output_file,"\nAll remove statements committed.")        
            cur.close()
            con.autocommit(True)
            
        except:
            mdb_exit("\nError: problem removing genomes with no replacements.\nNo remove actions have been implemented.")

        con.close()

    else:
        write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

else:
    write_out(output_file,"\nError: problem processing data list to remove genomes. Check input table format.\nNo remove actions have been implemented.")
    sys.exit(1)
    
write_out(output_file,"\nAll genome remove actions have been implemented.")
raw_input("Press ENTER to proceed to next import stage.")


  









#Add and replace actions implemented
#Now that the data list does not contain update-only or remove-only data, create a dictionary. This serves to verify that all genomes to be imported are unique, as well as to
#be able to quickly retrieve information based on the genbank file that is opened.
#Key = PhageID
#Value = Genome data list:
add_replace_data_dict = {}
for genome_data in add_replace_data_list:    
   
    #Verify there are no duplicate PhageIDs. This was checked before in a slightly different way when looking for duplicate actions.
    if genome_data[1] not in add_replace_data_dict:
        add_replace_data_dict[genome_data[1]] = genome_data
    
    else:
        write_out(output_file,"\nError: problem creating genome data dictionary. PhageID %s already in dictionary. Check input table format." % genome_data[1])
        write_out(output_file,"\nNo add/replace actions have been implemented.")
        sys.exit(1)


#Iterate over each genbank-formatted genome file in the directory
write_out(output_file,"\n\n\n\nAccessing genbank-formatted files for add/replace actions...")
files =  [X for X in os.listdir(phageListDir) if os.path.isfile(os.path.join(phageListDir,X))]
failed_genome_files = []
failed_actions = []
for filename in files:
    
    write_out(output_file,"\n\nProcessing file: %s" % filename)
    seqFile = phageListDir + filename
    
    #If the file extension is not admissible, then skip. Otherwise, proceed
    if(str(seqFile[-2:]) != "gb" and str(seqFile[-3:]) != "gbf" and str(seqFile[-3:]) != "gbk" and str(seqFile[-3:]) != "txt"):
        failed_genome_files.append(filename)
        continue
    
    #The file is parsed to grab the header information
    for seq_record in SeqIO.parse(seqFile, "genbank"):
        
        add_replace_statements = []
        record_errors = 0
        geneID_set = set()
        missing_phage_name_tally = 0
        
        
        #Create a list to hold summary info on the genome record:
        #Record ID
        #Record Definition
        #Record Source
        #Record Organism
        #Organism qualifier of Source feature
        #List of CDS info:
            #Assigned geneID (how the geneID will look when uploaded to the database)
            #Locus tag
            #Product
            #Note
            #Function
            #Translation table
            #Translation
            #Processed Description (how the descriptions will look when uploaded to the database)
        record_summary_header = [["Header Field","Data"]]      
        
        
        #Make sure the file header attributes are successfully retrieved
        try:
            #Name
            record_organism = seq_record.annotations["organism"]
            phageName = record_organism.split(' ')[-1]
            #phageName = seq_record.annotations["organism"].split(' ')[-1]
            if phageName == "Unclassified.":
                phageName = seq_record.annotations["organism"].split(' ')[-2]
            
            #Sequence, length, and GC
            phageSeq = seq_record.seq.upper()
            seqLength = len(phageSeq)
            seqGC = 100 * (float(phageSeq.count('G')) + float(phageSeq.count('C'))) / float(seqLength)
            
        except:
            record_organism = ""
            phageName = "ERROR"
            phageSeq = ""
            seqLength = 0
            seqGC = 0
            print "Record does not have record organism information."
            write_out(output_file,"\nProblem with header block in: %s" % filename)
            record_errors += 1


        try:
            record_id = str(seq_record.id)
        except:
            record_id = ""
            print "Record does not have record ID information."
            record_errors += question("\nError: problem with header info of file %s." % filename)


        try:
            record_def = str(seq_record.description)
        except:
            record_def = ""
            print "Record does not have record definition information."
            record_errors += question("\nError: problem with header info of file %s." % filename)


        try:
            record_source = str(seq_record.annotations["source"])
        except:
            record_source = ""
            print "Record does not have record source information."
            record_errors += question("\nError: problem with header info of file %s." % filename)
        
        
        
        










        #Nonessential Stuff. PhageNotes and Accessions
        try:
            phageNotes = str(seq_record.annotations["comment"])
        except:
            phageNotes = ""

        try:
            accessionNum = seq_record.annotations["accessions"][0]
        except:
            accessionNum = ""
        
        
        
        #Retrieve Host,Cluster,Status,Description,Drop-genome info
        matchedData = add_replace_data_dict.pop(phageName,"error")
        if matchedData != "error":
            write_out(output_file,"\nPreparing: " + str(matchedData))
            phageAction = matchedData[0]
            phageHost = matchedData[2]
            phageCluster = matchedData[3]
            phageStatus = matchedData[4]
            cdsQualifier = matchedData[5]
            genomeReplace = matchedData[6]
        else:
            write_out(output_file,"\nError: problem matching phage %s in file %s to genome data from table. This genome was not added. Check input table format." % (phageName,filename))            
            failed_genome_files.append(filename)
            continue
        
           
           
        #Check the database to make sure the genome sequence and name matches correctly.            
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()
        try:
            cur.execute("START TRANSACTION")
            cur.execute("""SELECT PhageID,status FROM phage WHERE Sequence = "%s" """ % phageSeq)
            query_results = cur.fetchall()
            cur.execute("COMMIT")
            cur.close()
            con.autocommit(True)
                  
        except:
            mdb_exit("\nError: retrieving genome information from database while processing file %s.\nNot all genbank files were processed." % filename)

        con.close()
        
        
        
        #Cross-check the import action against the current state of the database and create SQL statements
        
        #If adding a new genome, no genome sequence in database is expected to match the current genome sequence
        if (phageAction == "add" and len(query_results) > 0):
            record_errors += 1 
            write_out(output_file,"\nError: these genome(s) in the database currently contain the same genome sequence: %s.\nUnable to upload %s." % (query_results,phageName))


        #If replacing a genome, exactly one and only one genome in the database is expected to have the same sequence.                
        elif phageAction == "replace":
            
            if len(query_results) > 1:
                record_errors += 1 
                write_out(output_file,"\nError: the following genomes in the database currently contain the same genome sequence as %s: %s).\nUnable to perform replace action." % (phageName,query_results))

            elif len(query_results) == 0: 
                write_out(output_file,"\n%s appears to be a different genome sequence than %s. These genomes do not match." % (phageName,genomeReplace))
                record_errors += question("\nError: %s and %s have different genome sequences." % (phageName,genomeReplace))
                           
            elif len(query_results) == 1:

                #The genome to be replaced is not Draft.
                if query_results[0][1].lower() != "draft":            
                    print "The genome in the database with matching sequence, %s, is listed as %s status." % (query_results[0][0],query_results[0][1])
                    record_errors +=  question("\nError: the genome in the database with matching sequence was incorrect status")
                
                #The genome to be replaced does not match the genome name in the database with the same sequence.
                if query_results[0][0] != genomeReplace:
                    write_out(output_file,"\nError: the genome to be removed, %s, does not match the genome name in database, %s, that has the matching genome sequence to %s." % (genomeReplace,query_results[0][0],phageName))
                    record_errors += 1
            else:
                pass
            
            add_replace_statements.append("DELETE FROM phage WHERE PhageID = '" + genomeReplace + "';")

        else:
            pass
           

        #Create statements   
        add_replace_statements.append("""INSERT INTO phage (PhageID, Accession, Name, HostStrain, Sequence, SequenceLength, GC, status, DateLastModified) VALUES ("%s","%s","%s","%s","%s",%s,%s,"%s","%s")""" % (phageName, accessionNum, phageName, phageHost, phageSeq, seqLength, seqGC, phageStatus,date))
        add_replace_statements.append(create_cluster_statement(phageName,phageCluster))
        
        
        #Retrieve the Source Feature info
        #feature_source_info = ""
        #for feature in seq_record.features:
        #    if feature.type == "source":
        #        feature_source_info = str(feature.qualifiers["organism"][0])
        #record_summary_header.append(["Source Feature Organism Qualifier",feature_source_info])
                 

        #Before processing CDS info, if a genome is being replaced, remove all of its associated GeneIDs from the reference set of all GeneIDs.
        if phageAction == "replace":
            phageGene_set = phageGene_set - phageGene_dict.pop(genomeReplace)







        #Next each CDS feature is parsed from the file
        #The cdsCount will increment for each CDS processed, even if it does not pass the QC filters. This way, the genome still retains the info that another CDS was originally present.
        cdsCount = 0
        addCount = 0
        transl_table_set = set()
        missing_locus_tag_tally = 0
        feature_note_tally = 0
        feature_product_tally = 0
        feature_function_tally = 0
        feature_source_info = ""
        
        
        
        
        
        record_summary_cds = [["Locus Tag","Product","Note","Function","Translation Table","Translation","Assigned GeneID","Assigned Description"]]
        for feature in seq_record.features:
            if feature.type != "CDS":
            
                #Retrieve the Source Feature info
                if feature.type == "source":           
                    feature_source_info = str(feature.qualifiers["organism"][0])
            
                continue

            else:                
                cdsCount += 1
                typeID = feature.type
            
                 
            #GeneID
            #Feature_locus_tag is a record of the locus tag found in the file. GeneID is what will be assigned in the database.
            try:
                feature_locus_tag = feature.qualifiers["locus_tag"][0]                
                geneID = feature_locus_tag
            except:
                feature_locus_tag = ""
                missing_locus_tag_tally += 1
                geneID = phageName + "_" + str(cdsCount)

            if (geneID not in geneID_set and geneID not in phageGene_set):
                geneID_set.add(geneID)
            else:
                write_out(output_file,"\nError: feature %s of %s is a duplicate geneID." % (geneID,phageName))
                record_errors += 1
                continue










            #Name
            if (geneID.split('_')[-1].isdigit()):
                geneName = geneID.split('_')[-1]
            else:
                geneName = cdsCount                        

            #Gene boundary coordinates
            if str(feature.location)[:4] == "join":
                temp = str(feature.location).strip("join{}").split(",")[0]
                startCoord = int(temp.split(":")[0].split("[")[1])
                temp = str(feature.location).strip("join{}").split(",")[1]
                stopCoord = int(temp.split(":")[1].split("]")[0])

            else:

                strStart = str(feature.location.start)
                strStop = str(feature.location.end)

                if (strStart.isdigit() and strStop.isdigit()):
                    startCoord = int(strStart)
                    stopCoord = int(strStop)
                else:
                    write_out(output_file,"\n" + geneID + " " + strStart + " " + strStop + " are non-traditional coordinates. This CDS will be skipped, but processing of the other genes will continue.")
                    record_errors += question("\nError: feature %s of %s does not have correct coordinates." % (geneID,phageName))
                    continue
                    
            #Translation, Gene Length (via Translation)
            try:
                translation = feature.qualifiers["translation"][0]
                #translation_trunc = translation[:5]
                geneLen = (len(translation) * 3) + 3  #Add 3 for the stop codon...
                                
            except:
                translation = ""
                #translation_trunc = ""
                geneLen = 0
                write_out(output_file,"\nError: problem with %s translation in phage %s." % (geneID,phageName))
                record_errors += 1
                
       
       
       
       
       
                
            #Translation table used
            try:
                feature_transl_table = feature.qualifiers["transl_table"][0]
                transl_table_set.add(feature_transl_table)
            except:
                feature_transl_table = ""
                write_out(output_file,"\nError: problem with %s translation table in phage %s." % (geneID,phageName))
                record_errors += 1
               




            #Gene Description
            #Use the feature qualifier from the import file.
            #If it is a generic description, leave field empty.
            #For generic 'gp#' descriptions, make sure there is no other info in the product that is valuabe by checking for other words present. 


            #First retrieve gene function, note, and product fields.
            try:
                feature_product = retrieve_description[feature,"product"]
                #feature_product = feature.qualifiers["product"][0]
                if feature_product != "":
                    feature_product_tally += 1
                #if len(feature_product) > 20:
                #    feature_product = feature_product[:20]
            except:
                feature_product = ""

            try:
                feature_note = retrieve_description[feature,"note"]
                #feature_note = feature.qualifiers["note"][0]
                if feature_note != "":
                    feature_note_tally += 1

                #if len(feature_note) > 20:
                #    feature_note = feature_note[:20]                
            except:
                feature_note = ""

            try:
                feature_function = retrieve_description[feature,"function"]
                #feature_function = feature.qualifiers["function"][0]
                if feature_function != "":
                    feature_function_tally += 1
                #if len(feature_function) > 20:
                #    feature_function = feature_function[:20]                
            except:
                feature_function = ""
            
            
            
            #Now assign the appropriate description info to the featureNote variable, as indicated from the import table.
            try:
                
                if cdsQualifier == "product":
                    featureNote = feature_product 
                
                elif cdsQualifier == "function":
                    featureNote = feature_function 

                elif cdsQualifier == "note":
                    featureNote = feature_note
                    
                else:
                    featureNote = retrieve_description(feature,cdsQualifier)
                
                #description = feature.qualifiers[cdsQualifier][0].lower().strip()
                #split_description = description.split(' ')                
                #if description == "hypothetical protein":
                #    featureNote = ""
                #elif (split_description[0][:2] == "gp" and len(split_description) == 1):
                #    featureNote = ""
                #else:
                #    featureNote = feature.qualifiers[cdsQualifier][0].strip()
                
            except:
                featureNote = ""








            #Orientation
            if feature.strand == 1:
                orientation = "Forward"
            elif feature.strand == -1:
                orientation = "Reverse"
            #ssRNA phages
            elif feature.strand is None:
                orientation = "Forward"
            else:
                print "Feature %s of %s does not have a common orientation. This CDS will be skipped, but processing of the other genes will continue." % (geneID,phageName)
                record_errors += question("\nError: feature %s of %s does not have correct orientation." % (geneID,phageName))
                continue

            addCount+= 1
            add_replace_statements.append("""INSERT INTO gene (GeneID, PhageID, Start, Stop, Length, Name, TypeID, translation, Orientation, Notes) VALUES ("%s","%s",%s,%s,%s,"%s","%s","%s","%s","%s");""" % (geneID, phageName,startCoord, stopCoord, geneLen, geneName, typeID, translation, orientation[0], featureNote)) 


            #Retrieve summary info to verify quality of file
            
            #featureNote_trunc = featureNote
            #if len(featureNote_trunc) > 20:
            #    featureNote_trunc = featureNote_trunc[:20]            
            
            record_summary_cds.append([feature_locus_tag,feature_product[:20],feature_note[:20],feature_function[:20],feature_transl_table,translation[:5],geneID,featureNote[:20]])
            
            
            
            
            
            
            
            
        #Now that all CDS features processed, process the source and organism fields to look or problems


        #Print the summary of the header information        
        record_summary_header.append(["Record ID",record_id])
        record_summary_header.append(["Record Defintion",record_def])
        record_summary_header.append(["Record Source",record_source])
        record_summary_header.append(["Record Organism",record_organism])
        record_summary_header.append(["Source Feature Organism Qualifier",feature_source_info])
        print "\nHere is the header summary information for %s from file %s:\n" % (phageName,filename)
        print tabulate(record_summary_header,headers = "firstrow")
        print "Double-check genome name and host name information."
        record_errors += question("\nError: problem with header info of file %s." % filename)
        print "\n\n\n"
        



        
        #See if there are any phage name typos in the header block
        pattern1 = re.compile(phageName)
        search_result = pattern1.search(record_id)
        if search_result == False:
        
            print "Record ID does not have identical phage name as found in the record organism field."
            record_errors += question("\nError: problem with header info of file %s." % filename)

        search_result = pattern1.search(record_def)
        if search_result == False:
        
            print "Record definition does not have identical phage name as found in the record organism field."
            record_errors += question("\nError: problem with header info of file %s." % filename)

        search_result = pattern1.search(record_source)
        if search_result == False:
        
            print "Record source does not have identical phage name as found in the record organism field."
            record_errors += question("\nError: problem with header info of file %s." % filename)






        #See if there are any host name typos in the header block
        phageHost_trim = phageHost
        if phageHost_trim == "Mycobacterium":
            phageHost_trim = phageHost_trim[:-3]
            
        pattern2 = re.compile(lower(phageHost_trim))
         
        search_result = pattern2.search(lower(record_def))
        if search_result == False:
        
            print "Record definition does not appear to have same host data as found in import table."
            record_errors += question("\nError: problem with header info of file %s." % filename)

        search_result = pattern2.search(lower(record_source))
        if search_result == False:
        
            print "Record source does not appear to have same host data as found in import table."
            record_errors += question("\nError: problem with header info of file %s." % filename)
      
        search_result = pattern2.search(lower(feature_source_info))
        if search_result == False:
        
            print "Source feature does not appear to have same host data as found in import table."
            record_errors += question("\nError: problem with header info of file %s." % filename)
        
        
        
        
        
        
        
        #Print record summary for all CDS information for quality control
        print "\nHere is the CDS summary information for %s from file %s:\n" % (phageName,filename)
        print tabulate(record_summary_cds,headers = "firstrow")
        print "Double-check assigned geneID, gene descriptions, and translations."        
        #record_errors += question("\nError: problem with CDS features of file %s." % filename)
        print "\n"


        #Check locus tag info:
        if missing_locus_tag_tally > 0:
            print "Phage %s from file %s is missing %s CDS locus tags." % (phageName, filename, missing_locus_tag_tally)
            record_errors += question("\nError: problem with locus tags in file  %s." % filename)


        #Check the phage name spelling in the locus tags
        pattern3 = re.compile(lower(phageName))
        geneID_typo_tally = 0
        for geneID in geneID_set:
           search_result = pattern3.search(lower(geneID))            
           if search_result == False:
                geneID_typo_tally += 1

        if geneID_tally > 0:
            print "There are %s geneID(s) that do not have the identical phage name included." % geneID_typo_tally
            record_errors += question("\nError: problem with locus tags of file %s." % filename)





        #Check all translation table info:
        if len(transl_table_set) == 1:
            transl_table_list = list(transl_table_set)
            if trans_table_list[0] != '11':
                write_out("The translation table used for %s is: %s." % (phageName,transl_table_list[0]))
                record_errors += question("\nError: phage %s does not use correct translation table." % phageName)
        
        else:
            write_out(output_file,"\nError: more than one translation table used in file %s." % filename)
            record_errors += 1







        #Check to ensure the best gene description field was retained
        write_out(output_file,"\nNumber of gene product descriptions found for phage %s: %s" % (phageName, feature_product_tally)
        write_out(output_file,"\nNumber of gene function descriptions found for phage %s: %s" % (phageName, feature_function_tally) 
        write_out(output_file,"\nNumber of gene note descriptions found for phage %s: %s" % (phageName, feature_note_tally)       
        
        
        if (cdsQualifier == "product" or cdsQualifier == "note":
            if feature_function_tally > 0:
                print "There are %s gene functions found. These will be ignored." % feature_function_tally
                record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)
        
        if (cdsQualifier == "product" or cdsQualifier == "function":
            if feature_note_tally > 0:
                print "There are %s gene notes found. These will be ignored." % feature_note_tally
                record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)

        if (cdsQualifier == "function" or cdsQualifier == "note":
            if feature_product_tally > 0:
                print "There are %s gene products found. These will be ignored." % feature_product_tally
                record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)





                
        #If errors were encountered with the file parsing, do not add to the genome. Otherwise, proceed.
        if record_errors == 0:
            diffCount = cdsCount - addCount
            write_out(output_file,"\nNo errors encountered while parsing: %s." % filename)
            write_out(output_file,"\nNumber of %s CDS features that failed: %s." % (phageName,diffCount))
        else:
            write_out(output_file,"\n%s errors were encountered with file %s. %s genome was not added to the database." % (record_errors,filename,phageName))
            failed_genome_files.append(filename)
            failed_actions.append(matchedData)
            continue

      
        #Execute SQL transactions        
        if run_type == "production":
            con = mdb.connect(mysqlhost, username, password, database)
            con.autocommit(False)
            cur = con.cursor()
            try:
                cur.execute("START TRANSACTION")
                for statement in add_replace_statements:
                    cur.execute(statement)
                    write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
                cur.execute("COMMIT")
                write_out(output_file,"\nAll add/replace statements for %s committed." % matchedData)
                cur.close()
                con.autocommit(True)

            except:
                mdb_exit("\nError: problem importing the file %s with the following add/replace action: %s.\nNot all genbank files were processed." % (filename,matchedData))

            con.close()
            
        else:
            write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)
        
        raw_input("Press ENTER to proceed to next file.")           

      

write_out(output_file,"\nAll files have been iterated through.")
raw_input("Press ENTER to proceed to next import stage.")        

        
#Final verifications
write_out(output_file,"\n\n\n\nFinal checks...")
write_out(output_file,"\n\nAll update actions have been implemented.")
write_out(output_file,"\n\nAll remove actions have been implemented.")

#If there were failed actions during the genome upload stage, add it back to the add_replace_data_dictionary.
#The add_replace_data_dictionary now contains a list of failed actions as well as unaddressed actions (actions that did not have matching genome files)
if len(failed_actions) > 0:
    for genome_data in failed_actions:
        add_replace_data_dict[genome_data[1]] = genome_data


#Verify all add/replace actions from the import csv file were addressed.
if len(add_replace_data_dict) > 0:
    write_out(output_file,"\n\nThe following add/replace action(s) in the import table were NOT successfully implemented:")
    for key in add_replace_data_dict:
        write_out(output_file,"\n" + str(add_replace_data_dict[key]))
else:
    write_out(output_file,"\nAll add/replace actions have been implemented.")  


#Verify that none of the genbank files failed to be uploaded.
if len(failed_genome_files) > 0:
    write_out(output_file,"\n\nThe following genbank files were NOT successfully processed:")
    for filename in failed_genome_files:
        write_out(output_file,"\n" + filename)
else:
    write_out(output_file,"\n\nAll genbank files were successfully processed.")







#Retrieve the total number of genomes now in the updated database
try:
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password.\nImport script was not completed."
    sys.exit(1)

try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT count(*) FROM phage")
    final_tally = cur.fetchall()
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)

except:
    mdb_exit("\nUnable to access the database to retrieve genome information.\nImport script was not completed.")

write_out(output_file,"\n\nTotal phages in database after changes: " + str(final_tally[0][0]))


#Close script.
write_out(output_file,"\n\n\n\nImport script completed.")  
output_file.close()






