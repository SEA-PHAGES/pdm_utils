#!/usr/bin/env python
#Retrieve manually annotated Genbank flatfiles from phagesdb
#University of Pittsburgh
#Travis Mavrich
#20170206
#The purpose of this script is to retrieve manually annotated Genbank-formatted flatfiles from phagesdb

#Third-party libraries
import MySQLdb as mdb


#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib, urllib2




#Get the command line parameters
try:
    updateFileDir = sys.argv[1] #What is the directory into which the report should go
except:
    print "\n\n\
            This is a python script to retrieve manually annotated Genbank flatfiles to import into Phamerator.\n\
            It requires one arguments:\n\
            First argument: directory path to where the phage list is stored and where new genomes and associated import table should be made (csv-formatted).\n\
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




#Create output directories
date = time.strftime("%Y%m%d")

output_folder = '%s_new_flatfiles' % date
output_path = os.path.join(updateFileDir,output_folder)


try:
    os.mkdir(output_path)
except:
    print "\nUnable to create output folder: %s" %output_path
    sys.exit(1)

os.chdir(output_path)





#Retrieve list of genomes with new annotations available
#Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
phage_file = '/home/cbowman/Documents/PhameratorDB_Management/Updates/temp/new_phages.txt'
phage_file_handle = open(phage_file,'r')
phage_file_reader = csv.reader(phage_file_handle,delimiter='\t')




#Open file to create import table with changes that need to be implemented
import_table_file = open(os.path.join(updateFileDir,output_folder,date + "_flatfile_import_table.csv"), "w")
import_table_writer = csv.writer(import_table_file)



#Phagesdb API to retrieve genome information
api_prefix = "http://phagesdb.org/api/phages/"
api_suffix = "/?format=json"

retrieved_tally = 0
failed_tally = 0
retrieved_list = []
failed_list = []



for new_phage in phage_file_reader:


    print "Attempting to retrieve %s from phagesdb..." %new_phage[0]
    #Retrieve phage-specific data from phagesdb
    phage_url = api_prefix + new_phage + api_suffix
    online_data_json = urllib.urlopen(phage_url)
    online_data_dict = json.loads(online_data_json.read())


    #Check to see if there is a flatfile stored on phagesdb for this phage
    print online_data_dict["qced_genbank_file"]
    
    if online_data_dict["qced_genbank_file"] is not None:
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

    else:
        print "Error: no flatfile found for phage %s." %new_phage[0]
        failed_tally += 1
        failed_list.append(new_phage[0])




if retrieved_tally > 0:
    print "The following %s phage(s) were successfully retrieved:" %retrieved_tally
    for element in retrieved_list:
        print element
    
else:
    print "No new flatfiles available."


if failed_tally > 0:
    print "The following %s phage(s) failed to be retrieved:" %failed_tally
    for element in failed_list:
        print element

else:
    print "No phages failed to be retrieved."












#Close script.
print "\n\n\nFlatfile retrieval script completed."
import_table_file.close()












