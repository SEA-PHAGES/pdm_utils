#!/usr/bin/env python
#Retrieve auto-annotated genomes from PECAAN
#University of Pittsburgh
#Travis Mavrich
#20170206
#The purpose of this script is to retrieve auto-annotated draft genomes from PECAAN, as an alternative to DNA Master

#Third-party libraries
import MySQLdb as mdb


#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib, urllib2, time




#Get the command line parameters
try:
    updateFileDir = sys.argv[1] #What is the directory into which the report should go
except:
    print "\n\n\
            This is a python script to retrieve auto-annotated genomes from PECAAN to import into Phamerator.\n\
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




#Create output directories
date = time.strftime("%Y%m%d")

output_folder = '%s_new_draft_genomes' % date
output_path = os.path.join(updateFileDir,output_folder)


try:

    os.mkdir(output_path)
except:
    print "\nUnable to create output folder: %s" %output_path
    sys.exit(1)

os.chdir(output_path)

#Retrieve list of unphamerated genomes
#Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
phage_file = '/home/cbowman/Documents/PhameratorDB_Management/Updates/temp/new_phages.txt'


phage_file_handle = open(phage_file,'r')
phage_file_reader = csv.reader(phage_file_handle,delimiter='\t')




    
#Open file to create import table with changes that need to be implemented
import_table_file = open(os.path.join(updateFileDir,output_folder,date + "_draft_genome_import_table.csv"), "w")
import_table_writer = csv.writer(import_table_file)






#Retrieve auto-annotated genomes from PECAAN
pecaan_prefix = 'https://discoverdev.kbrinsgd.org/phameratoroutput/phage/'
retrieved_tally = 0
failed_tally = 0
retrieved_list = []
failed_list = []



for new_phage in phage_file_reader:

    #PECAAN should be able to generate any phage that is listed on phagesdb
    #print new_phage
    print "Attempting to retrieve %s from PECAAN..." %new_phage[0]
    pecaan_link = pecaan_prefix + new_phage[0]
    pecaan_file = new_phage[0] + "_Draft.txt"
    print pecaan_link
    try:
        response = urllib2.urlopen(pecaan_link)
        pecaan_file_handle = open(pecaan_file,'w')
        pecaan_file_handle.write(response.read())
        response.close()
        pecaan_file_handle.close()
        
        
        #Create the new import ticket
        #0 = Import action
        #1 = New phageID
        #2 = HostStrain
        #3 = Cluster
        #4 = Status
        #5 = Gene description field
        #6 = Phage to replace
        import_table_writer.writerow(["add",new_phage[0],"retrieve","retrieve","draft","product","none"])

        retrieved_tally += 1
        retrieved_list.append(new_phage[0])

    except:
        print "Error: unable to retrieve %s draft genome." %new_phage[0]
        failed_tally += 1
        failed_list.append(new_phage[0])




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


#Close script.
print "\n\n\nDraft genome retrieval script completed."
import_table_file.close()


















