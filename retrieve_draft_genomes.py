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
import json, urllib, urllib2




#Get the command line parameters
try:
    database = sys.argv[1] #What Phamerator database should be compared to phagesdb?
    updateFileDir = sys.argv[2] #What is the directory into which the report should go
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








#Retrieve list of unphamerated genomes
#Retrieved file should be tab-delimited text file, each row is a newly sequenced phage
phage_file = '/home/cbowman/Documents/Phamerator/Updates/temp/new_phages.txt'

phage_file_handle = open(phage_file,'r')
phage_file_reader = csv.reader(phage_file_handle,sep='\t')







#Open file to create import table with changes that need to be implemented
import_table_file = open(os.path.join(updateFileDir,output_folder,date + "_draft_genome_import_table.csv"), "w")
import_table_writer = csv.writer(import_table_file)






#Retrieve auto-annotated genomes from PECAAN
pecaan_prefix = 'https://discoverdev.kbrinsgd.org/phameratoroutput/phage/'

for new_phage in phage_file_reader:

    #PECAAN should be able to generate any phage that is listed on phagesdb
    pecaan_link = pecaan_prefix + '/' + new_phage
    pecaan_file = new_phage + "_Draft.txt'
    print pecaan_link
    try:
        urllib2.urlretrieve(pecaan_link,pecaan_file)

        #Create the new import ticket
        #0 = Import action
        #1 = New phageID
        #2 = HostStrain
        #3 = Cluster
        #4 = Status
        #5 = Gene description field
        #6 = Phage to replace
        import_table_writer.writerow(["add",new_phage,"retrieve","retrieve","draft","product","none"])



    except:
        print "Error: unable to retrieve %s draft genome." %new_phage










#Close script.
print "All new sequences downloaded.")  
import_table_file.close()


















