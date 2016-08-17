#!/usr/bin/env python
#Database Export Script
#Travis Mavrich
#University of Pittsburgh
#20160713

#After finishing updating the database, and ready to upload to webfactional...



#Import modules
import time, sys, os, getpass
import MySQLdb as mdb
import subprocess




#Get the command line parameters
try:
    database = sys.argv[1]
    main_dir = sys.argv[2]
    backup_dir = sys.argv[3]
except:
    print "Incorrect Parameters: ./export_database.py DATABASE MAIN_DIRECTORY BACKUP_DIRETCTORY"
    sys.exit(1)

#Set up MySQL parameters
mysqlhost = 'localhost'
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')


#Exits MySQL
def mdb_exit(message):
    print "\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`
    print "\nThe export script did not complete."
    print "\nExiting MySQL."
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print "\nExiting export script."
    sys.exit(1)




date = time.strftime("%Y%m%d")




#Verify connection to database
try:
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
    sys.exit(1)



#Change database version
try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT version FROM version")
    version_old = str(cur.fetchone()[0])
    print "Old database version: " + version_old
    version_new_int = int(version_old) + 1
    version_new = str(version_new_int)
    print "New database version: " + version_new
    statement = """UPDATE version SET version = %s;""" % version_new_int
    cur.execute(statement)
    cur.execute("COMMIT")
    print "Database version has been updated."
    
    
except:
    mdb_exit("\nError updating database version.\nNo changes have been made to the database.")



#Create a new version file
try:
    print "Creating version file..."
    versionfile="%s.version" % database
    versionfile_handle = open(main_dir + versionfile,'w')
    command_string = "echo %s" % version_new
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list,stdout=versionfile_handle)
    print "Version file has been created."

except:
    mdb_exit("\nError creating version file.")





#Export genome and gene data to file
#Filename formatting: DATE_DATABASE_VERSION_genes/genomes.csv
try:
    print "Exporting genome data..."
    filename1 = "%s_%s_v%s_genomes.csv" % (date,database,version_new)   
    statement1 = """SELECT phage.PhageID, phage.Name, phage.HostStrain, phage.Cluster, phage.status, phage.SequenceLength, phage.Accession, phage.DateLastModified FROM phage INTO OUTFILE '/tmp/%s' FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'""" % filename1
    cur.execute(statement1)


    print "Exporting gene data..."
    filename2 = "%s_%s_v%s_genes.csv" % (date,database,version_new)
    statement2 = """SELECT phage.PhageID, phage.Name, phage.HostStrain, phage.Cluster, phage.status, gene.GeneID, gene.Name, gene.Orientation, gene.Start, gene.Stop, gene.Notes, pham.name FROM gene JOIN phage on gene.PhageID = phage.PhageID JOIN pham on gene.GeneID = pham.GeneID INTO OUTFILE '/tmp/%s' FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'""" % filename2
    cur.execute(statement2)
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)
    print "Genome and gene data exported."
        
except:
    mdb_exit("\nError exporting genome or gene data to file.")

con.close()






#Now that the version has been updated, make a copy of the database in the MAIN directory, that will be uploaded to webfactional
print "Dumping new %s database to the Main directory..." % database
dumpfile1 = "%s.sql" % database
dumpfile1_handle = open(main_dir + dumpfile1,'w')
command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,database)
command_list = command_string.split(" ")
proc = subprocess.check_call(command_list,stdout=dumpfile1_handle)



#Also, create a backup of the update database in the BACKUP directory
print "Dumping copy of %s database to the Backup directory..." % database
dumpfile2 = "%s_v%s.sql" % (database,version_new)
dumpfile2_handle = open(backup_dir + dumpfile2,'w')
command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,database)
command_list = command_string.split(" ")
proc = subprocess.check_call(command_list,stdout=dumpfile2_handle)



#Close script.
print "\n\n\n\nExport script completed."




