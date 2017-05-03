#!/usr/bin/env python
#Database Export Script
#Travis Mavrich
#University of Pittsburgh
#20170503

#When it comes time to freeze a database for publication, run this script




#Import built-in modules
import time, sys, os, getpass
import subprocess

#Import third-party modules
try:
    import MySQLdb as mdb
    import paramiko #TODO may not need this module
except:
    print "\nUnable to import one or more of the following third-party modules: paramiko, MySQLdb."
    print "Install modules and try again.\n\n"
    sys.exit(1)




#Exits MySQL
def mdb_exit(message):
    print "\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`
    print "\nThe freeze database script did not complete."
    print "\nExiting MySQL."
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print "\nExiting freeze database script."
    sys.exit(1)


#Get the command line parameters
try:
    database = sys.argv[1]
    main_dir = sys.argv[2]
except:
    print "\n\n\
            This is a python script to create a frozen Phamerator database for publication.\n\
            It requires four arguments:\n\
            First argument: name of MySQL database that will be updated (e.g. 'Actino_Draft').\n\
            Second argument: directory path to the Main folder used to store the new database.\n\
    sys.exit(1)







#Expand home directory
home_dir = os.path.expanduser('~')


#Verify the main folder exists

#Expand the path if it references the home directory
if main_dir[0] == "~":
    main_dir = home_dir + main_dir[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
main_dir = os.path.abspath(main_dir)


if main_dir[-1] != "/":
    main_dir = main_dir + "/"

if os.path.isdir(main_dir) == False:
    print "\n\nInvalid main database directory path.\n\n"
    sys.exit(1)





#Set up MySQL parameters
mysqlhost = 'localhost'
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')



#Set up misc variables
date = time.strftime("%Y%m%d")











###Query MySQL database and retrieve information


#TODO not sure if the code block below is necessary for this script
#Verify connection to database
try:
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
    print "Exiting export script."
    sys.exit(1)




#Retrieve database version and change if requested
try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT version FROM version")
    version_current = str(cur.fetchone()[0])
    print "Old database version: " + version_current

    if version_change == "yes":


        version_new_int = int(version_current) + 1
        version_current = str(version_new_int)
        print "New database version: " + version_current
        statement = """UPDATE version SET version = %s;""" % version_new_int
        cur.execute(statement)
        cur.execute("COMMIT")
        print "Database version has been updated."

    else:
        print "Database version will not be updated."

except:
    mdb_exit("\nError retrieving database version.\nNo changes have been made to the database.")









###Create the new frozen db folder e.g. Actinobacteriophage_XYZ


#TODO variables below are not correct
versionfile="%s.version" % database
dumpfile1 = "%s.sql" % database




###Output database to new folder with new name


###Import new database into MySQL





###Drop all 'draft' databases




#Close script.
print "\n\n\n\Freeze database script completed."





















###Below is code from export database script



#Create a new version file
#It is okay if there is already a copy of this file present. Filenames like "Actino_Draft.sql" and "Actino_Draft.version" never change.
try:
    print "Creating version file..."
    versionfile_handle = open(main_dir + versionfile,'w')
    command_string = "echo %s" % version_current
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list,stdout=versionfile_handle)
    versionfile_handle.close()
    print "Version file has been created."

except:
    mdb_exit("\nError creating version file.")




#Now that the version has been updated, make a copy of the database in the MAIN directory, that will be uploaded to webfactional
#It is okay if there is already a copy of this file present. Filenames like "Actino_Draft.sql" and "Actino_Draft.version" never change.
print "Dumping new %s database to the Main directory..." % database

dumpfile1_handle = open(main_dir + dumpfile1,'w')
command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,database)
command_list = command_string.split(" ")
proc = subprocess.check_call(command_list,stdout=dumpfile1_handle)
dumpfile1_handle.close()


#Also, create a backup of the update database in the BACKUP directory
#The backup folder should contain files that serve as backup. Therefore, these files should not be overwritten by default, unlike the version and database files in the "Main" directory.
#First check if the path to the backup file already exists. If it does, this means that something may have gone wrong in the export process,
#or that this export has previously been completed.
dumpfile2 = "%s_v%s.sql" % (database,version_current)

if os.path.exists(backup_dir + dumpfile2) == True:
    print "\n\nThe backup database file already exists in the indicated backup directory."
    print "Backup file path: %s\n\n" %(backup_dir + dumpfile2)
    sys.exit(1)


print "Dumping copy of %s database to the Backup directory..." % database
dumpfile2_handle = open(backup_dir + dumpfile2,'w')
command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,database)
command_list = command_string.split(" ")
proc = subprocess.check_call(command_list,stdout=dumpfile2_handle)
dumpfile2_handle.close()



#Close MySQL connection
cur.close()
con.autocommit(True)
con.close()
