#!/usr/bin/env python
# Database Export Script
# Travis Mavrich
# University of Pittsburgh
# 20160713
# Updated by Christian Gauthier 20181115 - 20190114;
# paramiko version was old and had been deprecated, resulting in inability
# to use SFTP to export SQL and version files to phamerator server.
# Paramiko version was upgraded, then some minor changes were made to
# the script to comply with new paramiko standards.
# Updated by Travis 20191010; paramiko functions moved to server module.



#Import built-in modules
import time, sys, os
import subprocess

#Import third-party modules
try:
    import pymysql as pms
except:
    print("\nUnable to import one or more of the following third-party modules: pymysql.")
    print("Install modules and try again.\n\n")
    sys.exit(1)

from pdm_utils.functions import server
from pdm_utils.functions import basic

def main(unparsed_args_list):
    #Get the command line parameters
    try:
        database = unparsed_args_list[2]
        main_dir = unparsed_args_list[3]
        backup_dir = unparsed_args_list[4]
        mysql_query_final_dir = unparsed_args_list[5]
    except:
        print("\n\n\
                This is a python script to export Phamerator databases.\n\
                It requires four arguments:\n\
                First argument: name of MySQL database that will be updated (e.g. 'Actino_Draft').\n\
                Second argument: directory path to the Main folder used to store the new database.\n\
                Third argument: directory path to the Backup folder used to store frozen backup versions.\n\
                Fourth argument: directory path to the folder used to stored gene and genome data queried from the new database.\n")
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
        print("\n\nInvalid main database directory path.\n\n")
        sys.exit(1)



    #Verify the backup folder exists

    #Expand the path if it references the home directory
    if backup_dir[0] == "~":
        backup_dir = home_dir + backup_dir[1:]

    #Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
    backup_dir = os.path.abspath(backup_dir)



    if backup_dir[-1] != "/":
        backup_dir = backup_dir + "/"

    if os.path.isdir(backup_dir) == False:
        print("\n\nInvalid backup database directory path.\n\n")
        sys.exit(1)



    #Verify the query folder exists

    #Expand the path if it references the home directory
    if mysql_query_final_dir[0] == "~":
        mysql_query_final_dir = home_dir + mysql_query_final_dir[1:]

    #Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
    mysql_query_final_dir = os.path.abspath(mysql_query_final_dir)


    if mysql_query_final_dir[-1] != "/":
        mysql_query_final_dir = mysql_query_final_dir + "/"

    if os.path.isdir(mysql_query_final_dir) == False:
        print("\n\nInvalid mysql query directory path.\n\n")
        sys.exit(1)







    #MySQL has changed the way it outputs queries. By default, query files are stored in the directory below.
    #I am unable to figure out how to change this to a custom directory. So now the script outputs queries to files in this default directory, and the files get copied to a custom directory.
    mysql_query_default_dir = '/var/lib/mysql-files/'


    #Exits MySQL
    def mdb_exit(message):
        # print("\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
        print("\nError:")
        print("\nThe export script did not complete.")
        print("\nExiting MySQL.")
        cur.execute("ROLLBACK")
        cur.execute("SET autocommit = 1")
        cur.close()
        con.close()
        print("\nExiting export script.")
        sys.exit(1)



    #Set up misc variables
    date = time.strftime("%Y%m%d")
    versionfile="%s.version" % database
    dumpfile1 = "%s.sql" % database








    #See if the user wants to export the database from MySQL
    dump_database = "no"
    dump_database_valid = False
    while dump_database_valid == False:
        dump_database = input("\nDo you want to export the database from MySQL? ")

        if (dump_database.lower() == "yes" or dump_database.lower() == "y"):
            dump_database = "yes"
            dump_database_valid = True

        elif (dump_database.lower() == "no" or dump_database.lower() == "n"):
            dump_database = "no"
            dump_database_valid = True

        else:
            print("Invalid response.")




    if dump_database == "yes":


        #Set up MySQL parameters
        username, password = basic.get_user_pwd(user_prompt="MySQL username: ",
                                pwd_prompt="MySQL password: ")

        #Verify connection to database
        try:
            con = pms.connect("localhost", username, password, database)
            con.autocommit(False)
            cur = con.cursor()
        except pms.err.Error as err:
            print("Error connecting to MySQL database")
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)



        #Allow user to control whether the database version number is incremented or not.
        #This option allows this script to be used to re-export databases, if needed, when no changes to the database have been made (and thus no need to update the version).
        version_change = "no"
        version_change_valid = False
        while version_change_valid == False:
            version_change = input("\nDo you want to increment the database version number? ")

            if (version_change.lower() == "yes" or version_change.lower() == "y"):
                version_change = "yes"
                version_change_valid = True

            elif (version_change.lower() == "no" or version_change.lower() == "n"):
                version_change = "no"
                version_change_valid = True

            else:
                print("Invalid response.")


        #Retrieve database version and change if requested
        try:
            cur.execute("START TRANSACTION")
            cur.execute("SELECT Version FROM version")
            version_current = str(cur.fetchone()[0])
            print("Old database version: " + version_current)

            if version_change == "yes":


                version_new_int = int(version_current) + 1
                version_current = str(version_new_int)
                print("New database version: " + version_current)
                statement = """UPDATE version SET Version = %s;""" % version_new_int
                cur.execute(statement)
                cur.execute("COMMIT")
                print("Database version has been updated.")

            else:
                print("Database version will not be updated.")

        except:
            mdb_exit("\nError retrieving database version.\nNo changes have been made to the database.")



        #Create a new version file
        #It is okay if there is already a copy of this file present. Filenames like "Actino_Draft.sql" and "Actino_Draft.version" never change.
        try:
            print("Creating version file...")
            versionfile_handle = open(main_dir + versionfile,'w')
            command_string = "echo %s" % version_current
            command_list = command_string.split(" ")
            proc = subprocess.check_call(command_list,stdout=versionfile_handle)
            versionfile_handle.close()
            print("Version file has been created.")

        except:
            mdb_exit("\nError creating version file.")




        #Now that the version has been updated, make a copy of the database in the MAIN directory, that will be uploaded to webfactional
        #It is okay if there is already a copy of this file present. Filenames like "Actino_Draft.sql" and "Actino_Draft.version" never change.
        print("Dumping new %s database to the Main directory..." % database)

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
            print("\n\nThe backup database file already exists in the indicated backup directory.")
            print("Backup file path: %s\n\n" %(backup_dir + dumpfile2))
            sys.exit(1)


        print("Dumping copy of %s database to the Backup directory..." % database)
        dumpfile2_handle = open(backup_dir + dumpfile2,'w')
        command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,database)
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list,stdout=dumpfile2_handle)
        dumpfile2_handle.close()



        #Close MySQL connection
        cur.close()
        con.autocommit(True)
        con.close()







    #See if the user wants to create genome and gene data files
    print("\nData export option disabled.")
    query_database = "no"
    query_database_valid = True
    # query_database_valid = False
    while query_database_valid == False:
        query_database = input("\nDo you want to query the database for gene and genome data? ")

        if (query_database.lower() == "yes" or query_database.lower() == "y"):
            query_database = "yes"
            query_database_valid = True

        elif (query_database.lower() == "no" or query_database.lower() == "n"):
            query_database = "no"
            query_database_valid = True

        else:
            print("Invalid response.")






    #Code is structured such that each section is independent, so even though MySQL username and password
    # were provided for dumping the database, it needs to be provided again.
    if query_database == "yes":




        #Set up MySQL parameters
        username, password = basic.get_user_pwd(
                                user_prompt="MySQL username: ",
                                pwd_prompt="MySQL password: ")



        #Verify connection to database
        try:
            con = pms.connect("localhost", username, password, database)
            con.autocommit(False)
            cur = con.cursor()
        except pms.err.Error as err:
            print("Error connecting to MySQL database")
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)


        #If connection was successfull, retrieve the database version
        try:
            cur.execute("START TRANSACTION")
            cur.execute("SELECT Version FROM version")
            version_export = str(cur.fetchone()[0])
            cur.execute("COMMIT")
        except:
            mdb_exit("\nError retrieving database version.\nUnable to query the database for gene and genome data.")



        #Verify the query output folder exists
        #Expand the path if it references the home directory
        query_folder = os.path.join(mysql_query_final_dir,"%s_%s_v%s" %(date,database,version_export))

        if os.path.isdir(query_folder) == True:
            cur.close()
            con.close()
            print("\nError creating new mysql query output folder in specified output directory. Folder already exists.")
            sys.exit(1)



        #Export genome and gene data to file
        #Filename formatting: DATE_DATABASE_VERSION_genes/genomes.csv
        os.mkdir(query_folder)
        try:
            print("Exporting genome data...")
            filename1 = "%s_%s_v%s_genomes.csv" % (date,database,version_export)
            statement1 = """SELECT phage.PhageID, phage.Name, phage.HostStrain, phage.Cluster, phage.Cluster2, phage.Subcluster2,
                        phage.Status, phage.SequenceLength, phage.Accession, phage.DateLastModified, phage.AnnotationAuthor
                        FROM phage INTO OUTFILE '%s/%s' FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'""" \
                        % (mysql_query_default_dir,filename1)
            cur.execute(statement1)

            command_string = "sudo cp %s/%s %s" % (mysql_query_default_dir,filename1,query_folder)
            command_list = command_string.split(" ")
            proc = subprocess.check_call(command_list)



            print("Exporting gene data...")
            filename2 = "%s_%s_v%s_genes.csv" % (date,database,version_export)
            statement2 = """SELECT phage.PhageID, phage.Name, phage.HostStrain, phage.Cluster, phage.Cluster2, phage.Subcluster2, \
                        phage.Status, gene.GeneID, gene.Name, gene.Orientation, gene.Start, gene.Stop, gene.Notes, pham.Name \
                        FROM gene JOIN phage on gene.PhageID = phage.PhageID JOIN pham on gene.GeneID = pham.GeneID INTO OUTFILE '%s/%s' \
                        FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'""" % (mysql_query_default_dir,filename2)
            cur.execute(statement2)
            cur.execute("COMMIT")
            cur.close()
            con.autocommit(True)
            command_string = "sudo cp %s/%s %s" % (mysql_query_default_dir,filename2,query_folder)
            command_list = command_string.split(" ")
            proc = subprocess.check_call(command_list)

            print("Genome and gene data exported.")

        except:
            mdb_exit("\nError exporting genome or gene data to file.")



        #MySQL stores the gene and genome data query files in a specific directory.
        #After the script copies these two files from the original MySQL directory to the user-requested directory, delete the original files in the MySQL directory.
        command_string = "sudo rm %s/%s" % (mysql_query_default_dir,filename1)
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list)

        command_string = "sudo rm %s/%s" % (mysql_query_default_dir,filename2)
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list)



        con.close()





    #Allow user to control whether the new database should be uploaded to the server.
    server_upload = "no"
    server_upload_valid = False
    while server_upload_valid == False:
        server_upload = input("\nDo you want to upload the database to the server? ")

        if (server_upload.lower() == "yes" or server_upload.lower() == "y"):
            server_upload = "yes"
            server_upload_valid = True

        elif (server_upload.lower() == "no" or server_upload.lower() == "n"):
            server_upload = "no"
            server_upload_valid = True

        else:
            print("Invalid response.")



    if server_upload == "yes":

        server.set_log_file("/tmp/paramiko.log")

        #Verify the path to both files exists
        if os.path.exists(main_dir + versionfile) == False:
            print("\n\nProblem locating version file.")
            print("Version file: %s\n\n" %(main_dir + versionfile))
            sys.exit(1)

        if os.path.exists(main_dir + dumpfile1) == False:
            print("\n\nProblem locating database file.")
            print("Version file: %s\n\n" %(main_dir + dumpfile1))
            sys.exit(1)



        #Set up paramiko parameters
        host = 'phamerator.webfactional.com'
        transport = server.get_transport(host)
        if transport is None:
            print("Exiting export script.")
            sys.exit(1)

        sftp = server.setup_sftp_conn(transport)
        if sftp is None:
            print("Exiting export script.")
            sys.exit(1)

        remote_dir = "/home/phamerator/webapps/htdocs/databases_Hatfull/"

        #First upload the version file
        print("Uploading the version file...")
        local_version = main_dir + versionfile
        remote_version = remote_dir + versionfile
        result1 = server.upload_file(sftp, local_version, remote_version)
        if result1:
            print("Version file successfully uploaded.")
        else:
            print("Unable to upload version file.")

        #Second upload the sql file
        print("Uploading the sql database...")
        local_dumpfile = main_dir + dumpfile1
        remote_dumpfile = remote_dir + dumpfile1
        result2 = server.upload_file(sftp, local_dumpfile, remote_dumpfile)
        if result2:
            print("Database successfully uploaded.")
        else:
            print("Unable to upload database file.")

        # Exit if there was a problem uploading either file.
        if (result1 == False or result2 == False):
            print("\nThe export script did not complete.")
            print("\nExiting export script.")
            sys.exit(1)

        #Close the connections
        sftp.close()
        transport.close()

    #Close script.
    print("\n\n\n\nExport script completed.")

if __name__ == "__main__":
    main(sys.argv.insert(0, "empty"))
