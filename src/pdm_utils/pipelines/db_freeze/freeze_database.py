#!/usr/bin/env python
#Freeze Database Script
#Travis Mavrich
#University of Pittsburgh
#20170503

#When it comes time to freeze a database for publication, run this script to
#save database with a new name and drop all draft genomes




#Import built-in modules
import sys, os, getpass
import subprocess

#Import third-party modules
try:
    import pymysql as pms
    # import MySQLdb as mdb
except:
    print("\nUnable to import one or more of the following third-party modules: pymysql.")
    print("Install modules and try again.\n\n")
    sys.exit(1)




#Exits MySQL
def mdb_exit(message):
    print(message)
    print("Error")
    # print("\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
    print("\nThe freeze database script did not complete.")
    print("\nExiting MySQL.")
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print("\nExiting freeze database script.")
    sys.exit(1)


#Allows user to select specific options
def select_option(message,valid_response_set):

    response_valid = False
    while response_valid == False:
        response = input(message)
        if response.isdigit():
            response = int(response)
        else:
            response = response.lower()

        if response in valid_response_set:
            response_valid = True
            if response == 'y':
                response  = 'yes'
            elif response == 'n':
                response  = 'no'
        else:
            print('Invalid response.')
    return response




def main(unparsed_args_list):
    #Get the command line parameters
    try:
        old_database = unparsed_args_list[2]
        main_dir = unparsed_args_list[3]
    except:
        print("\n\n\
                This is a python script to create a frozen Phamerator database for publication.\n\
                It requires two arguments:\n\
                First argument: name of MySQL database that will be frozen (e.g. 'Actino_Draft').\n\
                Second argument: directory path to store the new database.\n")
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



    #Set up MySQL parameters
    mysqlhost = 'localhost'
    username = getpass.getpass(prompt='mySQL username:')
    password = getpass.getpass(prompt='mySQL password:')




    #Determine the new database name prefix
    database_prefix_options = [\
        'none',\
        'Actinobacteriophage',\
        'Bacteriophage',\
        'Cyanobacteriophage',\
        'Custom']
    print('\n\nDefault database prefix options:')
    print('1: ' + database_prefix_options[1])
    print('2: ' + database_prefix_options[2])
    print('3: ' + database_prefix_options[3])
    print('4: ' + database_prefix_options[4])
    database_prefix_selection = select_option(\
        "\nWhich database prefix should be used? ", \
        set([1,2,3,4]))


    if database_prefix_selection == 1:
        database_prefix = database_prefix_options[1]
    elif database_prefix_selection == 2:
        database_prefix = database_prefix_options[2]
    elif database_prefix_selection == 3:
        database_prefix = database_prefix_options[3]
    elif database_prefix_selection == 4:
        database_prefix = input("Provide the custom database prefix: ")






    #Query MySQL database and retrieve information
    #Verify connection to database
    try:
        con = pms.connect("localhost", username, password, old_database)
        # con = mdb.connect(mysqlhost, username, password, old_database)
        con.autocommit(False)
        cur = con.cursor()
    except pms.err.Error as err:
        print("Error connecting to MySQL database")
        print("Error {}: {}".format(err.args[0], err.args[1]))
        sys.exit(1)


    #Retrieve the number of non-draft phages in the database
    try:
        cur.execute("START TRANSACTION")
        cur.execute("SELECT count(*) FROM phage WHERE Status != 'draft'")
        phage_count = str(cur.fetchone()[0])
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)
    except:
        mdb_exit("\nError retrieving database information.")
    con.close()



    #Create the new frozen database folder e.g. Actinobacteriophage_XYZ
    new_database = "%s_%s" %(database_prefix,phage_count)
    dumpfile1_name = "%s.sql" %new_database
    new_database_directory = os.path.join(main_dir,new_database)


    if os.path.isdir(new_database_directory) == True:
        print("\n\nInvalid new database directory path.\n\n")
        sys.exit(1)


    new_database_current_directory = os.path.join(new_database_directory,"Current")
    os.mkdir(new_database_directory)
    os.mkdir(new_database_current_directory)
    os.mkdir(os.path.join(new_database_directory,"Backup"))
    os.mkdir(os.path.join(new_database_directory,"Update_history"))







    #Now create the new MySQL database
    #Verify connection to database
    try:
        con = pms.connect("localhost", username, password)
        # con = mdb.connect(mysqlhost, username, password)
        cur = con.cursor()
    except pms.err.Error as err:
        print("Error connecting to MySQL database")
        print("Error {}: {}".format(err.args[0], err.args[1]))
        sys.exit(1)


    try:
        #If there is already a database with the same name, an error is thrown
        #and the script will exit.
        cur.execute("CREATE DATABASE %s" %new_database)
        cur.close()
    except:
        mdb_exit("\nError creating new database.")
    con.close()






    #Output database to new folder with new name
    print("Creating new %s database file..." % new_database)

    dumpfile1 = os.path.join(new_database_current_directory,dumpfile1_name)
    output_handle = open(dumpfile1,'w')
    command_string = "mysqldump -u %s -p%s --skip-comments %s" % (username,password,old_database)
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list,stdout=output_handle)
    output_handle.close()


    #Import new database file into MySQL
    print("Importing new %s database file into MySQL..." % new_database)

    input_handle = open(dumpfile1,'r')
    command_string = "mysql -u %s -p%s %s" % (username,password,new_database)
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list,stdin=input_handle)
    input_handle.close()






    #Now drop all 'draft' databases
    #Verify connection to database
    try:
        con = pms.connect("localhost", username, password, new_database)
        # con = mdb.connect(mysqlhost, username, password, new_database)
        con.autocommit(False)
        cur = con.cursor()
    except pms.err.Error as err:
        print("Error connecting to MySQL database")
        print("Error {}: {}".format(err.args[0], err.args[1]))
        sys.exit(1)


    #Now delete all draft genomes
    try:
        cur.execute("START TRANSACTION")
        cur.execute("DELETE FROM phage WHERE Status = 'draft'")

        #Set version to 0 to indicate that the database needs re-phamerated
        cur.execute("UPDATE version SET Version = 0")

        #Close MySQL connection
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)

    except:
        mdb_exit("\nError removing draft genomes.\nNo changes have been made to the database.")
    con.close()


    #Delete the new database file created, since it still contains draft phagesdb
    #and incorrect version and pham information
    input("Press ENTER to proceed")
    os.remove(dumpfile1)



    #Close script.
    print("\n\n\nFreeze database script completed.")

if __name__ == "__main__":
    main(sys.argv.insert(0, "empty"))
