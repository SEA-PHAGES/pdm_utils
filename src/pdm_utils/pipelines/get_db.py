"""Pipeline to check for new versions of the database from the server."""

import argparse
import os
import pathlib
import subprocess
import sys
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils.constants import constants
from pdm_utils.functions import basic, phamerator



# TODO unittest.
def main(unparsed_args_list):
    """Run the get_db pipeline."""
    args = parse_args(unparsed_args_list)
    args.output_folder = basic.set_path(args.output_folder, kind="dir", expect=True)

    # curl website > output_file
    # Version file
    if args.filename is None:
        version_filename = args.database + ".version"
    else:
        version_filename = args.filename + ".version"
    version_url = constants.DB_WEBSITE + version_filename
    version_filepath = pathlib.Path(args.output_folder, version_filename)
    if args.download == True:
        if version_filepath.exists() == True:
            print(f"The file {version_filename} already exists.")
        else:
            command_string = f"curl {version_url}"
            command_list = command_string.split(" ")
            with version_filepath.open("w") as version_handle:
                try:
                    print("Downloading version file.")
                    subprocess.check_call(command_list, stdout=version_handle)
                except:
                    print(f"Unable to download {version_filename} from server.")

    # Database file
    if args.filename is None:
        db_filename = args.database + ".sql"
    else:
        db_filename = args.filename + ".sql"
    db_url = constants.DB_WEBSITE + db_filename
    db_filepath = pathlib.Path(args.output_folder, db_filename)
    if args.download == True:
        if db_filepath.exists() == True:
            print(f"The file {db_filename} already exists.")
        else:
            command_string2 = f"curl {db_url}"
            command_list2 = command_string2.split(" ")
            status = False
            # TODO if the file is not found on the server,
            # a file will still be created with text indicating
            # there was an error. So this step needs to catch that
            # error.
            with db_filepath.open("w") as db_handle:
                try:
                    print("Downloading sql file.")
                    subprocess.check_call(command_list2, stdout=db_handle)
                    status = True
                except:
                    print(f"Unable to download {db_filename} from server.")

    # Install new database
    if args.install == True:
        result1 = basic.verify_path2(db_filepath, kind="file", expect=True)
        if result1[0] == True:
            sql_handle = mch.MySQLConnectionHandler()
            sql_handle.open_connection()
            if sql_handle.credential_status:
                result2 = create_new_db(sql_handle, args.database)
                if result2 == 0:
                    sql_handle.database = args.database
                    sql_handle.open_connection()
                    if (sql_handle.credential_status == True and
                            sql_handle._database_status == True):
                        install_db(sql_handle, db_filepath)
                    else:
                        print(f"No connection to the {args.database} database due "
                              "to invalid credentials or database.")
                else:
                    print("Unable to create new, empty database.")
            else:
                print("Invalid MySQL credentials.")
        else:
            print("Unable to locate database file for installation.")
            print(result1[1])

    if args.remove == True:
        print("Removing downloaded data.")
        if version_filepath.exists() == True:
            os.remove(version_filepath)
        if db_filepath.exists() == True:
            os.remove(db_filepath)



# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting a new database."""

    UPDATE_DB_HELP = ("Pipeline to retrieve and install a new version of "
                   "a Phamerator MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the folder to create the folder for "
                          "downloading the database.")
    INSTALL_HELP = \
        ("Indicates if the new version should be installed.")
    DOWNLOAD_HELP = \
        ("Indicates if the new version should be downloaded.")
    FORCE_INSTALL_HELP = \
        ("Indicates if the downloaded version should overwrite the existing "
         "database or create a new database if it is not already present.")
    REMOVE_HELP = \
        ("Indicates if the downloaded file should be removed after installation.")
    ALL_HELP = \
        ("Indicates if the entire pipeline should run: "
         "download, install, and then remove the temp files.")
    FILENAME_HELP = \
        ("Indicates the name of the SQL file and verion file. "
         "By default, the filename is "
         "assumed to be the same as the database name. This option enables "
         "database to be created from a file of a different name.")

    parser = argparse.ArgumentParser(description=UPDATE_DB_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("output_folder", type=pathlib.Path,
        help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-d", "--download", action="store_true",
        default=False, help=DOWNLOAD_HELP)
    parser.add_argument("-i", "--install", action="store_true",
        default=False, help=INSTALL_HELP)
    parser.add_argument("-r", "--remove", action="store_true",
        default=False, help=REMOVE_HELP)
    parser.add_argument("-a", "--all_steps", action="store_true",
        default=False, help=ALL_HELP)
    parser.add_argument("-f", "--filename", type=str,
        help=FILENAME_HELP)


    # TODO implement this option.
    # parser.add_argument("-f", "--force_update", action="store_true",
    #     default=False, help=FORCE_INSTALL_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if args.all_steps == True:
        args.download = True
        args.install = True
        args.remove = True

    return args


# TODO unittest.
def create_new_db(sql_handle, database):
    """Creates a new, empty database."""
    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    query = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
             f"WHERE SCHEMA_NAME = '{database}'")
    result1 = sql_handle.execute_query(query)
    if len(result1) != 0:
        statement1 = [f"DROP DATABASE {database}"]
        result2 = sql_handle.execute_transaction(statement1)
    else:
        result2 = 0
    if result2 == 0:
        # Next, create the database within mysql.
        statement2 = [f"CREATE DATABASE {database}"]
        result2 = sql_handle.execute_transaction(statement2)
        sql_handle.close_connection()
    return result2

# TODO unittest.
def install_db(sql_handle, schema_filepath):
    """Install a MySQL file into the indicated database."""
    command_string = (f"mysql -u {sql_handle.username} "
                      f"-p{sql_handle.password} {sql_handle.database}")
    command_list = command_string.split(" ")
    with schema_filepath.open("r") as fh:
        try:
            print("Installing database...")
            subprocess.check_call(command_list, stdin=fh)
            print("Installation complete.")
        except:
            print(f"Unable to install {schema_filepath.name} in MySQL.")
