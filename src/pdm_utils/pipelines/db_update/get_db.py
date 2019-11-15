"""Pipeline to check for new versions of the database from the server."""

import argparse
import os
import pathlib
import subprocess
import sys
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils.constants import constants
from pdm_utils.functions import basic

# TODO not tested, but identical function in import_genome.py tested.
def setup_sql_handle(database):
    """Connect to a MySQL database."""
    sql_handle = mch.MySQLConnectionHandler()
    sql_handle.database = database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        print(f"No connection to the {database} database. "
            f"Valid credentials: {sql_handle.credential_status}. "
            f"Valid database: {sql_handle._database_status}")
        sys.exit(1)
    else:
        return sql_handle



# TODO not tested, but identical function in import_genome.py tested.
def set_path(path, kind=None, expect=True):
    """Confirm validity of path argument."""
    path = path.expanduser()
    path = path.resolve()
    result, msg = basic.verify_path2(path, kind=kind, expect=expect)
    if not result:
        print(msg)
        sys.exit(1)
    else:
        return path

# TODO in development
# TODO unittest.
def main(unparsed_args_list):
    """Run the get_db pipeline."""
    args = parse_args(unparsed_args_list)
    args.output_folder = set_path(args.output_folder, kind="dir", expect=True)

    if args.output_folder.exists() == False:
        print("The output ")
    # curl website > output_file
    # Version file
    version_filename = args.database + ".version"
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
    db_filename = args.database + ".sql"
    db_url = constants.DB_WEBSITE + db_filename
    db_filepath = pathlib.Path(args.output_folder, db_filename)
    if args.download == True:
        if db_filepath.exists() == True:
            print(f"The file {db_filename} already exists.")
        else:
            command_string2 = f"curl {db_url}"
            command_list2 = command_string2.split(" ")
            status = False
            with db_filepath.open("w") as db_handle:
                try:
                    print("Downloading sql file.")
                    subprocess.check_call(command_list2, stdout=db_handle)
                    status = True
                except:
                    print(f"Unable to download {db_filename} from server.")

    # Install new database
    # TODO this should first check if the database exists in MySQL.
    # If not, it should create the database.
    if args.install == True:
        result = basic.verify_path2(db_filepath, kind="file", expect=True)
        if result[0] == False:
            print("Unable to locate database file to install.")
            print(result[1])
        else:
            # TODO currently there really is no need to set up a sql_handle object.
            # But it will be used once the pipeline is improved.
            sql_handle = setup_sql_handle(args.database)
            command_string3 = (f"mysql -u {sql_handle.username} "
                               f"-p{sql_handle.password} {sql_handle.database}")
            command_list3 = command_string3.split(" ")
            with db_filepath.open("r") as db_handle:
                try:
                    print("Installing database...")
                    subprocess.check_call(command_list3, stdin=db_handle)
                    print("Installation complete.")
                except:
                    print(f"Unable to install {db_filename} in MySQL.")

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
