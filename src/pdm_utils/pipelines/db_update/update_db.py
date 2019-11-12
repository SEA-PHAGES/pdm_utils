"""Pipeline to check for new versions of the database from the server."""

import argparse
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
    # sql_handle = setup_sql_handle(args.database)

    download_folder = pathlib.Path(args.output_folder, "db_download")
    download_folder.mkdir()

    # curl website > output_file
    # Version file
    version_filename = args.database + ".version"
    version_url = constants.DB_WEBSITE + version_filename
    version_filepath = pathlib.Path(download_folder, version_url)
    command_string = f"curl {version_url}"
    command_list = command_string.split(" ")
    try:
        with version_filepath.open("w") as version_handle:
            proc = subprocess.check_call(command_list,stdout=version_handle)
    except:
        print(f"Unable to download {version_filename} from server.")

    # Database file
    db_filename = args.database + ".sql"
    db_url = constants.DB_WEBSITE + db_filename
    db_filepath = pathlib.Path(download_folder, db_filename)
    command_string2 = f"curl {db_url}"
    command_list2 = command_string2.split(" ")
    try:
        with db_filepath.open("w") as db_handle:
            proc = subprocess.check_call(command_list2,stdout=db_handle)
    except:
        print(f"Unable to download {db_filename} from server.")



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

    parser = argparse.ArgumentParser(description=UPDATE_DB_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path("/tmp/"), help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-d", "--download", action="store_true",
        default=False, help=DOWNLOAD_HELP)
    parser.add_argument("-i", "--install", action="store_true",
        default=False, help=INSTALL_HELP)
    parser.add_argument("-f", "--force_update", action="store_true",
        default=False, help=FORCE_INSTALL_HELP)
    parser.add_argument("-r", "--remove", action="store_true",
        default=False, help=REMOVE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    return args
