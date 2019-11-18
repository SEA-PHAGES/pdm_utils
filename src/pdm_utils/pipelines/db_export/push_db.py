"""Pipeline to push data to the server."""
import argparse
import os
import pathlib
import subprocess
import sys
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import server



# TODO unittest.
def main(unparsed_args_list):
    """Run the get_db pipeline."""
    args = parse_args(unparsed_args_list)
    args.input_folder = set_path(args.input_folder, kind="dir", expect=True)
    sql_handle = setup_sql_handle(args.database)





# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for uploading to the server."""
    PUSH_DB_HELP = ("Pipeline to retrieve and install a new version of "
                   "a Phamerator MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    INPUT_FOLDER_HELP = ("Path to the folder containing files for upload.")
    parser = argparse.ArgumentParser(description=PUSH_DB_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("input_folder", type=pathlib.Path,
        help=INPUT_FOLDER_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


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

###
