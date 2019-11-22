"""Pipeline to freeze a database for publication."""

import argparse
import sys
import os
import pathlib
import subprocess
from pdm_utils.functions import basic
from pdm_utils.functions import phamerator
from pdm_utils.classes import mysqlconnectionhandler as mch


# TODO not tested, but nearly identical function in import_genome.py tested.
def connect_to_db(database):
    """Connect to a MySQL database."""
    sql_handle, msg = phamerator.setup_sql_handle(database)
    if sql_handle is None:
        print(msg)
        sys.exit(1)
    else:
        return sql_handle


# TODO unittest.
def main(unparsed_args_list):
    """Run main freeze database pipeline."""
    args = parse_args(unparsed_args_list)
    args.output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    sql_handle = connect_to_db(args.database)

    # Get the number of draft genomes.
    query = "SELECT count(*) FROM phage WHERE Status != 'draft'"
    result = sql_handle.execute_query(query)
    phage_count = result[0]["count(*)"]

    # Create the new frozen database folder e.g. Actinobacteriophage_XYZ.
    prefix = get_prefix()
    new_database = f"{prefix}_{phage_count}"

    new_dir = pathlib.Path(args.output_folder, new_database)
    new_dir = basic.set_path(new_dir, kind="dir", expect=False)
    new_current = pathlib.Path(new_dir, "Current")
    new_backup = pathlib.Path(new_dir, "Backup")
    new_history = pathlib.Path(new_dir, "Update_history")
    new_dir.mkdir()
    new_current.mkdir()
    new_backup.mkdir()
    new_history.mkdir()

    # Create the new database
    statement = f"CREATE DATABASE {new_database}"
    result = sql_handle.execute_transaction([statement])
    if result == 0:

        #Output database to new folder with new name
        print(f"Creating new {new_database} database file...")
        filename = f"temp_db.sql"
        filepath = pathlib.Path(new_dir, filename)
        command_string1 = (f"mysqldump -u {sql_handle.username} "
                           f"-p{sql_handle.password} "
                           f"--skip-comments {sql_handle.database}")
        command_list1 = command_string1.split(" ")
        with filepath.open("w") as fh:
            subprocess.check_call(command_list1,stdout=fh)

        # Import new database file into MySQL
        print(f"Importing new {new_database} database file into MySQL...")
        command_string2 = (f"mysql -u {sql_handle.username} "
                          f"-p{sql_handle.password} {new_database}")
        command_list2 = command_string2.split(" ")
        with filepath.open("r") as fh:
            subprocess.check_call(command_list2, stdin=fh)

        # Drop all 'draft' databases and reset version.
        print(f"Deleting 'draft' genomes...")

        # TODO it would be better to instantiate a new handler object.
        # Pass previously verified username and password,
        # but use new database.
        statement2 = f"USE {new_database}"
        sql_handle.execute_transaction([statement2])
        statement3 = "DELETE FROM phage WHERE Status = 'draft'"
        sql_handle.execute_transaction([statement3])
        statement4 = "UPDATE version SET Version = 0"
        sql_handle.execute_transaction([statement4])
        os.remove(filepath)
    else:
        print(f"Error creating new database: {new_database}.")
    print("Freeze database script completed.")


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected."""

    PIPELINE_HELP = ("Pipeline to prepare a database for publication.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the directory for storing files.")
    parser = argparse.ArgumentParser(description=PIPELINE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("output_folder", type=pathlib.Path,
        help=OUTPUT_FOLDER_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args



# TODO unittest.
def get_prefix():
    """Allow user to select appropriate prefix for the new database."""
    OPTIONS = ["none",
               "Actinobacteriophage",
               "Bacteriophage",
               "Cyanobacteriophage",
               "Custom"]
    print("\n\nDefault database prefix options:"
        f"\n1: {OPTIONS[1]}"
        f"\n2: {OPTIONS[2]}"
        f"\n3: {OPTIONS[3]}"
        f"\n4: {OPTIONS[4]}")
    prompt = "Which database prefix should be used? "
    choice = basic.select_option(prompt, {1,2,3,4})
    if choice == 1:
        prefix = OPTIONS[1]
    elif choice == 2:
        prefix = OPTIONS[2]
    elif choice == 3:
        prefix = OPTIONS[3]
    elif choice == 4:
        prefix = input("Provide the custom database prefix: ")
    return prefix
