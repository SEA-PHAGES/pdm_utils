"""Pipeline to update specific fields in a MySQL database."""
import argparse
import sys
import csv
import pathlib
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from pdm_utils.classes.randomfieldupdatehandler import RandomFieldUpdateHandler
from pdm_utils.functions import basic, phamerator


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
def main(unparsed_args):
    """Runs the complete update pipeline."""
    args = parse_args(unparsed_args)

    # Establish the database connection using the MySQLConnectionHandler object
    mysql_handler = connect_to_db(args.database)

    if args.version == True:
        phamerator.change_version(mysql_handler)
        print("Database version updated.")

    if args.ticket_table is not None:
        update_table_path = basic.set_path(args.ticket_table,
                                kind="file", expect=True)

        # Iterate through the tickets and process them sequentially.
        list_of_update_tickets = []
        with update_table_path.open(mode='r') as f:
            file_reader = csv.DictReader(f)
            for dict in file_reader:
                list_of_update_tickets.append(dict)

        # Variables to be used for end summary
        processed = 0
        succeeded = 0
        failed = 0

        for dict in list_of_update_tickets:
            handler = RandomFieldUpdateHandler(mysql_handler.connection)
            handler.table = dict["table"]        # Which table will be updated?
            handler.field = dict["field"]       # Which field will be updated?
            handler.value = dict["value"]       # What value will be put in that field?
            handler.key_name = dict["key_name"]   # How will we know which row is the right one?
            handler.key_value = dict["key_value"]  # How will we know which row is the right one?
            handler.validate_ticket()   # Make sure all handler attributes are valid
            status = handler.execute_ticket()    # Do what was requested
            if status == 1:
                processed += 1
                succeeded += 1
            else:
                processed += 1
                failed += 1

        print("\nDone iterating through tickets.")
        if succeeded > 0:
            print(f"{succeeded} / {processed} tickets successfully handled.")
        if failed > 0:
            print(f"{failed} / {processed} tickets failed to be handled.")



# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""
    UPDATE_HELP = ("Pipeline to update specific fields in a "
                   "MySQL Phamerator database.")
    DATABASE_HELP = "Name of the MySQL database."
    TICKET_TABLE_HELP = """
            Path to the CSV-formatted table containing
            instructions to process each update.
            Structure of update ticket table:
                1. Table name
                2. Column name to be updated
                3. New value that will be added
                4. Conditional column name
                5. Conditional column value
            """
    VERSION_HELP = "Increment database version by 1."
    parser = argparse.ArgumentParser(description=UPDATE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-f", "--ticket_table", type=pathlib.Path,
        help=TICKET_TABLE_HELP)
    parser.add_argument("-v", "--version", action="store_true",
        default=False, help=VERSION_HELP)
    args = parser.parse_args(unparsed_args_list)
    return args
