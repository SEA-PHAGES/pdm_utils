# Import standard python modules
import argparse
import sys
import csv
import pathlib
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from pdm_utils.classes.randomfieldupdatehandler import RandomFieldUpdateHandler
from pdm_utils.functions import basic


# TODO not tested, but identical function in import_genome.py tested.
def setup_sql_handle(database):
    """Connect to a MySQL database."""
    sql_handle = MySQLConnectionHandler()
    sql_handle.database = database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        logger.info(f"No connection to the {database} database. "
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

# TODO unittest.
def main(unparsed_args):
    """Runs the complete update pipeline."""
    # Set up argparse
    script_description = """
    This is a script intended to handle specific, single-field updates to tables
    in a Phamerator database.
    """
    parser = argparse.ArgumentParser(description=script_description)
    parser.add_argument("database_name", metavar="db", type=str, nargs=1,
                        help="name of the Phamerator database to be updated")
    parser.add_argument("ticket_file", metavar="tf", type=pathlib.Path, nargs=1,
                        help="path to an update ticket")


    # Parse command line arguments.
    args = parser.parse_args(unparsed_args)
    database = args.database_name[0].split("=")[-1]
    # update_table = args.ticket_file[0].split("=")[-1]
    update_table = args.ticket_file[0]
    update_table_path = set_path(update_table, kind="file", expect=True)

    # Establish the database connection using the MySQLConnectionHandler object
    mysql_handler = setup_sql_handle(database)

    # Iterate through the tickets and process them sequentially.
    list_of_update_tickets = []
    with update_table_path.open(mode='r') as f:
        file_reader = csv.DictReader(f)
        for dict in file_reader:
            list_of_update_tickets.append(dict)

    # Variables to be used for end summary
    tickets_processed = 0
    tickets_succeeded = 0
    tickets_failed = 0

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
            tickets_processed += 1
            tickets_succeeded += 1
        else:
            tickets_processed += 1
            tickets_failed += 1

    print("\nDone iterating through tickets.")
    print(f"{tickets_succeeded} / {tickets_processed} "
          "tickets successfully handled.")
    if tickets_failed > 0:
        print(f"{tickets_failed} / {tickets_processed} "
              "tickets failed to be handled.")


if __name__ == "__main__":
    main(sys.argv)
