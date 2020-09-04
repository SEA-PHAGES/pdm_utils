"""Pipeline to update specific fields in a MySQL database."""
import argparse
import csv
import pathlib

from sqlalchemy import update

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import querying


# TODO unittest.
def main(unparsed_args):
    """Runs the complete update pipeline."""
    args = parse_args(unparsed_args[2:])

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]
    alchemist = AlchemyHandler(database=args.database,
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the update pipeline")

    if args.version is True:
        mysqldb.change_version(engine)
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
            status = update_field(alchemist, dict)

            if status == 1:
                processed += 1
                succeeded += 1
            else:
                processed += 1
                failed += 1

        engine.dispose()
        print("\nDone iterating through tickets.")
        if succeeded > 0:
            print(f"{succeeded} / {processed} tickets successfully handled.")
        if failed > 0:
            print(f"{failed} / {processed} tickets failed to be handled.")


def update_field(alchemist, update_ticket):
    """Attempts to update a field using information from an update_ticket.

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param update_ticket: Dictionary with instructions to update a field.
    :type update_ticket: dict
    """
    try:
        table_map = alchemist.mapper.classes[update_ticket["table"]]
    except:
        print(f"\nInvalid table '{update_ticket['table']}'")
        return 0

    table_obj = table_map.__table__

    try:
        field = table_obj.c[update_ticket["field"]]
    except:
        print(f"\nInvalid replacement field '{update_ticket['field']}' "
              f"for table '{table_obj.name}'")
        return 0

    try:
        key_field = table_obj.c[update_ticket["key_name"]]
    except:
        print(f"\nInvalid selection key '{update_ticket['key_name']}' "
              f"for table '{table_obj.name}'")
        return 0

    primary_keys = list(table_obj.primary_key.columns)
    if key_field not in primary_keys:
        print(f"\nInvalid selection key '{update_ticket['key_name']}' "
              f"for table '{table_obj.name}'")
        return 0

    key_value_clause = (key_field == update_ticket["key_value"])

    key_value_query = querying.build_count(alchemist.graph, key_field, where=key_value_clause)
    key_value_count = mysqldb_basic.scalar(alchemist.engine, key_value_query)

    if key_value_count != 1:
        print(f"\nInvalid selection value '{update_ticket['key_value']}' "
              f"for key '{key_field.name}' in table '{table_obj.name}'")
        return 0

    if update_ticket["value"] == "NULL":
        update_ticket["value"] = None

    statement = update(table_obj).where(key_value_clause).values(
                                    {field.name : update_ticket["value"]})
    alchemist.engine.execute(statement)
    return 1


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting updates."""
    update_help = "Pipeline to update specific fields in a MySQL database."
    database_help = "Name of the MySQL database."
    ticket_table_help = """
        Path to the CSV-formatted table containing
        instructions to process each update.
        Structure of update ticket table:
            1. Table name
            2. Column name to be updated
            3. New value that will be added
            4. Conditional column name
            5. Conditional column value
        """
    version_help = "Increment database version by 1."
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=update_help)
    parser.add_argument("database", type=str, help=database_help)
    parser.add_argument("-f", "--ticket_table", type=pathlib.Path,
                        help=ticket_table_help)
    parser.add_argument("-v", "--version", action="store_true",
                        default=False, help=version_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)
    args = parser.parse_args(unparsed_args_list)
    return args
