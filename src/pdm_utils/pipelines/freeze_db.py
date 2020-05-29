"""Pipeline to freeze a database."""

import argparse
import pathlib
import sys

from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb, mysqldb_basic
from pdm_utils.functions import parsing
from pdm_utils.functions import querying
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter

RESET_VERSION = "UPDATE version SET Version = 0"
TARGET_TABLE = "phage"

# TODO unittest.
def main(unparsed_args_list):
    """Run main freeze database pipeline."""
    args = parse_args(unparsed_args_list)
    ref_database = args.database
    reset = args.reset
    new_database = args.new_database_name
    prefix = args.prefix

    # Filters input: phage.Status=draft AND phage.HostGenus=Mycobacterium
    # Args structure: [['phage.Status=draft'], ['phage.HostGenus=Mycobacterium']]
    filters = args.filters

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    alchemist1 = AlchemyHandler(database=ref_database,
                                username=mysql_creds["user"],
                                password=mysql_creds["password"])
    alchemist1.connect(pipeline=True)
    engine1 = alchemist1.engine
    mysqldb.check_schema_compatibility(engine1, "the freeze pipeline")

    # Get SQLAlchemy metadata Table object
    # table_obj.primary_key.columns is a
    # SQLAlchemy ColumnCollection iterable object
    # Set primary key = 'phage.PhageID'
    alchemist1.build_metadata()
    table = querying.get_table(alchemist1.metadata, TARGET_TABLE)
    for column in table.primary_key.columns:
        primary_key = column

    # Create filter object and then add command line filter strings
    db_filter = Filter(alchemist=alchemist1, key=primary_key)
    db_filter.values = []

    # Attempt to add filters and exit if needed.
    add_filters(db_filter, filters)

    # Performs the query
    db_filter.update()

    # db_filter.values now contains list of PhageIDs that pass the filters.
    # Get the number of genomes that will be retained and build the
    # MYSQL DELETE statement.
    keep_set = set(db_filter.values)
    delete_stmt = construct_delete_stmt(TARGET_TABLE, primary_key, keep_set)
    count_query = construct_count_query(TARGET_TABLE, primary_key, keep_set)
    phage_count = mysqldb_basic.scalar(alchemist1.engine, count_query)

    # Determine the name of the new database.
    if new_database is None:
        if prefix is None:
            prefix = get_prefix()
        new_database = f"{prefix}_{phage_count}"

    # Create the new database, but prevent overwriting of current database.
    if engine1.url.database != new_database:
        result = mysqldb_basic.drop_create_db(engine1, new_database)
    else:
        print("Error: names of the reference and frozen databases are the same.")
        print("No database will be created.")
        result = 1

    # Copy database.
    if result == 0:
        print(f"Reference database: {ref_database}")
        print(f"New database: {new_database}")
        result = mysqldb_basic.copy_db(engine1, new_database)
        if result == 0:
            print(f"Deleting genomes...")
            alchemist2 = AlchemyHandler(database=new_database,
                                        username=engine1.url.username,
                                        password=engine1.url.password)
            alchemist2.connect(pipeline=True)
            engine2 = alchemist2.engine
            engine2.execute(delete_stmt)
            if reset:
                engine2.execute(RESET_VERSION)

            # Close up all connections in the connection pool.
            engine2.dispose()
        else:
            print("Unable to copy the database.")
        # Close up all connections in the connection pool.
        engine1.dispose()
    else:
        print(f"Error creating new database: {new_database}.")
    print("Freeze database script completed.")

# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected."""

    freeze_help = "Pipeline to prepare a database for publication."
    database_help = "Name of the MySQL database."
    filters_help = (
        "Indicates which genomes to retain in the new database, "
        "with each conditional formatted as 'table.Field=value'.")
    new_database_name_help = "The new name of the frozen database"
    prefix_help = "The prefix used in the new name of the frozen database"
    reset_help = "Reset version to 0 in new database."
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=freeze_help)
    parser.add_argument("database", type=str, help=database_help)
    parser.add_argument("-f", "--filters", nargs="?",
                        type=parsing.parse_cmd_string, help=filters_help,
                        default=[])
    parser.add_argument("-n", "--new_database_name", type=str,
        help=new_database_name_help)
    parser.add_argument("-p", "--prefix", type=str,
        help=prefix_help)
    parser.add_argument("-r", "--reset", action="store_true",
        default=False, help=reset_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)

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

# TODO test.
def construct_set_string(phage_id_set):
    """Convert set of phage_ids to string formatted for MySQL.

    e.g. set: {'Trixie', 'L5', 'D29'}
    returns: "('Trixie', 'L5', 'D29')""
    """
    string = "('" + "', '".join(phage_id_set) + "')"
    return string

# TODO test.
def construct_count_query(table, primary_key, phage_id_set):
    """Construct SQL query to determine count."""
    phage_id_string = construct_set_string(phage_id_set)
    query = (f"SELECT count(*) as count FROM {table} "
             f"WHERE {primary_key} IN {phage_id_string}")
    return query

# TODO test.
def construct_delete_stmt(table, primary_key, phage_id_set):
    """Construct SQL query to determine count."""
    phage_id_string = construct_set_string(phage_id_set)
    statement = (f"DELETE FROM {table} "
                 f"WHERE {primary_key} NOT IN {phage_id_string}")
    return statement

# TODO test.
def add_filters(filter_obj, filters):
    """Add filters from command line to filter object."""
    errors = 0
    for or_filters in filters:
        for filter in or_filters:
            # Catch the error if it is an invalid table.column
            try:
                filter_obj.add(filter)
            except:
                print(f"Invalid filter: {filter}")
                errors += 1
    if errors > 0:
        print("Unable to create new database.")
        sys.exit(1)
