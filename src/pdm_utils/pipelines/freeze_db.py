"""Pipeline to freeze a database for publication."""

import argparse
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb
from pdm_utils.functions import parsing
from pdm_utils.classes.filter import Filter

from pdm_utils.classes.alchemyhandler import AlchemyHandler


RESET_VERSION = "UPDATE version SET Version = 0"
DELETE_DRAFTS = "DELETE FROM phage WHERE Status = 'draft'"

# TODO unittest.
def main(unparsed_args_list):
    """Run main freeze database pipeline."""
    args = parse_args(unparsed_args_list)

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    alchemist = establish_database_connection(args.database)

    engine1 = alchemist.engine
    # engine1 = mysqldb.connect_to_db(args.database)
    mysqldb.check_schema_compatibility(engine1, "the freeze pipeline")


    # Get the number of draft genomes.
    query = "SELECT count(*) as count FROM phage WHERE Status != 'draft'"
    phage_count = alchemist.scalar(query)

    # result = engine1.execute(query).fetchall()
    # phage_count = result[0]["count"]

    # Create the new frozen database folder e.g. Actinobacteriophage_XYZ.

    if args.new_database_name is not None:
        new_database = args.new_database_name
    else:
        if args.prefix is not None:
            prefix = args.prefix
        else:
            prefix = get_prefix()
        new_database = f"{prefix}_{phage_count}"

    # Create the new database
    # Prevent overwriting of current database.
    if engine1.url.database != new_database:
        result = mysqldb.drop_create_db(engine1, new_database)
    else:
        print("Error: names of the reference and frozen databases are the same.")
        result == 1

    # Copy database.
    if result == 0:
        result = mysqldb.copy_db(engine1, new_database)
        if result == 0:
            print(f"Deleting 'draft' genomes...")
            engine2, msg = mysqldb.get_engine(
                                username=engine1.url.username,
                                password=engine1.url.password,
                                database=new_database)
            engine2.execute(DELETE_DRAFTS)
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


# TODO this may be moved elsewhere as a more generalized function.
# It may be best to add a sys.exit call.
def establish_database_connection(database: str):
    alchemist = AlchemyHandler(database=database)
    alchemist.connect()
    return alchemist


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected."""

    PIPELINE_HELP = ("Pipeline to prepare a database for publication.")
    DATABASE_HELP = "Name of the MySQL database."

    FILTERS_HELP = (
        "Genome selection option that indicates which genomes "
        " to remove or retain. "
        "Follow selection argument with formatted filter request: "
        "Table.Field=Value"
        )

    NEW_DATABASE_NAME_HELP = ("The new name of the frozen database")
    PREFIX_HELP = ("The prefix used in the new name of the frozen database")



    parser = argparse.ArgumentParser(description=PIPELINE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-f", "--filter", nargs="?",
                        type=parsing.parse_cmd_string, help=FILTERS_HELP)
    parser.add_argument("-n", "--new_database_name", type=str,
        help=NEW_DATABASE_NAME_HELP)
    parser.add_argument("-p", "--prefix", type=str,
        help=PREFIX_HELP)

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
