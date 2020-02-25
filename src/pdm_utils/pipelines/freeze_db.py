"""Pipeline to freeze a database for publication."""

import argparse
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb

# TODO unittest.
def main(unparsed_args_list):
    """Run main freeze database pipeline."""
    args = parse_args(unparsed_args_list)
    # engine1, msg = mysqldb.get_engine(database=args.database, echo=False)

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    engine1 = mysqldb.connect_to_db(args.database)
    mysqldb.check_schema_compatibility(engine1, "the freeze pipeline")


    # Get the number of draft genomes.
    query = "SELECT count(*) as count FROM phage WHERE Status != 'draft'"
    result = engine1.execute(query).fetchall()
    phage_count = result[0]["count"]

    # Create the new frozen database folder e.g. Actinobacteriophage_XYZ.
    prefix = get_prefix()
    new_database = f"{prefix}_{phage_count}"

    # Create the new database
    result = mysqldb.drop_create_db(engine1, new_database)

    # Copy database.
    if result == 0:
        result = mysqldb.copy_db(engine1, new_database)
        if result == 0:
            print(f"Deleting 'draft' genomes...")
            engine2, msg = mysqldb.get_engine(
                                username=engine1.url.username,
                                password=engine1.url.password,
                                database=new_database)

            statement3 = "DELETE FROM phage WHERE Status = 'draft'"
            engine2.execute(statement3)
            statement4 = "UPDATE version SET Version = 0"
            engine2.execute(statement4)
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

    PIPELINE_HELP = ("Pipeline to prepare a database for publication.")
    DATABASE_HELP = "Name of the MySQL database."
    parser = argparse.ArgumentParser(description=PIPELINE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)

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
