"""Pipeline to upgrade or downgrade the schema of a Phamerator database."""

import argparse
import pkgutil
from pdm_utils.functions import phamerator

# The schema version wasn't formally tracked until schema version = 3.
MIN_VERSION = 3
MAX_VERSION = 7
VERSIONS = list(range(MIN_VERSION, MAX_VERSION + 1))
CHOICES = set(VERSIONS)
MODULE = "pdm_utils"
SQL_SCRIPT_DIR = "resources/sql_scripts/"


# TODO unittest.
def get_script_filename(dir, step):
    """Generates the name of the script conversion filename."""
    if dir == "upgrade":
        first = step - 1
        second = step
    elif dir == "downgrade":
        first = step + 1
        second = step
    else:
        first = step
        second = step
    filename = f"{dir}_{first}_to_{second}.sql"
    return filename


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for converting database."""
    CONVERT_HELP = ("Pipeline to upgrade or downgrade the "
                    "schema of a Phamerator database.")
    DATABASE_HELP = "Name of the MySQL database."
    SCHEMA_VERSION_HELP = ("Database schema version to which "
                           "the database should be converted.")
    parser = argparse.ArgumentParser(description=CONVERT_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-s", "--schema_version", type=int,
        choices=list(CHOICES), default=MAX_VERSION,
        help=SCHEMA_VERSION_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


# TODO not tested, but nearly identical function in import_genome.py tested.
def connect_to_db(database):
    """Connect to a MySQL database."""
    sql_handle, msg = phamerator.setup_sql_handle(database)
    if sql_handle is None:
        print(msg)
        sys.exit(1)
    else:
        return sql_handle


# TODO not unittested, but identical to function in export_db.
def retrieve_database_version(sql_handle):
    """Helper function that queries a SQL database
    for the database version and schema version

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :returns:
        database_versions_list(dictionary) is a dictionary
        of size 2 that contains values tied to keys
        "Version" and "SchemaVersion"
    """
    result = phamerator.retrieve_data(sql_handle,
                                      query="SELECT * FROM version")
    return result[0]


# TODO unittest.
def get_schema_version(dict):
    """Attempt to retrieve the schema version."""
    if "schema_version" in dict.keys():
        schema_version = dict["schema_version"]
    elif "SchemaVersion" in dict.keys():
        schema_version = dict["SchemaVersion"]
    else:
        schema_version = None

    if schema_version is None:
        print("Database does not have a schema version. "
              "Unable to proceed with schema conversion.")
        sys.exit(1)
    else:
        return schema_version


# TODO unittest.
def main(unparsed_args_list):
    """Run main conversion pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    sql_handle = connect_to_db(args.database)
    version_data = retrieve_database_version(sql_handle)
    target = args.schema_version
    actual = get_schema_version(version_data)
    steps, dir = get_conversion_direction(actual, target)

    # Iterate through list of versions and implement SQL files.
    if dir == "none":
        print("No schema conversion is needed.")
    else:
        stop_step = convert_schema(sql_handle, actual, dir, steps)
        if stop_step == target:
            print(f"The database schema conversion was successful.")
        else:
            print(f"The database schema conversion was unsuccessful. "
                  f"Unable to proceed past schema version {stop_step}.")


# TODO unittest.
def get_conversion_direction(actual, target):
    """Determine needed conversion direction and steps."""
    # Compare actual version to desired selection to get list of versions.
    if actual == target:
        steps = []
        dir = "none"
    elif actual < target:
        dir = "upgrade"
        steps = list(range(actual + 1, target + 1))
    else:
        steps = list(reversed(range(target, actual)))
        dir = "downgrade"
    return steps, dir

# TODO unittest.
def convert_schema(sql_handle, actual, dir, steps):
    """Iterate through conversion steps and convert database schema."""
    index = 0
    convert = True
    stop_step = actual
    while (index < len(steps) and convert == True):
        step = steps[index]
        script = get_script_filename(dir, step)
        print(f"{dir[:-1].capitalize()}ing to schema version {step}...")
        script_path = SQL_SCRIPT_DIR + script
        data = pkgutil.get_data(MODULE, script_path).decode()
        commands = data.split(";")
        commands = [i for i in commands if i != ""]
        commands = [i for i in commands if i != "\n"]
        # Try to convert the schema.
        result = sql_handle.execute_transaction(commands)
        if result == 1:
            convert = False
            print("Error encountered while executing MySQL statements.")
        if convert == False:
            print(f"Unable to {dir} schema to version {step}.")
        else:
            stop_step = step
        index += 1
    return stop_step
