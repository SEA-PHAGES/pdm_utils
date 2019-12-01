"""Pipeline to upgrade or downgrade the schema of a Phamerator database."""

import argparse
import pkgutil
from pdm_utils.functions import phamerator


# The schema version wasn't formally tracked until schema version = 3.
MIN_VERSION = 3
MAX_VERSION = 7
VERSIONS = list(range(MIN_VERSION, MAX_VERSION + 1))
CHOICES = set([str(i) for i in VERSIONS]) | {"current"}



CONVERSION_DICT = {
    "none":{},
    "upgrade":{
                2:"upgrade_1_to_2.sql",
                3:"upgrade_2_to_3.sql",
                4:"upgrade_3_to_4.sql",
                5:"upgrade_4_to_5.sql",
                6:"upgrade_5_to_6.sql",
                7:"upgrade_6_to_7.sql"
                },
    "downgrade":{
                1:"downgrade_2_to_1.sql",
                2:"downgrade_3_to_2.sql",
                3:"downgrade_4_to_3.sql",
                4:"downgrade_5_to_4.sql",
                5:"downgrade_6_to_5.sql",
                6:"downgrade_7_to_6.sql"
                }
    }



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
    parser.add_argument("-s", "--schema_version", type=str.lower,
        choices=list(CHOICES), default="current",
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

    database_versions_list = phamerator.retrieve_data(
            sql_handle, query='SELECT * FROM version')
    return database_versions_list[0]



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

    if args.schema_version == "current":
        args.schema_version = MAX_VERSION
    sql_handle = connect_to_db(args.database)
    version_data = retrieve_database_version(sql_handle)

    # Get the schema versions.
    target = int(args.schema_version)
    actual = get_schema_version(version_data)

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

    # Get the list of sql scripts needed.
    script_dict = CONVERSION_DICT[dir]

    # Iterate through list of versions and implement SQL files.
    if dir == "none":
        print("No schema conversion is needed.")
    else:
        index = 0
        convert = True
        while (index < len(steps) and convert == True):
            step = steps[index]
            script = script_dict[step]
            print(f"{dir[:-1].capitalize()}ing to schema version {step}...")
            data = pkgutil.get_data("pdm_utils", "sql_scripts/" + script).decode()
            commands = data.split(";")
            commands = [i for i in commands if i != ""]
            commands = [i for i in commands if i != "\n"]
            # Try to convert the schema.
            try:
                sql_handle.execute_transaction(commands)
            except:
                convert = False
            index += 1


###
