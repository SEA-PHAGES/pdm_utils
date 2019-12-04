"""Pipeline to upgrade or downgrade the schema of a Phamerator database."""

import argparse
import pkgutil
import re
from pdm_utils.functions import phamerator

# The schema version wasn't formally tracked until schema version = 3.
MIN_VERSION = 0
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


# TODO unittest.
def get_schema_version(sql_handle):
    """Identify the schema version of the database_versions_list."""
    # If version table does not exist, schema_version = 0.
    # If no schema_version or SchemaVersion field,
    # it is either schema_version = 1 or 2.
    # If AnnotationAuthor, Program, AnnotationQC, and RetrieveRecord
    # columns are in phage table, schema_version = 2.
    db_tables = get_db_tables(sql_handle)
    if "version" in db_tables:
        version_table = True
    else:
        version_table = False

    if version_table == True:
        version_columns = retrieve_database_version(sql_handle)
        if "schema_version" in version_columns.keys():
            schema_version = version_columns["schema_version"]
        elif "SchemaVersion" in version_columns.keys():
            schema_version = version_columns["SchemaVersion"]
        else:
            phage_columns = get_phage_table_columns(sql_handle)
            expected = {"AnnotationAuthor", "Program",
                        "AnnotationQC", "RetrieveRecord"}
            diff = expected - phage_columns
            if len(diff) == 0:
                schema_version = 2
            else:
                schema_version = 1
    else:
        schema_version = 0
    return schema_version


# TODO unittest.
def get_db_tables(sql_handle):
    """Retrieve tables names from the database.

    Returns a set of table names."""
    query = ("SELECT table_name FROM information_schema.tables "
             f"WHERE table_schema = '{sql_handle.database}'")
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    db_tables = set()
    for dict in result_list:
        table = set(dict.values())
        db_tables = db_tables | table
    return db_tables


# TODO this should be abstracted to retrieve column names from any table.
# TODO unittest.
def get_phage_table_columns(sql_handle):
    """Retrieve columns names from the phage table.

    Returns a set of column names."""
    query = ("SELECT column_name FROM information_schema.columns WHERE "
              "table_schema = 'Actino_Draft' AND table_name = 'phage'")
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    columns = set()
    for result in result_list:
        columns = columns | set(result.values())
    return columns

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
def main(unparsed_args_list):
    """Run main conversion pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    sql_handle = connect_to_db(args.database)
    target = args.schema_version
    actual = get_schema_version(sql_handle)
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
        if (dir == "upgrade" and step == 3):
            commands = get_upgrade_2_to_3_commands(data)
        else:
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


# TODO unittest.
def get_upgrade_2_to_3_commands(data):
    """Parse commands to upgrade from 2 to 3."""
    # This upgrade step involves creating a stored procedure, which
    # needs to be parsed in a specific way.
    commands_step1 = data.split("DELIMITER //")

    # All statements prior to defined procedure.
    commands = commands_step1[0].split(";")

    # The statement defining the procedure.
    commands_step2 = commands_step1[1].split("DELIMITER;")
    commands_step3 = commands_step2[0].split("//")
    commands.append(commands_step3[0])

    # All statements following the defined procedure.
    commands.extend(commands_step2[1].split(";"))
    return commands
