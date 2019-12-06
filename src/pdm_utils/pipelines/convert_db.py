"""Pipeline to upgrade or downgrade the schema of a Phamerator database."""

import argparse
import pathlib
import subprocess
import sys
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils.functions import basic
from pdm_utils.functions import phamerator

MAX_VERSION = 7
CURRENT_VERSION = 6
VERSIONS = list(range(0, MAX_VERSION + 1))
CHOICES = set(VERSIONS)

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
    filename = pathlib.Path(f"{dir}_{first}_to_{second}.sql")
    return filename


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for converting database."""
    CONVERT_HELP = ("Pipeline to upgrade or downgrade the "
                    "schema of a Phamerator database.")
    DATABASE_HELP = "Name of the MySQL database."
    CONVERSION_SCRIPTS_FOLDER_HELP = \
        ("Path to the folder containing conversion SQL scripts.")
    SCHEMA_VERSION_HELP = ("Database schema version to which "
                           "the database should be converted.")
    NEW_DATABASE_NAME_HELP = ("The new name of the converted database"
                              "if different from the original database name.")
    parser = argparse.ArgumentParser(description=CONVERT_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("conversion_scripts_folder", type=pathlib.Path,
        help=CONVERSION_SCRIPTS_FOLDER_HELP)
    parser.add_argument("-s", "--schema_version", type=int,
        choices=list(CHOICES), default=CURRENT_VERSION,
        help=SCHEMA_VERSION_HELP)
    parser.add_argument("-n", "--new_database_name", type=str,
        help=NEW_DATABASE_NAME_HELP)



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
    db_tables = get_db_tables(sql_handle, sql_handle.database)
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
            phage_columns = get_table_columns(sql_handle, sql_handle.database,
                                              "phage")
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
def get_db_tables(sql_handle, database):
    """Retrieve tables names from the database.

    Returns a set of table names."""
    query = ("SELECT table_name FROM information_schema.tables "
             f"WHERE table_schema = '{database}'")
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    db_tables = get_values_from_dict_list(result_list)
    return db_tables


# TODO unittest.
def get_table_columns(sql_handle, database, table_name):
    """Retrieve columns names from the phage table.

    Returns a set of column names."""
    query = ("SELECT column_name FROM information_schema.columns WHERE "
              f"table_schema = '{database}' AND "
              f"table_name = '{table_name}'")
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    columns = get_values_from_dict_list(result_list)
    return columns

# TODO unittest.
def get_db_names(sql_handle):
    """Retrieve database names from MySQL.

    Returns a set of database names."""
    query = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA")
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    databases = get_values_from_dict_list(result_list)
    return databases



# TODO unittest.
def get_values_from_dict_list(list_of_dicts):
    """Convert a list of dictionaries to a set of the values."""
    output_set = set()
    for dict in list_of_dicts:
        output_set = output_set | set(dict.values())
    return output_set


# TODO not unittested, but identical to function in export_db.
def retrieve_database_version(sql_handle):
    """Gets the database version and schema version.
    """
    result = phamerator.retrieve_data(sql_handle,
                                      query="SELECT * FROM version")
    return result[0]


# TODO unittest.
def main(unparsed_args_list):
    """Run main conversion pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    args.conversion_scripts_folder = \
        basic.set_path(args.conversion_scripts_folder, kind="dir", expect=True)
    sql_handle1 = connect_to_db(args.database)
    target = args.schema_version
    actual = get_schema_version(sql_handle1)
    steps, dir = get_conversion_direction(actual, target)

    # Iterate through list of versions and implement SQL files.
    if dir == "none":
        print("No schema conversion is needed.")
        convert = False
    else:
        convert = True

    if convert == True:
        if (args.new_database_name is not None and
                args.new_database_name != args.database):
            result = create_new_db(sql_handle1, args.new_database_name)
            if result == 0:
                result = copy_db(sql_handle1, args.new_database_name)
                if result == 0:
                    # Create a new connection to the new database.
                    sql_handle2 = mch.MySQLConnectionHandler()
                    sql_handle2.username = sql_handle1.username
                    sql_handle2.password = sql_handle1.password
                    sql_handle2.database = args.new_database_name
                else:
                    print("Unable to copy the database for conversion.")
                    convert = False
            else:
                print("Unable to create the new database for conversion.")
                convert = False
        else:
            sql_handle2 = sql_handle1

        if convert == True:
            stop_step, summary = convert_schema(sql_handle2,
                                       args.conversion_scripts_folder,
                                       actual, dir, steps)
            if stop_step == target:
                print("\n\nThe database schema conversion was successful.")
            else:
                print("\n\nThe database schema conversion was not successful. "
                      f"Unable to proceed past schema version {stop_step}.")
            print_summary(summary)

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



# TODO similar to function in get_db.py
# TODO unittest.
def create_new_db(sql_handle, database):
    """Creates a new, empty database."""
    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a new database is installed.
    databases = get_db_names(sql_handle)
    if database in databases:
        statement1 = [f"DROP DATABASE {database}"]
        result = sql_handle.execute_transaction(statement1)
        sql_handle.close_connection()
    else:
        result = 0
    if result == 0:
        statement2 = [f"CREATE DATABASE {database}"]
        result = sql_handle.execute_transaction(statement2)
        sql_handle.close_connection()
    return result


# TODO unittest.
def copy_db(sql_handle, new_database):
    """Copies a database.

    The sql_handle contains pointer to the name of the database
    that will be copied into the new_database parameter.
    """
    #mysqldump -u root -pPWD database1 | mysql -u root -pPWD database2
    command_string1 = ("mysqldump "
                      f"-u {sql_handle.username} "
                      f"-p{sql_handle.password} "
                      f"{sql_handle.database}")
    command_string2 = ("mysql -u "
                      f"{sql_handle.username} "
                      f"-p{sql_handle.password} "
                      f"{new_database}")
    command_list1 = command_string1.split(" ")
    command_list2 = command_string2.split(" ")
    try:
        print("Copying database...")

        # Per subprocess documentation:
        # 1. For pipes, use Popen instead of check_call.
        # 2. Call p1.stdout.close() to allow p1 to receive a SIGPIPE if p2 exits.
        #    which gets called when used as a context manager.
        # communicate() waits for the process to complete.
        with subprocess.Popen(command_list1, stdout=subprocess.PIPE) as p1:
            with subprocess.Popen(command_list2, stdin=p1.stdout) as p2:
                p2.communicate()
        print("Copy complete.")
        result = 0
    except:
        print(f"Unable to copy {sql_handle.database} to {new_database} in MySQL.")
        result = 1
    return result

# TODO unittest.
def convert_schema(sql_handle, conversion_folder, actual, dir, steps):
    """Iterate through conversion steps and convert database schema."""
    summary = {} #Key = conversion step. Value = summary dictionary.
    index = 0
    convert = True
    stop_step = actual
    while (index < len(steps) and convert == True):
        step = steps[index]
        script = get_script_filename(dir, step)
        print(f"{dir[:-1].capitalize()}ing to schema version {step}...")
        data = get_script_data(conversion_folder, script)
        step_summary_dict = get_script_summary(data)
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
            summary[script] = step_summary_dict
        index += 1
    return stop_step, summary

# TODO unittest.
def get_script_data(conversion_folder, script):
    """Confirm path to SQL conversion script is correct."""
    script_path = pathlib.Path(conversion_folder, script)
    result, msg = basic.verify_path2(script_path, kind="file", expect=True)
    if not result:
        print(msg)
        print("Unable to proceed with database conversion.")
        sys.exit(1)
    else:
        with open(script_path, "r") as fh:
            data = fh.read()
        return data


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


# TODO unittest.
def get_script_summary(data):
    """Retrieve summary of which data could be lost or inaccurate."""
    summary_dict = {}
    split_data = data.split("### DATA_LOSS_SUMMARY")
    if len(split_data) == 2:
        summaries_string = split_data[1]
        summaries_list = summaries_string.split("\n")
        summaries_list = [i[2:].split(":") for i in summaries_list]
        summaries_list = [i for i in summaries_list if i != [""]]
        summaries_list = [i for i in summaries_list if i != ["\n"]]
        keys = set()
        for summary in summaries_list:
            keys.add(summary[0])
        for key in keys:
            new_list = []
            for summary in summaries_list:
                if summary[0] == key:
                    new_list.append(summary[1])
            summary_dict[key] = new_list
    return summary_dict



# TODO unittest.
def print_summary(summary):
    """Print summary of data that could be lost or inaccurate."""
    keys = summary.keys()
    if len(keys) > 0:
        for key in keys:
            step = summary[key]
            step_keys = step.keys()
            if len(step_keys) > 0:
                print(f"\n\n\nConversion step: {key.stem}")
            if "LOST_TABLE" in step_keys:
                print("\nWarning: data in the following tables have been lost:")
                for item in step["LOST_TABLE"]:
                    print(item)
            if "LOST_COLUMN" in step_keys:
                print("\nWarning: data in the following columns have been lost:")
                for item in step["LOST_COLUMN"]:
                    print(item)
            if "LOST_DATA" in step_keys:
                print("\nWarning: a subset of data in the "
                      "following columns have been lost:")
                for item in step["LOST_DATA"]:
                    print(item)
            if "INACCURATE_COLUMN" in step_keys:
                print("\nWarning: data in the following columns have "
                       "been auto-populated and may be inaccurate:")
                for item in step["INACCURATE_COLUMN"]:
                    print(item)
