"""Pipeline to upgrade or downgrade the schema of a MySQL database."""

import argparse
import pathlib
import sys
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb

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
                    "schema of a MySQL database.")
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


# TODO unittest.
def main(unparsed_args_list):
    """Run main conversion pipeline."""
    # Parse command line arguments
    args = parse_args(unparsed_args_list)
    args.conversion_scripts_folder = \
        basic.set_path(args.conversion_scripts_folder, kind="dir", expect=True)
    sql_handle1 = mysqldb.connect_to_db(args.database)
    target = args.schema_version
    actual = mysqldb.get_schema_version(sql_handle1)
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
            result = mysqldb.drop_create_db(sql_handle1, args.new_database_name)
            if result == 0:
                result = mysqldb.copy_db(sql_handle1, args.new_database_name)
                if result == 0:
                    # Create a new connection to the new database.
                    sql_handle2 = mysqldb.setup_sql_handle2(
                                    username=sql_handle1.username,
                                    password=sql_handle1.password,
                                    database=args.new_database_name)
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
