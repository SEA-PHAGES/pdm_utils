"""Pipeline to upgrade or downgrade the schema of a MySQL database."""

import argparse
import pathlib
import sys

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb, mysqldb_basic
from pdm_utils.constants import schema_conversions

MAX_VERSION = 10
CURRENT_VERSION = 10
VERSIONS = list(range(0, MAX_VERSION + 1))
CHOICES = set(VERSIONS)


# TODO unittest.
def get_step_name(dir, step):
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
    step_name = f"{dir}_{first}_to_{second}"
    return step_name


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for converting database."""
    convert_help = (
        "Pipeline to upgrade or downgrade the schema of a MySQL database.")
    database_help = "Name of the MySQL database."
    schema_version_help = (
        "Database schema version to which the database should be converted.")
    new_database_name_help = (
        "The new name of the converted database if different from "
        "the original database name.")
    verbose_help = (
        "Conversion progress will be printed.")
    config_file_help = (
        "Path to the file containing user-specific login data.")

    parser = argparse.ArgumentParser(description=convert_help)
    parser.add_argument("database", type=str, help=database_help)
    parser.add_argument("-s", "--schema_version", type=int,
        choices=list(CHOICES), default=CURRENT_VERSION,
        help=schema_version_help)
    parser.add_argument("-n", "--new_database_name", type=str,
        help=new_database_name_help)
    parser.add_argument("-v", "--verbose", action="store_true",
        default=False, help=verbose_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)

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
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]
    alchemist1 = AlchemyHandler(database=args.database,
                                username=mysql_creds["user"],
                                password=mysql_creds["password"])
    alchemist1.connect(pipeline=True)
    engine1 = alchemist1.engine


    target = args.schema_version
    actual = mysqldb.get_schema_version(engine1)
    steps, dir = get_conversion_direction(actual, target)

    # Iterate through list of versions and implement SQL files.
    if dir == "none":
        if args.verbose == True:
            print("No schema conversion is needed.")
        convert = False
    else:
        convert = True

    if convert == True:
        if (args.new_database_name is not None and
                args.new_database_name != args.database):
            result = mysqldb_basic.drop_create_db(engine1, args.new_database_name)
            if result == 0:
                result = mysqldb_basic.copy_db(engine1, args.new_database_name)
                if result == 0:
                    # Create a new connection to the new database.
                    alchemist2 = AlchemyHandler(database=args.new_database_name,
                                                username=engine1.url.username,
                                                password=engine1.url.password)
                    alchemist2.connect(pipeline=True)
                    engine2 = alchemist2.engine

                else:
                    print("Error: Unable to copy the database for conversion.")
                    convert = False
            else:
                print("Error: Unable to create the new database for conversion.")
                convert = False
        else:
            engine2 = engine1

        if convert == True:
            stop_step, summary = convert_schema(engine2, actual, dir,
                                                steps, verbose=args.verbose)
            engine2.dispose()
            if stop_step == target:
                if args.verbose == True:
                    print("\n\nThe database schema conversion was successful.")
            else:
                print("\n\nError: "
                      "The database schema conversion was not successful. "
                      f"Unable to proceed past schema version {stop_step}.")
            if args.verbose == True:
                print_summary(summary)
    engine1.dispose()

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
def convert_schema(engine, actual, dir, steps, verbose=False):
    """Iterate through conversion steps and convert database schema."""
    summary = {} #Key = conversion step. Value = summary dictionary.
    index = 0
    convert = True
    stop_step = actual
    while (index < len(steps) and convert == True):
        step = steps[index]
        if verbose == True:
            print(f"{dir[:-1].capitalize()}ing to schema version {step}...")
        step_name = get_step_name(dir, step)
        step_dict = get_step_data(step_name)
        commands = step_dict["statements"]
        commands = [i for i in commands if i != ""]
        commands = [i for i in commands if i != "\n"]
        # Try to convert the schema.
        result = mysqldb.execute_transaction(engine, commands)
        if result == 1:
            convert = False
            print("Error encountered while executing MySQL statements.")
        if convert == False:
            print(f"Error: Unable to {dir} schema to version {step}.")
        else:
            stop_step = step
            summary[step_name] = step_dict["step_summary_dict"]
        index += 1
    return stop_step, summary


# TODO unittest.
def get_step_data(step_name):
    """Get dictionary of conversion step data."""
    if step_name not in schema_conversions.CONVERSION_STEPS.keys():
        print(f"Error: The conversion step {step_name} is not available.")
        print("Unable to proceed with database conversion.")
        sys.exit(1)
    else:
        step_dict = schema_conversions.CONVERSION_STEPS[step_name]
        return step_dict


# TODO unittest.
def print_summary(summary):
    """Print summary of data that could be lost or inaccurate."""
    keys = summary.keys()
    if len(keys) > 0:
        for key in keys:
            step = summary[key]
            step_keys = step.keys()
            if len(step_keys) > 0:
                print(f"\n\n\nConversion step: {key}")
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
