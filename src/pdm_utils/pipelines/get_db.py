"""Pipeline to install a pdm_utils MySQL database."""

import argparse
from datetime import date
import pathlib
import shutil
import subprocess
import sys

from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import basic, mysqldb
from pdm_utils.pipelines import convert_db

DEFAULT_OUTPUT_FOLDER = "/tmp/"
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_db"

# TODO unittest.
def main(unparsed_args_list):
    """Run the get_db pipeline.

    The database data can be retrieved from three places:
    The server, in which is needs to be downloaded to a new folder.
    A local file, in which no download and no new folder are needed.
    The empty schema stored within pdm_utils, in which no download, new folder,
    or local file are needed.
    """
    args = parse_args(unparsed_args_list)

    # If filename was provided while also specifying to create new database
    # or download database, throw an error.
    if args.filename is not None:
        if (args.download == True or args.new == True):
            print(f"Conflicting arguments. If installing from local file, "
                  "do not also select to create a new database or "
                  "download database from server.")
            sys.exit(1)

    args.output_folder = basic.set_path(args.output_folder, kind="dir", expect=True)

    if args.filename is not None:
        args.filename = basic.set_path(args.filename, kind="file", expect=True)
        results_path = None
    elif args.new == True:
        results_path = None
    else:
        # If not installing from local file, create new dir to store
        # downloaded data.
        results_folder = pathlib.Path(RESULTS_FOLDER)
        results_path = basic.make_new_dir(args.output_folder,
                                          results_folder, attempt=50)
        if results_path is None:
            print("Unable to create results folder.")
            sys.exit(1)

    if args.new == True:
        schema_version = choose_schema_version()
    else:
        schema_version = None

    # curl website > output_file
    if args.download == True:
        version_filepath, status1 = get_data(results_path, constants.DB_WEBSITE,
                                             args.database, "version")
        db_filepath, status2 = get_data(results_path, constants.DB_WEBSITE,
                                        args.database, "sql")
    elif args.filename is not None:
        version_filepath = None
        db_filepath = args.filename
    else:
        version_filepath = None
        db_filename = None
        db_filepath = None

    if args.install == True:
        install_db(db_filepath, args, schema_version=schema_version)

    if args.remove == True:
        if results_path is not None:
            remove_data(results_path)




def install_db(db_filepath, args, schema_version=None):
    """Install database."""
    if args.new == True:
        result1 = [True, None]
    else:
        result1 = basic.verify_path2(db_filepath, kind="file", expect=True)
    if result1[0] == True:
        engine1, msg = mysqldb.get_engine(database="", echo=False)
        if engine1 is not None:
            result2 = mysqldb.drop_create_db(engine1, args.database)
            if result2 == 0:
                engine2, msg = mysqldb.get_engine(
                                    database=args.database,
                                    username=engine1.url.username,
                                    password=engine1.url.password,
                                    echo=False)
                if engine2 is not None:
                    if args.new == True:
                        mysqldb.execute_transaction(engine2, db_schema_0.STATEMENTS)
                        convert_args = ["pdm_utils.run",
                                        "convert",
                                        args.database,
                                        "-s",
                                        str(schema_version)]
                        convert_db.main(convert_args, engine2)
                    else:
                        mysqldb.install_db(engine2, db_filepath)
                    # Close up all connections in the connection pool.
                    engine2.dispose()
                else:
                    print(f"No connection to the {args.database} database due "
                          "to invalid credentials or database.")
            else:
                print("Unable to create new, empty database.")
            # Close up all connections in the connection pool.
            engine1.dispose()
        else:
            print("Invalid MySQL credentials.")
    else:
        print("Unable to locate database file for installation.")
        print(result1[1])


def remove_data(dir):
    """Remove the folder of downloaded files."""
    print("Removing downloaded data.")
    shutil.rmtree(dir)


def get_data(local_folder, url_folder, db_name, extension):
    """."""
    filename = ".".join([db_name, extension])
    url_path = url_folder + filename
    local_path = pathlib.Path(local_folder, filename)
    if local_path.exists() == True:
        print(f"The file {filename} already exists.")
    else:
        status = get_file(url_path, local_path)
    return local_path, status


# TODO if the file is not found on the server,
# a file will still be created with text indicating
# there was an error. So this step needs to catch that
# error.
def get_file(file_url, filepath):
    """Retrieve a file from the server."""
    print(f"Downloading {filepath.name} file.")
    command_string = f"curl {file_url}"
    command_list = command_string.split(" ")
    status = False
    with filepath.open("w") as fh:
        try:
            subprocess.check_call(command_list, stdout=fh)
            status = True
        except:
            print(f"Unable to download {filepath.name} from server.")
    return status


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting a new database."""

    UPDATE_DB_HELP = ("Pipeline to retrieve and install a new version of "
                   "a MySQL database.")
    DATABASE_HELP = "Name of the MySQL database."
    OUTPUT_FOLDER_HELP = ("Path to the folder to create the folder for "
                          "downloading the database.")
    INSTALL_HELP = \
        ("Indicates if the new version should be installed.")
    DOWNLOAD_HELP = \
        ("Indicates if the new version should be downloaded.")
    FORCE_INSTALL_HELP = \
        ("Indicates if the downloaded version should overwrite the existing "
         "database or create a new database if it is not already present.")
    REMOVE_HELP = \
        ("Indicates if the downloaded file should be removed after installation.")
    ALL_HELP = \
        ("Indicates if the entire pipeline should run: "
         "download, install, and then remove the temp files.")
    FILENAME_HELP = \
        ("Indicates the name of the SQL file and version file. "
         "By default, the filename is "
         "assumed to be the same as the database name. This option enables "
         "database to be created from a file of a different name.")
    NEW_HELP = \
        ("Indicates whether a new, empty database with the most current schema "
         "should be created. This option overrides all other options.")

    parser = argparse.ArgumentParser(description=UPDATE_DB_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER), help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-d", "--download", action="store_true",
        default=False, help=DOWNLOAD_HELP)
    parser.add_argument("-i", "--install", action="store_true",
        default=False, help=INSTALL_HELP)
    parser.add_argument("-r", "--remove", action="store_true",
        default=False, help=REMOVE_HELP)
    parser.add_argument("-a", "--all_steps", action="store_true",
        default=False, help=ALL_HELP)
    parser.add_argument("-f", "--filename", type=pathlib.Path,
        help=FILENAME_HELP)
    parser.add_argument("-n", "--new", action="store_true",
        default=False, help=NEW_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if args.new == True:
        args.all_steps == False
        args.download = False
        args.install = True
        args.remove = False
    elif args.all_steps == True:
        args.download = True
        args.install = True
        args.remove = True
    else:
        pass
    return args


# TODO unittest.
def choose_schema_version():
    prompt1 = "Do you want the current database schema? "
    response1 = basic.ask_yes_no(prompt=prompt1, response_attempt=3)
    if response1 is not None:
        if response1 == False:
            prompt2 = ("Select the schema version: "
                       f"\n{convert_db.VERSIONS} ")
            schema_version = basic.select_option(prompt2, convert_db.VERSIONS)
        else:
            schema_version = convert_db.CURRENT_VERSION
    else:
        print("The database will not be created.")
        sys.exit(1)
    return schema_version
