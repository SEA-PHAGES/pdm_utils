"""Pipeline to install a pdm_utils MySQL database."""

import argparse
from datetime import date
import pathlib
import shutil
import subprocess
import sys

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import basic, mysqldb, mysqldb_basic
from pdm_utils.pipelines import convert_db

DEFAULT_OUTPUT_FOLDER = "/tmp/"
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_db"

# TODO test.
def main(unparsed_args_list):
    """Run the get_db pipeline.

    The database data can be retrieved from three places:
    The server, which needs to be downloaded to a new folder.
    A local file, in which no download and no new folder are needed.
    The empty schema stored within pdm_utils, in which no download, new folder,
    or local file are needed.
    """
    args = parse_args(unparsed_args_list)

    # Set all values to get rid of args object and to set additional values.
    database = args.database
    option = args.option
    install = True
    schema_version = None
    db_filepath = None

    if option == "file":
        db_filepath = basic.set_path(args.filename, kind="file", expect=True)
    elif option == "new":
        schema_version = args.schema_version
    else:
        # option must be "server"
        output_folder = basic.set_path(args.output_folder, kind="dir", expect=True)
        download = True
        remove = True
        results_folder = pathlib.Path(RESULTS_FOLDER)
        results_path = basic.make_new_dir(output_folder, results_folder, attempt=50)
        if args.download_only == True:
            install = False
            remove = False

        if results_path is None:
            print("Unable to create results folder.")
            sys.exit(1)
        else:
            version_filepath, status1 = prepare_download(results_path,
                                            constants.DB_WEBSITE,
                                            args.database, "version")
            db_filepath, status2 = prepare_download(results_path,
                                            constants.DB_WEBSITE,
                                            args.database, "sql")
        if (status1 == False or status2 == False):
            print("Unable to download data from server.")
            sys.exit(1)

    # If downloading from server, user may have selected to not
    # install the database file.
    if install == True:
        install_db(database, db_filepath=db_filepath, schema_version=schema_version)

    # The output folder was only created for downloading from server.
    if option == "server":
        if remove == True:
            print("Removing downloaded data.")
            shutil.rmtree(results_path)

# TODO test.
def install_db(database, db_filepath=None, schema_version=None):
    """Install database. If database already exists, it is first removed."""
    # No need to specify database yet, since it needs to first check if the
    # database exists.

    alchemist1 = AlchemyHandler(database="")
    alchemist1.connect(pipeline=True)
    engine1 = alchemist1.engine
    result = mysqldb_basic.drop_create_db(engine1, database)
    if result != 0:
        print("Unable to create new, empty database.")
    else:
        alchemist2 = AlchemyHandler(database=database,
                                    username=engine1.url.username,
                                    password=engine1.url.password)
        alchemist2.connect(pipeline=True)
        engine2 = alchemist2.engine
        if engine2 is None:
            print(f"No connection to the {database} database due "
                  "to invalid credentials or database.")
        else:
            if db_filepath is not None:
                mysqldb_basic.install_db(engine2, db_filepath)
            else:
                mysqldb.execute_transaction(engine2, db_schema_0.STATEMENTS)
                convert_args = ["pdm_utils.run", "convert", database,
                                "-s", str(schema_version)]
                convert_db.main(convert_args, engine2)
            # Close up all connections in the connection pool.
            engine2.dispose()
    # Close up all connections in the connection pool.
    engine1.dispose()

# TODO test.
def prepare_download(local_folder, url_folder, db_name, extension):
    """Construct filepath and check if it already exists, then download."""
    filename = ".".join([db_name, extension])
    url_path = url_folder + filename
    local_path = pathlib.Path(local_folder, filename)
    if local_path.exists() == True:
        print(f"The file {filename} already exists.")
        status = False
    else:
        status = download_file(url_path, local_path)
    return local_path, status


# TODO if the file is not found on the server, a file will still be
# created with text indicating there was an error.
# So this step needs to catch that error.
# TODO test.
def download_file(file_url, filepath):
    """Retrieve a file from the server."""
    print(f"Downloading {filepath.name} file.")
    # Command line structure: curl website > output_file
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


# TODO test.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for getting a new database."""

    GET_DB_HELP = \
        "Pipeline to retrieve and install a new version of a MySQL database."
    DATABASE_HELP = \
        "Name of the MySQL database."
    OPTION_HELP = \
        "Source of data to create database."
    SERVER_HELP = \
        "Download database from server."
    OUTPUT_FOLDER_HELP = \
        ("Path to the folder to create the folder for downloading "
        f"the database. Default is {DEFAULT_OUTPUT_FOLDER}")
    DOWNLOAD_ONLY_HELP = \
        "The database should be downloaded but not installed locally."
    FILENAME_HELP = \
        "Name of the SQL file and version file."
    NEW_HELP = \
        "Indicates whether a new, empty database should be created."
    SCHEMA_VERSION_HELP = \
        "Database schema version to which the database should be converted."

    parser = argparse.ArgumentParser(description=GET_DB_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    subparsers = parser.add_subparsers(dest="option", help=OPTION_HELP)

    parser_a = subparsers.add_parser("server", help=SERVER_HELP)
    parser_a.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER), help=OUTPUT_FOLDER_HELP)
    parser_a.add_argument("-d", "--download_only", action="store_true",
        default=False, help=DOWNLOAD_ONLY_HELP)

    parser_b = subparsers.add_parser("file", help="Install from local file.")
    parser_b.add_argument("filename", type=pathlib.Path, help=FILENAME_HELP)

    parser_c = subparsers.add_parser("new", help="Create empty database.")
    parser_c.add_argument("-s", "--schema_version", type=int,
        choices=list(convert_db.CHOICES), default=convert_db.CURRENT_VERSION,
        help=SCHEMA_VERSION_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args
