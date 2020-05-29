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
from pdm_utils.functions import configfile
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

    # Set values that are shared between all three options.
    database = args.database
    option = args.option
    config_file = args.config_file

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(config_file)
    mysql_creds = config["mysql"]
    server_creds = config["download_server"]

    install = True
    schema_version = None
    db_filepath = None

    if option == "file":
        db_filepath = basic.set_path(args.filename, kind="file", expect=True)
    elif option == "new":
        schema_version = args.schema_version
    else:
        # option must be "server"
        # Give priority to config file to define url, although this is arbitrary.
        server_url = server_creds["url"]
        if server_url is None:
            server_url = args.url
        version_file = args.version
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
            # Only look for version file is selected.
            if version_file:
                version_filepath, status1 = prepare_download(results_path,
                                                server_url, database,
                                                "version")
            else:
                status1 = True

            db_filepath, status2 = prepare_download(results_path, server_url,
                                                    database, "sql")
        if (status1 == False or status2 == False):
            print("Unable to download data from server.")
            sys.exit(1)

    # If downloading from server, user may have selected to not
    # install the database file.
    if install == True:
        install_db(database, username=mysql_creds["user"],
                   password=mysql_creds["password"], db_filepath=db_filepath,
                   schema_version=schema_version, config_file=config_file)

    # The output folder was only created for downloading from server.
    if option == "server":
        if remove == True:
            print("Removing downloaded data.")
            shutil.rmtree(results_path)

# TODO test.
def install_db(database, username=None, password=None, db_filepath=None,
               schema_version=None, config_file=None):
    """Install database. If database already exists, it is first removed."""
    # No need to specify database yet, since it needs to first check if the
    # database exists.
    alchemist1 = AlchemyHandler(database="", username=username, password=password)
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
                if config_file is not None:
                    convert_args.extend(["-c", config_file])
                convert_db.main(convert_args)
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

    get_db_help = (
        "Pipeline to retrieve and install a new version of a MySQL database.")
    database_help = "Name of the MySQL database."
    option_help = "Source of data to create database."
    server_help = "Download database from server."
    url_help = "Server URL from which to retrieve files."
    version_help = "Indicates that a .version file should be downloaded."
    output_folder_help = (
        "Path to the folder to create the folder for downloading "
        f"the database. Default is {DEFAULT_OUTPUT_FOLDER}")
    download_only_help = (
        "The database should be downloaded but not installed locally.")
    file_help = "Install database from a SQL file."
    filename_help = "Name of the SQL file and version file."
    new_help = "Create a new empty database."
    schema_version_help = (
        "Database schema version to which the database should be converted.")
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=get_db_help)
    parser.add_argument("database", type=str, help=database_help)

    subparsers = parser.add_subparsers(dest="option", help=option_help)

    parser_a = subparsers.add_parser("server", help=server_help)
    parser_a.add_argument("-u", "--url", type=str,
        default=constants.DB_WEBSITE, help=url_help)
    parser_a.add_argument("-v", "--version", action="store_true",
        default=False, help=version_help)
    parser_a.add_argument("-o", "--output_folder", type=pathlib.Path,
        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER), help=output_folder_help)
    parser_a.add_argument("-d", "--download_only", action="store_true",
        default=False, help=download_only_help)

    parser_b = subparsers.add_parser("file", help=file_help)
    parser_b.add_argument("filename", type=pathlib.Path, help=filename_help)

    parser_c = subparsers.add_parser("new", help=new_help)
    parser_c.add_argument("-s", "--schema_version", type=int,
        choices=list(convert_db.CHOICES), default=convert_db.CURRENT_VERSION,
        help=schema_version_help)

    # Add config file option to all subparsers.
    # It could be added after database, but then the optional argument is
    # required to be placed before the required subparser option, which
    # doesn't make sense.
    for p in [parser_a, parser_b, parser_c]:
        p.add_argument("-c", "--config_file", type=pathlib.Path,
                       help=config_file_help, default=None)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args
