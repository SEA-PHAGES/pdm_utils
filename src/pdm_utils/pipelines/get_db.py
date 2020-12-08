"""Pipeline to install a pdm_utils MySQL database."""

import argparse
from datetime import date
import pathlib
import shutil
import sys
import urllib3
import re

from pdm_utils.classes.alchemyhandler import MySQLDatabaseError
from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import basic, mysqldb, mysqldb_basic
from pdm_utils.functions import configfile
from pdm_utils.functions import pipelines_basic
from pdm_utils.pipelines import convert_db

DEFAULT_OUTPUT_FOLDER = "/tmp/"
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_db"

DB_LINK = constants.DB_SERVER


# TODO test.
def main(unparsed_args_list):
    """
    Run the get_db pipeline.

    The database data can be retrieved from three places:
    The server, which needs to be downloaded to a new folder.
    A local file, in which no download and no new folder are needed.
    The empty schema stored within pdm_utils, in which no download, new folder,
    or local file are needed.

    :param unparsed_args_list: list of arguments to run the pipeline unparsed
    :type unparsed_args_list: list

    """
    args = parse_args(unparsed_args_list)

    # Set values that are shared between all three options.
    config = configfile.build_complete_config(args.config_file)
    alchemist = pipelines_basic.build_alchemist(None, config=config,
                                                ask_database=False)

    if args.option == "file":
        execute_get_file_db(alchemist, args.database, args.filename,
                            config_file=args.config_file,
                            schema_version=args.schema_version,
                            verbose=args.verbose)
    elif args.option == "new":
        execute_get_new_db(alchemist, args.database, args.schema_version,
                           config_file=args.config_file, verbose=args.verbose)
    else:
        execute_get_server_db(
                    alchemist, args.database, args.url,
                    folder_path=args.output_folder, db_name=args.db_name,
                    config_file=args.config_file, verbose=args.verbose,
                    subdirectory=args.remote_directory,
                    download_only=args.download_only,
                    get_version=args.get_version,
                    schema_version=args.schema_version)


# TODO test.
def parse_args(unparsed_args_list):
    """
    Verify the correct arguments are selected for getting a new database.
    :param unparsed_args_list: arguments in sys.argv format
    :type unparsed_args_list: list
    :returns: A parsed list of arguments
    :rtype: argparse.Namespace
    """

    get_db_help = ("Pipeline to retrieve and install a new version of a MySQL "
                   "database.")
    database_help = "Name of the MySQL database."
    option_help = "Source of data to create database."
    server_help = "Download database from server."
    url_help = "Server URL from which to retrieve files."
    VERBOSE_HELP = """
        Toggles get_db pipeline progress print statements.
        """
    REMOTE_DIRECTORY_HELP = """
        Remote directory at the server URL from which to retrieve files.
        """
    DB_NAME_HELP = """
        MySQL export option to allow renaming of the exported database.
            Follow selection argument with the name of the desired database.
        """
    get_version_help = "Indicates that a .version file should be downloaded."
    output_folder_help = f"Path to the folder to create the folder for " \
                         f"downloading the database. Default is " \
                         f"{DEFAULT_OUTPUT_FOLDER}"
    download_only_help = "The database should be downloaded but not " \
                         "installed locally."
    file_help = "Install database from a SQL file."
    filename_help = "Name of the SQL file and version file."
    new_help = "Create a new empty database."
    schema_version_help = "Database schema version to which the database " \
                          "should be converted."
    config_file_help = "Path to the file containing user-specific login data."

    # database optional for subparser a and required for b and c

    parser = argparse.ArgumentParser(description=get_db_help)

    subparsers = parser.add_subparsers(dest="option", help=option_help)

    # Command line structure
    # python3 -m pdm_utils get_db server -db <database>
    parser_a = subparsers.add_parser("server", help=server_help)
    parser_a.add_argument("-db", "--database", type=str, help=database_help,
                          default=None)
    parser_a.add_argument("-u", "--url", type=str,
                          default=constants.DB_SERVER, help=url_help)
    parser_a.add_argument("-gv", "--get_version", action="store_true",
                          default=False, help=get_version_help)
    parser_a.add_argument("-o", "--output_folder", type=pathlib.Path,
                          default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                          help=output_folder_help)
    parser_a.add_argument("-d", "--download_only", action="store_true",
                          default=False, help=download_only_help)
    parser_a.add_argument("-n", "--db_name", help=DB_NAME_HELP)
    parser_a.add_argument("-rd", "--remote_directory", type=pathlib.Path,
                          help=REMOTE_DIRECTORY_HELP, default=None)

    parser_b = subparsers.add_parser("file", help=file_help)
    parser_b.add_argument("database", type=str, help=database_help)
    parser_b.add_argument("filename", type=pathlib.Path, help=filename_help)

    parser_c = subparsers.add_parser("new", help=new_help)
    parser_c.add_argument("database", type=str, help=database_help)

    # Add config file option to all subparsers.
    # It could be added after database, but then the optional argument is
    # required to be placed before the required subparser option, which
    # doesn't make sense.
    for p in [parser_a, parser_b, parser_c]:
        p.add_argument("-c", "--config_file", type=pathlib.Path,
                       help=config_file_help, default=None)
        p.add_argument("-v", "--verbose", action="store_true",
                       help=VERBOSE_HELP)
        p.add_argument("-s", "--schema_version", type=int,
                       choices=list(convert_db.CHOICES),
                       default=convert_db.CURRENT_VERSION,
                       help=schema_version_help)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


def execute_get_server_db(
                alchemist, database, url, folder_path=None,
                folder_name=RESULTS_FOLDER, db_name=None,
                config_file=None, verbose=False,
                subdirectory=None, download_only=False, get_fastas=False,
                get_alns=False, get_version=False, schema_version=None):
    config = configfile.build_complete_config(config_file)

    if db_name is None:
        db_name = database

    if database is None:
        print("Interactive mode not supported at this time, please manually "
              "select a database using the -db flag.")
        sys.exit(1)
        database = interactive()

    # Give priority to config file to define url, although this is
    # arbitrary
    server_creds = config["download_server"]
    config_url = server_creds.get("url")
    if config_url is not None:
        url = config_url

    if subdirectory:
        url = "".join([url, str(subdirectory), "/"])

    pool = urllib3.PoolManager()
    response = pool_request(url, pool=pool, pipeline=True)
    directory_listing = get_url_listing(response)

    valid_db_request = False
    for data_dict in directory_listing:
        directory_name = data_dict.get("name")
        if directory_name is None:
            continue

        if not directory_name.endswith("/"):
            continue

        if database == directory_name.split("/")[0]:
            valid_db_request = True
            break

    if not valid_db_request:
        print("Requested database is not at the specified url.\n"
              "Please check the database availability.")
        return

    response.close()

    pkg_url = "".join([url, database, "/"])
    pkg_response = pool_request(pkg_url, pool=pool, pipeline=True)
    pkg_directory_listing = get_url_listing(pkg_response)

    database_filename = None
    version_filename = None
    for data_dict in pkg_directory_listing:
        name = data_dict.get("name", "")
        if name.endswith(".sql"):
            database_filename = name.split(".sql")[0]
        if name.endswith(".version"):
            version_filename = name.split(".version")[0]

    output_path = pipelines_basic.create_working_path(folder_path, folder_name)
    output_path.mkdir(exist_ok=True)

    if not output_path.is_dir():
        print("Unable to create results folder.\n Aborting pipeline.")
        return

    # Only look for version file is selected.
    if get_version:
        if version_filename is not None:
            version_filepath, status1 = prepare_download(
                                            output_path, pkg_url,
                                            version_filename,
                                            "version", verbose=verbose)
    else:
        status1 = True

    db_filepath, status2 = prepare_download(
                                    output_path, pkg_url, database_filename,
                                    "sql", verbose=verbose)
    if not status1 or not status2:
        print("Unable to download data from server.\n Aborting pipeline.")
        return

    # If downloading from server, user may have selected to not
    # install the database file.
    if (not download_only) and (not get_fastas) and (not get_alns):
        install_db(alchemist, db_name, db_filepath=db_filepath,
                   config_file=config_file, schema_version=schema_version,
                   verbose=verbose)

        # The output folder was only created for downloading from server.
        print("Removing downloaded data.")
        shutil.rmtree(output_path)


def execute_get_file_db(alchemist, database, filename, config_file=None,
                        schema_version=None, verbose=False):
    db_filepath = basic.set_path(filename, kind="file", expect=True)

    install_db(alchemist, database, db_filepath=db_filepath,
               config_file=config_file, schema_version=schema_version,
               verbose=verbose)


def execute_get_new_db(alchemist, database, schema_version,
                       config_file=None, verbose=False):
    install_db(alchemist, database, schema_version=schema_version,
               config_file=config_file, verbose=verbose)


# TODO test.
def install_db(alchemist, database, db_filepath=None, config_file=None,
               schema_version=None, verbose=False):
    """
    Install database. If database already exists, it is first removed.
    :param database: Name of the database to be installed
    :type database: str
    :param db_filepath: Directory for installation
    :type db_filepath: Path
    """
    # No need to specify database yet, since it needs to first check if the
    # database exists.
    engine = alchemist.engine
    result = mysqldb_basic.drop_create_db(engine, database)
    engine.dispose()
    if result != 0:
        print("Unable to create new, empty database.")
        sys.exit(1)

    alchemist.database = database
    try:
        alchemist.validate_database()
    except MySQLDatabaseError:
        print(f"No connection to database {database} due "
              "to invalid credentials or database.")
        sys.exit(1)

    alchemist.build_engine()
    engine = alchemist.engine

    if db_filepath is not None:
        mysqldb_basic.install_db(engine, db_filepath)
    else:
        mysqldb.execute_transaction(engine, db_schema_0.STATEMENTS)

    if schema_version is not None:
        curr_schema_version = mysqldb.get_schema_version(engine)

        if not curr_schema_version == schema_version:
            if verbose:
                print(f"Schema version {curr_schema_version} "
                      "database detected.\nBeginning database conversion to "
                      f"schema version {schema_version}...")
            convert_args = ["pdm_utils.run", "convert", database,
                            "-s", str(schema_version)]
            if verbose:
                convert_args.append("-v")

            if config_file is not None:
                convert_args.extend(["-c", config_file])

            convert_db.main(convert_args)

    engine.dispose()


# TODO test.
def prepare_download(local_folder, url_folder, db_name, extension,
                     verbose=False):
    """
    Construct filepath and check if it already exists, then download.
    :param local_folder: Working directory where the database is downloaded
    :type local_folder: Path
    :param url_folder: Base url where db_files are located.
    :type url_folder: str
    :param db_name: Name of the database to be downloaded
    :type db_name: str
    :param extension: file extension for the database
    :type extension: str
    :returns: Path to the destination directory and the status of the download
    :rtype: Path, bool
    """
    filename = ".".join([db_name, extension])
    url_path = url_folder + filename
    local_path = pathlib.Path(local_folder, filename)
    if local_path.exists():
        print(f"The file {filename} already exists.")
        status = False
    else:
        if verbose:
            print(" ".join(["Preparing to download:", str(db_name)]))
        status = download_file(url_path, local_path, verbose=verbose)
    return local_path, status


# TODO if the file is not found on the server, a file will still be
# created with text indicating there was an error.
# So this step needs to catch that error.
# TODO test.
def download_file(file_url, filepath, pool=None, chunk_size=512000,
                  preload=False, verbose=True):
    """
    Retrieve a file from the server.
    :param file_url:  URL for the desired file:
    :type file_url: str
    :param filepath: Local path where the file can be downloaded to.
    :type filepath: Path
    :returns: Status of the file retrieved from the server
    :rtype: bool
    """
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == 200:
        filesize = int(request.getheader("Content-Length"))
        output_file = filepath.open(mode="wb")
        for chunk in request.stream(chunk_size):
            output_file.write(chunk)
            progress = int(request._fp_bytes_read / filesize * 100)
            if verbose:
                print("\r[{}{}] {}%".format("#" * int(progress / 2),
                                            " " * (50 - int(progress / 2)),
                                            progress), end="")
    else:
        if verbose:
            print(" ".join(["ERROR.  HTTP Response", str(status)]))

    if verbose:
        print("")

    return status


def pool_request(url, pool=False,
                 pipeline=False, preload=False, expect_status=200):
    """
    Create a urllib3 PoolManager to access the url link for list of databases
    :returns An HTTPResponse object
    :rtype: urllib3.response.HTTPResponse
    """
    if pool is None:
        pool = urllib3.PoolManager()
    response = pool.request("GET", url, preload_content=preload)

    pool.clear()

    if pipeline:
        if response.status != expect_status:
            print("Received invalid response from server.\n"
                  "Aborting pipeline.")
            sys.exit(1)

    return response


def get_url_listing(response, get_dates=False):
    """
    Get list of databases from the link using a response object
    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :returns: A list of databases
    :rtype: list
    """
    databases = list()

    response_data = response.data.decode("utf-8")
    # print(response_data)

    # Get database names
    names_regex = """<a href="([-_\w\d./]+)">"""
    names_regex = re.compile(names_regex)
    names = names_regex.findall(response_data)

    if get_dates:
        # Get database dates
        dates_regex = ("""<td align="right">"""
                       "(\d+[-]\d+[-]\d+)\s+\d+[:]\d+\s+</td>")
        dates_regex = re.compile(dates_regex)
        dates = dates_regex.findall(response_data)

    response.close()

    # creating a list of dictionaries
    for i in range(len(names)):
        db_dict = dict()
        db_dict["num"] = i+1
        db_dict["name"] = names[i]
        if get_dates:
            db_dict["date"] = dates[i]

        databases.append(db_dict)

    return databases


def interactive():
    """
    Interactive mode to display all available databases at specified url link
    for download from server
    :returns: Name of the database for download
    :rtype: str
    """
    response = pool_request(DB_LINK)

    if response.status == 200:
        # get the names of the databases

        databases = get_url_listing(response, get_dates=True)

        print("Databases available at '", DB_LINK, "':\n")

        for i in databases:
            print("{:10}.\t{:<30} {}".format(i["num"], i["name"], i["date"]))

        prompt = f"\n\nWhich database would you like to download? (Enter " \
                 f"1-{len(databases)}) "
        db = input(prompt)

        selected = databases[int(db)-1]
        return selected["name"]

    else:
        response.close()
        return None
