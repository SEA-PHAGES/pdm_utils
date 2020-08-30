"""Pipeline to install a pdm_utils MySQL database."""

import argparse
from datetime import date
import pathlib
import shutil
import subprocess
import sys
import urllib3
import re

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import basic, mysqldb, mysqldb_basic
from pdm_utils.functions import configfile
from pdm_utils.pipelines import convert_db

DEFAULT_OUTPUT_FOLDER = "/tmp/"
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_db"

DB_LINK = "http://phamerator.webfactional.com/databases_Hatfull/"


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
        # Give priority to config file to define url, although this is
        # arbitrary

        if database is None:
            database = interactive()

        server_url = server_creds["url"]
        if server_url is None:
            server_url = args.url
        
        if args.remote_directory:
            server_url = "".join([server_url, str(args.remote_directory), "/"])

        version_file = args.version
        output_folder = basic.set_path(args.output_folder, kind="dir", expect=True)
        download = True
        remove = True
        results_folder = pathlib.Path(RESULTS_FOLDER)
        results_path = basic.make_new_dir(output_folder, results_folder, attempt=50)
        if args.download_only:
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
        if not status1 or not status2:
            print("Unable to download data from server.")
            sys.exit(1)

    # If downloading from server, user may have selected to not
    # install the database file.
    if install:
        install_db(database, username=mysql_creds["user"],
                   password=mysql_creds["password"], db_filepath=db_filepath,
                   schema_version=schema_version, config_file=config_file)

    # The output folder was only created for downloading from server.
    if option == "server":
        if remove:
            print("Removing downloaded data.")
            shutil.rmtree(results_path)


# TODO test.
def parse_args(unparsed_args_list):
    """
    Verify the correct arguments are selected for getting a new database.
    :param unparsed_args_list: arguments in sys.argv format
    :type unparsed_args_list: list
    :returns: A parsed list of arguments
    :rtype: argparse.Namespace
    """

    get_db_help = "Pipeline to retrieve and install a new version of a MySQL " \
                  "database."
    database_help = "Name of the MySQL database."
    option_help = "Source of data to create database."
    server_help = "Download database from server."
    url_help = "Server URL from which to retrieve files."
    REMOTE_DIRECTORY_HELP = """
        Remote directory at the server URL from which to retrieve files.
        """
    version_help = "Indicates that a .version file should be downloaded."
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
    # parser.add_argument("-db", "--database", type=str, help=database_help,
    # default=None)

    subparsers = parser.add_subparsers(dest="option", help=option_help)

    # Command line structure
    # python3 -m pdm_utils get_db server -db <database>
    parser_a = subparsers.add_parser("server", help=server_help)
    parser_a.add_argument("-db", "--database", type=str, help=database_help,
                          default=None)
    parser_a.add_argument("-u", "--url", type=str,
                          default=constants.DB_WEBSITE, help=url_help)
    parser_a.add_argument("-v", "--version", action="store_true",
                          default=False, help=version_help)
    parser_a.add_argument("-o", "--output_folder", type=pathlib.Path,
                          default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                          help=output_folder_help)
    parser_a.add_argument("-d", "--download_only", action="store_true",
                          default=False, help=download_only_help)
    parser_a.add_argument("-rd", "--remote_directory", type=pathlib.Path,
                                      help=REMOTE_DIRECTORY_HELP, default=None)

    parser_b = subparsers.add_parser("file", help=file_help)
    parser_b.add_argument("database", type=str, help=database_help)
    parser_b.add_argument("filename", type=pathlib.Path, help=filename_help)

    parser_c = subparsers.add_parser("new", help=new_help)
    parser_c.add_argument("database", type=str, help=database_help)
    parser_c.add_argument("-s", "--schema_version", type=int,
                          choices=list(convert_db.CHOICES),
                          default=convert_db.CURRENT_VERSION,
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


# TODO test.
def install_db(database, username=None, password=None, db_filepath=None,
               schema_version=None, config_file=None):
    """
    Install database. If database already exists, it is first removed.
    :param database: Name of the database to be installed
    :type database: str
    :param username: mySQL username
    :type username: str
    :param password: mySQL password
    :type password: str
    :param db_filepath: Directory for installation
    :type db_filepath: Path
    :param schema_version: Database schema version
    :type schema_version: int
    :param config_file: Config file with credentials available for pipeline use
    :type config_file: ConfigParser
    """

    # No need to specify database yet, since it needs to first check if the
    # database exists.
    alchemist1 = AlchemyHandler(database="", username=username,
                                password=password)
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
        status = download_file(url_path, local_path)
    return local_path, status


# TODO if the file is not found on the server, a file will still be
# created with text indicating there was an error.
# So this step needs to catch that error.
# TODO test.
def download_file(file_url, filepath):
    """
    Retrieve a file from the server.
    :param file_url: URL for database
    :type file_url: str
    :param filepath: Local path where the file is downloaded to.
    :type filepath: Path
    :returns: Status of the file retrieved from the server
    :rtype: bool
    """
    print(f"Downloading {filepath.name} file.")
    # Command line structure: curl website > output_file
    command_string = f"curl {file_url}"
    command_list = command_string.split(" ")
    status = False
    with filepath.open("w") as fh:
        try:
            subprocess.check_call(command_list, stdout=fh)
            status = True
        except subprocess.CalledProcessError:
            print(f"Unable to download {filepath.name} from server.")
    return status


def request_url():
    """
    Create a urllib3 PoolManager to access the url link for list of databases
    :returns An HTTPResponse object
    :rtype: urllib3.response.HTTPResponse
    """
    pool = urllib3.PoolManager()
    response = pool.request('GET', DB_LINK)

    pool.clear()

    return response


def get_database_list(response):
    """
    Get list of databases from the link using a response object
    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :returns: A list of databases
    :rtype: list
    """
    databases = list()

    # Get database names
    names_regex = """<a href="(\w+).sql">"""
    names_regex = re.compile(names_regex)
    names = names_regex.findall(response.data.decode('utf-8'))

    # Get database dates
    dates_regex = """<td align="right">(\d+[-]\d+[-]\d+)\s+\d+[:]\d+\s+</td>"""
    dates_regex = re.compile(dates_regex)
    dates = dates_regex.findall(response.data.decode('utf-8'))

    response.close()

    # creating a list of dictionaries
    for i in range(len(names)):
        db_dict = dict()
        db_dict["num"] = i+1
        db_dict["name"] = names[i]
        db_dict["date"] = dates[i*2]

        databases.append(db_dict)

    return databases


def interactive():
    """
    Interactive mode to display all available databases at specified url link for download from server
    :returns: Name of the database for download
    :rtype: str
    """
    response = request_url()

    if response.status == 200:
        # get the names of the databases

        databases = get_database_list(response)

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


