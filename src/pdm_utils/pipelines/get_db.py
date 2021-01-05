"""Pipeline to install a pdm_utils MySQL database."""

import argparse
from datetime import date
import pathlib
import shutil
import sys


from pdm_utils.classes.alchemyhandler import MySQLDatabaseError
from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import (basic, configfile, mysqldb, mysqldb_basic,
                                 pipelines_basic, pipeline_shells, url_basic)
from pdm_utils.pipelines import convert_db

DEFAULT_OUTPUT_FOLDER = "/tmp/"
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_get_db"

DEFAULT_SETTINGS = {"url": constants.DB_SERVER}


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
        url = args.url
        if url is None:
            server_creds = config["download_server"]
            url = server_creds.get("url")

        if url is None:
            url = DEFAULT_SETTINGS["url"]

        execute_get_server_db(
                    alchemist, args.database, url,
                    folder_path=args.output_folder, db_name=args.db_name,
                    config_file=args.config_file, verbose=args.verbose,
                    subdirectory=args.remote_directory,
                    download_only=args.download_only,
                    get_version=args.get_version, force_pull=args.force_pull,
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
        MySQL file option to allow renaming of the new database.
            Follow selection argument with the name of the desired database.
        """
    FORCE_PULL_HELP = """
        Get_db option to force a database pull and ignore database version
        concurrency.
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

    subparsers = parser.add_subparsers(dest="option", required=True,
                                       help=option_help)

    # Command line structure
    # python3 -m pdm_utils get_db server -db <database>
    parser_a = subparsers.add_parser("server", help=server_help)
    parser_a.add_argument("-db", "--database", type=str, help=database_help,
                          default=None)
    parser_a.add_argument("-u", "--url", type=str,
                          default=None, help=url_help)
    parser_a.add_argument("-gv", "--get_version", action="store_true",
                          default=False, help=get_version_help)
    parser_a.add_argument("-fp", "--force_pull", action="store_true",
                          help=FORCE_PULL_HELP)
    parser_a.add_argument("-o", "--output_folder", type=pathlib.Path,
                          default=None, help=output_folder_help)
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
                config_file=None, verbose=False, subdirectory=None,
                download_only=False, get_fastas=False, get_alns=False,
                force_pull=False, get_version=False, schema_version=None):
    

    if subdirectory:
        url = "".join([url, str(subdirectory), "/"])

    pool = url_basic.create_pool(pipeline=True)
    if database is None:
        print("Loading get_db interactive environment...")
        cmd = pipeline_shells.GetDBCMD(url, name=alchemist.username, pool=pool)
        cmd.cmdloop(intro=pipeline_shells.GET_DB_CMD_INTRO)

        if cmd.selected is None:
            return

        database = cmd.selected.name
        pkg_url = "".join([cmd.selected.get_abs_path(), "/"])
    else:
        response = url_basic.pool_request(url, pool=pool, pipeline=True)
        directory_listing = url_basic.get_url_listing_dirs(response)

        if database not in directory_listing:
            print("Requested database is not at the specified url.\n"
                  "Please check the database availability.")
            return

        response.close()
        pkg_url = "".join([url, database, "/"])

    if db_name is None:
        db_name = database

    pkg_response = url_basic.pool_request(pkg_url, pool=pool, pipeline=True)

    sql_file_listing = url_basic.get_url_listing_files(pkg_response, "sql")
    version_file_listing = url_basic.get_url_listing_files(
                                                       pkg_response, "version")

    if not sql_file_listing:
        print("Requested database file package does not have a SQL file.\n"
              "Please check SQL file availability at the specified url.")
        return
    database_filename = sql_file_listing[0]

    if not version_file_listing:
        if get_version:
            print("Requested database file package does not have"
                  "a version file.\nPlease check version file availability "
                  "at the specified url.")
            return
        else:
            version_filename = None
    else:
        version_filename = version_file_listing[0]

    if folder_path is None:
        output_path = pipelines_basic.create_working_path(
                                    pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                                    folder_name)
        if output_path.is_dir():
            shutil.rmtree(output_path)
        pipelines_basic.create_working_dir(output_path, force=True)
    else:
        output_path = pipelines_basic.create_working_path(
                                    folder_path, folder_name)
        pipelines_basic.create_working_dir(output_path)

    # Only look for version file is selected.
    if version_filename is not None:
        version_filepath, status1 = prepare_download(
                                            output_path, pkg_url,
                                            version_filename,
                                            "version")
        version_filehandle = version_filepath.open(mode="r")
        version = int(version_filehandle.readline().rstrip())
    else:
        status1 = True
        version = 0

    if (not force_pull) and (version > 0):
        if db_name in alchemist.databases:
            alchemist.database = db_name
            alchemist.build_engine()

            curr_schema_version = mysqldb.get_schema_version(alchemist.engine)
            if curr_schema_version > 2:
                curr_version_data = mysqldb_basic.get_first_row_data(
                                                alchemist.engine, "version")
                curr_version = int(curr_version_data.get("Version", 0))
                if curr_version >= version:
                    print(
                      f"Current database version of {db_name} "
                      "is greater than or equal to the database version "
                      "at the specified listing.\nPlease use "
                      "the --force_pull flag if you would like to "
                      "indiscriminately pull and install a database.")
                    return

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
# TODO move,
def install_db(alchemist, database, db_filepath=None, config_file=None,
               schema_version=None, verbose=False, pipeline=False):
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
        if pipeline:
            print("Unable to create new, empty database.\nPlease "
                  "check SQL service status and/or database availability.")
            sys.exit(1)

        raise OSError("Unable to create new, empty database.\nPlease "
                      "check SQL service status and/or database availability.")

    alchemist.database = database
    try:
        alchemist.validate_database()
    except MySQLDatabaseError:
        if pipeline:
            print(f"No connection to database {database} due "
                  "to invalid credentials and/or database.")
            sys.exit(1)

        raise MySQLDatabaseError(f"No connection to database {database} due "
                                 "to invalid credentials and/or database.")

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
        status = url_basic.download_file(url_path, local_path, verbose=verbose)
    return local_path, status
