"""Pipeline to check for new versions of the database from the server."""

import argparse
import os
import pathlib
import subprocess
import sys
from pdm_utils.constants import constants, db_schema_0
from pdm_utils.functions import basic, mysqldb
from pdm_utils.pipelines import convert_db

# TODO currently the 'output_folder' argument is not directly attached to
# other arguments. This causes a problem sometimes - if installing a database
# from a .sql file, the 'filename' argument points to the path to the file,
# so the 'output_folder' path is redundant. Further, if the 'output_folder'
# path is not the same as the 'filename' path, then an error is encountered.
# TODO unittest.
def main(unparsed_args_list):
    """Run the get_db pipeline."""
    args = parse_args(unparsed_args_list)
    args.output_folder = basic.set_path(args.output_folder, kind="dir", expect=True)

    if args.new == True:
        schema_version = choose_schema_version()

    # curl website > output_file
    # Version file
    if args.filename is None:
        version_filename = args.database + ".version"
    else:
        version_filename = args.filename.stem + ".version"
    version_url = constants.DB_WEBSITE + version_filename
    version_filepath = pathlib.Path(args.output_folder, version_filename)
    if args.download == True:
        if version_filepath.exists() == True:
            print(f"The file {version_filename} already exists.")
        else:
            command_string = f"curl {version_url}"
            command_list = command_string.split(" ")
            with version_filepath.open("w") as version_handle:
                try:
                    print("Downloading version file.")
                    subprocess.check_call(command_list, stdout=version_handle)
                except:
                    print(f"Unable to download {version_filename} from server.")

    # Database file
    if args.filename is None:
        db_filename = args.database + ".sql"
    else:
        db_filename = args.filename.name
    db_url = constants.DB_WEBSITE + db_filename
    db_filepath = pathlib.Path(args.output_folder, db_filename)
    if args.download == True:
        if db_filepath.exists() == True:
            print(f"The file {db_filename} already exists.")
        else:
            command_string2 = f"curl {db_url}"
            command_list2 = command_string2.split(" ")
            status = False
            # TODO if the file is not found on the server,
            # a file will still be created with text indicating
            # there was an error. So this step needs to catch that
            # error.
            with db_filepath.open("w") as db_handle:
                try:
                    print("Downloading sql file.")
                    subprocess.check_call(command_list2, stdout=db_handle)
                    status = True
                except:
                    print(f"Unable to download {db_filename} from server.")

    # Install new database
    if args.install == True:
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

    if args.remove == True:
        print("Removing downloaded data.")
        if version_filepath.exists() == True:
            os.remove(version_filepath)
        if db_filepath.exists() == True:
            os.remove(db_filepath)

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
    parser.add_argument("output_folder", type=pathlib.Path,
        help=OUTPUT_FOLDER_HELP)
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


    # TODO implement this option?
    # parser.add_argument("-f", "--force_update", action="store_true",
    #     default=False, help=FORCE_INSTALL_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])

    if args.new == True:
        args.all_steps == False
        args.download = False
        args.install = True
        args.remove = False
        args.filename = None
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
