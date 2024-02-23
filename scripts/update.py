import argparse
import multiprocessing
import pathlib
import shlex
from subprocess import Popen, PIPE
import sys
import time

from pdm_utils.functions import pipelines_basic

from commands import *
from email_submitters import email_submitters

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_phamerate"
CPUS = multiprocessing.cpu_count()


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--database", type=str, default="Actino_Draft")
    parser.add_argument("-c", "--config", type=pathlib.Path, default=None)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path)
    parser.add_argument("-m", "--folder_name", type=str,
                        default=DEFAULT_FOLDER_NAME)
    parser.add_argument("--cpus", type=int, default=CPUS)
    parser.add_argument("-e", "--email_submitters", action="store_true")

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    execute_update(args.database, args.config, args.output_folder,
                   args.folder_name, args.cpus, args.email_submitters)


def execute_update(database, config_file,
                   output_folder, folder_name, cpus, email=False):
    working_dir = pipelines_basic.create_working_path(
                                        output_folder, folder_name)
    pipelines_basic.create_working_dir(working_dir)

#    print(f"Retreiving latest version of {database}...")
#    get_db(database, config_file)

    print("Pulling data from remote servers...")
    get_data(database, config_file, working_dir)

    data_dir = working_dir.joinpath(f"{time.strftime('%Y%m%d')}_get_data")
    print(data_dir)
    if not data_dir.is_dir():
        print("No data was able to be retrieved.  Exiting early.")
        return

    import_map = dict()
    print("Importing retrieved data...")
    for source_dir in data_dir.iterdir():
        if source_dir.name == "updates":
            continue

        print(f"\tImporting data for {source_dir.name}...")
        genomes_dir = source_dir.joinpath("genomes")
        import_table = source_dir.joinpath("import_table.csv")
        output_folder = working_dir.joinpath(f"{source_dir.name}_import")
        output_folder.mkdir()

        import_gb(database, config_file, output_folder,
                  genomes_dir, import_table)

        import_map[source_dir.name] = output_folder

    update_file = data_dir.joinpath("updates/update_table.csv")
    if not update_file.is_file():
        update_file = None

    print("Processing database updates...")
    update(database, config_file, update_file)

    phagesdb_import = import_map.get("phagesdb")
    if email and phagesdb_import:
        print("Emailing genome submitters and staff...")

        import_dir = phagesdb_import.joinpath(
                                        f"{time.strftime('%Y%m%d')}_import")
        email_submitters(import_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
