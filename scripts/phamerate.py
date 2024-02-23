import argparse
import multiprocessing
import pathlib
import sys
import time

from commands import *

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_phamerate"
CPUS = multiprocessing.cpu_count()
RUN_MACHINE = f"remote"


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--database", type=str, default="Actino_Draft")
    parser.add_argument("-c", "--config", type=pathlib.Path, default=None)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path)
    parser.add_argument("--cpus", type=int, default=CPUS)
    parser.add_argument("--run_machine", type=str, default=RUN_MACHINE)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    phamerate(args.database, args.config, args.output_folder,
              args.cpus, args.run_machine)


def phamerate(database, config_file, output_folder, cpus, run_machine):
    output_folder.mkdir(exist_ok=True, parents=True)

    print("Phamerating phage gene clusters...")
    phammseqs(database, config_file, cpus)

    print("Searching for domains with rpsBLAST...")
    find_domains(database, config_file, output_folder, cpus)

    print("Searching for transmembrane domains with Deep-TMHMM")
    find_transmembrane(database, config_file, output_folder, run_machine)

    print("Exporting sql database...")
    export(database, config_file, output_folder, cpus, phams=True,
           name="Actino_Draft")

    # print("Downgrading to schema 10...")
    # convert(database, config_file, 10)

    # print("Exporting schema 10 database.")
    # export(f"{database}_v10", config_file, output_folder, cpus,
    #       name="Actino_Draft")

    # print("Downgrading to schema 8...")
    # convert(database, config_file, 8)

    # print("Exporting schema 8 database.")
    # export(f"{database}_v8", config_file, output_folder, cpus,
    #       name="Actino_Draft")

    print("Downgrading to schema 6...")
    convert(database, config_file, 6)

    print("Exporting schema 6 database.")
    export(f"{database}_v6", config_file, output_folder, cpus,
           name="Actino_Draft")


if __name__ == "__main__":
    main(sys.argv[1:])
