import argparse
import multiprocessing
import pathlib
import shlex
import sys
import time

import yaml

from commands import *

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_phamerate"
CPUS = multiprocessing.cpu_count()


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--indir", type=pathlib.Path, required=True)
    parser.add_argument("-n", "--database_name", type=str, required=True)
    parser.add_argument("-c", "--config", type=pathlib.Path, default=None)
    parser.add_argument("-k", "--permissions_file", type=pathlib.Path)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    upload(args.indir, args.database_name,
              args.config, args.permissions_file)


def upload(database_dir, database_name, config_file, pem_file):
    yaml_outfile = database_dir.joinpath(f"{database_name}.yml")
    version_file = database_dir.joinpath(f"{database_name}.version")

    curl_metadata(database_dir.name, yaml_outfile)

    increment_yaml(version_file, yaml_outfile)

    push(database_dir.stem, pem_file, config_file, database_dir)


def increment_yaml(version_file, yaml_file):
    with open(version_file, "r") as filehandle:
        lines = filehandle.readlines()

    version = int(lines[0].rstrip("\n"))

    with open(yaml_file, "r") as filehandle:
        lines = filehandle.readlines()

    lines[1] = f"version:    {version}\n"

    lines[2] = f"date:   {time.strftime('%Y-%m-%d')}\n"

    with open(yaml_file, "w") as filehandle:
        for line in lines:
            filehandle.write(line)


if __name__ == "__main__":
    main(sys.argv[1:])
