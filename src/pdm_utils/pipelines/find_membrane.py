import argparse
import shutil
import sys
import time
from pathlib import Path

from sqlalchemy import update

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb
from pdm_utils.functions import querying


def main(unparsed_args):
    """Runs the complete update pipeline."""
    args = parse_args(unparsed_args[2:])

    print("Connecting to the MySQL database...")

    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]
    alchemist = AlchemyHandler(database=args.database,
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    # mysqldb.check_schema_compatibility(engine, "")

    find_membrane(alchemist)


def find_membrane():
    pass



if __name__ == "__main__":
    main(sys.argv[1:])
