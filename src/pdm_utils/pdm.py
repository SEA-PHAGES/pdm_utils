"""Use this script to run all pipelines within the pipelines folder.
"""
# This avoids problems wtih import scope.
import sys
import argparse

# from pipelines.db_import import test
from pdm_utils.pipelines.db_import import test


parser = argparse.ArgumentParser(
    description="test.")
parser.add_argument("-db", "--database", type=str, default=False,
    help="Name of the MySQL database to import the genomes.")
parser.add_argument("-p", "--pipeline", type=str, default=False,
    help="Name of the pdm_utils pipeline to run.")
args = parser.parse_args()
print(args)


if args.pipeline == 'testrun':
    # pipelines.db_import.test.main()
    test.main(args)

print("Finished")
