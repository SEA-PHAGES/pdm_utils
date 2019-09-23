"""Use this script to run all pipelines within the pipelines folder.
It validates the command line arguments selected,
verifies the selected arguments match those required for the selected
pipeline, then passes that information to the main pipeline module.
"""
import sys
import os
import argparse
from pdm_utils.classes import mysqlconnectionhandler
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.pipelines.db_import import import_genome


VALID_PIPELINES = {
    "retrieve_data":"run_retrieve_data",
    "import":"run_import",
    "cdd":"run_cdd",
    "phamerate":"run_phamerate",
    "export":"run_export"
    }


PIPELINE_HELP = \
    ("Name of the pdm_utils pipeline to run: "
     "(retrieve_data, import, cdd, phamerate, export")

DATABASE_HELP = \
    "Name of the MySQL database to import the genomes."

GENOME_FOLDER_HELP = \
    ("Path to the folder containing GenBank-formatted "
     "flat files to be processed.")

IMPORT_TABLE_HELP = \
    """
    Path to the CSV-formatted table containing
    instructions to process each genome.
    Structure of import ticket table:
        1. Action to implement on the database (add, remove)
        2. PhageID
        3. Host genus
        4. Cluster
        5. Subcluster
        6. Annotation status (draft, final, unknown)
        7. Annotation authorship (hatfull, gbk)
        8. Gene description field (product, note, function)
        9. Accession
        10. Run mode
    """

FILENAME_HELP = \
    ("Indicates whether the filename should be used "
     "to identify the genome during import.")

TEST_RUN_HELP = \
    ("Indicates whether the import run is for ",
     "test or production purposes. "
     "A production run will implement all changes in the "
     "indicated database. A test run will not "
     "implement any changes")

RUN_MODE_HELP = \
    ("Indicates the evaluation configuration "
     "for importing genomes.")

DESCRIPTION_FIELD_HELP = \
    ("Indicates the field in CDS features that is expected "
     "to store the gene description.")




def create_parser():
    """Construct the args that need to be parsed."""

    parser = argparse.ArgumentParser(
        description="Top level command line script to call pdm_utils pipelines.")

    # Required for all pipelines.
    parser.add_argument("-p", "--pipeline", type=str, default=False,
        required=True, choices=list(VALID_PIPELINES.keys()),
        help=PIPELINE_HELP)

    # Commonly used args.
    parser.add_argument("-db", "--database", type=str,
        help=DATABASE_HELP)

    # Specific to import pipeline:
    parser.add_argument("-g", "--genome_folder", type=os.path.abspath,
        help=GENOME_FOLDER_HELP)
    parser.add_argument("-t", "--table", type=os.path.abspath,
        help=IMPORT_TABLE_HELP)
    parser.add_argument("-f", "--filename", action="store_true", default=False,
        help=FILENAME_HELP)
    parser.add_argument("-tr", "--testrun", action="store_false", default=True,
        help=TEST_RUN_HELP)
    parser.add_argument("-rm", "--run_mode", type=str.lower,
        choices=list(constants.RUN_MODES.keys()), default="phagesdb",
        help=RUN_MODE_HELP)
    parser.add_argument("-df", "--description_field", default="product",
        choices=list(constants.DESCRIPTION_FIELD_SET),
        help=DESCRIPTION_FIELD_HELP)
    return parser


def get_sql_handle(database):
    """Set up MySQl credentials."""
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        print("\nUnable to connect to the database")
        sys.exit(1)
    else:
        return sql_handle

def main():
    """Verify a valid pipeline is selected and arguments provided are valid."""

    parser = create_parser()
    args = parser.parse_args(sys.argv)


    # General validation of arguments:

    # --pipeline: validated automatically with choices.

    # --genome_folder: Confirm exits.
    if args.genome_folder is not None:
        if not basic.verify_path(args.genome_folder, "dir"):
            print("\n\nInvalid input for genome folder.\n\n")
            sys.exit(1)

    # --table: Confirm exists.
    if args.table is not None:
        if not basic.verify_path(args.table, "file"):
            print("\n\nInvalid input for import table file.\n\n")
            sys.exit(1)

    # --database: Confirm MySQL creds.
    if args.database is not False:
        sql_handle = get_sql_handle(args.database)

    # --run_mode: validated automatically with choices.
    # TODO This can also be a filename if custom?

    # Specific validation of arguments based on selected pipeline:
    called_pipeline = VALID_PIPELINES[args.pipeline]
    called_pipeline(args, sql_handle)
    print("Pipeline completed")

def run_import(args, sql_handle):
    """Verify the correct arguments are selected for import new genomes."""
    if args.genome_folder is None:
        print(GENOME_FOLDER_HELP)
        sys.exit(1)
    if args.table is None:
        print(IMPORT_TABLE_HELP)
        sys.exit(1)
    if args.database is None:
        print(DATABASE_HELP)
        sys.exit(1)

    # TODO build functionality to check run_mode option, and allow custom
    # run mode file to be used.
    if args.run_mode is None:
        print(RUN_MODE_HELP)
        sys.exit(1)


    # If everything checks out, pass args to the main import pipeline:
    import_genome.import_io(sql_handle=sql_handle,
        genome_folder=args.genome_folder, import_table_file=args.table,
        filename_flag=args.filename, test_run=args.testrun,
        description_field=args.description_field, run_mode=args.run_mode)

def run_retrieve_data(args, sql_handle):
    pass
def run_cdd(args, sql_handle):
    pass
def run_phamerate(args, sql_handle):
    pass
def run_export(args, sql_handle):
    pass



if __name__ == "__main__":
    main()
