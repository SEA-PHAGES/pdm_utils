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

VALID_PIPELINES = set([
    "retrieve_data",
    "import",
    "cdd",
    "phamerate",
    "export",
    "compare"])
RUN_HELP = "Command line script to call a pdm_utils pipeline."
PIPELINE_HELP = \
    ("Name of the pdm_utils pipeline to run: "
     "(retrieve_data, import, cdd, phamerate, export")
DATABASE_HELP = \
    "Name of the MySQL database to import the genomes."
INPUT_FOLDER_HELP = \
    ("Path to the folder containing files to be processed.")
OUTPUT_FOLDER_HELP = \
    ("Path to the folder where processed data and files can be stored.")
IMPORT_GENOME_FOLDER_HELP = \
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
FILENAME_FLAG_HELP = \
    ("Indicates whether the filename_flag should be used "
     "to identify the genome during import.")
TEST_RUN_HELP = \
    ("Indicates whether the script should make any changes to the database. "
     "If False, the production run will implement all changes in the "
     "indicated database. If True, the test run will not "
     "implement any changes")
RUN_MODE_HELP = \
    ("Indicates the evaluation configuration "
     "for importing genomes.")
DESCRIPTION_FIELD_HELP = \
    ("Indicates the field in CDS features that is expected "
     "to store the gene description.")

def create_parser():
    """Construct the args that need to be parsed."""
    parser = argparse.ArgumentParser(description=RUN_HELP)
    # Required for all pipelines:
    parser.add_argument("pipeline", type=str,
                        choices=list(VALID_PIPELINES),
                        help=PIPELINE_HELP)
    # Commonly used args:
    parser.add_argument("-db", "--database", type=str,
        help=DATABASE_HELP)
    parser.add_argument("-if", "--input_folder", type=os.path.abspath,
        help=INPUT_FOLDER_HELP)
    parser.add_argument("-of", "--output_folder", type=os.path.abspath,
        help=OUTPUT_FOLDER_HELP)
    # Specific to import pipeline:
    parser.add_argument("-it", "--import_table", type=os.path.abspath,
        help=IMPORT_TABLE_HELP)
    parser.add_argument("-ff", "--filename_flag", action="store_true",
        default=False, help=FILENAME_FLAG_HELP)
    parser.add_argument("-tr", "--test_run", action="store_false", default=True,
        help=TEST_RUN_HELP)
    parser.add_argument("-rm", "--run_mode", type=str.lower,
        choices=list(constants.RUN_MODES.keys()), default="phagesdb",
        help=RUN_MODE_HELP)
    parser.add_argument("-df", "--description_field", default="product",
        choices=list(constants.DESCRIPTION_FIELD_SET),
        help=DESCRIPTION_FIELD_HELP)
    return parser

def get_args(input_list):
    """Adds a 'sql_handle' arg attribute after processing user input.

    The sql_handle args attribute should not be directly accessible
    at the command line. It is only changed depending on
    get_sql_handle(), which is called based on the --database option.
    Arguments provided by the user are first processed, then the
    'sql_handle' argument is added.
    """
    parser = create_parser()
    args = parser.parse_args(args=input_list)
    parser.add_argument("--sql_handle")
    parser.parse_args(args=["--sql_handle", None], namespace=args)
    return args


# TODO this may need to be moved to phamerator module.
# TODO unittest.
def get_sql_handle(database):
    """Set up MySQL credentials."""
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = database
    sql_handle.open_connection()
    if (not sql_handle.credential_status or not sql_handle._database_status):
        return None
    else:
        return sql_handle

# TODO unittest.
def main():
    """Verify a valid pipeline is selected and arguments provided are valid.

    The command line arguments are parsed and performs several basic
    checks on the arguments. Then they are passed to sub-functions to
    specifically validate the arguments based on the selected pipeline.
    """
    args = get_args(sys.argv)

    # General validation of arguments:
    # pipeline: validated automatically with choices.

    # --database: Confirm MySQL creds.
    if args.database is not None:
        args.sql_handle = get_sql_handle(args.database)
        if args.sql_handle == None:
            print("No connection to the selected database.")
            sys.exit(1)

    # --input_folder: Confirm exits.
    if args.input_folder is not None:
        valid_folder = basic.verify_path(args.input_folder, "dir")
        if not valid_folder:
            print("Invalid input folder.")
            sys.exit(1)

    # --import_table: Confirm exists.
    if args.import_table is not None:
        valid_import_table = basic.verify_path(args.import_table, "file")
        if not valid_import_table:
            print("Invalid import table file.")
            sys.exit(1)

    # --run_mode: validated automatically with choices.
    # TODO This can also be a filename if custom?

    # Specific validation of arguments based on selected pipeline:
    if args.pipeline == "retrieve_data":
        run_retrieve_data(args)
    elif args.pipeline == "import":
        run_import(args)
    elif args.pipeline == "cdd":
        run_cdd(args)
    elif args.pipeline == "phamerate":
        run_phamerate(args)
    elif args.pipeline == "export":
        run_export(args)
    else:
        run_compare_dbs(args)
    print("Pipeline completed")

def run_import(args):
    """Verify the correct arguments are selected for import new genomes."""
    if args.input_folder is None:
        print(IMPORT_GENOME_FOLDER_HELP)
        sys.exit(1)
    if args.import_table is None:
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
    import_genome.import_io(sql_handle=args.sql_handle,
        genome_folder=args.input_folder, import_table_file=args.import_table,
        filename_flag=args.filename_flag, test_run=args.test_run,
        description_field=args.description_field, run_mode=args.run_mode)

# TODO implement.
def run_retrieve_data(args):
    """Validate arguments for retrieving new data for import."""
    pass

# TODO implement.
def run_cdd(args):
    """Validate arguments for evaluating conserved domain data."""
    pass

# TODO implement.
def run_phamerate(args):
    """Validate arguments for grouping gene products into phamilies."""
    pass

# TODO implement.
def run_export(args):
    """Validate arguments for exporting data from the database."""
    pass

# TODO implement.
def run_compare_dbs(args):
    """Validate arguments for comparing databases."""
    pass



if __name__ == "__main__":
    main()
