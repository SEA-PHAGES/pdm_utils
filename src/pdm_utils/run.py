"""Use this script to run all pipelines within the pipelines folder.
It validates the command line arguments selected,
verifies the selected arguments match those required for the selected
pipeline, then passes that information to the main pipeline module.
"""
import sys
import os
import argparse
from classes import mysqlconnectionhandler
from constants import constants
from functions import basic
from pipelines.db_import import import_main

def main():
    """Verify a valid pipeline is selected and arguments provided are valid."""
    parser = argparse.ArgumentParser(
        description="Top level command line script to call pdm_utils pipelines.")

    # Required for all pipelines.
    parser.add_argument("-p", "--pipeline", type=str, default=False,
        required=True,
        help=("Name of the pdm_utils pipeline to run:",
              " (retrieve_data, import, cdd, phamerate, export"))

    # Commonly used args.
    parser.add_argument("-db", "--database", type=str,
        help="Name of the MySQL database to import the genomes.")

    # Specific to import pipeline:
    parser.add_argument("-g", "--genome_folder", type=os.path.abspath,
        help=("Path to the folder containing GenBank-formatted "
              "flat files to be processed."))
    parser.add_argument("-t", "--table", type=os.path.abspath,
        help="Path to the CSV-formatted table containing instructions "
             "to process each genome.")
    parser.add_argument("-f", "--filename", action="store_true", default=False,
        help=("Indicates whether the filename should be used "
              "to identify the genome."))
    parser.add_argument("-tr", "--testrun", action="store_false", default=True,
        help=("Indicates whether the import run is for ",
              "test or production purposes. "
              "A production run will implement all changes in the "
              "indicated database. A test run will not implement any changes"))
    parser.add_argument("-rm", "--run_mode", type=str,
        help="Indicates the evaluation configuration for importing genomes.")
    parser.add_argument("-df", "--description_field", default="product",
        choices=list(constants.DESCRIPTION_FIELD_SET),
        help=("Indicates the field in CDS features that is expected "
              "to store the gene description."))

    args = parser.parse_args()


    # General validation of arguments:

    # --pipeline: Confirm valid.
    VALID_PIPELINES = set(["retrieve_data", "import", "cdd",
                           "phamerate", "export"])
    if args.pipeline not in VALID_PIPELINES:
        print("Invalid pipeline.")
        sys.exit(1)

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
        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.database = args.database
        sql_handle.open_connection()
        if (not sql_handle.credential_status or not sql_handle._database_status):
            print("\nUnable to connect to the database")
            sys.exit(1)


    # --run_mode: Confirm valid. This can be a filename or
    # one of several preset options.
    if args.run_mode is not None:
        if args.run_mode.lower() not in constants.RUN_MODES.keys():
            if not basic.verify_path(args.run_mode, "file"):
                print("\n\nInvalid input for run_mode file.\n\n")
                sys.exit(1)




    # Specific validation of arguments based on selected pipeline:

    # Import:
    if args.pipeline == "import":
        run_import(args, sql_handle)
    elif args.pipeline == "export":
        #run_export()
        pass
    else:
        print("Invalid pipeline")
        sys.exit(1)
    print("Pipeline completed")

def run_import(args, sql_handle):
    """Verify the correct arguments are selected for import new genomes."""
    if args.genome_folder is None:
        gf_description = """
            Use the -g option to indicate folder of genomes.
            The folder should contain GenBank-formatted flat files.
            """
        print(gf_description)
        sys.exit(1)
    if args.table is None:
        table_description = """
            Use the -t option to indicate csv-formatted import table.
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
        print(table_description)
        sys.exit(1)
    if args.database is None:
        print("Use the -db option to indicate the database.")
        sys.exit(1)

    # TODO build functionality to check run_mode option, and allow custom
    # run mode file to be used.
    if args.run_mode is None:
        print("Use the -rm option to indicate the run_mode.")
        sys.exit(1)


    # If everything checks out, pass args to the main import pipeline:
    pipelines.db_import.import_main.main(sql_handle=sql_handle,
        genome_folder=args.genome_folder, import_table=args.table,
        filename_flag=args.filename, test_run=args.testrun,
        description_field=args.description_field, run_mode=args.run_mode)



if __name__ == "__main__":
    main()
