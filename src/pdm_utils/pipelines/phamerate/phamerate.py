import argparse
import sys
import datetime
import os
from shutil import rmtree

try:
    import pymysql as pms
except ImportError:
    print("Failed to import pymysql. Please install it and try again.")
    sys.exit(1)

from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from pdm_utils.functions.phameration import *


def setup_argparser():
    """
    Builds argparse.ArgumentParser for this script
    :return:
    """
    # Pipeline description
    description = "Groups related CDS features into phamilies using MMseqs2 " \
                  "(default) or blastclust"

    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("db", type=str,
                        help="name of database to phamerate")
    parser.add_argument("--program", type=str, default="mmseqs", choices=[
        "mmseqs", "blast"], help="which program to use for clustering")
    parser.add_argument("--threads", default=1, help="number of threads to use")
    parser.add_argument("--identity", type=int, default=40,
                        help="percent identity threshold in range [0,100]")
    parser.add_argument("--coverage", type=int, default=80,
                        help="coverage threshold in range [0,100]")
    parser.add_argument("--steps", type=int, default=1,
                        help="number of clustering steps (mmseqs only)")
    parser.add_argument("--max_seqs", type=int, default=1000,
                        help="max number of targets per query per step "
                             "(mmseqs only)")
    parser.add_argument("--verbose", type=int, default=3,
                        help="verbosity of output in range [0, 3] (mmseqs "
                             "only)")
    parser.add_argument("--aln_mode", type=int, default=3,
                        help="alignment mode in range [0, 4] (mmseqs only)")
    parser.add_argument("--cov_mode", type=int, default=0,
                        help="coverage mode in range [0, 4] (mmseqs only)")
    parser.add_argument("--clu_mode", type=int, default=0,
                        help="cluster mode in range [0, 3] (mmseqs only)")
    return parser


def main(argument_list):
    # Set up the argument parser
    phamerate_parser = setup_argparser()

    # Parse arguments
    args = phamerate_parser.parse_args(argument_list)
    program = args.program

    # Initialize MySQLConnectionHandler with database provided at CLI
    mysql_handler = MySQLConnectionHandler()
    mysql_handler.database = args.db

    # Use open_connection() method to simultaneously get user credentials
    # and test database access
    mysql_handler.open_connection()

    # Handle possibility that connection failed (bad database or user/pass)
    if mysql_handler.connection_status() is not True:
        sys.exit(1)

    # If we made it past the above connection_status() check, database access
    # works (user at least has SELECT privileges on the indicated database).
    # We'll assume that they also have UPDATE, INSERT, and TRUNCATE privileges.

    # Record start time
    start = datetime.datetime.now()

    # Refresh temp_dir
    if os.path.exists("/tmp/phamerate"):
        try:
            rmtree("/tmp/phamerate")
        except OSError:
            print("Failed to delete existing '/tmp/phamerate'")
            sys.exit(1)
    try:
        os.makedirs("/tmp/phamerate")
    except OSError:
        print("Failed to create new '/tmp/phamerate'")
        sys.exit(1)

    # Get old pham data and un-phamerated genes
    old_phams, old_colors = read_existing_phams(mysql_handler)
    unphamerated = read_unphamerated_genes(mysql_handler)

    # Get GeneIDs & translations, and translation groups
    genes_and_trans, trans_groups = get_translations(mysql_handler)

    # Write input fasta file
    write_fasta(trans_groups)

    # Create clusterdb and perform clustering
    program_params = get_program_params(program, args)
    create_clusterdb(program)
    phamerate(program_params, program)

    # Convert output to parseable (mmseqs) and parse outputs
    convert_to_parseable(program)
    new_phams = parse_output(program)
    new_phams = reintroduce_duplicates(new_phams, trans_groups, genes_and_trans)

    # Preserve old pham names and colors
    new_phams, new_colors = preserve_phams(old_phams, new_phams, old_colors,
                                           unphamerated)

    # Early exit if we don't have new phams or new colors - avoids
    # overwriting the existing pham data with potentially incomplete new data
    if len(new_phams) == 0 or len(new_colors) == 0:
        print("Failed to parse new pham/color data... Terminating pipeline")
        return

    # If we got past the early exit, we are probably safe to truncate the
    # pham and pham_color tables, and insert the new pham data
    # Clear old pham data - auto commits at end of transaction
    commands = ["TRUNCATE TABLE pham", "TRUNCATE TABLE pham_color"]
    mysql_handler.execute_transaction(commands)

    # Insert new pham/color data
    reinsert_pham_data(new_phams, new_colors, mysql_handler)

    # Fix miscolored phams/orphams
    fix_miscolored_phams(mysql_handler)

    # Record stop time
    stop = datetime.datetime.now()

    # Report phameration elapsed time
    print("Elapsed time: {}".format(str(stop - start)))


if __name__ == "__main__":
    main(sys.argv[1:])
