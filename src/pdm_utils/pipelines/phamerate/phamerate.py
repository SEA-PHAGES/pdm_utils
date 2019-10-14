"""This script allows user to "phamerate" a MySQL database compliant
with the current Phamerator database schema using either blastclust OR
mmseqs2 (default).
"""
import argparse
import sys
import datetime

try:
    import pymysql as pms
except ImportError:
    print("Failed to import pymysql. Please install it and try again.")
    sys.exit(1)

from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from pdm_utils.classes.phameratorhandler import PhameratorHandler


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
    parser.add_argument("--use_blast", action='store_true', default=False,
                        help="use blastclust instead of MMseqs2")
    parser.add_argument("--threads", default=1, help="number of threads to use")
    parser.add_argument("--identity", type=int, default=40,
                        help="percent identity threshold in range [0,100]")
    parser.add_argument("--coverage", type=int, default=80,
                        help="coverage threshold in range [0,100]")
    parser.add_argument("--steps", type=int, default=1,
                        help="number of clustering steps (mmseqs only)")
    parser.add_argument("--max_seqs", type=int, default=500,
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
    parser.add_argument("--no_filter_dups", action='store_true', default=False,
                        help="override filtering of duplicate genes - slower")
    return parser


def main(argument_list):
    # Set up the argument parser
    phamerate_parser = setup_argparser()

    # Parse arguments
    args = phamerate_parser.parse_args(argument_list)

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

    # Initialize PhameratorHandler with our mysql_handler and run main
    pham_handler = PhameratorHandler(mysql_handler, args)

    # Filter redundant translations from phameration (insignificant
    # improvement for MMseqs2 but HUGE improvement for blastclust) - default
    # is to do it
    if args.no_filter_dups is True:
        pham_handler.filter_redundant = False

    pham_handler.main()

    # Record stop time
    stop = datetime.datetime.now()

    # Report phameration elapsed time
    print("Elapsed time: {}".format(str(stop - start)))


if __name__ == "__main__":
    main(sys.argv[1:])
