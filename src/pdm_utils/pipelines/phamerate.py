"""Primary pipeline to group related gene products into phamilies."""

import argparse
import datetime
import os
import shutil
from pdm_utils.functions.phameration import *
from pdm_utils.functions import mysqldb

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
    parser.add_argument("--identity", type=float, default=32.5,
                        help="percent identity threshold in range [0,100]")
    parser.add_argument("--coverage", type=float, default=65,
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
    parser.add_argument("--temp_dir", type=str, default="/tmp/phamerate",
                        help="temporary directory for phameration file I/O")
    return parser


def main(argument_list):
    # Set up the argument parser
    phamerate_parser = setup_argparser()

    # Parse arguments
    args = phamerate_parser.parse_args(argument_list)
    program = args.program
    temp_dir = args.temp_dir

    # Initialize SQLAlchemy engine with database provided at CLI
    engine = mysqldb.connect_to_db(args.db)

    # If we made it past the above connection_status() check, database access
    # works (user at least has SELECT privileges on the indicated database).
    # We'll assume that they also have UPDATE, INSERT, and TRUNCATE privileges.

    # Record start time
    start = datetime.datetime.now()

    # Refresh temp_dir
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except OSError:
            print(f"Failed to delete existing temp directory '{temp_dir}'")
            return
    try:
        os.makedirs(temp_dir)
    except OSError:
        print(f"Failed to create new temp directory '{temp_dir}")
        return

    # Get old pham data and un-phamerated genes
    old_phams = get_pham_geneids(engine)
    old_colors = get_pham_colors(engine)
    unphamerated = get_new_geneids(engine)

    # Get GeneIDs & translations, and translation groups
    genes_and_trans = map_geneids_to_translations(engine)
    translation_groups = map_translations_to_geneids(engine)

    # Write input fasta file
    write_fasta(translation_groups, temp_dir)

    # Create clusterdb and perform clustering
    program_params = get_program_params(program, args)
    create_clusterdb(program, temp_dir)
    phamerate(program_params, program, temp_dir)

    # Parse phameration output
    new_phams = parse_output(program, temp_dir)
    new_phams = reintroduce_duplicates(new_phams, translation_groups, genes_and_trans)

    # Preserve old pham names and colors
    new_phams, new_colors = preserve_phams(old_phams, new_phams, old_colors,
                                           unphamerated)

    # Early exit if we don't have new phams or new colors - avoids
    # overwriting the existing pham data with potentially incomplete new data
    if len(new_phams) == 0 or len(new_colors) == 0:
        print("Failed to parse new pham/color data... Terminating pipeline")
        return

    # If we got past the early exit, we are probably safe to truncate the
    # pham table, and insert the new pham data
    # Clear old pham data - auto commits at end of transaction - this will also
    # set all PhamID values in gene table to NULL
    commands = ["DELETE FROM pham"]
    mysqldb.execute_transaction(engine, commands)


    # Insert new pham/color data
    reinsert_pham_data(new_phams, new_colors, engine)

    # Fix miscolored phams/orphams
    fix_miscolored_phams(engine)

    # Close all connections in the connection pool.
    engine.dispose()

    # Record stop time
    stop = datetime.datetime.now()

    # Report phameration elapsed time
    print("Elapsed time: {}".format(str(stop - start)))
