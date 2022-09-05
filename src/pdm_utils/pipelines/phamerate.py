"""Program to group related gene products into phamilies using MMseqs2
for similarity search and clustering."""

import argparse
from datetime import datetime
import pathlib

from phammseqs import assemble_phams, SequenceDB
from phammseqs.multiprocess import CPUS

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions.configfile import *
from pdm_utils.functions.phameration import *
from pdm_utils.functions.multiprocess import *

DATE = datetime.now().strftime("%d_%b_%Y")
OUT_DIR = pathlib.Path().cwd().joinpath(f"phamerate__{DATE}")

# Sequence-sequence
SEQ_MIN_SEQ_ID = 35.0
SEQ_COVERAGE = 80.0
SEQ_EVALUE = 0.001
SEQ_SENSITIVITY = 7
SEQ_CLUSTER_MODE = 0
SEQ_CLUSTER_STEPS = 1

HMM_MIN_SEQ_ID = 15.0
HMM_COVERAGE = 70.0
HMM_EVALUE = 0.001
HMM_SENSITIVITY = 7
HMM_CLUSTER_MODE = 0
HMM_CLUSTER_STEPS = 3

EPILOG = """
Steinegger M. and Söding J. MMseqs2 enables sensitive protein
sequence searching for the analysis of massive data sets. Nature
Biotechnology, 2017. doi: 10.1038/nbt.3988"""

DB_SUMMARY = """
Database summary:
=============================
 {} total phams
 {} orphams
 {} genes in phams
 {} genes not in phams
 {} total genes
 {} non-redundant genes
=============================
"""


def setup_argparser():
    """Parse command line arguments."""
    p = argparse.ArgumentParser(description=__doc__, epilog=EPILOG)

    p.add_argument("db", type=str, help="name of database to phamerate")

    s = p.add_argument_group("MMseqs2 sequence-sequence clustering arguments")
    s.add_argument("--identity",
                   type=float, default=SEQ_MIN_SEQ_ID, metavar='',
                   help=f"percent identity for sequence-sequence clustering "
                        f"[default: {SEQ_MIN_SEQ_ID}%%]")
    s.add_argument("--coverage",
                   type=float, default=SEQ_COVERAGE, metavar='',
                   help=f"percent coverage for sequence-sequence clustering "
                        f"[default: {SEQ_COVERAGE}%%]")
    s.add_argument("--evalue",
                   type=float, default=SEQ_EVALUE, metavar='',
                   help=f"E-value threshold for sequence-sequence clustering "
                        f"[default: {SEQ_EVALUE}]")
    s.add_argument("--sensitivity",
                   type=float, default=SEQ_SENSITIVITY, metavar='',
                   help=f"sensitivity: 1 favors speed, 7 favors "
                        f"sensitivity [default: {SEQ_SENSITIVITY}]")
    s.add_argument("--cluster-mode",
                   type=int, default=SEQ_CLUSTER_MODE, metavar='',
                   help=f"clustering algorithm [default: {SEQ_CLUSTER_MODE}]")
    s.add_argument("--cluster-steps",
                   type=int, default=SEQ_CLUSTER_STEPS, metavar='',
                   help=f"number of steps for sequence-sequence clustering "
                        f"to proceed in [default: {SEQ_CLUSTER_STEPS}]")

    h = p.add_argument_group("MMseqs2 profile-sequence clustering arguments")
    h.add_argument("--hmm-identity",
                   type=float, default=HMM_MIN_SEQ_ID, metavar='',
                   help=f"percent identity for profile-consensus clustering "
                        f"[default: {HMM_MIN_SEQ_ID}%%]")
    h.add_argument("--hmm-coverage",
                   type=float, default=HMM_COVERAGE, metavar='',
                   help=f"percent coverage for profile-consensus clustering "
                        f"[default: {HMM_COVERAGE}%%]")
    h.add_argument("--hmm-evalue",
                   type=float, default=HMM_EVALUE, metavar='',
                   help=f"E-value threshold for profile-consensus clustering "
                        f"[default: {HMM_EVALUE}]")
    h.add_argument("--hmm-sensitivity",
                   type=float, default=HMM_SENSITIVITY, metavar='',
                   help=f"sensitivity: 1 favors speed, 7 favors "
                        f"sensitivity [default: {HMM_SENSITIVITY}]")
    h.add_argument("--hmm-cluster-mode",
                   type=int, default=HMM_CLUSTER_MODE, metavar='',
                   help=f"clustering algorithm [default: {HMM_CLUSTER_MODE}]")
    h.add_argument("--hmm-cluster-steps",
                   type=int, default=HMM_CLUSTER_STEPS, metavar='',
                   help=f"number of steps for profile-consensus clustering "
                        f"to proceed in [default: {HMM_CLUSTER_STEPS}]")
    h.add_argument("--skip-hmm", action="store_true",
                   help="do not perform profile-consensus clustering")

    p.add_argument("-v", "--verbose", action="store_true",
                   help="print progress messages to the console")
    p.add_argument("-d", "--debug", action="store_true",
                   help="run in debug mode")
    p.add_argument("-x", "--reset", action="store_true",
                   help="reset existing pham data before running")
    p.add_argument("-t", "--threads", type=int, default=CPUS,
                   help=f"number of threads to use [default: {CPUS}]")
    p.add_argument("-c", "--config-file", type=pathlib.Path, default=None,
                   help="path to file containing MySQL login details")

    return p


def setup_mmseqs_params(args):
    """Parse MMseqs2-related arguments into dictionaries of clustering
    parameters.

    :param args: dictionary of commandline arguments
    :type args: dict[str, int or float]
    :return: seq_params, hmm_params
    """
    seq_params = {"identity":         args["identity"],
                  "coverage":         args["coverage"],
                  "evalue":           args["evalue"],
                  "sensitivity":      args["sensitivity"],
                  "cluster_mode":     args["cluster_mode"],
                  "cluster_steps":    args["cluster_steps"]}

    if args["skip_hmm"]:
        hmm_params = None
    else:
        hmm_params = {"identity":      args["hmm_identity"],
                      "coverage":      args["hmm_coverage"],
                      "evalue":        args["hmm_evalue"],
                      "sensitivity":   args["sensitivity"],
                      "cluster_mode":  args["hmm_cluster_mode"],
                      "cluster_steps": args["hmm_cluster_steps"]}

    return seq_params, hmm_params


def main(argument_list):
    """Commandline entrypoint to this module."""
    # Set up the argument parser
    parser = setup_argparser()

    # Parse arguments
    args = vars(parser.parse_args(argument_list[2:]))

    # Create config object with data obtained from file and/or defaults.
    config = build_complete_config(args["config_file"])
    mysql_creds = config["mysql"]

    # Parse MMseqs2 parameters
    seq_params, hmm_params = setup_mmseqs_params(args)

    # Other runtime params
    verbose = args["verbose"]
    debug = args["debug"]
    threads = args["threads"]
    reset = args["reset"]

    # Record start time
    start_time = datetime.now()

    # Initialize SQLAlchemy engine with database provided at CLI
    alchemist = AlchemyHandler(database=args["db"],
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(login_attempts=5, pipeline=True)
    engine = alchemist.engine

    # If reset, execute pham reset statements
    if reset:
        engine.execute("UPDATE gene SET PhamID = NULL")
        engine.execute("DELETE FROM pham")

    # Get old pham data and un-phamerated genes
    old_phams = get_pham_geneids(engine)
    old_colors = get_pham_colors(engine)
    new_genes = get_new_geneids(engine)

    # Create a SequenceDB
    database = SequenceDB()
    for geneid, translation in get_geneids_and_translations(engine).items():
        database.add_gene(geneid, translation)

    # Print initial snapshot of MySQL database
    if verbose:
        initial_summary = \
            DB_SUMMARY.format(len(old_phams),
                              sum([len(x) == 1 for x in old_phams.values()]),
                              sum([len(x) for x in old_phams.values()]),
                              len(new_genes), len(database),
                              len(database.translations))
        print(initial_summary)

    phams = assemble_phams(database, seq_params, hmm_params,
                           cpus=threads, verbose=verbose, debug=debug)

    new_phams = dict()
    for i, pham in enumerate(phams):
        pham_geneids = set(pham.geneids)
        new_phams[i+1] = pham_geneids

    # Preserve old pham names and colors
    print("Preserving old phamily names/colors where possible...")
    new_phams, new_colors = preserve_phams(old_phams, new_phams,
                                           old_colors, new_genes)

    # Early exit if we don't have new phams or new colors - avoids
    # overwriting the existing pham data with incomplete new data
    if len(new_phams) == 0 or len(new_colors) == 0:
        print("Failed to parse new pham/color data properly... Terminating "
              "pipeline")
        return

    # Update gene/pham tables with new pham data. Pham colors need to be done
    # first, because gene.PhamID is a foreign key to pham.PhamID.
    print("Updating pham data in database...")
    update_pham_table(new_colors, engine)
    update_gene_table(new_phams, engine)

    # Fix miscolored phams/orphams
    print("Phixing phalsely phlagged orphams...", end=" ")
    fix_white_phams(engine)
    print("Phixing phalsely hued phams...", end=" ")
    fix_colored_orphams(engine)

    # Close all connections in the connection pool.
    engine.dispose()

    # Print final snapshot of MySQL database
    if verbose:
        final_summary = \
            DB_SUMMARY.format(len(phams),
                              sum([len(x) == 1 for x in phams]),
                              sum([len(x) for x in phams]), 0,
                              len(database), len(database.translations))
        print(final_summary)

    # Record stop time
    stop_time = datetime.now()
    elapsed_time = str(stop_time - start_time)

    # Report phameration elapsed time
    print(f"Elapsed time: {elapsed_time}")


if __name__ == "__main__":
    main(sys.argv[1:])
