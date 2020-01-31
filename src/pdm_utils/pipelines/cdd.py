import argparse
import os
import platform
import shlex
from subprocess import Popen, PIPE

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbirpsblastCommandline

from pdm_utils.functions import mysqldb
from pdm_utils.functions.basic import expand_path
from pdm_utils.functions.parallelize import *

# SQL QUERIES
GET_GENES_FOR_CDD = (
    "SELECT GeneID, Translation FROM gene WHERE DomainStatus = 0")
GET_UNIQUE_HIT_IDS = "SELECT HitID FROM domain"

# SQL COMMANDS
INSERT_INTO_DOMAIN = """INSERT IGNORE INTO domain (HitID, DomainID, Name, Description) VALUES ("{}", "{}", "{}", "{}")"""
INSERT_INTO_GENE_DOMAIN = """INSERT IGNORE INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("{}", "{}", {}, {}, {})"""
UPDATE_GENE = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"


def setup_argparser():
    """
    Builds argparse.ArgumentParser for this script
    :return:
    """
    # Pipeline description
    description = (
        "Uses rpsblast to search the NCBI conserved domain database "
        "for significant domain hits in all new proteins of a "
        "MySQL database.")

    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("db", type=str,
                        help="name of database to phamerate")
    parser.add_argument("dir", type=str,
                        help="path to local directory housing cdd database")
    parser.add_argument("--threads", default=mp.cpu_count(), type=int,
                        help="number of concurrent cdd searches to run")
    parser.add_argument("--evalue", default=0.001, type=float,
                        help="evalue cutoff for rpsblast hits")
    parser.add_argument("--tmp_dir", default="/tmp/cdd", type=str,
                        help="path to temporary directory for file I/O")
    parser.add_argument("--rpsblast", default="", type=str,
                        help="path to rpsblast(+) binary")
    return parser


def make_tempdir(tmp_dir):
    """
    Uses pdm_utils.functions.basic.expand_path to expand TMP_DIR; then
    checks whether tmpdir exists - if it doesn't, uses os.makedirs to
    make it recursively.
    :param tmp_dir: location where I/O should take place
    :return:
    """
    try:
        path = expand_path(tmp_dir)
        # If the path doesn't exist yet
        if not os.path.exists(path):
            # Make it recursively
            os.makedirs(path)
    except OSError as err:
        print(f"Error {err.args[0]}: {err.args[1]}")


def search_and_process(rpsblast, cdd_name, tmp_dir, evalue,
                       geneid, translation):
    """
    Uses rpsblast to search indicated gene against the indicated cdd
    :param rpsblast: path to rpsblast binary
    :param cdd_name: cdd database path/name
    :param tmp_dir: path to directory where I/O will take place
    :param evalue: evalue cutoff for rpsblast
    :param geneid: name of the gene to query
    :param translation: protein sequence for gene to query
    :return: results
    """
    # Setup I/O variables
    i = "{}/{}.txt".format(tmp_dir, geneid)
    o = "{}/{}.xml".format(tmp_dir, geneid)

    # Write the input file
    with open(i, "w") as fh:
        fh.write(">{}\n{}".format(geneid, translation))

    # Setup and run the rpsblast command
    rps_command = NcbirpsblastCommandline(cmd=rpsblast, db=cdd_name,
                                          query=i, out=o, outfmt=5,
                                          evalue=evalue)
    rps_command()

    # Process results into a single list
    results = []

    with open(o, "r") as fh:
        for record in NCBIXML.parse(fh):
            # Only need to process if there are record alignments
            if record.alignments:
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect <= evalue:
                            align.hit_def = align.hit_def.replace("\"", "\'")

                            des_list = align.hit_def.split(",")
                            if len(des_list) == 1:
                                description = des_list[0].strip()
                                domain_id = None
                                name = None
                            elif len(des_list) == 2:
                                domain_id = des_list[0].strip()
                                description = des_list[1].strip()
                                name = None
                            else:
                                domain_id = des_list[0].strip()
                                name = des_list[1].strip()
                                description = ",".join(des_list[2:]).strip()

                            # Try to put domain into domain table
                            results.append(INSERT_INTO_DOMAIN.format(
                                align.hit_id, domain_id, name, description))

                            # Try to put this hit into gene_domain table
                            results.append(INSERT_INTO_GENE_DOMAIN.format(
                                geneid, align.hit_id, float(hsp.expect),
                                int(hsp.query_start), int(hsp.query_end)))

    # Update this gene's DomainStatus to 1
    results.append(UPDATE_GENE.format(geneid))
    return results


def learn_cdd_name(cdd_dir):
    cdd_files = os.listdir(cdd_dir)
    cdd_files = [os.path.join(cdd_dir, x.split(".")[0]) for x in cdd_files]
    file_set = set(cdd_files)
    if len(file_set) == 1:
        cdd_name = cdd_files[0]
    else:
        cdd_name = ""
    return cdd_name


def main(argument_list):
    """
    :param argument_list:
    :return:
    """
    # Setup argument parser
    cdd_parser = setup_argparser()

    # Use argument parser to parse argument_list
    args = cdd_parser.parse_args(argument_list)

    # Store arguments in more easily accessible variables
    database = args.db
    cdd_dir = expand_path(args.dir)
    cdd_name = learn_cdd_name(cdd_dir)
    threads = args.threads
    evalue = args.evalue
    rpsblast = args.rpsblast
    tmp_dir = args.tmp_dir

    # Early exit if either 1) cdd_name == "" or 2) no rpsblast given and we are
    # unable to find one
    if cdd_name == "":
        print(f"Unable to learn cdd database name. Make sure the files in "
              f"{cdd_dir} all have the same basename.")
        return

    if rpsblast == "":
        # See if we're running on a Mac
        if platform.system() == "Darwin":
            print("Detected MacOS operating system...")
            command = shlex.split("which rpsblast")
        # Otherwise see if we're on a Linux machine
        elif platform.system() == "Linux":
            print("Detected Linux operating system...")
            command = shlex.split("which rpsblast+")
        # Windows or others - unsupported, leave early
        else:
            print(f"Unsupported system '{platform.system()}'; cannot run "
                  f"cdd pipeline.")
            return

        # If we didn't return, we have a command.
        # Run it, and PIPE stdout into rpsblast_path
        with Popen(args=command, stdout=PIPE) as proc:
            rpsblast = proc.stdout.read().decode("utf-8").rstrip("\n")

        # If empty string, rpsblast not found in globally
        # available executables, otherwise proceed with value.
        if rpsblast == "":
            print("No rpsblast binary found. "
                  "If you have rpsblast on your machine, please try "
                  "again with the '--rpsblast' flag and provide the "
                  "full path to your rpsblast binary.")
            return
        else:
            print(f"Found rpsblast binary at '{rpsblast}'...")

    engine = mysqldb.connect_to_db(database)
    result = engine.execute(GET_GENES_FOR_CDD)
    cdd_genes = []
    for row in result:
        row_as_dict = dict(row)
        cdd_genes.append(row_as_dict)
    engine.dispose()

    # Print number of genes to process
    print(f"{len(cdd_genes)} genes to search for conserved domains...")

    # Only run the pipeline if there are genes returned that need it
    if len(cdd_genes) > 0:
        # Create temp_dir
        make_tempdir(tmp_dir)

        # Build jobs list
        jobs = []
        for cdd_gene in cdd_genes:
            jobs.append((rpsblast, cdd_name, tmp_dir, evalue,
                         cdd_gene["GeneID"], cdd_gene["Translation"]))

        results = parallelize(jobs, threads, search_and_process)
        print("\n")

        for result in results:
            mysqldb.execute_transaction(engine, result)
        engine.dispose()

    return
