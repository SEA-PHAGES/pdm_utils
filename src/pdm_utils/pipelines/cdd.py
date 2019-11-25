import argparse
import sys
import os
import multiprocessing as mp

from Bio.Blast.Applications import NcbirpsblastCommandline
from Bio.Blast import NCBIXML

from pdm_utils.functions.basic import expand_path
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler

# CHANGE THESE VARIABLES TO MODIFY BEHAVIOR
RPSBLAST_BINARY = "/usr/bin/rpsblast+"
TMP_DIR = "/tmp/cdd"

# SQL QUERIES
GET_GENES_FOR_CDD = "SELECT GeneID, Translation FROM gene WHERE DomainStatus = 0"
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
    description = "Uses rpsblast to search the NCBI conserved domain database" \
                  "for significant domain hits in all new proteins of a " \
                  "Phamerator database."

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
    return parser


def make_tempdir():
    """
    Uses pdm_utils.functions.basic.expand_path to expand TMP_DIR; then
    checks whether TMP_DIR exists - if it doesn't, uses os.makedirs to
    make it recursively.
    :return:
    """
    try:
        path = expand_path(TMP_DIR)
        # If the path doesn't exist yet
        if not os.path.exists(path):
            # Make it recursively
            os.makedirs(path)
    except OSError as err:
        print(f"Error {err.args[0]}: {err.args[1]}")


def search_and_process(datum, database, evalue,
                       domain_queue, gene_domain_queue,
                       gene_queue, tmp="/tmp/cdd"):
    """
    Uses blast to search indicated gene against the indicated database
    :param database: path to target database
    :param datum: dictionary containing GeneID and Translation for the
    gene being processed
    :param evalue: evalue cutoff for rpsblast hits
    :param domain_queue: Queue to store 'insert ignore into domain'
    statements
    :param gene_domain_queue: Queue to store 'insert into gene domain'
    statements
    :param gene_queue: Queue to store 'update gene set DomainStatus'
    statements
    :param tmp: path to directory where I/O files will be stored
    :return:
    """
    # Setup I/O variables
    i = "{}/{}.txt".format(tmp, datum["GeneID"])
    o = "{}/{}.xml".format(tmp, datum["GeneID"])

    # Write the input file
    f = open(i, "w")
    f.write(">{}\n{}".format(datum["GeneID"], datum["Translation"]))
    f.close()

    # Setup and run the rpsblast command
    rps_command = NcbirpsblastCommandline(cmd=RPSBLAST_BINARY, db=database,
                                          query=i, out=o, outfmt=5,
                                          evalue=evalue)
    rps_command()

    """
    command = "{} -db {} -query {} -out {} -outfmt 5 -evalue {}".format(
        RPSBLAST_BINARY, database, i, o, evalue)
    os.system(command)
    """
    with open(o, "r") as oh:
        for record in NCBIXML.parse(oh):
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
                            domain_queue.put(INSERT_INTO_DOMAIN.format(
                                    align.hit_id, domain_id, name, description))

                            # Try to put this hit into gene_domain table
                            gene_domain_queue.put(INSERT_INTO_GENE_DOMAIN.format(
                                    datum["GeneID"], align.hit_id,
                                    float(hsp.expect), int(hsp.query_start),
                                    int(hsp.query_end)))

    # Update this gene's DomainStatus to 1
    gene_queue.put(UPDATE_GENE.format(datum["GeneID"]))


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
    db = args.db
    cdd_dir = expand_path(args.dir)
    cdd_db = os.path.join(cdd_dir, cdd_dir.split("/")[-1])
    threads = args.threads
    evalue = args.evalue

    # Use MySQLConnectionHandler to query for translations that need to be
    # put through this pipeline
    mysql_handler = MySQLConnectionHandler(database=db)
    mysql_handler.open_connection()     # Gets user/pass automatically

    if mysql_handler.connection_status() is not True:
        sys.exit(1)                     # If user/pass wrong too many times

    results = mysql_handler.execute_query(GET_GENES_FOR_CDD)
    mysql_handler.close_connection()

    # Print number of genes to process
    print(f"{len(results)} genes to search for conserved domains...")

    # Only run the pipeline if there are genes returned that need it
    if len(results) > 0:
        # Create temp_dir
        make_tempdir()

        # Create threadsafe queues for workers to put resulting SQL commands
        # into
        domain_queue = mp.Manager().Queue()
        gene_domain_queue = mp.Manager().Queue()
        gene_queue = mp.Manager().Queue()

        # create list of job tuples - each tuple contains all the arguments
        # for the search_and_process function
        jobs = [(gene, cdd_db, evalue, domain_queue, gene_domain_queue,
                 gene_queue) for gene in results]

        # Create worker pool and use starmap to run concurrent jobs with
        # multiple args
        with mp.Pool(threads) as pool:
            pool.starmap(search_and_process, jobs)

        mysql_handler.open_connection()
        # Insert into domain commands
        while not domain_queue.empty():
            mysql_handler.execute_transaction([domain_queue.get()])

        # Insert into gene_domain commands
        while not gene_domain_queue.empty():
            mysql_handler.execute_transaction([gene_domain_queue.get()])

        # Update gene commands
        while not gene_queue.empty():
            mysql_handler.execute_transaction([gene_queue.get()])

    return
