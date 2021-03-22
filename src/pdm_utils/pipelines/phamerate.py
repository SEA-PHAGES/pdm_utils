"""Program to group related gene products into phamilies using either MMseqs2
for both similarity search and clustering, or blastp for similarity search
and mcl for clustering."""

import argparse
from datetime import datetime
import os
import shutil

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions.configfile import *
from pdm_utils.functions.phameration import *
from pdm_utils.functions.parallelize import *

MMSEQS_DESCRIPTION = """
Assort protein sequences into phamilies using MMseqs2.

[1] Steinegger M and Soeding J. MMseqs2 enables sensitive protein 
sequence searching for the analysis of massive data sets. Nature 
Biotechnology, 2017. doi: 10.1038/nbt.3988).
"""

BLAST_DESCRIPTION = """
Assort protein sequences into phamilies using blastp for homology 
search and mcl to extract groups of homologs from the blastp output.

[1] Altschul et al. Basic Local Alignment Search Tool. Journal of 
Molecular Biology, 1990. doi: 10.1016/S0022-2836(05)80360-2
[2] Stijn van Dongen and Cei Abreu-Goodger. Using MCL to Extract 
Clusters From Networks. Methods in Molecular Biology, 2012. doi: 
10.1007/978-1-61779-361-5_15.
"""


def setup_argparser():
    """
    Builds argparse.ArgumentParser for this script
    :return:
    """
    # Initialize parser and add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()

    # Create sub-parser for MMseqs2 invocation
    mmseqs_parser = subparsers.add_parser("mmseqs",
                                          help="pham assembly uses MMseqs2")
    mmseqs_parser.description = MMSEQS_DESCRIPTION
    mmseqs_parser.add_argument("db", type=str,
                               help="name of database to phamerate")
    mmseqs_parser.add_argument("--identity", type=float, default=0.3,
                               help="pct identity to enter pre-pham [0, 1]")
    mmseqs_parser.add_argument("--coverage", type=float, default=0.85,
                               help="pct coverage to enter pre-pham [0, 1]")
    mmseqs_parser.add_argument("--e-value", type=float, default=0.001,
                               help="E-value to enter pre-pham HMM")
    mmseqs_parser.add_argument("--sens", type=float, default=4,
                               help="sensitivity in range [1, 8.5]")
    mmseqs_parser.add_argument("--steps", type=int, default=1,
                               help="number of steps to build pre-pham HMMs")
    mmseqs_parser.add_argument("--aln-mode", type=int, default=3,
                               help="alignment mode in range [0, 4]")
    mmseqs_parser.add_argument("--cov-mode", type=int, default=0,
                               help="coverage mode in range [0, 4]")
    mmseqs_parser.add_argument("--clu-mode", type=int, default=1,
                               help="cluster mode in range [0, 3]")
    mmseqs_parser.add_argument("--skip-hmm", action="store_true",
                               help="skip the HMM clustering step")
    mmseqs_parser.add_argument("--hmmident", type=float, default=0.25,
                               help="pct identity to cluster HMMs [0, 1]")
    mmseqs_parser.add_argument("--hmmcover", type=float, default=0.50,
                               help="pct coverage to cluster HMMs [0, 1]")
    mmseqs_parser.add_argument("--hmm-eval", type=float, default=0.001,
                               help="E-value to cluster HMMs")
    mmseqs_parser.add_argument("--threads", type=int, default=mp.cpu_count(),
                               help="number of threads to use")
    mmseqs_parser.add_argument("--tmp-dir", type=str, default="/tmp/phamerate",
                               help="temporary directory for file I/O")
    mmseqs_parser.add_argument("-c", "--config_file", type=pathlib.Path, default=None,
                               help="path to file containing login details")
    mmseqs_parser.formatter_class = argparse.RawTextHelpFormatter

    # Create sub-parser for blast-mcl invocation
    blast_parser = subparsers.add_parser("blast-mcl",
                                         help="pham assembly uses blastp+mcl")
    blast_parser.description = BLAST_DESCRIPTION
    blast_parser.add_argument("db", type=str,
                              help="name of database to phamerate")
    blast_parser.add_argument("--e-value", type=float, default=0.001,
                              help="blastp e-value to store result")
    blast_parser.add_argument("--query-cov", type=float, default=0.50,
                              help="blastp query coverage to keep HSPs [0, 1]")
    blast_parser.add_argument("--inflate", type=float, default=5.0,
                              help="MCL inflation parameter")
    blast_parser.add_argument("--threads", type=int, default=mp.cpu_count(),
                              help="blastp instances to run in parallel")
    blast_parser.add_argument("--tmp-dir", type=str, default="/tmp/phamerate",
                              help="temporary directory for file I/O")
    blast_parser.add_argument("-c", "--config_file", type=pathlib.Path, default=None,
                              help="path to file containing login details")
    blast_parser.formatter_class = argparse.RawTextHelpFormatter
    return parser


def refresh_tempdir(tmpdir):
    """
    Recursively deletes tmpdir if it exists, otherwise makes it
    :param tmpdir: directory to refresh
    :return:
    """
    if os.path.exists(tmpdir):
        try:
            shutil.rmtree(tmpdir)
        except OSError:
            print(f"Failed to delete existing temp directory '{tmpdir}'")
            return
    try:
        os.makedirs(tmpdir)
    except OSError:
        print(f"Failed to create new temp directory '{tmpdir}")
        return


def main(argument_list):
    # Set up the argument parser
    parser = setup_argparser()

    # Parse arguments
    args = vars(parser.parse_args(argument_list[2:]))

    # Temporary directory gets its own variable because we'll use it a lot
    tmp = args["tmp_dir"]

    # Create config object with data obtained from file and/or defaults.
    config = build_complete_config(args["config_file"])
    mysql_creds = config["mysql"]

    # Make a note of which workflow we're using based on len(args)
    if len(args) > 10:
        program = "mmseqs"
    else:
        program = "blast-mcl"

    # Record start time
    start_time = datetime.now()

    # Initialize SQLAlchemy engine with database provided at CLI
    alchemist = AlchemyHandler(database=args["db"],
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(login_attempts=5, pipeline=True)
    engine = alchemist.engine

    # Refresh temp_dir
    refresh_tempdir(tmp)

    # Get old pham data and un-phamerated genes
    old_phams = get_pham_geneids(engine)
    old_colors = get_pham_colors(engine)
    new_genes = get_new_geneids(engine)

    # Get GeneIDs & translations, and translation groups
    # gene_x: translation_x
    genes_and_translations = get_geneids_and_translations(engine)
    # translation_x: [gene_x, ..., gene_z]
    translation_groups = get_translation_groups(engine)

    # Print initial state
    initial_summary = f"""
Initial database summary:
=============================
 {len(old_phams)} total phams
 {sum([len(x) == 1 for x in old_phams.values()])} orphams
 {sum([len(x) for x in old_phams.values()])} genes in phams
 {len(new_genes)} genes not in phams
 {len(genes_and_translations)} total genes
 {len(translation_groups)} non-redundant genes
=============================
"""
    print(initial_summary)

    # Write input fasta file
    print("Writing non-redundant sequences to input fasta...")
    infile = f"{tmp}/input.fasta"
    write_fasta(translation_groups, infile)

    # Here is where the workflow selection comes into play
    if program == "mmseqs":
        seq_db = f"{tmp}/sequenceDB"            # MMseqs2 sequence database
        clu_db = f"{tmp}/clusterDB"             # MMseqs2 cluster database
        psf_db = f"{tmp}/seqfileDB"             # pre-pham seqfile database
        p_out = f"{tmp}/pre_out.fasta"          # pre-pham output (FASTA)

        print("Creating MMseqs2 sequence database...")
        mmseqs_createdb(infile, seq_db)

        print("Clustering sequence database...")
        mmseqs_cluster(seq_db, clu_db, args)

        print("Storing sequence-based phamilies...")
        mmseqs_createseqfiledb(seq_db, clu_db, psf_db)
        mmseqs_result2flat(seq_db, seq_db, psf_db, p_out)
        pre_phams = parse_mmseqs_output(p_out)      # Parse pre-pham output

        # Proceed with profile clustering, if allowed
        if not args["skip_hmm"]:
            con_lookup = dict()
            for name, geneids in pre_phams.items():
                for geneid in geneids:
                    con_lookup[geneid] = name

            pro_db = f"{tmp}/profileDB"         # MMseqs2 profile database
            con_db = f"{tmp}/consensusDB"       # Consensus sequence database
            aln_db = f"{tmp}/alignDB"           # Alignment database
            res_db = f"{tmp}/resultDB"          # Cluster database
            hsf_db = f"{tmp}/hmmSeqfileDB"      # hmm-pham seqfile database
            h_out = f"{tmp}/hmm_out.fasta"      # hmm-pham output (FASTA)

            print("Creating HMM profiles from sequence-based phamilies...")
            mmseqs_result2profile(seq_db, clu_db, pro_db)

            print("Extracting consensus sequences from HMM profiles...")
            mmseqs_profile2consensus(pro_db, con_db)

            print("Searching for profile-profile hits...")
            mmseqs_search(pro_db, con_db, aln_db, args)

            print("Clustering based on profile-profile alignments...")
            mmseqs_clust(con_db, aln_db, res_db)

            print("Storing profile-based phamilies...")
            mmseqs_createseqfiledb(seq_db, res_db, hsf_db)
            mmseqs_result2flat(pro_db, con_db, hsf_db, h_out)
            hmm_phams = parse_mmseqs_output(h_out)

            print("Merging sequence and profile-based phamilies...")
            new_phams = merge_pre_and_hmm_phams(
                hmm_phams, pre_phams, con_lookup)
        else:
            new_phams = pre_phams
    else:
        blast_db = "blastdb"
        blast_path = f"{tmp}/{blast_db}"

        print("Creating blast protein database...")
        create_blastdb(infile, blast_db, blast_path)

        print("Splitting non-redundant sequences into multiple blastp query "
              "files...")
        chunks = chunk_translations(translation_groups)

        jobs = []
        for key, chunk in chunks.items():
            jobs.append((key, chunk, tmp, blast_path, 
                         args["e_value"], args["query_cov"]))

        print("Running blastp...")
        parallelize(jobs, args["threads"], blastp)

        print("Converting blastp output into adjacency matrix for mcl...")
        results = [x for x in os.listdir(tmp) if x.endswith(".tsv")]
        adjacency = f"{tmp}/blast_adjacency.abc"
        with open(adjacency, "w") as fh:
            for result in results:
                f = open(f"{tmp}/{result}", "r")
                for line in f:
                    fh.write(line)
                f.close()

        print("Running mcl on adjacency matrix...")
        outfile = markov_cluster(adjacency, args["inflate"], tmp)

        print("Storing blast-mcl phamilies...")
        new_phams = parse_mcl_output(outfile)

        # Some proteins don't have even self-hits in blastp - take a
        # census of who is missing, and add them as "orphams"
        mcl_genes = set()
        for name, pham in new_phams.items():
            for gene in pham:
                mcl_genes.add(genes_and_translations[gene])

        all_trans = set(translation_groups.keys())

        # Some genes don't have blast hits, even to themselves. These are
        # not in the blast output and need to be re-inserted as orphams.
        missing = all_trans - mcl_genes
        for translation in missing:
            new_phams[len(new_phams) + 1] = \
                [translation_groups[translation][0]]

    # Reintroduce duplicates
    print("Propagating phamily assignments to duplicate genes...")
    new_phams = reintroduce_duplicates(new_phams, translation_groups,
                                       genes_and_translations)

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

    # Print final state
    final_summary = f"""
Final database summary:
=============================
 {len(new_phams)} total phams
 {sum([len(x) == 1 for x in new_phams.values()])} orphams
 {sum([len(x) for x in new_phams.values()])} genes in phams
 {len(genes_and_translations) - sum([len(x) for x in new_phams.values()])} genes not in phams
 {len(genes_and_translations)} total genes
 {len(translation_groups)} non-redundant genes
=============================
"""
    print(final_summary)

    # Record stop time
    stop_time = datetime.now()
    elapsed_time = str(stop_time - start_time)

    # Report phameration elapsed time
    print(f"Elapsed time: {elapsed_time}")


if __name__ == "__main__":
    main(sys.argv[1:])
