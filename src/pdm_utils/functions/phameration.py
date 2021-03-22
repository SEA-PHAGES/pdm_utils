"""Functions that are used in the phameration pipeline"""

import shlex
from subprocess import Popen, PIPE
import random
import colorsys

from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic


# DATABASE FUNCTIONS
def get_pham_geneids(engine):
    """
    Queries the database for those genes that are already phamerated.
    :param engine: the Engine allowing access to the database
    :return: pham_geneids
    """
    pham_geneids = dict()

    geneid_query = "SELECT GeneID, PhamID FROM gene WHERE PhamID IS NOT NULL"
    geneid_results = mysqldb_basic.query_dict_list(engine, geneid_query)

    for dictionary in geneid_results:
        pham_id = dictionary["PhamID"]
        geneid = dictionary["GeneID"]

        if pham_id in pham_geneids.keys():
            pham_geneids[pham_id] = pham_geneids[pham_id] | {geneid}
        else:
            pham_geneids[pham_id] = {geneid}

    return pham_geneids


def get_pham_colors(engine):
    """
    Queries the database for the colors of existing phams
    :param engine: the Engine allowing access to the database
    :return: pham_colors
    """
    pham_colors = dict()

    color_query = "SELECT PhamID, Color FROM pham"
    color_results = mysqldb_basic.query_dict_list(engine, color_query)

    for dictionary in color_results:
        pham_id = dictionary["PhamID"]
        color = dictionary["Color"]

        pham_colors[pham_id] = color.upper()

    return pham_colors


def get_new_geneids(engine):
    """
    Queries the database for those genes that are not yet phamerated.
    :param engine: the Engine allowing access to the database
    :return: new_geneids
    """
    new_geneids = list()

    gene_query = "SELECT GeneID FROM gene WHERE PhamID IS NULL"
    gene_results = mysqldb_basic.query_dict_list(engine, gene_query)

    # At scale, much cheaper to convert list to set than to build the
    # set one gene at a time
    for dictionary in gene_results:
        geneid = dictionary["GeneID"]
        new_geneids.append(geneid)

    return set(new_geneids)


def get_geneids_and_translations(engine):
    """
    Constructs a dictionary mapping all geneids to their translations.
    :param engine: the Engine allowing access to the database
    :return: gs_to_ts
    """
    gs_to_ts = dict()

    query = ("SELECT GeneID, CONVERT(Translation USING utf8) as Translation "
             "FROM gene")
    results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        translation = dictionary["Translation"]
        gs_to_ts[geneid] = translation

    return gs_to_ts


def get_translation_groups(engine):
    """
    Constructs a dictionary mapping all unique translations to their
    groups of geneids that share them
    :param engine: the Engine allowing access to the database
    :return: ts_to_gs
    """
    ts_to_gs = dict()

    query = ("SELECT GeneID, CONVERT(Translation USING utf8) as Translation "
             "FROM gene")
    results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        trans = dictionary["Translation"]
        geneids = ts_to_gs.get(trans, [])
        geneids.append(geneid)
        ts_to_gs[trans] = geneids

    return ts_to_gs


def update_pham_table(colors, engine):
    """
    Populates the pham table with the new PhamIDs and their colors.
    :param colors: new pham color data
    :type colors: dict
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    # First command needs to clear the pham table
    commands = ["DELETE FROM pham"]

    # Then we need to issue insert commands for each pham
    for key in colors.keys():
        commands.append(f"INSERT INTO pham (PhamID, Color) VALUES ({key}, "
                        f"'{colors[key]}')")

    mysqldb.execute_transaction(engine, commands)


def update_gene_table(phams, engine):
    """
    Updates the gene table with new pham data
    :param phams: new pham gene data
    :type phams: dict
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    commands = []

    # We need to issue an update command for each gene in each pham
    for key in phams.keys():
        for gene in phams[key]:
            commands.append(f"UPDATE gene SET PhamID = {key} WHERE GeneID = '{gene}'")

    mysqldb.execute_transaction(engine, commands)


def fix_white_phams(engine):
    """
    Find any phams with 2+ members which are colored as though they are
    orphams (#FFFFFF in pham.Color).
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    query = "SELECT c.PhamID FROM (SELECT g.PhamID, COUNT(GeneID) AS count, "\
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID " \
            "= p.PhamID GROUP BY PhamID) AS c WHERE Color = '#FFFFFF' "\
            "AND count > 1"

    results = mysqldb_basic.query_dict_list(engine, query)
    print(f"Found {len(results)} white phams...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        h = s = v = 0
        while h <= 0:
            h = random.random()
        while s < 0.5:
            s = random.random()
        while v < 0.8:
            v = random.random()
        rgb = colorsys.hsv_to_rgb(h, s, v)
        rgb = (rgb[0] * 255, rgb[1] * 255, rgb[2] * 255)
        hexrgb = "#{:02x}{:02x}{:02x}".format(int(rgb[0]), int(rgb[1]),
                                              int(rgb[2]))
        new_color = hexrgb.upper()
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE "
                        f"PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)


def fix_colored_orphams(engine):
    """
    Find any single-member phams which are colored as though they are
    multi-member phams (not #FFFFFF in pham.Color).
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    query = "SELECT * FROM (SELECT g.PhamID, COUNT(GeneID) AS count, " \
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID " \
            "=p.PhamID GROUP BY PhamID) AS c WHERE Color != '#FFFFFF' " \
            "AND count = 1"

    results = mysqldb_basic.query_dict_list(engine, query)
    print(f"Found {len(results)} non-white orphams...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        count = dictionary["count"]
        color = dictionary["Color"]
        new_color = "#FFFFFF"
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE "
                        f"PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)


# FILE I/O FUNCTIONS
def write_fasta(translation_groups, outfile):
    """
    Writes a FASTA file of the non-redundant protein sequences to be
    assorted into phamilies.
    :param translation_groups: groups of genes that share a translation
    :type translation_groups: dict
    :param outfile: FASTA filename
    :type outfile: str
    :return:
    """
    fasta = open(f"{outfile}", "w")
    for translation in translation_groups.keys():
        fasta.write(f">{translation_groups[translation][0]}\n{translation}\n")
    fasta.close()


def parse_mmseqs_output(outfile):
    """
    Parses the indicated MMseqs2 FASTA-like file into a dictionary of
    integer-named phams.
    :param outfile: FASTA-like parseable output
    :type outfile: str
    :return: phams
    :rtype: dict
    """
    phams = {}
    pham_geneids = list()
    pham_name = 0

    with open(f"{outfile}", "r") as fh:
        prior = fh.readline()
        current = fh.readline()

        # While loop to interate until EOF
        while current:
            # If current line is a header line
            if current.startswith(">"):
                # If the prior line was also a header line
                if prior.startswith(">"):
                    # We've reached a new pham block
                    try:
                        # We need to remove prior line's geneid before dumping
                        pham_geneids.pop(-1)
                    except IndexError:
                        # This should happen for the first pham dump
                        pass
                    phams[pham_name] = pham_geneids
                    pham_name += 1
                    pham_geneids = [current.lstrip(">").rstrip()]
                # Otherwise, we're still in a pham block
                else:
                    pham_geneids.append(current.lstrip(">").rstrip())
            # If the current line is a translation line
            else:
                # We don't care about it - do nothing and move on
                pass
            prior, current = current, fh.readline()
        # Need to dump the last pham into the dictionary
        phams[pham_name] = pham_geneids
    # Pham 0 is a placeholder
    phams.pop(0)

    return phams


def parse_mcl_output(outfile):
    """
    Parse the mci output into phams
    :param outfile: mci output file
    :type outfile: str
    :return: phams
    :rtype: dict
    """
    phams = dict()

    counter = 0
    with open(outfile, "r") as fh:
        for line in fh:
            counter += 1
            pham = line.rstrip().split()
            for i in range(len(pham)):
                gene = pham[i]
                if "|" in gene:
                    gene = gene.split("|")
                    gene = "_".join(gene[1:])
                pham[i] = gene
            phams[counter] = pham

    return phams


# PHAM MANIPULATION FUNCTIONS
def merge_pre_and_hmm_phams(hmm_phams, pre_phams, consensus_lookup):
    """
    Merges the pre-pham sequences (which contain all nr sequences) with
    the hmm phams (which contain only hmm consensus sequences) into the
    full hmm-based clustering output. Uses consensus_lookup dictionary
    to find the pre-pham that each consensus belongs to, and then adds
    each pre-pham geneid to a full pham based on the hmm phams.
    :param hmm_phams: clustered consensus sequences
    :type hmm_phams: dict
    :param pre_phams: clustered sequences (used to generate hmms)
    :type pre_phams: dict
    :param consensus_lookup: reverse-mapped pre_phams
    :type consensus_lookup: dict
    :return: phams
    :rtype: dict
    """
    phams = dict()

    for name, members in hmm_phams.items():
        full_pham = list()
        for member in members:
            pre_pham_key = consensus_lookup[member]
            for geneid in pre_phams[pre_pham_key]:
                full_pham.append(geneid)
        phams[name] = full_pham

    return phams


def reintroduce_duplicates(new_phams, trans_groups, genes_and_trans):
    """
    Reintroduces into each pham ALL GeneIDs that map onto the set of
    translations in the pham.
    :param new_phams: the pham dictionary for which duplicates are
    to be reintroduced
    :param trans_groups: the dictionary that maps translations to the
    GeneIDs that share them
    :param genes_and_trans: the dictionary that maps GeneIDs to their
    translations
    :return:
    """
    for key in new_phams.keys():
        geneids = new_phams[key]
        dup_geneids = list()
        for geneid in geneids:
            translation = genes_and_trans[geneid]
            dup_group = trans_groups[translation]
            for gene in dup_group:
                dup_geneids.append(gene)
        new_phams[key] = set(dup_geneids)
    return new_phams


def preserve_phams(old_phams, new_phams, old_colors, new_genes):
    """
    Attempts to keep pham numbers consistent from one round of pham
    building to the next
    :param old_phams: the dictionary that maps old phams to their genes
    :param new_phams: the dictionary that maps new phams to their genes
    :param old_colors: the dictionary that maps old phams to colors
    :param new_genes: the set of previously unphamerated genes
    :return:
    """
    final_phams = dict()
    final_colors = dict()

    outcount = 0
    total = len(old_phams)
    new_phams_copy = new_phams.copy()

    # Iterate through old and new phams
    for old_key in old_phams.keys():
        outcount += 1
        # print("Pham Name Conservation: {} / {}".format(outcount, total))

        old_pham = old_phams[old_key]

        if old_key in final_phams.keys():
            continue

        for new_key in new_phams_copy.keys():
            new_pham = new_phams_copy[new_key]

            if old_pham & new_pham == set():
                continue

            # Case 1 + 5 (Identity and Subtraction)
            if old_pham == new_pham:
                final_phams[old_key] = new_pham
                final_colors[old_key] = old_colors[old_key]
                new_phams_copy.pop(new_key)
                break

            # Case 2 and 4 (Addition and Join) - PHAM GREW
            elif new_pham - old_pham != set():

                # Case 2 and 4 (Addition and Join)
                if new_pham & new_genes != set():

                    # Case 4 - Join with new gene
                    if (new_pham - (new_pham & new_genes)) - old_pham != set():
                        break

                    # Case 2 - Addition with new gene
                    final_phams[old_key] = new_pham
                    final_colors[old_key] = old_colors[old_key]
                    new_phams_copy.pop(new_key)
                    break

                # Case 4 - Join without new gene
                else:
                    break

            # Case 3 - split - PHAM SHRANK, BUT NOT BY REMOVAL
            elif old_pham - new_pham != set():
                break

    final_phams[0] = "placeholder"
    highest_pham = max(map(int, final_phams.keys())) + 1
    final_phams.pop(0)

    # Reassign data for split or joined phams
    for key in new_phams_copy.keys():
        new_key = highest_pham
        highest_pham += 1

        final_phams[new_key] = new_phams_copy[key]

        # Multi-member pham gets non-white color
        if len(new_phams_copy[key]) > 1:
            h = s = v = 0
            while h <= 0:
                h = random.random()
            while s < 0.5:
                s = random.random()
            while v < 0.8:
                v = random.random()
            rgb = colorsys.hsv_to_rgb(h, s, v)
            rgb = (rgb[0] * 255, rgb[1] * 255, rgb[2] * 255)
            hexrgb = "#{:02x}{:02x}{:02x}".format(int(rgb[0]), int(rgb[1]),
                                                  int(rgb[2]))
            final_colors[new_key] = hexrgb.upper()
        # Orpham gets white
        else:
            final_colors[new_key] = '#FFFFFF'

    return final_phams, final_colors


# MMSEQS2 CLUSTERING FUNCTIONS
def mmseqs_createdb(fasta, sequence_db):
    """
    Runs 'mmseqs createdb' to convert a FASTA file into an MMseqs2
    sequence database.
    :param fasta: path to the FASTA file to convert
    :type fasta: str
    :param sequence_db: MMseqs2 sequence database
    :type sequence_db: str
    """
    command = f"mmseqs createdb {fasta} {sequence_db} -v 3"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_cluster(sequence_db, cluster_db, args):
    """
    Runs 'mmseqs cluster' to cluster an MMseqs2 sequence database.
    :param sequence_db: MMseqs2 sequence database
    :type sequence_db: str
    :param cluster_db: MMseqs2 clustered database
    :type cluster_db: str
    :param args: parsed command line arguments
    :type args: dict
    """
    command = f"mmseqs cluster {sequence_db} {cluster_db} {args['tmp_dir']}" \
              f" -v 3 --min-seq-id {args['identity']} -c {args['coverage']} " \
              f"-e {args['e_value']} -s {args['sens']} --max-seqs 1000 " \
              f"--cluster-steps {args['steps']} --threads {args['threads']} " \
              f"--alignment-mode {args['aln_mode']} --cov-mode " \
              f"{args['cov_mode']} --cluster-mode {args['clu_mode']} --cluster-reassign"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_result2profile(sequence_db, cluster_db, profile_db):
    """
    Runs 'mmseqs result2profile' to convert clusters from one MMseqs2
    clustered database into a profile database.
    :param sequence_db: MMseqs2 sequence database
    :type sequence_db: str
    :param cluster_db: MMseqs2 clustered database
    :type cluster_db: str
    :param profile_db: MMseqs2 profile database
    :type profile_db: str
    """
    command = f"mmseqs result2profile {sequence_db} {sequence_db} " \
              f"{cluster_db} {profile_db} -v 3"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_profile2consensus(profile_db, consensus_db):
    """
    Runs 'mmseqs profile2consensus' to extract consensus sequences from
    an MMseqs2 profile database, and creates an MMseqs2 sequence
    database from the consensuses.
    :param profile_db: MMseqs2 profile database
    :type profile_db: str
    :param consensus_db: MMseqs2 sequence database
    :type consensus_db: str
    """
    command = f"mmseqs profile2consensus {profile_db} {consensus_db} -v 3"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_search(profile_db, consensus_db, align_db, args):
    """
    Runs 'mmseqs search' to search profiles against their consensus
    sequences and save the alignment results to an MMseqs2 alignment
    database. The profile_db and consensus_db MUST be the same size.
    :param profile_db: MMseqs2 profile database
    :type profile_db: str
    :param consensus_db: MMseqs2 sequence database
    :type consensus_db: str
    :param align_db: MMseqs2 alignment database
    :type align_db: str
    :param args: parsed command line arguments
    :type args: dict
    """
    command = f"mmseqs search {profile_db} {consensus_db} {align_db} " \
              f"{args['tmp_dir']} --min-seq-id {args['hmmident']} -c " \
              f"{args['hmmcover']} --e-profile {args['hmm_eval']} -v 3 " \
              f"--add-self-matches --max-seqs 1000 -e {args['hmm_eval']}"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_clust(consensus_db, align_db, cluster_db):
    """
    Runs 'mmseqs clust' to cluster an MMseqs2 consensus database
    using an MMseqs2 alignment database, with results being saved to
    an MMseqs2 cluster database.
    :param consensus_db: MMseqs sequence database
    :type consensus_db: str
    :param align_db: MMseqs2 alignment database
    :type align_db: str
    :param cluster_db: MMseqs2 cluster database
    :type cluster_db: str
    """
    command = f"mmseqs clust {consensus_db} {align_db} {cluster_db}"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_createseqfiledb(sequence_db, cluster_db, seqfile_db):
    """
    Runs 'mmseqs createseqfiledb' to create the intermediate to the
    FASTA-like parseable output.
    :param sequence_db: MMseqs2 sequence database
    :type sequence_db: str
    :param cluster_db: MMseqs2 clustered database
    :type cluster_db: str
    :param seqfile_db: MMseqs2 seqfile database
    :type seqfile_db: str
    """
    command = f"mmseqs createseqfiledb {sequence_db} {cluster_db} " \
              f"{seqfile_db} -v 3"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def mmseqs_result2flat(query_db, target_db, seqfile_db, outfile):
    """
    Runs 'mmseqs result2flat' to create FASTA-like parseable output.
    :param query_db: MMseqs2 sequence or profile database
    :type query_db: str
    :param target_db: MMseqs2 sequence database
    :type target_db: str
    :param seqfile_db: MMseqs2 seqfile database
    :type seqfile_db: str
    :param outfile: FASTA-like parseable output
    :type outfile: str
    """
    command = f"mmseqs result2flat {query_db} {target_db} {seqfile_db} " \
              f"{outfile} -v 3"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


# BLAST-MCL CLUSTERING FUNCTIONS
def create_blastdb(fasta, db_name, db_path):
    """
    Runs 'makeblastdb' to create a BLAST-searchable database.
    :param fasta: FASTA-formatted input file
    :type fasta: str
    :param db_name: BLAST sequence database
    :type db_name: str
    :param db_path: BLAST sequence database path
    :type db_path: str
    """
    command = f"makeblastdb -in {fasta} -dbtype prot -title {db_name} " \
              f"-parse_seqids -out {db_path}"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()


def chunk_translations(translation_groups, chunksize=500):
    """
    Break translation_groups into a dictionary of chunksize-tuples of
    2-tuples where each 2-tuple is a translation and its corresponding
    geneid.
    :param translation_groups: translations and their geneids
    :type translation_groups: dict
    :param chunksize: how many translations will be in a chunk?
    :type chunksize: int
    :return: chunks
    :rtype: dict
    """
    chunks = dict()
    keys = list(translation_groups.keys())
    num_chunks = len(keys) // chunksize

    index = 0
    for i in range(num_chunks):
        temp_chunk = tuple((x, translation_groups[x][0]) for x in keys[index:index+chunksize])
        chunks[i] = temp_chunk
        index += chunksize

    # Add any leftovers to one last (smaller) chunk
    temp_chunk = tuple((x, translation_groups[x][0]) for x in keys[index:])
    chunks[num_chunks + 1] = temp_chunk

    return chunks


def blastp(index, chunk, tmp, db_path, evalue, query_cov):
    """
    Runs 'blastp' using the given chunk as the input gene set. The
    blast output is an adjacency matrix for this chunk.
    :param index: chunk index being run
    :type index: int
    :param chunk: the translations to run right now
    :type chunk: tuple of 2-tuples
    :param tmp: path where I/O can go on
    :type tmp: str
    :param db_path: path to the target blast database
    :type db_path: str
    :param evalue: e-value cutoff to report hits
    :type evalue: float
    """
    in_name = f"{tmp}/input{index}.fasta"
    out_name = f"{tmp}/output{index}.tsv"

    with open(in_name, "w") as fh:
        for t in chunk:
            fh.write(f">{t[1]}\n{t[0]}\n")

    command = f"blastp -query {in_name} -db {db_path} -out {out_name} " \
              f"-outfmt '6 qseqid sseqid evalue' -max_target_seqs " \
              f"10000 -num_threads 1 -use_sw_tback -evalue {evalue} " \
              f"-qcov_hsp_perc {int(100*query_cov)} -max_hsps 1"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        stdout = process.stdout.read().decode("utf-8")
        stderr = process.stderr.read().decode("utf-8")
        if stdout != "":
            print(stdout)
        if stderr != "":
            print(stderr)

    return [out_name]


def markov_cluster(adj_mat_file, inflation, tmp_dir):
    """
    Run 'mcl' on an adjacency matrix to cluster the blastp results.
    :param adj_mat_file: 3-column file with blastp resultant
    queries, subjects, and evalues
    :type adj_mat_file: str
    :param inflation: mcl inflation parameter
    :type inflation: float
    :param tmp_dir: file I/O directory
    :type tmp_dir: str
    :return: outfile
    :rtype: str
    """
    outfile = f"{tmp_dir}/mcl_clusters.txt"
    command = f"mcl {adj_mat_file} -I {inflation} --abc -o {outfile} " \
              f"-abc-tf 'ceil(200)' --abc-neg-log10"
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        # print(process.stdout.read().decode("utf-8"))
        # print(process.stderr.read().decode("utf-8"))
        process.wait()

    return outfile
