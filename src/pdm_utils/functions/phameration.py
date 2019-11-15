"""Functions that are used in the phameration pipeline"""

import shlex
from subprocess import Popen, PIPE
import csv
import random
import colorsys

from pdm_utils.constants.constants import BLASTCLUST_PATH


def get_program_params(program, args):
    if program == "blast":
        params = {"-S": args.identity,
                  "-L": float(args.coverage)/100,
                  "-a": args.threads}
    elif program == "mmseqs":
        params = {"--threads": args.threads,
                  "-v": args.verbose,
                  "--cluster-steps": args.steps,
                  "--max-seqs": args.max_seqs,
                  "--min-seq-id": float(args.identity)/100,
                  "-c": float(args.coverage)/100,
                  "--alignment-mode": args.aln_mode,
                  "--cov-mode": args.cov_mode,
                  "--cluster-mode": args.clu_mode}
    else:
        print(f"Unknown program {program}")
        return

    return params


def read_existing_phams(mysql_handler):
    """
    Uses the provided mysql_handler to query the database for existing
    pham names and colors, and returns two dictionaries: one mapping
    phams to geneids, and one mapping phams to colors
    :param mysql_handler:
    :return: old_phams, old_colors
    """
    old_phams = dict()
    old_colors = dict()

    query = "SELECT * FROM (SELECT a.Name, a.GeneID, b.Color FROM pham AS a " \
            "INNER JOIN pham_color AS b on a.Name = b.Name) AS c ORDER BY " \
            "c.Name ASC"
    results = mysql_handler.execute_query(query)

    print("{} old genes".format(len(results)))

    for dictionary in results:
        name = dictionary["Name"]
        geneid = dictionary["GeneID"]
        color = dictionary["Color"]
        if name in old_phams.keys():
            old_phams[name] = old_phams[name] | {geneid}
        else:
            old_phams[name] = {geneid}
            old_colors[name] = color

    return old_phams, old_colors


def read_unphamerated_genes(mysql_handler):
    """
    Uses the provided mysql_handler to query the database for genes
    that are not represented in the pham table (genes whose phages
    were updated in the current round of database updates). Returns
    the list of GeneIDs not yet in a pham
    :return:
    """
    new_genes = set()

    query = "SELECT GeneID FROM gene WHERE GeneID NOT IN (SELECT " \
            "GeneID from pham)"
    results = mysql_handler.execute_query(query)

    print("{} new genes".format(len(results)))

    for dictionary in results:
        geneid = dictionary["GeneID"]
        new_genes = new_genes | {geneid}

    return new_genes


def get_translations(mysql_handler):
    """
    Uses the provided mysql_handler to query the database for all
    GeneIDs and translations.  Returns two dictionaries: one mapping
    GeneIDs to translations, and the other mapping translations to
    the list of GeneIDs associated with them.
    :param mysql_handler:
    :return:
    """
    geneids_to_translations = dict()
    translation_groups = dict()

    query = "SELECT GeneID, Translation FROM gene"
    results = mysql_handler.execute_query(query)

    print(f"{len(results)} genes in the database")

    for dictionary in results:
        geneid = dictionary["GeneID"]
        translation = dictionary["Translation"]
        # map each geneid to its translation
        geneids_to_translations[geneid] = translation

        # map translations to their corresponding geneids
        translation_group = translation_groups.get(translation, [])
        translation_group.append(geneid)
        translation_groups[translation] = translation_group

    return geneids_to_translations, translation_groups


def write_fasta(trans_groups, dir="/tmp/phamerate"):
    """
    Writes the translations in `trans_dict` to a fasta_file in `dir`.
    :param trans_groups: dictionary mapping translations to their
    associated GeneIDs
    :param dir: the directory to write input.fasta into
    :return:
    """
    print("Begin writing genes to fasta...")
    fasta = open(f"{dir}/input.fasta", "w")
    for translation in trans_groups.keys():
        fasta.write(f">{trans_groups[translation][0]}\n{translation}\n")
    fasta.close()


def create_clusterdb(program="mmseqs", dir="/tmp/phamerate"):
    """
    Decides which program is being used for clustering (mmseqs by
    default) and builds a database from input.fasta appropriate for
    use with the indicated program
    :param program: the program to be used for clustering (currently
    supported options are "mmseqs" or "blast")
    :param dir: temporary directory where all I/O will take place
    :return:
    """
    # If program is blast, make blastdb
    if program == "blast":
        command = f"{BLASTCLUST_PATH}/formatdb -i {dir}/input.fasta " \
                  f"-p T -o T -n {dir}/sequenceDB"

    # Default action is mmseqs
    elif program == "mmseqs":
        command = f"mmseqs createdb {dir}/input.fasta {dir}/sequenceDB"

    # Otherwise (e.g. kClust or other programs not currently supported)
    else:
        print(f"Unknown program {program}")
        return

    print(command)
    # Run the command - assumes mmseqs and makeblastdb are either at the
    # same scope as toplevel that calls these functions, or in $PATH
    with Popen(args=shlex.split(command), stdout=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")


def phamerate(params, program="mmseqs", dir="/tmp/phamerate"):
    """
    Uses inputs to build the command string, then runs it with Popen
    :param params: dictionary of program parameters
    :param program: the program to be used for clustering (currently
    supported options are "mmseqs" or "blast")
    :param dir: temporary directory where all I/O will take place
    :return:
    """
    # If program is blast, create blastclust command string
    if program == "blast":
        base = f"{BLASTCLUST_PATH}/blastclust -i {dir}/input.fasta -d " \
               f"{dir}/sequenceDB -o {dir}/output.txt"

    # Default action is mmseqs command string
    elif program == "mmseqs":
        base = f"mmseqs cluster {dir}/sequenceDB {dir}/clusterDB {dir}"

    # Otherwise (e.g. kClust or other programs not currently supported)
    else:
        print(f"Unknown program {program}")
        return

    # Now add arguments to the base command:
    for param in params.keys():
        base += f" {param} {params[param]}"

    print("Beginning clustering...")
    with Popen(args=shlex.split(base), stdout=PIPE, stderr=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")
        print(process.stderr.read().decode("utf-8") + "\n\n")

    print("Finished clustering...")


def convert_to_parseable(program="mmseqs", dir="/tmp/phamerate"):
    """
    If program is not mmseqs, does nothing. Otherwise it converts
    clusterDB to more parseable output.fasta
    :param program: the program used for clustering (currently
    supported options are "mmseqs" or "blast")
    :param dir: temporary directory where all I/O takes place
    :return:
    """
    # Do nothing if blast - output is already easily parseable
    if program == "blast":
        return

    # For mmseqs we want to convert to ~fasta~ format for easier parsing
    elif program == "mmseqs":
        command = f"mmseqs createseqfiledb {dir}/sequenceDB {dir}/clusterDB " \
                  f"{dir}/output.txt"

    # Otherwise (e.g. kClust or other programs not currently supported)
    else:
        print(f"Unknown program {program}")
        return

    # Run the command if we got here
    with Popen(args=shlex.split(command), stdout=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")


def parse_output(program="mmseqs", dir="/tmp/phamerate"):
    """

    :param program:
    :param dir:
    :return:
    """
    # Variables to store the new pham data
    phams = dict()
    name = 1
    geneids = list()

    print(f"Beginning to parse {program} output...")

    if program == "blast":
        with open(f"{dir}/output.txt", "r") as fh:
            in_reader = csv.reader(fh, delimiter=" ")
            for i, row in enumerate(in_reader):
                geneids = [geneid for geneid in row[:-1]]   # last token is ''
                phams[i + 1] = geneids

    elif program == "mmseqs":
        with open(f"{dir}/output.txt", "r") as fh:
            line = fh.readline()
            # EOF marked by null bit
            while line != "\x00":
                if line.startswith(">"):
                    geneids.append(line.lstrip(">").rstrip("\n").rstrip(" "))
                # If line starts with a null bit, we hit a new pham
                elif line.startswith("\x00"):
                    phams[name] = geneids
                    name += 1
                    geneids = [line.lstrip("\x00").lstrip(">").rstrip(
                        "\n").rstrip(" ")]
                # Otherwise we're on a translation - skip
                else:
                    pass
                line = fh.readline()
            # Dump the last pham
            phams[name] = geneids

    else:
        print(f"Unknown program {program}...")
        return {}

    print(f"Done parsing {program} output.")
    return phams


def reintroduce_duplicates(new_phams, trans_groups, genes_and_trans):
    """

    :param new_phams:
    :param trans_groups:
    :param genes_and_trans:
    :return:
    """
    for key in new_phams.keys():
        geneids = new_phams[key]
        dup_geneids = list()
        for gene in geneids:
            translation = genes_and_trans[gene]
            dup_group = trans_groups[translation]
            for geneid in dup_group:
                dup_geneids.append(geneid)
        new_phams[key] = set(dup_geneids)
    return new_phams


def preserve_phams(old_phams, new_phams, old_colors, new_genes):
    final_phams = dict()
    new_colors = dict()

    outcount = 0
    total = len(old_phams)
    new_phams_1 = new_phams.copy()

    # Iterate through old and new phams
    for old_key in old_phams.keys():
        outcount += 1
        print("Pham Name Conservation: {} / {}".format(outcount, total))

        old_pham = old_phams[old_key]

        if old_key in final_phams.keys():
            continue

        for new_key in new_phams_1:
            new_pham = new_phams_1[new_key]

            if old_pham & new_pham == set():
                continue

            # Case 1 + 5 (Identity and Subtraction)
            if old_pham == new_pham:
                final_phams[old_key] = new_pham
                new_colors[old_key] = old_colors[old_key]
                new_phams_1.pop(new_key, None)
                break

            # Case 2 and 4 (Addition and Join) - PHAM GREW
            elif new_pham - old_pham != set():

                # Case 2 and 4 (Addition and Join)
                if new_pham & new_genes != set():

                    # Case 4 - Join with new gene
                    if (new_pham - (new_pham & new_genes)) - \
                            old_pham \
                            != set():
                        break

                    # Case 2 - Addition with new gene
                    final_phams[old_key] = new_pham
                    new_colors[old_key] = old_colors[old_key]
                    new_phams_1.pop(new_key, None)
                    break

                # Case 4 - Join without new gene
                else:
                    break

            # Case 3 - split - PHAM SHRANK, BUT NOT BY REMOVAL
            elif old_pham - new_pham != set():
                break

    final_phams[0] = "placeholder"
    highest_pham = max(map(int, final_phams.keys())) + 1
    del final_phams[0]

    # Reassign data for split or joined phams
    for key in new_phams_1:
        new_key = highest_pham
        highest_pham += 1

        final_phams[new_key] = new_phams_1[key]

        if len(new_phams_1[key]) > 1:
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
            new_colors[new_key] = hexrgb
        else:
            new_colors[new_key] = '#FFFFFF'

    return final_phams, new_colors


def reinsert_pham_data(new_phams, new_colors, mysql_handler):
    commands = []
    for key in new_phams.keys():
        for gene in new_phams[key]:
            commands.append(f"INSERT INTO pham (GeneID, Name) VALUES ('"
                            f"{gene}', {key})")

    mysql_handler.execute_transaction(commands)

    commands = []
    for key in new_colors.keys():
        commands.append(f"INSERT INTO pham_color (Name, Color) VALUES ({key}, "
                        f"'{new_colors[key]}')")

    mysql_handler.execute_transaction(commands)


def fix_miscolored_phams(mysql_handler):
    print("Phixing Phalsely Hued Phams...")

    query = "SELECT * FROM (SELECT b.ID, COUNT(GeneID) AS count, " \
            "a.Name, b.Color FROM pham AS a INNER JOIN pham_color AS " \
            "b ON a.Name = b.Name GROUP BY a.Name, b.ID) AS c WHERE " \
            "Color = '#FFFFFF' AND count > 1"

    results = mysql_handler.execute_query(query)

    print(f"Found {len(results)} miscolored phams to fix")

    commands = []
    for dictionary in results:
        pham_id = dictionary["ID"]
        count = dictionary["count"]
        name = dictionary["Name"]
        color = dictionary["Color"]
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
        new_color = hexrgb

        commands.append("UPDATE pham_color SET Color = '{}' WHERE "
                        "ID = '{}'".format(new_color, pham_id))

    mysql_handler.execute_transaction(commands)

    print("Phixing Phalsely Phlagged Orphams...")

    query = "SELECT * FROM (SELECT b.ID, COUNT(GeneID) AS count," \
            "a.Name, b.Color FROM pham AS a INNER JOIN pham_color AS " \
            "b ON a.Name = b.Name GROUP BY a.Name, b.ID) AS c WHERE " \
            "Color != '#FFFFFF' AND count = 1"

    results = mysql_handler.execute_query(query)

    print(f"Found {len(results)} miscolored orphams to fix...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["ID"]
        count = dictionary["count"]
        name = dictionary["Name"]
        color = dictionary["Color"]
        new_color = "#FFFFFF"

        commands.append("UPDATE pham_color SET Color = '{}' WHERE "
                        "ID = '{}'".format(new_color, pham_id))

    mysql_handler.execute_transaction(commands)
