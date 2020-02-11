"""Functions that are used in the phameration pipeline"""

import shlex
from subprocess import Popen, PIPE
import csv
import random
import colorsys

from pdm_utils.constants.constants import BLASTCLUST_PATH
from pdm_utils.functions import mysqldb

def get_program_params(program, args):
    program_params = dict()

    if program == "blast":
        program_params = get_blast_params(args)
    elif program == "mmseqs":
        program_params = get_mmseqs_params(args)
    else:
        print(f"Unknown program {program}")

    return program_params


def get_blast_params(args):
    """
    Uses parsed arguments from phamerate.main() to build the blastclust
    commandline argument dictionary
    :param args: dictionary from Argparse.ArgumentParser.parse_args()
    :return: dictionary of blast CLI parameters
    """
    blast_params = dict()
    try:
        blast_params["-S"] = args.identity
        blast_params["-L"] = float(args.coverage)/100
        blast_params["-a"] = args.threads
    except AttributeError as err:
        print(f"Error {err.args[0]}: {err.args[1]}")
    return blast_params


def get_mmseqs_params(args):
    """
    Uses parsed arguments from phamerate.main() to build the MMseqs2
    commandline argument dictionary
    :param args: dictionary from Argparse.ArgumentParser.parse_args()
    :return: dictionary of MMseqs2 CLI parameters
    """
    mmseqs_params = dict()
    try:
        mmseqs_params["--threads"] = args.threads
        mmseqs_params["-v"] = args.verbose
        mmseqs_params["--cluster-steps"] = args.steps
        mmseqs_params["--max-seqs"] = args.max_seqs
        mmseqs_params["--min-seq-id"] = float(args.identity)/100
        mmseqs_params["-c"] = float(args.coverage)/100
        mmseqs_params["--alignment-mode"] = args.aln_mode
        mmseqs_params["--cov-mode"] = args.cov_mode
        mmseqs_params["--cluster-mode"] = args.clu_mode
    except AttributeError as err:
        print(f"Error {err.args[0]}: {err.args[1]}")
    return mmseqs_params


def get_pham_geneids(engine):
    """
    Queries the database for those genes that are already phamerated.
    :param engine: the Engine allowing access to the database
    :return: pham_geneids
    """
    pham_geneids = dict()

    geneid_query = "SELECT GeneID, PhamID FROM gene WHERE PhamID IS NOT NULL"
    geneid_results = mysqldb.query_dict_list(engine, geneid_query)

    print(f"Found {len(geneid_results)} genes in phams...")

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
    color_results = mysqldb.query_dict_list(engine, color_query)

    print(f"Found colors for {len(color_results)} phams...")

    for dictionary in color_results:
        pham_id = dictionary["PhamID"]
        color = dictionary["Color"]

        pham_colors[pham_id] = color

    return pham_colors


def get_new_geneids(engine):
    """
    Queries the database for those genes that are not yet phamerated.
    :param engine: the Engine allowing access to the database
    :return: new_geneids
    """
    new_geneids = set()

    gene_query = "SELECT GeneID FROM gene WHERE PhamID IS NULL"
    gene_results = mysqldb.query_dict_list(engine, gene_query)

    for dictionary in gene_results:
        geneid = dictionary["GeneID"]
        new_geneids = new_geneids | {geneid}

    print(f"Found {len(new_geneids)} genes not in phams...")

    return new_geneids


def map_geneids_to_translations(engine):
    """
    Constructs a dictionary mapping all geneids to their translations.
    :param engine: the Engine allowing access to the database
    :return: gs_to_ts
    """
    gs_to_ts = dict()

    query = "SELECT GeneID, Translation FROM gene"
    results = mysqldb.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        trans = dictionary["Translation"]
        gs_to_ts[geneid] = trans

    print(f"Found {len(results)} genes in the database...")

    return gs_to_ts


def map_translations_to_geneids(engine):
    """
    Constructs a dictionary mapping all unique translations to their
    groups of geneids
    :param engine: the Engine allowing access to the database
    :return: ts_to_gs
    """
    ts_to_gs = dict()

    query = "SELECT GeneID, Translation FROM gene"
    results = mysqldb.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        trans = dictionary["Translation"]
        geneids = ts_to_gs.get(trans, [])
        geneids.append(geneid)
        ts_to_gs[trans] = geneids

    print(f"Found {len(ts_to_gs)} unique translations in the database...")

    return ts_to_gs


def write_fasta(trans_groups, wd):
    """
    Writes the translations in `trans_dict` to a fasta_file in `wd`.
    :param trans_groups: dictionary mapping translations to their
    associated GeneIDs
    :param wd: the temporary working directory for phameration
    :return:
    """
    print("Begin writing genes to fasta...")
    fasta = open(f"{wd}/input.fasta", "w")
    for translation in trans_groups.keys():
        fasta.write(f">{trans_groups[translation][0]}\n{translation}\n")
    fasta.close()


def create_clusterdb(program, wd):
    """
    Decides which program is being used for clustering and calls
    the database constructor method for that program to build a
    database from input.fasta
    :param program: the program to be used for clustering (currently
    supported options are "mmseqs" or "blast")
    :param wd: the temporary working directory for phameration
    :return:
    """
    # If program is MMseqs2 (default), make MMseqs2 database
    if program == "mmseqs":
        command = mmseqsdb_command(wd)
    # If program is blast, make blastdb
    elif program == "blast":
        command = blastdb_command(wd)
    # Otherwise (e.g. kClust or other programs not currently supported)
    else:
        print(f"Unknown program {program}")
        command = "echo No database construction command..."

    # Run the command - assumes mmseqs and makeblastdb are either at the
    # same scope as toplevel that calls these functions, or in $PATH
    with Popen(args=shlex.split(command), stdout=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")


def blastdb_command(wd):
    """
    Builds a blastp database to use for blastclust phameration
    :param wd: the temporary working directory for phameration
    :return:
    """
    # -i: fasta file to be converted to blastdb <filepath>
    # -p: is it a protein? <boolean T/F>
    # -o: parse/index sequence IDs? <boolean T/F>
    # -n: database name </path/to/name>
    # -l: log file <filepath>
    command = f"{BLASTCLUST_PATH}/formatdb -i {wd}/input.fasta " \
              f"-p T -o T -n {wd}/sequenceDB -l {wd}/formatdb.log"

    print("Build blast protein database for blastclust...")

    return command


def mmseqsdb_command(wd):
    """
    Builds an MMseqs2 database to use for MMseqs2 phameration
    :param wd: the temporary working directory for phameration
    :return:
    """
    command = f"mmseqs createdb {wd}/input.fasta {wd}/sequenceDB"

    print("Build MMseqs2 protein database for mmseqs cluster...")

    return command


def phamerate(params, program, wd):
    """
    Phamerates unique translations using specified program
    :param params: dictionary of program parameters
    :param program: the program to be used for clustering (currently
    supported options are "mmseqs" or "blast")
    :param wd: the temporary working directory for phameration
    :return:
    """
    # If program is MMseqs2 (default), create mmseqs command string
    if program == "mmseqs":
        command = mmseqs_phamerate_command(params, wd)
    # If program is blast, create blastclust command string
    elif program == "blast":
        command = blast_phamerate_command(params, wd)
    # Otherwise (e.g. kClust or other programs not currently supported)
    else:
        print(f"Unknown program {program}")
        command = "echo No phameration command..."

    print(f"Begin clustering with {program}...")
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")
        print(process.stderr.read().decode("utf-8") + "\n\n")

    print(f"Finish clustering with {program}...")


def blast_phamerate_command(parameters, wd):
    """
    Builds blast phameration command (base command + user args)
    :param parameters: dictionary of blastclust parameters
    :param wd: the temporary working directory for phameration
    :return: command
    """
    # -d: database to cluster </path/to/name>
    # -o: output file <filepath>
    command = f"{BLASTCLUST_PATH}/blastclust -d {wd}/sequenceDB " \
              f"-o {wd}/output.txt"

    for parameter in parameters.keys():
        command += f" {parameter} {parameters[parameter]}"

    return command


def mmseqs_phamerate_command(parameters, wd):
    """
    Builds MMseqs2 phameration command (base command + user args)
    :param parameters: dictionary of MMseqs2 parameters
    :param wd: the temporary working directory for phameration
    :return: command
    """
    command = f"mmseqs cluster {wd}/sequenceDB {wd}/clusterDB {wd}"

    for parameter in parameters.keys():
        command += f" {parameter} {parameters[parameter]}"

    return command


def parse_output(program, wd):
    """
    Runs the parser appropriate to the specified phameration program
    :param program: program used for clustering
    :param wd: the temporary working directory for phameration
    :return:
    """
    # Variables to store the new pham data
    parsed_phams = dict()

    if program == "mmseqs":
        parsed_phams = parse_mmseqs(wd)
    elif program == "blast":
        parsed_phams = parse_blast(wd)
    else:
        print(f"Unknown program {program}...")

    return parsed_phams


def parse_mmseqs(wd):
    """
    Parses MMseqs2 clustering output by running:
    - mmseqs createseqfiledb sequenceDB clusterDB clusterSeqs
    - mmseqs result2flat sequenceDB sequenceDB clusterSeqs output.txt
    and then parsing output.txt
    :param wd: the temporary working directory for phameration
    :return: parsed_phams (dictionary)
    """
    filename = f"{wd}/output.txt"

    print("Convert MMseqs2 output into parseable format...")

    command = f"mmseqs createseqfiledb {wd}/sequenceDB {wd}/clusterDB " \
              f"{wd}/clusterSeqs"

    with Popen(args=shlex.split(command), stdout=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")

    command = f"mmseqs result2flat {wd}/sequenceDB {wd}/sequenceDB " \
              f"{wd}/clusterSeqs {filename}"

    with Popen(args=shlex.split(command), stdout=PIPE) as process:
        print(process.stdout.read().decode("utf-8") + "\n\n")

    print("Begin parsing MMseqs2 output...")

    parsed_phams = dict()
    pham_geneids = list()
    pham_name = 0

    with open(filename, "r") as fh:
        # Latest MMseqs2 output format indicates the start of a new
        # pham by repeating the pham representative's identifier in
        # two adjacent lines - need references to the prior line as
        # well as the current one.
        prior_line = fh.readline()
        current_line = fh.readline()

        # While not EOF, iterate through lines
        while current_line:
            if current_line.startswith(">"):
                if prior_line.startswith(">"):
                    try:
                        pham_geneids.pop(-1)
                    except IndexError:
                        # First pham should fail because pham_geneids is empty
                        pass
                    parsed_phams[pham_name] = pham_geneids
                    pham_name += 1
                    pham_geneids = [current_line.lstrip(">").rstrip()]
                else:
                    pham_geneids.append(current_line.lstrip(">").rstrip())
            else:
                # Translation, skip
                pass
            # prior gets current's value, current gets new line
            prior_line, current_line = current_line, fh.readline()

        # Dump the last working pham into the dictionary
        parsed_phams[pham_name] = pham_geneids

    # Pham 0 is a placeholder - remove it and then return parsed_phams
    parsed_phams.pop(0)

    print("Finish parsing MMseqs2 output...")
    print(f"Genes were sorted into {len(parsed_phams)} phams...")

    return parsed_phams


def parse_blast(wd):
    """
    Parses blastclust output (space-delimited values)
    :param wd: the temporary working directory for phameration
    :return: parsed_phams (dictionary)
    """
    parsed_phams = dict()

    filename = f"{wd}/output.txt"

    print("Begin parsing blastclust output...")

    with open(filename, "r") as fh:
        in_reader = csv.reader(fh, delimiter=" ")
        for i, row in enumerate(in_reader):
            geneids = [geneid for geneid in row[:-1]]   # last token is ''
            parsed_phams[i + 1] = geneids

    print("Finish parsing blastclust output...")
    print(f"Genes were sorted into {len(parsed_phams)} phams...")

    return parsed_phams


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
        print("Pham Name Conservation: {} / {}".format(outcount, total))

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
            final_colors[new_key] = hexrgb
        else:
            final_colors[new_key] = '#FFFFFF'

    return final_phams, final_colors


def reinsert_pham_data(new_phams, new_colors, engine):
    """
    Puts pham data back into the database
    :param new_phams:
    :param new_colors:
    :param engine:
    :return:
    """
    # Colors have to go first, since PhamID column in gene table references
    # PhamID in pham table
    commands = []
    for key in new_colors.keys():
        commands.append(f"INSERT INTO pham (PhamID, Color) VALUES ({key}, "
                        f"'{new_colors[key]}')")

    mysqldb.execute_transaction(engine, commands)



    commands = []
    for key in new_phams.keys():
        for gene in new_phams[key]:
            commands.append(f"UPDATE gene SET PhamID = {key} WHERE GeneID = '{gene}'")

    mysqldb.execute_transaction(engine, commands)


def fix_miscolored_phams(engine):
    print("Phixing Phalsely Hued Phams...")
    # Phams which are colored as though they are orphams, when really
    # they are multi-member phams
    query = "SELECT * FROM (SELECT g.PhamID, COUNT(GeneID) AS count, "\
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID " \
            "= p.PhamID GROUP BY PhamID) AS c WHERE Color = '#FFFFFF' "\
            "AND count > 1"

    results = mysqldb.query_dict_list(engine, query)


    print(f"Found {len(results)} miscolored phams to fix")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        count = dictionary["count"]
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
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)


    print("Phixing Phalsely Phlagged Orphams...")
    # Phams which are colored as though they are multi-member phams
    # when really they are orphams
    query = "SELECT * FROM (SELECT g.PhamID, COUNT(GeneID) AS count, "\
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID "\
            "=p.PhamID GROUP BY PhamID) AS c WHERE Color != '#FFFFFF' "\
            "AND count = 1"

    results = mysqldb.query_dict_list(engine, query)
    print(f"Found {len(results)} miscolored orphams to fix...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        count = dictionary["count"]
        color = dictionary["Color"]
        new_color = "#FFFFFF"
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)
