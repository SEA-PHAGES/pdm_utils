import os
import pymysql as pms
import sys
import random
import colorsys

from misc_functions import *


def reset_mmseqs_tempdir(tempdir="/tmp/MMseqs2/"):
    """
    This function recursively deletes existing /tmp/MMseqs2 directory
    if it exists, and creates a new one.
    :param tempdir: /tmp/MMseqs2 unless otherwise specified
    :return tempdir: refreshed (new, empty) tempdir
    """
    if os.path.exists(tempdir) is True:
        try:
            os.system("rm -r {}".format(tempdir))
            os.mkdir(tempdir)
            return tempdir
        except OSError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)
    else:
        try:
            os.mkdir(tempdir)
            return tempdir
        except OSError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)


def get_phamerated_data(connection):
    """
    This function creates a cursor for the indicated connection, and
    queries the database for pham names, geneids, and colors associated
    with already phamerated genes in the database.
    :param connection: the MySQL connection to use
    :return tuples: tuple of tuples containing query results
    """
    try:
        cur = connection.cursor()
        cur.execute("SELECT a.name, a.GeneID, b.color FROM (SELECT p.name, "
                    "g.GeneID FROM gene g INNER JOIN pham p ON g.GeneID = "
                    "p.GeneID) AS a INNER JOIN pham_color AS b ON a.name = "
                    "b.name ORDER BY a.Name ASC")
        tuples = cur.fetchall()
        cur.close()
        return tuples
    except pms.err as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)


def get_unphamerated_data(connection):
    """
    This function creates a cursor for the indicated connection, and
    queries the database for the geneids associated with un-phamerated
    genes in the database.
    :param connection: the MySQL connection to use
    :return tuples: tuple of tuples containing query results
    """
    try:
        cur = connection.cursor()
        cur.execute("SELECT GeneID FROM gene WHERE GeneID NOT IN (SELECT "
                    "g.GeneID FROM gene AS g INNER JOIN pham AS p ON "
                    "g.GeneID =  p.GeneID)")
        tuples = cur.fetchall()
        cur.close()
        return tuples
    except pms.err as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)


def get_geneids_and_translations(connection):
    """
    This function creates a cursor for the indicated connection, and
    uses it to clear old pham data, then query for all geneids and
    translations.
    :param connection: the MySQL connection to use
    :return tuples: tuple of tuples containing query results
    """
    try:
        cur = connection.cursor()
        # Clear previous pham data.
        cur.execute("TRUNCATE TABLE pham")
        cur.execute("TRUNCATE TABLE pham_color")
        cur.execute("COMMIT")

        # Now get geneids and translations.
        cur.execute("SELECT GeneID, translation FROM gene")
        tuples = cur.fetchall()
        return tuples
    except pms.err as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)


def write_mmseqs_input(tuple_of_tuples):
    """
    This function takes an input dictionary of geneids and translations
    and filters out duplicate protein sequences.  It writes the non-
    redundant protein sequences to a fasta-formatted file in
    /tmp/MMseqs2, choosing a random geneid from the list of genes with
    identical translations.  It returns the file handle for the written
    file, and a dictionary associating duplicated translations with all
    their relevant geneids for post-phameration propagation of pham
    numbers to all their geneids.
    :param tuple_of_tuples: a tuple containing tuples of geneids and
    translations
    :return fh: file handle for the MMseqs2 input file
    :return duplicate_groups: dictionary whose keys are duplicated
    translations, and values are lists of geneids sharing a translation
    """
    # Create file to be written.
    try:
        fh = "/tmp/MMseqs2/tempquery.txt"
        f = open(fh, "w")
    except IOError as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)

    # Write geneids and translations to file.
    try:
        # Run without filtering duplicate protein sequences.
        # for t in tuples:
        #   f.write(">" + t[0] + '\n')
        #   f.write(t[1].replace('-','M') + '\n')

        # Filter out duplicate sequences, then write to file.
        translations_and_geneids = {}
        for t in tuple_of_tuples:
            geneids = translations_and_geneids.get(t[1], set())
            geneids.add(t[0])
            translations_and_geneids[t[1]] = geneids

        duplicate_groups = {}

        for trans in translations_and_geneids.keys():
            geneid = random.sample(translations_and_geneids[trans],1)[0]
            duplicate_groups[geneid] = translations_and_geneids[trans]
            f.write(">{}\n{}\n".format(geneid, trans.replace('-', 'M')))
        f.close()
    except IndexError as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)
    except IOError as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)
    return fh, duplicate_groups


def convert_to_mmseqs(fasta_handle, mmseqs_handle="/tmp/MMseqs2/sequenceDB"):
    """
    This function converts an existing fasta-formatted file to MMseqs2
    internal database format so that we can cluster the database.
    :param fasta_handle: the path to the MMseqs2 input fasta file
    :return mmseqs_path: the path to the MMseqs2 formatted database
    file
    """
    # Verify paths:
    fasta_exists = verify_path(fasta_handle)
    mmseqs_exists = verify_path(mmseqs_handle)

    if fasta_exists is True and mmseqs_exists is False:
        command = "mmseqs createdb {} {}".format(fasta_handle, mmseqs_handle)
        try:
            os.system(command)
        except OSError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)
    elif fasta_exists is False and mmseqs_exists is False:
        print("{} does not exist...".format(fasta_handle))
        sys.exit(1)
    elif fasta_exists is True and mmseqs_exists is True:
        print("{} already exists, can't proceed".format(mmseqs_handle))
        sys.exit(1)
    elif fasta_exists is False and mmseqs_exists is True:
        print("{} doesn't exist and {} already exists".format(
            fasta_handle, mmseqs_handle))
        sys.exit(1)
    return mmseqs_handle


def cluster_mmseqs_database(mmseqs_input,
                            mmseqs_output="/tmp/MMseqs2/clusterDB",
                            tmp_dir="/tmp/MMseqs2/",
                            threads=4,
                            verbosity=3,
                            single_step=True,
                            max_seqs=1000,
                            min_seq_id=0.4,
                            cov=0.8,
                            alignment_mode=3,
                            cov_mode=0,
                            cluster_mode=0):
    try:
        if single_step is True:
            command = "mmseqs cluster {} {} {} --remove-tmp-files --threads " \
                      "{} -v {} --single-step-clustering --max-seqs {} " \
                      "--min-seq-id {} -c {} --alignment-mode {} --cov-mode " \
                      "{} --cluster-mode {}".format(mmseqs_input,
                                                    mmseqs_output, tmp_dir,
                                                    threads, verbosity,
                                                    max_seqs, min_seq_id,
                                                    cov, alignment_mode,
                                                    cov_mode, cluster_mode)
        else:
            command = "mmseqs cluster {} {} {} --remove-tmp-files --threads " \
                      "{} -v {} --max-seqs {} --min-seq-id {} -c {} " \
                      "--alignment-mode {} --cov-mode {} --cluster-mode {}" \
                      "".format(mmseqs_input, mmseqs_output, tmp_dir,
                                threads, verbosity, max_seqs, min_seq_id,
                                cov, alignment_mode, cov_mode, cluster_mode)
        os.system(command)
        return mmseqs_output
    except OSError as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)


def convert_output(mmseqs_input, mmseqs_output,
                   readable_output="/tmp/MMseqs2/output.txt"):
    try:
        command = "mmseqs createseqfiledb {} {} {}".format(mmseqs_input,
                                                           mmseqs_output,
                                                           readable_output)
        os.system(command)
        return readable_output
    except OSError as e:
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)


def parse_output(readable_output):
    f = open(readable_output, "r")
    line = f.readline()

    phams = {}
    number = 1
    geneids = []

    while line != "\x00":
        if line[0] == ">":
            geneids.append(line[1:].rstrip("\n").rstrip(" "))
        elif "\x00" in line:
            phams[number] = set(geneids)
            number += 1
            geneids = [line[2:].rstrip("\n").rstrip(" ")]
        else:
            pass
        line = f.readline()
    phams[number] = set(geneids)
    f.close()

    return phams


def reintroduce_duplicates(phams, duplicates):
    for key in phams.keys():
        geneids = phams[key]
        new_geneids = []
        for gene in geneids:
            group = duplicates[gene]
            for geneid in group:
                new_geneids.append(geneid)
        phams[key] = set(new_geneids)

    return phams


def conserve_pham_names(new_phams, old_phams, old_colors, new_set):
    phams_final = {}
    phams_colors = {}

    new_phams_temp = new_phams.copy()
    outcount = 0
    total = len(old_phams)

    # Iterate through old and new phams
    for old_key in old_phams:
        outcount += 1
        print("Pham Name Conservation Search: {} / {}".format(outcount, total))

        old_pham = old_phams[old_key]

        if old_key in phams_final.keys():
            continue

        for new_key in new_phams_temp:
            new_pham = new_phams_temp[new_key]

            if old_pham & new_pham == set():
                continue

            # Case 1 + 5 (Identity and Subtraction)
            if old_pham == new_pham:
                phams_final[old_key] = new_pham
                phams_colors[old_key] = old_colors[old_key]
                new_phams.pop(new_key, None)
                break

            # Case 2 and 4 (Addition and Join) - PHAM GREW
            elif new_pham - old_pham != set():

                # Case 2 and 4 (Addition and Join)
                if new_pham & new_set != set():

                    # Case 4 - Join with new gene
                    if (new_pham - (new_pham & new_set)) - old_pham != set():
                        break

                    # Case 2 - Addition with new gene
                    phams_final[old_key] = new_pham
                    phams_colors[old_key] = old_colors[old_key]
                    new_phams.pop(new_key, None)
                    break

                # Case 4 - Join without new gene
                else:
                    break

            # Case 3 - split - PHAM SHRANK, BUT NOT BY REMOVAL
            elif old_pham - new_pham != set():
                break

    phams_final[0] = "placeholder"
    highest_pham = max(map(int, phams_final.keys())) + 1
    del phams_final[0]

    # Reassign data for split or joined phams
    for key in new_phams:
        new_key = highest_pham
        highest_pham += 1

        phams_final[new_key] = new_phams[key]

        if len(new_phams[key]) > 1:
            h = s = v = 0
            while h <= 0:
                h = random.random()
            while s < 0.5:
                s = random.random()
            while v < 0.8:
                v = random.random()
            rgb = colorsys.hsv_to_rgb(h, s, v)
            rgb = (rgb[0] * 255, rgb[1] * 255, rgb[2] * 255)
            hexrgb = '#{}02x{}02x{}02x'.format(int(rgb[0]), int(rgb[1]),
                                               int(rgb[2]))
            phams_colors[new_key] = hexrgb
        else:
            phams_colors[new_key] = '#FFFFFF'

    return phams_final, phams_colors


def reinsert_pham_data(connection, pham_data, color_data):
    cur = connection.cursor()
    try:
        print("Inserting Pham Data")

        for key in pham_data:
            for gene in pham_data[key]:
                cur.execute("INSERT INTO pham(GeneID, Name) VALUES ({}, "
                            "{})".format(gene, key))
                cur.execute("COMMIT")

        print("Inserting Color Data")

        for key in color_data.keys():
            cur.execute("INSERT IGNORE INTO pham_color(Name, color) VALUES ("
                        "{}, {})".format(key, color_data[key]))
            cur.execute("COMMIT")
    except pms.err as e:
        cur.close()
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)
    cur.close()
    return True


def fix_miscolored_phams(connection):
    cur = connection.cursor()
    try:
        print("Phixing Phalsely Hued Phams...")
        cur.execute("SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count, "
                    "a.Name, b.color FROM pham AS a INNER JOIN pham_color AS "
                    "b ON a.Name = b.Name GROUP BY a.Name, b.id) AS c WHERE "
                    "color = '#FFFFFF' AND count > 1")
        tuples = cur.fetchall()
        print("Found {} miscolored phams to fix...".format(len(tuples)))

        for t in tuples:
            pham_id, count, name, color = t
            h = s = v = 0
            while h <= 0:
                h = random.random()
            while s < 0.5:
                s = random.random()
            while v < 0.8:
                v = random.random()
            rgb = colorsys.hsv_to_rgb(h, s, v)
            rgb = (rgb[0]*255, rgb[1]*255, rgb[2]*255)
            hexrgb = "#%02x%02x%02x" % rgb
            new_color = hexrgb

            cur.execute("UPDATE pham_color SET color = '{}' WHERE id = "
                        "'{}'".format(new_color, pham_id))
    except pms.err as e:
        cur.close()
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)
    cur.close()
    return


def fix_miscolored_orphams(connection):
    cur = connection.cursor()
    try:
        print("Phixing Phalsely Phlagged Orphams")
        cur.execute("SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count, "
                    "a.Name, b.color FROM pham AS a INNER JOIN pham_color AS "
                    "b ON a.Name = b.Name GROUP BY a.Name, b.id) AS c WHERE "
                    "color != '#FFFFFF' AND count = 1")
        tuples = cur.fetchall()
        print("Found {} miscolored phams to fix...")

        for t in tuples:
            pham_id, count, name, color = t
            new_color = "#FFFFFF"

            cur.execute("UPDATE pham_color SET color = '{}' WHERE id = "
                        "'{}'".format(new_color, pham_id))
    except pms.err as e:
        cur.close()
        print("Error {}: {}".format(e.args[0], e.args[1]))
        sys.exit(1)
    cur.close()
    return
