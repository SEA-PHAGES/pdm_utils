"""
This class handles start-to-end phameration of the translations present
in a Phamerator database. Phameration can proceed using either MMseqs2
(default) or blastclust. MMseqs2 is always faster than blastclust by at
least an order of magnitude, but clustering specificity, sensitivity,
and accuracy vary based on chosen parameters.
"""

import os
import shutil
import sys
from subprocess import Popen, PIPE
import shlex
import pymysql as pms
import random
import colorsys
import csv

from pdm_utils.constants.constants import BLASTCLUST_PATH


class PhameratorHandler:
    def __init__(self, mysql_handler, parsed_args):
        """
        :param mysql_handler: a MySQLConnectionHandler object
        """
        # MySQL connection object to be used for querying the database
        self.mysql_handler = mysql_handler

        # Flag to denote whether redundant translations should be filtered
        self.filter_redundant = True

        # Flag to direct program flow at stages where MMseqs2 and blastclust
        # differ in their requirements
        self.use_blast = parsed_args.use_blast

        # I/O attributes
        self.temp_dir = "/tmp/phamerate"
        self.query_fasta = "phamerate_input.fasta"

        self.mmseqs_in = "sequenceDB"
        self.mmseqs_out = "clusterDB"
        self.mmseqs_parseable_out = "phamerate_output.fasta"

        self.blast_db = "targetDB"
        self.blast_out = "phamerate_output.txt"

        # Clustering parameters (method populates them based on use_blast flag)
        self.cluster_params = self.get_default_params(parsed_args)

        # Store current database status for use at the end
        self.old_phams = dict()
        self.old_colors = dict()
        self.geneids_to_translations = dict()
        self.duplicate_translations = dict()
        self.new_genes = set()

        # Store new database data for injection back into the database
        self.new_phams = dict()
        self.new_colors = dict()

        # Log file to write output to
        self.log = open("/tmp/MMseqsPhamerationLog.txt", "w")

    def main(self):
        """
        Main workflow for phameration (doesn't care which parameters
        or program are to be used).
        :return:
        """
        # Refresh temp dir (essential for MMseqs2 - otherwise it will use old
        # outputs as "new" and run super fast but also potentially with the
        # wrong output
        self.refresh_temp_dir()

        # Get old pham data
        self.read_existing_phams()
        self.read_unphamerated_genes()

        # Write fasta input and create appropriate cluster database
        self.write_fasta()
        self.create_cluster_db()

        # Phamerate
        self.phamerate()

        # Read output, reintroduce duplicates
        self.read_new_phams()
        self.reintroduce_duplicates()

        # Preserve pham names/colors and reinsert
        self.preserve_pham_names()
        self.reinsert_pham_data()

        # Fix miscolored phams
        self.fix_miscolored_phams()

    @staticmethod
    def get_default_params(parsed_args):
        """
        Determines which default parameters to use depending on the
        use_blast flag's value. MMseqs2 and blastclust have different
        parameters (some value overlap but names are all different).
        :param parsed_args: parsed command line arguments list
        :type parsed_args: Namespace object from argparse
        :return: params
        """
        if parsed_args.use_blast:
            params = {"-S": parsed_args.identity,
                      "-L": float(parsed_args.coverage)/100,
                      "-a": parsed_args.threads}
        else:
            params = {"--threads": parsed_args.threads,
                      "-v": parsed_args.verbose,
                      "--cluster-steps": parsed_args.steps,
                      "--max-seqs": parsed_args.max_seqs,
                      "--min-seq-id": float(parsed_args.identity)/100,
                      "-c": float(parsed_args.coverage)/100,
                      "--alignment-mode": parsed_args.aln_mode,
                      "--cov-mode": parsed_args.cov_mode,
                      "--cluster-mode": parsed_args.clu_mode}
        return params

    def refresh_temp_dir(self):
        """
        This function checks to see if there is already a directory
        at self.temp_dir.  If there is, it is from an earlier run
        and should be recursively removed so as to prevent MMseqs2
        from using an earlier output (which is its default action).
        If self.temp_dir does not exist, it creates the directory.
        :return:
        """
        if os.path.exists(self.temp_dir) is True:
            try:
                shutil.rmtree(self.temp_dir)
                os.mkdir(self.temp_dir)
                return
            except OSError as e:
                print("Error {}: {}".format(e.args[0], e.args[1]))
                sys.exit(1)
        else:
            try:
                os.mkdir(self.temp_dir)
                return
            except OSError as e:
                print("Error {}: {}".format(e.args[0], e.args[1]))
                sys.exit(1)

    def read_existing_phams(self):
        """
        Retrieves pham names, GeneIDs, and colors for genes present in
        the pham/pham_color tables (those genes whose phages were NOT
        updated in the current set of database updates).
        :return:
        """
        try:
            query = "SELECT * FROM (SELECT a.name, a.GeneID, b.color FROM " \
                    "pham AS a INNER JOIN pham_color AS b on a.name = " \
                    "b.name) AS c ORDER BY c.name ASC"
            result_list = self.mysql_handler.execute_query(query)
        except pms.err.InternalError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)

        print("{} old genes".format(len(result_list)))

        for dictionary in result_list:
            name = dictionary["name"]
            geneid = dictionary["GeneID"]
            color = dictionary["color"]
            if name in self.old_phams.keys():
                self.old_phams[name] = self.old_phams[name] | {geneid}
            else:
                self.old_phams[name] = {geneid}
                self.old_colors[name] = color

        return

    def read_unphamerated_genes(self):
        """
        Retrieves GeneIDs for those genes present in the gene table but
        not the pham table (genes whose phages WERE updated in the
        current set of database updates.
        :return:
        """
        try:
            query = "SELECT GeneID FROM gene WHERE GeneID NOT IN (SELECT " \
                    "GeneID from pham)"
            result_list = self.mysql_handler.execute_query(query)
        except pms.err.InternalError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)

        print("{} new genes".format(len(result_list)))

        for dictionary in result_list:
            geneid = dictionary["GeneID"]
            self.new_genes = self.new_genes | {geneid}

        return

    def write_fasta(self):
        """
        Retrieves all GeneIDs and translations in the gene table.
        Filters duplicates if it's supposed to, then writes phameration
        input FASTA file. Clears old pham/color data from their tables.
        :return:
        """
        try:
            # Clear old pham data - auto commits at end of transaction
            commands = ["TRUNCATE TABLE pham", "TRUNCATE TABLE pham_color"]
            self.mysql_handler.execute_transaction(commands)

            query = "SELECT GeneID, translation FROM gene"
            result_list = self.mysql_handler.execute_query(query)
        except pms.err.Error as err:
            print("Error: {}".format(err))
            sys.exit(1)

        # Parse geneids and translations - sort into duplicate translation
        # groups
        for dictionary in result_list:
            geneid = dictionary["GeneID"]
            translation = dictionary["translation"].replace("-", "M")
            self.geneids_to_translations[geneid] = translation

            duplicate_group = self.duplicate_translations.get(translation, [])
            duplicate_group.append(geneid)

            self.duplicate_translations[translation] = duplicate_group

        # Write GeneIDs and translations to file
        print("Begin write genes to fasta...")
        try:
            f = open("{}/{}".format(self.temp_dir, self.query_fasta), "w")
            # Write all genes if not filtering redundant
            if not self.filter_redundant:
                for geneid in self.geneids_to_translations.keys():
                    f.write(">{}\n{}\n".format(
                        geneid, self.geneids_to_translations[geneid]))
            else:
                for translation in self.duplicate_translations.keys():
                    f.write(">{}\n{}\n".format(
                        self.duplicate_translations[translation][0],
                        translation))
        except IOError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)
        except IndexError as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)

        print("Finish write genes to fasta...")

        return

    def create_cluster_db(self):
        """
        Converts query_fasta to the appropriate database format for
        the clustering program to be used.
        :return:
        """
        # If we'll be running blastclust, make a blastdb
        if len(self.cluster_params) == 3:
            command = "makeblastdb -in {} -out {} -dbtype prot " \
                      "-parse_seqids".format(self.query_fasta, self.blast_db)

        # Else, make an mmseqsdb
        else:
            command = "mmseqs createdb {}/{} {}/{}".format(
                self.temp_dir, self.query_fasta, self.temp_dir, self.mmseqs_in)

        with Popen(args=shlex.split(command), stdout=PIPE) as process:
            self.log.write(process.stdout.read().decode("utf-8") + "\n\n")

        return

    def phamerate(self):
        """
        Runs the appropriate clustering program with user-defined
        parameters.
        :return:
        """
        # If only 3 parameters, use blastclust
        if len(self.cluster_params) == 3:
            base = "{}/bin/blastclust -i {}/{} -d {}/{} -o {}/{}".format(
                str(BLASTCLUST_PATH), self.temp_dir, self.query_fasta,
                self.temp_dir, self.blast_db, self.temp_dir, self.blast_out)
        # Otherwise use mmseqs2
        else:
            base = "mmseqs cluster {}/{} {}/{} {}".format(
                self.temp_dir, self.mmseqs_in, self.temp_dir,
                self.mmseqs_out, self.temp_dir)

        # Parameters are already sorted out properly so we can just add them
        # to the base command string
        for param in self.cluster_params.keys():
            base += " {} {}".format(param, self.cluster_params[param])

        # Print starting clustering
        print("Beginning clustering...")

        with Popen(args=shlex.split(base), stdout=PIPE, stderr=PIPE) as process:
            self.log.write(process.stdout.read().decode("utf-8") + "\n\n")
            self.log.write(process.stderr.read().decode("utf-8") + "\n\n")

        print("Finish clustering...")

        return

    def read_new_phams(self):
        """
        Reads phameration output (for mmseqs it converts output to more
        parseable format first).
        :return:
        """
        # Variables to temporarily store the new pham data
        phams = {}
        name = 1
        geneids = []

        # blastclust output reading is easy, try that first
        if len(self.cluster_params) == 3:
            print("Begin reading blastclust phams...")
            with open(self.blast_out, "r") as fh:
                in_reader = csv.reader(fh, delimiter=" ")
                for i, row in enumerate(in_reader):
                    geneids = [geneid for geneid in row[:-1]]
                    phams[i + 1] = geneids
            print("Finish reading blastclust phams...")

        # mmseqs more complex - convert output first
        else:
            print("Begin converting mmseqs output...")
            command = "mmseqs createseqfiledb {}/{} {}/{} {}/{}".format(
                self.temp_dir, self.mmseqs_in, self.temp_dir,
                self.mmseqs_out, self.temp_dir, self.mmseqs_parseable_out)

            with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) \
                    as process:
                self.log.write(process.stdout.read().decode("utf-8") + "\n\n")
                self.log.write(process.stderr.read().decode("utf-8") + "\n\n")

            print("Finish converting mmseqs output...")
            print("Begin reading mmseqs phams...")

            with open("{}/{}".format(
                    self.temp_dir, self.mmseqs_parseable_out, "r")) as fh:
                line = fh.readline()
                # EOF is marked by null bit
                while line != "\x00":
                    # If line is FASTA header, grab GeneID
                    if line.startswith(">"):
                        geneids.append(line.lstrip(">").rstrip("\n").rstrip(" "))
                    # If line starts with null bit, new pham encountered
                    elif line.startswith("\x00"):
                        phams[name] = set(geneids)
                        name += 1
                        geneids = [line.lstrip("\x00").lstrip(">").rstrip(
                            "\n").rstrip(" ")]
                    # Otherwise we're on a translation - skip
                    else:
                        pass
                    line = fh.readline()
                # Dump the last pham
                phams[name] = set(geneids)

            print("Finish reading mmseqs phams...")

        self.new_phams = phams

        return

    def reintroduce_duplicates(self):
        """
        This function overwrites each pham in self.new_phams with the
        set of all GeneIDs sharing translations with any member of the
        pham in question
        :return:
        """
        # Only have something to do if filter_redundant is set to True
        if self.filter_redundant is True:
            for key in self.new_phams.keys():
                geneids = self.new_phams[key]
                new_geneids = []
                for gene in geneids:
                    # Look up the translation and duplicate group
                    translation = self.geneids_to_translations[gene]
                    duplicate_group = self.duplicate_translations[translation]
                    for geneid in duplicate_group:
                        new_geneids.append(geneid)
                # Overwrite the pham with expanded geneid set
                self.new_phams[key] = set(new_geneids)

        return

    def preserve_pham_names(self):
        """
        This function overwrites self.new_phams with the conserved pham
        names, once they're figured out.  It also populates
        self.new_colors with conserved colorimetric data for matched
        phams, and new colors for phams that couldn't be conserved.
        :return:
        """
        final_phams = {}
        phams_colors = {}

        new_phams = self.new_phams.copy()
        outcount = 0
        total = len(self.old_phams)

        # Iterate through old and new phams
        for old_key in self.old_phams.keys():
            outcount += 1
            print("Pham Name Conservation: {} / {}".format(outcount, total))

            old_pham = self.old_phams[old_key]

            if old_key in final_phams.keys():
                continue

            for new_key in new_phams:
                new_pham = new_phams[new_key]

                if old_pham & new_pham == set():
                    continue

                # Case 1 + 5 (Identity and Subtraction)
                if old_pham == new_pham:
                    final_phams[old_key] = new_pham
                    phams_colors[old_key] = self.old_colors[old_key]
                    new_phams.pop(new_key, None)
                    break

                # Case 2 and 4 (Addition and Join) - PHAM GREW
                elif new_pham - old_pham != set():

                    # Case 2 and 4 (Addition and Join)
                    if new_pham & self.new_genes != set():

                        # Case 4 - Join with new gene
                        if (new_pham - (new_pham & self.new_genes)) - \
                                old_pham \
                                != set():
                            break

                        # Case 2 - Addition with new gene
                        final_phams[old_key] = new_pham
                        phams_colors[old_key] = self.old_colors[old_key]
                        new_phams.pop(new_key, None)
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
        for key in new_phams:
            new_key = highest_pham
            highest_pham += 1

            final_phams[new_key] = new_phams[key]

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
                hexrgb = "#{:02x}{:02x}{:02x}".format(int(rgb[0]), int(rgb[1]),
                                                      int(rgb[2]))
                phams_colors[new_key] = hexrgb
            else:
                phams_colors[new_key] = '#FFFFFF'

        self.new_phams = final_phams
        self.new_colors = phams_colors
        return

    def reinsert_pham_data(self):
        """
        This function reinserts pham name and GeneID pairs into the
        pham table, and also puts pham name and color pairs into the
        pham_color table.
        :return:
        """
        try:
            commands = []
            print("Inserting Pham Data")
            for key in self.new_phams.keys():
                for gene in self.new_phams[key]:
                    commands.append("INSERT INTO pham (GeneID, name) VALUES "
                                    "('{}', {})".format(gene, key))
            self.mysql_handler.execute_transaction(commands)

            commands = []
            print("Inserting Color Data")
            for key in self.new_colors.keys():
                commands.append("INSERT INTO pham_color (name, color) VALUES "
                                "({}, '{}')".format(key, self.new_colors[key]))
            self.mysql_handler.execute_transaction(commands)
        except pms.err.Error as err:
            print("Error: {}".format(err))
            sys.exit(1)
        return

    def fix_miscolored_phams(self):
        """
        This function calculates new colors for phams where the color
        is white but the pham is not an orpham, or for phams where the
        color is not white but the pham is an orpham.  The new colors
        are reinserted into the pham_color table, if any miscolored
        phams and orphams were found.
        :return:
        """
        try:
            print("Phixing Phalsely Hued Phams...")
            query = "SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count, " \
                    "a.name, b.color FROM pham AS a INNER JOIN pham_color AS " \
                    "b ON a.name = b.name GROUP BY a.name, b.id) AS c WHERE " \
                    "color = '#FFFFFF' AND count > 1"
            result_list = self.mysql_handler.execute_query(query)
            print("Found {} miscolored phams to fix".format(len(result_list)))

            commands = []
            for dictionary in result_list:
                pham_id = dictionary["id"]
                count = dictionary["count"]
                name = dictionary["name"]
                color = dictionary["color"]
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

                commands.append("UPDATE pham_color SET color = '{}' WHERE "
                                "id = '{}'".format(new_color, pham_id))
            self.mysql_handler.execute_transaction(commands)

            print("Phixing Phalsely Phlagged Orphams...")
            query = "SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count," \
                    "a.name, b.color FROM pham AS a INNER JOIN pham_color AS " \
                    "b ON a.name = b.name GROUP BY a.name, b.id) AS c WHERE " \
                    "color != '#FFFFFF' AND count = 1"
            result_list = self.mysql_handler.execute_query(query)
            print("Found {} miscolored orphams to fix...".format(len(result_list)))

            commands = []
            for dictionary in result_list:
                pham_id = dictionary["id"]
                count = dictionary["count"]
                name = dictionary["name"]
                color = dictionary["color"]
                new_color = "#FFFFFF"

                commands.append("UPDATE pham_color SET color = '{}' WHERE "
                                "id = '{}'".format(new_color, pham_id))
            self.mysql_handler.execute_transaction(commands)

        except pms.err.Error as err:
            print("Error: {}".format(err))
