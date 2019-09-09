import os
import shutil
import sys
from subprocess import *
import pymysql as pms
import random
import colorsys


class MMseqsPhameratorHandler:
    def __init__(self, mysql_handler):
        """
        This object is intended to handle the phameration of genes in
        a Phamerator database using MMseqs2.  It has methods and
        attributes suitable for: storing extant pham data to conserve
        pham numbers and colors where possible; retrieving geneid and
        translation data and writing them to a file; filtering out
        redundant translations before writing the genes to a file;
        re-associating redundant geneids after phameration; converting
        the FASTA input into the MMseqs2 database format; running
        MMseqs2 with a host of different parameters; parsing MMseqs2
        output; reinserting MMseqs2 pham data into the database, and
        fixing miscolored phams and orphams.
        :param mysql_connection: a pymysql connection object
        """
        # MySQL connection object to be used for querying the database
        self.mysql_handler = mysql_handler

        # Flag to filter genes with duplicate translations, or not
        self.filter_redundant = True

        # File and directory names for use in several methods
        self.temp_dir = "/tmp/MMseqs2"
        self.fasta_input = "temp_input.fasta"
        self.mmseqs_input = "sequenceDB"
        self.mmseqs_output = "clusterDB"
        self.parseable_output = "temp_output.fasta"

        # Store current database status for use at the end
        self.old_phams = {}
        self.old_colors = {}
        self.duplicate_genes = {}
        self.new_genes = set()

        # Store new database data for injection back into the database
        self.new_phams = {}
        self.new_colors = {}

        # Set the MMseqs2 clustering parameter defaults.
        self.threads = 4
        self.verbosity = 3
        self.single_step = False
        self.max_seqs = 250
        self.min_seq_id = 0.8
        self.coverage = 0.8
        self.alignment_mode = 3
        self.coverage_mode = 0
        self.cluster_mode = 0

        # Log file to write output to
        self.log = open("/tmp/MMseqsPhamerationLog.txt", "w")

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

    def get_existing_pham_data(self):
        """
        This function retrieves pham names, GeneIDs, and colors for
        those genes present in the pham/pham_color tables (i.e. genes
        whose phages weren't updated in the present round of updates).
        It then iterates through the returned tuple of tuples and, for
        each pham name (keys in self.old_phams and self.old_colors), if
        the name is present in the existing keys for old_phams, add the
        geneid to that pham.  Otherwise, it's a pham we haven't stored
        information for yet, and we need to store both the current gene
        and the pham color.
        :return:
        """
        try:
            query = "SELECT a.name, a.GeneID, b.color FROM (SELECT p.name," \
                    "g.GeneID FROM gene g INNER JOIN pham p ON g.GeneID = " \
                    "p.GeneID) AS a INNER JOIN pham_color AS b ON a.name = " \
                    "b.name ORDER BY a.Name ASC"
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
                self.old_phams[name] = self.old_phams[name] | set([geneid])
            else:
                self.old_phams[name] = set([geneid])
                self.old_colors[name] = color
        return

    def get_unphamerated_geneids(self):
        """
        This function retrieves the GeneIDs for those genes present in
        the gene table, but not the pham table (i.e. genes whose phages
        were updated in the present round of updates).  It iterates
        through the returned tuple of tuples and adds each GeneID to
        the set of unphamerated genes.
        :return:
        """
        try:
            query = "SELECT GeneID FROM gene WHERE GeneID NOT IN (SELECT " \
                    "g.GeneID FROM gene AS g INNER JOIN pham AS p ON " \
                    "g.GeneID = p.GeneID)"
            result_list = self.mysql_handler.execute_query(query)
        except pms.err.InternalError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)

        print("{} new genes".format(len(result_list)))

        for dictionary in result_list:
            geneid = dictionary["GeneID"]
            self.new_genes = self.new_genes | set([geneid])
        return

    def write_geneids_and_translations_to_fasta(self):
        """
        This function retrieves all GeneIDs and translations in the
        gene table.  It either writes them all to self.fasta_input or
        filters out GeneIDs with duplicated translations and stores
        the duplicate groups in self.duplicate_genes and writes the
        non-redundant translations to self.fasta_input.  It also clears
        the old pham and pham_color data from their respective tables.
        :return:
        """
        # Interact with database.
        try:
            # Clear old pham data - auto commits at end of transaction
            commands = ["TRUNCATE TABLE pham", "TRUNCATE TABLE pham_color"]
            self.mysql_handler.execute_transaction(commands)

            query = "SELECT GeneID, translation FROM gene"
            result_list = self.mysql_handler.execute_query(query)
        except pms.err.Error as err:
            print("Error: {}".format(err))
            sys.exit(1)

        # Write GeneIDs and translations to file
        try:
            print("Writing genes to fasta file")
            f = open("{}/{}".format(self.temp_dir, self.fasta_input), "w")
        except IOError as e:
            print("Error {}: {}".format(e.args[0], e.args[1]))
            sys.exit(1)

        # Write all gene data if self.filter_redundant is set to False
        if self.filter_redundant is False:
            try:
                for dictionary in result_list:
                    f.write(">{}\n".format(dictionary["GeneID"]))
                    f.write("{}\n".format(dictionary["translation"].replace(
                        "-", "M")))
            except IOError as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                sys.exit(1)
        # Otherwise, filter redundant translations, group their geneids
        # and pick a random geneid from the group for the FASTA file.
        else:
            try:
                trans_and_geneids = {}
                for dictionary in result_list:
                    geneids = trans_and_geneids.get(dictionary["translation"],
                                                    set())
                    geneids.add(dictionary["GeneID"])
                    trans_and_geneids[dictionary["translation"]] = geneids

                for trans in trans_and_geneids.keys():
                    geneid = random.sample(trans_and_geneids[trans], 1)[0]
                    self.duplicate_genes[geneid] = trans_and_geneids[trans]
                    f.write(">{}\n".format(geneid))
                    f.write("{}\n".format(trans.replace("-", "M")))
                f.close()
            except IndexError as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                sys.exit(1)
            except IOError as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                sys.exit(1)
        return

    def convert_fasta_to_mmseqsdb(self):
        """
        This function converts the file at self.fasta_input to the db
        format used internally by MMseqs2 (self.mmseqs_input).
        :return:
        """
        try:
            args = ["mmseqs", "createdb", "{}/{}".format(self.temp_dir,
                                                         self.fasta_input),
                    "{}/{}".format(self.temp_dir, self.mmseqs_input)]
            output, errors = Popen(args=args, stdin=PIPE, stdout=PIPE,
                                   stderr=PIPE).communicate()
            self.log.write(str(output))
        except OSError as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)
        return

    def cluster_database(self):
        """
        This function runs MMseqs2 with the defined parameters.
        :return:
        """
        if self.single_step is True:
            try:
                args = ["mmseqs", "cluster", "{}/{}".format(self.temp_dir,
                                                           self.mmseqs_input),
                        "{}/{}".format(self.temp_dir, self.mmseqs_output),
                        self.temp_dir, "--remove-tmp-files",
                        "--single-step-clustering", "--threads",
                        str(self.threads), "-v", str(self.verbosity),
                        "--max-seqs", str(self.max_seqs), "--min-seq-id",
                        str(self.min_seq_id), "-c", str(self.coverage),
                        "--alignment-mode", str(self.alignment_mode),
                        "--cov-mode", str(self.coverage_mode),
                        "--cluster-mode", str(self.cluster_mode)]
                output, errors = Popen(args=args, stdin=PIPE, stdout=PIPE,
                                       stderr=PIPE).communicate()
                self.log.write(str(output))
            except OSError as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                sys.exit(1)
        else:
            try:
                args = ["mmseqs", "cluster", "{}/{}".format(self.temp_dir,
                                                           self.mmseqs_input),
                        "{}/{}".format(self.temp_dir, self.mmseqs_output),
                        self.temp_dir, "--remove-tmp-files", "--threads",
                        str(self.threads), "-v", str(self.verbosity),
                        "--max-seqs", str(self.max_seqs), "--min-seq-id",
                        str(self.min_seq_id), "-c", str(self.coverage),
                        "--alignment-mode", str(self.alignment_mode),
                        "--cov-mode", str(self.coverage_mode),
                        "--cluster-mode", str(self.cluster_mode)]
                output, errors = Popen(args=args, stdin=PIPE, stdout=PIPE,
                                       stderr=PIPE).communicate()
                self.log.write(str(output))
            except OSError as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                sys.exit(1)
        return

    def convert_and_parse_output(self):
        """
        This function converts the MMseqs2 output into a format similar
        to the Pearson FASTA format.  It then parses this output, with
        the knowledge that null bits (\x00) signify new clusters or if
        found alone in a line, that the end of file is here.
        :return:
        """
        # Convert the output
        try:
            args = ["mmseqs", "createseqfiledb", "{}/{}".format(
                self.temp_dir, self.mmseqs_input), "{}/{}".format(
                self.temp_dir, self.mmseqs_output), "{}/{}".format(
                self.temp_dir, self.parseable_output)]
            output, errors = Popen(args=args, stdin=PIPE, stdout=PIPE,
                                   stderr=PIPE).communicate()
            self.log.write(str(output))
        except OSError as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)

        # Parse the output
        try:
            f = open("{}/{}".format(self.temp_dir, self.parseable_output), "r")
            line = f.readline()

            phams = {}
            number = 1
            geneids = []

            while line != "\x00":  # While end of file hasn't been reached
                if line[0] == ">":  # FASTA header, grab GeneID
                    geneids.append(line[1:].rstrip("\n").rstrip(" "))
                elif "\x00" in line:  # New pham
                    phams[number] = set(geneids)  # Dump current geneids
                    number += 1  # Increment pham number
                    geneids = [line[2:].rstrip("\n").rstrip(" ")]  # GeneID
                else:
                    pass  # Skip translation lines
                line = f.readline()  # Get next line
            phams[number] = set(geneids)  # Dump the last pham
            f.close()
        except IOError as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            sys.exit(1)

        # Write the new_phams attribute with these new pham compositions
        self.new_phams = phams

        return

    def reintroduce_duplicates(self):
        """
        This function overwrites each pham in self.new_phams with the
        set of all GeneIDs sharing translations with any member of the
        pham in question
        :return:
        """
        if self.filter_redundant is True:
            for key in self.new_phams.keys():
                geneids = self.new_phams[key]
                new_geneids = []
                for gene in geneids:
                    duplicate_group = self.duplicate_genes[gene]
                    for geneid in duplicate_group:
                        new_geneids.append(geneid)
                self.new_phams[key] = set(new_geneids)
        else:
            pass
        return

    def conserve_pham_names(self):
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
