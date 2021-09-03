import os
import shlex
from subprocess import Popen
import re

from pdm_utils.functions import basic


class AragornHandler:
    def __init__(self, identifier, sequence):
        self.id = identifier
        self.sequence = sequence

        # I/O attributes
        self.temp_dir = "/tmp/aragorn"
        self.input = os.path.join(self.temp_dir, self.id) + ".fasta"
        self.output = os.path.join(self.temp_dir, self.id) + ".out"
        self.out_str = ""

        self.trna_tally = 0
        self.tmrna_tally = 0

        self.trnas = list()
        self.tmrnas = list()

    def write_fasta(self):
        """
        Writes the search sequence to input file in FASTA format.
        :return:
        """
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
        with open(self.input, "w") as fh:
            fh.write(f">{self.id}\n{self.sequence}\n")

    def run_aragorn(self, c=False, d=True, m=False, t=True):
        """
        Set up Aragorn command, then run it. Default arguments will
        assume linear sequence to be scanned on both strands for tRNAs
        only (no tmRNAs).
        :param c: treat sequence as circular
        :type c: bool
        :param d: search both strands of DNA
        :type d: bool
        :param m: search for tmRNAs
        :type m: bool
        :param t: search for tRNAs
        :type t: bool
        :return:
        """
        command = f"aragorn -gcbact -br -wa "
        if c is False:
            command += "-l "
        if d is False:
            command += "-s "
        if m is True:
            command += "-m "
        if t is True:
            command += "-t "
        command += f"-o {self.output} {self.input}"

        Popen(shlex.split(command)).wait()

    def read_output(self):
        """
        Reads the Aragorn output file and joins the lines into a single
        string which it populates into `out_str`.
        :return:
        """
        with open(self.output, "r") as fh:
            self.out_str = "".join(fh.readlines())

    def parse_tmrnas(self):
        """
        Searches `out_str` for matches to a regular expression for
        Aragorn tmRNAs.
        :return:
        """
        # values:        orient,  start, stop,            peptide tag
        # indices:          0       1     2                    3
        re_str = "tmRNA\s+(\w+)?\[(-?\d+),(\d+)\]\s+\d+,\d+\s+([\w|*]*)"
        regex = re.compile(re_str, re.MULTILINE | re.DOTALL)

        tmrnas = regex.findall(self.out_str)

        # Iterate through tmRNAs and create dictionaries for each
        for tmrna in tmrnas:
            tmrna_data = dict()

            # Increment tally
            self.tmrna_tally += 1

            # If no complement character found, forward orientation
            if tmrna[0] == "":
                orientation = "forward"
            else:
                orientation = "reverse"

            # Coordinates need to be converted from 1-based closed to
            # 0-based half-open.
            start, stop = basic.reformat_coordinates(
                int(tmrna[1]), int(tmrna[2]), "1_closed", "0_half_open")

            tmrna_data["Orientation"] = orientation
            tmrna_data["Start"] = start
            tmrna_data["Stop"] = stop
            tmrna_data["PeptideTag"] = tmrna[3]

            self.tmrnas.append(tmrna_data)

    def parse_trnas(self):
        """
        Calls two helper methods to parse determinate and indeterminate
        tRNAs.
        :return:
        """
        self.parse_determinate_trnas()
        self.parse_indeterminate_trnas()

    def parse_determinate_trnas(self):
        """
        Searches `out_str` for matches to a regular expression for
        Aragorn tRNAs of determinate isotype.
        :return:
        """
        # values:       aa., orient, start, stop,          anticodon,
        # indices:       0      1      2     3                 4
        re_str = "tRNA-(\w+)\s+(c)?\[(-?\d+),(\d+)\]\s+\d+\s+\((\w+)\)\s+" \
                 "([.]*?\w+[.]*?)\s+([( )dAtv.]*)"
        # values:     sequence,      structure
        # indices:       5               6
        regex = re.compile(re_str, re.MULTILINE | re.DOTALL)

        trnas = regex.findall(self.out_str)

        # Iterate through tRNAs and create dictionaries for each
        for trna in trnas:
            trna_data = dict()

            # Increment tally
            self.trna_tally += 1

            amino_acid = trna[0]

            # If no complement character found, forward orientation
            if trna[1] == "":
                orientation = "forward"
            else:
                orientation = "reverse"

            # If start coordinate is negative, add 1 because Aragorn indexing
            # goes from -1 to 1 (skips 0)
            start, stop = int(trna[2]), int(trna[3])
            if start < 0:
                start += 1

            # Coordinates need to be converted from 1-based closed to
            # 0-based half open
            start, stop = basic.reformat_coordinates(
                start, stop, "1_closed", "0_half_open")

            anticodon = trna[4]
            sequence = trna[5]
            structure = trna[6]

            # Convert structure to true dot-bracket notation
            for char in [" ", "d", "A", "t", "v"]:
                structure = structure.replace(char, ".")

            while len(structure) < len(sequence):
                structure += "."

            trna_data["Orientation"] = orientation
            trna_data["Start"] = start
            trna_data["Stop"] = stop
            trna_data["AminoAcid"] = amino_acid
            trna_data["Anticodon"] = anticodon.lower()
            trna_data["Sequence"] = sequence
            trna_data["Structure"] = structure

            self.trnas.append(trna_data)

    def parse_indeterminate_trnas(self):
        """
        Searches `out_str` for matches to a regular expression for
        Aragorn tRNAs of indeterminate isotype.
        :return:
        """
        # values: possible amino acids., orient, start, stop,
        # indices:          0      1        2      3     4
        re_str = "tRNA-\?\((\w+)\|(\w+)\)\s+(c)?\[(-?\d+),(\d+)\]\s+\d+\s+" \
                 "\((\w+)\)\s+([.]*?\w+[.]*?)\s+([( )dAtv.]*)"
        # values: anticodon,      sequence,       structure
        # indices:    5              6                7
        regex = re.compile(re_str, re.MULTILINE | re.DOTALL)

        trnas = regex.findall(self.out_str)

        # Iterate through tRNAs and create dictionaries for each
        for trna in trnas:
            trna_data = dict()

            # Increment tally
            self.trna_tally += 1

            amino_acid = trna[0] + "|" + trna[1]

            # If no complement character found, forward orientation
            if trna[2] == "":
                orientation = "forward"
            else:
                orientation = "reverse"

            # If start coordinate is negative, add 1 because Aragorn indexing
            # goes from -1 to 1 (skips 0)
            start, stop = int(trna[3]), int(trna[4])
            if start < 0:
                start += 1

            # Coordinates need to be converted from 1-based closed to
            # 0-based half open
            start, stop = basic.reformat_coordinates(
                start, stop, "1_closed", "0_half_open")

            anticodon = trna[5]
            sequence = trna[6]
            structure = trna[7]

            # Convert structure to true dot-bracket notation
            for char in [" ", "d", "A", "t", "v"]:
                structure = structure.replace(char, ".")

            while len(structure) < len(sequence):
                structure += "."

            trna_data["Orientation"] = orientation
            trna_data["Start"] = start
            trna_data["Stop"] = stop
            trna_data["AminoAcid"] = amino_acid
            trna_data["Anticodon"] = anticodon.lower()
            trna_data["Sequence"] = sequence
            trna_data["Structure"] = structure

            self.trnas.append(trna_data)
