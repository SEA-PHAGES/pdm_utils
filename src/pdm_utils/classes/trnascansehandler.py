import os
import shlex
from subprocess import Popen
import re

from pdm_utils.functions import basic


class TRNAscanSEHandler:
    def __init__(self, identifier, sequence):
        self.id = identifier
        self.sequence = sequence

        # I/O attributes
        self.temp_dir = "/tmp/trnascanse"
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

    def run_trnascanse(self, x=10):
        """
        Set up tRNAscan-SE command, then run it. Explanation of
        arguments:
        :param x: score cutoff for tRNAscan-SE
        :type x: int
        :return:
        """
        command = f"tRNAscan-SE -B -H -qQ --detail -X {x} -o /dev/null "
        command += f"-f {self.output} {self.input}"
        Popen(shlex.split(command)).wait()

    def read_output(self):
        """
        Reads the Aragorn output file and joins the lines into a single
        string which it populates into `out_str`.
        :return:
        """
        with open(self.output, "r") as fh:
            self.out_str = "".join(fh.readlines())

    def parse_trnas(self):
        """
        Searches `out_str` for matches to a regular expression for
        tRNAscan-SE tRNAs.
        :return:
        """
        # values:   start, stop,           length,           amino acid,
        # indices:    0     1                 2                   3
        re_str = "\((\d+)-(\d+)\)\s+Length: (\d+)\s+\w+\s+Type: (\w+)\s+" \
                 "Anticodon: (\w+).*?\s+Seq: (\w+)\s+Str: ([>.<]*)"
        # values:           anticodon,      sequence,     structure
        # indices:             4               5             6
        regex = re.compile(re_str, re.MULTILINE | re.DOTALL)

        trnas = regex.findall(self.out_str)

        # Iterate through tRNAs and create dictionaries for each
        for trna in trnas:
            trna_data = dict()

            # Increment tally
            self.trna_tally += 1

            start, stop = int(trna[0]), int(trna[1])
            # If stop > start, forward orientation
            if start < stop:
                orientation = "forward"
            else:
                orientation = "reverse"
                start, stop = stop, start

            # Coordinates need to be converted from 1-based closed to
            # 0-based half open
            start, stop = basic.reformat_coordinates(
                start, stop, "1_closed", "0_half_open")

            amino_acid = trna[3]
            anticodon = trna[4]
            sequence = trna[5]

            # Convert structure to true dot-bracket notation
            structure = trna[6].replace(">", "(").replace("<", ")")
            # Pad structure with periods to ensure full coverage of sequence
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
