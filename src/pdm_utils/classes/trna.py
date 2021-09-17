"""Represents a collection of data about a tRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import re
from collections import OrderedDict

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.classes import evaluation
from pdm_utils.classes.aragornhandler import AragornHandler
from pdm_utils.classes.trnascansehandler import TRNAscanSEHandler
from pdm_utils.functions import basic

# Amino acids that we allow in the database
MYSQL_AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "fMet", "Gln", "Glu",
                     "Gly", "His", "Ile", "Ile2", "Leu", "Lys", "Met", "Phe",
                     "Pro", "Pyl", "SeC", "Ser", "Stop", "Thr", "Trp", "Tyr",
                     "Val", "OTHER"}

# Amino acids that Genbank allows in the product field
GENBANK_AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
                       "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
                       "Thr", "Trp", "Tyr", "Val", "OTHER"}

# Additional amino acids that SEA-PHAGES may annotate in the note field
SPECIAL_AMINO_ACIDS = {"fMet": "fMet", "Ile2": "Ile2", "Pyl": "Pyl",
                       "SeC": "SeC", "Sec": "SeC", "Term": "Stop",
                       "Stop": "Stop", "Sup": "Stop", "Undet": "OTHER"}

# Extracts amino acid from product field
PRODUCT_REGEX = re.compile("tRNA-(\w+)")

# Extracts amino acid and anticodon from note field for Aragorn-determinate
# or tRNAscan-SE- determinate or indeterminate tRNAs
NOTE_STANDARD_REGEX = re.compile("tRNA-(\w+)\s?\((\w+)\)")

# Extracts amino acid possibilities and anticodon from note field for
# Aragorn-indeterminate tRNAs
NOTE_SPECIAL_REGEX = re.compile("tRNA-\?\((\w+)\|(\w+)\)\s?\((\w+)\)")


class Trna:
    def __init__(self):
        """
        Constructor method for a tRNA object.
        """
        # The only attribute a tRNA really needs to have is a sequence
        self.seq = Seq("", IUPAC.ambiguous_dna)
        self.length = 0

        # Information about the tRNA with respect to its parent genome
        self.genome_id = ""     # Identifier for the parent genome
        self.genome_length = -1
        self.start = -1         # Start coord in parent genome, 0-indexed
        self.stop = -1          # Stop coord in parent genome, 0-indexed
        self.coordinate_format = ""
        self.orientation = ""   # "F" or "R", etc. - relative to parent genome

        # MySQL tRNAs will also need these:
        self.id = ""            # Identifier for this gene
        self.name = ""          # tRNA number in the parent genome
        self.amino_acid = "OTHER"   # Which aa. is this tRNA charged with?
        self.anticodon = "nnn"  # Which anticodon is at bottom of the C-loop?
        self.structure = ""     # What is the predicted secondary structure?

        # Genbank-formatted flat files should also have these attributes:
        self.seqfeature = None  # BioPython SeqFeature object for this tRNA
        self.locus_tag = ""     # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = ""
        self._gene_num = ""
        self.parts = 0          # How many regions does the tRNA have
        self.gene = ""          # Gene number parsed from feature
        self.product = ""       # Raw product field
        self.note = ""          # Raw note field

        # Aragorn data
        self.aragorn_data = None

        # tRNAscan-SE data
        self.trnascanse_data = None

        # Which program(s) support the annotation?
        self.sources = set()

        # Which program to check for valid anticodon?
        self.use = None

        # Useful for processing data from various sources:
        self.evaluations = list()
        self.type = ""
        self._start_stop_orient_id = tuple()
        self._end_orient_id = tuple()
        self._start_end_id = tuple()

    # TODO: create base feature class
    def set_locus_tag(self, tag="", delimiter="_", check_value=None):
        """
        Populate the locus_tag and parse the locus_tag number.
        :param tag: Input locus_tag data
        :type tag: str
        :param delimiter: Value used to split locus_tag data
        :type delimiter: str
        :param check_value: Genome name or other value that will be
        used to parse the locus_tag to identify the feature number
        :type check_value: str
        """
        # List of delimiters commonly found in locus tags. "_" is most common
        delims = ["_", "-", "."]

        # If tag is None, locus_tag and _locus_tag_num remain default ""
        if tag is None:
            return

        # Else, we'll set the tag directly, then decide how to parse number
        self.locus_tag = tag

        # If delimiter is None, choose based on what appears in locus_tag
        if delimiter is None:
            delimiter = basic.choose_most_common(tag, delims)

        # If no check_value was given, use self.genome_id
        if check_value is None:
            check_value = self.genome_id

        # Compile a search pattern and split the incoming tag on delimiter
        pattern = re.compile(check_value.lower())
        parts = tag.split(delimiter)

        # Attempt to parse locus tag number from the locus tag
        parsed_value = None
        for i in range(len(parts)):
            # Search each part of the locus tag for check_value
            search_result = pattern.search(parts[i].lower())
            if search_result is not None:
                if i == len(parts) - 1:
                    # If the check_value was found in last part of locus tag,
                    # value is probably everything after check_value
                    parsed_value = parts[i][search_result.end():]
                else:
                    # If the check_value was found anywhere else in locus tag,
                    # value is probably found in the next index
                    parsed_value = parts[-1]
                # No need to keep searching if we found a match
                break

        if parsed_value is None:
            # Last-ditch effort to find locus tag number
            parsed_value = parts[-1]

        # Remove generic 'gp' or 'orf' prefixes
        if parsed_value.lower().startswith("gp"):
            parsed_value = parsed_value[2:]
        elif parsed_value.lower().startswith("orf"):
            parsed_value = parsed_value[3:]
        else:
            pass

        # If the value we parsed is a digit, populate _locus_tag_num
        left, right = basic.split_string(parsed_value)
        if right != "":
            self._locus_tag_num = right
        else:
            self._locus_tag_num = parsed_value

    # TODO: create base feature class
    def set_name(self, value=None):
        """
        Set the feature name. Ideally, the name of the CDS will be an
        integer. This information can be stored in multiple fields in
        the GenBank-formatted flat file. The name is derived from one
        of several qualifiers.
        :param value: Indicates a value that should be used to
        directly set the name regardless of the 'gene' and
        '_locus_tag_num' attributes.
        :type value: str
        """
        # The CDS feature number may be stored in several places in the
        # gene, product, note, function, and locus_tag qualifiers of a
        # flat file, but it is unpredictable. Not all of these qualifiers
        # are always present.
        # 1. PECAAN Draft and new SEA-PHAGES Final:
        #    The 'gene' qualifier should be present and contain an integer.
        # 2. SEA-PHAGES Final in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.
        # 3. Non-SEA-PHAGES in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.
        if value is None:
            name = ""
            list1 = ["_locus_tag_num", "_gene_num"]

            # First see if any num attributes have a float.
            x = 0
            while name == "" and x < len(list1):
                value = getattr(self, list1[x])
                if basic.is_float(value):
                    name = value
                x += 1

            # Second see if any num attributes have non-empty values.
            # At this point, it's very unpredictable. Values could be like
            # 'terL' in the gene or like '10a' in the locus_tag.
            if name == "":
                list2 = ["gene"]
                list2.extend(list1)
                y = 0
                while name == "" and y < len(list2):
                    value = getattr(self, list2[y])
                    if value != "":
                        name = value
                    y += 1
        self.name = value

    # TODO: create base feature class
    def set_gene(self, value, delimiter=None, prefix_set=None):
        """
        Set the gene attribute.
        :param value: Gene data to parse. Also passed to set_num().
        :type value: str
        :param delimiter: Passed to set_num().
        :type delimiter: str
        :param prefix_set: Passed to set_num().
        :type prefix_set: set
        """
        value = value.strip()
        self.gene = value
        self._gene_num = value
        self.set_num("_gene_num", value,
                     delimiter=delimiter, prefix_set=prefix_set)

    # TODO: create base feature class
    def set_num(self, attr, description, delimiter=None, prefix_set=None):
        """
        Set a number attribute from a description.
        :param attr: Attribute to set the number.
        :type attr: str
        :param description: Description data from which to parse the number.
        :type description: str
        :param delimiter: Value used to split the description data.
        :type delimiter: str
        :param prefix_set: Valid possible delimiters in the description.
        :type prefix_set: set
        """
        # Used to set product_num, function_num, note_num, gene_num
        # Sometimes feature number is stored within the description
        # (e.g. 'gp10; terminase' or 'terminase; gp10')
        # List of delimiters found in description fields
        # ";" is probably the most common.
        delimiters = [";", ","]
        if delimiter is None:
            delimiter = basic.choose_most_common(description, delimiters)
        if prefix_set is None:
            prefix_set = {"gp", "orf", ""}

        # Iterate through the list of strings, and select the first
        # string that looks like a gene number.
        split_value = description.split(delimiter)
        num = ""
        x = 0
        while (num == "" and x < len(split_value)):
            string = split_value[x]
            string = string.strip()
            left, right = basic.split_string(string)
            # Possible returns:
            # 1. left is alpha and right is number or float (gp10 = gp, 10),
            # 2. string is alpha (terminase = terminase, ""),
            # 3. string is number (10 = "", 10).
            if left.lower() in prefix_set:
                num = right
            x += 1
        setattr(self, attr, num)

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's more Pythonic.
    def get_begin_end(self):
        """
        Accesses feature coordinates in transcription begin-end format.
        :return: (begin, end)
        """
        # Get a copy of the orientation in fr_long format:
        orientation = basic.reformat_strand(self.orientation, "fr_long")

        if orientation == "forward":
            # Rightward transcribed gene
            begin, end = self.start, self.stop
        elif orientation == "reverse":
            # Leftward transcribed gene
            begin, end = self.stop, self.start
        else:
            # Unexpected orientation
            begin = end = -1

        return begin, end

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's documented more cleanly and avoids using
    #  keyword 'format'.
    def set_orientation(self, value, fmt, case=False):
        """
        Set the orientation based on the indicated format.
        :param value: orientation value
        :type value: int or str
        :param fmt: how orientation should be formatted
        :type fmt: str
        :param case: whether to capitalize the first letter of
        orientation
        :type case: bool
        :return:
        """
        self.orientation = basic.reformat_strand(value, fmt, case)

    # TODO: create base feature class
    def set_location_id(self):
        """
        Create identifier tuples containing feature location data. For
        this method we only care about gene boundaries and will ignore
        any multi-part elements to the gene.
        :return:
        """
        # Tuple with ALL relevant information
        self._start_stop_orient_id = (self.start, self.stop, self.orientation)

        begin, end = self.get_begin_end()

        # If only end coordinate is used, we need to know what orientation the
        # gene is in - two genes could theoretically be transcribed in opposed
        # directions and have the same stop position
        self._end_orient_id = (end, self.orientation)

        # Or we can have transcription begin/end coordinates, in which case
        # the direction can be gleaned from these values
        self._start_end_id = (begin, end)

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it has slightly improved logic.
    def reformat_start_and_stop(self, fmt):
        """
        Convert existing start and stop coordinates to the indicated
        new format; also updates the coordinate format attribute to
        reflect any change.
        :param fmt: the new desired coordinate format
        :type fmt: str
        :return:
        """
        # Only need to do work if the new coordinate format differs from
        # the existing coordinate format
        if fmt != self.coordinate_format:
            self.start, self.stop = basic.reformat_coordinates(
                self.start, self.stop, self.coordinate_format, fmt)
            self.coordinate_format = fmt

    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented differently.
    def set_nucleotide_length(self, use_seq=False):
        """
        Set the nucleotide length of this gene feature.
        :param use_seq: whether to use the Seq feature to calculate
        nucleotide length of this feature
        :type use_seq: bool
        :return:
        """
        if use_seq:
            self.length = len(self.seq)
        else:
            if self.coordinate_format == "0_half_open":
                # Python coordinate format
                self.length = self.stop - self.start
            elif self.coordinate_format == "1_closed":
                # Genbank coordinate format
                self.length = self.stop - self.start + 1
            else:
                # Unknown coordinate format
                self.length = -1

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because the logic is slightly better and it raises an
    #  Error if something is wrong with given seq.
    def set_nucleotide_sequence(self, value=None, parent_genome_seq=None):
        """
        Set this feature's nucleotide sequence
        :param value: sequence
        :type value: str or Seq
        :param parent_genome_seq: parent genome sequence
        :type parent_genome_seq: Seq
        :raise: ValueError
        :return:
        """
        # If seq is given we'll try to use that
        if value is not None:
            # If seq is a BioPython Seq object, use it directly
            if isinstance(value, Seq):
                self.seq = value.upper()
            # If seq is a string, we can coerce it to a BioPython Seq object
            elif isinstance(value, str):
                self.seq = Seq(value.upper(), IUPAC.ambiguous_dna)
            # If seq is something else entirely, we cannot use it
            else:
                pass
                # raise ValueError(f"tRNA.seq cannot use type '{type(value)}'")
        # If instead a parent_genome_seq is given, we'll try to use that
        elif parent_genome_seq is not None and self.seqfeature is not None:
            try:
                self.seq = self.seqfeature.extract(parent_genome_seq)
            # TODO: tighten exception clause, OR let the error pass to a higher
            #  level in the call stack
            except:
                # Leave as default
                pass
        # If neither is given, leave self.seq as default
        else:
            pass

    # TODO: create base feature class - fully equivalent to the version in Cds,
    #  but better documented.
    def set_seqfeature(self, type=None):
        """
        Create a SeqFeature object with which to populate the
        `seqfeature` attribute.
        :return:
        """
        if type is None:
            type = self.type

        # SeqFeature coordinates are 0-based half-open
        start, stop = basic.reformat_coordinates(
            self.start, self.stop, self.coordinate_format, "0_half_open")

        # SeqFeature orientation is (-1, 1) instead of ("R", "F")
        strand = basic.reformat_strand(self.orientation, "numeric")

        self.seqfeature = self.create_seqfeature(type, start, stop, strand)

    def create_seqfeature(self, type, start, stop, strand):
        # Standard genes will have start < stop
        if start <= stop:
            seqfeature = SeqFeature(FeatureLocation(start, stop),
                                         strand=strand, type=type)
        # Wrap-around genes will have stop < start
        else:
            seqfeature = SeqFeature(CompoundLocation(
                [FeatureLocation(start, self.genome_length),
                 FeatureLocation(0, stop)]), strand=strand, type=type)
        # Add feature qualifiers
        seqfeature.qualifiers = self.get_qualifiers(type)

        return seqfeature

    def get_qualifiers(self, type):
        """
        Helper function that uses tRNA data to populate the qualifiers
        attribute of `seqfeature`.
        :return: qualifiers OrderedDict()
        """
        qualifiers = OrderedDict()
        if type == "tRNA":
            qualifiers["gene"] = [self.name]
            if self.locus_tag != "":
                qualifiers["locus_tag"] = [self.locus_tag]

            if self.amino_acid not in GENBANK_AMINO_ACIDS:
                amino_acid = "OTHER"
            else:
                amino_acid = self.amino_acid

            qualifiers["product"] = [f"tRNA-{amino_acid}"]
            qualifiers["note"] = [f"tRNA-{amino_acid} ({self.anticodon})"]
        if type == "gene":
            qualifiers["gene"] = [self.name]
            qualifiers["locus_tag"] = [self.locus_tag]

        return qualifiers

    def run_aragorn(self):
        """
        Uses an AragornHandler object to negotiate the flow of
        information between this object and Aragorn.
        :return:
        """
        if self.id is not None and self.id != "":
            identifier = self.id
        else:
            identifier = "aragorn"

        if len(self.seq) > 0:
            ah = AragornHandler(identifier, str(self.seq))
            ah.write_fasta()
            ah.run_aragorn()            # search (linear) self.seq for tRNAs
            ah.read_output()
            ah.parse_trnas()
            if ah.trna_tally == 1:
                self.aragorn_data = ah.trnas[0]
                self.sources.add("aragorn")
            # else:
                # print(f"Aragorn found {ah.trna_tally} tRNAs in this region.")
        else:
            print("Cannot run Aragorn on 0-length sequence.")

    def run_trnascanse(self):
        """
        Uses a TRNAscanSEHandler object to negotiate the flow of
        information between this object and tRNAscan-SE.
        :return:
        """
        if self.id is not None and self.id != "":
            identifier = self.id
        else:
            identifier = "trnascanse"

        if len(self.seq) > 0:
            th = TRNAscanSEHandler(identifier, str(self.seq))
            th.write_fasta()
            th.run_trnascanse()
            th.read_output()
            th.parse_trnas()
            if th.trna_tally == 1:
                self.trnascanse_data = th.trnas[0]
                self.sources.add("trnascan")
            # else:
                # print(f"tRNAscan-SE found {th.trna_tally} tRNAs in this region.")
        else:
            print("Cannot run tRNAscan-SE on 0-length sequence.")

    def set_amino_acid(self, value):
        """
        Sets the `amino_acid` attribute using the indicated value.
        :param value: the Amino acid to be used
        :type value: str
        :raise: ValueError
        :return:
        """
        self.amino_acid = value

    def parse_amino_acid(self):
        """
        Attempts to parse the `amino_acid` attribute from the `product`
        and `note` attributes.
        :return:
        """
        results = PRODUCT_REGEX.findall(self.product)
        # If we got a single item (e.g. the amino acid)
        if len(results) == 1:
            # Begin using it as the amino acid
            amino_acid = results[0]
            # If amino acid we got is "OTHER", we need to check the note
            # field for a more specific amino acid
            if amino_acid == "OTHER":
                results = NOTE_STANDARD_REGEX.findall(self.note)
                # If we got two items (e.g. the amino acid and anticodon)
                # and the amino acid is in our allowed special amino acids
                if len(results) != 0 and results[0][0] in SPECIAL_AMINO_ACIDS:
                    # Use the re-formatted special amino acid
                    amino_acid = SPECIAL_AMINO_ACIDS[results[0][0]]
                # Else, we will leave amino acid as "OTHER"
            # Else, there is nothing else to do but set the attribute
            self.amino_acid = amino_acid
        # Else, the product is not formatted properly and will fail checks

    def set_anticodon(self, value):
        """
        Sets the `anticodon` attribute using the indicated value.
        :param value: the anticodon to use for this tRNA
        :type value: str
        :return:
        """
        self.anticodon = value.lower()

    def parse_anticodon(self):
        """
        Attempts to parse the `anticodon` attribute from the `note`
        attribute.
        :return:
        """
        results = NOTE_STANDARD_REGEX.findall(self.note)
        # If we got two items (e.g. the amino acid and anticodon)
        if len(results) == 1:
            self.anticodon = results[0][-1].lower().replace('u', 't')
        else:
            results = NOTE_SPECIAL_REGEX.findall(self.note)
            # If we got three items (e.g. the two amino acid choices and
            # part of anticodon)
            if len(results) == 1:
                self.anticodon = results[0][-1].lower().replace('u', 't')
        # If the regex fails to parse an anticodon using either regex, it
        # will be left as "nnn"

    def set_structure(self, value):
        """
        Set the secondary structure string so downstream users can
        easily display the predicted fold of this tRNA.
        :param value: the string to use as the secondary structure
        :type value: str
        :return:
        """
        self.structure = value

    # Evaluations
    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented slightly differently.
    def set_eval(self, eval_id, definition, result, status):
        """
        Constructs and adds and Evaluation object to this feature's list of
        evaluations.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param definition: description of the evaluation
        :type definition: str
        :param result: description of the evaluation outcome
        :type result: str
        :param status: overall outcome of the evaluation
        :type status: str
        :return:
        """
        evl = evaluation.Evaluation(eval_id, definition, result, status)
        self.evaluations.append(evl)

    # TODO: create base feature class - nearly identical to Cds version,
    #  but exception clause is appropriately narrow.
    def check_attribute(self, attribute, check_set, expect=False, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """
        Checks whether the indicated feature attribute is present in
        the given check_set. Uses expect to determine whether the
        presence (or lack thereof) is an error, or correct.
        :param attribute: the gene feature attribute to evaluate
        :type attribute: str
        :param check_set: set of reverence values
        :type check_set: set
        :param expect: whether the attribute's value is expected to be
        in the reference set
        :type expect: bool
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        try:
            value = getattr(self, attribute)
        except AttributeError:
            value = None

        if value is not None:
            value_short = basic.truncate_value(str(value), 30, "...")
            result = f"The '{attribute}' value '{value_short}' is "

            # Returns a boolean
            outcome = basic.check_value_expected_in_set(
                value, check_set, expect)

            if outcome:
                result += "valid."
                status = success
            else:
                result += "not valid."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"

        definition = f"Check the value of the '{attribute}' attribute for " \
                     f"{self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class - Travis' version assumes value can't
    #  be None; if wrong, that method will fail because None can't be <=> x.
    def check_magnitude(self, attribute, expect, ref_value, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """
        Check that the magnitude of a numerical attribute meets
        expectations.
        :param attribute: the gene feature attribute to evaluate
        :type attribute: str
        :param expect: symbol designating direction of magnitude (>=<)
        :type expect: str
        :param ref_value: numerical value for comparison
        :type ref_value: int, float, datetime
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        try:
            value = getattr(self, attribute)
        except AttributeError:
            value = None

        if value is not None:
            result = f"The '{attribute}' value '{value}' is "
            if value > ref_value:
                compare = ">"
                result += "greater than "
            elif value == ref_value:
                compare = "="
                result += "equal to "
            else:
                compare = "<"
                result += "less than "
            result += f"{ref_value}, which is "
            if compare == expect:
                result += "expected."
                status = success
            else:
                result += "not expected."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"

        definition = f"Check the magnitude of the '{attribute}' attribute " \
                     f"for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_sources(self, eval_id=None, success="correct", fail="error",
                      eval_def=None):
        """
        Check that this tRNA's DNA sequence can successfully turn up
        a tRNA when run through Aragorn and tRNAscan-SE.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        self.run_aragorn()
        self.run_trnascanse()

        result = f"This tRNA gene's DNA sequence "

        if self.sources == set():
            # Neither program found a tRNA in this region
            self.use = None
            result += "does not appear to correspond to a tRNA."
            status = fail
        elif self.sources == {"aragorn"}:
            # Aragorn found a tRNA in this region but tRNAscan-SE did not
            self.use = "aragorn"
            result += "corresponds to a tRNA in Aragorn."
            status = success
        elif self.sources == {"trnascan"}:
            # tRNAscan-SE found a tRNA in this region but Aragorn did not
            self.use = "trnascan"
            result += "corresponds to a tRNA in tRNAscan-SE."
            status = success
        elif self.sources == {"aragorn", "trnascan"}:
            # Both programs found a tRNA in this region
            self.use = "both"
            result += "corresponds to a tRNA in both Aragorn and tRNAscan-SE."
            status = success
        else:
            # Something unexpected happened...
            self.use = None
            result += "gave rise to something unexpected."
            status = fail
            print(f"Invalid tRNA validation source in ({self.sources}).")

        definition = f"Check that there is bioinformatic support for the " \
                     f"existence of {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    # TODO: create base feature class - this version avoids overloading
    #  built-in format()
    def check_orientation(self, fmt="fr_short", case=True, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """
        Check that the orientation is set appropriately.
        :param fmt: indicates how coordinates should be formatted
        :type fmt: str
        :param case: indicates whether orientation data should be cased
        :type case: bool
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        expect = basic.reformat_strand(self.orientation, fmt, case)
        result = f"The orientation is {self.orientation}, and it is formatted "

        if self.orientation == expect:
            result += "correctly."
            status = success
        else:
            result += "incorrectly."
            status = fail

        definition = f"Check if the orientation is set correctly for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_orientation_correct(self, fmt="fr_short", case=True,
                                  eval_id=None, success="correct",
                                  fail="error", eval_def=None):
        """
        Check that the orientation agrees with the Aragorn and/or
        tRNAscan-SE predicted orientation. If Aragorn/tRNAscan-SE
        report a forward orientation, it means they agree with the
        annotated orientation. If they report reverse orientation,
        they think the annotation is backwards.
        :param fmt: indicates how coordinates should be formatted
        :type fmt: str
        :param case: indicates whether orientation data should be cased
        :type case: bool
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        orientation = basic.reformat_strand(self.orientation, fmt, case)
        result = f"The annotated orientation '{orientation}' "

        if self.use == "both":
            a_orient = basic.reformat_strand(
                self.aragorn_data["Orientation"], fmt, case)
            t_orient = basic.reformat_strand(
                self.trnascanse_data["Orientation"], fmt, case)
            if a_orient == t_orient == "F":
                result += "matches Aragorn and tRNAscan-SE predictions."
                status = success
            elif a_orient == t_orient == "R":
                result += "is backwards relative to Aragorn and tRNAscan-SE " \
                          "predictions."
                status = fail
            else:
                result += "cannot be sensibly checked, as Aragorn and " \
                          "tRNAscan-SE disagree."
                status = "unchecked"
        elif self.use == "aragorn":
            a_orient = basic.reformat_strand(
                self.aragorn_data["Orientation"], fmt, case)
            if a_orient == "F":
                result += "matches the Aragorn prediction."
                status = success
            elif a_orient == "R":
                result += "is backwards relative to Aragorn prediction."
                status = fail
            else:
                result += "cannot be sensibly checked, as Aragorn did" \
                          "something unexpected"
                status = "unchecked"
        elif self.use == "trnascan":
            t_orient = basic.reformat_strand(
                self.trnascanse_data["Orientation"], fmt, case)
            if t_orient == "F":
                result += "matches the tRNAscan-SE prediction."
                status = success
            elif t_orient == "R":
                result += "is backwards relative to tRNAscan-SE prediction."
                status = fail
            else:
                result += "cannot be sensibly checked, as tRNAscan-SE did" \
                          "something unexpected"
                status = "unchecked"
        # If neither predict a tRNA, we presume the annotation is wrong.
        else:
            result += "is at odds with expectations (no tRNA)."
            status = fail

        definition = f"Check whether the annotated orientation for {self.locus_tag} ({self.id})" \
                     f" matches the Aragorn/tRNAscan-SE prediction(s)."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_amino_acid_valid(self, eval_id=None, success="correct",
                               fail="error", eval_def=None):
        """
        Checks that the amino acid that has been annotated for this
        tRNA is in the set of amino acids that we have opted to allow
        in the MySQL database.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The annotated amino acid '{self.amino_acid}' is "

        if self.amino_acid in MYSQL_AMINO_ACIDS:
            result += "valid."
            status = success
        else:
            result += "invalid."
            status = fail

        definition = f"Check that the annotated amino acid is in the allowed" \
                     f" set for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_amino_acid_correct(self, eval_id=None, success="correct",
                                 fail="error", eval_def=None):
        """
        Checks that the amino acid that has been annotated for this
        tRNA agrees with the Aragorn and/or tRNAscan-SE prediction(s).
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The annotated isotype ({self.amino_acid}) "

        # If both programs found a tRNA here
        if self.use == "both":
            # All Aragorn isotypes are valid in our database except undet
            a_isotype = self.aragorn_data["AminoAcid"]
            if "|" in a_isotype:
                a_isotype = "OTHER"
            # Some tRNAscan-SE isotypes need to be converted to allowed value
            t_isotype = self.trnascanse_data["AminoAcid"]
            if t_isotype not in MYSQL_AMINO_ACIDS:
                if t_isotype in SPECIAL_AMINO_ACIDS:
                    t_isotype = SPECIAL_AMINO_ACIDS[t_isotype]
                else:
                    t_isotype = "OTHER"
            # if Aragorn and tRNAscan-SE are the same
            if a_isotype == t_isotype == self.amino_acid:
                result += f"is consistent with Aragorn ({a_isotype}) and " \
                          f"tRNAscan-SE ({t_isotype})."
                status = success
            # if Aragorn is Met and tRNAscan-se is consistent
            elif a_isotype == "Met" and t_isotype in ("Met", "fMet", "Ile2"):
                if self.amino_acid in (a_isotype, t_isotype):
                    result += f"is consistent with Aragorn ({a_isotype}) " \
                              f"and tRNAscan-SE ({t_isotype})."
                    status = success
                else:
                    result += f"is inconsistent with Aragorn ({a_isotype}) " \
                              f"and tRNAscan-SE ({t_isotype})."
                    status = fail
            # if tRNAscan-SE is Sup and Aragorn is consistent
            elif a_isotype in ("Pyl", "Stop") and t_isotype == "Sup":
                if self.amino_acid in (a_isotype, t_isotype):
                    result += f"is consistent with Aragorn ({a_isotype}) " \
                              f"and tRNAscan-SE ({t_isotype})."
                    status = success
                else:
                    result += f"is inconsistent with Aragorn ({a_isotype}) " \
                              f"and tRNAscan-SE ({t_isotype})."
                    status = fail
            else:
                result += f"cannot be sensibly checked because Aragorn " \
                          f"({a_isotype}) prediction differs from " \
                          f"tRNAscan-SE ({t_isotype})."
                status = "unchecked"
        # If only Aragorn found a tRNA here
        elif self.use == "aragorn":
            a_isotype = self.aragorn_data["AminoAcid"]
            if "|" in a_isotype:
                a_isotype = "OTHER"
            t_isotype = "no tRNA"
            if a_isotype == self.amino_acid:
                result += f"is consistent with Aragorn ({a_isotype}) and " \
                          f"tRNAscan-SE ({t_isotype})."
                status = success
            else:
                result += f"is inconsistent with Aragorn ({a_isotype}) " \
                          f"and tRNAscan-SE ({t_isotype})."
                status = fail
        # If only tRNAscan-SE found a tRNA here
        elif self.use == "trnascan":
            a_isotype = "no tRNA"
            t_isotype = self.trnascanse_data["AminoAcid"]
            if t_isotype not in MYSQL_AMINO_ACIDS:
                if t_isotype in SPECIAL_AMINO_ACIDS:
                    t_isotype = SPECIAL_AMINO_ACIDS[t_isotype]
                else:
                    t_isotype = "OTHER"
            if t_isotype == self.amino_acid:
                result += f"is consistent with Aragorn ({a_isotype}) and " \
                          f"tRNAscan-SE ({t_isotype})."
                status = success
            else:
                result += f"is inconsistent with Aragorn ({a_isotype}) " \
                          f"and tRNAscan-SE ({t_isotype})."
                status = fail
        # If neither program found a tRNA here
        else:
            result += "is inconsistent with Aragorn (no tRNA) and " \
                      "tRNAscan-SE (no tRNA)."
            status = fail

        definition = f"Check that the annotated amino acid is consistent " \
                     f"with the prediction(s) made by Aragorn and/or " \
                     f"tRNAscan-SE for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_anticodon_valid(self, eval_id=None, success="correct",
                              fail="error", eval_def=None):
        """
        Checks that the anticodon conforms to the expected length
        (2-4) and alphabet ("a", "c", "g", "t") or is "nnn".
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The anticodon '{self.anticodon}' is "
        expect = ["a", "c", "g", "t"]
        valid = [letter in expect for letter in self.anticodon]
        if len(self.anticodon) == 3:
            if sum(valid) == 3:
                result += "valid."
                status = success
            elif self.anticodon == "nnn":
                result += "valid."
                status = success
            else:
                result += "invalid."
                status = fail
        elif len(self.anticodon) in [2, 4]:
            if sum(valid) == len(self.anticodon):
                result += "valid."
                status = success
            else:
                result += "invalid."
                status = fail
        else:
            result += "invalid."
            status = fail

        definition = f"Check that the anticodon is valid for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_anticodon_correct(self, eval_id=None, success="correct",
                              fail="error", eval_def=None):
        """
        Checks that the annotated anticodon agrees with the prediction
        by Aragorn or tRNAscan-SE.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The anticodon '{self.anticodon}' "

        # If both programs found a tRNA here
        if self.use == "both":
            a_anti = self.aragorn_data["Anticodon"]
            t_anti = self.trnascanse_data["Anticodon"]
            if a_anti == t_anti == self.anticodon:
                result += f"is consistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = success
            elif a_anti == t_anti != self.anticodon:
                result += f"is inconsistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = fail
            else:
                result += f"cannot be sensibly checked because Aragorn " \
                          f"({a_anti}) prediction differs from tRNAscan-SE " \
                          f"({t_anti})."
                status = "unchecked"
        # If only Aragorn found a tRNA here
        elif self.use == "aragorn":
            a_anti = self.aragorn_data["Anticodon"]
            t_anti = "no tRNA"
            if a_anti == self.anticodon:
                result += f"is consistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = success
            else:
                result += f"is inconsistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = fail
        # If only tRNAscan-SE found a tRNA here
        elif self.use == "trnascan":
            a_anti = "no tRNA"
            t_anti = self.trnascanse_data["Anticodon"]
            if t_anti == self.anticodon:
                result += f"is consistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = success
            else:
                result += f"is inconsistent with Aragorn ({a_anti}) and " \
                          f"tRNAscan-SE ({t_anti})."
                status = fail
        # If neither program found a tRNA here
        else:
            result += "is inconsistent with Aragorn (no tRNA) and " \
                      "tRNAscan-SE (no tRNA)."
            status = fail

        definition = f"Check that the annotated anticodon agrees with the " \
                     f"Aragorn or tRNAscan-SE anticodon for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class
    def check_locus_tag_structure(self, check_value=None, only_typo=False,
                                  prefix_set=set(), case=True, eval_id=None,
                                  success="correct", fail="error", eval_def=None):
        """Check if the locus_tag is structured correctly.

        :param check_value:
            Indicates the genome id that is expected to be present.
            If None, the 'genome_id' parameter is used.
        :type check_value: str
        :param only_typo:
            Indicates if only the genome id spelling should be evaluated.
        :type only_typo: bool
        :param prefix_set:
            Indicates valid common prefixes, if a prefix is expected.
        :type prefix_set: set
        :param case:
            Indicates whether the locus_tag is expected to be capitalized.
        :type case: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        if check_value is None:
            check_value = self.genome_id

        results = []
        if only_typo:
            pattern = re.compile(check_value.lower())
            search_result = pattern.search(self.locus_tag.lower())
            if search_result == None:
                results.append("The genome ID is missing.")
        else:
            # Expected structure: SEA_TRIXIE_20
            parts = self.locus_tag.split("_")
            if case:
                if self.locus_tag != self.locus_tag.upper():
                    results.append("The capitalization is incorrect.")
            if len(parts) == 3:
                if prefix_set is not None:
                    if parts[0].upper() not in prefix_set:
                        results.append("The prefix is missing.")
                if parts[1].upper() != check_value.upper():
                    results.append("The genome ID is missing.")
                if not parts[2].isdigit():
                    results.append("The feature number is missing.")
            else:
                results.append("The locus_tag does not contain three parts.")
        result = f"The locus_tag qualifier is {self.locus_tag}. It is "
        if len(results) == 0:
            result = result + "structured correctly."
            status = success
        else:
            result = result + "not structured correctly. " + " ".join(results)
            status = fail
        definition = f"Check if the locus_tag qualifier is structured correctly for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class - fully equivalent to the version in Cds,
    #  but exception clause is appropriately narrow.
    def check_gene_structure(self, eval_id=None, success="correct",
                             fail="error", eval_def=None):
        """
        Check that the gene qualifier contains an integer.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        try:
            value = int(self.gene)
        except ValueError:
            value = self.gene

        result = f"The gene qualifier is {self.gene}. It"
        if isinstance(value, int):
            result += "contains an integer, as expected."
            status = success
        else:
            result += "does ont contain an integer, which is not expected."
            status = fail

        definition = f"Check if the gene qualifier contains an integer for " \
                     f"{self.locus_tag}({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class
    def check_compatible_gene_and_locus_tag(self, eval_id=None,
                                            success="correct", fail="error",
                                            eval_def=None):
        """
        Check that gene and locus_tag attributes contain identical
        numbers.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The numbers in the gene ({self.gene}) and " \
                 f"locus_tag ({self.locus_tag}) qualifiers are "
        if self.gene == self._locus_tag_num:
            result += "consistent."
            status = success
        else:
            result += "not consistent."
            status = fail

        definition = f"Check if the gene and locus_tag numbers are " \
                     f"consistent for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_length(self, eval_id=None, success="correct", fail="error",
                     eval_def=None):
        """
        Checks that the tRNA is in the expected range of lengths. The
        average tRNA gene is 70-90bp in length, but it is not uncommon
        to identify well-scoring tRNAs in the 60-100bp range.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The tRNA length ({self.length}) is "
        if self.length in range(60, 101):
            result += "in the expected range."
            status = success
        elif self.length < 60:
            result += "shorter than expected."
            status = fail
        else:
            result += "longer than expected."
            status = fail

        definition = f"Check if the tRNA length is in the expected range " \
                     f"for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_coordinates(self, eval_id=None, success="correct", fail="error",
                          eval_def=None):
        """

        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The tRNA coordinates ({self.start}, {self.stop}) "

        # We'll use Aragorn by default when available
        if self.use == "both" or self.use == "aragorn":
            start_offset = self.aragorn_data["Start"]
            aragorn_stop = self.aragorn_data["Stop"]
            stop_offset = self.length - aragorn_stop
            structure = self.aragorn_data["Structure"]
        # Use tRNAscan-SE if that's all we've got
        elif self.use == "trnascan":
            start_offset = self.trnascanse_data["Start"]
            trnascanse_stop = self.trnascanse_data["Stop"]
            stop_offset = self.length - trnascanse_stop
            structure = self.trnascanse_data["Structure"]
        else:
            start_offset = -99999     # arbitrary number that will never happen
            stop_offset = -99999      # arbitrary number that will never happen
            structure = ""

        # If tRNA is not forward-oriented, the start/stop offsets are reversed
        # because we put the real tRNA DNA sequence as the input, NOT the top
        # strand DNA sequence
        orientation = basic.reformat_strand(
            self.orientation, format="fr_short", case=True)
        if orientation == "R":
            start_offset, stop_offset = stop_offset, start_offset

        # Now check the offsets and adjust the structure as appropriate
        if start_offset == stop_offset == 0:
            result += "appear to be correct."
            status = success
        elif start_offset == stop_offset == -99999:
            result += f"should be investigated - neither Aragorn nor " \
                      f"tRNAscan-SE identify a tRNA in this position."
            status = fail
        else:
            result += "appear to be incorrect. Left coordinate "
            status = fail

            if start_offset == 0:
                result += "looks correct. "
            # Else we need to modify the structure to accommodate extra or
            # missing bases
            elif start_offset < 0:
                # tRNA is annotated shorter than Aragorn says it should be
                structure = structure[abs(start_offset):]
                result += f"should be moved left by {abs(start_offset)} " \
                          f"base(s). "
            else:
                # tRNA is annotated longer than Aragorn thinks it should be
                structure = ("." * start_offset) + structure
                result += f"should be moved right by {start_offset} base(s). "

            result += "Right coordinate "
            if stop_offset == 0:
                result += "looks correct."
            elif stop_offset < 0:
                # tRNA is annotated shorter than Aragorn says it should be
                structure = structure[:self.length]
                result += f"should be moved right by {abs(stop_offset)} " \
                          f"base(s)."
            else:
                # tRNA is annotated longer than Aragorn thinks it should be
                structure = structure + ("." * stop_offset)
                result += f"should be moved left by {stop_offset} base(s)."

        self.structure = structure

        definition = f"Check if the tRNA coordinates appear to match the " \
                     f"Aragorn or tRNAscan-SE prediction(s) for " \
                     f"{self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_terminal_nucleotides(self, eval_id=None, success="correct",
                                   fail="warning", eval_def=None):
        """
        Checks that the tRNA ends with "CCA" or "CC" or "C".
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        if len(self.seq) > 4:
            result = f"The tRNA ends with '...{self.seq[-4:]}'. "
            # tRNAs must end in "NCCA" to function - some are this way naturally
            # The others must end in "NCC" or "NC" or N
            if self.seq.endswith("CCA"):
                result += "This is expected."
                status = success
            elif self.seq.endswith("CC"):
                result += "This is expected."
                status = success
            elif self.seq.endswith("C"):
                result += "This is expected."
                status = success
            else:
                result += " This is not expected."
                status = fail
        else:
            result = f"Cannot check terminal nucleotides on a sequence of " \
                     f"length {len(self.seq)}."
            status = "unchecked"

        definition = f"Check that the correct terminal nucleotide(s) are " \
                     f"present for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_product_structure(self, eval_id=None, success="correct",
                                fail="error", eval_def=None):
        """
        Checks that the product field is formatted properly, and that
        the annotated amino acid is valid.
        Genbank enforces that all tRNA annotations that have a product
        field have it annotated as either 'tRNA-Xxx', where Xxx is one
        of the 20 standard amino acids, or 'tRNA-OTHER' for those tRNAs
        which decode a non-standard amino acid (e.g. SeC, Pyl, fMet).
        SEA-PHAGES may also append the anticodon parenthetically for a
        product field such as 'tRNA-Xxx(nnn)'.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The tRNA product is "
        if self.product != "":
            result += f"present ('{self.product}'). "

            # This regex will ignore the (nnn) anticodon if it's there; it
            # parses the amino acid prediction and makes sure it's valid.
            results = PRODUCT_REGEX.findall(self.product)

            # If the product field is formatted correctly, there should only
            # be a single result from the regex search
            if results is not None and len(results) != 1:
                result += f"Product is formatted properly. "
                status = success

            # If there is more or less than a single result, the product field
            # is not properly formatted
            else:
                result += "Product is not formatted properly."
                status = fail
        # If there is no product field, it's not formatted properly.
        else:
            result += "missing (''). "
            result += "Product is not formatted properly."
            status = fail

        definition = f"Check that the product is formatted properly for " \
                     f"{self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_note_structure(self, eval_id=None, success="correct",
                             fail="error", eval_def=None):
        """
        Checks that the note field is formatted properly.

        Genbank does not enforce any standard for the note field.
        This means that a note does not have to exist.

        SEA-PHAGES note fields should look like 'tRNA-Xxx(nnn)'.
        :param eval_id: unique identifier for the evaluation
        :type eval_id: str
        :param success: status if the outcome is successful
        :type success: str
        :param fail: status if the outcome is unsuccessful
        :type fail: str
        :param eval_def: description of the evaluation
        :type eval_def: str
        :return:
        """
        result = f"The tRNA note is "
        if self.note != "":
            result += f"present ('{self.note}'). "

            results1 = NOTE_STANDARD_REGEX.findall(self.note)
            results2 = NOTE_SPECIAL_REGEX.findall(self.note)

            if results1 is not None and len(results1) == 2:
                result += f"Note is formatted properly."
                status = success
            elif results2 is not None and len(results2) == 3:
                result += f"Note is formatted properly."
                status = success
            else:
                result += f"Note is not formatted properly."
                status = fail
        else:
            result += "missing (''). "
            result += "Note cannot be checked."
            status = "unchecked"

        definition = f"Check that the note is formatted properly for " \
                     f"{self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
