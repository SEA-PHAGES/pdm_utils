"""Represents a collection of data about a tmRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import re
from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC

from pdm_utils.classes import evaluation
from pdm_utils.classes.aragornhandler import AragornHandler
from pdm_utils.functions import basic

# Extracts peptide tag from note field acid and anticodon from note field for Aragorn-determinate
# or tRNAscan-SE- determinate or indeterminate tRNAs
NOTE_STANDARD_REGEX = re.compile("Tag peptide:\s+([\w|*]*)")


class Tmrna:
    def __init__(self):
        """
        Constructor method for a tmRNA object.
        """
        # The only attribute a tmRNA really needs to have is a sequence
        self.seq = Seq("", IUPAC.ambiguous_dna)
        self.length = 0

        # Information about the tmRNA with respect to its parent genome
        self.genome_id = ""     # Identifier for the parent genome
        self.genome_length = -1
        self.start = -1         # Start coord in parent genome, 0-indexed
        self.stop = -1          # Stop coord in parent genome, 0-indexed
        self.coordinate_format = ""
        self.orientation = ""   # "F" or "R", etc. - relative to parent genome

        # MySQL tmRNAs will also need these:
        self.id = ""            # Identifier for this gene
        self.name = ""          # tRNA number in the parent genome
        self.peptide_tag = ""   # Degradation tag

        # Genbank-formatted flat files should also have these attributes:
        self.seqfeature = None  # BioPython SeqFeature object for this tRNA
        self.locus_tag = ""     # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = ""
        self._gene_num = ""
        self.parts = 0          # Number of regions this gene encompasses
        self.gene = ""          # Gene number parsed from feature
        self.note = ""          # Raw note field

        # Aragorn data
        self.aragorn_run = False
        self.aragorn_data = None

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

    def set_gene(self, value, delimiter=None, prefix_set=None):
        """Set the gene attribute.

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
                # raise ValueError(f"tRNA.seq cannot use type '{type(seq)}'")
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
    def set_seqfeature(self):
        """
        Create a SeqFeature object with which to populate the
        `seqfeature` attribute.
        :return:
        """
        # SeqFeature coordinates are 0-based half-open
        start, stop = basic.reformat_coordinates(
            self.start, self.stop, self.coordinate_format, "0_half_open")

        # SeqFeature orientation is (-1, 1) instead of ("R", "F")
        strand = basic.reformat_strand(self.orientation, "numeric")

        # Standard genes will have start < stop
        if self.start <= self.stop:
            self.seqfeature = SeqFeature(FeatureLocation(start, stop),
                                         strand=strand, type=self.type)
        # Wrap-around genes will have stop < start
        else:
            self.seqfeature = SeqFeature(CompoundLocation(
                [FeatureLocation(start, self.genome_length),
                 FeatureLocation(0, stop)]), strand=strand, type=self.type)
        # Add feature qualifiers
        self.seqfeature.qualifiers = self.get_qualifiers()

    def get_qualifiers(self):
        """
        Helper function that uses tRNA data to populate the qualifiers
        attribute of `seqfeature`.
        :return: qualifiers OrderedDict()
        """
        qualifiers = OrderedDict()
        qualifiers["gene"] = [self.name]
        if self.locus_tag != "":
            qualifiers["locus_tag"] = [self.locus_tag]
        if self.peptide_tag != "":
            qualifiers["note"] = [f"Tag peptide: {self.peptide_tag}"]

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
            ah = AragornHandler(self.id, str(self.seq))
            ah.write_fasta()
            ah.run_aragorn(m=True, t=False)  # search (linear) self.seq for tmRNAs
            ah.read_output()
            ah.parse_tmrnas()
            if ah.tmrna_tally == 1:
                self.aragorn_data = ah.tmrnas[0]
            # else:
                # print(f"Aragorn found {ah.tmrna_tally} tmRNAs in this region.")
            self.aragorn_run = True
        else:
            print("Cannot run Aragorn on 0-length sequence.")

    def parse_peptide_tag(self):
        """
        Parse the `peptide_tag` attribute out of the note field.
        :return:
        """
        # For SEA annotations, the note should contain the peptide tag
        results = NOTE_STANDARD_REGEX.findall(self.note)
        if len(results) == 1:
            self.peptide_tag = results[0]

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
                     f"{self.locus_tag}({self.id})."
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
        # Make sure self.orientation is "F"/"R"
        orientation = basic.reformat_strand(self.orientation, fmt, case)

        # If Aragorn predicts a tmRNA, get orientation in "F"/"R" format
        if self.aragorn_data is not None:
            aragorn_orient = basic.reformat_strand(
                self.aragorn_data["Orientation"], fmt, case)
        else:
            aragorn_orient = None

        result = f"The annotated orientation '{orientation}' "

        # If only Aragorn predicts a tRNA
        if aragorn_orient is not None:
            # If Aragorn predicts "F" it agrees with annotated orientation
            if aragorn_orient == "F":
                result += "matches the Aragorn prediction."
                status = success
            else:
                result += "is backwards relative to the Aragorn prediction."
                status = fail
        # If Aragorn doesn't predict a tmRNA, we presume annotation is wrong.
        else:
            result += "is at odds with expectations (NO tmRNA HERE)."
            status = fail

        definition = f"Check whether the annotated orientation for {self.locus_tag} ({self.id})" \
                     f" matches the Aragorn prediction."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class - I think this is exactly the same as
    #  in Cds
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

        result = f"The gene qualifier is {self.gene}. It "
        if isinstance(value, int):
            result += "contains an integer, as expected."
            status = success
        else:
            result += "does not contain an integer, which is not expected."
            status = fail

        definition = f"Check if the gene qualifier contains an integer for " \
                     f"{self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class
    def check_compatible_gene_and_locus_tag(self, eval_id=None,
                                            success="correct", fail="error",
                                            eval_def=None):
        """
        Check that gene and locus_tag attributes contain identical numbers
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

    def check_peptide_tag_valid(self, eval_id=None, success="correct",
                                fail="error", eval_def=None):
        """
        Checks whether the annotated peptide tag contains any letters
        not strictly within the protein alphabet.
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
        result = f"The annotated peptide tag '{self.peptide_tag}' is "
        if self.peptide_tag != "" and self.peptide_tag[-1] == "*":
            letters = set(self.peptide_tag[:-1])
            bad_letters = letters.difference(IUPAC.IUPACProtein.letters)
            if not bad_letters:
                result += "formatted correctly, and uses the correct " \
                          "alphabet."
                status = success
            else:
                result += f"formatted correctly, but has " \
                          f"{len(bad_letters)} characters that fall outside " \
                          f"the protein alphabet."
                status = fail
        elif self.peptide_tag != "" and self.peptide_tag[-1] != "*":
            letters = set(self.peptide_tag)
            bad_letters = letters.difference(IUPAC.IUPACProtein.letters)
            if not bad_letters:
                result += "using the correct alphabet, but doesn't end in '*'."
                status = fail
            else:
                result += f"using {len(bad_letters)} characters outside the " \
                          f"protein alphabet, and doesn't end in '*'."
                status = fail
        else:
            result += "does not exist, so could not be checked."
            status = "unchecked"

        definition = f"Check that the peptide tag appears to be structured " \
                     f"correctly and fits the protein alphabet for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_peptide_tag_correct(self, eval_id=None, success="correct",
                                  fail="error", eval_def=None):
        """
        Checks whether the annotated peptide tag matches the Aragorn
        output.
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
        result = f"The annotated peptide tag '{self.peptide_tag}' "

        if self.aragorn_run and self.aragorn_data is not None:
            aragorn_tag = self.aragorn_data["PeptideTag"]

            if self.peptide_tag == aragorn_tag:
                result += "matches the Aragorn output."
                status = success
            else:
                result += f"does not match the Aragorn output ({aragorn_tag})."
                status = fail
        elif self.aragorn_run and self.aragorn_data is None:
            result += "does not match the Aragorn output (no tmRNA)."
            status = fail
        else:
            result += "could not be checked, because Aragorn wasn't run."
            status = "unchecked"

        definition = f"Check that the annotated peptide tag matches the " \
                     f"Aragorn output for {self.locus_tag} {(self.id)}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_parts(self, eval_id=None, success="correct", fail="error",
                    eval_def=None):
        """
        Makes sure only one region exists for this tRNA.
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
        result = f"The tRNA has '{self.parts}' parts. This is "
        if self.parts == 1:
            result += "expected."
            status = success
        else:
            result += "unexpected."
            status = fail

        definition = f"Make sure there is only one region for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
