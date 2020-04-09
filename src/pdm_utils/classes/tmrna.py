"""Represents a collection of data about a tmRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import re
from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC

from pdm_utils.classes import eval
from pdm_utils.functions import basic


class TmrnaFeature:
    def __init__(self):
        """
        Constructor method for a tmRNA object.
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
        self.peptide_tag = ""   # Degradation tag

        # Aragorn data
        self.aragorn_data = None

        # Genbank-formatted flat files should also have these attributes:
        self.seqfeature = None  # BioPython SeqFeature object for this tRNA
        self.locus_tag = ""     # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = ""
        self.parts = 0          # Number of regions this gene encompasses
        self.gene = ""          # Gene number parsed from feature
        self.note = ""          # Raw note field

        # Useful for processing data from various sources:
        self.evaluations = list()
        self.type = "tmRNA"
        self._start_stop_orient_id = tuple()
        self._end_orient_id = tuple()
        self._start_end_id = tuple()

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's simpler and better documented.
    def set_locus_tag(self, tag="", delimiter="_", check_value=None):
        """
        Populate the locus tag attribute.
        :param tag: Input locus_tag data
        :type tag: str
        :param delimiter: Value used to split locus_tag data
        :type delimiter: str
        :param check_value: Genome name or other value that will be
        used to parse the locus_tag to identify the feature number
        :type check_value: str
        """
        # Set locus tag
        self.locus_tag = tag

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
                    parsed_value = parts[i + 1]
                # No need to keep searching if we found a match
                break

        if parsed_value is None:
            # Last-ditch effort to find locus tag number
            parsed_value = parts[-1]

        # If the value we parsed is a digit, populate _locus_tag_num
        if parsed_value.isdigit():
            self._locus_tag_num = parsed_value

    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented differently.
    def set_name(self, value=None):
        """
        Set the feature name.
        :param value: Value to use for `name` regardless of `gene`
        and `_locus_tag_num`
        :type value: str
        :return:
        """
        # If a value was given, we'll use this by default
        if value is not None:
            self.name = value
        # If no value was given we'll use the following hierarchy:
        # 1. self.gene
        # 2. self._locus_tag_num
        # 3. leave as default ""
        elif self.gene != "":
            self.name = self.gene
        elif self._locus_tag_num != "":
            self.name = self._locus_tag_num
        else:
            pass

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's documented more cleanly and avoids using
    #  keyword 'format'.
    def set_orientation(self, value, fmt, capitalize=False):
        """
        Set the orientation based on the indicated format.
        :param value: orientation value
        :type value: int or str
        :param fmt: how orientation should be formatted
        :type fmt: str
        :param capitalize: whether to capitalize the first letter of
        orientation
        :type capitalize: bool
        :return:
        """
        self.orientation = basic.reformat_strand(value, fmt, capitalize)

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

        return begin, end                       # tuple format is implicit

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

    def set_peptide_tag(self, value):
        """
        Set the `amino_acid` attribute.
        :param value: what to use as the amino acid
        :type value: str
        :raise: ValueError
        :return:
        """
        if isinstance(value, str):
            self.peptide_tag = value
        else:
            raise ValueError(f"Invalid type '{type(value)}' for peptide tag.")

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
            qualifiers["note"] = [f"Peptide tag: {self.peptide_tag}"]

        return qualifiers

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because the logic is slightly better and it raises an
    #  Error if something is wrong with given seq.
    def set_nucleotide_sequence(self, seq=None, parent_genome_seq=None):
        """
        Set this feature's nucleotide sequence
        :param seq: sequence
        :type seq: str or Seq
        :param parent_genome_seq: parent genome sequence
        :type parent_genome_seq: Seq
        :raise: ValueError
        :return:
        """
        # If seq is given we'll try to use that
        if seq is not None:
            # If seq is a BioPython Seq object, use it directly
            if isinstance(seq, Seq):
                self.seq = seq.upper()
            # If seq is a string, we can coerce it to a BioPython Seq object
            elif isinstance(seq, str):
                self.seq = Seq(seq, IUPAC.ambiguous_dna)
            # If seq is something else entirely, we cannot use it
            else:
                raise ValueError(f"tRNA.seq cannot use type '{type(seq)}'")
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

    # Evaluations
    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented slightly differently.
    def set_eval(self, eval_id, definition, result, status):
        """
        Constructs and adds and Eval object to this feature's list of
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
        evl = eval.Eval(eval_id, definition, result, status)
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
                     f"{self.id}."
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
                     f"for {self.id}."
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

        definition = f"Check if the orientation is set correctly for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO: create base feature class - I think this is exactly the same as
    #  in Cds
    def check_locus_tag_structure(self, check_value=None, only_typo=False,
                                  prefix_set=set(), case=True, eval_id=None,
                                  success="correct", fail="error",
                                  eval_def=None):
        """
        Check that the locus tag is structured correctly.
        :param check_value: expected genome id
        :type check_value: str
        :param only_typo: only genomeid spelling is evaluated
        :type only_typo: bool
        :param prefix_set: valid common prefixes, if one is expected
        :type prefix_set: set
        :param case: whether locus_tag is expected to be capitalized
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
        # Use genome_id from the object if check_value not given
        if prefix_set is None:
            prefix_set = {}
        if check_value is None:
            check_value = self.genome_id

        results = list()
        # If only check is to see that genome id is spelled correclty
        if only_typo:
            pattern = re.compile(check_value.lower())
            search_result = pattern.search(self.locus_tag.lower())
            if search_result is None:
                results.append("The genome ID is missing.")
        # Else if we're doing all checks
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
            result += "structured correctly."
            status = success
        else:
            result += "not structured correctly. " + " ".join(results)
            status = fail

        definition = f"Check if the locus_tag qualifier is structured " \
                     f"correctly for {self.id}."
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
                     f"{self.id}."
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
                     f"consistent for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_length(self, eval_id=None, success="correct", fail="error",
                     eval_def=None):
        """
        Checks that the tRNA is in the expected range of lengths. The
        average tRNA gene is 70-90bp in length, but it is not uncommon
        to identify high-scoring tRNAs in the 60-100bp range.
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
                     f"for {self.id}."
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

            # CHECK NOTE
        else:
            result += "missing (''). "
            result += "Note cannot be checked."
            status = "unchecked"

        definition = f"Check that the note is formatted properly for " \
                     f"{self.id}."
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

        definition = f"Make sure there is only one region for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
