"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from pdm_utils.functions import basic
from pdm_utils.constants import constants
from pdm_utils.classes import eval
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import re
from collections import OrderedDict


class Cds:
    """Class to hold data about a CDS feature."""

    def __init__(self):

        # The following attributes are common to any CDS.
        self.id = "" # Gene ID
        self.name = "" # Tends to be an integer for SEA-PHAGES.
        self.seqfeature = None # Biopython SeqFeature object.
        self.start = -1 # Genomic position
        self.stop = -1 # Genomic position
        self.coordinate_format = "" # Indexing format used for coordinates.
        self.orientation = "" #'forward', 'reverse', 'top', 'bottom', etc.
        self.parts = 0 # Number of regions that define the feature
        self.translation_table = 0
        self.translation = Seq("", IUPAC.protein) # Biopython amino acid Seq object.
        self.translation_length = 0
        self.seq = Seq("", IUPAC.ambiguous_dna) # Biopython nucleotide Seq object.
        self.length = 0

        # Information about the genome from which the feature is derived.
        self.genome_id = ""
        self.genome_length = -1


        # The following attributes are common to MySQL database.
        self.pham_id = 0
        self.domain_status = -1
        self.raw_description = "" # Raw gene description
        self.description = "" # Non-generic gene descriptions


        # The following attributes are common to
        # GenBank-formatted flat file records.
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = "" # Should be digit, but keep as string.
        self.gene = "" # Tends to be an integer, but not guaranteed.
        self.raw_product = ""
        self.raw_function = ""
        self.raw_note = ""
        self.product = ""
        self.function = ""
        self.note = ""

        # The following attributes are usefule for processing data
        # from various data sources.
        self.evaluations = []
        self.type = ""
        self._start_stop_orient_id = ()
        self._end_orient_id = ()
        self._start_end_id = ()




    def set_locus_tag(self, tag="", delimiter="_", check_value=None):
        """Set locus tag and parse the locus_tag feature number.

        :param tag:
            Input locus_tag data.
        :type tag: str
        :param delimiter:
            Value used to split locus_tag data.
        :type delimiter: str
        :param check_value:
            Indicates genome name or other value that will be used to parse
            the locus_tag to identify the feature number. If no check_value
            is provided, the genome_id attribute is used.
        :type check_value: str
        """

        # MySQL database-output format
        if tag is None:
            tag = ""
        self.locus_tag = tag
        if check_value is None:
            check_value = self.genome_id
        pattern = re.compile(check_value.lower())
        parts = tag.split(delimiter)

        index = 0
        found = False
        while (index < len(parts) and not found):
            search_result = pattern.search(parts[index].lower())
            if search_result != None:
                found = True
            else:
                index += 1
        if found:
            if index == len(parts) - 1:
                value = parts[index][search_result.end():]
            else:
                value = parts[-1]
        else:
            value = parts[-1]
        if value.isdigit():
            self._locus_tag_num = value



    def set_name(self, value=None):
        """Set the feature name.

        Ideally, the name of the CDS will be an integer. This information
        can be stored in multiple fields in the GenBank-formatted flat file.
        The name is first derived from the 'gene' qualifier, then
        the 'locus_tag' qualifier, and finally left empty.

        :param value:
            Indicates a value that should be used to directly set
            the name regardless of the 'gene' and '_locus_tag_num'
            attributes.
        :type value: str
        """

        # 1. PECAAN Draft:
        #    The 'gene' qualifier should be present and contain an integer.
        # 2. New SEA-PHAGES Final:
        #    The 'gene' qualifier should be present and contain an integer.
        # 3. SEA-PHAGES Final in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.
        # 4. Non-SEA-PHAGES in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.

        if value is not None:
            self.name = value
        elif self.gene != "":
            self.name = self.gene
        elif self._locus_tag_num != "":
            self.name = self._locus_tag_num
        else:
            self.name = ""


    def set_description(self, value):
        """Set the raw and processed description attributes.

        :param value:
            Indicates which reference attributes are used
            to set the attributes ('product', 'function', 'note').
        :type value: str
        """

        if value == "product":
            self.raw_description = self.raw_product
            self.description = self.product
        elif value == "function":
            self.raw_description = self.raw_function
            self.description = self.function
        elif value == "note":
            self.raw_description = self.raw_note
            self.description = self.note
        else:
            pass


    def get_begin_end(self):
        """Get feature coordinates in transcription begin-end format.

        :returns:
            (Begin, End) Start and stop coordinates ordered by which coordinate
            indicates the transcriptional beginning and end of the feature.
        :rtype: tuple
        """

        # Ensure format of orientation info.
        orientation = basic.reformat_strand(self.orientation, "fr_long")

        if orientation == "forward":
            begin = self.start
            end = self.stop
        elif orientation == "reverse":
            begin = self.stop
            end = self.start
        else:
            begin = -1
            end = -1

        return (begin, end)

    def translate_seq(self):
        """Translate the CDS nucleotide sequence.

        Use Biopython to translate the nucleotide sequece.
        The method expects the nucleotide sequence to be a valid CDS
        sequence in which:

            1. it begins with a valid start codon,
            2. it ends with a stop codon,
            3. it contains only one stop codon,
            4. its length is divisible by 3,
            5. it translates non-standard start codons to methionine.

        If these criteria are not met, an empty Seq object is returned.

        :returns: Amino acid sequence
        :rtype: Seq
        """

        try:
            translation = self.seq.translate(table=self.translation_table,
                                             cds=True)
        except:
            translation = Seq("", IUPAC.protein)
        return translation


    def set_translation(self, value=None, translate=False):
        """Set translation and its length.

        The translation is coerced into a Biopython Seq object.
        If no input translation value is provided, the translation is
        generated from the parent genome nucleotide sequence.
        If an input translation value is provided, the 'translate'
        parameter has no impact.

        :param value: Amino acid sequence
        :type value: str or Seq
        :param translate:
            Indicates whether the translation should be generated
            from the parent genome nucleotide sequence.
        :type translate: bool
        """

        if isinstance(value, Seq):
            self.translation = value.upper()
        elif value is not None:
            try:
                self.translation = Seq(value.upper(), IUPAC.protein)
            except:
                self.translation = Seq("", IUPAC.protein)
        elif translate:
            self.translation = self.translate_seq()
        else:
            self.translation = Seq("", IUPAC.protein)
        self.translation_length = len(self.translation)


    def set_translation_table(self, value):
        """Set translation table integer.

        :param value:
            Translation table that should be used to generate the translation.
        :type value: int
        """

        try:
            self.translation_table = int(value)
        except:
            self.translation_table = -1


    def set_orientation(self, value, format, case=False):
        """Sets orientation based on indicated format.

        Relies on the  `reformat_strand` function to manage orientation data.

        :param value: Input orientation value.
        :type value: misc.
        :param format: Indicates how the orientation data should be formatted.
        :type format: str
        :param case: Indicates whether the output orientation data should be cased.
        :type case: bool
        """
        self.orientation = basic.reformat_strand(value, format, case)


    def set_location_id(self):
        """Create a tuple of feature location data.

        For start and stop coordinates of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "start" and "stop" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates.
        If only the feature transcription "end" coordinate is used,
        orientation information is required.
        If transcription "begin" and "end" coordinates are used
        instead of "start" and "stop" coordinates,
        no orientation information is required.
        """
        self._start_stop_orient_id = (self.start, self.stop, self.orientation)

        begin, end = self.get_begin_end()
        self._end_orient_id = (end, self.orientation)
        self._start_end_id = (begin, end)


    def reformat_start_and_stop(self, new_format):
        """Convert start and stop coordinates to new coordinate format.
        This also updates the coordinate format attribute to reflect
        change.

        Relies on the `reformat_coordinates` function.

        :param new_format: Indicates how coordinates should be formatted.
        :type new_format: str
        """
        new_start, new_stop = basic.reformat_coordinates(
            self.start, self.stop, self.coordinate_format, new_format)
        if (new_format != self.coordinate_format):
            self.start = new_start
            self.stop = new_stop
            self.coordinate_format = new_format


    def set_nucleotide_length(self, seq=False):
        """Set the length of the nucleotide sequence.

        :param seq:
            Nucleotide length can be computed multiple ways.
            If the 'seq' parameter is False, the 'start' and 'stop'
            coordinates are used to determine the length (and take into
            account the 'coordinate_format' attribute. However, for
            compound features, this value may not be accurate.
            If the 'seq' parameter is True, the nucleotide sequence is used
            to determine the length. When the sequence reflects the
            entire feature (including compound features), the length
            will be accurate.
        :type seq: bool
        """

        if seq:
            self.length = len(self.seq)
            pass
        else:
            if self.coordinate_format == "0_half_open":
                self.length = self.stop - self.start
            elif self.coordinate_format == "1_closed":
                self.length = self.stop - self.start + 1
            else:
                self.length = -1


    def set_nucleotide_sequence(self, value=None, parent_genome_seq=None):
        """Set the nucleotide sequence of the feature.

        This method can directly set the attribute from a supplied 'value',
        or it can retrieve the sequence from the parent genome using
        Biopython. In this latter case, it relies on a Biopython SeqFeature
        object for the sequence extraction method and coordinates.
        If this object was generated from a Biopython-parsed
        GenBank-formatted flat file, the coordinates are by default
        '0-based half-open', the object contains coordinates for every
        part of the feature (e.g. if it is a compound feature) and
        fuzzy locations. As a result, the length of the retrieved sequence
        may not exactly match the length indicated from the 'start' and 'stop'
        coordinates.
        If the nucleotide sequence 'value' is provided, the
        'parent_genome_seq' does not impact the result.

        :param value: Input nucleotide sequence
        :type value: str of Seq
        :param parent_genome_seq: Input parent genome nucleotide sequence.
        :type parent_genome_seq: Seq
        """

        if isinstance(value, Seq):
            self.seq = value.upper()
        elif value is not None:
            try:
                self.seq = Seq(value.upper(), IUPAC.ambiguous_dna)
            except:
                self.seq = Seq("", IUPAC.ambiguous_dna)
        elif parent_genome_seq is not None:
            try:
                self.seq = self.seqfeature.extract(parent_genome_seq)
            except:
                # TODO not sure if this is the proper alphabet to use.
                self.seq = Seq("", IUPAC.ambiguous_dna)
        else:
            self.seq = Seq("", IUPAC.ambiguous_dna)


    # TODO Owen unittest for added steps.
    def set_seqfeature(self):
        """Set the 'seqfeature' attribute.

        The 'seqfeature' attribute stores a Biopython SeqFeature object,
        which contains methods valuable to extracting sequence data
        relevant to the feature.
        """

        # SeqFeature methods rely on coordinates in 0-based half-open
        # format and orientation to be numeric.
        new_start, new_stop = \
            basic.reformat_coordinates(self.start,
                                       self.stop,
                                       self.coordinate_format,
                                       "0_half_open")

        new_strand = basic.reformat_strand(self.orientation, "numeric")
        if self.start <= self.stop:
            self.seqfeature = SeqFeature(FeatureLocation(new_start, new_stop),
                                     strand=new_strand, type="CDS")
        else:
            self.seqfeature = SeqFeature(CompoundLocation
                    ([FeatureLocation(new_start, self.genome_length),
                    FeatureLocation(0, new_stop)]),
                    strand=new_strand, type="CDS")
        self.seqfeature.qualifiers = self.get_qualifiers()

    # TODO Owen unittest.
    def get_qualifiers(self):
        """Helper function that uses cds data to populate
        the qualifiers SeqFeature attribute

        :returns:
            qualifiers(dictionary) is a dictionary with the
            formating of BioPython's SeqFeature qualifiers
            attribute.
        """

        qualifiers = OrderedDict()
        qualifiers["gene"] = [self.name]
        if self.locus_tag != "":
            qualifiers["locus_tag"] = [self.locus_tag]
        qualifiers["note"] = ["gp{}".format(self.name)]
        qualifiers["codon_start"] = ["1"]
        qualifiers["transl_table"] = ["11"]
        if self.raw_description != "":
            qualifiers["product"] = [self.raw_description]
        qualifiers["id"] = [self.id]
        qualifiers["translation"] = [self.translation]

        return qualifiers




    # Evaluations.
    def set_eval(self, eval_id, definition, result, status):
        """Constructs and adds an Eval object to the evaluations list.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param definition: Description of the evaluation.
        :type definition: str
        :param result: Description of the outcome of the evaluation.
        :type result: str
        :param status: Outcome of the evaluation.
        :type status: str
        """
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_attribute(self, attribute, check_set, expect=False, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check that the attribute value is valid.

        :param attribute: Name of the CDS object attribute to evaluate.
        :type attribute: str
        :param check_set:
            Set of reference ids.
        :type check_set: set
        :param expect:
            Indicates whether the attribute value is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param success: Default status if the outcome is a success.
        :type success: str
        :param fail: Default status if the outcome is not a success.
        :type fail: str
        :param eval_def: Description of the evaluation.
        :type eval_def: str
        """
        try:
            test = True
            value1 = getattr(self, attribute)
        except:
            test = False
            value1 = None
        if test:
            value1_short = basic.truncate_value(str(value1), 30, "...")
            result = f"The {attribute} value '{value1_short}' is "

            value2 = basic.check_value_expected_in_set(
                        value1, check_set, expect)
            if value2:
                result = result + "valid."
                status = success
            else:
                result = result + "not valid."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"
        definition = f"Check the value of the '{attribute}' attribute for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_magnitude(self, attribute, expect, ref_value, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check that the magnitude of a numerical attribute is valid.

        :param attribute: same as for check_attribute().
        :param expect:
            Comparison symbol indicating direction of magnitude (>, =, <).
        :type expect: str
        :param ref_value: Numerical value for comparison.
        :type ref_value: int, float, datetime
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        try:
            test = True
            query_value = getattr(self, attribute)
        except:
            test = False
            query_value = None
        if test:
            result = f"The {attribute} value {query_value} is "
            if query_value > ref_value:
                compare = ">"
                result = result + "greater than "
            elif query_value == ref_value:
                compare = "="
                result = result + "equal to "
            else:
                compare = "<"
                result = result + "less than "
            result = result + f"{ref_value}, which is "
            if compare == expect:
                result = result + "expected."
                status = success
            else:
                result = result + "not expected."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"
        definition = f"Check the magnitude of the '{attribute}' attribute for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_translation(self, eval_id=None, success="correct",
                          fail="error", eval_def=None):
        """Check that the current and expected translations match.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """

        translation = self.translate_seq()
        exp_len = len(translation)
        result = f"The translation length ({self.translation_length}) "
        if self.translation_length < exp_len:
            result = result + f"is shorter than expected ({exp_len})."
            status = fail
        elif self.translation_length > exp_len:
            result = result + f"is longer than expected ({exp_len})."
            status = fail
        elif self.translation != translation:
            result = result + (f"is as expected, but the "
                               f"translation itself ({self.translation}) is "
                               f"different than expected ({translation}).")
            status = fail
        else:
            result = result + "and sequence are correct."
            status = success
        definition = f"Check that {self.id} contains the expected translation."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_amino_acids(self, check_set=set(), eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Check whether all amino acids in the translation are valid.

        :param check_set: Set of valid amino acids.
        :type check_set: set
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - check_set
        result = "The translation contains "
        if len(amino_acid_error_set) > 0:
            aae_string = ", ".join(amino_acid_error_set)
            result = result + f"the following unexpected amino acids: {aae_string}."
            status = fail
        else:
            result = result + "no unexpected amino acids."
            status = success

        definition = f"Check if all amino acids in the translation are expected for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_orientation(self, format="fr_short", case=True, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Check if orientation is set appropriately.

        Relies on the `reformat_strand` function to manage orientation data.

        :param format: Indicates how coordinates should be formatted.
        :type format: str
        :param case: Indicates whether the orientation data should be cased.
        :type case: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        expected_orient = basic.reformat_strand(self.orientation,
                                                format=format,
                                                case=case)
        result = f"The orientation is {self.orientation}, and it is formatted "
        if self.orientation == expected_orient:
            result = result + "correctly."
            status = success
        else:
            result = result + "incorrectly."
            status = fail
        definition = f"Check if the orientation is set correctly for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


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
        definition = f"Check if the locus_tag qualifier is structured correctly for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_gene_structure(self, eval_id=None, success="correct",
                             fail="error", eval_def=None):
        """Check if the gene qualifier contains an integer.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        try:
            value = int(self.gene)
        except:
            value = self.gene

        result = f"The gene qualifier is {self.gene}. It "
        if isinstance(value, int):
            result = result + "contains an integer, as expected."
            status = success
        else:
            result = result + "does not contain an integer, which is not expected."
            status = fail
        definition = f"Check if the gene qualifier contains an integer for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_compatible_gene_and_locus_tag(self, eval_id=None,
                                            success="correct", fail="error",
                                            eval_def=None):
        """Check if gene and locus_tag attributes contain identical numbers.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        result = (f"The numbers in the gene ({self.gene}) and "
                  f"locus_tag ({self._locus_tag_num}) qualifiers are ")
        if self.gene == self._locus_tag_num:
            result = result + "consistent."
            status = success
        else:
            result = result + "not consistent."
            status = fail
        definition = f"Check if the gene and locus_tag numbers are consistent for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_description_field(self, attribute="product", eval_id=None,
                                success="correct", fail="error", eval_def=None):
        """Check if there are CDS descriptions in unexpected fields.

        Evaluates whether the indicated attribute is empty or generic,
        and other fields contain non-generic data.

        :param attribute:
            Indicates the reference attribute for the evaluation
            ('product', 'function', 'note').
        :type attribute: str
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """

        description = ""
        description_set = set()

        if attribute == "product":
            description = self.product
        else:
            if self.product != "":
                description_set.add(self.product)

        if attribute == "function":
            description = self.function
        else:
            if self.function != "":
                description_set.add(self.function)

        if attribute == "note":
            description = self.note
        else:
            if self.note != "":
                description_set.add(self.note)

        result = f"The CDS description is '{description}', "
        if (description == "" and len(description_set) > 0):
            d_string = ", ".join(description_set)
            result = result + ("but there is a non-generic description in "
                               "at least one other qualifier: "
                               f"'{d_string}'. This is not expected.")
            status = fail
        else:
            result = result + "as expected."
            status = success
        definition = ("Check if there is a discrepancy "
                      f"between description fields for {self.id}.")
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO this should probably be implemented at the genome level,
    # calculating how many CDS features contain non-generic data.
    def check_generic_data(self, attribute=None, eval_id=None,
                           success="correct", fail="error", eval_def=None):
        """Check if the indicated attribute contains generic data.

        :param attribute:
            Indicates the attribute for the evaluation
            ('product', 'function', 'note').
        :type attribute: str
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """

        if attribute == "product":
            original = self.raw_product
            processed = self.product
        elif attribute == "function":
            original = self.raw_function
            processed = self.function
        elif attribute == "note":
            original = self.raw_note
            processed = self.note
        else:
            original = ""
            processed = ""

        if original == processed:
            result = f"The '{attribute}' field is correct."
            status = success
        else:
            result = f"The '{attribute}' field is not correct."
            status = fail
        definition = f"Check if the '{attribute}' field contains generic data for {self.id}."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)



    # TODO can be implemented once there are pre-defined lists available
    # of acceptable and non-acceptable gene descriptions.
    # TODO unittest.
    # def check_valid_description(self, check_set=None, attribute=None,
    #                             eval_id=None):
    #     """Check if the CDS description in the indicated attribute is valid."""
    #     pass
