"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from collections import OrderedDict
import re

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from pdm_utils.constants import constants
from pdm_utils.classes import evaluation
from pdm_utils.functions import basic




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
        self._product_num = ""
        self._function_num = ""
        self._note_num = ""
        self._gene_num = ""

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
        # List of delimiters found in locus tags
        # "_" is the most common. The others are rare.
        delimiters = ["_", "-", "."]

        # MySQL database-output format
        if tag is None:
            tag = ""
        self.locus_tag = tag
        if delimiter is None:
            delimiter = basic.choose_most_common(tag, delimiters)
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

        # Remove generic 'gp' prefix if present (e.g. TRIXIE_gp10)
        # Sometimes locus tags contain number with character (e.g. TRIXIE_10A)
        if value.lower().startswith("gp"):
            value = value[2:]
        elif value.lower().startswith("orf"):
            value = value[3:]
        else:
            pass

        # If a numeric value can be identified, take it, other wise
        # take the entire value
        left, right = basic.split_string(value)
        if right != "":
            self._locus_tag_num = right
        else:
            self._locus_tag_num = value


    def set_name(self, value=None):
        """Set the feature name.

        Ideally, the name of the CDS will be an integer. This information
        can be stored in multiple fields in the GenBank-formatted flat file.
        The name is derived from one of several qualifiers.

        :param value:
            Indicates a value that should be used to directly set
            the name regardless of the 'gene' and '_locus_tag_num'
            attributes.
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
            value = ""
            list1 = ["_locus_tag_num", "_gene_num", "_product_num",
                     "_note_num", "_function_num"]

            # First see if any num attributes have a float.
            x = 0
            while (value == "" and x < len(list1)):
                name = getattr(self, list1[x])
                if basic.is_float(name):
                    value = name
                x += 1

            # Second see if any num attributes have non-empty values.
            # At this point, it's very unpredictable. Values could be like
            # 'terL' in the gene or like '10a' in the locus_tag.
            if value == "":
                list2 = ["gene"]
                list2.extend(list1)
                y = 0
                while (value == "" and y < len(list2)):
                    name = getattr(self, list2[y])
                    if name != "":
                        value = name
                    y += 1
        self.name = value


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
        self.set_num("_gene_num", value,
                     delimiter=delimiter, prefix_set=prefix_set)


    def set_description_field(self, attr, description,
                              delimiter=None, prefix_set=None):
        """Set a description attribute parsed from a description.

        :param attr: Attribute to set the description.
        :type attr: str
        :param description: Description data to parse. Also passed to set_num().
        :type description: str
        :param delimiter: Passed to set_num().
        :type delimiter: str
        :param prefix_set: Passed to set_num().
        :type prefix_set: set
        """
        # Used to set product/raw_product, function/raw_function, note/raw_note.
        # attr  # e.g. self.product
        raw_attr = "raw_" + attr  # e.g. self.raw_product
        num_attr = "_" + attr + "_num"  # e.g. self._product_num
        raw, processed = basic.reformat_description(description)
        setattr(self, attr, processed)
        setattr(self, raw_attr, raw)
        self.set_num(num_attr, description,
                     delimiter=delimiter, prefix_set=prefix_set)


    def set_num(self, attr, description, delimiter=None, prefix_set=None):
        """Set a number attribute from a description.

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


    def set_description(self, value):
        """Set the primary raw and processed description attributes.

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


    def set_nucleotide_length(self, seq=False, translation=False):
        """Set the length of the nucleotide sequence.

        Nucleotide length can be computed several different ways, including
        from the difference of the start and stop coordinates, the length
        of the transcribed nucleotide sequence, or the length of
        the translation. For compound features, using either the nucleotide or
        translation sequence is the accurate way to determine the
        true length of the feature, but 'length' may mean different things
        in different contexts.

        :param seq:
            Use the nucleotide sequence from the 'seq' attribute to
            compute the length.
        :type seq: bool
        :param translation:
            Use the translation sequence from the 'translation' attribute to
            compute the length.
        :type translation: bool
        """

        if seq:
            self.length = len(self.seq)
        elif translation:
            self.length = (len(self.translation) * 3) + 3 # include stop codon
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



#Needs unittests, however:
#Seqfeature retrieval and generation is clunky probably requires some
#over arching seqfeature generation.
#This cannot be in flat_files due to the need to import Cds which would
#create circular dependancies.
#May delay unittests until structure is revamped
    def set_seqfeature(self, type="CDS"):
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

        strand = basic.reformat_strand(self.orientation, "numeric")

        self.seqfeature = self.create_seqfeature(type, new_start, new_stop,
                                                                  strand)

    def create_seqfeature(self, type, start, stop, strand):
        if start <= stop:
            seqfeature = SeqFeature(FeatureLocation(start, stop),
                                            strand=strand, type=type)
        else:
            seqfeature = SeqFeature(CompoundLocation
                                  ([FeatureLocation(start, self.genome_length),
                                    FeatureLocation(0, stop)]),
                                            strand=strand, type=type)

        seqfeature.qualifiers = self.get_qualifiers(type)

        return seqfeature

    def get_qualifiers(self, type):
        """Helper function that uses cds data to populate
        the qualifiers SeqFeature attribute

        :returns:
            qualifiers(dictionary) is a dictionary with the
            formating of BioPython's SeqFeature qualifiers
            attribute.
        """
        qualifiers = OrderedDict()
        if self.description == "":
            product = "hypothetical protein"
        else:
            product = self.description

        if type == "CDS":
            qualifiers["gene"] = [self.name]
            if self.locus_tag != "":
                qualifiers["locus_tag"] = [self.locus_tag]
            qualifiers["note"] = ["".join(["gp", self.name])]
            qualifiers["codon_start"] = ["1"]
            qualifiers["transl_table"] = [self.translation_table]
            qualifiers["product"] = [product]
            qualifiers["translation"] = [self.translation]


        elif type == "Protein":
            qualifiers["product"] = [product]
        elif type == "gene":
           qualifiers["gene"] = [self.name]
           qualifiers["locus_tag"] = [self.locus_tag]

        return qualifiers




    # Evaluations.
    def set_eval(self, eval_id, definition, result, status):
        """Constructs and adds an Evaluation object to the evaluations list.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param definition: Description of the evaluation.
        :type definition: str
        :param result: Description of the outcome of the evaluation.
        :type result: str
        :param status: Outcome of the evaluation.
        :type status: str
        """
        evl = evaluation.Evaluation(eval_id, definition, result, status)
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
        definition = f"Check the value of the '{attribute}' attribute for {self.locus_tag} ({self.id})."
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
        definition = f"Check the magnitude of the '{attribute}' attribute for {self.locus_tag} ({self.id})."
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
        definition = f"Check that {self.locus_tag} ({self.id}) contains the expected translation."
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

        definition = f"Check if all amino acids in the translation are expected for {self.locus_tag} ({self.id})."
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
        definition = f"Check if the orientation is set correctly for {self.locus_tag} ({self.id})."
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
        definition = f"Check if the locus_tag qualifier is structured correctly for {self.locus_tag} ({self.id})."
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
        definition = f"Check if the gene qualifier contains an integer for {self.locus_tag} ({self.id})."
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
        definition = f"Check if the gene and locus_tag numbers are consistent for {self.locus_tag} ({self.id})."
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
                      f"between description fields for {self.locus_tag} ({self.id}).")
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
        definition = f"Check if the '{attribute}' field contains generic data for {self.locus_tag} ({self.id})."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)



    # TODO can be implemented once there are pre-defined lists available
    # of acceptable and non-acceptable gene descriptions.
    # TODO unittest.
    # def check_valid_description(self, check_set=None, attribute=None,
    #                             eval_id=None):
    #     """Check if the CDS description in the indicated attribute is valid."""
    #     pass
