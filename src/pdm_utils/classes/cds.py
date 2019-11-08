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
        self.left = -1 # Genomic position
        self.right = -1 # Genomic position
        self.coordinate_format = "" # Indexing format used for coordinates.
        self.strand = "" #'forward', 'reverse', 'top', 'bottom', etc.
        self.parts = 0 # Number of regions that define the feature
        self.translation_table = 0
        self.translation = Seq("", IUPAC.protein) # Biopython amino acid Seq object.
        self.translation_length = 0
        self.seq = Seq("", IUPAC.ambiguous_dna) # Biopython nucleotide Seq object.
        self.length = 0

        # Information about the genome from which the feature is derived.
        self.genome_id = ""
        self.genome_length = -1


        # The following attributes are common to PhameratorDB.
        self.pham = 0 # TODO build method to implement this.
        self.description = "" # Raw gene description
        self.processed_description = "" # Non-generic gene descriptions


        # The following attributes are common to
        # GenBank-formatted flat file records.
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = "" # Should be digit, but keep as string.
        self.gene = "" # Tends to be an integer, but not guaranteed.
        self.product = ""
        self.function = ""
        self.note = ""
        self.processed_product = ""
        self.processed_function = ""
        self.processed_note = ""

        # The following attributes are usefule for processing data
        # from various data sources.
        self.evaluations = []
        self.type = ""
        self._left_right_strand_id = ()
        self._end_strand_id = ()
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
        """Set the description and processed_description attributes.

        :param value:
            Indicates which reference attributes are used
            to set the attributes ('product', 'function', 'note').
        :type value: str
        """

        if value == "product":
            self.description = self.product
            self.processed_description = self.processed_product
        elif value == "function":
            self.description = self.function
            self.processed_description = self.processed_function
        elif value == "note":
            self.description = self.note
            self.processed_description = self.processed_note
        else:
            pass


    def get_start_end(self):
        """Return the coordinates in start-end format."""

        # Ensure format of strand info.
        strand = basic.reformat_strand(self.strand, "fr_long")

        if strand == "forward":
            start = self.left
            end = self.right
        elif strand == "reverse":
            start = self.right
            end = self.left
        else:
            start = -1
            end = -1

        return (start, end)

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


    def set_strand(self, value, format, case=False):
        """Sets strand based on indicated format.

        Relies on the  `reformat_strand` function to manage strand data.

        :param value: Input strand value.
        :type value: misc.
        :param format: Indicates how the strand data should be formatted.
        :type format: str
        :param case: Indicates whether the output strand data should be cased.
        :type case: bool
        """
        self.strand = basic.reformat_strand(value, format, case)


    def set_location_id(self):
        """Create a tuple of feature location data.

        For left and right coordinates of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "left" and "right" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates.
        If only the feature "end" coordinate is used, strand information is
        required.
        If "start" and "end" coordinates are used instead of "left" and "right"
        coordinates, no strand information is required.
        """
        self._left_right_strand_id = (self.left, self.right, self.strand)

        start, end = self.get_start_end()
        self._end_strand_id = (end, self.strand)
        self._start_end_id = (start, end)


    def reformat_left_and_right(self, new_format):
        """Convert left and right coordinates to new coordinate format.
        This also updates the coordinate format attribute to reflect
        change.

        Relies on the `reformat_coordinates` function.

        :param new_format: Indicates how coordinates should be formatted.
        :type new_format: str
        """
        new_left, new_right = basic.reformat_coordinates(
            self.left, self.right, self.coordinate_format, new_format)
        if (new_format != self.coordinate_format):
            self.left = new_left
            self.right = new_right
            self.coordinate_format = new_format


    def set_nucleotide_length(self, seq=False):
        """Set the length of the nucleotide sequence.

        :param seq:
            Nucleotide length can be computed multiple ways.
            If the 'seq' parameter is False, the 'left' and 'right'
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
                self.length = \
                    self.right - self.left
            elif self.coordinate_format == "1_closed":
                self.length = \
                    self.right - self.left + 1
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
        may not exactly match the length indicated from the 'left' and 'right'
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
        # format and strand to be numeric.
        new_left, new_right = \
            basic.reformat_coordinates(self.left, \
                                       self.right, \
                                       self.coordinate_format, \
                                       "0_half_open")

        new_strand = basic.reformat_strand(self.strand, "numeric")
        if self.left <= self.right:
            self.seqfeature = SeqFeature(FeatureLocation(new_left, new_right),\
                                     strand=new_strand, type="CDS")
        else:
            self.seqfeature = SeqFeature(CompoundLocation\
                    ([FeatureLocation(new_left, self.genome_length),\
                    FeatureLocation(0, new_right)]),\
                    strand = new_strand, type = "CDS")
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
        if self.description != "":
            qualifiers["product"] = [self.description]
        qualifiers["id"] = [self.id]
        qualifiers["translation"] = [self.translation]

        return qualifiers



    # Evaluations.

    def check_translation_table(self, check_table=11, eval_id=None):
        """Check that the translation table is correct.

        :param check_table: Translation table used to check the translation.
        :type check_table: int
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.translation_table == check_table:
            result = "The translation table is correct."
            status = "correct"
        else:
            result = "The translation table is not correct."
            status = "error"
        definition = "Check that the translation table is correct."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_translation_present(self, eval_id=None):
        """Confirm that a translation is present.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if self.translation_length < 1:
            result = "A translation is not present."
            status = "error"
        else:
            result = "A translation is present."
            status = "correct"

        definition = "Check that there is a translation present."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_translation(self, eval_id=None):
        """Check that the current and expected translations match.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        translation = self.translate_seq()
        exp_len = len(translation)
        if self.translation_length < exp_len:
            result = (f"The translation length ({self.translation_length}) "
                     f"is shorter than expected ({exp_len}).")
            status = "error"
        elif self.translation_length > exp_len:
            result = (f"The translation length ({self.translation_length}) "
                     f"is longer than expected ({exp_len}).")
            status = "error"
        elif self.translation != translation:
            result = (f"The translation ({self.translation})is different "
                     f"than expected ({translation}).")
            status = "error"
        else:
            result = "The translation is correct."
            status = "correct"
        definition = "Check that the feature contains the expected translation."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_amino_acids(self, check_set=set(), eval_id=None):
        """Check whether all amino acids in the translation are valid.

        :param check_set: Set of valid amino acids.
        :type check_set: set
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - check_set

        if len(amino_acid_error_set) > 0:
            result = "There are unexpected amino acids in the translation: " \
                + str(amino_acid_error_set)
            status = "error"
        else:
            result = "There are no unexpected amino acid residues."
            status = "correct"

        definition = "Check validity of amino acid residues."
        evl = eval.Eval(eval_id, definition, result, status = status)
        self.evaluations.append(evl)


    def check_strand(self, format="fr_short", case=True, eval_id=None):
        """Check if strand is set appropriately.

        Relies on the `reformat_strand` function to manage strand data.

        :param format: Indicates how coordinates should be formatted.
        :type format: str
        :param case: Indicates whether the strand data should be cased.
        :type case: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        expected_strand = basic.reformat_strand(self.strand,
                                                format=format,
                                                case=case)
        if self.strand == expected_strand:
            result = "The feature strand is correct."
            status = "correct"
        else:
            result = "The feature strand is not correct."
            status = "error"
        definition = "Check if the strand is set appropriately."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_coordinates(self, eval_id=None):
        """Check if coordinates are exact.

        This method assumes that if the coordinates are not exact, they
        have been set to -1 or are not integers.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if not (isinstance(self.left, int) and isinstance(self.right, int)):
            result = ("The feature coordinates are not integers: "
                      + str((self.left, self.right)))
            status = "error"
        elif (self.left == -1 or self.right == -1):
            result = ("The feature coordinates are not determined: "
                      + str((self.left, self.right)))
            status = "error"
        else:
            result = "Feature coordinates are exact."
            status = "correct"

        definition = ("Check if the left and right boundary coordinates "
                      "are exact or fuzzy.")
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_locus_tag_present(self, expect=True, eval_id=None):
        """Check if status of locus tag matches expectations.

        :param expect:
            Indicates whether the locus_tag is expected to be present.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.locus_tag != "":
            present = True
        else:
            present = False

        if expect:
            if present:
                result = "The locus_tag qualifier is present."
                status = "correct"
            else:
                result = "The locus_tag qualifier is not present."
                status = "error"
        else:
            if present:
                result = "The locus_tag qualifier is present."
                status = "error"
            else:
                result = "The locus_tag qualifier is not present."
                status = "correct"
        definition = "Check if the locus_tag qualifier status is expected."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_locus_tag_structure(self, check_value=None, only_typo=False,
                                  prefix_set=set(), case=True, eval_id=None):
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
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
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
        if len(results) == 0:
            result = "The locus_tag qualifier is structured correctly."
            status = "correct"
        else:
            result = "The locus_tag qualifier is not structured correctly." \
                     + " ".join(results)
            status = "error"
        definition = "Check if the locus_tag qualifier is structured correctly."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    # TODO is this needed? CDS ids are not directly parsed from file.
    # They are now automatically generated. So there doesn't seem to be
    # a need to check that it is spelled correctly.
    # TODO implement.
    # TODO unittest
    # def check_id_typo(self, check_value=None, eval_id=None):
    #     """Check if the id contains potential typos."""
    #
    #     if check_value is None:
    #         check_value = self.id
    #
    #     pattern = re.compile(check_value.lower())
    #     search_result = pattern.search(self.id.lower())
    #
    #     if search_result == None:
    #         result = "The id has a typo."
    #         status = "error"
    #     else:
    #         result = "The id is correct."
    #         status = "correct"
    #     definition = "Check if the id contains a typo."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)


    def check_gene_present(self, expect=True, eval_id=None):
        """Check if the status of gene matches expectations.

        :param expect:
            Indicates whether the gene qualifier is expected to be present.
        :type expect: bool
        :param eval_id:
            Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.gene != "":
            present = True
        else:
            present = False

        if expect:
            if present:
                result = "The gene qualifier is present."
                status = "correct"
            else:
                result = "The gene qualifier is not present."
                status = "error"
        else:
            if present:
                result = "The gene qualifier is present."
                status = "error"
            else:
                result = "The gene qualifier is not present."
                status = "correct"
        definition = "Check if the gene status is expected."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_gene_structure(self, eval_id=None):
        """Check if the gene qualifier contains an integer.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        try:
            value = int(self.gene)
        except:
            value = self.gene

        if isinstance(value, int):
            result = "The gene qualifier contains an integer."
            status = "correct"
        else:
            result = "The gene qualifier does not contain an integer."
            status = "error"
        definition = "Check if the gene qualifier contains an integer."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_compatible_gene_and_locus_tag(self, eval_id=None):
        """Check if the gene and locus_tag attributes contain the same
        gene number.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if self.gene == self._locus_tag_num:
            result = "The gene and locus_tag numbers are consistent."
            status = "correct"
        else:
            result = "The gene and locus_tag numbers are not consistent."
            status = "error"
        definition = "Check if the gene and locus_tag numbers are consistent."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_description_field(self, attribute="product", eval_id=None):
        """Check if there are CDS descriptions in unexpected fields.

        Evaluates whether the indicated attribute is empty or generic,
        and other fields contain non-generic data.

        :param attribute:
            Indicates the reference attribute for the evaluation
            ('product', 'function', 'note').
        :type attribute: str
        :param eval_id:
            Unique identifier for the evaluation.
        :type eval_id: str
        """

        description = ""
        description_set = set()

        if attribute == "product":
            description = self.processed_product
        else:
            if self.processed_product != "":
                description_set.add(self.processed_product)

        if attribute == "function":
            description = self.processed_function
        else:
            if self.processed_function != "":
                description_set.add(self.processed_function)

        if attribute == "note":
            description = self.processed_note
        else:
            if self.processed_note != "":
                description_set.add(self.processed_note)

        if (description == "" and len(description_set) > 0):
            result = "The description is not in the expected field."
            status = "error"
        else:
            result = "The description is in the expected field."
            status = "correct"

        definition = \
            "Check if there is a discrepancy between description fields."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_generic_data(self, attribute=None, eval_id=None):
        """Check if the indicated attribute contains generic data.

        :param attribute:
            Indicates the attribute for the evaluation
            ('product', 'function', 'note').
        :type attribute: str
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if attribute == "product":
            original = self.product
            processed = self.processed_product
        elif attribute == "function":
            original = self.function
            processed = self.processed_function
        elif attribute == "note":
            original = self.note
            processed = self.processed_note
        else:
            original = ""
            processed = ""

        if original == processed:
            result = f"The '{attribute}' field is correct."
            status = "correct"
        else:
            result = f"The '{attribute}' field is not correct."
            status = "error"
        definition = f"Check if the '{attribute}' field contains generic data."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)



    # TODO can be implemented once there are pre-defined lists available
    # of acceptable and non-acceptable gene descriptions.
    # TODO unittest.
    # def check_valid_description(self, check_set=None, attribute=None,
    #                             eval_id=None):
    #     """Check if the CDS description in the indicated attribute is valid."""
    #     pass




###
