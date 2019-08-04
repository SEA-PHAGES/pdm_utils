"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from functions import basic
from classes import Eval
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re




class Cds:
    """Class to hold data about a CDS feature."""

    def __init__(self):

        # The following attributes are common to any CDS.

        # TODO: eventually change how id is computed.
        self.id = "" # Gene ID comprised of PhageID and Gene name
        self.name = ""
        self.genome_id = "" # Genome from which CDS feature is derived.
        self.left = "" # Genomic position
        self.right = "" # Genomic position
        self.start = "" # Genomic position
        self.end = "" # Genomic position
        self.strand = "" #'forward', 'reverse', 'top', 'bottom', etc.
        self.compound_parts = 0 # Number of regions that define the feature
        self.translation_table = ""
        self.translation = "" # Biopython amino acid Seq object.
        self._translation_length = 0
        self.seq = "" # Biopython nucleotide Seq object.
        self._length = 0 # Replaces gene_length, stored in Phamerator?
        self._left_right_strand_id = ()
        self._end_strand_id = ()
        self._start_end_id = ()

        # This enables several QC checks
        # that utilize Biopython, such as retrieving the nucleotide
        # sequence from the parent genome and re-translating the CDS.
        self.seqfeature = None # Biopython SeqFeature object.

        # Indexing format for coordinates. Current valid formats:
        # 0-based half open (stored in Phamerator), 1-based closed.
        self.coordinate_format = ""
        self.description = ""
        self.processed_description = "" # Non-generic gene descriptions


        # The following attributes are common to
        # GenBank-formatted flat file records.
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self.gene_number = ""
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



    def choose_description(self, value):
        """Set the primary description and processed primary description."""

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
        """

        try:
            translation = self.seq.translate(table=self.translation_table,
                                             cds=True)
        except:
            translation = Seq("", IUPAC.protein)
        return translation



    def set_translation(self, value=None, translate=False):
        """Set translation and its length."""

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
        self._translation_length = len(self.translation)


    def set_translation_table(self, value):
        """Set translation table integer."""

        try:
            self.translation_table = int(value)
        except:
            self.translation_table = -1


    def set_strand(self, value, format, case=False):
        """Sets strand based on indicated format."""
        self.strand = basic.reformat_strand(value, format, case)

    def set_start_end(self):
        """Determines which boundary coordinate is the start and end of
        gene based on the strand.
        """

        # Ensure format of strand info.
        strand = basic.reformat_strand(self.strand, "fr_long")

        if strand == "forward":
            self.start = self.left
            self.end = self.right
        elif strand == "reverse":
            self.start = self.right
            self.end = self.left
        else:
            pass

    def set_location_id(self):
        """Create a tuple of feature location data.

        For left and right boundaries of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "left" and "right" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates.
        If only the feature "end" coordinate is used, strand information is
        required.
        If "start" and "end" coordinates are used instead of "left" and "right"
        coordinates, no strand information is required.
        """
        self._left_right_strand_id = (self.left, \
                                      self.right, \
                                      self.strand)
        self._end_strand_id = (self.end, self.strand)
        self._start_end_id = (self.start, self.end)


    def set_nucleotide_length(self):
        """From the set coordinates, determine the length of the
        nucleotide sequence. This method is not correct for
        non-compound features.
        """

        if self.coordinate_format == "0_half_open":
            self._length = \
                self.right - self.left
        elif self.coordinate_format == "1_closed":
            self._length = \
                self.right - self.left + 1
        else:
            self._length = -1


    def reformat_left_and_right_boundaries(self, new_format):
        """Convert left and right boundaries to new coordinate format.
        This also updates the coordinate format attribute to reflect
        change. However, it does not update start and end attributes."""

        new_left, new_right = \
            basic.reformat_coordinates(self.left, \
                                       self.right, \
                                       self.coordinate_format, \
                                       new_format)

        if (new_left != "" and new_right != ""):
            self.left = new_left
            self.right = new_right
            self.coordinate_format = new_format


    def set_nucleotide_sequence(self, value=None, parent_genome_seq=None):
        """Set the nucleotide sequence of the feature.

        The method can directly set the attribute from a supplied 'value',
        or it can retrieve the sequence from the parent genome."""

        if isinstance(value, Seq):
            self.seq = value.upper()
        elif value is not None:
            try:
                self.seq = Seq(value.upper(), IUPAC.ambiguous_dna)
            except:
                self.seq = Seq("", IUPAC.ambiguous_dna)
        elif parent_genome_seq is not None:
            # The seqfeature object contains the entire
            # coordinates from the flat file record, including all compound
            # parts and fuzzy coordinates. So the retrieved sequence may
            # not exactly match the length indicated from the
            # self.left and self.right attributes.
            try:
                self.seq = self.seqfeature.extract(parent_genome_seq)
            except:
                # TODO not sure if this is the proper alphabet to use.
                self.seq = Seq("", IUPAC.ambiguous_dna)
        else:
            self.seq = Seq("", IUPAC.ambiguous_dna)

        # TODO consider whether _length should be automatically set,
        # parallel to genome method.
        # self._length = len(self.seq)






    # Evaluations.

    def check_amino_acids(self, protein_alphabet_set):
        """Check whether all amino acids in the translation are valid."""
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set

        if len(amino_acid_error_set) > 0:
            result = "There are unexpected amino acids in the translation: " \
                + str(amino_acid_error_set)
            status = "error"
        else:
            result = "There are no unexpected amino acid residues."
            status = "correct"

        definition = "Check validity of amino acid residues."
        eval = Eval.Eval("CDS0001", definition, result, status = status)
        self.evaluations.append(eval)

    def check_translation_length(self):
        """Confirm that a translation is present."""
        if self._translation_length < 1:
            result = "There is no translation."
            status = "error"
        else:
            result = "Translation is identified."
            status = "correct"

        definition = "Confirm there is a translation."
        eval = Eval.Eval("CDS0002", definition, result, status)
        self.evaluations.append(eval)


    # TODO this method can be improved by taking into account coordinate
    # indexing format. The current implementation assumes only one format.
    # Also, instead of computing the nucleotide length, this can now
    # reference the self._length attribute.
    def check_lengths(self):
        """Confirm coordinates match translation length.
        This method can only be used on non-compound features."""

        if self.compound_parts == 1:
            length1 = self.right - self.left + 1
            length2 = (self._translation_length * 3) + 3

            if length1 != length2:
                result = "The translation length and nucleotide length " + \
                            "do not match."
                status = "error"
            else:
                result = "Translation and nucleotide lengths match."
                status = "correct"

        else:
            result = "Translation and nucleotide lengths not compared."
            status = "untested"

        definition = "Confirm the translation and nucleotide lengths."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)



    # TODO unit test.
    def check_strand(self):
        """Check if strand is set appropriately."""

        strand = basic.reformat_strand(self.strand, format = "numeric")

        if (strand == 1 or self.strand == -1):
            result = "The feature strand is not determined: " \
                + str((self.left, self.right))
            status = "error"

        else:
            result = "Feature strand is correct."
            status = "correct"

        definition = "Check if the strand is set appropriately."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_boundaries(self):
        """Check if start and end coordinates are exact.
        This method assumes that if the coordinates are not exact, they
        have been set to -1 or are not integers."""

        if not (str(self.left).isdigit() and \
            str(self.right).isdigit()):

            result = "The feature boundaries are not determined: " \
                + str((self.left, self.right))
            status = "error"

        elif (self.left == -1 or \
            self.right == -1):

            # TODO unit test this elif clause.
            result = "The feature boundaries are not determined: " \
                + str((self.left, self.right))
            status = "error"

        else:
            result = "Feature boundaries are exact."
            status = "correct"

        definition = "Check if the left and right boundary coordinates " + \
                        "are exact or fuzzy."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_locus_tag_present(self, expectation):
        """Check if status of locus tag matches expectations."""
        if expectation == "present":
            if self.locus_tag == "":
                result = "The feature has no locus tag."
                status = "error"
            else:
                result = "The locus_tag is as expected."
                status = "correct"

        elif expectation == "absent":
            if self.locus_tag != "":
                result = "The feature has a locus tag."
                status = "error"
            else:
                result = "The locus_tag is as expected."
                status = "correct"

        # TODO unit test.
        else:
            result = "The presence/absence of the locus_tag was not evaluated."
            status = "untested"

        definition = "Check if the locus_tag status is expected."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)




    def check_locus_tag_typo(self, value):
        """Check if the locus tag contains potential typos."""
        pattern = re.compile(value.lower())
        search_result = pattern.search(self.locus_tag.lower())

        if search_result == None:

            result = "The feature locus tag has a typo."
            status = "error"

        else:
            result = "The feature locus tag is correct."
            status = "correct"

        definition = "Check if the locus_tag contains a typo."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_description(self):
        """If the product description is empty or generic, and the
        function or note descriptions are not, there is an error."""
        if self.processed_product == '' and \
            (self.processed_function != '' or \
            self.processed_note != ''):

            result = "The feature is missing a product description."
            status = "error"

        else:
            result = "The feature description is correct."
            status = "correct"

        definition = "Check if there is a discrepancy between description fields."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_translation_table_present(self):
        """Check that translation table data is present."""
        if self.translation_table == "":
            result = "The feature is missing a translation table."
            status = "error"
        else:
            result = "The feature contains a translation table."
            status = "correct"

        definition = "Check that translation table data is present."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)

    def check_translation_table_typo(self, parent_trans_table=11):
        """Check that translation table data matches data from parent Genome."""

        if self.translation_table != parent_trans_table:
            result = "The feature contains a translation table that " + \
                        "is different from the parent Genome translation table."
            status = "error"

        else:
            result = "The feature contains the expected translation table."
            status = "correct"

        definition = "Check that feature contains the expected translation table."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_translation(self):
        """Check that the current and expected translations match."""

        translation = self.translate_seq()
        if self._translation_length < len(translation):
            result = "The translation length is shorter than expected."
            status = "error"
        elif self._translation_length > len(translation):
            result = "The translation length is longer than expected."
            status = "error"
        elif self.translation != translation:
            result = "The translation is different than expected."
            status = "error"
        else:
            result = "The translation is correct."
            status = "correct"
        definition = "Check that the feature contains the expected translation."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)







###
