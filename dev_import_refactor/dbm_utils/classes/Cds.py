"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from functions import basic
from classes import Eval
import re




class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        # Datafields from Phamerator database.

        # TODO: eventually change how id is computed.
        self.id = "" # Gene ID comprised of PhageID and Gene name
        self.name = ""
        self.type = "CDS"
        self.left_boundary = "" # Genomic position
        self.right_boundary = "" # Genomic position
        self.start = "" # Genomic position
        self.end = "" # Genomic position
        self.strand = "" #'forward', 'reverse', or 'NA'
        self.compound_parts = 0 # Number of regions that form the feature
        self.translation = "" # Biopython Seq object with protein alphabet.
        self.translation_table = ""

        # Indexing format for coordinates. Current valid formats:
        # 0-based half open (stored in Phamerator), 1-based closed.
        self.coordinate_format = ""


        # Common to Phamerator.
        self.parent_genome_id = ""
        self.parent_translation_table = ""


        self.primary_description = ""
        self.processed_primary_description = "" # Non-generic gene descriptions



        # Common to NCBI.


        self.seqfeature = None # Biopython SeqFeature object.
        # Thiis enables several QC checks
        # that utilize Biopython, such as retrieving the nucleotide
        # sequence from the parent genome and re-translating the CDS.
        self.seq = "" # Biopython Seq object containing nucleotide seq.

        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self.gene_number = ""
        self.product_description = ""
        self.function_description = ""
        self.note_description = ""
        self.processed_product_description = ""
        self.processed_function_description = ""
        self.processed_note_description = ""

        # Common to all.
        self.evaluations = []



        # Computed attributes:

        #Common to all. Computed internally.
        self._translation_length = 0
        self._nucleotide_length = 0 # Replaces gene_length, stored in Phamerator?
        self._left_right_strand_id = ()
        self._end_strand_id = ()
        self._start_end_id = ()



        # TODO may need to move these attributes to an inherited class.
        # # Computed error attributes.
        # self._amino_acid_errors = False
        # self._boundary_error = False
        # self._total_errors = 0
        #
        #
        # self._unmatched_error = False # Indicates if it is matched or not.
        #
        #
        #
        # #Common to NCBI
        # self._locus_tag_missing = False
        # self._locus_tag_typo = False # Requies phage name in Genome object.
        # self._description_field_error = False
        # self._compound_parts_error = False #TODO implement this error












    # Define all attribute setters:


    # def set_evaluation(self, type, message1 = None, message2 = None):
    #     """Creates an EvalResult object and adds it to the list of all
    #     evaluations.
    #     """
    #     if type == "warning":
    #         eval_object = Eval.construct_warning(message1, message2)
    #     elif type == "error":
    #         eval_object = Eval.construct_error(message1)
    #     else:
    #         eval_object = Eval.EvalResult()
    #     self.evaluations.append(eval_object)

    def set_parent_genome_id(self, value):
        """Set the parent_genome_id."""

        self.parent_genome_id = value

    def choose_description(self, value):
        """Set the primary description and processed primary description."""

        if value == "product":
            self.primary_description = self.product_description
            self.processed_primary_description = \
                self.processed_primary_description

        elif value == "function":
            self.primary_description = self.function_description
            self.processed_primary_description = \
                self.processed_function_description

        elif value == "note":
            self.primary_description = self.note_description
            self.processed_primary_description = \
                self.processed_note_description

        else:
            pass




    def set_translation(self, value):
        """Sets translation and determines length of translation.
        """
        self.translation = value.upper()
        self._translation_length = len(self.translation)

    def set_strand(self, value, format, case = False):
        """Sets strand based on indicated format.
        """
        self.strand = basic.reformat_strand(value, format, case)

    def set_start_end(self):
        """Determines which boundary coordinate is the start and end of
        gene based on the strand.
        """

        # Ensure format of strand info.
        strand = basic.reformat_strand(self.strand, "fr_long")

        if strand == "forward":
            self.start = self.left_boundary
            self.end = self.right_boundary
        elif strand == "reverse":
            self.start = self.right_boundary
            self.end = self.left_boundary
        else:
            pass

    def set_location_id(self):
        """ Create a tuple of feature location data.
        For left and right boundaries of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "left" and "right" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates.
        If only the feature "end" coordinate is used, strand information is
        required.
        If "start" and "end" coordinates are used instead of "left" and "right"
        coordinates, no strand information is required."""
        self._left_right_strand_id = (self.left_boundary, \
                                    self.right_boundary, \
                                    self.strand)
        self._end_strand_id = (self.end, self.strand)
        self._start_end_id = (self.start, self.end)


    def set_nucleotide_length(self):
        """From the set coordinates, determine the length of the
        nucleotide sequence. This method is not correct for
        non-compound features.
        """

        if self.coordinate_format == "0_half_open":
            self._nucleotide_length = \
                self.right_boundary - self.left_boundary

        elif self.coordinate_format == "1_closed":
            self._nucleotide_length = \
                self.right_boundary - self.left_boundary + 1

        else:
            self._nucleotide_length = -1



    def reformat_left_and_right_boundaries(self, new_format):
        """Convert left and right boundaries to new coordinate format.
        This also updates the coordinate format attribute to reflect
        change. However, it does not update start and end attributes."""

        new_left, new_right = \
            basic.reformat_coordinates(self.left_boundary, \
                                        self.right_boundary, \
                                        self.coordinate_format, \
                                        new_format)

        if (new_left != "" and new_right != ""):
            self.left_boundary = new_left
            self.right_boundary = new_right
            self.coordinate_format = new_format




    # TODO implement.
    # TODO unit test.
    def set_nucleotide_sequence(self, parent_genome):
        """Retrieves the nucleotide sequence from the parent genome."""
        # TODO pass the nucleotide sequence of the parent genome.
        # Use self.seqfeature object - this has methods to retrieve the
        # nucleotide sequence. The seqfeature object contains the entire
        # coordinates from the flat file record, including all compound
        # parts and fuzzy coordinates. So the retrieved sequence may
        # not exactly match the length indicated from the
        # self.left_boundary and self.right_boundary attributes.
        pass









    # Evaluations.

    def check_amino_acids(self, protein_alphabet_set):
        """Check whether all amino acids in the translation are valid.
        """
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
        eval = Eval.Eval(id = "CDS0001", \
                        definition = definition, \
                        result = result, \
                        status = status)
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
        eval = Eval.Eval(id = "CDS0002", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    # TODO this method can be improved by taking account coordinate
    # indexing format. The current implementation assumes only one format.
    # Also, instead of computing the nucleotide length, this can now
    # reference the self._nucleotide_length attribute.
    def check_lengths(self):
        """Confirm coordinates match translation length.
        This method can only be used on non-compound features."""

        if self.compound_parts == 1:
            length1 = self.right_boundary - self.left_boundary + 1
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
        eval = Eval.Eval(id = "CDS0003", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)



    # TODO unit test.
    def check_strand(self):
        """Check if strand is set appropriately."""

        strand = basic.reformat_strand(self.strand, format = "numeric")

        if (strand == 1 or self.strand == -1):
            result = "The feature strand is not determined: " \
                + str((self.left_boundary, self.right_boundary))
            status = "error"

        else:
            result = "Feature strand is correct."
            status = "correct"

        definition = "Check if the strand is set appropriately."
        eval = Eval.Eval(id = "CDS0004", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_boundaries(self):
        """Check if start and end coordinates are exact.
        This method assumes that if the coordinates are not exact, they
        have been set to -1 or are not integers."""

        if not (str(self.left_boundary).isdigit() and \
            str(self.right_boundary).isdigit()):

            result = "The feature boundaries are not determined: " \
                + str((self.left_boundary, self.right_boundary))
            status = "error"

        elif (self.left_boundary == -1 or \
            self.right_boundary == -1):

            # TODO unit test this elif clause.
            result = "The feature boundaries are not determined: " \
                + str((self.left_boundary, self.right_boundary))
            status = "error"

        else:
            result = "Feature boundaries are exact."
            status = "correct"

        definition = "Check if the left and right boundary coordinates " + \
                        "are exact or fuzzy."
        eval = Eval.Eval(id = "CDS0005", \
                        definition = definition, \
                        result = result, \
                        status = status)
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
        eval = Eval.Eval(id = "CDS0006", \
                        definition = definition, \
                        result = result, \
                        status = status)
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
        eval = Eval.Eval(id = "CDS0007", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_description(self):
        """If the product description is empty or generic, and the
        function or note descriptions are not, there is an error."""
        if self.processed_product_description == '' and \
            (self.processed_function_description != '' or \
            self.processed_note_description != ''):

            result = "The feature is missing a product description."
            status = "error"

        else:
            result = "The feature description is correct."
            status = "correct"

        definition = "Check if there is a discrepancy between description fields."
        eval = Eval.Eval(id = "CDS0008", \
                        definition = definition, \
                        result = result, \
                        status = status)
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
        eval = Eval.Eval(id = "CDS0009", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_translation_table_typo(self):
        """Check that translation table data matches data from parent Genome."""

        if self.translation_table != self.parent_translation_table:
            result = "The feature contains a translation table that " + \
                        "is different from the parent Genome translation table."
            status = "error"

        else:
            result = "The feature contains the expected translation table."
            status = "correct"

        definition = "Check that feature contains the exected translation table."
        eval = Eval.Eval(id = "CDS0010", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    # TODO implement.
    # TODO unit test.
    def check_translation(self):
        """Check that the current translation matches the expected
        translation."""
        # Once the nucleotide sequence is retrieved using the SeqFeature
        # methods, use Biopython to translate the nucleotide sequence.
        # This method can confirm that the CDS contains an initial start
        # codon, a final stop codon, no stop codons in the middle,
        # and that the translated product matches the translated product
        # stored in self.translation.
        pass






###
