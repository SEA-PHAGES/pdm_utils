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
        self.type_id = "" # Feature type: CDS
        self.left_boundary = "" # Genomic position, 0-indexed
        self.right_boundary = "" # Genomic position, 0-indexed
        self.start = "" # Genomic position, 0-indexed
        self.end = "" # Genomic position, 0-indexed
        self.strand = "" #'forward', 'reverse', or 'NA'
        self.compound_parts = 0 # Number of regions that form the feature
        self.translation = ""
        self.translation_table = ""

        # Indexing format for coordinates. Current valid formats:
        # 0-based half open (stored in Phamerator), 1-based closed.
        self.coordinate_format = ""

        # TODO: create coordinate indexing attribute, to be able to switch
        # coordinates between different types of indexing strategies.


        # Common to Phamerator.
        self.phage_id = ""

        # TODO: eventually change how Gene_ID is computed.
        self.gene_id = "" # Gene ID comprised of PhageID and Gene name
        self.gene_name = ""
        self.primary_description = ""
        self.processed_primary_description = "" # Non-generic gene descriptions

        # TODO: I don't think I need these anymore, since they have
        # been replaced by 'primary_description'
        # self.notes = ""
        # self.search_notes = "" # Non-generic gene descriptions


        # Common to NCBI.
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
        self._search_id = ""
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


    def set_evaluation(self, type, message1 = None, message2 = None):
        """Creates an EvalResult object and adds it to the list of all
        evaluations.
        """
        if type == "warning":
            eval_object = Eval.construct_warning(message1, message2)
        elif type == "error":
            eval_object = Eval.construct_error(message1)
        else:
            eval_object = Eval.EvalResult()
        self.evaluations.append(eval_object)

    def set_phage_id(self, value):
        """Set the phage_id and search_id at the same time, removing
        '_Draft' suffix if present."""

        self.phage_id = value
        self._search_id = basic.edit_draft_suffix(self.phage_id, "remove")


    # TODO unit test.
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



    # Evaluations.

    def check_amino_acids(self, protein_alphabet_set):
        """Check whether all amino acids in the translation are valid.
        """
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set

        if len(amino_acid_error_set) > 0:
            message = "There are unexpected amino acids in the translation: " \
                + str(amino_acid_error_set)
            self.set_evaluation("error", message)

    def check_translation_length(self):
        """Confirm that a translation is present."""
        if self._translation_length < 1:
            message = "There is no translation."
            self.set_evaluation("error", message)


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

                message = "The translation length and nucleotide length " + \
                            "do not match."
                self.set_evaluation("error", message)


    def check_boundaries(self):
        """Check if start and end coordinates are fuzzy."""
        if not (str(self.left_boundary).isdigit() and \
            str(self.right_boundary).isdigit()):

            message = "The feature boundaries are not determined: " \
                + str((self.left_boundary, self.right_boundary))
            self.set_evaluation("error", message)


    def check_locus_tag_present(self, expectation):
        """Check is status of locus tag matches expectations."""
        if expectation == "present":
            if self.locus_tag == "":
                message = "The feature has no locus tag."
                self.set_evaluation("error", message)

        elif expectation == "absent":
            if self.locus_tag != "":
                message = "The feature has a locus tag."
                self.set_evaluation("error", message)

        else:
            pass




    def check_locus_tag_typo(self, value):
        """Check is the locus tag contains potential typos."""
        pattern = re.compile(value.lower())
        search_result = pattern.search(self.locus_tag.lower())

        if search_result == None:

            message = "The feature locus tag has a typo."
            self.set_evaluation("error", message)



    def check_description(self):
        """If the product description is empty or generic, and the
        function or note descriptions are not, there is an error."""
        if self.processed_product_description == '' and \
            (self.processed_function_description != '' or \
            self.processed_note_description != ''):

            message = "The feature is missing a product description."
            self.set_evaluation("error", message)



    # TODO implement.
    # TODO unittest.
    def check_translation_table(self):
        """Check that translation table data is present."""
        pass





    # TODO implement.
    # TODO unit test.
    def check_translation(self):
        """Check that the current translation matches the expected
        translation."""
        # Using Biopython, retrieve the nucleotide sequence, translate the
        # sequence, and compare to the current value stored in the
        # translation attribute.
        pass






    #TODO this function cannot be implemented easily within the CDS feature,
    #since several of these require parameters that are determined by other
    #parts of the script, such as the amino acid alphabet and the _locus_tag
    #reference value.
    #TODO Unit test
    # def check_feature(self):
    #     self.check_translation(alphabet)
    #     self.check_boundaries()
    #     self.check_locus_tag_present()
    #     self.check_description()
    #     self.check_locus_tag_typo(value)








        # TODO: I will need to re-implement this for database_comparison script.
        # if self.get_unmatched_error():
        #     self._total_errors += 1










###
