"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import functions_general
import Eval




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
        self.gene_length = "" # Value stored in Phamerator
        self.translation_table = ""


        # Common to Phamerator.
        self.phage_id = ""
        self.gene_id = "" # Gene ID comprised of PhageID and Gene name = change this probably
        self.gene_name = ""
        self.notes = ""
        self.search_notes = "" # Non-generic gene descriptions


        # Common to NCBI.
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self.gene_number = ""
        self.product_description = ""
        self.function_description = ""
        self.note_description = ""
        self.processed_product_description = ""
        self.processed__function_description = ""
        self.processed__note_description = ""

        # Common to all.
        self.evaluations = []



        # Computed attributes:

        #Common to all.
        self._search_id = ""
        self._translation_length = "" # Computed internally
        self._left_right_strand_id = ()
        self._end_strand_id = ()



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

    def set_strand(self, value, format):
        """Sets strand based on indicated format.
        """
        self.strand = functions_general.reformat_strand(value, format)

    def set_translation(self, value):
        """Sets translation and determines length of translation.
        """
        self.translation = value.upper()
        self._translation_length = len(self.translation)




    def set_start_end(self):
        """Determines which boundary coordinate is the start and end of
        gene based on the strand.
        """

        # Ensure format of strand info.
        strand = functions_general.reformat_strand(self.strand, "long")

        if strand == "forward":

            self.start = self.left_boundary
            self.end = self.right_boundary

        elif strand == "reverse":

            self.start = self.right_boundary
            self.end = self.left_boundary

        else:
            pass





    def check_translation(self, protein_alphabet_set):
        """Check whether all amino acids in the translation are valid.
        """
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set

        if len(amino_acid_error_set) > 0:
            message = "There are unexpected amino acids in the translation: " \
                + str(amino_acid_error_set)
            self.set_evaluation("error", "message")

    def set_location_id(self):
        """ Create a tuple of feature location data.
        For left and right boundaries of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "left" and "right" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates."""
        self._left_right_strand_id = (self.left_boundary, \
                                    self.right_boundary, \
                                    self.strand)

        self._end_strand_id = (self.end, self.strand)




    #TODO Unit test
    def compute_boundary_error(self):
        #Check if start and end coordinates are fuzzy
        if not (str(self.left_boundary).isdigit() and str(self.right_boundary).isdigit()):
            self._boundary_error = True






    #Common to Phamerator
    #TODO Unit test
    def set_phage_id(self,value):
        self.phage_id = value
        self._search_id = remove_draft_suffix(self.phage_id)


    #TODO Unit test
    def set_notes(self,value1,value2):
        self.notes = value1
        self.search_notes = value2





    # Compute errors.

    #TODO Revamp.
    #TODO Unit test
    def set_locus_tag(self,value):
        self.locus_tag = value
        if self.locus_tag == '':
            self._locus_tag_missing = True









    # Methods that need revamped to create eval.




    #TODO revamp
    #TODO Unit test
    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self._total_errors += 1
        if self.get_boundary_error():
            self._total_errors += 1
        if self.get_unmatched_error():
            self._total_errors += 1


    #TODO revamp
    #TODO Unit test
    def compute_description_error(self):

        # If the product description is empty or generic, and the
        # function or note descriptions are not, there is an error
        if self.__search_product_description == '' and \
            (self.__search_function_description != '' or \
            self.__search_note_description != ''):

            self._description_field_error = True


    #TODO revamp
    #TODO Unit test
    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self._total_errors += 1
        if self.get_boundary_error():
            self._total_errors += 1
        if self._description_field_error:
            self._total_errors += 1
        if self._locus_tag_missing:
            self._total_errors += 1
        if self._locus_tag_typo:
            self._total_errors += 1
        if self.get_unmatched_error():
            self._total_errors += 1








###
