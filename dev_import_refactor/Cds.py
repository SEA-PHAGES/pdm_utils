"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import functions_general




class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        self.type_id = "" #Feature type: CDS, GenomeBoundary,or tRNA
        self.left_boundary = "" #Position of left boundary of feature, 0-indexed
        self.right_boundary = "" #Position of right boundary of feature, 0-indexed
        self.strand = "" #'forward', 'reverse', or 'NA'
        self.compound_parts = 0 # number of regions that form the feature
        self.translation = ""
        self.gene_length = "" #stored in Phamerator
        self.translation_table = ""


        #Common to Phamerator
        self.phage_id = ""
        self.gene_id = "" #Gene ID comprised of PhageID and Gene name
        self.gene_name = ""
        self.notes = ""
        self.search_notes = "" #non-generic gene descriptions




        #Common to NCBI
        self.locus_tag = "" #Gene ID comprised of PhageID and Gene name
        self.gene_number = ""
        self.product_description = ""
        self.function_description = ""
        self.note_description = ""
        self.processed_product_description = ""
        self.processed__function_description = ""
        self.processed__note_description = ""






        # Computed datafields

        #Common to all
        self._translation_length = "" #computed internally
        self._amino_acid_errors = False
        self._start_end_strand_id = ""
        self._end_strand_id = ""
        self._boundary_error = False
        self._unmatched_error = False #keeps track if it is contains a match or not

        self._search_id = ""
        self._total_errors = 0





        #Common to NCBi
        self._locus_tag_missing = False
        self._locus_tag_typo = False #This can only be computed at the genome level since it uses the phage name
        self._description_field_error = False
        self._compound_parts_error = False #TODO implement this error












    # Define all attribute setters:



    #Common to all
    def set_strand(self,value):
        self.strand = functions_general.reformat_strand(value, "long")


    def set_translation(self,value):
        self.translation = value.upper()
        self._translation_length = len(self.translation)


    def set_type_id(self,value):
        self.__type_id = value


    def compute_amino_acid_errors(self,protein_alphabet_set):
        amino_acid_set = set(self.__translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set
        if len(amino_acid_error_set) > 0:
            self._amino_acid_errors = True


    def set_start_end_strand_id(self):
        #Create a tuple of feature location data.
        #For start and end of feature, it doesn't matter whether the feature is complex with a translational
        #frameshift or not. Retrieving the "start" and "end" attributes return the very beginning and end of
        #the feature, disregarding the inner "join" coordinates.
        self._start_end_strand_id = (str(self.__left_boundary),str(self.__right_boundary),self.__strand)

        #Since this id matched genes with different start sites,
        #the strand impacts whether the left or right boundary is used
        if self.__strand == 'forward':
            self._end_strand_id = (str(self.__right_boundary),self.__strand)
        elif self.__strand == 'reverse':
            self._end_strand_id = (str(self.__left_boundary),self.__strand)
        else:
            pass
    def set_unmatched_error(self):
        self._unmatched_error = True
    def compute_boundary_error(self):
        #Check if start and end coordinates are fuzzy
        if not (str(self.__left_boundary).isdigit() and str(self.__right_boundary).isdigit()):
            self._boundary_error = True






    #Common to Phamerator
    def set_phage_id(self,value):
        self.__phage_id = value
        self._search_id = remove_draft_suffix(self.__phage_id)
    def set_gene_id(self,value):
        self.__gene_id = value
    def set_gene_name(self,name):
        self.__gene_name = name
    def set_notes(self,value1,value2):
        self.__notes = value1
        self.__search_notes = value2
    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self._total_errors += 1
        if self.get_boundary_error():
            self._total_errors += 1
        if self.get_unmatched_error():
            self._total_errors += 1




    #Common to NCBI
    def set_locus_tag(self,value):
        self.locus_tag = value
        if self.locus_tag == '':
            self._locus_tag_missing = True
    def set_gene_number(self,value):
        self.__gene_number = value
    def set_locus_tag_typo(self):
        self._locus_tag_typo = True

    def compute_description_error(self):

        #If the product description is empty or generic, and the function or note descriptions are not, there is an error
        if self.__search_product_description == '' and \
            (self.__search_function_description != '' or \
            self.__search_note_description != ''):

            self._description_field_error = True

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











    # Define all attribute getters:


    #Common to all
    def get_left_boundary(self):
        return self.__left_boundary
    def get_right_boundary(self):
        return self.__right_boundary
    def get_type_id(self):
        return self.__type_id
    def get_strand(self):
        return self.__strand
    def get_amino_acid_errors(self):
        return self._amino_acid_errors
    def get_translation(self):
        return self.__translation
    def get_translation_length(self):
        return self.__translation_length
    def get_start_end_strand_id(self):
        return self._start_end_strand_id
    def get_end_strand_id(self):
        return self._end_strand_id
    def get_unmatched_error(self):
        return self._unmatched_error
    def get_boundary_error(self):
        return self._boundary_error


    #Common to Phamerator
    def get_gene_id(self):
        return self.__gene_id
    def get_gene_name(self):
        return self.__gene_name
    def get_notes(self):
        return self.__notes
    def get_search_notes(self):
        return self.__search_notes
    def get_phage_id(self):
        return self.__phage_id
    def get_search_id(self):
        return self._search_id
    def get_total_errors(self):
        return self._total_errors





    #Common to NCBI
    def get_locus_tag(self):
        return self.__locus_tag
    def get_gene_number(self):
        return self.__gene_number
    def get_product_description(self):
        return self.__product_description
    def get_function_description(self):
        return self.__function_description
    def get_note_description(self):
        return self.__note_description
    def get_search_product_description(self):
        return self.__search_product_description
    def get_search_function_description(self):
        return self.__search_function_description
    def get_search_note_description(self):
        return self.__search_note_description
    def get_locus_tag_missing(self):
        return self._locus_tag_missing
    def get_locus_tag_typo(self):
        return self._locus_tag_typo
    def get_description_field_error(self):
        return self._description_field_error
    def get_total_errors(self):
        return self._total_errors
