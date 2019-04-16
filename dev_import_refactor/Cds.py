"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""






class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        self.__type_id = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.__left_boundary = '' #Position of left boundary of feature, 0-indexed
        self.__right_boundary = '' #Position of right boundary of feature, 0-indexed
        self.__strand = '' #'forward', 'reverse', or 'NA'
        self.__translation = ''
        self.__translation_length = ''


        #Common to Phamerator
        self.__phage_id = ''
        self.__gene_id = '' #Gene ID comprised of PhageID and Gene name
        self.__gene_name = ''
        self.__notes = ''
        self.__search_notes = '' #non-generic gene descriptions




        #Common to NCBI
        self.__locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.__gene_number = ''
        self.__product_description = ''
        self.__function_description = ''
        self.__note_description = ''
        self.__search_product_description = ''
        self.__search_function_description = ''
        self.__search_note_description = ''






        # Computed datafields

        #Common to all
        self.__amino_acid_errors = False
        self.__start_end_strand_id = ''
        self.__end_strand_id = ''
        self.__boundary_error = False
        self.__unmatched_error = False #keeps track if it is contains a match or not

        self.__search_id = ''
        self.__total_errors = 0






        #Common to NCBi
        self.__locus_tag_missing = False
        self.__locus_tag_typo = False #This can only be computed at the genome level since it uses the phage name
        self.__description_field_error = False
        self.__total_errors = 0












    # Define all attribute setters:



    #Common to all
    def set_left_boundary(self,value):
        self.__left_boundary = value
    def set_right_boundary(self,value):
        self.__right_boundary = value
    def set_strand(self,value):
        self.__strand = parse_strand(value)
    def set_translation(self,value):
        self.__translation = value.upper()
        self.__translation_length = len(self.__translation)
    def set_type_id(self,value):
        self.__type_id = value
    def compute_amino_acid_errors(self,protein_alphabet_set):
        amino_acid_set = set(self.__translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set
        if len(amino_acid_error_set) > 0:
            self.__amino_acid_errors = True
    def set_start_end_strand_id(self):
        #Create a tuple of feature location data.
        #For start and end of feature, it doesn't matter whether the feature is complex with a translational
        #frameshift or not. Retrieving the "start" and "end" attributes return the very beginning and end of
        #the feature, disregarding the inner "join" coordinates.
        self.__start_end_strand_id = (str(self.__left_boundary),str(self.__right_boundary),self.__strand)

        #Since this id matched genes with different start sites,
        #the strand impacts whether the left or right boundary is used
        if self.__strand == 'forward':
            self.__end_strand_id = (str(self.__right_boundary),self.__strand)
        elif self.__strand == 'reverse':
            self.__end_strand_id = (str(self.__left_boundary),self.__strand)
        else:
            pass
    def set_unmatched_error(self):
        self.__unmatched_error = True
    def compute_boundary_error(self):
        #Check if start and end coordinates are fuzzy
        if not (str(self.__left_boundary).isdigit() and str(self.__right_boundary).isdigit()):
            self.__boundary_error = True






    #Common to Phamerator
    def set_phage_id(self,value):
        self.__phage_id = value
        self.__search_id = remove_draft_suffix(self.__phage_id)
    def set_gene_id(self,value):
        self.__gene_id = value
    def set_gene_name(self,name):
        self.__gene_name = name
    def set_notes(self,value1,value2):
        self.__notes = value1
        self.__search_notes = value2
    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self.__total_errors += 1
        if self.get_boundary_error():
            self.__total_errors += 1
        if self.get_unmatched_error():
            self.__total_errors += 1




    #Common to NCBI
    def set_locus_tag(self,value):
        self.__locus_tag = value
        if self.__locus_tag == '':
            self.__locus_tag_missing = True
    def set_gene_number(self,value):
        self.__gene_number = value
    def set_product_description(self,value1,value2):
        self.__product_description = value1
        self.__search_product_description = value2
    def set_function_description(self,value1,value2):
        self.__function_description = value1
        self.__search_function_description = value2
    def set_note_description(self,value1,value2):
        self.__note_description = value1
        self.__search_note_description = value2
    def set_locus_tag_typo(self):
        self.__locus_tag_typo = True

    def compute_description_error(self):

        #If the product description is empty or generic, and the function or note descriptions are not, there is an error
        if self.__search_product_description == '' and \
            (self.__search_function_description != '' or \
            self.__search_note_description != ''):

            self.__description_field_error = True

    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self.__total_errors += 1
        if self.get_boundary_error():
            self.__total_errors += 1
        if self.__description_field_error:
            self.__total_errors += 1
        if self.__locus_tag_missing:
            self.__total_errors += 1
        if self.__locus_tag_typo:
            self.__total_errors += 1
        if self.get_unmatched_error():
            self.__total_errors += 1











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
        return self.__amino_acid_errors
    def get_translation(self):
        return self.__translation
    def get_translation_length(self):
        return self.__translation_length
    def get_start_end_strand_id(self):
        return self.__start_end_strand_id
    def get_end_strand_id(self):
        return self.__end_strand_id
    def get_unmatched_error(self):
        return self.__unmatched_error
    def get_boundary_error(self):
        return self.__boundary_error


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
        return self.__search_id
    def get_total_errors(self):
        return self.__total_errors





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
        return self.__locus_tag_missing
    def get_locus_tag_typo(self):
        return self.__locus_tag_typo
    def get_description_field_error(self):
        return self.__description_field_error
    def get_total_errors(self):
        return self.__total_errors
