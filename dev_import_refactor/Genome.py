"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""

class Genome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields:

        #Common to all genomes
        self.phage_name = '' #Case sensitive and contains "_Draft"
        self.host = ''
        self.sequence = '' #Todo this should be Biopython Seq object
        self.accession = ''




        #Common to Phamerator and PhagesDB
        self.phage_id = ''
        self.status = '' #Final, Draft, Gbk version of genome data
        self.cluster = ''
        self.subcluster = ''



        #Common to Phamerator
        self.cluster_subcluster = '' #Combined cluster_subcluster data from Phamerator
        self.ncbi_update_flag = ''
        self.date_last_modified = ''
        self.annotation_author = '' #1 (Hatfull), 0 (Genbank)



        #Common to NCBI records
        self.record_name = '' #Todo might not need this anymore
        self.record_id = '' #Todo might not need this anymore
        self.record_accession = ''
        self.record_description = ''
        self.record_source = ''
        self.record_organism = ''
        self.source_feature_organism = ''
        self.source_feature_host = ''
        self.source_feature_lab_host = ''
        self.record_authors = ''








        #Computed datafields


        # Computed datafields: common to all genomes
        self.search_name = '' # Lowercase phage name void of "_Draft"
        self.__length = 0
        self.__nucleotide_errors = False


        #Common to annotated genomes
        self.__cds_features = []
        self.__cds_features_tally = 0
        self.__cds_features_with_translation_error_tally = 0
        self.__cds_features_boundary_error_tally = 0


        #Common to Phamerator
        self.search_id = '' #Lowercase phage ID void of "_Draft" and converted
        self.__status_accession_error = False
        self.__status_description_error = False
        self.__description_tally = 0
        self.__genes_with_errors_tally = 0 #How many genes in this genome have at least one error?


        # Computed datafields: common to NCBI records
        self.__function_descriptions_tally = 0
        self.__product_descriptions_tally = 0
        self.__note_descriptions_tally = 0
        self.__missing_locus_tags_tally = 0
        self.__locus_tag_typos_tally = 0
        self.__description_field_error_tally = 0
















    # Define all attribute setters:

    #common to phamerator
    def set_phage_name(self,value):
        self.phage_name = value
        self.search_name = remove_draft_suffix(self.phage_name) #Todo import remove_draft_suffix func
    def set_host(self,value):
        self.host = value #Todo improve this to parse genus species?



    def set_sequence(self,value):
        self.sequence = value.upper() #Todo should be biopython object?
        self.__length = len(self.sequence)

    def set_accession(self,value):
        if value is None or value.strip() == '':
            self.accession = ''
        else:
            value = value.strip()
            self.accession = value.split('.')[0]

    def compute_nucleotide_errors(self,dna_alphabet_set):
        nucleotide_set = set(self.sequence)
        nucleotide_error_set = nucleotide_set - dna_alphabet_set
        if len(nucleotide_error_set) > 0:
            self.__nucleotide_errors = True


    def set_cds_features(self,value):
        self.__cds_features = value #Should be a list
        self.__cds_features_tally = len(self.__cds_features)

    def compute_cds_feature_errors(self):
        for cds_feature in self.__cds_features:
            if cds_feature.get_amino_acid_errors():
                self.__cds_features_with_translation_error_tally += 1
            if cds_feature.get_boundary_error():
                self.__cds_features_boundary_error_tally += 1




    def set_phage_id(self,value):
        self.phage_id = value
        self.search_id = remove_draft_suffix(self.phage_id) #Todo make sure this func is imported
    def set_status(self,value):
        self.status = value

        #Be sure to first set the accession attribute before the status attribute,
        #else this will throw an error.
        #Now that the AnnotationAuthor field contains authorship data, the
        #'gbk' annotation status now reflects an 'unknown' annotation (in
        #regards to if it was auto-annotated or manually annotated).
        #So for the status-accession error, if the status is 'gbk' ('unkown'),
        #there is no reason to assume whether there should be an accession or not.
        #Only for 'final' (manually annotated) genomes should there be an accession.
        #Old code, using 'gbk':
        # if (self.__status == 'final' or self.__status == 'gbk') and self.get_accession() == '':
        #     self.__status_accession_error = True
        if self.status == 'final' and self.get_accession() == '':
            self.__status_accession_error = True #Todo probably need to change this logic. No need for accessor method



    def set_cluster_subcluster(self,value):
        if value is None:
            self.cluster_subcluster = 'Singleton'
        elif value == 'UNK':
            self.cluster_subcluster = ''
        else:
            self.cluster_subcluster = value

    def set_ncbi_update_flag(self,value):
        self.ncbi_update_flag = value

    def set_date_last_modified(self,value):
        if value is None:
            self.date_last_modified = ''
        else:
            self.date_last_modified = value

    def set_annotation_author(self,value):
        self.annotation_author = value

    def compute_status_description_error(self):
        #Iterate through all CDS features, see if they have descriptions, then compare to the status
        for feature in self.get_cds_features():
            if feature.get_notes() != '':
                self.__description_tally += 1
        if self.status == 'draft' and self.__description_tally > 0:
            self.__status_description_error = True
        elif self.status == 'final' and self.__description_tally == 0:
            self.__status_description_error = True
        else:
            pass


    #Even though this method iterates through the cds features like the compute_status_description_error does,
    #it has to be kept separate, since you need to wait to run this method after all genome and gene matching
    #is completed.
    def compute_genes_with_errors_tally(self):
        for feature in self.get_cds_features():
            #Need to first compute the number of errors per gene
            feature.compute_total_cds_errors()
            if feature.get_total_errors() > 0:
                self.__genes_with_errors_tally += 1


    #Common to Phamerator and phagesdb
    def set_cluster(self,value):
        self.cluster = value #Todo need to handle singletons better?
    def set_subcluster(self,value):
        self.subcluster = value #Todo need to handle non-subclustered data better


    #Common to NCBI records
    def compute_ncbi_cds_feature_errors(self):
        for cds_feature in self.get_cds_features():

            #counting descriptions should skip if it is blank or "hypothetical protein"
            if cds_feature.get_search_product_description() != '':
                self.__product_descriptions_tally += 1

            if cds_feature.get_search_function_description() != '':
                self.__function_descriptions_tally += 1

            if cds_feature.get_search_note_description() != '':
                self.__note_descriptions_tally += 1



            if cds_feature.get_locus_tag_missing():
                self.__missing_locus_tags_tally += 1
            else:
                pattern4 = re.compile(self.get_search_name())
                search_result = pattern4.search(cds_feature.get_locus_tag().lower())

                if search_result == None:
                    self.__locus_tag_typos_tally += 1
                    cds_feature.set_locus_tag_typo() #Sets this attribute to True

            if cds_feature.get_description_field_error():
                self.__description_field_error_tally += 1












    # Define all attribute getters:

    #common to all
    def get_sequence(self):
        return self.sequence #Todo convert to string?


    def get_nucleotide_errors(self):
        return self.__nucleotide_errors
    def get_cds_features(self):
        return self.__cds_features
    def get_cds_features_tally(self):
        return self.__cds_features_tally
    def get_cds_features_with_translation_error_tally(self):
        return self.__cds_features_with_translation_error_tally
    def get_cds_features_boundary_error_tally(self):
        return self.__cds_features_boundary_error_tally



    #common to phamerator
    def get_status_description_error(self):
        return self.__status_description_error
    def get_status_accession_error(self):
        return self.__status_accession_error
    def get_description_tally(self):
        return self.__description_tally
    def get_genes_with_errors_tally(self):
        return self.__genes_with_errors_tally





    #common to NCBI
    def get_function_descriptions_tally(self):
        return self.__function_descriptions_tally
    def get_product_descriptions_tally(self):
        return self.__product_descriptions_tally
    def get_note_descriptions_tally(self):
        return self.__note_descriptions_tally
    def get_missing_locus_tags_tally(self):
        return self.__missing_locus_tags_tally
    def get_locus_tag_typos_tally(self):
        return self.__locus_tag_typos_tally
    def get_description_field_error_tally(self):
        return self.__description_field_error_tally
