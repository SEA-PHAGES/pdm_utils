"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""
import functions_general
import Eval


class Genome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields:

        # Common to all genomes
        self.phage_name = "" # Case sensitive and contains "_Draft"
        self.host = ""
        self.sequence = "" # TODO should this be a Biopython Seq object?
        self.accession = ""




        # Common to Phamerator and PhagesDB
        self.phage_id = "" # Case sensitive and does not contain "_Draft"
        self.status = "" # Final, Draft, Gbk version of genome data
        self.cluster = ""
        self.subcluster = ""



        # Common to Phamerator
        self.cluster_subcluster = "" # Combined cluster_subcluster data.
        self.ncbi_update_flag = ""
        self.date_last_modified = ""
        self.annotation_author = "" # 1 (Hatfull), 0 (Genbank)
        self.annotation_qc = "" # 1 (reliable), 0, (not reliable)
        self.retrieve_record = "" # 1 (auto update), 0 (do not auto update)



        # Common to NCBI records
        self.record_name = ""
        self.record_id = "" # Todo might not need this anymore
        self.record_accession = ""
        self.record_description = ""
        self.record_source = ""
        self.record_organism = ""
        self.source_feature_organism = ""
        self.source_feature_host = ""
        self.source_feature_lab_host = ""
        self.record_authors = ""
        self.parsed_record = [] # Holds parsed flat file record. TODO do I need this?
        self.filename = "" # Case sensitive file name from which the record is derived








        #Computed datafields


        # Computed datafields: common to all genomes
        self.search_name = '' # Lowercase phage_name void of "_draft"
        self._length = 0
        self.evaluations = []


        #Common to annotated genomes
        self.cds_features = []
        self._cds_features_tally = 0

        #TODO not sure if I need these.
        # self.__cds_features_with_translation_error_tally = 0
        # self.__cds_features_boundary_error_tally = 0


        self.trna_features = []
        self._trna_features_tally = 0


        self.source_features = []
        self._source_features_tally = 0



        #Common to Phamerator
        self.search_id = '' #Lowercase phage_id void of "_draft"

        self._description_tally = 0



        #TODO I don't think I need these.
        # self.__status_accession_error = False
        # self.__status_description_error = False
        # self.__genes_with_errors_tally = 0 #How many genes in this genome have at least one error?


        # Computed datafields: common to NCBI records
        self.search_filename = "" # Lowercase file name
        self._function_descriptions_tally = 0
        self._product_descriptions_tally = 0
        self._note_descriptions_tally = 0
        self._missing_locus_tags_tally = 0
        self._locus_tag_typos_tally = 0
        self._description_field_error_tally = 0









    # TODO create method to choose how to set phage name from NCBI record.
    # Can choose between Organism, Description, etc.


    # TODO set_CDS_features, set_tRNA_features, set_source_features in which
    # it adds a list of features and automatically computes how many there are.




#HERE


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




    # Define all attribute setters:


    def set_filename(self,value):
        """Set the filename. Discard the path and file extension."""
        split_filepath = value.split('/')
        filename = split_filepath[-1]
        filename = filename.split('.')[0]
        self.filename = filename
        self.search_filename = filename.lower()


    #common to phamerator

    def set_phage_id(self,value):
        self.phage_id = functions_general.remove_draft_suffix(value)
        self.search_id = self.phage_id.lower() #TODO make sure this func is imported


    # TODO I don't think I need this anymore.
    # Original script removed draft suffix.
    # Not sure if "_draft" should remain or not.
    # def set_phage_name(self,value):
    #     self.phage_name = value
    #
    #     self.search_name = functions_general.remove_draft_suffix(self.phage_name).lower()


    def set_host(self,value):
        """Set the host genus and discard the species."""
        value = value.strip()
        if value != "none":
            self.host = value.split(" ")[0]
        else:
            self.host = ""

    def set_sequence(self,value):
        """Set the nucleotide sequence and compute the length."""
        self.sequence = value.upper() #TODO should this be biopython object?
        self._length = len(self.sequence)


    def set_accession(self,value,strategy):
        """Set the accession. The Accession field in Phamerator defaults to "".
        Some flat file accessions have the version number suffix, so discard
        the version number. Strategy indicates how an empty accession value
        should be stored."""

        if value is None:
            value = ""

        value = value.strip()
        value = value.split(".")[0]

        if value == "":

            if strategy == "empty":
                self.accession = ""
            elif strategy == "filled":
                self.accession = "none"
            elif strategy == "null":
                self.accession = None
            else:
                pass

        else:
            self.accession = value


    def set_cds_features(self,value):
        """Set the CDS features and compute the number of CDS features."""
        self._cds_features = value #Should be a list
        self._cds_features_tally = len(self._cds_features)








    # HERE
    # TODO unit test
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


    # TODO unit test
    def set_cluster_subcluster(self,value):
        if value is None:
            self.cluster_subcluster = 'Singleton'
        elif value == 'UNK':
            self.cluster_subcluster = ''
        else:
            self.cluster_subcluster = value





    # TODO unit test
    def set_date_last_modified_empty(self,value):
        if value is None:
            self.date_last_modified = ''
        else:
            self.date_last_modified = value


	# If there is no date in the DateLastModified field,
    # set it to a very early date
    # TODO unit test
    def set_date_last_modified_filled(self,value):
        if value is None:
            self.date_last_modified = datetime.strptime('1/1/1900','%m/%d/%Y')
        else:
            self.date_last_modified = value
























    # Evaluations.

    def check_nucleotides(self,dna_alphabet_set):
        """Check if all nucleotides in the sequence are expected."""
        nucleotide_set = set(self.sequence)
        nucleotide_error_set = nucleotide_set - dna_alphabet_set
        if len(nucleotide_error_set) > 0:

            message = "There are unexpected nucleotides in the genome: " \
                + str(nucleotide_error_set)
            self.set_evaluation("error", message)









    # TODO unit test
    def compute_cds_feature_errors(self):
        for cds_feature in self._cds_features:
            if cds_feature.get_amino_acid_errors():
                self.__cds_features_with_translation_error_tally += 1
            if cds_feature.get_boundary_error():
                self.__cds_features_boundary_error_tally += 1












    # TODO unit test
    def compute_status_description_error(self):
        #Iterate through all CDS features, see if they have descriptions, then compare to the status
        for feature in self.get_cds_features():
            if feature.get_notes() != '':
                self._description_tally += 1
        if self.status == 'draft' and self._description_tally > 0:
            self.__status_description_error = True
        elif self.status == 'final' and self._description_tally == 0:
            self.__status_description_error = True
        else:
            pass




    # TODO unit test
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

    # TODO unit test
    def set_cluster(self,value):

        if value.lower() == "singleton":
            self.cluster = value.lower()
        else:
            self.cluster = value




    #Common to NCBI records

    # TODO unit test
    def compute_ncbi_cds_feature_errors(self):
        for cds_feature in self.get_cds_features():

            #counting descriptions should skip if it is blank or "hypothetical protein"
            if cds_feature.get_search_product_description() != '':
                self._product_descriptions_tally += 1

            if cds_feature.get_search_function_description() != '':
                self._function_descriptions_tally += 1

            if cds_feature.get_search_note_description() != '':
                self._note_descriptions_tally += 1



            if cds_feature.get_locus_tag_missing():
                self._missing_locus_tags_tally += 1
            else:
                pattern4 = re.compile(self.get_search_name())
                search_result = pattern4.search(cds_feature.get_locus_tag().lower())

                if search_result == None:
                    self._locus_tag_typos_tally += 1
                    cds_feature.set_locus_tag_typo() #Sets this attribute to True

            if cds_feature.get_description_field_error():
                self._description_field_error_tally += 1











    # TODO I don't think I need the following getters.
    # # Define all attribute getters:
    #
    # #common to all
    # def get_sequence(self):
    #     return self.sequence #Todo convert to string?
    #
    #
    # def get_nucleotide_errors(self):
    #     return self.__nucleotide_errors
    # def get_cds_features(self):
    #     return self._cds_features
    # def get_cds_features_tally(self):
    #     return self.__cds_features_tally
    # def get_cds_features_with_translation_error_tally(self):
    #     return self.__cds_features_with_translation_error_tally
    # def get_cds_features_boundary_error_tally(self):
    #     return self.__cds_features_boundary_error_tally
    #
    #
    #
    # #common to phamerator
    # def get_status_description_error(self):
    #     return self.__status_description_error
    # def get_status_accession_error(self):
    #     return self.__status_accession_error
    # def get_description_tally(self):
    #     return self.__description_tally
    # def get_genes_with_errors_tally(self):
    #     return self.__genes_with_errors_tally
    #
    #
    #
    #
    #
    # #common to NCBI
    # def get_function_descriptions_tally(self):
    #     return self.__function_descriptions_tally
    # def get_product_descriptions_tally(self):
    #     return self.__product_descriptions_tally
    # def get_note_descriptions_tally(self):
    #     return self.__note_descriptions_tally
    # def get_missing_locus_tags_tally(self):
    #     return self.__missing_locus_tags_tally
    # def get_locus_tag_typos_tally(self):
    #     return self.__locus_tag_typos_tally
    # def get_description_field_error_tally(self):
    #     return self.__description_field_error_tally
