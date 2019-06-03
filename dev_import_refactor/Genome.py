"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""
import FunctionsSimple
import Eval
from datetime import datetime


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



        # Common to GenBank-formatted flat file (NCBI) records
        self.record_name = ""
        self.record_id = "" # TODO might not need this anymore
        self.record_accession = ""
        self.record_description = ""
        self.record_source = ""
        self.record_organism = ""
        self.source_feature_organism = ""
        self.source_feature_host = ""
        self.source_feature_lab_host = ""
        self.record_authors = ""
        self.parsed_record = [] # Holds parsed flat file record.
        self.filename = "" # The file name from which the record is derived








        #Computed datafields


        # Computed datafields: common to all genomes
        self.search_name = '' # Lowercase phage_name void of "_draft"
        self._length = 0 # Size of the nucleotide sequence
        self.evaluations = [] # List of warnings and errors about the genome


        #Common to annotated genomes
        self.cds_features = [] # List of all parsed CDS features
        self._cds_features_tally = 0
        self._cds_processed_primary_descriptions_tally = 0


        self.trna_features = []
        self._trna_features_tally = 0


        self.source_features = []
        self._source_features_tally = 0



        #Common to Phamerator
        self.search_id = '' #Lowercase phage_id void of "_draft"


        # Computed datafields: common to flat file (NCBI) records
        self.search_filename = "" # Lowercase file name

        self._cds_processed_product_descriptions_tally = 0
        self._cds_processed_function_descriptions_tally = 0
        self._cds_processed_note_descriptions_tally = 0




        # TODO create functions to compute these values. Not sure if they
        # need distinct attributes instead of simply EvalResult objects.
        self.record_header_fields_phage_name_error = False
        self.record_header_fields_host_error = False

        self.cds_features_unique_ids = set()
        self.cds_features_duplicate_ids = set()









    def set_evaluation(self, type, message1 = None, message2 = None):
        """Creates an EvalResult object and adds it to the list of all
        evaluations."""

        if type == "warning":
            eval_object = Eval.construct_warning(message1, message2)
        elif type == "error":
            eval_object = Eval.construct_error(message1)
        else:
            eval_object = Eval.EvalResult()
        self.evaluations.append(eval_object)


    def set_filename(self, value):
        """Set the filename. Discard the path and file extension."""

        split_filepath = value.split('/')
        filename = split_filepath[-1]
        filename = filename.split('.')[0]
        self.filename = filename
        self.search_filename = filename.lower()


    # Common to Phamerator

    def set_phage_id(self, value):
        """Set the phage_id and search_id."""

        self.phage_id = FunctionsSimple.remove_draft_suffix(value)
        self.search_id = self.phage_id.lower()




    # TODO create method to choose how to set phage name from NCBI record.
    # Can choose between Organism, Description, etc.




    # TODO I don't think I need this anymore.
    # Original script removed draft suffix.
    # Not sure if "_draft" should remain or not.
    # def set_phage_name(self,value):
    #     self.phage_name = value
    #
    #     self.search_name = \
    #           FunctionsSimple.remove_draft_suffix(self.phage_name).lower()


    def set_host(self, value):
        """Set the host genus and discard the species."""

        value = value.strip()
        if value != "none":
            self.host = value.split(" ")[0]
        else:
            self.host = ""

    def set_sequence(self, value):
        """Set the nucleotide sequence and compute the length."""

        self.sequence = value.upper() #TODO should this be biopython object?
        self._length = len(self.sequence)


    def set_accession(self, value, strategy):
        """Set the accession. The Accession field in Phamerator defaults to "".
        Some flat file accessions have the version number suffix, so discard
        the version number. Strategy indicates how an empty accession value
        should be stored:
        'empty_string' = empty string.
        'none_string' = 'none' string.
        'none_object' = None object."""

        if value is None:
            value = ""

        value = value.strip()
        value = value.split(".")[0]

        if value == "":

            if strategy == "empty_string":
                self.accession = ""
            elif strategy == "none_string":
                self.accession = "none"
            elif strategy == "none_object":
                self.accession = None
            else:
                pass

        else:
            self.accession = value


    def set_cds_features(self, value):
        """Set and tally the CDS features."""

        self.cds_features = value # Should be a list
        self._cds_features_tally = len(self.cds_features)


    def set_trna_features(self, value):
        """Set and tally the tRNA features."""

        self.trna_features = value # Should be a list
        self._trna_features_tally = len(self.trna_features)


    def set_source_features(self, value):
        """Set and tally the source features."""

        self.source_features = value # Should be a list
        self._source_features_tally = len(self.source_features)





    #Common to Phamerator and phagesdb

    def set_cluster(self,value):
        """Set the cluster, and modify singleton if needed."""

        if value.lower() == "singleton":
            self.cluster = value.lower()
        else:
            self.cluster = value

    def set_cluster_subcluster(self, cluster, subcluster):
        """Set the combined Cluster-Subcluster field."""

        if subcluster == "":
            if cluster is None:
                self.cluster_subcluster = "singleton"
            else:
                self.cluster_subcluster = cluster
        else:
            self.cluster_subcluster = subcluster



    # TODO this method may need to be improved. May be better to structure
    # it more similarly to the general strand format conversion function
    # that utilizes a dictionary. Also, if input value is None, it gets
    # converted first to empty string, which may not be ideal. And if an
    # invalid strategy is chosen, it converts to an empty string, which
    # may not be ideal.
    def set_date_last_modified(self, value, strategy):
        """Set the date_last_modified field. Originally this field in
        Phamerator was not used, so for many genomes this field is empty.
        Strategy indicates how an empty accession value
        should be stored:
        'empty_string' = empty string.
        'empty_datetime_obj' = datetime object with an arbitrary early date.
        'none_object' = None object.
        'none_string' = 'none' string."""

        if value is None:
            value = ""

        if value == "":

            if strategy == "empty_string":
                self.date_last_modified = ""
            elif strategy == "empty_datetime_obj":
                self.date_last_modified = \
                    datetime.strptime('1/1/1900', '%m/%d/%Y')
            elif strategy == "none_object":
                self.date_last_modified = None
            elif strategy == "none_string":
                self.date_last_modified = "none"
            else:
                self.date_last_modified = ""

        else:
            self.date_last_modified = value


    def tally_descriptions(self):
        """Iterate through all CDS features and determine how many
        non-generic descriptions are present."""

        for cds in self.cds_features:

            if cds.processed_primary_description != "":
                self._cds_processed_primary_descriptions_tally += 1

            if cds.processed_product_description != "":
                self._cds_processed_product_descriptions_tally += 1

            if cds.processed_function_description != "":
                self._cds_processed_function_descriptions_tally += 1

            if cds.processed_note_description != "":
                self._cds_processed_note_descriptions_tally += 1








    # Evaluations.

    def check_nucleotides(self,dna_alphabet_set):
        """Check if all nucleotides in the sequence are expected."""

        nucleotide_set = set(self.sequence)
        nucleotide_error_set = nucleotide_set - dna_alphabet_set
        if len(nucleotide_error_set) > 0:

            message = "There are unexpected nucleotides in the genome: " \
                + str(nucleotide_error_set)
            self.set_evaluation("error", message)

    def check_status_accession(self):
        """Compare genome status and accession.
        Now that the AnnotationAuthor field contains authorship data, the
        'gbk' annotation status now reflects an 'unknown' annotation (in
        regards to if it was auto-annotated or manually annotated).
        So for the status-accession error, if the status is 'gbk' ('unkown'),
        there is no reason to assume whether there should be an accession
        or not. Only for 'final' (manually annotated) genomes should
        there be an accession."""

        if self.status == 'final' and self.accession == '':
            message = "The genome is final but does not have an accession. "
            self.set_evaluation("error", message)

    def check_status_descriptions(self):
        """Depending on the status of the genome, CDS features are expected
        to contain or not contain descriptions.
        Draft genomes should not have any descriptions.
        Final genomes should not have any descriptions.
        There are no expectations for other types of genomes."""

        if self.status == 'draft' and \
            self._cds_processed_primary_descriptions_tally > 0:

            message = "The genome is draft status " + \
                        "but contains CDS descriptions."
            self.set_evaluation("error", message)

        elif self.status == 'final' and \
            self._cds_processed_primary_descriptions_tally == 0:

            message = "The genome is final status " + \
                        "but does not contain any CDS descriptions."
            self.set_evaluation("error", message)

        else:
            pass





    # TODO need to implement this function. Pasted from original
    # MatchedGenomes object.
    def check_header_phage_name_typos_xxxxx(self):

        #Compare phage names
        pattern1 = re.compile('^' + ph_genome.get_phage_name() + '$')
        pattern2 = re.compile('^' + ph_genome.get_phage_name())

        if find_name(pattern2,ncbi_genome.get_record_description().split(' ')) == 0 or \
            find_name(pattern1,ncbi_genome.get_record_source().split(' ')) == 0 or \
            find_name(pattern1,ncbi_genome.get_record_organism().split(' ')) == 0 or \
            find_name(pattern1,ncbi_genome.get_source_feature_organism().split(' ')) == 0:

            self.__ncbi_record_header_fields_phage_name_mismatch = True






    # TODO need to implement this function. Pasted from original
    # MatchedGenomes object.
    def check_header_host_name_typos_xxxxx(self):

        #Compare host data
        search_host = ph_genome.get_host()
        if search_host == 'Mycobacterium':
            search_host = search_host[:-3]
        pattern3 = re.compile('^' + search_host)

        if (find_name(pattern3,ncbi_genome.get_record_description().split(' ')) == 0 or \
            find_name(pattern3,ncbi_genome.get_record_source().split(' ')) == 0 or \
            find_name(pattern3,ncbi_genome.get_record_organism().split(' ')) == 0 or \
            find_name(pattern3,ncbi_genome.get_source_feature_organism().split(' ')) == 0) or \
            (ncbi_genome.get_source_feature_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_host().split(' ')) == 0) or \
            (ncbi_genome.get_source_feature_lab_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_lab_host().split(' ')) == 0):

            self.__ncbi_host_mismatch = True






    # TODO need to implement this function. Pasted from original
    # MatchedGenomes object.
    def check_author_xxxxx(self):

        #Check author list for errors
        #For genomes with AnnotationAuthor = 1 (Hatfull), Graham is expected
        #to be an author.
        #For genomes with AnnotationAuthor = 0 (non-Hatfull/Genbank), Graham
        #is NOT expected to be an author.
        pattern5 = re.compile('hatfull')
        search_result = pattern5.search(ncbi_genome.get_record_authors().lower())
        if ph_genome.get_annotation_author() == 1 and search_result == None:
            self.__ph_ncbi_author_error = True
        elif ph_genome.get_annotation_author() == 0 and search_result != None:
            self.__ph_ncbi_author_error = True
        else:
            #Any other combination of phamerator and ncbi author can be skipped
            pass





    # TODO create methods to compute this:
    # TODO unit test
    def compare_feature_identifiers(self):

        # This set will contain all unique feature identifiers created
        # by concatenating the start coordinate,
        # stop coordinate, and strand into a tuple.
        feature_set1 = set()


        # This will will contain all feature identifiers that do not
        # have a unique tuple identifier.
        feature_set1_duplicate_set = set()


        for feature in self.cds_features:

            if feature._left_right_strand_id not in feature_set1:
                feature_set1.add(feature._left_right_strand_id)
            else:
                feature_set1_duplicate_set.add(feature._left_right_strand_id)


        # Remove the duplicate end_strand ids from the main id_set
        feature_set1 = feature_set1 - feature_set1_duplicate_set

        self.cds_features_unique_ids = feature_set1
        self.cds_features_duplicate_ids = feature_set1_duplicate_set










    # TODO decide how to assess whether CDS features contain any issues and
    # how to report that at the genome level.
    # TODO unit test
    # def compute_cds_feature_errors(self):
    #     for cds_feature in self._cds_features:
    #         if cds_feature.get_amino_acid_errors():
    #             self.__cds_features_with_translation_error_tally += 1
    #         if cds_feature.get_boundary_error():
    #             self.__cds_features_boundary_error_tally += 1
