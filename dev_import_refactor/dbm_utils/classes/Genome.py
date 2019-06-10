"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""
from functions import basic
from classes import Eval
from datetime import datetime
import re


class Genome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields:



        # Common to all genomes
        self.phage_id = "" # Unique identifier. Case sensitive, no "_Draft".
        self.phage_name = "" # Case sensitive and contains "_Draft".
        self.host = ""
        self.sequence = "" # TODO should this be a Biopython Seq object?
        self.accession = ""
        self.author = ""



        # Common to Phamerator and PhagesDB
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
        self.search_id = '' # Lowercase phage_id
        self.search_name = '' # Lowercase phage_name
        self._length = 0 # Size of the nucleotide sequence
        self.evaluations = [] # List of warnings and errors about the genome


        #Common to annotated genomes
        self.cds_features = [] # List of all parsed CDS features
        self._cds_features_tally = 0
        self._cds_start_end_ids = []
        self._cds_end_strand_ids = []
        self._cds_processed_primary_descriptions_tally = 0


        self.trna_features = []
        self._trna_features_tally = 0


        self.source_features = []
        self._source_features_tally = 0





        # Computed datafields: common to flat file (NCBI) records
        self.search_filename = "" # Lowercase file name


        self._record_description_phage_name = ""
        self._record_source_phage_name = ""
        self._record_organism_phage_name = ""
        self._source_feature_organism_phage_name = ""

        self._record_description_host_name = ""
        self._record_source_host_name = ""
        self._record_organism_host_name = ""
        self._source_feature_organism_host_name = ""
        self._source_feature_host_host_name = ""
        self._source_feature_lab_host_host_name = ""


        self._cds_processed_product_descriptions_tally = 0
        self._cds_processed_function_descriptions_tally = 0
        self._cds_processed_note_descriptions_tally = 0


        self._cds_unique_start_end_ids = set()
        self._cds_duplicate_start_end_ids = set()
        self._cds_unique_end_strand_ids = set()
        self._cds_duplicate_end_strand_ids = set()




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

        self.phage_id = basic.remove_draft_suffix(value)
        self.search_id = self.phage_id.lower()







    def parse_record_description(self):
        """Retrieve the phage name and host name from the record's
        'description' field."""
        string = self.record_description
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._record_description_phage_name = phage_name
        self._record_description_host_name = host_name

    def parse_record_source(self):
        """Retrieve the phage name and host name from the record's
        'source' field."""
        string = self.record_source
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._record_source_phage_name = phage_name
        self._record_source_host_name = host_name

    def parse_record_organism(self):
        """Retrieve the phage name and host name from the record's
        'organism' field."""
        string = self.record_organism
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._record_organism_phage_name = phage_name
        self._record_organism_host_name = host_name

    def parse_source_feature_organism(self):
        """Retrieve the phage name and host name from the 'organism'
        field in the record's annotated 'source' feature."""
        string = self.source_feature_organism
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._source_feature_organism_phage_name = phage_name
        self._source_feature_organism_host_name = host_name

    def parse_source_feature_host(self):
        """Retrieve the host name from the 'host'
        field in the record's annotated 'source' feature."""
        string = self.source_feature_host
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._source_feature_host_host_name = host_name
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.

    def parse_source_feature_lab_host(self):
        """Retrieve the host name from the 'lab_host'
        field in the record's annotated 'source' feature."""
        string = self.source_feature_lab_host
        phage_name, host_name = \
            basic.parse_names_from_record_field(string)
        self._source_feature_lab_host_host_name = host_name
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.















    # TODO I don't think I need this anymore.
    # Original script removed draft suffix.
    # Not sure if "_draft" should remain or not.
    # def set_phage_name(self,value):
    #     self.phage_name = value
    #
    #     self.search_name = \
    #           basic.remove_draft_suffix(self.phage_name).lower()


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
        """Set the CDS features. Tally the CDS features."""

        self.cds_features = value # Should be a list
        self._cds_features_tally = len(self.cds_features)


    def set_cds_ids(self):
        """
        Create a list of CDS feature identifiers using the
        start and end coordinates.
        Create a list of CDS feature identifiers using the
        end coordinate and strand information.
        """

        start_end_id_list = []
        end_strand_id_list = []

        for cds in self.cds_features:
            start_end_id_list.append(cds._start_end_id)
            end_strand_id_list.append(cds._end_strand_id)

        self._cds_start_end_ids = start_end_id_list
        self._cds_end_strand_ids = end_strand_id_list




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

        value = value.strip()

        if value.lower() == "singleton":
            self.cluster = value.lower()
        else:
            self.cluster = value


    def set_subcluster(self, value, strategy):
        """Set the subcluster. Strategy indicates how an empty subcluster value
        should be stored:
        'empty_string' = empty string.
        'none_string' = 'none' string.
        'none_object' = None object."""

        if (value is None or value.lower() == "none"):
            value = ""

        value = value.strip()

        if value == "":

            if strategy == "empty_string":
                self.subcluster = ""
            elif strategy == "none_string":
                self.subcluster = "none"
            elif strategy == "none_object":
                self.subcluster = None
            else:
                pass

        else:
            self.subcluster = value




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


    def identify_unique_cds_start_end_ids(self):
        """Identify which CDS features contain unique start-end
        coordinates."""

        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_start_end_ids)

        self._cds_unique_start_end_ids = set(unique_id_tuples)
        self._cds_duplicate_start_end_ids = set(duplicate_id_tuples)


    def identify_unique_cds_end_strand_ids(self):
        """Identify which CDS features contain unique end-strand
        coordinates."""

        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_end_strand_ids)

        self._cds_unique_end_strand_ids = set(unique_id_tuples)
        self._cds_duplicate_end_strand_ids = set(duplicate_id_tuples)







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




    # TODO now that the parse_field methods are created,
    # I can improve or refine this method.
    def check_phage_name_typos(self, phage_name):
        """Check phage name spelling in various fields."""

        pattern1 = re.compile("^" + phage_name + "$")
        pattern2 = re.compile("^" + phage_name)

        split_description = self.record_description.split(" ")
        split_source = self.record_source.split(" ")
        split_organism1 = self.record_organism.split(" ")
        split_organism2 = self.source_feature_organism.split(" ")

        if basic.find_expression(pattern2, split_description) == 0 \
            or \
            basic.find_expression(pattern1, split_source) == 0 \
            or \
            basic.find_expression(pattern1, split_organism1) == 0 \
            or \
            basic.find_expression(pattern1, split_organism2) == 0:

            message1 = "There appears to be a phage name discrepancy."
            message2 = "There is a phage name discrepancy."
            self.set_evaluation("warning", message1, message2)


    # TODO now that the parse_field methods are created,
    # I can improve or refine this method.
    def check_host_name_typos(self, host_name):
        """Check host name spelling in various fields."""

        if host_name == 'Mycobacterium':
            host_name = host_name[:-3]
        pattern = re.compile('^' + host_name)

        split_description = self.record_description.split(" ")
        split_source = self.record_source.split(" ")
        split_organism1 = self.record_organism.split(" ")
        split_organism2 = self.source_feature_organism.split(" ")
        split_host1 = self.source_feature_host.split(" ")
        split_host2 = self.source_feature_lab_host.split(" ")


        if (basic.find_expression(pattern,split_description) == 0 \
            or \
            basic.find_expression(pattern,split_source) == 0 \
            or \
            basic.find_expression(pattern,split_organism1) == 0 \
            or \
            basic.find_expression(pattern,split_organism2) == 0) \
            or \
            (self.source_feature_host != "" and \
                basic.find_expression(pattern,split_host1) == 0) \
            or \
            (self.source_feature_lab_host != "" and \
                basic.find_expression(pattern,split_host2) == 0):


            message1 = "There appears to be a host name discrepancy."
            message2 = "There is a host name discrepancy."
            self.set_evaluation("warning", message1, message2)


    def check_author(self):
        """Check author name spelling.
        When AnnotationAuthor is set to 1, it will expect to find the
        provided author in the list of authors.
        When AnnotationAuthor is set to 0, it will expect to NOT find the
        provided author in the list of authors."""

        authors = self.record_authors.lower()
        pattern = re.compile(self.author.lower())
        search_result = pattern.search(authors)

        if self.annotation_author == 1 and search_result == None:

            message1 = "The expected author is not listed."
            self.set_evaluation("warning", message1, message1)

        elif self.annotation_author == 0 and search_result != None:

            message2 = "The author is not expected to be present."
            self.set_evaluation("warning", message2, message2)

        else:
            pass


    def check_cds_start_end_ids(self):
        """Check if there are any duplicate start-end coordinates.
        Duplicated start-end coordinates may represent
        unintentional duplicate CDS features."""

        if len(self._cds_duplicate_start_end_ids) > 0:

            message = "There are multiple CDS features with the same " + \
                "start and end coordinates."
            self.set_evaluation("warning", message, message)


    def check_cds_end_strand_ids(self):
        """Check if there are any duplicate end-strand coordinates.
        Duplicated end-strand coordinates may represent
        unintentional duplicate CDS features with slightly
        different start coordinates."""

        if len(self._cds_duplicate_end_strand_ids) > 0:

            message = "There are multiple CDS features with the same " + \
                "end coordinate and strand."
            self.set_evaluation("warning", message, message)









###
