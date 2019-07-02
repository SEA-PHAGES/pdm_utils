"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""
from functions import basic
from classes import Eval
from datetime import datetime
from Bio.SeqUtils import GC
import re


class Genome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields:



        # Common to all genomes
        self.phage_id = "" # Unique identifier. Case sensitive, no "_Draft".
        self.phage_name = "" # Case sensitive and contains "_Draft".
        self.host = ""
        self.sequence = "" # Biopython Seq object
        self.accession = ""
        self.author = "" # TODO do I need this in addition to annotation_author?



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
        self.record_authors = ""
        self.record_date = ""


        self.filename = "" # The file name from which the record is derived
        self.translation_table = ""

        # TODO necessary to retain this?
        self.seqrecord = [] # Holds parsed Biopython SeqRecord object.









        # Computed datafields


        # Computed datafields: common to all genomes
        self.search_id = '' # Lowercase phage_id
        self.search_name = '' # Lowercase phage_name
        self._length = 0 # Size of the nucleotide sequence
        self._gc = 0 # %GC content
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
        self._record_description_host_name = ""
        self._record_source_host_name = ""
        self._record_organism_host_name = ""


        self._cds_processed_product_descriptions_tally = 0
        self._cds_processed_function_descriptions_tally = 0
        self._cds_processed_note_descriptions_tally = 0

        self._cds_unique_start_end_ids = set()
        self._cds_duplicate_start_end_ids = set()
        self._cds_unique_end_strand_ids = set()
        self._cds_duplicate_end_strand_ids = set()




    # def set_evaluation(self, type, message1 = None, message2 = None):
    #     """Creates an EvalResult object and adds it to the list of all
    #     evaluations."""
    #
    #     if type == "warning":
    #         eval_object = Eval.construct_warning(message1, message2)
    #     elif type == "error":
    #         eval_object = Eval.construct_error(message1)
    #     else:
    #         eval_object = Eval.EvalResult()
    #     self.evaluations.append(eval_object)


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

        self.phage_id = basic.edit_draft_suffix(value, "remove")
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














    # TODO I don't think I need this anymore.
    # Original script removed draft suffix.
    # Not sure if "_draft" should remain or not.
    # def set_phage_name(self,value):
    #     self.phage_name = value
    #
    #     self.search_name = \
    #           basic.edit_draft_suffix(self.phage_name, "remove").lower()


    def set_host(self, value, format = "empty_string"):
        """Set the host genus and discard the species."""

        if isinstance(value, str):
            value = value.strip()

        # The host value may need to be split. But don't split until
        # it is determined if the value is a null value.
        value = basic.convert_empty(value, "empty_string")

        if value != "":
            self.host = value.split(" ")[0]

        else:
            self.host = basic.convert_empty(value, format)




    def set_sequence(self, value):
        """Set the nucleotide sequence and compute the length."""

        self.sequence = value.upper() # Biopython Seq object
        self._length = len(self.sequence)

        if self._length > 0:
            self._gc = round(GC(self.sequence), 4)
        else:
            self._gc = -1


    def set_accession(self, value, format = "empty_string"):
        """Set the accession. The Accession field in Phamerator defaults to "".
        Some flat file accessions have the version number suffix, so discard
        the version number."""

        if isinstance(value, str):
            value = value.strip()
            value = value.split(".")[0]

        self.accession = basic.convert_empty(value, format)




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


    def set_subcluster(self, value, format = "empty_string"):
        """Set the subcluster."""

        if isinstance(value, str):
            value = value.strip()
        self.subcluster = basic.convert_empty(value, format)


    def set_cluster_subcluster(self):
        """Set the combined Cluster-Subcluster field using the Cluster
        and Subcluster designations."""

        if (self.subcluster is None or \
            self.subcluster.lower() == "none" or \
            self.subcluster == ""):

            if self.cluster is None:
                self.cluster_subcluster = "singleton"
            else:
                self.cluster_subcluster = self.cluster
        else:
            self.cluster_subcluster = self.subcluster


    def split_cluster_subcluster(self, format = "none_string"):
        """From the combined cluster_subcluster value,
        set the Cluster and Subcluster fields.
        If the combined cluster_subcluster value is None, "none", or "",
        no changes are implemented to the current cluster and subcluster
        attributes."""

        if (self.cluster_subcluster is None or \
            self.cluster_subcluster.lower() == "none" or \
            self.cluster_subcluster == ""):

            pass

        else:

            left, right = basic.split_string(self.cluster_subcluster)

            # If splitting produces only an alphabetic string,
            # set only the cluster.
            if self.cluster_subcluster == left:
                self.cluster = left
                self.subcluster = ""

            # If splitting produces an alphabetic string and a
            # numberic string, set the cluster and subcluster.
            elif self.cluster_subcluster == left + right:
                self.cluster = left
                self.subcluster = left + right

            # If the value can't be split into alphabetic and numeric strings,
            # cluster and subcluster should be set to empty.
            else:
                self.cluster == ""
                self.subcluster == ""

        # Now format the cluster and subclusters if they are empty.
        self.cluster = basic.convert_empty(self.cluster, format)
        self.subcluster = basic.convert_empty(self.subcluster, format)


    def set_date_last_modified(self, value, format = "empty_datetime_obj"):
        """Set the date_last_modified field. Originally this field in
        Phamerator was not used, so for many genomes this field is empty.
        """

        self.date_last_modified = basic.convert_empty(value, format)


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
        # When Biopython SeqIO parses the GenBank record, it automatically
        # determines that it is a DNA sequence. It assigns the Seq object
        # alphabet as IUPACAmbiguousDNA. The alphabet could be coerced
        # to a different alphabet, and then tested using the
        # Bio.Alphabet._verify_alphabet() function. Since this is a
        # function though, it is not clear how stable this function is.
        # Instead, Bio.Alphabet.IUPAC.unambiguous_dna alphabet can be passed
        # to the check_nucleotides method.

        nucleotide_set = set(self.sequence)
        nucleotide_error_set = nucleotide_set - dna_alphabet_set

        if len(nucleotide_error_set) > 0:
            result = \
                "There are unexpected nucleotides in the sequence: " \
                + str(nucleotide_error_set)
            status = "error"

        else:
            result = "There are no unexpected nucleotides in the sequence."
            status = "correct"

        definition = "Check if all nucleotides in the sequence are expected."
        eval = Eval.Eval(id = "GENOME0001", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_status_accession(self):
        """Compare genome status and accession.
        Now that the AnnotationAuthor field contains authorship data, the
        'gbk' annotation status now reflects an 'unknown' annotation (in
        regards to if it was auto-annotated or manually annotated).
        So for the status-accession error, if the status is 'gbk' ('unkown'),
        there is no reason to assume whether there should be an accession
        or not. Only for 'final' (manually annotated) genomes should
        there be an accession."""


        if (self.status == 'final' and self.accession == ''):
            result = "The genome is final but does not have an accession. "
            status = "error"

        else:
            result = "No conflict between status and accession."
            status = "correct"

        definition = "Compare the status and accession."
        eval = Eval.Eval(id = "GENOME0002", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_status_descriptions(self):
        """Depending on the status of the genome, CDS features are expected
        to contain or not contain descriptions.
        Draft genomes should not have any descriptions.
        Final genomes should not have any descriptions.
        There are no expectations for other types of genomes."""


        if self.status == 'draft' and \
            self._cds_processed_primary_descriptions_tally > 0:

            result = "The genome is draft status " + \
                            "but contains CDS descriptions."
            status = "error"

        elif self.status == 'final' and \
            self._cds_processed_primary_descriptions_tally == 0:

            result = "The genome is final status " + \
                            "but does not contain any CDS descriptions."
            status = "error"

        else:
            result = "There is no conflict between status and " + \
                            "CDS descriptions."
            status = "correct"

        definition = "Compare the status and presence of CDS descriptions."
        eval = Eval.Eval(id = "GENOME0003", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_record_description_phage_name(self):
        """Check phage name spelling in the record description field."""


        if self.phage_id != self._record_description_phage_name:
            result = "The phage name in the record_description field " + \
                        "does not match the phage_id."
            status = "error"

        else:
            result = "Record descriptions field contains the phage name."
            status = "correct"

        definition = "Check phage name spelling in the " + \
                            "record description field."
        eval = Eval.Eval(id = "GENOME0004", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_record_source_phage_name(self):
        """Check phage name spelling in the record source field."""


        if self.phage_id != self._record_source_phage_name:
            result = "The phage name in the record_source field " + \
                        "does not match the phage_id."
            status = "error"

        else:
            result = "Record source field contains the phage name."
            status = "correct"

        definition = "Check phage name spelling in the record source field."
        eval = Eval.Eval(id = "GENOME0005", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_record_organism_phage_name(self):
        """Check phage name spelling in the record organism field."""


        if self.phage_id != self._record_organism_phage_name:
            result = "The phage name in the record_organism field " + \
                        "does not match the phage_id."
            status = "error"

        else:
            result = "Record organism field contains the phage name."
            status = "correct"

        definition = "Check phage name spelling in the record organism field."
        eval = Eval.Eval(id = "GENOME0006", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_record_description_host_name(self):
        """Check host name spelling in the record description field."""


        if self.host != self._record_description_host_name:
            result = "The host name in the record_description field " + \
                        "does not match the host."
            status = "error"

        else:
            result = "Record description field contains the host name."
            status = "correct"

        definition = "Check host name spelling in the record description field."
        eval = Eval.Eval(id = "GENOME0007", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_record_source_host_name(self):
        """Check host name spelling in the record source field."""


        if self.host != self._record_source_host_name:
            result = "The host name in the record_source field " + \
                        "does not match the host."
            status = "error"

        else:
            result = "Record source field contains the host name."
            status = "correct"

        definition = "Check host name spelling in the record source field."
        eval = Eval.Eval(id = "GENOME0008", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_record_organism_host_name(self):
        """Check host name spelling in the record organism field."""


        if self.host != self._record_organism_host_name:
            result = "The host name in the record_organism field " + \
                        "does not match the host."
            status = "error"

        else:
            result = "Record organism field contains the host name."
            status = "correct"

        definition = "Check host name spelling in the record organism field."
        eval = Eval.Eval(id = "GENOME0009", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)



    # TODO this could be improved by simply creating a set of authors
    # from the long string of authors, then searching whether a specific
    # author is present or not. It could be improved even further by
    # allowing a list of author names to be provided so that more
    # than one author can be searched for in the long author string
    # parsed from the record.
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

            result = "The expected author is not listed."
            status = "error"

        elif self.annotation_author == 0 and search_result != None:

            result = "The author is not expected to be present."
            status = "error"

        else:
            result = "The authorship is as expected."
            status = "correct"

        definition = "Check authorship."
        eval = Eval.Eval(id = "GENOME0010", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_cds_feature_tally(self):
        """Check to confirm that CDS features have been parsed."""


        if self._cds_features_tally == 0:
            result = "There are no CDS features for this genome."
            status = "error"

        else:
            result = "CDS features were annotated."
            status = "correct"

        definition = "Check if CDS features are annotated."
        eval = Eval.Eval(id = "GENOME0011", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    def check_cds_start_end_ids(self):
        """Check if there are any duplicate start-end coordinates.
        Duplicated start-end coordinates may represent
        unintentional duplicate CDS features."""


        if len(self._cds_duplicate_start_end_ids) > 0:
            result = "There are multiple CDS features with the same " + \
                "start and end coordinates."
            status = "error"
        else:
            result = "All CDS features contain unique start and " + \
                            "end coordinate information."
            status = "correct"

        definition = "Check whether CDS features can be uniquely " + \
                        "identified by their start and end coordinates."
        eval = Eval.Eval(id = "GENOME0012", \
                        definition = definition, \
                        result = result, \
                        status = status)

        self.evaluations.append(eval)


    def check_cds_end_strand_ids(self):
        """Check if there are any duplicate end-strand coordinates.
        Duplicated end-strand coordinates may represent
        unintentional duplicate CDS features with slightly
        different start coordinates."""

        if len(self._cds_duplicate_end_strand_ids) > 0:
            result = "There are multiple CDS features with the same " + \
                        "end coordinate and strand."
            status = "error"
        else:
            result = "All CDS features contain unique strand and " + \
                            "end coordinate information."
            status = "correct"

        definition = "Check whether CDS features can be uniquely " + \
                            "identified by their strand and end coordinate."
        eval = Eval.Eval(id = "GENOME0013", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)








###
