"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""

from constants import constants
from functions import basic
from classes import Eval
from datetime import datetime
from Bio.SeqUtils import GC
import re


class Genome:

    # Initialize all attributes:
    def __init__(self):



        # Non-computed datafields:

        # Genome classification
        self.nucleic_acid_type = "" # dsDNA, ssDNA, etc.
        self.order = "" # Caudovirales, Nidovirales, etc.
        self.family = "" # Siphoviridae, Myoviridae, etc.
        self.cluster = "" # A, B, C, Singleton, etc.
        self.subcluster = "" #A1, A2, etc.
        self.phage_id = "" # Unique identifier. Case sensitive, no "_Draft".
        self.phage_name = "" # Case sensitive and contains "_Draft".
        self.sequence = "" # Biopython Seq object



        # Common to all genomes
        self.host_genus = ""
        self.accession = ""
        self.author = "" # TODO do I need this in addition to annotation_author?



        # Common to Phamerator
        self.cluster_subcluster = "" # Combined cluster_subcluster data.
        self.annotation_status = "" # Final, Draft, Gbk version of genome data
        self.date_last_modified = ""
        self.annotation_author = "" # 1 (Hatfull), 0 (Genbank)
        self.annotation_qc = "" # 1 (reliable), 0, (not reliable)
        self.retrieve_record = "" # 1 (auto update), 0 (do not auto update)


        self.translation_table = ""
        self.type = "" # Describes how this genome is used
                        # (e.g. import, phamerator, phagesdb, etc.)



        # Common to GenBank-formatted flat file (NCBI) records

        # TODO necessary to retain this?
        self.record = "" # Holds parsed Biopython SeqRecord object.



        self.record_filename = "" # The file name from which the record is derived
        self.record_name = ""
        self.record_id = "" # TODO might not need this anymore
        self.record_accession = ""
        self.record_description = ""
        self.record_source = ""
        self.record_organism = ""
        self.record_authors = ""
        self.record_date = ""







        # Computed datafields


        # Computed datafields: common to all genomes
        self.search_id = '' # Lowercase phage_id
        self.search_name = '' # Lowercase phage_name
        self._length = 0 # Size of the nucleotide sequence
        self._gc = 0 # %GC content
        self.evaluations = [] # List of warnings and errors about the genome
        self._empty_fields = False


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
        self.search_record_filename = "" # Lowercase file name


        self._record_description_phage_name = ""
        self._record_source_phage_name = ""
        self._record_organism_phage_name = ""
        self._record_description_host_genus = ""
        self._record_source_host_genus = ""
        self._record_organism_host_genus = ""


        self._cds_processed_product_descriptions_tally = 0
        self._cds_processed_function_descriptions_tally = 0
        self._cds_processed_note_descriptions_tally = 0

        self._cds_unique_start_end_ids = set()
        self._cds_duplicate_start_end_ids = set()
        self._cds_unique_end_strand_ids = set()
        self._cds_duplicate_end_strand_ids = set()



    def set_record_filename(self, value):
        """Set the record_filename. Discard the path and file extension."""

        split_filepath = value.split('/')
        record_filename = split_filepath[-1]
        record_filename = record_filename.split('.')[0]
        self.record_filename = record_filename
        self.search_record_filename = record_filename.lower()


    # Common to Phamerator

    def set_phage_id(self, value):
        """Set the phage_id and search_id."""

        self.phage_id = basic.edit_suffix(value, "remove")
        self.search_id = self.phage_id.lower()

    def set_phage_id_from_field(self, value):
        """Set the phage_id from a value parsed from the indicated field."""

        if value == "phage_name":
            self.set_phage_id(self.phage_name)
        elif value == "accession":
            self.set_phage_id(self.accession)
        elif value == "record_name":
            self.set_phage_id(self.record_name)
        elif value == "record_id":
            self.set_phage_id(self.record_id)
        elif value == "record_accession":
            self.set_phage_id(self.record_accession)
        elif value == "record_description":
            self.set_phage_id(self.record_description)
        elif value == "record_source":
            self.set_phage_id(self.record_source)
        elif value == "record_organism":
            self.set_phage_id(self.record_organism)
        elif value == "record_filename":
            self.set_phage_id(self.record_filename)
        elif value == "record_description_phage_name":
            self.set_phage_id(self._record_description_phage_name)
        elif value == "record_source_phage_name":
            self.set_phage_id(self._record_source_phage_name)
        elif value == "record_organism_phage_name":
            self.set_phage_id(self._record_organism_phage_name)
        else:
            self.set_phage_id("")

    def set_host_genus_from_field(self, value):
        """Set the host_genus from a value parsed from the indicated field."""

        if value == "record_name":
            self.set_host(self.record_name)
        elif value =="record_id":
            self.set_host(self.record_id)
        elif value =="record_accession":
            self.set_host(self.record_accession)
        elif value =="record_description":
            self.set_host(self.record_description)
        elif value =="record_source":
            self.set_host(self.record_source)
        elif value =="record_organism":
            self.set_host(self.record_organism)
        elif value =="record_filename":
            self.set_host(self.record_filename)
        elif value =="record_description_host_genus":
            self.set_host(self._record_description_host_genus)
        elif value =="record_source_host_genus":
            self.set_host(self._record_source_host_genus)
        elif value =="record_organism_host_genus":
            self.set_host(self._record_organism_host_genus)
        else:
            self.set_host("")

    def parse_record_description(self):
        """Retrieve the phage name and host_genus name from the record's
        'description' field."""
        string = self.record_description
        phage_name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._record_description_phage_name = phage_name
        self._record_description_host_genus = host_genus

    def parse_record_source(self):
        """Retrieve the phage name and host_genus name from the record's
        'source' field."""
        string = self.record_source
        phage_name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._record_source_phage_name = phage_name
        self._record_source_host_genus = host_genus

    def parse_record_organism(self):
        """Retrieve the phage name and host_genus name from the record's
        'organism' field."""
        string = self.record_organism
        phage_name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._record_organism_phage_name = phage_name
        self._record_organism_host_genus = host_genus


    def set_host(self, value, format = "empty_string"):
        """Set the host_genus and discard the species."""

        if isinstance(value, str):
            value = value.strip()

        # The host_genus value may need to be split. But don't split until
        # it is determined if the value is a null value.
        value = basic.convert_empty(value, "empty_string")

        if value != "":
            self.host_genus = value.split(" ")[0]

        else:
            self.host_genus = basic.convert_empty(value, format)


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


    def set_annotation_author(self,value):
        """Convert author name listed in ticket to binary value if needed."""

        self.annotation_author = basic.convert_author(value)



    def set_annotation_qc(self):
        """Set annotation_qc."""

        # TODO not sure if this is needed.
        # if self.annotation_status == 'final':
        #     self.annotation_qc = 1
        # else:
        #     self.annotation_qc = 0
        pass

    def set_retrieve_record(self):
        """Set retrieve_record."""

        # TODO not sure if this is needed.
        # if self.annotation_author == 1:
        #     self.retrieve_record = 1
        # else:
        #     self.retrieve_record = 0
        pass


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


    def set_empty_fields(self, value):
        """Search all attributes for a specified 'value'.
        The 'expect' parameter indicates whether any attributes are
        expected to contain 'value'."""

        if value in vars(self).values():
            self._empty_fields = True
        else:
            self._empty_fields = False












    # Evaluations.
    def check_phage_id(self, phage_id_set, expect = False):
        """Check that the phage_id is valid."""

        value = basic.check_value_expected_in_set(self.phage_id,
                phage_id_set, expect)
        if value:
            result = "The phage_id is valid."
            status = "correct"
        else:
            result = "The phage_id is not valid."
            status = "error"

        definition = "Check that the phage_id is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_phage_name(self, phage_name_set, expect = False):
        """Check that the phage_name is valid."""

        value = basic.check_value_expected_in_set(self.phage_name,
                phage_name_set, expect)
        if value:
            result = "The phage_name is valid."
            status = "correct"
        else:
            result = "The phage_name is not valid."
            status = "error"

        definition = "Check that the phage_name is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_annotation_status(self,
                    status_set = constants.ANNOTATION_STATUS_SET,
                    expect = False):
        """Check that the annotation status is valid."""

        value = basic.check_value_expected_in_set(self.annotation_status,
                status_set, expect)
        if value:
            result = "The annotation status is valid."
            status = "correct"
        else:
            result = "The annotation status is not valid."
            status = "error"

        definition = "Check that the annotation status is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_host_genus(self, host_set, expect = False):
        """Check that the host_genus is valid."""

        value = basic.check_value_expected_in_set(self.host_genus,
                host_set, expect)
        if value:
            result = "The host_genus is valid."
            status = "correct"
        else:
            result = "The host_genus is not valid."
            status = "error"

        definition = "Check that the host_genus is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_cluster(self, cluster_set, expect = False):
        """Check that the cluster is valid."""

        value = basic.check_value_expected_in_set(self.cluster,
                cluster_set, expect)
        if value:
            result = "The cluster is valid."
            status = "correct"
        else:
            result = "The cluster is not valid."
            status = "error"

        definition = "Check that the cluster is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster(self, subcluster_set, expect = False):
        """Check that the subcluster is valid."""


        value = basic.check_value_expected_in_set(self.subcluster,
                subcluster_set, expect)
        if value:
            result = "The subcluster is valid."
            status = "correct"
        else:
            result = "The subcluster is not valid."
            status = "error"

        definition = "Check that the subcluster is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    # TODO this may need to be updated to account for Seq alphabets better.
    def check_sequence(self, seq_set, expect = False):
        """Check that the phage_name is valid."""

        value = basic.check_value_expected_in_set(self.sequence,
                seq_set, expect)
        if value:
            result = "The sequence is valid."
            status = "correct"
        else:
            result = "The sequence is not valid."
            status = "error"

        definition = "Check that the sequence is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, accession_set, expect = False):
        """Check that the accession is valid."""

        value = basic.check_value_expected_in_set(self.accession,
                accession_set, expect)
        if value:
            result = "The accession is valid."
            status = "correct"
        else:
            result = "The accession is not valid."
            status = "error"

        definition = "Check that the accession is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    # TODO implement.
    # TODO unit test.
    def check_date_last_modified(self):
        """Check that the date_last_modified is valid."""

        # definition = "Check that the date_last_modified is valid."
        # eval = Eval.Eval("GENOME", definition, result, status)
        # self.evaluations.append(eval)
        pass


    # TODO improve this to accept an author list?
    def check_annotation_author(self):
        """Check that the annotation_author is valid."""

        if (self.annotation_author == 0 or self.annotation_author == 1):
            result = "The annotation_author is valid."
            status = "correct"
        else:
            result = "The annotation_author is not valid."
            status = "error"

        definition = "Check that the annotation_author is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_qc(self):
        """Check that the annotation_qc is valid."""

        if (self.annotation_qc == 0 or self.annotation_qc == 1):
            result = "The annotation_qc is valid."
            status = "correct"
        else:
            result = "The annotation_qc is not valid."
            status = "error"

        definition = "Check that the annotation_qc is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_retrieve_record(self):
        """Check that the retrieve_record is valid."""

        if (self.retrieve_record == 0 or self.retrieve_record == 1):
            result = "The retrieve_record is valid."
            status = "correct"
        else:
            result = "The retrieve_record is not valid."
            status = "error"

        definition = "Check that the retrieve_record is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_record_filename(self, filename_set, expect = False):
        """Check that the record_filename is valid."""

        value = basic.check_value_expected_in_set(self.record_filename,
                filename_set, expect)
        if value:
            result = "The record_filename is valid."
            status = "correct"
        else:
            result = "The record_filename is not valid."
            status = "error"

        definition = "Check that the record_filename is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    # TODO is this needed?
    def check_record(self, record_set, expect = False):
        """Check that the record is valid."""

        value = basic.check_value_expected_in_set(self.record,
                record_set, expect)
        if value:
            result = "The record is valid."
            status = "correct"
        else:
            result = "The record is not valid."
            status = "error"

        definition = "Check that the record is valid."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)




    def check_subcluster_structure(self):
        """Check whether the cluster field is structured appropriately."""

        if self.subcluster != "none":
            left, right = basic.split_string(self.subcluster)

            if (left.isalpha() == False or right.isdigit() == False):
                result = "Subcluster is not structured correctly."
                status = "error"
            else:
                result = "Subcluster is structured correctly."
                status = "correct"

        else:
            result = "Subcluster is empty."
            status = "not_evaluated"

        definition = "Check if subcluster field is structured correctly."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_cluster_structure(self):
        """Check whether the cluster field is structured appropriately."""

        if self.cluster != "none":
            left, right = basic.split_string(self.cluster)

            if (right != "" or left.isalpha() == False):
                result = "Cluster is not structured correctly."
                status = "error"
            else:
                result = "Cluster is structured correctly."
                status = "correct"
        else:
            result = "Cluster is empty."
            status = "not_evaluated"


        definition = "Check if cluster field is structured correctly."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def compare_cluster_subcluster_structure(self):
        """Check whether the cluster and subcluster fields are
        compatible."""

        output = basic.compare_cluster_subcluster(self.cluster, self.subcluster)
        if not output:
            result = "Cluster and Subcluster designations are not compatible."
            status = "error"
        else:
            result = "Cluster and Subcluster designations are compatible."
            status = "correct"

        definition = "Check for compatibility between Cluster and Subcluster."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


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
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_status_accession(self):
        """Compare genome annotation_status and accession.
        Now that the AnnotationAuthor field contains authorship data, the
        'gbk' annotation status now reflects an 'unknown' annotation (in
        regards to if it was auto-annotated or manually annotated).
        So for the annotation_status-accession error,
        if the annotation_status is 'gbk' ('unkown'),
        there is no reason to assume whether there should be an accession
        or not. Only for 'final' (manually annotated) genomes should
        there be an accession."""


        if (self.annotation_status == 'final' and self.accession == ''):
            result = "The genome is final but does not have an accession. "
            status = "error"

        else:
            result = "No conflict between annotation_status and accession."
            status = "correct"

        definition = "Compare the annotation_status and accession."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_status_descriptions(self):
        """Depending on the annotation_status of the genome,
        CDS features are expected to contain or not contain descriptions.
        Draft genomes should not have any descriptions.
        Final genomes should not have any descriptions.
        There are no expectations for other types of genomes."""


        if self.annotation_status == 'draft' and \
            self._cds_processed_primary_descriptions_tally > 0:

            result = "The genome is draft status " + \
                            "but contains CDS descriptions."
            status = "error"

        elif self.annotation_status == 'final' and \
            self._cds_processed_primary_descriptions_tally == 0:

            result = "The genome is final status " + \
                            "but does not contain any CDS descriptions."
            status = "error"

        else:
            result = "There is no conflict between status and " + \
                            "CDS descriptions."
            status = "correct"

        definition = "Compare the status and presence of CDS descriptions."
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_record_description_host_genus(self):
        """Check host_genus name spelling in the record description field."""


        if self.host_genus != self._record_description_host_genus:
            result = "The host_genus name in the record_description field " + \
                        "does not match the host_genus."
            status = "error"

        else:
            result = "Record description field contains the host_genus name."
            status = "correct"

        definition = "Check host_genus name spelling in the record description field."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_record_source_host_genus(self):
        """Check host_genus name spelling in the record source field."""


        if self.host_genus != self._record_source_host_genus:
            result = "The host_genus name in the record_source field " + \
                        "does not match the host_genus."
            status = "error"

        else:
            result = "Record source field contains the host_genus name."
            status = "correct"

        definition = "Check host_genus name spelling in the record source field."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)

    def check_record_organism_host_genus(self):
        """Check host_genus name spelling in the record organism field."""


        if self.host_genus != self._record_organism_host_genus:
            result = "The host_genus name in the record_organism field " + \
                        "does not match the host_genus."
            status = "error"

        else:
            result = "Record organism field contains the host_genus name."
            status = "correct"

        definition = "Check host_genus name spelling in the record organism field."
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
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
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)


    def check_empty_fields(self, expect = False):
        """Check if there are any fields that are not populated as expected."""

        if self._empty_fields:
            if expect:
                result = "All fields are populated."
                status = "correct"
            else:
                result = "Some fields are not populated."
                status = "error"
        else:
            if not expect:
                result = "All fields are populated."
                status = "correct"
            else:
                result = "Some fields are not populated."
                status = "error"

        definition = "Check if there are any fields that are set to %s."
        eval = Eval.Eval("GENOME", definition, result, status)
        self.evaluations.append(eval)




###
