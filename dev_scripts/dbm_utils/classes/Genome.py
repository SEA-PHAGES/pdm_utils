"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""

from functions import basic
from classes import Eval
from datetime import datetime
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import re
from operator import attrgetter
from classes import Cds, Trna


class Genome:
    """Class to hold data about a phage genome."""

    def __init__(self):

        # The following attributes are common to any genome.
        self.nucleic_acid_type = "" # dsDNA, ssDNA, etc.
        self.order = "" # Caudovirales, Nidovirales, etc.
        self.family = "" # Siphoviridae, Myoviridae, etc.
        self.cluster = "" # A, B, C, Singleton, etc.
        self.subcluster = "" #A1, A2, etc.
        self.id = "" # Unique identifier. Case sensitive, no "_Draft".
        self.name = "" # Case sensitive and contains "_Draft".
        self.seq = "" # Biopython Seq object
        self._length = 0 # Size of the nucleotide sequence
        self._gc = 0 # %GC content
        self.host_genus = ""
        self.accession = ""
        self.lifestyle = "" # E.g. temperate, obligately lytic, unknown, etc.
        self.translation_table = ""

        # The following attributes are common to PhameratorDB.
        self.cluster_subcluster = "" # Combined cluster/subcluster data.
        self.annotation_status = "" # Final, Draft, Gbk version of genome data
        self.annotation_author = "" # 1 (can be changed), 0 (can not be changed)
        self.annotation_qc = "" # 1 (reliable), 0, (not reliable)
        self.retrieve_record = "" # 1 (auto update), 0 (do not auto update)
        self.date = "" # Used for the DateLastModified field.

        # The following attributes are common to
        # GenBank-formatted flat file records.
        self.description = ""
        self._description_name = ""
        self._description_host_genus = ""
        self.source = ""
        self._source_name = ""
        self._source_host_genus = ""
        self.organism = ""
        self._organism_name = ""
        self._organism_host_genus = ""
        self.authors = "" # Compiled list of all authors named in the record.

        # The following attributes are computed datafields that are
        # common to annotated genomes.
        self.cds_features = [] # List of all parsed CDS features
        self._cds_features_tally = 0
        self._cds_start_end_ids = []
        self._cds_end_strand_ids = []
        self._cds_processed_descriptions_tally = 0
        self._cds_processed_products_tally = 0
        self._cds_processed_functions_tally = 0
        self._cds_processed_notes_tally = 0
        self._cds_unique_start_end_ids = set() # TODO still in development.
        self._cds_duplicate_start_end_ids = set() # TODO still in development.
        self._cds_unique_end_strand_ids = set() # TODO still in development.
        self._cds_duplicate_end_strand_ids = set() # TODO still in development.
        self.trna_features = []
        self._trna_features_tally = 0
        self.tmrna_features = []
        self._tmrna_features_tally = 0
        self.source_features = []
        self._source_features_tally = 0

        # The following attributes are usefule for processing data
        # from various data sources.
        self.filename = "" # The file name from which the data is derived
        self.type = "" # Identifier to describes how this genome is used
                       # (e.g. import, phamerator, phagesdb, etc.)
        self.evaluations = [] # List of warnings and errors about the data
        self._value_flag = False


    def set_filename(self, value):
        """Set the filename. Discard the path and file extension."""
        split_filepath = value.split('/')
        filename = split_filepath[-1]
        filename = filename.split('.')[0]
        self.filename = filename


    def set_id(self, value=None, attribute=None):
        """Set the id from either an input value or an indicated attribute."""
        if value is None:
            if attribute == "name":
                self.id = basic.edit_suffix(self.name, "remove")
            elif attribute == "accession":
                self.id = basic.edit_suffix(self.accession, "remove")
            elif attribute == "description":
                self.id = basic.edit_suffix(self.description, "remove")
            elif attribute == "source":
                self.id = basic.edit_suffix(self.source, "remove")
            elif attribute == "organism":
                self.id = basic.edit_suffix(self.organism, "remove")
            elif attribute == "filename":
                self.id = basic.edit_suffix(self.filename, "remove")
            elif attribute == "description_name":
                self.id = basic.edit_suffix(self._description_name, "remove")
            elif attribute == "source_name":
                self.id = basic.edit_suffix(self._source_name, "remove")
            elif attribute == "organism_name":
                self.id = basic.edit_suffix(self._organism_name, "remove")
            else:
                self.id = basic.edit_suffix("", "remove")
        else:
            self.id = basic.edit_suffix(value, "remove")


    def set_host_genus(self, value=None, attribute=None, format="empty_string"):
        """Set the host_genus from a value parsed from the indicated attribute.

        The input data is split into multiple parts, and the first word is
        used to set host_genus.

        :param value: the host genus of the phage genome
        :type value: str
        :param attribute: the name of the genome attribute from which the
            host_genus attribute will be set
        :type attribute: str
        :param format: the default format if the input is an empty/null value.
        :type format: str
        """

        if value is None:
            if attribute == "name":
                value = self.name
            elif attribute == "description":
                value = self.description
            elif attribute == "source":
                value = self.source
            elif attribute == "organism":
                value = self.organism
            elif attribute == "filename":
                value = self.filename
            elif attribute == "description_host_genus":
                value = self._description_host_genus
            elif attribute == "source_host_genus":
                value = self._source_host_genus
            elif attribute == "organism_host_genus":
                value = self._organism_host_genus
            else:
                value = ""

        if isinstance(value, str):
            value = value.strip()

        # The host_genus value may need to be split. But don't split until
        # it is determined if the value is a null value.
        value = basic.convert_empty(value, "empty_string")
        if value != "":
            self.host_genus = value.split(" ")[0]
        else:
            self.host_genus = basic.convert_empty(value, format)


    def parse_description(self):
        """Retrieve the name and host_genus from the 'description' attribute."""
        string = self.description
        name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._description_name = name
        self._description_host_genus = host_genus


    def parse_source(self):
        """Retrieve the name and host_genus from the 'source' attribute."""
        string = self.source
        name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._source_name = name
        self._source_host_genus = host_genus


    def parse_organism(self):
        """Retrieve the name and host_genus from the 'organism' attribute."""
        string = self.organism
        name, host_genus = \
            basic.parse_names_from_record_field(string)
        self._organism_name = name
        self._organism_host_genus = host_genus


    def set_sequence(self, value):
        """Set the nucleotide sequence and compute the length.

        This method coerces sequences into a Biopython Seq object.
        """
        if not isinstance(value, Seq):
            self.seq = Seq(value).upper()
        else:
            self.seq = value.upper()
        self._length = len(self.seq)
        if self._length > 0:
            self._gc = round(GC(self.seq), 4)
        else:
            self._gc = -1


    def set_accession(self, value, format="empty_string"):
        """Set the accession.

        The Accession field in Phamerator defaults to "".
        Some flat file accessions have the version number suffix, so discard
        the version number.
        """
        if isinstance(value, str):
            value = value.strip()
            value = value.split(".")[0]
        self.accession = basic.convert_empty(value, format)


    def set_cds_features(self, value):
        """Set and tally the CDS features."""
        self.cds_features = value # Should be a list.
        self._cds_features_tally = len(self.cds_features)


    def set_cds_id_list(self):
        """Creates lists of CDS feature identifiers.

        The first identifier is derived from the start and end coordinates.
        The second identifier is derived from the end coordinate and strand.
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


    def set_cluster(self, value):
        """Set the cluster and modify singleton if needed."""
        if isinstance(value, str):
            value = value.strip()
            if value.lower() == "singleton":
                self.cluster = value.lower()
            else:
                self.cluster = value
        if value is None:
            self.cluster = "singleton"


    def set_subcluster(self, value, format="empty_string"):
        """Set the subcluster."""
        if isinstance(value, str):
            value = value.strip()
        self.subcluster = basic.convert_empty(value, format)


    def set_cluster_subcluster(self, value="internal"):
        """Set the combined Cluster-Subcluster attribute.

        If the value is set to 'internal', it is determined from the Cluster
        and Subcluster designations. Otherwise, the value is directly
        used to populate this attribute.
        """
        if value is "internal":
            if (self.subcluster is None or \
                self.subcluster.lower() == "none" or \
                self.subcluster == ""):
                if self.cluster is None:
                    self.cluster_subcluster = "singleton"
                else:
                    self.cluster_subcluster = self.cluster
            else:
                self.cluster_subcluster = self.subcluster
        elif value is None:
            self.cluster_subcluster = "singleton"
        elif value.lower() == "singleton":
            self.cluster_subcluster = value.lower()
        else:
            self.cluster_subcluster = value


    def split_cluster_subcluster(self, format="none_string"):
        """Split the combined cluster_subcluster data.

        Sets the 'cluster' and 'subcluster' attributes from the
        'cluster_subcluster' attribute. If the combined 'cluster_subcluster'
        attribute is None, "none", or "", no changes are implemented
        to the current cluster and subcluster attributes.
        """

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


    def set_date(self, value, format="empty_datetime_obj"):
        """Set the date attribute."""
        self.date = basic.convert_empty(value, format)



    # TODO this may no longer be needed.
    def set_annotation_author(self,value):
        """Convert author name listed in ticket to binary value if needed."""
        # self.annotation_author = basic.convert_author(value)
        self.annotation_author = value


    # TODO implement.
    # TODO unit test.
    def set_annotation_qc(self):
        """Set annotation_qc."""

        # TODO not sure if this is needed.
        # if self.annotation_status == 'final':
        #     self.annotation_qc = 1
        # else:
        #     self.annotation_qc = 0
        pass

    # TODO implement.
    # TODO unit test.
    def set_retrieve_record(self):
        """Set retrieve_record."""

        # TODO not sure if this is needed.
        # if self.annotation_author == 1:
        #     self.retrieve_record = 1
        # else:
        #     self.retrieve_record = 0
        pass


    def tally_descriptions(self):
        """Tally the non-generic CDS descriptions."""
        for cds in self.cds_features:
            if cds.processed_description != "":
                self._cds_processed_descriptions_tally += 1
            if cds.processed_product != "":
                self._cds_processed_products_tally += 1
            if cds.processed_function != "":
                self._cds_processed_functions_tally += 1
            if cds.processed_note != "":
                self._cds_processed_notes_tally += 1


    def set_unique_cds_start_end_ids(self):
        """Identify CDS features contain unique start-end coordinates."""
        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_start_end_ids)
        self._cds_unique_start_end_ids = set(unique_id_tuples)
        self._cds_duplicate_start_end_ids = set(duplicate_id_tuples)


    def set_unique_cds_end_strand_ids(self):
        """Identify CDS features contain unique end-strand coordinates."""
        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_end_strand_ids)
        self._cds_unique_end_strand_ids = set(unique_id_tuples)
        self._cds_duplicate_end_strand_ids = set(duplicate_id_tuples)




    def set_value_flag(self, value):
        """Sets the flag if any attributes contain the specified 'value'."""
        if value in vars(self).values():
            self._value_flag = True
        else:
            self._value_flag = False


    def set_feature_ids(self, use_type=False, use_cds=False,
                        use_trna=False, use_tmrna=False):
        """Sets the id of each feature.

        Lists of features can be added to this method. The method assumes
        that all elements in all lists contain 'id', 'left', and 'right'
        attributes. This feature attribute is processed within
        the Genome object because and not within the feature itself since
        the method sorts all features and generates systematic IDs based on
        feature order in the genome.
        """

        # Both coordinates are used to control the order of features
        # that may share one, but not both, coordinates (e.g. tail
        # assembly chaperone).
        list_to_sort = []
        if use_cds:
            list_to_sort.extend(self.cds_features)
        if use_trna:
            list_to_sort.extend(self.trna_features)

        # TODO unit test after tmrna_features are implemented.
        if use_tmrna:
            list_to_sort.extend(self.tmrna_features)

        sorted_list = sorted(list_to_sort, key=attrgetter("left", "right"))
        index = 0
        while index < len(sorted_list):
            if use_type:
                if isinstance(sorted_list[index], Cds.Cds):
                    delimiter = "_CDS_"

                # TODO unit test after tRNA class implemented.
                elif isinstance(sorted_list[index], Trna.TrnaFeature):
                    delimiter = "_TRNA_"

                # TODO unit test after tmRNA class implemented.
                elif isinstance(sorted_list[index], Tmrna.Tmrna):
                    delimiter = "_TMRNA_"
                else:
                    delimiter = "_"
            else:
                delimiter = "_"

            sorted_list[index].id = self.id + delimiter + str(index + 1)
            index += 1






    # Evaluations.

    def check_id(self, check_set, expect=False, eval_id=None):
        """Check that the id is valid."""
        value = basic.check_value_expected_in_set(self.id,
                check_set, expect)
        if value:
            result = "The id is valid."
            status = "correct"
        else:
            result = "The id is not valid."
            status = "error"
        definition = "Check that the id is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_name(self, check_set, expect=False, eval_id=None):
        """Check that the name is valid."""
        value = basic.check_value_expected_in_set(self.name,
                check_set, expect)
        if value:
            result = "The name is valid."
            status = "correct"
        else:
            result = "The name is not valid."
            status = "error"
        definition = "Check that the name is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_status(
            self, check_set=set(), expect=False, eval_id=None):
        """Check that the annotation_status is valid."""
        value = basic.check_value_expected_in_set(self.annotation_status,
                check_set, expect)
        if value:
            result = "The annotation status is valid."
            status = "correct"
        else:
            result = "The annotation status is not valid."
            status = "error"
        definition = "Check that the annotation status is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_host_genus(self, check_set, expect=False, eval_id=None):
        """Check that the host_genus is valid."""
        value = basic.check_value_expected_in_set(self.host_genus,
                check_set, expect)
        if value:
            result = "The host_genus is valid."
            status = "correct"
        else:
            result = "The host_genus is not valid."
            status = "error"
        definition = "Check that the host_genus is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cluster(self, check_set, expect=False, eval_id=None):
        """Check that the cluster is valid."""
        value = basic.check_value_expected_in_set(self.cluster,
                check_set, expect)
        if value:
            result = "The cluster is valid."
            status = "correct"
        else:
            result = "The cluster is not valid."
            status = "error"
        definition = "Check that the cluster is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cluster_structure(self, eval_id=None):
        """Check whether the cluster attribute is structured appropriately."""
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
        definition = "Check if cluster attribute is structured correctly."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster(self, check_set, expect=False, eval_id=None):
        """Check that the subcluster is valid."""
        value = basic.check_value_expected_in_set(self.subcluster,
                check_set, expect)
        if value:
            result = "The subcluster is valid."
            status = "correct"
        else:
            result = "The subcluster is not valid."
            status = "error"
        definition = "Check that the subcluster is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster_structure(self, eval_id=None):
        """Check whether the subcluster attribute is structured appropriately."""
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
        definition = "Check if subcluster attribute is structured correctly."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_compatible_cluster_and_subcluster(self, eval_id=None):
        """Check compatibility of cluster and subcluster attributes."""
        output = basic.compare_cluster_subcluster(self.cluster, self.subcluster)
        if not output:
            result = "Cluster and Subcluster designations are not compatible."
            status = "error"
        else:
            result = "Cluster and Subcluster designations are compatible."
            status = "correct"
        definition = "Check for compatibility between Cluster and Subcluster."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_sequence(self, check_set, expect=False, eval_id=None):
        """Check that the sequence is valid."""
        value = basic.check_value_expected_in_set(self.seq,
                check_set, expect)
        if value:
            result = "The sequence is valid."
            status = "correct"
        else:
            result = "The sequence is not valid."
            status = "error"

        definition = "Check that the sequence is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_nucleotides(self, dna_alphabet_set=set(), eval_id=None):
        """Check if all nucleotides in the sequence are expected."""
        # When Biopython SeqIO parses the GenBank record, it automatically
        # determines that it is a DNA sequence. It assigns the Seq object
        # alphabet as IUPACAmbiguousDNA. The alphabet could be coerced
        # to a different alphabet, and then tested using the
        # Bio.Alphabet._verify_alphabet() function. Since this is a private
        # function though, it is not clear how stable/reliable it is.
        # Instead, Bio.Alphabet.IUPAC.unambiguous_dna alphabet can be passed
        # to the check_nucleotides method.

        nucleotide_set = set(self.seq)
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, check_set, expect=False, eval_id=None):
        """Check that the accession is valid."""
        value = basic.check_value_expected_in_set(self.accession,
                check_set, expect)
        if value:
            result = "The accession is valid."
            status = "correct"
        else:
            result = "The accession is not valid."
            status = "error"
        definition = "Check that the accession is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_authors(self, check_set=set(), expect=True, eval_id=None):
        """Check author list.

        The method will compare the list of authors provided in the
        'check_set' parameter with the list of authors stored in the
        'authors' attribute. It will then evaluate whether the presence
        of at least one author matches what is expected, as indicated
        by the 'expect' parameter.
        """
        authors_list = self.authors.lower().split(";")
        authors_list = [x.strip() for x in authors_list]
        authors_set = set(authors_list)
        mutual_authors_set = authors_set & check_set
        if len(mutual_authors_set) == 0:
            if expect:
                result = "The expected authors are not listed."
                status = "error"
            else:
                result = "The authorship is as expected."
                status = "correct"
        else:
            if expect:
                result = "The authorship is as expected."
                status = "correct"
            else:
                result = "The authors are not expected to be present."
                status = "error"
        definition = "Check authorship."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_author(self, check_set=set(), eval_id=None):
        """Check that the annotation_author is valid."""
        if self.annotation_author in check_set:
            result = "The annotation_author is valid."
            status = "correct"
        else:
            result = "The annotation_author is not valid."
            status = "error"
        definition = "Check that the annotation_author is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_qc(self, check_set=set(), eval_id=None):
        """Check that the annotation_qc is valid."""
        if self.annotation_qc in check_set:
            result = "The annotation_qc is valid."
            status = "correct"
        else:
            result = "The annotation_qc is not valid."
            status = "error"
        definition = "Check that the annotation_qc is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_retrieve_record(self, check_set=set(), eval_id=None):
        """Check that the retrieve_record is valid."""
        if self.retrieve_record in check_set:
            result = "The retrieve_record is valid."
            status = "correct"
        else:
            result = "The retrieve_record is not valid."
            status = "error"
        definition = "Check that the retrieve_record is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_filename(self, check_set, expect=False, eval_id=None):
        """Check that the filename is valid."""
        value = basic.check_value_expected_in_set(self.filename,
                check_set, expect)
        if value:
            result = "The filename is valid."
            status = "correct"
        else:
            result = "The filename is not valid."
            status = "error"
        definition = "Check that the filename is valid."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_compatible_status_and_accession(self, eval_id=None):
        """Compare genome annotation_status and accession."""
        # Now that the AnnotationAuthor attribute contains authorship data, the
        # 'gbk' annotation status now reflects an 'unknown' annotation (in
        # regards to if it was auto-annotated or manually annotated).
        # So for the annotation_status-accession error,
        # if the annotation_status is 'gbk' ('unkown'),
        # there is no reason to assume whether there should be an accession
        # or not. Only for 'final' (manually annotated) genomes should
        # there be an accession.

        if (self.annotation_status == 'final' and self.accession == ''):
            result = "The genome is final but does not have an accession. "
            status = "error"
        else:
            result = "No conflict between annotation_status and accession."
            status = "correct"
        definition = "Compare the annotation_status and accession."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_compatible_status_and_descriptions(self, eval_id=None):
        """Compare annotation_status and description tally."""

        # Depending on the annotation_status of the genome,
        # CDS features are expected to contain or not contain descriptions.
        # Draft genomes should not have any descriptions.
        # Final genomes should not have any descriptions.
        # There are no expectations for other types of genomes.

        if (self.annotation_status == 'draft' and \
                self._cds_processed_descriptions_tally > 0):
            result = "The genome is draft status " + \
                     "but contains CDS descriptions."
            status = "error"
        elif (self.annotation_status == 'final' and \
                self._cds_processed_descriptions_tally == 0):
            result = "The genome is final status " + \
                     "but does not contain any CDS descriptions."
            status = "error"
        else:
            result = "There is no conflict between status and " + \
                     "CDS descriptions."
            status = "correct"
        definition = "Compare the status and presence of CDS descriptions."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_description_name(self, eval_id=None):
        """Check genome id spelling in the description attribute."""
        if self.id != self._description_name:
            result = "The name in the description attribute " + \
                     "does not match the genome's id."
            status = "error"
        else:
            result = "Record descriptions attribute contains the genome id."
            status = "correct"
        definition = "Check genome id spelling in the " + \
                     "description attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_source_name(self, eval_id=None):
        """Check genome id spelling in the source attribute."""
        if self.id != self._source_name:
            result = "The name in the source attribute " + \
                     "does not match the genome's id."
            status = "error"
        else:
            result = "The source attribute contains the genome id."
            status = "correct"
        definition = "Check genome id spelling in the source attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_organism_name(self, eval_id=None):
        """Check genome id in the organism attribute."""
        if self.id != self._organism_name:
            result = "The name in the organism attribute " + \
                     "does not match the genome's id."
            status = "error"
        else:
            result = "The organism attribute contains the genome id."
            status = "correct"
        definition = "Check genome id spelling in the organism attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_description_host_genus(self, eval_id=None):
        """Check host_genus spelling in the description attribute."""
        if self.host_genus != self._description_host_genus:
            result = "The host_genus in the description attribute " + \
                     "does not match the genome's host_genus."
            status = "error"
        else:
            result = "The description attribute contains the host_genus."
            status = "correct"
        definition = "Check host_genus spelling in the description attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_source_host_genus(self, eval_id=None):
        """Check host_genus spelling in the source attribute."""
        if self.host_genus != self._source_host_genus:
            result = "The host_genus name in the source attribute " + \
                     "does not match the genome's host_genus."
            status = "error"
        else:
            result = "The source attribute contains the host_genus."
            status = "correct"
        definition = "Check host_genus spelling in the source attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_organism_host_genus(self, eval_id=None):
        """Check host_genus spelling in the organism attribute."""
        if self.host_genus != self._organism_host_genus:
            result = "The host_genus in the organism attribute " + \
                     "does not match the genome's host_genus."
            status = "error"
        else:
            result = "The organism attribute contains the host_genus."
            status = "correct"
        definition = "Check host_genus spelling in the organism attribute."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cds_feature_tally(self, eval_id=None):
        """Check to confirm that CDS features have been parsed."""
        if self._cds_features_tally == 0:
            result = "There are no CDS features for this genome."
            status = "error"
        else:
            result = "CDS features were annotated."
            status = "correct"
        definition = "Check if CDS features are annotated."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cds_start_end_ids(self, eval_id=None):
        """Check if there are any duplicate start-end coordinates.

        Duplicated start-end coordinates may represent
        unintentional duplicate CDS features.
        """

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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cds_end_strand_ids(self, eval_id=None):
        """Check if there are any duplicate end-strand coordinates.

        Duplicated end-strand coordinates may represent
        unintentional duplicate CDS features with slightly
        different start coordinates.
        """

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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_value_flag(self, expect=False, eval_id=None):
        """Check if there are any attributes that are not populated as expected."""
        if self._value_flag:
            if expect:
                result = "All attributes are populated."
                status = "correct"
            else:
                result = "Some attributes are not populated."
                status = "error"
        else:
            if not expect:
                result = "All attributes are populated."
                status = "correct"
            else:
                result = "Some attributes are not populated."
                status = "error"
        definition = "Check if there are any attributes that are set to %s."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_feature_ids(self, cds=False, trna=False, tmrna=False,
                          other=None, strand=False, eval_id=None):
        """Identify overlapping, duplicated, or partially-duplicated
        features."""
        unsorted_feature_lists = []
        unsorted_features = []
        if cds:
            unsorted_features.extend(self.cds_features)
        if trna:
            unsorted_features.extend(self.trna_features)
        if tmrna:
            unsorted_features.extend(self.tmrna_features)
        if other is not None:
            unsorted_features.extend(other)
        if strand:
            unsorted_f_features = [] # Forward strand
            unsorted_r_features = [] # Reverse strand
            index = 0
            while index < len(unsorted_features):
                feature = unsorted_features[index]
                strand = basic.reformat_strand(feature.strand,
                                               format="fr_short")
                if strand == "f":
                    unsorted_f_features.append(feature)
                else:
                    unsorted_r_features.append(feature)
                index += 1
            unsorted_feature_lists.append(unsorted_f_features)
            unsorted_feature_lists.append(unsorted_r_features)
        else:
            unsorted_feature_lists.append(unsorted_features)
        msgs = ["There are one or more errors with the feature coordinates."]
        for unsorted_features in unsorted_feature_lists:
            sorted_features = sorted(unsorted_features,
                                     key=attrgetter("left", "right"))
            index = 0
            while index < len(sorted_features) - 1:
                current = sorted_features[index]
                next = sorted_features[index + 1]
                if (current.left == next.left and current.right == next.right):
                    msgs.append("Features contain identical left and " \
                                + "right coordinates.")
                elif (current.left < next.left and current.right > next.right):
                    msgs.append("Features are nested.")
                elif (current.left == next.left and \
                        basic.reformat_strand(current.strand,
                            format="fr_short") == "r" and \
                        basic.reformat_strand(next.strand,
                            format="fr_short") == "r"):
                    msgs.append("Features contain the same stop coordinate.")
                elif (current.right == next.right and \
                        basic.reformat_strand(current.strand,
                            format="fr_short") == "f" and \
                        basic.reformat_strand(next.strand,
                            format="fr_short") == "f"):
                    msgs.append("Features contain the same stop coordinate.")
                else:
                    pass
                index += 1
        if len(msgs) > 1:
            result = " ".join(msgs)
            status = "error"
        else:
            result = "The feature coordinates are correct."
            status = "correct"
        definition = "Check if there are any errors with the " \
                     + "genome's feature coordinates."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)





###
