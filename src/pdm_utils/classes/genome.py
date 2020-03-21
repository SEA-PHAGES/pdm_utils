"""Represents a collection of data about a genome that are commonly used to
maintain and update SEA-PHAGES phage genomics data.
"""
from operator import attrgetter

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqUtils import GC

from pdm_utils.classes import eval, cds, trna, source
from pdm_utils.constants import constants
from pdm_utils.functions import basic

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
        self.seq = Seq("", IUPAC.ambiguous_dna) # Biopython Seq object
        self.length = 0 # Size of the nucleotide sequence
        self.gc = -1 # %GC content
        self.host_genus = ""
        self.accession = ""
        self.lifestyle = "" # E.g. temperate, obligately lytic, unknown, etc.
        self.translation_table = 0

        # The following attributes are common to MySQL database.
        self.annotation_status = "" # Final, Draft, Unknown version of genome data
        self.annotation_author = -1 # 1 (can be changed), 0 (can not be changed)
        self.retrieve_record = -1 # 1 (auto update), 0 (do not auto update)
        self.date = constants.EMPTY_DATE # Used for the DateLastModified field.

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
        self.authors = "" # Compiled string of all authors named in the record.

        # The following attributes are computed datafields that are
        # common to annotated genomes.
        self.cds_features = [] # List of all parsed CDS features
        self._cds_features_tally = 0
        self._cds_start_end_ids = []
        self._cds_end_orient_ids = []
        self._cds_descriptions_tally = 0
        self._cds_products_tally = 0
        self._cds_functions_tally = 0
        self._cds_notes_tally = 0
        self._cds_unique_start_end_ids = set() # TODO still in development.
        self._cds_duplicate_start_end_ids = set() # TODO still in development.
        self._cds_unique_end_orient_ids = set() # TODO still in development.
        self._cds_duplicate_end_orient_ids = set() # TODO still in development.
        self.trna_features = []
        self._trna_features_tally = 0
        self.tmrna_features = []
        self._tmrna_features_tally = 0
        self.source_features = []
        self._source_features_tally = 0

        # The following attributes are useful for processing data
        # from various data sources.
        self.filename = "" # The file name from which the data is derived
        self.type = "" # Identifier to describes how this genome is used
                       # (e.g. import, MySQL database, PhagesDB, etc.)
        self.evaluations = [] # List of warnings and errors about the data
        self.misc = None # Unstructured attribute to store misc. data not
                         # applicable to any of the other attributes and that
                         # may be used differently for downstream applications.

    def __str__(self):
        str_list = []
        if self.id != "":
            str_list.append(f"ID: {self.id}")
        if self.type != "":
            str_list.append(f"Type: {self.type}")
        if self.filename != "":
            str_list.append(f"Filename: {self.filename}")
        return ", ".join(str_list)


    def set_filename(self, filepath):
        """Set the filename. Discard the path and file extension.

        :param filepath: name of the file reference.
        :type filepath: Path
        """
        self.filename = filepath.stem
        # split_filepath = value.split('/')
        # filename = split_filepath[-1]
        # filename = filename.split('.')[0]
        # self.filename = filename

    def set_id(self, value=None, attribute=None):
        """Set the id from either an input value or an indicated attribute.

        :param value:
            unique identifier for the genome.
        :type value: str
        :param attribute:
            name of a genome object attribute that stores
            a unique identifier for the genome.
        :type attribute: str
        """

        if value is None:
            if (attribute is not None and hasattr(self, attribute)):
                value = getattr(self, attribute)
            else:
                value = ""
        self.id = basic.edit_suffix(value, "remove")


    def set_host_genus(self, value=None, attribute=None, format="empty_string"):
        """Set the host_genus from a value parsed from the indicated attribute.

        The input data is split into multiple parts, and the first word is
        used to set host_genus.

        :param value:
            the host genus of the phage genome
        :type value: str
        :param attribute:
            the name of the genome attribute from which the
            host_genus attribute will be set
        :type attribute: str
        :param format:
            the default format if the input is an empty/null value.
        :type format: str
        """

        if value is None:
            if (attribute is not None and hasattr(self, attribute)):
                value = getattr(self, attribute)
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

        :param value: the genome's nucleotide sequence.
        :type value: str or Seq
        """
        if not isinstance(value, Seq):
            self.seq = Seq(value).upper()
        else:
            self.seq = value.upper()
        self.length = len(self.seq)
        if self.length > 0:
            self.gc = round(GC(self.seq), 4)
        else:
            self.gc = -1


    def set_accession(self, value, format="empty_string"):
        """Set the accession.

        The Accession field in the MySQL database defaults to ''.
        Some flat file accessions have the version number suffix, so discard
        the version number.

        :param value:
            GenBank accession number.
        :type value: str
        :param format:
            indicates the format of the data if it is not a
            valid accession. Default is ''.
        :type format: misc.
        """
        if isinstance(value, str):
            value = value.strip()
            value = value.split(".")[0]
        self.accession = basic.convert_empty(value, format)


    def set_cds_features(self, value):
        """Set and tally the CDS features.

        :param value: list of Cds objects.
        :type value: list
        """
        self.cds_features = value # Should be a list.
        self._cds_features_tally = len(self.cds_features)


    def set_cds_id_list(self):
        """Creates lists of CDS feature identifiers.

        The first identifier is derived from the start and end coordinates.
        The second identifier is derived from the transcription end
        coordinate and orientation.
        """
        start_end_id_list = []
        end_orient_id_list = []
        for cds_ftr in self.cds_features:
            start_end_id_list.append(cds_ftr._start_end_id)
            end_orient_id_list.append(cds_ftr._end_orient_id)
        self._cds_start_end_ids = start_end_id_list
        self._cds_end_orient_ids = end_orient_id_list


    def set_trna_features(self, value):
        """Set and tally the tRNA features.

        :param value: list of Trna objects.
        :type value: list
        """
        self.trna_features = value # Should be a list
        self._trna_features_tally = len(self.trna_features)


    def set_source_features(self, value):
        """Set and tally the source features.

        :param value: list of Source objects.
        :type value: list
        """
        self.source_features = value # Should be a list
        self._source_features_tally = len(self.source_features)


    def set_cluster(self, value):
        """Set the cluster and modify singleton if needed.

        :param value: Cluster designation of the genome.
        :type value: str
        """
        singleton = "Singleton"
        if isinstance(value, str):
            value = value.strip()

            # PhagesDB-output format.
            if value.capitalize() == singleton:
                self.cluster = singleton
            else:
                self.cluster = value

        # MySQL database-output format
        if value is None:
            self.cluster = singleton


    def set_subcluster(self, value):
        """Set the subcluster.

        :param value: Subcluster designation of the genome.
        :type value: str
        """
        if isinstance(value, str):
            self.subcluster = value.strip()

        # MySQL database-output format
        if value is None:
            self.subcluster = "none"


    def set_date(self, value, format="empty_datetime_obj"):
        """Set the date attribute.

        :param value: Date
        :type value: misc
        :param format: Indicates the format if the value is empty.
        :type format: str
        """
        self.date = basic.convert_empty(value, format)


    def set_annotation_author(self, value):
        """Convert annotation_author to integer value if possible.

        :param value: Numeric value.
        :type value: str, int
        """
        try:
            self.annotation_author = int(value)
        except:
            self.annotation_author = value


    def set_retrieve_record(self, value):
        """Convert retrieve_record to integer value if possible.

        :param value: Numeric value.
        :type value: str, int
        """
        try:
            self.retrieve_record = int(value)
        except:
            self.retrieve_record = value


    def set_cds_descriptions(self, value):
        """Set each CDS processed description as indicated.

        :param value: Name of the description field.
        :type value: str
        """
        for x in range(len(self.cds_features)):
            cds_ftr = self.cds_features[x]
            cds_ftr.set_description(value)


    def tally_cds_descriptions(self):
        """Tally the non-generic CDS descriptions."""
        self._cds_descriptions_tally = 0
        self._cds_products_tally = 0
        self._cds_functions_tally = 0
        self._cds_notes_tally = 0
        for cds_ftr in self.cds_features:
            if cds_ftr.description != "":
                self._cds_descriptions_tally += 1
            if cds_ftr.product != "":
                self._cds_products_tally += 1
            if cds_ftr.function != "":
                self._cds_functions_tally += 1
            if cds_ftr.note != "":
                self._cds_notes_tally += 1


    def set_unique_cds_start_end_ids(self):
        """Identify CDS features contain unique start-end coordinates."""
        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_start_end_ids)
        self._cds_unique_start_end_ids = set(unique_id_tuples)
        self._cds_duplicate_start_end_ids = set(duplicate_id_tuples)


    def set_unique_cds_end_orient_ids(self):
        """Identify CDS features contain unique transcription end-orientation coordinates."""
        unique_id_tuples, duplicate_id_tuples = \
            basic.identify_unique_items(self._cds_end_orient_ids)
        self._cds_unique_end_orient_ids = set(unique_id_tuples)
        self._cds_duplicate_end_orient_ids = set(duplicate_id_tuples)


    def set_feature_ids(self, use_type=False, use_cds=False,
                        use_trna=False, use_tmrna=False, use_source=False):
        """Sets the id of each feature.

        Lists of features can be added to this method. The method assumes
        that all elements in all lists contain 'id', 'start', and 'stop'
        attributes. This feature attribute is processed within
        the Genome object because and not within the feature itself since
        the method sorts all features and generates systematic IDs based on
        feature order in the genome.

        :param use_type:
            Indicates whether the type of object should be
            added to the feature id.
        :type use_type: bool
        :param use_cds:
            Indicates whether ids for CDS features should be generated.
        :type use_cds: bool
        :param use_trna:
            Indicates whether ids for tRNA features should be generated.
        :type use_trna: bool
        :param use_tmrna:
            Indicates whether ids for tmRNA features should be generated.
        :type use_tmrna: bool
        :param use_source:
            Indicates whether ids for source features should be generated.
        :type use_source: bool
        """

        # Both coordinates are used to control the order of features
        # that may share one, but not both, coordinates (e.g. tail
        # assembly chaperone).
        list_to_sort = []
        if use_cds:
            list_to_sort.extend(self.cds_features)
        if use_source:
            list_to_sort.extend(self.source_features)
        if use_trna:
            list_to_sort.extend(self.trna_features)

        # TODO unit test after tmrna_features are implemented.
        if use_tmrna:
            list_to_sort.extend(self.tmrna_features)

        sorted_list = sorted(list_to_sort, key=attrgetter("start", "stop"))

        for index in range(len(sorted_list)):
            if use_type:
                if isinstance(sorted_list[index], cds.Cds):
                    delimiter = "_CDS_"
                elif isinstance(sorted_list[index], source.Source):
                    delimiter = "_SRC_"
                elif isinstance(sorted_list[index], trna.TrnaFeature):
                    delimiter = "_TRNA_"

                # TODO unit test after tmRNA class implemented.
                # elif isinstance(sorted_list[index], tmrna.Tmrna):
                #     delimiter = "_TMRNA_"
                else:
                    delimiter = "_"
            else:
                delimiter = "_"

            sorted_list[index].id = self.id + delimiter + str(index + 1)

    # TODO add parameters to specify which feature types (e.g. cds=True, trna=True, ...)
    def clear_locus_tags(self):
        """Resets locus_tags to empty string."""
        x = 0
        while x < len(self.cds_features):
            self.cds_features[x].locus_tag = ""
            x += 1

        # TODO implement for tRNA and tmRNA feature lists.
        # y = 0
        # while y < len(self.trna_features):
        #     self.trna_features[y].locus_tag = ""
        #     y += 1
        #
        # z = 0
        # while z < len(self.tmrna_features):
        #     self.tmrna_features[z].locus_tag = ""
        #     z += 1



    # TODO unittest
    def set_feature_genome_ids(self, feature_type, value=None):
        """Sets the genome_id of each feature.

        :param feature_type:
            Type of features to set genome_id for (CDS, tRNA, etc.)
        :type feature_type: str
        :param value: Genome identifier.
        :type value: str
        """
        if value is None:
            value = self.id

        if feature_type.lower() == "cds":
            feature_list = self.cds_features
        elif feature_type.lower() == "source":
            feature_list = self.source_features
        # TODO implement.
        # elif feature_type.lower() == "trna":
        #     feature_list = self.trna_features
        # elif feature_type.lower() == "tmrna":
        #     feature_list = self.source_features
        else:
            feature_list = []

        for feature in feature_list:
            feature.genome_id = value


    # Evaluations.
    def set_eval(self, eval_id, definition, result, status):
        """Constructs and adds an Eval object to the evaluations list.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param definition: Description of the evaluation.
        :type definition: str
        :param result: Description of the outcome of the evaluation.
        :type result: str
        :param status: Outcome of the evaluation.
        :type status: str
        """
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_attribute(self, attribute, check_set, expect=False, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check that the attribute value is valid.

        :param attribute: Name of the Genome object attribute to evaluate.
        :type attribute: str
        :param check_set: Set of reference ids.
        :type check_set: set
        :param expect:
            Indicates whether the attribute value is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param success: Default status if the outcome is a success.
        :type success: str
        :param fail: Default status if the outcome is not a success.
        :type fail: str
        :param eval_def: Description of the evaluation.
        :type eval_def: str
        """
        try:
            test = True
            value1 = getattr(self, attribute)
        except:
            test = False
            value1 = None
        if test:
            value1_short = basic.truncate_value(str(value1), 30, "...")
            result = f"The {attribute} value '{value1_short}' is "

            value2 = basic.check_value_expected_in_set(
                        value1, check_set, expect)
            if value2:
                result = result + "valid."
                status = success
            else:
                result = result + "not valid."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"
        definition = f"Check the value of the '{attribute}' attribute."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def compare_two_attributes(self, attribute1, attribute2,
                               expect_same=False, eval_id=None,
                               success="correct", fail="error", eval_def=None):
        """Determine if two attributes are the same.

        :param attribute1: First attribute to compare.
        :type attribute1: str
        :param attribute2: Second attribute to compare.
        :type attribute2: str
        :param expect_same:
            Indicates whether the two attribute values are expected to be
            the same.
        :type expect_same: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        try:
            test = True
            value1 = getattr(self, attribute1)
            value2 = getattr(self, attribute2)
        except:
            test = False
            value1 = None
            value2 = None
        if test:
            if value1 == value2:
                actual_same = True
            else:
                actual_same = False

            v1_short = basic.truncate_value(str(value1), 30, "...")
            v2_short = basic.truncate_value(str(value2), 30, "...")
            result = (f"The '{attribute1}' attribute contains: '{v1_short}'. "
                      f"The '{attribute2}' attribute contains: '{v2_short}'. "
                      "These two values are ")

            if actual_same:
                result = result + "identical, "
            else:
                result = result + "different, "

            if actual_same and expect_same:
                result = result + "as expected."
                status = success
            elif not actual_same and not expect_same:
                result = result + "as expected."
                status = success
            else:
                result = result + "which is not expected."
                status = fail
        else:
            result = (f"'{attribute1}' and/or '{attribute2}' is "
                      "not a valid field to be compared.")
            status = "untested"
        definition = f"Compare values of '{attribute1}' and '{attribute2}' attributes."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_cluster_structure(self, eval_id=None, success="correct",
                                fail="error", eval_def=None):
        """Check whether the cluster attribute is structured appropriately.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        result = f"The Cluster is '{self.cluster}'. "
        if self.cluster != "none":
            result = result + "It is structured "
            left, right = basic.split_string(self.cluster)

            if (right != "" or left.isalpha() == False):
                result = result + "incorrectly."
                status = fail
            else:
                result = result + "correctly."
                status = success
        else:
            result = result + "It is empty."
            status = "untested"
        definition = "Check if the Cluster attribute is structured correctly."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_subcluster_structure(self, eval_id=None, success="correct",
                                   fail="error", eval_def=None):
        """Check whether the subcluster attribute is structured appropriately.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        result = f"The Subcluster is '{self.subcluster}'. "
        if (self.subcluster != "none" and self.subcluster != ""):
            result = result + "It is structured "
            left, right = basic.split_string(self.subcluster)
            if (left.isalpha() == False or right.isdigit() == False):
                result = result + "incorrectly."
                status = fail
            else:
                result = result + "correctly."
                status = success
        else:
            result = result + "It is empty."
            status = "untested"
        definition = "Check if the Subcluster attribute is structured correctly."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_compatible_cluster_and_subcluster(self, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check compatibility of cluster and subcluster attributes.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        result = (f"The Cluster is '{self.cluster}', "
                  f"the Subcluster is '{self.subcluster}', and they are ")
        output = basic.compare_cluster_subcluster(self.cluster, self.subcluster)
        if not output:
            result = result + "not compatible."
            status = fail
        else:
            result = result + "compatible."
            status = success
        definition = "Check for compatibility between Cluster and Subcluster."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_nucleotides(self, check_set=set(), eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Check if all nucleotides in the sequence are expected.

        :param check_set: Set of reference nucleotides.
        :type check_set: set
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        # When Biopython SeqIO parses the GenBank record, it automatically
        # determines that it is a DNA sequence. It assigns the Seq object
        # alphabet as IUPACAmbiguousDNA. The alphabet could be coerced
        # to a different alphabet, and then tested using the
        # Bio.Alphabet._verify_alphabet() function. Since this is a private
        # function though, it is not clear how stable/reliable it is.
        # Instead, Bio.Alphabet.IUPAC.unambiguous_dna alphabet can be passed
        # to the check_nucleotides method.
        nucleotide_set = set(self.seq)
        nucleotide_error_set = nucleotide_set - check_set
        if len(nucleotide_error_set) > 0:
            nes_string = basic.join_strings(nucleotide_error_set, delimiter=", ")
            result = ("There are unexpected nucleotides in the sequence: "
                      f"{nes_string}.")
            status = fail
        else:
            result = "There are no unexpected nucleotides in the sequence."
            status = success
        definition = "Check if all nucleotides in the sequence are expected."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_authors(self, check_set=set(), expect=True, eval_id=None,
                      success="correct", fail="error", eval_def=None):
        """Check author list.

        Evaluates whether at least one author in the in the list of
        authors is present in a set of reference authors.

        :param check_set: Set of reference authors.
        :type check_set: set
        :param expect:
            Indicates whether at least one author in the
            list of authors is expected to be present in the check set.
        :type expect: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        authors = self.authors.lower()
        authors = authors.replace(";", ",")
        authors = authors.replace(" ", ",")
        authors_list = authors.split(",")
        authors_list = [x.strip() for x in authors_list]
        authors_set = set(authors_list)
        mutual_authors_set = authors_set & check_set
        if len(mutual_authors_set) == 0:
            missing_set = check_set - authors_set
            missing_string = basic.join_strings(missing_set, delimiter=", ")
            result = ("The following authors are not "
                      f"listed: {missing_string}. This is ")
            if expect:
                result = result + "not expected."
                status = fail
            else:
                result = result + "expected."
                status = success
        else:
            mas_string = basic.join_strings(mutual_authors_set, delimiter=", ")
            result = f"The following authors are listed: {mas_string}. This is "
            if expect:
                result = result + "expected."
                status = success
            else:
                result = result + "not expected."
                status = fail
        definition = "Check authorship."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_magnitude(self, attribute, expect, ref_value, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check that the magnitude of a numerical attribute is valid.

        :param attribute: same as for check_attribute().
        :param expect:
            Comparison symbol indicating direction of magnitude (>, =, <).
        :type expect: str
        :param ref_value: Numerical value for comparison.
        :type ref_value: int, float, datetime
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        try:
            test = True
            query_value = getattr(self, attribute)
        except:
            test = False
            query_value = None
        if test:
            result = f"The {attribute} value {query_value} is "
            if query_value > ref_value:
                compare = ">"
                result = result + "greater than "
            elif query_value == ref_value:
                compare = "="
                result = result + "equal to "
            else:
                compare = "<"
                result = result + "less than "
            result = result + f"{ref_value}, which is "
            if compare == expect:
                result = result + "expected."
                status = success
            else:
                result = result + "not expected."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"
        definition = f"Check the magnitude of the '{attribute}' attribute."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO this may no longer be needed.
    def check_cds_start_end_ids(self, eval_id=None, success="correct",
                                fail="error", eval_def=None):
        """Check if there are any duplicate start-end coordinates.

        Duplicated start-end coordinates may represent
        unintentional duplicate CDS features.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """

        if len(self._cds_duplicate_start_end_ids) > 0:
            result = ("There are multiple CDS features with the same "
                      "start and end coordinates.")
            status = fail
        else:
            result = ("All CDS features contain unique start and "
                      "end coordinate information.")
            status = success
        definition = ("Check whether CDS features can be uniquely "
                      "identified by their start and end coordinates.")
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    # TODO this may no longer be needed.
    def check_cds_end_orient_ids(self, eval_id=None, success="correct",
                                 fail="error", eval_def=None):
        """Check if there are any duplicate transcription end-orientation coordinates.

        Duplicated transcription end-orientation coordinates may represent
        unintentional duplicate CDS features with slightly
        different start coordinates.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """

        if len(self._cds_duplicate_end_orient_ids) > 0:
            result = ("There are multiple CDS features with the same "
                      "transcription end coordinate and orientation.")
            status = fail
        else:
            result = ("All CDS features contain unique orientation and "
                      "transcription end coordinate information.")
            status = success
        definition = ("Check whether CDS features can be uniquely "
                      "identified by their orientation and transcription end "
                      "coordinate.")
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_feature_coordinates(self, cds_ftr=False, trna_ftr=False,
                tmrna_ftr=False, other=None, strand=False, eval_id=None,
                success="correct", fail="error", eval_def=None):
        """Identify nested, duplicated, or partially-duplicated
        features.

        :param cds_ftr: Indicates whether ids of CDS features should be included.
        :type cds_ftr: bool
        :param trna_ftr: Indicates whether ids of tRNA features should be included.
        :type trna_ftr: bool
        :param tmrna_ftr: Indicates whether ids of tmRNA features should be included.
        :type tmrna_ftr: bool
        :param other: List of features that should be included.
        :type other: list
        :param strand: Indicates if feature orientation should be included.
        :type strand: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        unsorted_feature_lists = []
        unsorted_features = []
        ftr_types = set()
        if cds_ftr:
            ftr_types.add("cds")
            unsorted_features.extend(self.cds_features)
        if trna_ftr:
            ftr_types.add("trna")
            unsorted_features.extend(self.trna_features)
        if tmrna_ftr:
            ftr_types.add("tmrna")
            unsorted_features.extend(self.tmrna_features)
        if other is not None:
            ftr_types.add("other")
            unsorted_features.extend(other)
        if strand:
            s_info = "were"
            unsorted_f_features = [] # Forward orientation
            unsorted_r_features = [] # Reverse orientation
            for index in range(len(unsorted_features)):
                feature = unsorted_features[index]
                strand = basic.reformat_strand(feature.orientation,
                                               format="fr_short")
                if strand == "f":
                    unsorted_f_features.append(feature)
                else:
                    unsorted_r_features.append(feature)
            unsorted_feature_lists.append(unsorted_f_features)
            unsorted_feature_lists.append(unsorted_r_features)
        else:
            s_info = "were not"
            unsorted_feature_lists.append(unsorted_features)
        ft_string = basic.join_strings(ftr_types, delimiter=", ")
        result = (f"The following types of features were evaluated: {ft_string}. "
                  f"Features {ft_string} separately grouped "
                  "by orientation for evaluation. ")

        # result = (f"The following types of features were evaluated {s_info} "
        #           f"regard to feature orientation: {ft_string}. ")

        msgs = ["There are one or more errors with the feature coordinates."]
        for unsorted_features in unsorted_feature_lists:
            sorted_features = sorted(unsorted_features,
                                     key=attrgetter("start", "stop"))
            index = 0
            while index < len(sorted_features) - 1:
                current = sorted_features[index]
                next = sorted_features[index + 1]
                ftrs = (f"Feature1 ID: {current.id}, "
                        f"start coordinate: {current.start}, "
                        f"stop coordinate: {current.stop}, "
                        f"orientation: {current.orientation}. "
                        f"Feature2 ID: {next.id}, "
                        f"start coordinate: {next.start}, "
                        f"stop coordinate: {next.stop}, "
                        f"orientation: {next.orientation}. ")

                if (current.start == next.start and current.stop == next.stop):
                    msgs.append(ftrs)
                    msgs.append("Feature1 and Feature2 contain identical "
                                "start and stop coordinates.")

                # To identify nested features, the following tests
                # avoid false errors due to genes that may wrap around the
                # genome.
                elif (current.start < next.start and
                      current.start < next.stop and
                      current.stop > next.start and
                      current.stop > next.stop):
                    msgs.append(ftrs)
                    msgs.append("Feature2 is nested within Feature1.")
                elif (current.start == next.start and
                        basic.reformat_strand(current.orientation,
                            format="fr_short") == "r" and
                        basic.reformat_strand(next.orientation,
                            format="fr_short") == "r"):
                    msgs.append(ftrs)
                    msgs.append(("Feature1 and Feature2 contain "
                                 "identical stop coordinates."))
                elif (current.stop == next.stop and
                        basic.reformat_strand(current.orientation,
                            format="fr_short") == "f" and
                        basic.reformat_strand(next.orientation,
                            format="fr_short") == "f"):
                    msgs.append(ftrs)
                    msgs.append(("Feature1 and Feature2 contain "
                                 "identical stop coordinates."))
                else:
                    pass
                index += 1
        if len(msgs) > 1:
            result = result + " ".join(msgs)
            status = fail
        else:
            result = result + "The feature coordinates are correct."
            status = success
        definition = ("Check if there are any feature coordinate conflicts.")
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
