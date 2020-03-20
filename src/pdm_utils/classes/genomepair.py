"""Represents a structure to pair two Genome objects and
perform comparisons between them to identify inconsistencies."""

from pdm_utils.classes import eval
from pdm_utils.functions import basic




class GenomePair:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.genome1 = None
        self.genome2 = None


        # Initialize calculated attributes
        self.evaluations = []

        self.matched_cds_list = [] # List of MatchedCds objects. TODO this contains perfect and imperfect matches
        self.unmatched_cds_genome1 = [] # List of Cds objects
        self.unmatched_cds_genome2 = [] # List of Cds objects



        self.perfect_matched_features_tally = 0
        self.imperfect_matched_features_tally = 0
        self.genome1_features_unmatched_tally = 0
        self.genome2_features_unmatched_tally = 0




    def copy_data(self, attr, first, second, keyword=None):
        """Copy data between paired genomes.

        :param attr:
            Indicates the Genome object attribute used to identify direction
            of copy.
        :type attr: str
        :param first:
            Indicates the value of the indicated attribute to identify
            the Genome object from which data will be copied.
        :type first: str
        :param second:
            Indicates the value of the indicated attribute to identify
            the Genome object to which data will be copied.
        :type second: str
        :param keyword:
            Indicates the Genome attribute from which data should be copied.
            If None, all attributes are copied.
        :type keyword: str
        """



        # Decide the direction of copy.
        direction = True
        if first == second:
            # If unable to identify direction of copy, set to False.
            direction = False
        elif attr == "type":
            if (self.genome1.type == first and self.genome2.type == second):
                first = self.genome1
                second = self.genome2
            elif (self.genome1.type == second and self.genome2.type == first):
                first = self.genome2
                second = self.genome1
            else:
                # If unable to identify direction of copy, set to False.
                direction = False
        elif attr == "id":
            if (self.genome1.id == first and self.genome2.id == second):
                first = self.genome1
                second = self.genome2
            elif (self.genome1.id == second and self.genome2.id == first):
                first = self.genome2
                second = self.genome1
            else:
                # If unable to identify direction of copy, set to False.
                direction = False
        else:
            direction = False


        # Now copy the data.
        if direction:

            # Do not copy over the 'type' attribute that is used to
            # differentiate the two genomes.
            if not attr == "id":
                if (second.id == keyword or keyword is None):
                    second.id = first.id
            if not attr == "type":
                if (second.type == keyword or keyword is None):
                    second.type = first.type

            # All other fields can be copied as needed.
            if (second.name == keyword or keyword is None):
                second.name = first.name
            if (second.host_genus == keyword or keyword is None):
                second.host_genus = first.host_genus
            if (second.seq == keyword or keyword is None):
                second.seq = first.seq
            if (second.accession == keyword or keyword is None):
                second.accession = first.accession
            if (second.annotation_status == keyword or keyword is None):
                second.annotation_status = first.annotation_status
            if (second.cluster == keyword or keyword is None):
                second.cluster = first.cluster
            if (second.subcluster == keyword or keyword is None):
                second.subcluster = first.subcluster
            if (second.date == keyword or keyword is None):
                second.date = first.date
            if (second.annotation_author == keyword or keyword is None):
                second.annotation_author = first.annotation_author
            if (second.retrieve_record == keyword or keyword is None):
                second.retrieve_record = first.retrieve_record
            if (second.description == keyword or keyword is None):
                second.description = first.description
            if (second.source == keyword or keyword is None):
                second.source = first.source
            if (second.organism == keyword or keyword is None):
                second.organism = first.organism
            if (second.authors == keyword or keyword is None):
                second.authors = first.authors
            if (second.filename == keyword or keyword is None):
                second.filename = first.filename
            if (second.translation_table == keyword or keyword is None):
                second.translation_table = first.translation_table
            if (second.length == keyword or keyword is None):
                second.length = first.length
            if (second.gc == keyword or keyword is None):
                second.gc = first.gc
            if (second.evaluations == keyword or keyword is None):
                second.evaluations = first.evaluations
            if (second.cds_features == keyword or keyword is None):
                second.cds_features = first.cds_features
            if (second._cds_features_tally == keyword or keyword is None):
                second._cds_features_tally = first._cds_features_tally
            if (second._cds_start_end_ids == keyword or keyword is None):
                second._cds_start_end_ids = first._cds_start_end_ids
            if (second._cds_end_orient_ids == keyword or keyword is None):
                second._cds_end_orient_ids = first._cds_end_orient_ids
            if (second._cds_descriptions_tally == keyword
                    or keyword is None):
                second._cds_descriptions_tally = \
                    first._cds_descriptions_tally
            if (second.trna_features == keyword or keyword is None):
                second.trna_features = first.trna_features
            if (second._trna_features_tally == keyword or keyword is None):
                second._trna_features_tally = first._trna_features_tally
            if (second.source_features == keyword or keyword is None):
                second.source_features = first.source_features
            if (second._source_features_tally == keyword or keyword is None):
                second._source_features_tally = first._source_features_tally
            if (second._description_name == keyword or keyword is None):
                second._description_name = first._description_name
            if (second._source_name == keyword or keyword is None):
                second._source_name = first._source_name
            if (second._organism_name == keyword or keyword is None):
                second._organism_name = first._organism_name
            if (second._description_host_genus == keyword or keyword is None):
                second._description_host_genus = first._description_host_genus
            if (second._source_host_genus == keyword or keyword is None):
                second._source_host_genus = first._source_host_genus
            if (second._organism_host_genus == keyword or keyword is None):
                second._organism_host_genus = first._organism_host_genus
            if (second._cds_products_tally == keyword
                    or keyword is None):
                second._cds_products_tally = first._cds_products_tally
            if (second._cds_functions_tally == keyword
                    or keyword is None):
                second._cds_functions_tally = first._cds_functions_tally
            if (second._cds_notes_tally == keyword
                    or keyword is None):
                second._cds_notes_tally = first._cds_notes_tally
            if (second._cds_unique_start_end_ids == keyword
                    or keyword is None):
                second._cds_unique_start_end_ids = \
                    first._cds_unique_start_end_ids
            if (second._cds_duplicate_start_end_ids == keyword
                    or keyword is None):
                second._cds_duplicate_start_end_ids = \
                    first._cds_duplicate_start_end_ids
            if (second._cds_unique_end_orient_ids == keyword
                    or keyword is None):
                second._cds_unique_end_orient_ids = \
                    first._cds_unique_end_orient_ids
            if (second._cds_duplicate_end_orient_ids == keyword
                    or keyword is None):
                second._cds_duplicate_end_orient_ids = \
                    first._cds_duplicate_end_orient_ids


    # TODO finish revamping code for matching features.
    # TODO unit test.
    def match_cds_start_stop_ids(self):
        """Match annotated features in each genome with the same start and
        stop coordinates and the same orientation (perfect match)."""


        # Create the perfect matched and unmatched sets.
        start_stop_ids_1 = self.genome1.get_cds_start_stop_ids()
        start_stop_ids_2 = self.genome2.get_cds_start_stop_ids()

        ids_matched, ids_1_unmatched, ids_2_unmatched = \
            FunctionsSimple.match_items(start_stop_ids_1, start_stop_ids_2)

        self.matched_cds_start_stop_ids.add(ids_matched)
        self.unmatched_cds_start_stop_ids_genome1 += ids_1_unmatched
        self.unmatched_cds_start_stop_ids_genome2 += ids_2_unmatched




    # TODO finish revamping code for matching features.
    # TODO unit test.
    def match_cds_stop_ids(self, list1, list2):
        """Match annotated features in each genome with the same stop
        coordinates and same orientation but different start coordinates
        (imperfect match)."""

        # From the unmatched sets, created second round of
        # end-orientation id sets.
        stop_ids_1 = []
        stop_ids_2 = []

        for start_stop_id in list1:
            stop_ids_1.append(start_stop_id[1])

        for start_stop_id in list2:
            stop_ids_2.append(start_stop_id[1])

        # Create the imperfect matched set
        ids_matched, ids_1_unmatched, ids_2_unmatched = \
            FunctionsSimple.match_items(start_stop_ids_1, start_stop_ids_2)


        self.matched_cds_stop_ids.add(ids_matched)
        self.unmatched_cds_stop_ids_genome1 += ids_1_unmatched
        self.unmatched_cds_stop_ids_genome2 += ids_2_unmatched






    # TODO finish revamping code for matching features.
    # TODO unit test.
    def match_cds_features(self):
        """Match annotated features in each genome with the same stop
        coordinates and same orientation but different start coordinates
        (imperfect match)."""

        matched_start_stop_set = set(self.matched_cds_start_stop_ids)
        matched_stop_set = set(self.matched_cds_stop_ids)

        g1_feature_list = self.genome1.cds_features
        g2_feature_list = self.genome2.cds_features


        # Iterate through all cds features of each genome and assign
        # to the appropriate dictionary or list.
        g1_matched_cds_start_stop_dict = {}
        g1_matched_cds_stop_dict = {}
        g1_unmatched_cds_list = []


        for cds_ftr in g1_feature_list:
            if cds_ftr._start_stop_id in matched_start_stop_set:
                g1_matched_cds_start_stop_dict[cds_ftr._start_stop_id] = cds_ftr
            elif cds_ftr._stop_id in matched_stop_set:
                g1_matched_cds_stop_dict[cds_ftr._stop_id] = cds_ftr
            else:
                g1_unmatched_cds_list.append(cds_ftr)


        g2_matched_cds_start_stop_dict = {}
        g2_matched_cds_stop_dict = {}
        g2_unmatched_cds_list = []


        for cds_ftr in g2_feature_list:
            if cds_ftr._start_stop_id in matched_start_stop_set:
                g2_matched_cds_start_stop_dict[cds_ftr._start_stop_id] = cds_ftr
            elif cds_ftr._stop_id in matched_stop_set:
                g2_matched_cds_stop_dict[cds_ftr._stop_id] = cds_ftr
            else:
                g2_unmatched_cds_list.append(cds_ftr)


        # Create MatchedCdsFeatures objects
        for start_stop_tup in matched_start_stop_set:

            matched_cds_object = MatchedCdsFeatures()
            matched_cds_object.feature1 = \
                g1_matched_cds_start_stop_dict[start_stop_tup]
            matched_cds_object.feature2 = \
                g2_matched_cds_start_stop_dict[start_stop_tup]
            # TODO add a step to run all MatchedCDS evaluation methods?

            self.matched_start_stop_features.append(matched_cds_object)


        # Imperfectly matched features
        for end_orient_tup in imperfect_matched_cds_id_set:

            matched_cds_object = MatchedCdsFeatures()
            matched_cds_object.feature1 = \
                g1_matched_cds_stop_dict[end_orient_tup]
            matched_cds_object.feature2 = \
                ncbi_imperfect_matched_cds_dict[end_orient_tup]
            # TODO add a step to run all MatchedCDS evaluation methods?

            self.matched_stop_features.append(matched_cds_object)


        # Compute unmatched error and gene total errors for
        # all unmatched features.
        for cds_ftr in g1_unmatched_cds_list:
            cds_ftr.set_unmatched_error()
            cds_ftr.compute_total_cds_errors()
            if cds_ftr.get_total_errors() > 0:
                self.__total_number_genes_with_errors += 1

        for cds_ftr in ncbi_unmatched_cds_list:
            cds_ftr.set_unmatched_error()
            cds_ftr.compute_total_cds_errors()
            if cds_ftr.get_total_errors() > 0:
                self.__total_number_genes_with_errors += 1


        # Set unmatched cds lists
        self.genome1_unmatched_features = g1_unmatched_cds_list
        self.genome2_unmatched_features = g2_unmatched_cds_list


        # TODO not sure if I need this:
        # # Now compute the number of features in each category
        # self.__phamerator_ncbi_perfect_matched_features_tally = len(self.__phamerator_ncbi_perfect_matched_features)
        # self.__phamerator_ncbi_imperfect_matched_features_tally = len(self.__phamerator_ncbi_imperfect_matched_features)
        # self.__phamerator_features_unmatched_in_ncbi_tally = len(self.__phamerator_features_unmatched_in_ncbi)
        # self.__ncbi_features_unmatched_in_phamerator_tally = len(self.__ncbi_features_unmatched_in_phamerator)


        # TODO make sure this below is accounted for. If there is no
        # matching genome...
        # #If there is no matching NCBI genome, assign all MySQL genes to Unmatched
        # else:
        #
        #     #Set unmatched cds lists, but do NOT count them in the unmatched tally.
        #     #The unmatched tally should reflect unmatched genes if there is actually a metching NCBI genome.
        #     self.__phamerator_features_unmatched_in_ncbi = g1_feature_list

        # TODO above code is still in development.


    # Evaluations
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


    def compare_attribute(self, attribute, expect_same=False, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Compare values of the specified attribute in each genome.

        :param attribute: Name of the GenomePair object attribute to evaluate.
        :type attribute: str
        :param expect_same:
            Indicates whether the two attribute values are expected to be
            the same.
        :type expect_same: bool
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
            value1 = getattr(self.genome1, attribute)
            value2 = getattr(self.genome2, attribute)
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
            result = (f"The first genome is ID: {self.genome1.id}, "
                      f"Type: {self.genome1.type}. The '{attribute}' attribute "
                      f" contains: '{v1_short}'. "
                      f"The second genome is ID: {self.genome2.id}, "
                      f"Type: {self.genome2.type}. The '{attribute}' attribute "
                      f" contains: '{v2_short}'. These two values are ")

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
            result = f"'{attribute}' is not a valid field to be compared."
            status = "untested"
        definition = f"Compare values of the '{attribute}' attribute in each genome."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def compare_date(self, expect, eval_id=None, success="correct",
                     fail="error", eval_def=None):
        """Compare the date of each genome.

        :param expect:
            Is the first genome expected to be "newer", "equal", or "older"
            than the second genome.
        :type expect: str
        :param eval_id: same as for compare_attribute().
        :param success: same as for compare_attribute().
        :param fail: same as for compare_attribute().
        :param eval_def: same as for compare_attribute().
        """

        if expect in set(["newer", "equal", "older"]):
            if self.genome1.date > self.genome2.date:
                actual = "newer"
                actual2 = actual + "than"
            elif self.genome1.date == self.genome2.date:
                actual = "equal"
                actual2 = actual + "to"
            else:
                actual = "older"
                actual2 = actual + "than"

            msg = (f"The query genome '{self.genome1.id}' date "
                   f"is '{self.genome1.date}'."
                   f"The reference genome '{self.genome2.id}' date "
                   f"is '{self.genome2.date}'."
                   f"The date of query genome is {actual2} the "
                   "date of the reference genome, which is ")

            if actual == expect:
                result = msg + "expected."
                status = success
            else:
                result = msg + "not expected."
                status = fail
        else:
            result = f"'{expect}' is an invalid comparison."
            status = "untested"
        definition = "Compare the date of both genomes."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
