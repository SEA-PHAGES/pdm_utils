"""Represents a structure to pair two Genome objects and
perform comparisons between them to identify inconsistencies."""

from classes import Eval




class GenomePair:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.genome1 = None
        self.genome2 = None


        # Initialize calculated attributes
        self.evaluations = []

        self.matched_cds = [] # List of MatchedCds objects. TODO this contains perfect and imperfect matches
        self.unmatched_cds_genome1 = [] # List of CdsFeature objects
        self.unmatched_cds_genome2 = [] # List of CdsFeature objects



        self.perfect_matched_features_tally = 0
        self.imperfect_matched_features_tally = 0
        self.genome1_features_unmatched_tally = 0
        self.genome2_features_unmatched_tally = 0







    def set_evaluation(self, type, message1 = None, message2 = None):

        if type == "warning":
            eval_object = Eval.construct_warning(message1, message2)

        elif type == "error":
            eval_object = Eval.construct_error(message1)

        else:
            eval_object = Eval.EvalResult()

        self.evaluations.append(eval_object)


    def compare_genome_sequence(self):
        """Compare the sequence of each genome."""

        if self.genome1.sequence != self.genome2.sequence:
            message = "The two genomes have different sequences."
            self.set_evaluation("warning", message, message)


    def compare_genome_length(self):
        """Compare the sequence length of each genome."""

        if self.genome1._length != self.genome2._length:
            message = "The two genomes have different sequence lengths."
            self.set_evaluation("warning", message, message)


    def compare_cluster(self):
        """Compare the cluster of each genome."""

        if self.genome1.cluster != self.genome2.cluster:
            message = "The two genomes are assigned to different clusters."
            self.set_evaluation("warning", message, message)


    def compare_subcluster(self):
        """Compare the subcluster of each genome."""

        if self.genome1.subcluster != self.genome2.subcluster:
            message = "The two genomes are assigned to different subclusters."
            self.set_evaluation("warning", message, message)


    def compare_accession(self):
        """Compare the accession of each genome."""

        if self.genome1.accession != self.genome2.accession:
            message = "The two genomes have different accessions."
            self.set_evaluation("warning", message, message)


    def compare_host(self):
        """Compare the host of each genome."""

        if self.genome1.host != self.genome2.host:
            message = "The two genomes have different hosts."
            self.set_evaluation("warning", message, message)


    # TODO implement this method. Since authorship is not as straightforward
    # as other fields, it is tricky.
    def compare_author(self):
        """Compare the authorship of each genome."""
        pass














    # TODO finish revamping code for matching features.
    # TODO unit test.
    def match_cds_start_stop_ids(self):
        """Match annotated features in each genome with the same start and
        stop coordinates and the same strand (perfect match)."""


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
        coordinates and same strand but different start coordinates
        (imperfect match)."""

        # From the unmatched sets, created second round of
        # end-strand id sets.
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
        coordinates and same strand but different start coordinates
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


        for cds in g1_feature_list:
            if cds._start_stop_id in matched_start_stop_set:
                g1_matched_cds_start_stop_dict[cds._start_stop_id] = cds
            elif cds._stop_id in matched_stop_set:
                g1_matched_cds_stop_dict[cds._stop_id] = cds
            else:
                g1_unmatched_cds_list.append(cds)


        g2_matched_cds_start_stop_dict = {}
        g2_matched_cds_stop_dict = {}
        g2_unmatched_cds_list = []


        for cds in g2_feature_list:
            if cds._start_stop_id in matched_start_stop_set:
                g2_matched_cds_start_stop_dict[cds._start_stop_id] = cds
            elif cds._stop_id in matched_stop_set:
                g2_matched_cds_stop_dict[cds._stop_id] = cds
            else:
                g2_unmatched_cds_list.append(cds)


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
        for end_strand_tup in imperfect_matched_cds_id_set:

            matched_cds_object = MatchedCdsFeatures()
            matched_cds_object.feature1 = \
                g1_matched_cds_stop_dict[end_strand_tup]
            matched_cds_object.feature2 = \
                ncbi_imperfect_matched_cds_dict[end_strand_tup]
            # TODO add a step to run all MatchedCDS evaluation methods?

            self.matched_stop_features.append(matched_cds_object)


        # Compute unmatched error and gene total errors for
        # all unmatched features.
        for cds in g1_unmatched_cds_list:
            cds.set_unmatched_error()
            cds.compute_total_cds_errors()
            if cds.get_total_errors() > 0:
                self.__total_number_genes_with_errors += 1

        for cds in ncbi_unmatched_cds_list:
            cds.set_unmatched_error()
            cds.compute_total_cds_errors()
            if cds.get_total_errors() > 0:
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
        # #If there is no matching NCBI genome, assign all Phamerator genes to Unmatched
        # else:
        #
        #     #Set unmatched cds lists, but do NOT count them in the unmatched tally.
        #     #The unmatched tally should reflect unmatched genes if there is actually a metching NCBI genome.
        #     self.__phamerator_features_unmatched_in_ncbi = g1_feature_list




###
