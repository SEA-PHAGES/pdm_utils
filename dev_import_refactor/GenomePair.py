"""Represents a structure to pair two Genome objects and
perform comparisons between them to identify inconsistencies."""

import Eval




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








        # TODO instead of setting attributes, should probably
        # just create EvalResult objects.
        self._host_mismatch = False
        self._accession_mismatch = False
        self._cluster_mismatch = False
        self._subcluster_mismatch = False
        self._author_mismatch = False















    def set_evaluation(self, type, message1 = None, message2 = None):

        if type == "warning":
            eval_object = Eval.construct_warning(message1, message2)

        elif type == "error":
            eval_object = Eval.construct_error(message1)

        else:
            eval_object = Eval.EvalResult()

        self.evaluations.append(eval_object)











    # TODO unit test
    def compare_genome_sequence(self):
        """Compare the sequence of each genome."""

        if self.genome1.sequence != self.genome2.sequence:

            message = "The two genomes have different sequences."
            self.set_evaluation("error", message)

    # TODO unit test
    def compare_genome_length(self):
        """Compare the sequence length of each genome."""

        if self.genome1._length != self.genome2._length:

            message = "The two genomes have different sequence lengths."
            self.set_evaluation("error", message)

    # TODO unit test
    def compare_cluster(self):
        """Compare the cluster of each genome."""

        if self.genome1.cluster != self.genome2.cluster:

            message = "The two genomes are assigned to different clusters."
            self.set_evaluation("error", message)

    # TODO unit test
    def compare_subcluster(self):
        """Compare the subcluster of each genome."""


        if self.genome1.subcluster != self.genome2.subcluster:

            message = "The two genomes are assigned to different subclusters."
            self.set_evaluation("error", message)


    # TODO unit test
    def compare_accession(self):
        """Compare the accession of each genome."""

        if self.genome1.accession != self.genome2.accession:

            message = "The two genomes have different accessions."
            self.set_evaluation("error", message)


    # TODO unit test
    def compare_host(self):
        """Compare the host of each genome."""

        if self.genome1.host != self.genome2.host:

            message = "The two genomes have different hosts."
            self.set_evaluation("error", message)
















#Below: old python2 code



    self.matched_cds
    self.unmatched_cds_genome1
    self.unmatched_cds_genome2
    #self._end_strand_id = ()



    # TODO revamping code.
    # TODO unit test.
    def match_features_by_start_end_strand(self):
        """Match annotated features in each genome with the same start and
        stop coordinates and the same strand (perfect match)."""



        if isinstance(self.genome1,Genome.Genome) and isinstance(self.genome2,Genome.Genome):

            g1_feature_set1 = self.genome1.cds_features_unique_ids
            g1_feature_set1_duplicate_set = self.genome1_cds_features_duplicate_ids

            g2_feature_set1 = self.genome2.cds_features_unique_ids
            g2_feature_set1_duplicate_set = self.genome2_cds_features_duplicate_ids




            # Create the perfect matched and unmatched sets
            feature_set1_matched = g1_feature_set1 & g2_feature_set1
            g1_feature_set1_unmatched = g1_feature_set1 - g2_feature_set1
            g2_feature_set1_unmatched = g2_feature_set1 - g1_feature_set1


            # TODO set matched and unmatched id sets to self.




    # TODO revamping code.
    # TODO unit test.
    def match_features_by_end_strand(self):
        """Match annotated features in each genome with the same stop
        coordinates and same strand but different start coordinates
        (imperfect match)."""

        if isinstance(self.genome1,Genome.Genome) and isinstance(self.genome2,Genome.Genome):



            # From the unmatched sets, created second round of
            # end-strand id sets.
            g1_feature_set = set()
            g2_feature_set = set()

            g1_feature_list = self.unmatched_cds_genome1
            g2_feature_list = self.unmatched_cds_genome2



            # HERE. I may be able to create a method in the Genome class
            # to create a set of unique end_strand identifiers.

            # Phamerator end_strand data

            ph_end_strand_id_set = set()
            ph_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in g1_feature_list:
                if cds.get_start_end_strand_id() in g1_feature_set1_unmatched:
                    if cds.get_end_strand_id() not in ph_end_strand_id_set:
                        ph_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ph_end_strand_duplicate_id_set.add(cds.get_end_strand_id())

            #Remove the duplicate end_strand ids from the main id_set
            ph_end_strand_id_set = ph_end_strand_id_set - ph_end_strand_duplicate_id_set






            ncbi_end_strand_id_set = set()
            ncbi_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in g2_feature_list:
                if cds.get_start_end_strand_id() in g2_feature_set1_unmatched:
                    if cds.get_end_strand_id() not in ncbi_end_strand_id_set:
                        ncbi_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ncbi_end_strand_duplicate_id_set.add(cds.get_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ncbi_end_strand_id_set = ncbi_end_strand_id_set - ncbi_end_strand_duplicate_id_set





            #Create the imperfect matched set
            imperfect_matched_cds_id_set = ph_end_strand_id_set & ncbi_end_strand_id_set


            #Now go back through all cds features and assign to the appropriate dictionary or list
            ph_perfect_matched_cds_dict = {}
            ph_imperfect_matched_cds_dict = {}
            ph_unmatched_cds_list = []


            for cds in g1_feature_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ph_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ph_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ph_unmatched_cds_list.append(cds)


            ncbi_perfect_matched_cds_dict = {}
            ncbi_imperfect_matched_cds_dict = {}
            ncbi_unmatched_cds_list = []

            for cds in g2_feature_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ncbi_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ncbi_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ncbi_unmatched_cds_list.append(cds)


            #Create MatchedCdsFeatures objects
            #Compute matched gene errors and tallies
            #Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1



                if matched_cds_object.get_phamerator_ncbi_different_translations():
                    self.__phamerator_ncbi_different_translations_tally += 1
                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_perfect_matched_features.append(matched_cds_object)


            #Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_imperfect_matched_features.append(matched_cds_object)







            #Compute unmatched error and gene total errors for all unmatched features
            for cds in ph_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

            for cds in ncbi_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1















            #Set unmatched cds lists
            self.__phamerator_features_unmatched_in_ncbi = ph_unmatched_cds_list
            self.__ncbi_features_unmatched_in_phamerator = ncbi_unmatched_cds_list

            #Now compute the number of features in each category
            self.__phamerator_ncbi_perfect_matched_features_tally = len(self.__phamerator_ncbi_perfect_matched_features)
            self.__phamerator_ncbi_imperfect_matched_features_tally = len(self.__phamerator_ncbi_imperfect_matched_features)
            self.__phamerator_features_unmatched_in_ncbi_tally = len(self.__phamerator_features_unmatched_in_ncbi)
            self.__ncbi_features_unmatched_in_phamerator_tally = len(self.__ncbi_features_unmatched_in_phamerator)


        #If there is no matching NCBI genome, assign all Phamerator genes to Unmatched
        else:

            #Set unmatched cds lists, but do NOT count them in the unmatched tally.
            #The unmatched tally should reflect unmatched genes if there is actually a metching NCBI genome.
            self.__phamerator_features_unmatched_in_ncbi = g1_feature_list
































###
