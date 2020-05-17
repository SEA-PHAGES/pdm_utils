"""Represents a structure to group three Genome objects and
perform comparisons between them to identify inconsistencies."""


# TODO this class needs to be refactored, with attributes and methods
# simplified. The class is used in the compare pipeline, which has only
# been partially refactored since integrating into pdm_utils.
# It relies on extra Genome object attributes that are NOT present
# in the base Genome class and that are added during compare pipeline.
# Eventually, GenomeTriad should be generalized and/or replaced with
# implementing multiple GenomePair objects to represent all pairwise
# comparisons.
# So do NOT use this class for anything other than in the compare pipeline
# until it has been properly refactored.

# Variables are prefixed to indicate genome type:
# GenBank =  "gbk", "g"
# MySQL = "mysql", "m"
# PhagesDB = "pdb", "p"

import re

from pdm_utils.functions import basic
from pdm_utils.classes import cdspair

# TODO refactor and test.
class GenomeTriad:
    """Stores three Genome objects."""
    # TODO: similar to GenomePair, but it has three Genome slots instead of two.
    # Many methods here have been improved elsewhere. Ultimately, this class
    # should be removed, and the compare pipeline should rely on GenomePair
    # objects. However, all GenomeTriad methods need to be accounted for,
    # and the compare pipeline would need to run specific checks in specific
    # GenomePair objects, based on the type of paired genomes
    # (MySQL-PhagesDB, MySQL-Genbank, etc.)

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.m_genome = ""
        self.p_genome = ""
        self.g_genome = ""

        # MySQL and GenBank matched data comparison results
        self._m_g_seq_mismatch = False
        self._m_g_seq_length_mismatch = False
        self._g_header_fields_name_mismatch = False
        self._g_host_mismatch = False

        # List of MatchedCdsFeature objects
        self._m_g_perfect_matched_ftrs = []
        self._m_g_imperfect_matched_ftrs = []

        # List of Cds objects
        self._m_ftrs_unmatched_in_g = []
        self._g_ftrs_unmatched_in_m = []

        self._m_g_perfect_matched_ftrs_tally = 0
        self._m_g_imperfect_matched_ftrs_tally = 0
        self._m_ftrs_unmatched_in_g_tally = 0
        self._g_ftrs_unmatched_in_m_tally = 0
        self._m_g_different_descriptions_tally = 0
        self._m_g_different_translations_tally = 0
        self._m_g_author_error = False

        # MySQL and PhagesDB matched data comparison results
        self._m_p_seq_mismatch = False
        self._m_p_seq_length_mismatch = False
        self._m_p_host_mismatch = False
        self._m_p_accession_mismatch = False
        self._m_p_cluster_mismatch = False


        # PhagesDB and GenBank matched data comparison results
        self._p_g_seq_mismatch = False
        self._p_g_seq_length_mismatch = False

        # Total errors summary
        self._contains_errors = False
        self._total_number_genes_with_errors = 0


    # Define all attribute setters:


    # TODO refactor and test.
    def compare_mysql_gbk_genomes(self, gnm_mysql, gnm_gbk, cdspair_mysql_gbk):

        # verify that there is a MySQL and GenBank genome in
        # the matched genome object.
        m_gnm = self.m_genome
        g_gnm = self.g_genome
        m_cds_lst = m_gnm.cds_features

        if m_gnm.type == gnm_mysql and g_gnm.type == gnm_gbk:

            if m_gnm.seq != g_gnm.seq:
                self._m_g_seq_mismatch = True
            if m_gnm.length != g_gnm.length:
                self._m_g_seq_length_mismatch = True

            # Compare phage names
            pattern1 = re.compile("^" + m_gnm.name + "$")
            pattern2 = re.compile("^" + m_gnm.name)

            if (basic.find_expression(pattern2,g_gnm.description.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.source.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.organism.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.source_feature_organism.split(" ")) == 0):

                self._g_header_fields_name_mismatch = True

            # Compare host_genus data
            search_host = m_gnm.host_genus
            if search_host == "Mycobacterium":
                search_host = search_host[:-3]
            pattern3 = re.compile("^" + search_host)

            if ((basic.find_expression(pattern3,g_gnm.description.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.source.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.organism.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.source_feature_organism.split(" ")) == 0) or
                    (g_gnm.source_feature_host != "" and basic.find_expression(pattern3,g_gnm.source_feature_host.split(" ")) == 0) or
                    (g_gnm.source_feature_lab_host != "" and basic.find_expression(pattern3,g_gnm.source_feature_lab_host.split(" ")) == 0)):

                self._g_host_mismatch = True

            # Check author list for errors
            # For genomes with AnnotationAuthor = 1 (Hatfull), Graham is
            # expected to be an author.
            # For genomes with AnnotationAuthor = 0 (non-Hatfull/GenBank),
            # Graham is NOT expected to be an author.
            pattern5 = re.compile("hatfull")
            search_result = pattern5.search(g_gnm.authors.lower())
            if m_gnm.annotation_author == 1 and search_result == None:
                self._m_g_author_error = True
            elif m_gnm.annotation_author == 0 and search_result != None:
                self._m_g_author_error = True
            else:
                # Any other combination of MySQL and GenBank
                # author can be skipped.
                pass


            # Compare CDS features

            # First find all unique start-end-orientation CDS identifiers
            # for MySQL and GenBank genomes.
            m_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            m_start_end_strand_duplicate_id_set = set()

            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id not in m_start_end_strand_id_set:
                    m_start_end_strand_id_set.add(cds_ftr._start_end_strand_id)
                else:
                    m_start_end_strand_duplicate_id_set.add(cds_ftr._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            m_start_end_strand_id_set = \
                m_start_end_strand_id_set - m_start_end_strand_duplicate_id_set


            g_cds_list = g_gnm.cds_features
            g_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            g_start_end_strand_duplicate_id_set = set()
            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id not in g_start_end_strand_id_set:
                    g_start_end_strand_id_set.add(cds_ftr._start_end_strand_id)

                else:
                    g_start_end_strand_duplicate_id_set.add(cds_ftr._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            g_start_end_strand_id_set = \
                g_start_end_strand_id_set - g_start_end_strand_duplicate_id_set

            # Create the perfect matched and unmatched sets
            m_unmatched_start_end_strand_id_set = \
                        m_start_end_strand_id_set - g_start_end_strand_id_set
            g_unmatched_start_end_strand_id_set = \
                        g_start_end_strand_id_set - m_start_end_strand_id_set
            perfect_matched_cds_id_set = \
                        m_start_end_strand_id_set & g_start_end_strand_id_set


            # From the unmatched sets, created second round of
            # end-orientation id sets MySQL end_strand data
            m_end_strand_id_set = set()

            # All end_strand ids that are not unique
            m_end_strand_duplicate_id_set = set()
            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id in m_unmatched_start_end_strand_id_set:
                    if cds_ftr._end_strand_id not in m_end_strand_id_set:
                        m_end_strand_id_set.add(cds_ftr._end_strand_id)
                    else:
                        m_end_strand_duplicate_id_set.add(cds_ftr._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            m_end_strand_id_set = \
                        m_end_strand_id_set - m_end_strand_duplicate_id_set


            g_end_strand_id_set = set()

            # All end_strand ids that are not unique
            g_end_strand_duplicates = set()
            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id in g_unmatched_start_end_strand_id_set:
                    if cds_ftr._end_strand_id not in g_end_strand_id_set:
                        g_end_strand_id_set.add(cds_ftr._end_strand_id)
                    else:
                        g_end_strand_duplicates.add(cds_ftr._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            g_end_strand_id_set = \
                        g_end_strand_id_set - g_end_strand_duplicates


            # Create the imperfect matched set
            imperfect_matched_cds_id_set = \
                        m_end_strand_id_set & g_end_strand_id_set


            # Now go back through all CDS features and assign
            # to the appropriate dictionary or list.
            m_perfect_matched_cds_dict = {}
            m_imperfect_matched_cds_dict = {}
            m_unmatched_cds_list = []


            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id in perfect_matched_cds_id_set:
                    m_perfect_matched_cds_dict[cds_ftr._start_end_strand_id] = cds_ftr
                elif cds_ftr._end_strand_id in imperfect_matched_cds_id_set:
                    m_imperfect_matched_cds_dict[cds_ftr._end_strand_id] = cds_ftr
                else:
                    m_unmatched_cds_list.append(cds_ftr)


            g_perfect_matched_cds_dict = {}
            g_imperfect_matched_cds_dict = {}
            g_unmatched_cds_list = []

            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id in perfect_matched_cds_id_set:
                    g_perfect_matched_cds_dict[cds_ftr._start_end_strand_id] = cds_ftr
                elif cds_ftr._end_strand_id in imperfect_matched_cds_id_set:
                    g_imperfect_matched_cds_dict[cds_ftr._end_strand_id] = cds_ftr
                else:
                    g_unmatched_cds_list.append(cds_ftr)


            # Create CdsPair objects
            # Compute matched gene errors and tallies
            # Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = cdspair.CdsPair()
                matched_cds_object.type = cdspair_mysql_gbk
                matched_cds_object.cds1 = \
                            m_perfect_matched_cds_dict[start_end_strand_tup]
                matched_cds_object.cds2 = \
                            g_perfect_matched_cds_dict[start_end_strand_tup]
                matched_cds_object.compare_cds()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1



                if matched_cds_object.different_translation:
                    self._m_g_different_translations_tally += 1
                if matched_cds_object.different_description:
                    self._m_g_different_descriptions_tally += 1
                self._m_g_perfect_matched_ftrs.append(matched_cds_object)


            # Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = cdspair.CdsPair()
                matched_cds_object.type = cdspair_mysql_gbk
                matched_cds_object.cds1 = \
                                m_imperfect_matched_cds_dict[end_strand_tup]
                matched_cds_object.cds2 = \
                                g_imperfect_matched_cds_dict[end_strand_tup]
                matched_cds_object.compare_cds()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1

                if matched_cds_object.different_description:
                    self._m_g_different_descriptions_tally += 1
                self._m_g_imperfect_matched_ftrs.append(matched_cds_object)


            # Compute unmatched error and gene total errors for
            # all unmatched features.
            for cds_ftr in m_unmatched_cds_list:
                cds_ftr._unmatched_error = True
                cds_ftr.check_for_errors()
                if cds_ftr._total_errors > 0:
                    self._total_number_genes_with_errors += 1

            for cds_ftr in g_unmatched_cds_list:
                cds_ftr._unmatched_error = True
                cds_ftr.check_for_errors()
                if cds_ftr._total_errors > 0:
                    self._total_number_genes_with_errors += 1


            # Set unmatched CDS lists
            self._m_ftrs_unmatched_in_g = m_unmatched_cds_list
            self._g_ftrs_unmatched_in_m = g_unmatched_cds_list

            # Now compute the number of features in each category
            self._m_g_perfect_matched_ftrs_tally = \
                                        len(self._m_g_perfect_matched_ftrs)
            self._m_g_imperfect_matched_ftrs_tally = \
                                        len(self._m_g_imperfect_matched_ftrs)
            self._m_ftrs_unmatched_in_g_tally = len(self._m_ftrs_unmatched_in_g)
            self._g_ftrs_unmatched_in_m_tally = len(self._g_ftrs_unmatched_in_m)


        # If there is no matching GenBank genome,
        # assign all MySQL genes to Unmatched.
        else:

            # Set unmatched CDS lists, but do NOT count them in unmatched tally.
            # Unmatched tally should reflect unmatched genes if there is
            # actually a matching GenBank genome.
            self._m_ftrs_unmatched_in_g = m_cds_lst

            # Now that all errors have been computed for each gene,
            # compute how many genes have errors
            # If there is no matching GenBank genome,
            # gene error tallies are only computed for the MySQL genome.
            m_gnm.compute_genes_with_errors_tally()
            self._total_number_genes_with_errors = m_gnm._genes_with_errors_tally


    # TODO refactor and test.
    def compare_mysql_phagesdb_genomes(self, gnm_mysql, gnm_pdb):

        # verify that there is a MySQL and PhagesDB genome
        # in the matched genome object.
        m_gnm = self.m_genome
        p_gnm = self.p_genome

        if m_gnm.type == gnm_mysql and p_gnm.type == gnm_pdb:

            if m_gnm.seq != p_gnm.seq:
                self._m_p_seq_mismatch = True
            if m_gnm.length != p_gnm.length:
                self._m_p_seq_length_mismatch = True
            if m_gnm.accession != p_gnm.accession:
                self._m_p_accession_mismatch = True
            if m_gnm.host_genus != p_gnm.host_genus:
                self._m_p_host_mismatch = True
            if m_gnm.cluster != p_gnm.cluster:
                self._m_p_cluster_mismatch = True


    # TODO refactor and test.
    def compare_phagesdb_gbk_genomes(self, gnm_pdb, gnm_gbk):

        # verify that there is a PhagesDB and GenBank genome
        # in the matched genome object.
        p_gnm = self.p_genome
        g_gnm = self.g_genome

        if p_gnm.type == gnm_pdb and g_gnm.type == gnm_gbk:

            if p_gnm.seq != g_gnm.seq:
                self._p_g_seq_mismatch = True
            if p_gnm.length != g_gnm.length:
                self._p_g_seq_length_mismatch = True



    # TODO refactor and test.
    def compute_total_genome_errors(self, gnm_gbk, gnm_pdb):
        m_gnm = self.m_genome
        g_gnm = self.g_genome
        p_gnm = self.p_genome

        if self._total_number_genes_with_errors > 0:
            self._contains_errors = True

        # MySQL genome
        if m_gnm._nucleotide_errors:
            self._contains_errors = True
        if m_gnm._status_description_error:
            self._contains_errors = True
        if m_gnm._status_accession_error:
            self._contains_errors = True

        # GenBank genome
        if g_gnm.type == gnm_gbk:
            if g_gnm._nucleotide_errors:
                self._contains_errors = True

        # PhagesDB genome
        if p_gnm.type == gnm_pdb:
            if p_gnm._nucleotide_errors:
                self._contains_errors = True

        # MySQL-GenBank
        if self._m_g_seq_mismatch:
            self._contains_errors = True
        if self._m_g_seq_length_mismatch:
            self._contains_errors = True
        if self._g_header_fields_name_mismatch:
            self._contains_errors = True
        if self._g_host_mismatch:
            self._contains_errors = True
        if self._m_g_imperfect_matched_ftrs_tally > 0:
            self._contains_errors = True
        if self._m_ftrs_unmatched_in_g_tally > 0:
            self._contains_errors = True
        if self._g_ftrs_unmatched_in_m_tally > 0:
            self._contains_errors = True
        if self._m_g_different_descriptions_tally > 0:
            self._contains_errors = True
        if self._m_g_different_translations_tally > 0:
            self._contains_errors = True
        if self._m_g_author_error:
            self._contains_errors = True

        # MySQL-PhagesDB
        if self._m_p_seq_mismatch:
            self._contains_errors = True
        if self._m_p_seq_length_mismatch:
            self._contains_errors = True
        if self._m_p_host_mismatch:
            self._contains_errors = True
        if self._m_p_accession_mismatch:
            self._contains_errors = True
        if self._m_p_cluster_mismatch:
            self._contains_errors = True

        # PhagesDB-GenBank
        if self._p_g_seq_mismatch:
            self._contains_errors = True
        if self._p_g_seq_length_mismatch:
            self._contains_errors = True
