"""Represents a collection of data about how databases storing the same data
differ from each other.
"""

# TODO this class needs to be refactored, with attributes and methods
# simplified. The class is used in the compare pipeline, which has only
# been partially refactored since integrating into pdm_utils.
# So do NOT use this class for anything other than in the compare pipeline
# until it has been properly refactored.


# TODO refactor and test.
class DbCompareSummary:

    # Initialize all attributes:
    def __init__(self,matched_genomes_list):

        # Initialize all non-calculated attributes:
        self._matched_genomes_list = matched_genomes_list

        # Initialize all calculated attributes:

        # MySQL data
        # General genome data
        self._m_g_update_flag_tally = 0

        # Genome data checks
        self._m_gnms_with_nucleotide_errors_tally = 0
        self._m_gnms_with_translation_errors_tally = 0
        self._m_gnms_with_boundary_errors_tally = 0
        self._m_gnms_with_status_accession_error_tally = 0
        self._m_gnms_with_status_description_error_tally = 0

        # PhagesDB data
        # Genome data checks
        self._p_gnms_with_nucleotide_errors_tally = 0

        # GenBank data
        # Genome data checks
        self._g_gnms_with_nucleotide_errors_tally = 0
        self._g_gnms_with_translation_errors_tally = 0
        self._g_gnms_with_boundary_errors_tally = 0
        self._g_gnms_with_missing_locus_tags_tally = 0
        self._g_gnms_with_locus_tag_typos_tally = 0
        self._g_gnms_with_description_field_errors_tally = 0

        # MySQL-PhagesDB checks
        self._m_p_seq_mismatch_tally = 0
        self._m_p_seq_length_mismatch_tally = 0
        self._m_p_cluster_mismatch_tally = 0
        self._m_p_accession_mismatch_tally = 0
        self._m_p_host_mismatch_tally = 0

        # MySQL-GenBank checks
        self._m_g_seq_mismatch_tally = 0
        self._m_g_seq_length_mismatch_tally = 0
        self._m_g_header_phage_mismatch_tally = 0
        self._m_g_header_host_mismatch_tally = 0
        self._m_g_gnms_with_imperfectly_matched_ftrs_tally = 0
        self._m_g_gnms_with_unmatched_m_ftrs_tally = 0
        self._m_g_gnms_with_unmatched_g_ftrs_tally = 0
        self._m_g_gnms_with_different_descriptions_tally = 0
        self._m_g_gnms_with_different_translations_tally = 0
        self._m_g_gnms_with_author_errors_tally = 0

        # PhagesDB-GenBank checks
        self._p_g_seq_mismatch_tally = 0
        self._p_g_seq_length_mismatch_tally = 0


        # MySQL feature
        # Gene data checks
        self._m_translation_errors_tally = 0
        self._m_boundary_errors_tally = 0


        # GenBank feature
        # Gene data checks
        self._g_translation_errors_tally = 0
        self._g_boundary_errors_tally = 0
        self._g_missing_locus_tags_tally = 0
        self._g_locus_tag_typos_tally = 0
        self._g_description_field_errors_tally = 0

        # MySQL-GenBank checks
        self._m_g_different_descriptions_tally = 0
        self._m_g_different_start_sites_tally = 0
        self._m_g_different_translation_tally = 0
        self._m_g_unmatched_m_ftrs_tally = 0
        self._m_g_unmatched_g_ftrs_tally = 0


        # Calculate summary metrics
        self._m_total_gnms_analyzed = 0
        self._m_gnms_unmatched_to_p_tally = 0
        self._m_gnms_unmatched_to_g_tally = 0
        self._total_genomes_with_errors = 0


    # TODO refactor and test.
    def compute_summary(self, gnm_mysql, gnm_pdb, gnm_gbk):
        for gnms in self._matched_genomes_list:
            self.compute_matched_genomes_summary(gnms, gnm_mysql, gnm_pdb,
                                                 gnm_gbk)


    # TODO refactor and test.
    def  compute_matched_genomes_summary(self, gnms, gnm_mysql,
                                         gnm_pdb, gnm_gbk):
        """Check errors within matched genomes."""

        self._m_total_gnms_analyzed += 1

        if gnms._contains_errors:
            self._total_genomes_with_errors += 1

        m_gnm = gnms.m_genome
        p_gnm = gnms.p_genome
        g_gnm = gnms.g_genome

        if m_gnm.type == gnm_mysql:
            self._m_g_update_flag_tally += m_gnm.retrieve_record
            self.compute_mysql_gnm_summary(m_gnm)

        if p_gnm.type == gnm_pdb:
            self.compute_pdb_gnm_summary(p_gnm)
        else:
            self._m_gnms_unmatched_to_p_tally += 1

        if g_gnm.type == gnm_gbk:
            self.compute_gbk_gnm_summary(g_gnm)
        else:
            self._m_gnms_unmatched_to_g_tally += 1

        self.compute_mysql_pdb_summary(gnms)
        self.compute_mysql_gbk_summary(gnms)
        self.compute_pdb_gbk_summary(gnms)


    # TODO refactor and test.
    def  compute_mysql_gnm_summary(self, gnm):
        """Check errors within MySQL genome."""

        if gnm._nucleotide_errors:
            self._m_gnms_with_nucleotide_errors_tally += 1
        if gnm._cds_features_with_translation_error_tally > 0:
            self._m_gnms_with_translation_errors_tally += 1
            self._m_translation_errors_tally += \
                                gnm._cds_features_with_translation_error_tally
        if gnm._cds_features_boundary_error_tally > 0:
            self._m_gnms_with_boundary_errors_tally += 1
            self._m_boundary_errors_tally += \
                                gnm._cds_features_boundary_error_tally
        if gnm._status_accession_error:
            self._m_gnms_with_status_accession_error_tally += 1
        if gnm._status_description_error:
            self._m_gnms_with_status_description_error_tally += 1

    # TODO refactor and test.
    def  compute_pdb_gnm_summary(self, gnm):
        """Check errors within PhagesDB genome."""
        if gnm._nucleotide_errors:
            self._p_gnms_with_nucleotide_errors_tally += 1

    # TODO refactor and test.
    def  compute_gbk_gnm_summary(self, gnm):
        """Check errors within GenBank genome."""

        if gnm._nucleotide_errors:
            self._g_gnms_with_nucleotide_errors_tally += 1
        if gnm._cds_features_with_translation_error_tally > 0:
            self._g_gnms_with_translation_errors_tally += 1
            self._g_translation_errors_tally += \
                                gnm._cds_features_with_translation_error_tally
        if gnm._cds_features_boundary_error_tally > 0:
            self._g_gnms_with_boundary_errors_tally += 1
            self._g_boundary_errors_tally += \
                                    gnm._cds_features_boundary_error_tally
        if gnm._missing_locus_tags_tally > 0:
            self._g_gnms_with_missing_locus_tags_tally += 1
            self._g_missing_locus_tags_tally += gnm._missing_locus_tags_tally
        if gnm._locus_tag_typos_tally > 0:
            self._g_gnms_with_locus_tag_typos_tally += 1
            self._g_locus_tag_typos_tally += gnm._locus_tag_typos_tally
        if gnm._description_field_error_tally > 0:
            self._g_gnms_with_description_field_errors_tally += 1
            self._g_description_field_errors_tally += \
                                    gnm._description_field_error_tally


    # TODO refactor and test.
    def  compute_mysql_pdb_summary(self, gnms):
        """Check differences between MySQL and PhagesDB genomes."""

        if gnms._m_p_seq_mismatch:
            self._m_p_seq_mismatch_tally += 1
        if gnms._m_p_seq_length_mismatch:
            self._m_p_seq_length_mismatch_tally += 1
        if gnms._m_p_cluster_mismatch:
            self._m_p_cluster_mismatch_tally += 1
        if gnms._m_p_accession_mismatch:
            self._m_p_accession_mismatch_tally += 1
        if gnms._m_p_host_mismatch:
            self._m_p_host_mismatch_tally += 1

    # TODO refactor and test.
    def  compute_mysql_gbk_summary(self, gnms):
        """Check differences between MySQL and GenBank genomes."""

        if gnms._m_g_seq_mismatch:
            self._m_g_seq_mismatch_tally += 1
        if gnms._m_g_seq_length_mismatch:
            self._m_g_seq_length_mismatch_tally += 1
        if gnms._g_header_fields_name_mismatch:
            self._m_g_header_phage_mismatch_tally += 1
        if gnms._g_host_mismatch:
            self._m_g_header_host_mismatch_tally += 1
        if gnms._m_g_imperfect_matched_ftrs_tally > 0:
            self._m_g_gnms_with_imperfectly_matched_ftrs_tally += 1
            self._m_g_different_start_sites_tally += \
                                    gnms._m_g_imperfect_matched_ftrs_tally
        if gnms._m_ftrs_unmatched_in_g_tally > 0:
            self._m_g_gnms_with_unmatched_m_ftrs_tally += 1
            self._m_g_unmatched_m_ftrs_tally += \
                                    gnms._m_ftrs_unmatched_in_g_tally
        if gnms._g_ftrs_unmatched_in_m_tally > 0:
            self._m_g_gnms_with_unmatched_g_ftrs_tally += 1
            self._m_g_unmatched_g_ftrs_tally += \
                                    gnms._g_ftrs_unmatched_in_m_tally
        if gnms._m_g_different_descriptions_tally > 0:
            self._m_g_gnms_with_different_descriptions_tally += 1
            self._m_g_different_descriptions_tally += \
                                    gnms._m_g_different_descriptions_tally
        if gnms._m_g_different_translations_tally > 0:
            self._m_g_gnms_with_different_translations_tally += 1
            self._m_g_different_translation_tally += \
                                    gnms._m_g_different_translations_tally
        if gnms._m_g_author_error:
            self._m_g_gnms_with_author_errors_tally += 1

    # TODO refactor and test.
    def  compute_pdb_gbk_summary(self, gnms):
        """Check differences between PhagesDB and GenBank genomes."""

        if gnms._p_g_seq_mismatch:
            self._p_g_seq_mismatch_tally += 1
        if gnms._p_g_seq_length_mismatch:
            self._p_g_seq_length_mismatch_tally += 1
