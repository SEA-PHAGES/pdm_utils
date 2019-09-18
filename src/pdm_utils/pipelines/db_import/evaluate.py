"""Groups of evaluation/check functions."""


from pdm_utils.constants import constants
from pdm_utils.functions import basic

def check_bundle_for_import(bndl):
    """Check a Bundle for errors.

    Evaluate whether all genomes have been successfully grouped,
    and whether all genomes have been paired, as expected.
    Based on the ticket type, there are expected to be certain
    types of genomes and pairs of genomes in the bundle.

    :param bndl: A pdm_utils Bundle object.
    :type bndl: Bundle
    """
    bndl.check_ticket(eval_id="BNDL_001")
    if bndl.ticket is not None:
        tkt = bndl.ticket
        bndl.check_genome_dict("add", eval_id="BNDL_002")
        bndl.check_genome_dict("flat_file", eval_id="BNDL_003")
        bndl.check_genome_pair_dict("flat_file_add", eval_id="BNDL_004")

        # There may or may not be data retrieved from PhagesDB.
        tkt.set_value_flag("retrieve")
        if tkt._value_flag:
            bndl.check_genome_dict("phagesdb", eval_id="BNDL_005")
            bndl.check_genome_pair_dict("flat_file_phagesdb",
                                        eval_id="BNDL_006")

        if tkt.type == "replace":
            bndl.check_genome_dict("phamerator", eval_id="BNDL_007")

            # There may or may not be a genome_pair to retain some data.
            tkt.set_value_flag("retain")
            if tkt._value_flag:
                bndl.check_genome_pair_dict("add_phamerator",
                                            eval_id="BNDL_008")

            # There should be a genome_pair between the current phamerator
            # genome and the new flat_file genome.
            bndl.check_genome_pair_dict("flat_file_phamerator",
                                        eval_id="BNDL_009")


def check_ticket_structure(tkt, type_set=set(), description_field_set=set(),
        null_set=set(), run_mode_set=set(), id_dupe_set=set(),
        phage_id_dupe_set=set(), accession_dupe_set=set()):
    """Evaluate a ticket to confirm it is structured appropriately.
    The assumptions for how each field is populated varies depending on
    the type of ticket.

    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param description_field_set: Valid description_field options.
    :type description_field_set: set
    :param null_set: Values that represent an empty field.
    :type null_set: set
    :param run_mode_set: Valid run mode options.
    :type run_mode_set: set
    :param id_dupe_set: Predetermined duplicate ticket ids.
    :type id_dupe_set: set
    :param phage_id_dupe_set: Predetermined duplicate PhageIDs.
    :type phage_id_dupe_set: set
    :param accession_dupe_set: Predetermined duplicate accessions.
    :type accession_dupe_set: set
    """
    # This function simply evaluates whether there is data in the
    # appropriate ticket attributes given the type of ticket.
    # It confirms that ticket attributes 'type', 'run_mode', and
    # 'description_field' are populated with specific values.
    # But it does not evaluate the quality of the data itself for
    # the other fields, since those are genome-specific fields and
    # can be checked within Genome objects.

    # Check for duplicated values.
    tkt.check_duplicate_id(id_dupe_set, eval_id="TKT_001")
    tkt.check_duplicate_phage_id(phage_id_dupe_set, eval_id="TKT_002")
    tkt.check_duplicate_accession(accession_dupe_set, eval_id="TKT_003")

    # Check these fields for specific values.
    tkt.check_type(type_set, True, eval_id="TKT_004")
    tkt.check_description_field(description_field_set, True, eval_id="TKT_005")
    tkt.check_run_mode(run_mode_set, True, eval_id="TKT_006")

    # For these fields, simply check that they are not empty.
    tkt.check_phage_id(null_set, False, eval_id="TKT_007")
    tkt.check_host_genus(null_set, False, eval_id="TKT_008")
    tkt.check_cluster(null_set, False, eval_id="TKT_009")
    tkt.check_annotation_status(null_set, False, eval_id="TKT_010")
    tkt.check_annotation_author(null_set, False, eval_id="TKT_011")
    tkt.check_retrieve_record(null_set, False, eval_id="TKT_012")

    # No need to evaluate the Accession and Subcluster fields
    # since they may or may not be populated.

    # Check if certain combinations of fields make sense.
    tkt.check_compatible_type_and_annotation_status(eval_id="TKT_013")

# TODO this may no longer be needed.
def check_phagesdb_genome(gnm, null_set):
    """Check a Genome object for specific errors when it has been
    parsed from PhagesDB data in preparation for completing import tickets.
    """
    gnm.check_id(null_set, False)
    gnm.check_name(null_set, False)
    gnm.check_host_genus(null_set, False)
    gnm.check_cluster(null_set, False)
    gnm.check_subcluster(null_set, False)
    gnm.check_accession(null_set, False)
    gnm.check_filename(null_set, False)
    gnm.check_sequence(null_set, False)


def check_genome_for_import(gnm, tkt, null_set=set(), phage_id_set=set(),
                           seq_set=set(), host_set=set(),
                           cluster_set=set(), subcluster_set=set(),
                           accession_set=set()):
    """Check a Genome object for errors.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt: A pdm_utils Ticket object.
    :type tkt: Ticket
    :param null_set: A set of values representing empty or null data.
    :type null_set: set
    :param phage_id_set: A set PhageIDs.
    :type phage_id_set: set
    :param seq_set: A set of genome sequences.
    :type seq_set: set
    :param host_set: A set of host genera.
    :type host_set: set
    :param cluster_set: A set of clusters.
    :type cluster_set: set
    :param subcluster_set: A set of subclusters.
    :type subcluster_set: set
    :param accession_set: A set of accessions.
    :type accession_set: set
    """

    if tkt.type == "add":
        gnm.check_id(phage_id_set | null_set, False, eval_id="GNM_001")
        gnm.check_name(phage_id_set | null_set, False, eval_id="GNM_002")
        gnm.check_sequence(seq_set | null_set, False, eval_id="GNM_003")

        # It is unusual for 'final' status genomes to be on 'add' tickets.
        gnm.check_annotation_status(check_set=set(["final"]), expect=False,
                                    eval_id="GNM_004")

        # If the genome is being added, and if it has an accession,
        # no other genome is expected to have an identical accession.
        # If the genome is being replaced, and if it has an accession,
        # the prior version of the genome may or may not have had
        # accession data, so no need to check for 'replace' tickets.
        if gnm.accession != "":
            gnm.check_accession(check_set=accession_set, expect=False,
                                eval_id="GNM_005")

    # 'replace' ticket checks.
    else:
        gnm.check_id(phage_id_set, True, eval_id="GNM_006")
        gnm.check_name(phage_id_set, True, eval_id="GNM_007")
        gnm.check_sequence(seq_set, True, eval_id="GNM_008")

        # It is unusual for 'draft' status genomes to be on 'replace' tickets.
        gnm.check_annotation_status(check_set=set(["draft"]), expect=False,
                                    eval_id="GNM_009")
    gnm.check_annotation_status(check_set=constants.ANNOTATION_STATUS_SET,
                                expect=True, eval_id="GNM_010")
    gnm.check_annotation_author(check_set=constants.ANNOTATION_AUTHOR_SET,
                                eval_id="GNM_011")
    gnm.check_retrieve_record(check_set=constants.RETRIEVE_RECORD_SET,
                              eval_id="GNM_012")
    gnm.check_cluster(cluster_set, True, eval_id="GNM_013")
    gnm.check_subcluster(subcluster_set, True, eval_id="GNM_014")
    gnm.check_subcluster_structure(eval_id="GNM_015")
    gnm.check_cluster_structure(eval_id="GNM_016")
    gnm.check_compatible_cluster_and_subcluster(eval_id="GNM_017")
    if tkt.eval_flags["check_seq"]:
        gnm.check_nucleotides(check_set=constants.DNA_ALPHABET,
                              eval_id="GNM_018")
    gnm.check_compatible_status_and_accession(eval_id="GNM_019")
    gnm.check_compatible_status_and_descriptions(eval_id="GNM_020")
    if tkt.eval_flags["check_id_typo"]:
        gnm.check_description_name(eval_id="GNM_021")
        gnm.check_source_name(eval_id="GNM_022")
        gnm.check_organism_name(eval_id="GNM_023")
    if tkt.eval_flags["check_host_typo"]:
        gnm.check_host_genus(host_set, True, eval_id="GNM_024")
        gnm.check_description_host_genus(eval_id="GNM_025")
        gnm.check_source_host_genus(eval_id="GNM_026")
        gnm.check_organism_host_genus(eval_id="GNM_027")
    if tkt.eval_flags["check_author"]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET,
                              eval_id="GNM_028")
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                                 expect=False, eval_id="GNM_029")
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET, expect=False,
                              eval_id="GNM_030")

    gnm.check_cds_feature_tally(eval_id="GNM_031")
    gnm.check_feature_ids(cds_ftr=True, trna_ftr=True, tmrna=True,
                          eval_id="GNM_032")

    # TODO confirm that these check_value_flag() are needed here.
    # Currently all "copy_data_to/from" functions run the check method
    # to throw an error if not all data was copied.
    gnm.set_value_flag("retrieve")
    gnm.check_value_flag(eval_id="GNM_033")
    gnm.set_value_flag("retain")
    gnm.check_value_flag(eval_id="GNM_034")


def check_source_for_import(src_ftr, check_id_typo=True, check_host_typo=True):
    """Check a Source object for errors."""
    if check_id_typo:
        src_ftr.check_organism_name(eval_id="SRC_001")
    if check_host_typo:
        src_ftr.check_organism_host_genus(eval_id="SRC_002")
        src_ftr.check_host_host_genus(eval_id="SRC_003")
        src_ftr.check_lab_host_host_genus(eval_id="SRC_004")


def check_cds_for_import(cds_ftr, check_locus_tag=True,
                         check_gene=True, check_description=True,
                         check_description_field=True):
    """Check a Cds object for errors."""
    cds_ftr.check_amino_acids(eval_id="CDS_001")
    cds_ftr.check_translation(eval_id="CDS_002")
    cds_ftr.check_translation_length(eval_id="CDS_003")
    cds_ftr.check_translation_table(eval_id="CDS_004")
    cds_ftr.check_coordinates(eval_id="CDS_005")
    cds_ftr.check_strand(eval_id="CDS_006")

    # These evaluations vary by genome type, stage of import, etc.
    if check_locus_tag:
        cds_ftr.check_locus_tag_present(eval_id="CDS_007")
        cds_ftr.check_locus_tag_structure(eval_id="CDS_008")
    if check_gene:
        cds_ftr.check_gene_present(eval_id="CDS_009")
        cds_ftr.check_gene_structure(eval_id="CDS_010")
    if check_locus_tag and check_gene:
        cds_ftr.check_compatible_gene_and_locus_tag(eval_id="CDS_011")
    if check_description:
        cds_ftr.check_generic_data(eval_id="CDS_012")
        cds_ftr.check_valid_description(eval_id="CDS_013")
    if check_description_field:
        cds_ftr.check_description_field(eval_id="CDS_014")


def compare_genomes(genome_pair, check_replace=True):
    """Compare two genomes to identify discrepancies."""
    genome_pair.compare_genome_sequence(eval_id="GP_001")
    genome_pair.compare_genome_length(eval_id="GP_002")
    genome_pair.compare_cluster(eval_id="GP_003")
    genome_pair.compare_subcluster(eval_id="GP_004")
    genome_pair.compare_accession(eval_id="GP_005")
    genome_pair.compare_host_genus(eval_id="GP_006")
    genome_pair.compare_author(eval_id="GP_007")
    if check_replace:
        genome_pair.compare_annotation_status("type","phamerator",
            "flat_file","draft","final", eval_id="GP_008")


# TODO implement.
# TODO unit test.
def check_trna_for_import(trna_obj):
    """Check a TrnaFeature object for errors."""

    pass
