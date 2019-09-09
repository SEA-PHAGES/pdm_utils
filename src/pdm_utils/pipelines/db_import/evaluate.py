"""Groups of evaluation/check functions."""


from pdm_utils.constants import constants




#
# def check_import_ticket_structure(
#         tkt,
#         type_set=constants.IMPORT_TICKET_TYPE_SET,
#         annotation_status_set=constants.ANNOTATION_STATUS_SET,
#         description_field_set=constants.DESCRIPTION_FIELD_SET,
#         annotation_author_set=constants.ANNOTATION_AUTHOR_SET,
#         run_mode_set=constants.RUN_MODE_SET,
#         null_set=constants.EMPTY_SET):
#     """Evaluate a ticket to confirm it is structured appropriately.
#     The assumptions for how each field is populated varies depending on
#     the type of ticket."""
#
#     # This function simply evaluates whether there is data in the
#     # appropriate ticket attributes given the type of ticket.
#     # It confirms that ticket attributes 'type', 'run_mode', and
#     # 'description_field' are populated with specific values.
#     # But it does not evaluate the quality of the data itself for
#     # the other fields, since those are genome-specific fields and
#     # can be checked within Genome objects.
#
#     #null_set = set(["none"])
#
#     # This is the only evaluation that is not dependent on the ticket type.
#     tkt.check_type(type_set, True)
#     tkt.check_phage_id(null_set, False)
#     tkt.check_host_genus(null_set, False)
#     tkt.check_cluster(null_set, False)
#     tkt.check_annotation_status(annotation_status_set)
#     tkt.check_description_field(description_field_set)
#     tkt.check_annotation_author(annotation_author_set)
#     tkt.check_run_mode(run_mode_set)
#
#     # No need to evaluate the Accession and Subcluster fields
#     # since they may or may not be populated.
#



# TODO this may no longer be needed.
def check_ticket_structure(tkt,type_set,null_set, run_mode_set):
    """Evaluate a ticket to confirm it is structured appropriately.
    The assumptions for how each field is populated varies depending on
    the type of ticket."""

    # This function simply evaluates whether there is data in the
    # appropriate ticket attributes given the type of ticket.
    # It confirms that ticket attributes 'type', 'run_mode', and
    # 'description_field' are populated with specific values.
    # But it does not evaluate the quality of the data itself for
    # the other fields, since those are genome-specific fields and
    # can be checked within Genome objects.

    #null_set = set(["none"])

    # This is the only evaluation that is not dependent on the ticket type.
    tkt.check_type(type_set, True)

    if (tkt.type == "add" or tkt.type == "replace"):
        tkt.check_phage_id(null_set, False)
        tkt.check_host_genus(null_set, False)
        tkt.check_cluster(null_set, False)
        tkt.check_annotation_status(null_set, False)
        tkt.check_description_field(null_set, False)
        tkt.check_annotation_author(null_set, False)
        tkt.check_run_mode(run_mode_set, True)

        # No need to evaluate the Accession and Subcluster fields
        # since they may or may not be populated.

        # TODO no more secondary_phage_id attribute, so this may need
        # to evaluate something else.
        # if tkt.type == "replace":
        #     tkt.check_secondary_phage_id(null_set, False)
        # else:
        #     tkt.check_secondary_phage_id(null_set, True)


    # TODO this may be deleted.
    # TODO unit test.
    elif tkt.type == "update":
        tkt.check_phage_id(null_set, False)
        tkt.check_host_genus(null_set, False)
        tkt.check_cluster(null_set, False)
        tkt.check_annotation_status(null_set, False)
        tkt.check_description_field(null_set, False)
        tkt.check_annotation_author(null_set, False)
        tkt.check_run_mode(null_set, True)

        # No need to evaluate the Accession and Subcluster fields
        # since they may or may not be populated.


    # TODO this may be deleted.
    # TODO unit test.
    elif tkt.type == "remove":

        # Everything except the primary phage_id field should be 'none'
        tkt.check_phage_id(null_set, False)
        tkt.check_host_genus(null_set, True)
        tkt.check_subcluster(null_set, True)
        tkt.check_cluster(null_set, True)
        tkt.check_annotation_status(null_set, True)
        tkt.check_description_field(null_set, True)
        tkt.check_accession(null_set, True)
        tkt.check_annotation_author(null_set, True)
        tkt.check_run_mode(null_set, True)

    else:
        pass






def check_phagesdb_genome(gnm, null_set):
    """Check a Genome object for specific errors when it has been
    parsed from PhagesDB data in preparation for completing import tickets."""

    gnm.check_id(null_set, False)
    gnm.check_name(null_set, False)
    gnm.check_host_genus(null_set, False)
    gnm.check_cluster(null_set, False)
    gnm.check_subcluster(null_set, False)
    gnm.check_accession(null_set, False)
    gnm.check_filename(null_set, False)
    gnm.check_sequence(null_set, False)







def check_source_for_import(src_ftr, check_id_typo=True, check_host_typo=True):
    """Check a Source object for errors."""
    if check_id_typo:
        src_ftr.check_organism_name()
    if check_host_typo:
        src_ftr.check_organism_host_genus()
        src_ftr.check_host_host_genus()
        src_ftr.check_lab_host_host_genus()


def check_cds_for_import(cds_ftr, check_locus_tag=True,
                         check_gene=True, check_description=True,
                         check_description_field=True):
    """Check a Cds object for errors."""
    cds_ftr.check_amino_acids()
    cds_ftr.check_translation()
    cds_ftr.check_translation_length()
    cds_ftr.check_translation_table()
    cds_ftr.check_coordinates()
    cds_ftr.check_strand()

    # These evaluations vary by genome type, stage of import, etc.
    if check_locus_tag:
        cds_ftr.check_locus_tag_present()
        cds_ftr.check_locus_tag_structure()
    if check_gene:
        cds_ftr.check_gene_present()
        cds_ftr.check_gene_structure()
    if check_locus_tag and check_gene:
        cds_ftr.check_compatible_gene_and_locus_tag()
    if check_description:
        cds_ftr.check_generic_data()
        cds_ftr.check_valid_description()
    if check_description_field:
        cds_ftr.check_description_field()



def compare_genomes(genome_pair, check_replace=True):
    """Compare two genomes to identify discrepancies."""
    genome_pair.compare_genome_sequence()
    genome_pair.compare_genome_length()
    genome_pair.compare_cluster()
    genome_pair.compare_subcluster()
    genome_pair.compare_accession()
    genome_pair.compare_host_genus()
    genome_pair.compare_author()

    if check_replace:
        genome_pair.compare_annotation_status("type","phamerator",
            "flat_file","draft","final")






















# TODO unit test below....


# TODO implement.
# TODO unit test.
def check_bundle_for_import(bndl):
    """Check a Bundle for errors."""

    tkt = bndl.ticket


    # First, evaluate whether all genomes have been successfully grouped,
    # and whether all genomes have been paired, as expected.
    # Based on the ticket type, there are expected to be certain
    # types of genomes and pairs of genomes in the bundle.

    if tkt.type == "add" or tkt.type == "replace":


        bndl.check_genome_dict("add")
        bndl.check_genome_dict("flat_file")
        bndl.check_genome_pair_dict("flat_file_add")

        # There may or may not be data retrieved from PhagesDB.
        tkt.set_value_flag("retrieve")
        if tkt._value_flag:
            bndl.check_genome_dict("phagesdb")
            try:
                check_phagesdb_genome(bndl.genome_dict["phagesdb"])
            except:
                pass
            try:
                bndl.check_genome_pair_dict("flat_file_phagesdb")
            except:
                pass

    # TODO this may need to be moved elsewhere.
    tkt.check_compatible_type_and_annotation_status()

    if tkt.type == "replace":

        # There should be a "phamerator" genome.
        bndl.check_genome_dict("phamerator")

        # There may or may not be a genome_pair to retain some data.
        tkt.set_value_flag("retain")
        if tkt._value_flag:
            bndl.check_genome_pair_dict("add_phamerator")

        # There should be a genome_pair between the current phamerator
        # genome and the new flat_file genome.
        bndl.check_genome_pair_dict("flat_file_phamerator")
        try:
            compare_genomes_for_replace(
                bndl.genome_pair_dict["flat_file_phamerator"])
        except:
            pass


    # Second, evaluate each genome or pair of genomes as needed.
    check_genome_to_import(bndl.genome_dict["flat_file"], tkt.type)








# TODO implement.
# TODO unit test.
def check_genome_to_import(gnm, tkt, null_set, phage_id_set,
                           seq_set, host_set, cluster_set, subcluster_set):
    """Check a Genome object for errors."""

    if tkt.type == "add":
        gnm.check_id(phage_id_set | null_set, False)
        gnm.check_name(phage_id_set | null_set, False)
        gnm.check_sequence(seq_set | null_set, False)

    # 'replace' ticket checks.
    else:
        gnm.check_id(phage_id_set, True)
        gnm.check_name(phage_id_set, True)
        gnm.check_sequence(seq_set, True)


    # TODO if these values are set at the command line, these checks
    # may no longer be needed.
    # tkt = tkt.check_description_field(description_field_set)
    # tkt = tkt.check_run_mode(run_mode_set)

    gnm.check_annotation_status(check_set=constants.ANNOTATION_STATUS_SET,
                                   expect=True)



    # TODO This may no longer be needed, if cluster data not used
    # during genome import.
    gnm.check_cluster(cluster_set, True)
    gnm.check_subcluster(subcluster_set, True)
    gnm.check_subcluster_structure()
    gnm.check_cluster_structure()
    gnm.check_compatible_cluster_and_subcluster()

    # TODO This may no longer be needed, if this database field is removed.
    # gnm.check_annotation_qc()

    # TODO This may no longer be needed.
    # gnm.check_filename()

    # TODO not sure if this is needed.
    # gnm.check_accession(accession_set, False)


    gnm.check_annotation_author()
    gnm.check_retrieve_record()



    if tkt.run_mode[check_seq]:
        gnm.check_nucleotides(dna_alphabet_set=constants.DNA_ALPHABET)
    gnm.check_compatible_status_and_accession()
    gnm.check_compatible_status_and_descriptions()

    if tkt.run_mode[check_id_typo]:
        gnm.check_description_name()
        gnm.check_source_name()
        gnm.check_organism_name()

    if tkt.run_mode[check_host_typo]:
        gnm.check_host_genus(host_set, True)
        gnm.check_description_host_genus()
        gnm.check_source_host_genus()
        gnm.check_organism_host_genus()


    if tkt.run_mode[check_author]:
        if gnm.annotation_author == 1:
            gnm.check_authors(check_set=constants.AUTHOR_SET)
            gnm.check_authors(check_set=set(["lastname", "firstname"]),
                                 expect=False)
        else:
            gnm.check_authors(check_set=constants.AUTHOR_SET, expect=False)



    gnm.check_cds_feature_tally()
    gnm.check_feature_ids(cds_ftr=True, trna_ftr=True, tmrna=True)


    # TODO not sure if these are needed now that check_feature_ids()
    # is implemented.
    # gnm.check_cds_start_end_ids()
    # gnm.check_cds_end_strand_ids()

    # TODO confirm that these check_value_flag() are needed here.
    gnm.set_value_flag("retrieve")
    gnm.check_value_flag()
    gnm.set_value_flag("retain")
    gnm.check_value_flag()

    # Check all CDS features
    index1 = 0
    while index1 < len(gnm.cds_features):
        check_cds_for_import(gnm.cds_features[index1])
        index1 += 1

    # Check all tRNA features
    if tkt.run_mode[check_trna]:
        index2 = 0
        while index2 < len(gnm.trna_features):
            check_trna_for_import(gnm.trna_features[index2])
            index2 += 1

    # Check all Source features
    index3 = 0
    while index3 < len(gnm.source_features):
        check_source_for_import(gnm.source_features[index3])
        index3 += 1






# TODO implement.
# TODO unit test.
def check_trna_for_import(trna_obj):
    """Check a TrnaFeature object for errors."""

    pass














def check_replace_tickets(bndl):
    """Check several aspects about a genome only if it is being replaced."""

    if len(bndl.genome_pair_dict.keys()) == 0:

        # TODO throw an error - there should be a matched genome_pair object
        # since this is a check_replace function.

        #
        #
        # # If the genome to be added is not spelled the same as the genome
        # # to be removed, the new genome needs to have a unique name.
        # if self.phage_id != self.secondary_phage_id:
        #     tkt.check_phage_id(phage_id_set, False)
        #
        # # No need to evaluate the following fields:
        # # Accession = it will either be an accession or it will be "none"
        # # Subcluster = it will either be a Subcluster or it will be "none"

        pass
    else:
        for key in bndl.genome_pair_dict.keys():
            compare_genomes(bndl.genome_pair_dict[key])


    pass














# TODO complete function. It should create SQL statements to update data.
# TODO implement.
# TODO unit test.
def check_update_tickets(list_of_update_objects):
    """."""
    pass
    # index = 0
    # while index < len(list_of_update_objects):
    #
    #     bndl = list_of_update_objects[index]
    #
    #     if len(bndl.genome_pairs_dict.keys()) == 0:
    #         # TODO throw an error if there is no matched Phamerator genome?
    #         pass
    #
    #     for key in bndl.genome_pairs_dict.keys():
    #
    #         genome_pair = bndl.genome_pairs_dict[key]
    #
    #         # TODO check for conflicting hosts. It is not common to
    #         # change hosts.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicting status. It is not common to
    #         # change status unless the current phamerator is draft status.
    #         # The only common change is from draft to final.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicting accession.
    #         # It is not common to change from real accession to another
    #         # real accession. But it is common to change from 'none' to
    #         # real accession.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicint authorship.
    #         # It is not common to change authorships.
    #         genome_pair.check_xyz()
    #
    #     index += 1




# TODO complete function. It should create SQL statements to remove data.
# TODO implement.
# TODO unit test.
def check_remove_tickets(list_of_remove_objects, genome_type):
    """."""
    pass
    #
    #
    # index = 0
    # while index < len(list_of_remove_objects):
    #
    #     bndl = list_of_update_objects[index]
    #
    #     try:
    #         gnm = bndl.genomes_dict[genome_type]
    #     except:
    #         # TODO throw an error if there is no matched Phamerator genome?
    #         continue
    #
    #
    #     # TODO list of evaluations for remove ticket.
    #
    #     # TODO check the status of the removing genome.
    #     # It is not common to remove anything but a 'draft' genome.
    #
    #
    #     index += 1








































# TODO implement.
# TODO unit test.
# Cds object now contains a method to reset the primary description based
# on a user-selected choice.
#If other CDS fields contain descriptions, they can be chosen to
#replace the default import_cds_qualifier descriptions.
#Then provide option to verify changes.
#This block is skipped if user selects to do so.
# def check_description_field_choice():
#
#     if ignore_description_field_check != 'yes':
#
#         changed = ""
#         if (import_cds_qualifier != "product" and feature_product_tally > 0):
#            print "\nThere are %s CDS products found." % feature_product_tally
#            change_descriptions()
#
#            if question("\nCDS products will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[10]
#                 changed = "product"
#
#         if (import_cds_qualifier != "function" and feature_function_tally > 0):
#             print "\nThere are %s CDS functions found." % feature_function_tally
#             change_descriptions()
#
#             if question("\nCDS functions will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[11]
#                 changed = "function"
#         if (import_cds_qualifier != "note" and feature_note_tally > 0):
#
#             print "\nThere are %s CDS notes found." % feature_note_tally
#             change_descriptions()
#
#             if question("\nCDS notes will be used for phage %s in file %s." % (phageName,filename)) == 1:
#                 for feature in all_features_data_list:
#                     feature[9] = feature[12]
#                 changed = "note"
#
#         if changed != "":
#             record_warnings += 1
#             write_out(output_file,"\nWarning: CDS descriptions only from the %s field will be retained." % changed)
#             record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)


###
