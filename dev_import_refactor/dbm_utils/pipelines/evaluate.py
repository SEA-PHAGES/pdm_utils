"""Groups of evaluation/check functions."""




def check_ticket_structure(ticket, type_set, null_set, run_mode_set):
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
    ticket.check_type(type_set, True)

    if (ticket.type == "add" or ticket.type == "replace"):
        ticket.check_primary_phage_id(null_set, False)
        ticket.check_host(null_set, False)
        ticket.check_cluster(null_set, False)
        ticket.check_status(null_set, False)
        ticket.check_description_field(null_set, False)
        ticket.check_annotation_author(null_set, False)
        ticket.check_run_mode(run_mode_set, True)

        # No need to evaluate the Accession and Subcluster fields
        # since they may or may not be populated.

        if ticket.type == "replace":
            ticket.check_secondary_phage_id(null_set, False)
        else:
            ticket.check_secondary_phage_id(null_set, True)


    # TODO this may be deleted.
    # TODO unit test.
    elif ticket.type == "update":
        ticket.check_primary_phage_id(null_set, False)
        ticket.check_host(null_set, False)
        ticket.check_cluster(null_set, False)
        ticket.check_status(null_set, False)
        ticket.check_description_field(null_set, False)
        ticket.check_annotation_author(null_set, False)
        ticket.check_run_mode(null_set, True)
        ticket.check_secondary_phage_id(null_set, True)

        # No need to evaluate the Accession and Subcluster fields
        # since they may or may not be populated.


    # TODO this may be deleted.
    # TODO unit test.
    elif ticket.type == "remove":

        # Everything except the primary phage_id field should be 'none'
        ticket.check_primary_phage_id(null_set, False)
        ticket.check_secondary_phage_id(null_set, True)
        ticket.check_host(null_set, True)
        ticket.check_subcluster(null_set, True)
        ticket.check_cluster(null_set, True)
        ticket.check_status(null_set, True)
        ticket.check_description_field(null_set, True)
        ticket.check_accession(null_set, True)
        ticket.check_annotation_author(null_set, True)
        ticket.check_run_mode(null_set, True)

    else:
        pass







# TODO this function needs to be completely revamped. It should
# be implemented within the Genome object.
# TODO unit test.
# def check_ticket_structure(ticket,
#                             null_set,
#                             phage_id_set,
#                             host_set,
#                             cluster_set,
#                             status_set,
#                             author_set,
#                             description_field_set,
#                             run_mode_set):
#     """Evaluate a ticket to confirm it is structured appropriately.
#     The assumptions for how each field is populated varies depending on
#     the type of ticket."""
#
#     #null_set = set(["none"])
#
#     # This is the only evaluation that is not dependent on the ticket type.
#     ticket.check_type(type_set)
#
#     if ticket.type == "update":
#         ticket.check_primary_phage_id(phage_id_set)
#         ticket.check_host(host_set)
#         ticket.check_cluster(cluster_set)
#         ticket.check_subcluster_structure()
#         ticket.check_cluster_structure()
#         ticket.check_cluster_subcluster()
#         ticket.check_status(status_set)
#         ticket.check_description_field(null_set)
#         ticket.check_annotation_author(author_set)
#         ticket.check_run_mode(null_set)
#         ticket.check_secondary_phage_id(null_set)
#
#         # No need to evaluate the following fields:
#         # Accession = it will either be an accession or it will be "none"
#         # Subcluster = it will either be a Subcluster or it will be "none"
#
#
#     elif ticket.type == "add":
#         # TODO make sure it checks that the primary_phage_id
#         # is not 'none' as well.
#         ticket.check_primary_phage_id(phage_id_set, False)
#         ticket.check_secondary_phage_id(null_set)
#         ticket.check_host(host_set)
#         ticket.check_cluster(cluster_set)
#         ticket.check_subcluster_structure()
#         ticket.check_cluster_structure()
#         ticket.check_cluster_subcluster()
#         ticket.check_status(status_set)
#         ticket.check_description_field(description_field_set)
#         ticket.check_annotation_author(author_set)
#         ticket.check_run_mode(run_mode_set)
#
#         # No need to evaluate the following fields:
#         # Accession = it will either be an accession or it will be "none"
#         # Subcluster = it will either be a Subcluster or it will be "none"
#
#     elif ticket.type == "remove":
#
#         # Everything except the secondary phage_id field should be 'none'
#         ticket.check_primary_phage_id(null_set)
#         ticket.check_secondary_phage_id(phage_id_set)
#         ticket.check_host(null_set)
#         ticket.check_subcluster(null_set)
#         ticket.check_cluster(null_set)
#         ticket.check_status(null_set)
#         ticket.check_description_field(null_set)
#         ticket.check_accession(null_set)
#         ticket.check_annotation_author(null_set)
#         ticket.check_run_mode(null_set)
#
#     elif ticket.type == "replace":
#         ticket.check_secondary_phage_id(phage_id_set)
#         ticket.check_host(host_set)
#         ticket.check_cluster(cluster_set)
#         ticket.check_subcluster_structure()
#         ticket.check_cluster_structure()
#         ticket.check_cluster_subcluster()
#         ticket.check_status(status_set)
#         ticket.check_description_field(description_field_set)
#         ticket.check_annotation_author(author_set)
#         ticket.check_run_mode(run_mode_set)
#         ticket.check_primary_secondary_phage_ids()
#
#
#         # If the genome to be added is not spelled the same as the genome
#         # to be removed, the new genome needs to have a unique name.
#         if self.primary_phage_id != self.secondary_phage_id:
#             ticket.check_primary_phage_id(phage_id_set, False)
#
#         # No need to evaluate the following fields:
#         # Accession = it will either be an accession or it will be "none"
#         # Subcluster = it will either be a Subcluster or it will be "none"
#     else:
#         pass










# TODO implement.
# TODO unit test.
def check_add_replace_tickets(list_of_matched_objects, genome_type):

    index = 0
    while index < len(list_of_matched_objects):

        matched_object = list_of_matched_objects[index]

        # Since a DataGroup object can hold multiple genomes, the specific
        # genome that needs to be evaluated for import needs to be specified
        # by 'genome_type'.

        genome = matched_object.matched_genomes_dict[genome_type]
        check_genome(genome)

        # TODO is this needed?
        compare_ticket_to_flat_file(matched_object)

        if matched_object.ticket.type == "add":
            check_add_tickets(genome)
        elif matched_object.ticket.type == "replace":
            check_replace_tickets(matched_object)
        else:
            pass

        index += 1










# TODO implement.
# TODO unit test.
def check_genome(genome_obj):
    """Check a Genome object for errors."""

    # TODO decide how to implement alphabet
    genome_obj.check_nucleotides(alphabet = alphabet)
    genome_obj.check_status_accession()
    genome_obj.check_status_descriptions()
    genome_obj.check_record_description_phage_name()
    genome_obj.check_record_source_phage_name()
    genome_obj.check_record_organism_phage_name()
    genome_obj.check_record_description_host_name()
    genome_obj.check_record_source_host_name()
    genome_obj.check_record_organism_host_name()
    genome_obj.check_author()
    genome_obj.check_cds_feature_tally()


    # TODO decide how to evaluate duplicate feature coordinates.
    genome_obj.check_cds_start_end_ids()
    genome_obj.check_cds_end_strand_ids()


    # Check all CDS features
    index1 = 0
    while index1 < len(genome_obj.cds_features):
        check_cds(genome_obj.cds_features[index1])
        index1 += 1

    # Check all tRNA features
    index2 = 0
    while index2 < len(genome_obj.trna_features):
        check_trna(genome_obj.trna_features[index2])
        index2 += 1


    # Check all Source features
    index3 = 0
    while index3 < len(genome_obj.source_features):
        check_source(genome_obj.source_features[index3])
        index3 += 1





# TODO implement.
# TODO unit test.
def check_cds(cds_obj):
    """Check a CdsFeature object for errors."""

    # TODO decide how to implement alphabet
    cds_obj.check_amino_acids(alphabet = alphabet)
    cds_obj.check_translation_length()
    cds_obj.check_lengths()
    cds_obj.check_strand()
    cds_obj.check_boundaries()
    cds_obj.check_locus_tag_present()
    cds_obj.check_locus_tag_typo()
    cds_obj.check_description()
    cds_obj.check_translation_table_present()
    cds_obj.check_translation_table_typo()
    cds_obj.check_translation()





# TODO implement.
# TODO unit test.
def check_trna(trna_obj):
    """Check a TrnaFeature object for errors."""

    pass





# TODO implement.
# TODO unit test.
def check_source(src_obj):
    """Check a SourceFeature object for errors."""

    src_obj.check_organism_phage_name()
    src_obj.check_organism_host_name()
    src_obj.check_host_host_name()
    src_obj.check_lab_host_host_name()






def compare_genomes(genome_pair_obj):
    """Compare two genomes to identify discrepancies."""

    genome_pair_obj.compare_genome_sequence()
    genome_pair_obj.compare_genome_length()
    genome_pair_obj.compare_cluster()
    genome_pair_obj.compare_subcluster()
    genome_pair_obj.compare_accession()
    genome_pair_obj.compare_host()
    genome_pair_obj.compare_author()



    # TODO at this stage check the status of the genome. If it is a final,
    # and there is no other paired genome, it should throw an error. This was
    # moved from the ticket evaluation stage.
    # if self.status == "final":
    #     result6 = "The phage %s to be added is listed " + \
    #             "as Final status, but no Draft (or other) genome " + \
    #             " is listed to be removed."
    #     status6 = "error"
    # else:
    #     result6 = ""
    #     status6 = "correct"



def check_add_tickets(genome_obj):
    """Check several aspects about a genome only if it is being added."""

    # TODO determine other assumptions about what is/isn't in the database
    # if this is a new genome.

    # TODO verify there is no other identical genome sequence.
    # TODO verify there is no other identical phage id.

    pass




def check_replace_tickets(matched_object):
    """Check several aspects about a genome only if it is being replaced."""

    if len(matched_object.genome_pair_dict.keys()) == 0:
        # TODO throw an error - there should be a matched genome_pair object
        # since this is a check_replace function.
        pass
    else:
        for key in matched_object.genome_pair_dict.keys():
            compare_genomes(matched_object.genome_pair_dict[key])

    # TODO determine other assumptions about what is/isn't in the database
    # if this is a replacement.
    # TODO verify there is only one identical genome sequence.
    # TODO confirm that phage_ids of both genomes are the same.

    pass









# TODO complete this function.
# TODO unit test.
def compare_ticket_to_flat_file(matched_object):

    pass








# TODO complete function. It should create SQL statements to update data.
# TODO implement.
# TODO unit test.
def check_update_tickets(list_of_update_objects):
    """."""

    index = 0
    while index < len(list_of_update_objects):

        matched_object = list_of_update_objects[index]

        if len(matched_object.genome_pairs_dict.keys()) == 0:
            # TODO throw an error if there is no matched Phamerator genome?
            pass

        for key in matched_object.genome_pairs_dict.keys():

            genome_pair = matched_object.genome_pairs_dict[key]

            # TODO check for conflicting hosts. It is not common to
            # change hosts.
            genome_pair.check_xyz()

            # TODO check for conflicting status. It is not common to
            # change status unless the current phamerator is draft status.
            # The only common change is from draft to final.
            genome_pair.check_xyz()

            # TODO check for conflicting accession.
            # It is not common to change from real accession to another
            # real accession. But it is common to change from 'none' to
            # real accession.
            genome_pair.check_xyz()

            # TODO check for conflicint authorship.
            # It is not common to change authorships.
            genome_pair.check_xyz()

        index += 1




# TODO complete function. It should create SQL statements to remove data.
# TODO implement.
# TODO unit test.
def check_remove_tickets(list_of_remove_objects, genome_type):
    """."""



    index = 0
    while index < len(list_of_remove_objects):

        matched_object = list_of_update_objects[index]

        try:
            genome = matched_object.genomes_dict[genome_type]
        except:
            # TODO throw an error if there is no matched Phamerator genome?
            continue


        # TODO list of evaluations for remove ticket.

        # TODO check the status of the removing genome.
        # It is not common to remove anything but a 'draft' genome.


        index += 1








































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
