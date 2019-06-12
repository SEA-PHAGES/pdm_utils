



def check_update_tickets(list_of_update_objects):


    # TODO complete function. It should create SQL statements to update data.
    pass
    return list_of_update_objects




def check_remove_tickets(list_of_remove_objects):

    # TODO complete function. It should create SQL statements to remove data.
    pass
    return list_of_remove_objects









def check_add_replace_tickets(list_of_matched_objects):



    # Flat file specific evaluations
    for matched_object in list_of_matched_objects:

        genome = matched_object.matched_genomes_dict["flat_file"]

        #TODO need to return the genome object?
        evaluate_genome(genome)



        # Compare flat file genome to ticket data
        compare_ticket_to_flat_file(matched_object)



        # Compare flat file genome to phamerator genome
        compare_ticket_to_phamerator(matched_object)



        # TODO other types of comparisons?





    # TODO complete this function.
    pass


    #TODO need to return object?
    return list_of_matched_objects












#TODO work on this function
# Cds object now contains a method to reset the primary description based
# on a user-selected choice.
#If other CDS fields contain descriptions, they can be chosen to
#replace the default import_cds_qualifier descriptions.
#Then provide option to verify changes.
#This block is skipped if user selects to do so.
def check_description_field_choice():

    if ignore_description_field_check != 'yes':

        changed = ""
        if (import_cds_qualifier != "product" and feature_product_tally > 0):
           print "\nThere are %s CDS products found." % feature_product_tally
           change_descriptions()

           if question("\nCDS products will be used for phage %s in file %s." % (phageName,filename)) == 1:
                for feature in all_features_data_list:
                    feature[9] = feature[10]
                changed = "product"

        if (import_cds_qualifier != "function" and feature_function_tally > 0):
            print "\nThere are %s CDS functions found." % feature_function_tally
            change_descriptions()

            if question("\nCDS functions will be used for phage %s in file %s." % (phageName,filename)) == 1:
                for feature in all_features_data_list:
                    feature[9] = feature[11]
                changed = "function"
        if (import_cds_qualifier != "note" and feature_note_tally > 0):

            print "\nThere are %s CDS notes found." % feature_note_tally
            change_descriptions()

            if question("\nCDS notes will be used for phage %s in file %s." % (phageName,filename)) == 1:
                for feature in all_features_data_list:
                    feature[9] = feature[12]
                changed = "note"

        if changed != "":
            record_warnings += 1
            write_out(output_file,"\nWarning: CDS descriptions only from the %s field will be retained." % changed)
            record_errors += question("\nError: problem with CDS descriptions of file %s." % filename)


###
