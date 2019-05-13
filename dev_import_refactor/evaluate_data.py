



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
###
