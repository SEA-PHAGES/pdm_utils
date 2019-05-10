




def parse_import_tickets(list_of_data):

	list_of_tickets = []
	for input_row in list_of_data:

		#Verify the row of information has the correct number of fields to parse.
		if len(input_row) != 11:

			#TODO error handling
			write_out(output_file,"\nRow in import table is not formatted correctly: " + str(input_row))
			table_errors += 1
			continue


		ticket = ImportTicket()
		ticket.type = input_row[0] #Import action
		ticket.primary_phage_id = input_row[1] #New PhageID
		ticket.host = input_row[2] #Host
		ticket.cluster = input_row[3] #Cluster
		ticket.subcluster = input_row[4] #Subcluster
		ticket.status = input_row[5] #Status
		ticket.description_field = input_row[7] #Feature field
		ticket.accession = input_row[8] #Accession
		ticket.annotation_author = input_row[6] #AnnotationAuthor
		ticket.secondary_phage_id = input_row[10] #PhageID to be removed
		ticket.run_mode = input_row[9] #Run mode

		list_of_tickets.append(ticket)

		return(list_of_tickets)













#TODO need to complete
# Verify there are no duplicate actions

def validate_tickets(list_of_tickets):


    # Initialize all calculated attributes:
    add_set = set()
    remove_set = set()
    action_add_set = set()
    action_remove_set = set()
    action_add_remove_set = set()
    import_accession_set = set()



    #Create each set and do initial checks for duplications.
    #If the Add name or Remove name is "none", skip that because there are
    #probably duplicates of those.
    for ticket in list_of_tickets:
    	current_add = (ticket.primary_phage_id,)
    	current_remove = (ticket.secondary_phage_id,)
    	current_action_add = (ticket.type,ticket.primary_phage_id)
    	current_action_remove = (ticket.type,ticket.secondary_phage_id)
    	current_action_add_remove = (ticket.type,ticket.primary_phage_id,ticket.secondary_phage_id)


    	#First check the one-field and two-field combinations
    	if current_add[0] != "none":
    		if current_add in add_set:
    			print ticket.primary_phage_id + " appears to be involved in more than one step."

                #TODO error handling
                table_errors += question("\nError: %s is duplicated" % str(current_add))

    		else:
    			add_set.add(current_add)

    		if current_action_add in action_add_set:

                #TODO error handling
    			write_out(output_file,"\nError: %s is duplicated" % str(current_action_add))
    			table_errors += 1

    		else:
    			action_add_set.add(current_action_add)

    	if current_remove[0] != "none":

    		if current_remove in remove_set:

                #TODO error handling
    			print ticket.secondary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_remove))

    		else:
    			remove_set.add(current_remove)

    		if current_action_remove in action_remove_set:

                #TODO error handling
    			write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
    			table_errors += 1

    		else:
    			action_remove_set.add(current_action_remove)

    	#Now check the three-field combinations
    	if current_action_add_remove in action_add_remove_set:

            #TODO error handling
    		write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
    		table_errors += 1

    	else:
    		action_add_remove_set.add(current_action_add_remove)


    #Once the sets are created, also check if genomes to be removed are
    #found in the Add field and vice versa.
    for ticket in list_of_tickets:
    	current_add = (ticket.primary_phage_id,)
    	current_remove = (ticket.secondary_phage_id,)

    	#If the phage name is not replacing itself, the Add name is not expected
    	#to be in the Remove set and vice versa.
    	if current_add != current_remove:
    		if (current_add in remove_set and current_add != "none"):

                #TODO error handling
    			print ticket.primary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_add))


    		if (current_remove in add_set and current_remove != "none"):

                #TODO error handling
    			print ticket.secondary_phage_id + " appears to be involved in more than one step."
    			table_errors += question("\nError: %s is duplicated" % str(current_remove))




    #Verify there are no duplicate accessions in the import table
    for ticket in list_of_tickets:

    	if ticket.accession != "none":
    		if ticket.accession in import_accession_set:

                #TODO error handling
    			write_out(output_file,"\nError: Accession %s is duplicated in the import table." %ticket.accession)
    			table_errors += 1
    		else:
    			import_accession_set.add(ticket.accession)




	#TODO return proper information
	return pass
