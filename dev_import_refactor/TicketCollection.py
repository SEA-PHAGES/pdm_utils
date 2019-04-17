"""Represents a structure to contain groups of tickets for importing
genomes into Phamerator."""






class TicketCollection:

    # Initialize all attributes:
    def __init__(self):


        # Initialize all non-calculated attributes:
        self.list_of_tickets = []




        # Initialize all calculated attributes:
        self.add_set = set()
        self.remove_set = set()
        self.action_add_set = set()
        self.action_remove_set = set()
        self.action_add_remove_set = set()
        self.import_accession_set = set()
        self.list_of_update_tickets = []
        self.list_of_remove_tickets = []
        self.list_of_add_replace_tickets = []







    def check_tickets():


        #Now that all rows have been added to the list, verify there are no duplicate actions

        #Create each set and do initial checks for duplications.
        #If the Add name or Remove name is "none", skip that because there are
        #probably duplicates of those.
        for ticket in self.list_of_tickets:
        	current_add = (ticket.primary_phage_id,)
        	current_remove = (ticket.secondary_phage_id,)
        	current_action_add = (ticket.type,ticket.primary_phage_id)
        	current_action_remove = (ticket.type,ticket.secondary_phage_id)
        	current_action_add_remove = (ticket.type,ticket.primary_phage_id,ticket.secondary_phage_id)


        	#First check the one-field and two-field combinations
        	if current_add[0] != "none":
        		if current_add in self.add_set:
        			print ticket.primary_phage_id + " appears to be involved in more than one step."

                    #TODO error handling
                    table_errors += question("\nError: %s is duplicated" % str(current_add))

        		else:
        			self.add_set.add(current_add)

        		if current_action_add in self.action_add_set:

                    #TODO error handling
        			write_out(output_file,"\nError: %s is duplicated" % str(current_action_add))
        			table_errors += 1

        		else:
        			self.action_add_set.add(current_action_add)

        	if current_remove[0] != "none":

        		if current_remove in self.remove_set:

                    #TODO error handling
        			print ticket.secondary_phage_id + " appears to be involved in more than one step."
        			table_errors += question("\nError: %s is duplicated" % str(current_remove))

        		else:
        			self.remove_set.add(current_remove)

        		if current_action_remove in self.action_remove_set:

                    #TODO error handling
        			write_out(output_file,"\nError: %s is duplicated" % str(current_action_remove))
        			table_errors += 1

        		else:
        			self.action_remove_set.add(current_action_remove)

        	#Now check the three-field combinations
        	if current_action_add_remove in self.action_add_remove_set:

                #TODO error handling
        		write_out(output_file,"\nError: %s is duplicated" % str(current_action_add_remove))
        		table_errors += 1

        	else:
        		self.action_add_remove_set.add(current_action_add_remove)


        #Once the sets are created, also check if genomes to be removed are
        #found in the Add field and vice versa.
        for ticket in self.list_of_tickets:
        	current_add = (ticket.primary_phage_id,)
        	current_remove = (ticket.secondary_phage_id,)

        	#If the phage name is not replacing itself, the Add name is not expected
        	#to be in the Remove set and vice versa.
        	if current_add != current_remove:
        		if (current_add in self.remove_set and current_add != "none"):

                    #TODO error handling
        			print ticket.primary_phage_id + " appears to be involved in more than one step."
        			table_errors += question("\nError: %s is duplicated" % str(current_add))


        		if (current_remove in self.add_set and current_remove != "none"):

                    #TODO error handling
        			print ticket.secondary_phage_id + " appears to be involved in more than one step."
        			table_errors += question("\nError: %s is duplicated" % str(current_remove))




        #Verify there are no duplicate accessions in the import table
        for ticket in self.list_of_tickets:

        	if ticket.accession != "none":
        		if ticket.accession in self.import_accession_set:

                    #TODO error handling
        			write_out(output_file,"\nError: Accession %s is duplicated in the import table." %ticket.accession)
        			table_errors += 1
        		else:
        			self.import_accession_set.add(ticket.accession)










    def split_tickets(self):

        #Create separate lists of ticket based on the indicated action: update, add/replace, remove
        for ticket in self.list_of_tickets:
        	if ticket.type == "update":
        		self.list_of_update_tickets.append(ticket)
        	elif ticket.type == "remove":
        		self.list_of_remove_tickets.append(ticket)
        	elif (ticket.type == "add" or ticket.type == "replace"):
        		self.list_of_add_replace_tickets.append(ticket)
        	else:
        		write_out(output_file,"\nError: during parsing of actions.")
        		table_errors += 1














	#Count # of ticket types
	#TODO not sure when I should re-implement this
	if row[0] in action_set:
		if row[0] == "add":
			add_total += 1
		elif row[0] == "remove":
			remove_total += 1
		elif row[0] == "replace":
			replace_total += 1
		elif row[0] == "update":
			update_total += 1
		else:
			pass




###
