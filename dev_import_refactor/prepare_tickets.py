"""Parse import table and parse into ticket objects.
"""



#Receive filename.


#TODO complete

#Validate filename.




# Retrieve import data
#List of ticket data.


#Retrieve import info from indicated import table file and read all lines into a list and verify contents are correctly populated.
#0 = Type of database action to be performed (add, remove, replace, update)
#1 = New PhageID that will be added to database
#2 = Host of new phage
#3 = Cluster of new phage (singletons should be reported as "singleton")
#4 = Subcluster of new phage (no subcluster should be reported as "none")
#5 = Annotation status of new phage
#6 = Annotation author of the new phage
#7 = Feature field containing gene descriptions of new phage
#8 = Accession
#9 = Run mode
#10 = PhageID of genome to be removed from the database



write_out(output_file,"\n\n\n\nRetrieving import info from table in file...")

file_object = open(updateFile,'r')
file_reader = csv.reader(file_object)
import_table_data_list = []
for input_row in file_reader:
    import_table_data_list.append(input_row)
file_object.close()














#TODO not sure if I need these counters any more.
# table_errors = 0
# add_total = 0
# remove_total = 0
# replace_total = 0
# update_total = 0
# run_mode_custom_total = 0




#Parse list of data and construct tickets.

#Convert data from import file into ticket objects
ticket_list = parse_import_tickets(import_table_data_list)



#Verify all data is cased appropriately.
for ticket_list in ticket_list:
	ticket.check_case()



# For each ticket, retrieve any data from PhagesDB genome if necessary.
for ticket in ticket_list:

    if (ticket.host == "retrieve" or \
    	ticket.cluster == "retrieve" or \
    	ticket.subcluster == "retrieve" or \
    	ticket.accession == "retrieve"):

    	phagesdb_genome = Genome()

    	#TODO make sure api_prefix and api_suffix variables are set
    	phage_url = api_prefix + ticket.primary_phage_id + api_suffix


    	try:
    		online_data_json = urllib.urlopen(phage_url)
    		online_data_dict = json.loads(online_data_json.read())

            #Returns a genome object
            phagesdb_genome = parse_phagesdb_data(phagesdb_genome,online_data_dict)

    	except:
    		online_data_json = ""
    		online_data_dict = {}

            #TODO handle error better
    		write_out(output_file,"\nError: unable to retrieve Host, Cluster, Subcluster, or Accession data for phage %s from phagesdb." %row[1])


    	if ticket.host == "retrieve":
    		ticket.host = phagesdb_genome.host
    	if ticket.cluster == "retrieve":
    		ticket.cluster = phagesdb_genome.cluster
    	if ticket.subcluster == "retrieve":
    		ticket.subcluster = phagesdb_genome.subcluster
    	if ticket.accession == "retrieve":
    		ticket.accession = phagesdb_genome.accession


# Each ticket should be complete, now that data from PhagesDB has been
# retrieved. Validate each ticket by checking each field in the ticket
# that it is populated correctly.

#TODO not sure if I should pass a list of valid types to this function.
for ticket in ticket_list:
    validate(ticket)



# Now that individual tickets have been validated,
# validate the entire group of tickets.

#TODO this should return information
validate_tickets(temp_list_of_tickets)












###
