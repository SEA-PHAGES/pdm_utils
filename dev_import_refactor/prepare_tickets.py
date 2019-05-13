"""Parse import table and parse into ticket objects.
"""










def main(ticket_file):




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


    # Some data may need to be retrieved from PhagesDB.
    ticket_list = retrieve_online_data(ticket_list)



    # Each ticket should be complete, now that data from PhagesDB has been
    # retrieved. Validate each ticket by checking each field in the ticket
    # that it is populated correctly.

    #TODO not sure if I should pass a list of valid types to this function.
    for ticket in ticket_list:
        validate(ticket)



    # Now that individual tickets have been validated,
    # validate the entire group of tickets.

    #TODO this should return information
    validate_tickets(ticket_list)


    return ticket_list


###
