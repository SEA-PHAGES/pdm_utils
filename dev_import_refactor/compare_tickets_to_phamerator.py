"""Parse import table and parse into ticket objects.
"""








#TODO below: old python2 code that needs to be re-factored. It requires
# that tickets are matched to Phamerator data.


#Check to see if there are any inconsistencies with the
#update data compared to current phamerator data
for genome_data in update_data_list:

	#Initialize variable
	matched_phamerator_data = ''

	#Now that the Draft suffix is no longer appended to the import ticket,
	#this is less complications with matching to Phamerator PhageIDs
	try:
		matched_phamerator_data = phamerator_data_dict[genome_data[1]]
	except:
		write_out(output_file,"\nError: unable to retrieve phamerator data for %s." %genome_data[1])
		table_errors += 1
		continue

	#Host data check
	if genome_data[2] != matched_phamerator_data[2]:

		print "\n\nThere is conflicting host data for genome %s" % genome_data[1]
		print "Phamerator host: %s" % matched_phamerator_data[2]
		print "Import ticket host: %s" % genome_data[2]
		print "The new host data will be imported."
		table_errors += question("\nError: incorrect host data for %s." % genome_data[1])


	#Status data check
	if genome_data[4] != matched_phamerator_data[4]:

		#It is not common to change from 'gbk' or 'final' to anything else
		if matched_phamerator_data[4] != 'draft':

			print "\n\nThere is conflicting status data for genome %s" % genome_data[1]
			print "Phamerator status: %s" % matched_phamerator_data[4]
			print "Import ticket status: %s" % genome_data[4]
			print "The new status data will be imported."
			table_errors += question("\nError: incorrect status data for %s." % genome_data[1])

		#It is common to change status from 'draft' to 'final', but not anything else
		elif genome_data[4] != "final":

			print "\n\nThere is conflicting status data for genome %s" % genome_data[1]
			print "Phamerator status: %s" % matched_phamerator_data[4]
			print "Import ticket status: %s" % genome_data[4]
			print "The new status data will be imported."
			table_errors += question("\nError: incorrect status data for %s." % genome_data[1])


	#Accession data check
	if genome_data[7] == "none" and matched_phamerator_data[7] != "none":

		print "\n\nThere is conflicting accession data for genome %s" % genome_data[1]
		print "Phamerator accession: %s" % matched_phamerator_data[7]
		print "Import ticket accession: %s" % genome_data[7]
		print "The new accession data will be imported."
		table_errors += question("\nError: incorrect accession data for %s." % genome_data[1])

	elif genome_data[7] != "none" and matched_phamerator_data[7] != "none" and genome_data[7] != matched_phamerator_data[7]:

		print "\n\nThere is conflicting accession data for genome %s" % genome_data[1]
		print "Phamerator accession: %s" % matched_phamerator_data[7]
		print "Import ticket accession: %s" % genome_data[7]
		print "The new accession data will be imported."
		table_errors += question("\nError: incorrect accession data for %s." % genome_data[1])

	#Cluster, Subcluster check = no need to check this, as this data may be
	#more frequently updated than other fields.

	#Author
	#It is not common for authorship to change
	if genome_data[9] != matched_phamerator_data[9]:

		print "\n\nThere is conflicting author data for genome %s" % genome_data[1]
		print "Phamerator author: %s" % author_dictionary[matched_phamerator_data[9]]
		print "Import ticket author: %s" % author_dictionary[genome_data[9]]
		print "The new author data will be imported."
		table_errors += question("\nError: incorrect author data for %s." % genome_data[1])


#Check to see if genomes to be removed are the correct status
for genome_data in remove_data_list:

	#The next QC check relies on the PhageID being present in the Phamerator.
	#If it is not, this error will have already been identified, and a table_error
	#will already have been added. However, the code will crash without
	#here without adding an exception. So no need to add another table_error
	#in the except clause.
	try:
		matched_phamerator_data = phamerator_data_dict[genome_data[6]]
		if matched_phamerator_data[4] != "draft":
			print "The genome %s to be removed is currently %s status." % \
					(genome_data[6],matched_phamerator_data[4])
			table_errors += question("\nError: %s is %s status and should not be removed." % \
					(genome_data[6],matched_phamerator_data[4]))
	except:
		pass


#If no errors encountered, print list of action items to be
#implemented and continue. Otherwise, exit the script.
if table_errors == 0:
	write_out(output_file,"\nImport table verified with 0 errors.")
	write_out(output_file,"\nList of all actions to be implemented:")
	for genome_data in genome_data_list:
		write_out(output_file,"\n" + str(genome_data))
	raw_input("\nPress ENTER to proceed to next import stage.")

else:
	write_out(output_file,"\n%s error(s) encountered with import file.\nNo changes have been made to the database." % table_errors)
	write_out(output_file,"\nExiting import script.")
	output_file.close()
	sys.exit(1)
