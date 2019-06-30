#!/usr/bin/env python
#Phage Import Script
#Based on original script by Charles Bowman and Damani Brown on 20150216
#University of Pittsburgh
#Revamped by Travis Mavrich on 20160701
#Now it is designed to be a single script to upload genomes, update host/cluster/status data, remove genomes, and check integrity of files.
#WORKS FOR .gb, .gbf, .gbk, or .txt files




#Flow of the import process
#1 Import python modules, set up variables and functions
#2 Retrieve current database information
#3 Retrieve import information from file
#4 Update genome information for genomes already in the database
#5 Remove genomes from database that have no replacements
#6 Add or replace genomes

#The strategy for handling MySQL errors:
#The database is accessed in Steps 2, 4, 5, and 6. Connections are opened and closed at each Step. This ensures that any unforeseen errors encountered elsewhere in the script
#that force the script to exit do not cause ROLLBACK errors or connection errors.









#Create output file to store successful actions implemented
success_action_file = '%s_successful_import_table_actions.csv' % date
success_action_file_handle = open(os.path.join(phageListDir,success_folder,success_action_file),"w")
success_action_file_writer = csv.writer(success_action_file_handle)



#Update actions implemented
#Prepare to update fields that are not accompanied by adding or replacing a genome.
write_out(output_file,"\n\n\n\nUpdating Host, Cluster, Subcluster, and Status fields...")
updated = 0





















#Add and replace actions implemented
#Now that the data list does not contain update-only or remove-only data,
#create a dictionary. This serves to verify that all genomes
#to be imported are unique, as well as to
#be able to quickly retrieve information based on the genbank file that is opened.
#Key = PhageID
#Value = Genome data list:
add_replace_data_dict = {}
for genome_data in add_replace_data_list:

	#Verify there are no duplicate PhageIDs. This was checked before in a
	#slightly different way when looking for duplicate actions.
	if genome_data[1] not in add_replace_data_dict:
		add_replace_data_dict[genome_data[1]] = genome_data

	else:
		write_out(output_file,"\nError: problem creating genome data dictionary. PhageID %s already in dictionary. Check input table format." % genome_data[1])
		write_out(output_file,"\nNo add/replace actions have been implemented.")
		output_file.close()
		success_action_file_handle.close()
		sys.exit(1)

















#Iterate over each genbank-formatted genome file in the directory and parse all the needed data
write_out(output_file,"\nA total of %s file(s) containing Genbank-formatted records will be parsed." % len(genbank_files))
for filename in genbank_files:

	file_tally += 1
	write_out(output_file,"\n\nProcessing file %s: %s" % (file_tally,filename))

	basename = filename.split('.')[0]

	#The file is parsed to grab the header information
	for seq_record in SeqIO.parse(os.path.join(phageListDir,filename), "genbank"):

		matched_by_basename = ''
		add_replace_statements = []
		record_errors = 0
		record_warnings = 0
		geneID_set = set()
		phage_data_list = []



		#Create a list to hold summary info on the genome record:
		#Record Name
		#Record ID
		#Record Definition
		#Record Source
		#Record Organism
		#Organism qualifier of Source feature
		#List of CDS info:
			#Assigned geneID (how the geneID will look when uploaded to the database)
			#Locus tag
			#Product
			#Note
			#Function
			#Translation table
			#Translation
			#Processed Description (how the descriptions will look when uploaded to the database)
		record_summary_header = [["Header Field","Data"]]



		#Sequence, length, and GC
		# try:
			# phageSeq = seq_record.seq.upper()
			# seqLength = len(phageSeq)
			# seqGC = 100 * (float(phageSeq.count('G')) + float(phageSeq.count('C'))) / float(seqLength)

		# except:
		# 	write_out(output_file,"\nError: problem retrieving DNA sequence information in file %s. This file will not be processed." % filename)
		# 	record_errors += 1
		# 	failed_genome_files.append(filename)
		# 	script_warnings += record_warnings
		# 	script_errors += record_errors
		# 	raw_input("\nPress ENTER to proceed to next file.")
		# 	continue

		# #Check DNA sequence for possible errors
		# nucleotide_set = set(phageSeq)
		# nucleotide_error_set = nucleotide_set - dna_alphabet_set
		# if len(nucleotide_error_set) > 0:
		# 	record_warnings += 1
		# 	write_out(output_file,"\nWarning: phage %s contains unexpected nucleotide(s): %s" % (phageName,str(nucleotide_error_set)))
		# 	for element in nucleotide_error_set:
		# 		print "\nThere are %s unexpected %s nucleotides in %s." % (phageSeq.count(element),element,phageName)
		# 	record_errors += question("\nError: problem with DNA sequence in phage %s." % phageName)




		#File header fields are retrieved to be able to check phageName and HostStrain typos
		#The Accession field, with the appended version number, is stored as the record.id
		#The Locus name at the top of the file is stored as the record.name
		#The base accession number, without the version, is stored in the 'accession' annotation list
		# try:
		# 	record_name = str(seq_record.name)
		# except:
		# 	record_name = ""
		# 	print "\nRecord does not have record Locus information."
		# 	record_errors += question("\nError: problem with header info of file %s." % filename)


		# try:
		# 	record_id = str(seq_record.id)
		# except:
		# 	record_id = ""
		# 	print "\nRecord does not have record ID information."
		# 	record_errors += question("\nError: problem with header info of file %s." % filename)


		# try:
		# 	record_def = str(seq_record.description)
		# except:
		# 	record_def = ""
		# 	print "\nRecord does not have record definition information."
		# 	record_errors += question("\nError: problem with header info of file %s." % filename)


		# try:
		# 	record_source = str(seq_record.annotations["source"])
		# except:
		# 	record_source = ""
		# 	print "\nRecord does not have record source information."
		# 	record_errors += question("\nError: problem with header info of file %s." % filename)


		# #Date of the record
		# try:
		# 	seq_record_date = seq_record.annotations["date"]
		# 	seq_record_date = datetime.strptime(seq_record_date,'%d-%b-%Y')
        #
		# except:
		# 	write_out(output_file,"\nError: problem retrieving date in file %s. This file will not be processed." % filename)
		# 	record_errors += 1
		# 	failed_genome_files.append(filename)
		# 	script_warnings += record_warnings
		# 	script_errors += record_errors
		# 	raw_input("\nPress ENTER to proceed to next file.")
		# 	continue


		# #PhageNotes
		# try:
		# 	phageNotes = str(seq_record.annotations["comment"])
		# except:
		# 	phageNotes = ""


		# #Author list
		# try:
		# 	#The retrieved authors can be stored in multiple Reference elements
		# 	record_references_author_list = []
		# 	for reference in seq_record.annotations["references"]:
		# 		if reference.authors != "":
		# 			record_references_author_list.append(reference.authors)
		# 	if len(record_references_author_list) > 0:
		# 		record_author_string = ";".join(record_references_author_list)
		# 	else:
		# 		record_author_string = ""
		# except:
		# 	record_author_string = ""

		#Accession
		#This initiates this variable. Since different things affect this variable
		#depending on whether it is an add or replace action, it is easiest to
		#initiate it in advance to avoid throwing an error.
		accession_to_upload = ""

        #
		# try:
		# 	#There may be a list of accessions associated with this file.
		# 	#I think the first accession in the list is the most recent.
		# 	#Discard the version suffix if it is present in the
		# 	#Accession field (it might not be present).
		# 	#If an Accession is present, the script will interpret this to mean
		# 	#the genome record is derived from NCBI. An empty Accession will
		# 	#be interpreted as a manual annotation not retrieved from NCBI.
		# 	parsed_accession = seq_record.annotations["accessions"][0]
		# 	parsed_accession = parsed_accession.split('.')[0]
		# except:
		# 	parsed_accession = "none"


		#Old code used to match up to the import ticket data
		# if use_basename == "yes":
		#     matchedData = add_replace_data_dict.pop(basename,"error")
		# else:
		#     matchedData = add_replace_data_dict.pop(phageName,"error")




		#Since the Genbank record and the import ticket hasn't been matched yet
		#it is not clear whether it should be matched by filename or PhageID
		#If the basename is the same as the PhageID, it doesn't matter.
		#Otherwise it first checks by PhageID, then the basename.
		#Once matched, it confirms later that the matching strategy indicated in the
		#import ticket's run mode was the same that was actually used.
		if basename == phageName:
			matchedData = add_replace_data_dict.pop(basename,"error")
			matched_by_basename = 'both'

		elif phageName in add_replace_data_dict.keys():
			matchedData = add_replace_data_dict.pop(phageName,"error")
			matched_by_basename = 'no'

		elif basename in add_replace_data_dict.keys():
			matchedData = add_replace_data_dict.pop(basename,"error")
			matched_by_basename = 'yes'

		else:
			matchedData = "error"
			matched_by_basename = 'neither'


		if matchedData != "error":
			write_out(output_file,"\nPreparing: " + str(matchedData))
			import_action = matchedData[0]
			import_host = matchedData[2]
			import_cluster = matchedData[3]
			import_status = matchedData[4]
			import_cds_qualifier = matchedData[5]
			import_genome_replace = matchedData[6]
			import_accession = matchedData[7]
			import_subcluster = matchedData[8]
			import_author = matchedData[9]
			import_run_mode_dict = run_mode_options_dict[matchedData[10]]



			#Assign run mode parameters according to import ticket.
			use_basename = import_run_mode_dict['use_basename']
			custom_gene_id = import_run_mode_dict['custom_gene_id']
			ignore_gene_id_typo = import_run_mode_dict['ignore_gene_id_typo']
			ignore_description_field_check = import_run_mode_dict['ignore_description_field_check']
			ignore_replace_warning = import_run_mode_dict['ignore_replace_warning']
			ignore_trna_check = import_run_mode_dict['ignore_trna_check']
			ignore_locus_tag_import = import_run_mode_dict['ignore_locus_tag_import']
			ignore_phage_name_typos = import_run_mode_dict['ignore_phage_name_typos']
			ignore_host_typos = import_run_mode_dict['ignore_host_typos']
			ignore_generic_author = import_run_mode_dict['ignore_generic_author']
			ignore_description_check = import_run_mode_dict['ignore_description_check']


		else:
			write_out(output_file,"\nError: problem matching phage %s in file %s to genome data from table. This genome was not added. Check input table format." % (phageName,filename))
			record_errors += 1
			failed_genome_files.append(filename)
			script_warnings += record_warnings
			script_errors += record_errors
			raw_input("\nPress ENTER to proceed to next file.")
			continue



		#Now that the import ticket is matched, verify that the matching strategy
		#indicated in the tickets is equivalent to the strategy actually used
		if use_basename != matched_by_basename and matched_by_basename != 'both':
			record_errors += 1
			write_out(output_file,\
			"\nError: the genome in file %s was not matched to the import ticket as expected in the indicated run mode." \
			% (filename))



		#Check the database to make sure the genome sequence and name matches correctly.
		con = mdb.connect(mysqlhost, username, password, database)
		con.autocommit(False)
		cur = con.cursor()
		try:
			cur.execute("START TRANSACTION")
			cur.execute("""SELECT PhageID,status FROM phage WHERE Sequence = "%s" """ % phageSeq)
			query_results = cur.fetchall()
			cur.execute("SELECT GeneID,PhageID FROM gene")
			current_gene_data_tuples = cur.fetchall()
			cur.execute("COMMIT")
			cur.close()
			con.autocommit(True)


		except:
			success_action_file_handle.close()
			mdb_exit("\nError: retrieving genome information from database while processing file %s.\nNot all genbank files were processed." % filename)

		con.close()


		#Create a set of GeneIDs. If a genome will be replaced,
		#do not add those GeneIDs to the set.
		all_GeneID_set = set()
		for gene_tuple in current_gene_data_tuples:
			if (import_action == "replace" and gene_tuple[1] == import_genome_replace):
				continue
			all_GeneID_set.add(gene_tuple[0])



		#Cross-check the import action against the current state of
		#the database and create SQL statements

		#If adding a new genome, no genome sequence in database is expected
		#to match the current genome sequence
		if (import_action == "add" and len(query_results) > 0):
			record_errors += 1
			write_out(output_file,\
			"\nError: these genome(s) in the database currently contain the same genome sequence as %s: %s." \
			% (phageName,query_results))


		#If replacing a genome:
		elif import_action == "replace":

			#Retrieve current phamerator data (if it is a 'replace' ticket)
			try:
				matched_phamerator_data = phamerator_data_dict[import_genome_replace]
			except:
				matched_phamerator_data = "error"


			if matched_phamerator_data != "error":
				write_out(output_file,"\nRetrieved phamerator data for phage: %s" % phageName)

				phamerator_host = matched_phamerator_data[2]
				phamerator_cluster = matched_phamerator_data[5]
				phamerator_status = matched_phamerator_data[4]
				phamerator_datelastmod = matched_phamerator_data[6]
				phamerator_accession = matched_phamerator_data[7]
				phamerator_subcluster = matched_phamerator_data[8]
				phamerator_author = matched_phamerator_data[9]
				phamerator_annotation_qc = matched_phamerator_data[10]
				phamerator_retrieve_record = matched_phamerator_data[11]

			else:
				write_out(output_file,"\nError: problem matching phage %s in file %s to phamerator data. This genome was not added. Check input table format." % (phageName,filename))
				record_errors += 1
				failed_genome_files.append(filename)
				script_warnings += record_warnings
				script_errors += record_errors
				raw_input("\nPress ENTER to proceed to next file.")
				continue


			#Import and Phamerator host data check
			if import_host != phamerator_host:

				record_warnings += 1
				write_out(output_file,"\nWarning: There is conflicting host data for genome %s" % phageName)
				print "Phamerator host: %s" % phamerator_host
				print "Import ticket host: %s" % import_host
				print "The new host data will be imported."
				record_errors += question("\nError: incorrect host data for %s." % phageName)


			#Import and Phamerator status data check
			if import_status != phamerator_status:

				#It is not common to change from 'gbk' or 'final' to anything else
				if phamerator_status != 'draft':

					record_warnings += 1
					write_out(output_file,"\nWarning: There is conflicting status data for genome %s" % phageName)
					print "Phamerator status: %s" % phamerator_status
					print "Import ticket status: %s" % import_status
					print "The new status data will be imported."
					record_errors += question("\nError: incorrect status data for %s." % phageName)

				#It is common to change status from 'draft' to 'final', but not anything else
				elif import_status != "final":

					record_warnings += 1
					write_out(output_file,"\nWarning: There is conflicting status data for genome %s" % phageName)
					print "Phamerator status: %s" % phamerator_status
					print "Import ticket status: %s" % import_status
					print "The new status data will be imported."
					record_errors += question("\nError: incorrect status data for %s." % phageName)


			#Verify there are no accession conflicts.
			#Phamerator may or may not have accession data.
			#Import table may or may not have accession data.
			#Genbank-formatted file may or may not have accession data.

			accession_comparison_set = set()
			if phamerator_accession != "none":
				accession_comparison_set.add(phamerator_accession)
			if import_accession != "none":
				accession_comparison_set.add(import_accession)
			if parsed_accession != "none":
				accession_comparison_set.add(parsed_accession)


			if len(accession_comparison_set) == 1:
				accession_to_upload = list(accession_comparison_set)[0]
			elif len(accession_comparison_set) > 1:
				record_warnings += 1
				write_out(output_file,"\nWarning: There is conflicting accession data for genome %s" % phageName)
				print "Phamerator accession: %s" % phamerator_accession
				print "Import ticket accession: %s" % import_accession
				print "Parsed accession from file: %s" % parsed_accession
				print "If the parsed accession is not None, it will be imported."
				print "If the parsed accession is None, but the import ticket accession is not None, it will be imported."
				record_errors += question("\nError: incorrect accession data for %s." % phageName)

				if parsed_accession != "none":
					accession_to_upload = parsed_accession
				elif import_accession != "none":
					accession_to_upload = import_accession
				elif phamerator_accession != "none":
					accession_to_upload = phamerator_accession
				else:
					accession_to_upload = ""
			else:
				accession_to_upload = ""



			#Author check
			#It is not common for authorship to change
			if import_author != phamerator_author:

				record_warnings += 1
				write_out(output_file,"\nWarning: there is conflicting author data for genome %s." % phageName)
				print "Phamerator author: %s" % author_dictionary[phamerator_author]
				print "Import ticket author: %s" % author_dictionary[import_author]
				print "The new author data will be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)



			#Exactly one and only one genome in the database is
			#expected to have the same sequence.
			if len(query_results) > 1:
				record_errors += 1
				write_out(output_file,"\nError: the following genomes in the database currently contain the same genome sequence as %s: %s).\nUnable to perform replace action." % (phageName,query_results))

			elif len(query_results) == 0:

				record_warnings += 1
				write_out(output_file,"\nWarning: %s appears to be a different genome sequence than %s. These genomes do not match." % (phageName,import_genome_replace))
				print "The genome will still be replaced."
				record_errors += question("\nError: %s and %s have different genome sequences." % (phageName,import_genome_replace))

			elif len(query_results) == 1:

				#If the new genome is a Final or Gbk, this code block
				#alerts the user, unless this warning is turned off.
				if query_results[0][1].lower() != "draft" and ignore_replace_warning != 'yes':

					record_warnings += 1
					write_out(output_file,"\nWarning: The genome in the database with matching sequence, %s, is listed as %s status." % (query_results[0][0],query_results[0][1]))
					print "The genome will still be replaced."
					record_errors +=  question("\nError: the genome to be removed, %s, was incorrect status." %import_genome_replace)

				#The genome to be replaced does not match the genome
				#name in the database with the same sequence.
				if query_results[0][0] != import_genome_replace:
					write_out(output_file,"\nError: the genome to be removed, %s, does not match the genome name in the database, %s, that has the matching genome sequence to %s." % (import_genome_replace,query_results[0][0],phageName))
					record_errors += 1
			else:
				pass

			#Check to see if the date in the new record is more recent
			#than when the old record was uploaded into Phamerator
			#(stored in DateLastModified)
			if not seq_record_date > phamerator_datelastmod:
				record_warnings += 1
				write_out(output_file,"\nWarning: The date %s in file %s is not more recent than the Phamerator date %s." %(seq_record_date,filename,phamerator_datelastmod))
				print 'Despite it being an older record, the phage %s will continue to be imported.' % phageName
				record_errors +=  question("\nError: the date %s in file %s is not more recent than the Phamerator date %s." %(seq_record_date,filename,phamerator_datelastmod))

			#Create the DELETE command
			add_replace_statements.append("DELETE FROM phage WHERE PhageID = '" + import_genome_replace + "';")

		else:
			#At this point, if the genome is being added, simply assign the accession data
			#from the import table to the accession_to_upload variable. The assignment
			#logic is more complex if the genome is being replaced.
			if import_accession != "none":
				accession_to_upload = import_accession






		#Author list check
		#For annotation author = hatfull
		#If annotation status is draft, author field can be missing Hatfull
		#If annotation status is final, author field should have Hatfull
		#For annotation author = gbk, author should NOT be in either field
		pattern5 = re.compile("hatfull")
		search_result = pattern5.search(record_author_string.lower())


		#For Hatfull authored final annotations, Hatfull is an expected author
		if import_author == '1' and import_status == 'final':

			if search_result == None:
				record_warnings += 1
				write_out(output_file,"\nWarning: Graham Hatfull is not a listed author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)


		#For Hatfull authored draft annotations, it doesn't matter whether
		#there are authors, since they will be added at the Final stage.
		elif import_author == '1' and import_status == 'draft':
			pass

		#For non-Hatfull authored annotations of any kind (draft, final, gbk),
		#Graham may or may not be an author
		else:
			if search_result != None:
				record_warnings += 1
				write_out(output_file,"\nWarning: Graham Hatfull is a listed author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)





		#Check for generic author that sometimes gets added by DNA Master
		if ignore_generic_author != 'yes':
			pattern6 = re.compile("lastname")
			search_result6 = pattern6.search(record_author_string.lower())

			pattern7 = re.compile("firstname")
			search_result7 = pattern7.search(record_author_string.lower())


			if search_result6 != None or search_result7 != None:
				record_warnings += 1
				write_out(output_file,"\nWarning: the author list appears to contain a generic Lastname, Firstname author for genome %s" % phageName)
				print "The genome will continue to be imported."
				record_errors += question("\nError: incorrect author data for %s." % phageName)






		#Determine AnnotationQC and RetrieveRecord settings
		#AnnotationQC and RetrieveRecord settings can be carried over
		#from the previous settings in Phamerator under specific
		#circumstances.
		#This enables manual updates made to the database for these two fields
		#to be carried over as the genome is manually or automatically updated.

		#Since a change in AnnotatioQC or RetrieveRecord from default settings
		#is expected to be rare, these two fields can be updated without
		#requiring two additional fields in the import table.

		#For genome replacements, RetrieveRecord is determined based on
		#previous setting. AnnotationQC is determined from previous setting
		#UNLESS it is a Final genome replacing a Draft, in which the AnnotationQC
		#must be changed from 0 to 1.
		if import_action == "replace":

			ncbi_update_status = phamerator_retrieve_record

			#Only under the specific event of a Final replacing a Draft
			#should the AnnotationQC be set to 1. All other events should simply
			#carry it over from the previous setting.
			if phamerator_status == 'draft' and import_status == 'final':
				annotation_qc = '1'

			else:
				annotation_qc = phamerator_annotation_qc

		#For adding a new genome, AnnotationQC is determined by the status,
		#and the RetrieveRecord is determined by the author.
		else:

			if import_author == '1':
				ncbi_update_status = '1'
			else:
				ncbi_update_status = '0'

			if import_status == 'draft' or import_status == 'gbk':
				annotation_qc = '0'
			else:
				annotation_qc = '1'








		#Next each CDS feature is parsed from the file
		#The cdsCount will increment for each CDS processed, even if it does not pass the QC filters. This way, the genome still retains the info that another CDS was originally present.
		cdsCount = 0
		addCount = 0
		transl_table_set = set()
		missing_transl_table_tally = 0
		missing_locus_tag_tally = 0
		assigned_description_tally = 0
		feature_note_tally = 0
		feature_product_tally = 0
		feature_function_tally = 0
		feature_source_organism = ""
		feature_source_lab_host = ""
		feature_source_host = ""
		all_features_data_list = []
		all_coordinates_set = set()
		record_summary_cds = [["Locus Tag",\
								"Product",\
								"Function",\
								"Note",\
								"Translation Table",\
								"Translation",\
								"Assigned GeneID",\
								"Assigned Description"]]

		for feature in seq_record.features:






















###Below: error code from import script that needs to be implemented


		#Not sure where to put this now:



		#Some phage names are spelled differently in phagesdb and Phamerator
		#compared to the Genbank record. Convert the parsed name to the
		#expected name in these databases before proceeding.
		if phageName in phage_name_typo_dict.keys():
			incorrect_phageName = phageName
			phageName = phage_name_typo_dict[phageName]
			record_warnings += 1
			write_out(output_file,"\nWarning: parsed phage name %s converted to %s." \
				% (incorrect_phageName,phageName))




    	genome_object.compute_nucleotide_errors(dna_alphabet_set)




		#CDS errors from compare_datases.py
            #Compute other fields
            gene_object.compute_amino_acid_errors(protein_alphabet_set)
            gene_object.set_start_end_strand_id()
            gene_object.compute_boundary_error()
            gene_object.compute_description_error()


		#Genome errors after all features are parsed
	    genome_object.compute_cds_feature_errors()
	    genome_object.compute_ncbi_cds_feature_errors()








		#TODO decide when to re-assign gene_id if needed

			#If user selected customized GeneIDs, OR if there is no locus tag,
			#then create a concatenated GeneID
			if custom_gene_id == 'yes' or feature_locus_tag == '':
				if use_basename == 'yes':
					geneID = basename.upper() + "_" + str(cdsCount)
				else:
					geneID = phageName.upper() + "_" + str(cdsCount)
			else:
				geneID = feature_locus_tag

			#See if the geneID is already in the database
			duplicate = False
			if geneID in geneID_set:
				duplicate = True
				record_warnings += 1
				write_out(output_file,"\nWarning: there is a duplicate geneID %s in phage %s." % (geneID,phageName))

			elif geneID in all_GeneID_set:
				duplicate = True
				record_warnings += 1
				write_out(output_file,"\nWarning: there is a duplicate geneID %s in the current database." % geneID)
			else:
				geneID_set.add(geneID)

			#If there is a geneID duplication conflict, try to resolve it
			old_ID = geneID
			dupe_value = 1
			while duplicate == True and dupe_value < 99:

				geneID = old_ID + "_duplicateID" + str(dupe_value)
				record_warnings += 1
				#Check to see if the new geneID is found in the set of all geneIDs
				if (geneID not in geneID_set and geneID not in all_GeneID_set):
					duplicate = False
					write_out(output_file,"\nGeneID %s duplication has been automatically resolved by renaming ID to %s." % (old_ID,geneID))
					duplicate_answer = question("\nError: feature %s of %s is a duplicate geneID." % (old_ID,phageName))

					#If user indicates the feature with the new geneID should not be added, add to record_errors. Otherwise, assign new geneID to the geneID_set
					if duplicate_answer == 1:
						record_errors += 1
					else:
						geneID_set.add(geneID)
				dupe_value += 1

			#Once the while loop exits, check if the duplication was resolved
			if duplicate == True:
				record_warnings += 1
				write_out(output_file,"\nWarning: unable to resolve duplicate geneID %s conflict. This CDS will be skipped, but processing of the other genes will continue." % old_ID)
				record_errors += question("\nError: feature %s of phage %s is a duplicate geneID and cannot be renamed to %s." % (old_ID,phageName,geneID))
				continue


			#Name
			if (geneID.split('_')[-1].isdigit()):
				geneName = geneID.split('_')[-1]
			else:
				geneName = cdsCount




			#Orientation
			if feature.strand == 1:
				orientation = "Forward"
			elif feature.strand == -1:
				orientation = "Reverse"
			#ssRNA phages
			elif feature.strand is None:
				orientation = "Forward"
			else:
				record_warnings += 1
				write_out(output_file,"\nWarning: feature %s of %s does not have a common orientation. This CDS will be skipped, but processing of the other genes will continue." % (geneID,phageName))
				record_errors += question("\nError: feature %s of %s does not have correct orientation." % (geneID,phageName))
				continue

            #
			# #Now that start and stop have been parsed, check if coordinates are fuzzy or not
			# if (strStart.isdigit() and strStop.isdigit()):
			# 	startCoord = int(strStart)
			# 	stopCoord = int(strStop)
			# else:
			# 	record_warnings += 1
			# 	write_out(output_file,"\nWarning: gene %s start %s and stop %s are non-traditional coordinates. This CDS will be skipped, but processing of the other genes will continue." % (geneID,strStart,strStop))
			# 	record_errors += question("\nError: feature %s of %s does not have correct coordinates." % (geneID,phageName))
			# 	continue

			# #Test if there is a gene with the same coordinates already parsed.
			# coordinate_tuple = tuple([startCoord,stopCoord,orientation])
			# if coordinate_tuple not in all_coordinates_set:
			# 	all_coordinates_set.add(coordinate_tuple)
			# else:
			# 	record_warnings += 1
			# 	write_out(output_file,"\nWarning: multiple genes have coordinates %s. This is likely a gene feature duplication." % str(coordinate_tuple))
			# 	record_errors += question("\nError: gene coordinates %s are duplicated in this genome." % str(coordinate_tuple))



			# #Translation, Gene Length (via Translation)
			# try:
			# 	translation = feature.qualifiers["translation"][0].upper()
			# 	geneLen = (len(translation) * 3) + 3  #Add 3 for the stop codon...
			# except:
			# 	translation = ""
			# 	geneLen = 0
			# 	record_warnings += 1
			# 	write_out(output_file,"\nWarning: gene %s has no translation. This CDS will be skipped, but processing of the other genes will continue." % geneID)
			# 	record_errors += question("\nError: problem with %s translation in phage %s." % (geneID,phageName))
			# 	continue




			#TODO error to implement - different translation tables used?



			# #Check translation for possible errors
			# amino_acid_set = set(translation)
			# amino_acid_error_set = amino_acid_set - protein_alphabet_set
			# if len(amino_acid_error_set) > 0:
			# 	record_warnings += 1
			# 	write_out(output_file,"\nWarning: feature %s of %s appears to have unexpected amino acid(s)." % (geneID,phageName))
			# 	print "Unexpected amino acids: " + str(amino_acid_error_set)
			# 	record_errors += question("\nError: problem with %s translation in phage %s." % (geneID,phageName))

			#TODO compute description tally for feature, product, and note




			#Now assign the appropriate description info to the
			#assigned_description variable, as indicated from the import table.
			try:

				if import_cds_qualifier == "product":
					assigned_description = feature_product

				elif import_cds_qualifier == "function":
					assigned_description = feature_function

				elif import_cds_qualifier == "note":
					assigned_description = feature_note

				#This clause allows the user to specify an uncommon
				#feature qualifier to retrieve the gene description from.
				else:
					assigned_description = retrieve_description(feature,import_cds_qualifier)

			except:
				assigned_description = ""

			if assigned_description != "":
				assigned_description_tally += 1



			#Check to verify misc details about the assigned gene description field
			#This is distinct from the ignore_description_field_check, which
			#is used to determine which field (product, function, note) is parsed.

			if ignore_description_check != 'yes':
				if assigned_description[-1:] == '.':

					record_warnings += 1
					write_out(output_file,"\nWarning: the description for feature %s of %s contains a period." % (geneID,phageName))
					print assigned_description
					record_errors += question("\nError: problem with description for feature %s of %s." % (geneID,phageName))



		#Check to see if there are any CDS features processed. If not, then the genbank record does not have any called genes.
		#The record_summary_cds list contains the column headers, so at minimum, it is length == 1
		# if len(record_summary_cds) == 1:
		# 	print "\nNo CDS features were found in this record. The genome will still be added to the database."
		# 	record_errors += question("\nError: no CDS features found in %s." % filename)


		#See if there are any phage name typos in the header block
		if ignore_phage_name_typos != 'yes':


			pattern1 = re.compile('^' + phageName + '$')
			pattern2 = re.compile('^' + phageName)


			#REVIEW This QC check may not be necessary. It really doesn't matter
			#for Draft genomes, since the fields are auto-populated. For NCBI genomes,
			#the field should be an Accession. For Final genomes yet to be submitted
			#to NCBI, even if there is a typo, it should not be retained since this
			#field gets populated with the Accession.
			# if find_name(pattern1,record_name.split(' ')) == 0:
			#
			#     if record_name.split('.')[0] != parsed_accession:
			#         print "\nRecord name does not have the accession number or the identical phage name as found in the record organism field."
			#         record_errors += question("\nError: problem with header info of file %s." % filename)


			#Record definition QC
			#It can contain ", complete genome." or "." at the end,
			#so remove this before doing search.
			# if record_def[-1:].lower() == '.':
			# 	if record_def[-18:].lower() == ', complete genome.':
			# 		record_def_trimmed = record_def[:-18]
			# 	else:
			# 		record_def_trimmed = record_def[:-1]
			# else:
			# 	record_def_trimmed = record_def

            #
			# if find_name(pattern1,record_def_trimmed.split(' ')) == 0:
			# 	print "\nRecord definition does not have identical phage name as found in the record organism field."
			# 	record_errors += question("\nError: problem with header info of file %s." % filename)
            #
            #
			# #Record source QC
			# if find_name(pattern1,record_source.split(' ')) == 0:
			# 	print "\nRecord source does not have identical phage name as found in the record organism field."
			# 	record_errors += question("\nError: problem with header info of file %s." % filename)
            #
			# #Source feature organism QC
			# if find_name(pattern1,feature_source_organism.split(' ')) == 0:
			# 	print "\nSource feature organism does not have identical phage name as found in the record organism field."
			# 	record_errors += question("\nError: problem with header info of file %s." % filename)



		#See if there are any host name typos in the header block.
		#Skip this step if it is a Draft genome, because it won't correctly have this information.
		# if import_status != 'draft' and ignore_host_typos != 'yes':
		# 	import_host_trim = import_host
		# 	if import_host_trim == "Mycobacterium":
		# 		import_host_trim = import_host_trim[:-3]
        #
		# 	pattern3 = re.compile('^' + import_host_trim)
        #
        #
		# 	if (find_name(pattern3,record_def.split(' ')) == 0 and record_def.split(' ')[0].lower() not in host_ignore):
		# 		print "\nRecord definition does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)
        #
		# 	if (find_name(pattern3,record_source.split(' ')) == 0 and record_source.split(' ')[0].lower() not in host_ignore):
        #
		# 		print "\nRecord source does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)
        #
		# 	if (find_name(pattern3,record_organism.split(' ')) == 0 and record_organism.split(' ')[0].lower() not in host_ignore):
        #
		# 		print "\nRecord organism does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)
        #
		# 	if (find_name(pattern3,feature_source_organism.split(' ')) == 0 and feature_source_organism.split(' ')[0].lower() not in host_ignore):
        #
		# 		print "\nSource feature organism does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)
        #
		# 	#Host and Lab_Host data may not have been present, so skip if it is blank
		# 	if (feature_source_host != "" and find_name(pattern3,feature_source_host.split(' ')) == 0 and feature_source_host.split(' ')[0].lower() not in host_ignore):
        #
		# 		print "\nSource feature host does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)
        #
		# 	if (feature_source_lab_host != "" and find_name(pattern3,feature_source_lab_host.split(' ')) == 0 and feature_source_lab_host.split(' ')[0].lower() not in host_ignore):
        #
		# 		print "\nSource feature lab host does not appear to have same host data as found in import table."
		# 		record_errors += question("\nError: problem with header info of file %s." % filename)




		# #Check locus tag info:
		# if missing_locus_tag_tally > 0:
		# 	record_warnings += 1
		# 	write_out(output_file,"\nWarning: phage %s from file %s is missing %s CDS locus tag(s)." % (phageName, filename, missing_locus_tag_tally))
		# 	record_errors += question("\nError: problem with locus tags in file  %s." % filename)


		# #Check the phage name spelling in the locus tags.
		# pattern4 = re.compile(phageName.lower())
		# geneID_typo_tally = 0
		# geneID_typo_list = []
        #
		# if ignore_gene_id_typo != "yes":
		# 	for geneID in geneID_set:
        #
		# 	   search_result = pattern4.search(geneID.lower())
		# 	   if search_result == None:
		# 			geneID_typo_tally += 1
		# 			geneID_typo_list.append(geneID)
        #
		# 	if geneID_typo_tally > 0:
		# 		record_warnings += 1
		# 		write_out(output_file,"\nWarning: there are %s geneID(s) that do not have the identical phage name included." % geneID_typo_tally)
		# 		print geneID_typo_list
		# 		record_errors += question("\nError: problem with locus tags of file %s." % filename)



		# #Check all translation table info:
		# if len(transl_table_set) > 1:
		# 	write_out(output_file,"\nError: more than one translation table used in file %s." % filename)
		# 	record_errors += 1
        #
		# elif len(transl_table_set) == 1:
		# 	transl_table_list = list(transl_table_set)
		# 	if transl_table_list[0] != '11':
		# 		write_out(output_file,"\nThe translation table used for %s is: %s." % (phageName,transl_table_list[0]))
		# 		record_errors += question("\nError: phage %s does not use correct translation table." % phageName)
		# else:
		# 	pass
        #
		# if missing_transl_table_tally > 0:
		# 	record_warnings += 1
		# 	write_out(output_file,"\nWarning: there are %s genes with no translation table for phage %s." % (missing_transl_table_tally,phageName))
		# 	record_errors += question("\nError: phage %s has missing translation table information." % phageName)




# TODO:
# Now that the flat file data has been retrieved and
# FlatFile genome objects created, match them to ticket data

















#TODO compile all gene feature data for import

#Now that it has acquired all gene feature info, create list of
#gene data and append to list of all gene feature data
#0 = geneID
#1 = PhageID or basename
#2 = startCoord
#3 = stopCoord
#4 = geneLen
#5 = geneName
#6 = typeID
#7 = translation
#8 = orientation[0]
#9 = assigned_description
#10 = feature_product
#11 = feature_function
#12 = feature_note
#13 = feature_locus_tag
addCount+= 1

feature_data_list.append(geneID) #0


if use_basename == "yes":
	feature_data_list.append(basename) #1
else:
	feature_data_list.append(phageName) #1

feature_data_list.append(startCoord) #2
feature_data_list.append(stopCoord) #3
feature_data_list.append(geneLen) #4
feature_data_list.append(geneName) #5
feature_data_list.append(typeID) #6
feature_data_list.append(translation) #7
feature_data_list.append(orientation[0]) #8
feature_data_list.append(assigned_description) #9
feature_data_list.append(feature_product) #10
feature_data_list.append(feature_function) #11
feature_data_list.append(feature_note) #12

#If there is a parsed accession, the file is interpreted to have
#been retrieved from NCBI, which means the locus tags are the official
#locus tags. Otherwise, do not retain the locus tags in the record.
if ignore_locus_tag_import != 'yes':
	feature_data_list.append(feature_locus_tag) #13
else:
	feature_data_list.append('') #13

all_features_data_list.append(feature_data_list)








#If errors were encountered with the file parsing, do not add to the genome. Otherwise, proceed.
if record_errors == 0:
	write_out(output_file,"\nNo errors encountered while parsing: %s." % filename)
	diffCount = cdsCount - addCount
	if record_warnings > 0:
		write_out(output_file,"\nWarning summary: there are %s warning(s) with phage %s and %s CDS feature(s) were not added." % (record_warnings,phageName,diffCount))
else:
	write_out(output_file,"\n%s errors were encountered with file %s. %s genome was not added to the database." % (record_errors,filename,phageName))
	failed_genome_files.append(filename)
	failed_actions.append(matchedData)
	script_warnings += record_warnings
	script_errors += record_errors
	raw_input("\nPress ENTER to proceed to next file.")
	continue


#Execute SQL transactions
if run_type == "production":
	con = mdb.connect(mysqlhost, username, password, database)
	con.autocommit(False)
	cur = con.cursor()
	try:
		cur.execute("START TRANSACTION")
		for statement in add_replace_statements:
			cur.execute(statement)
			write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
		cur.execute("COMMIT")
		write_out(output_file,"\nAll add/replace statements for %s committed." % matchedData)
		cur.close()
		con.autocommit(True)

	except:
		success_action_file_handle.close()
		mdb_exit("\nError: problem importing the file %s with the following add/replace action: %s.\nNot all genbank files were processed." % (filename,matchedData))

	con.close()

else:
	write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)

#Now that genome has been successfully uploaded, proceed
try:
	shutil.move(os.path.join(phageListDir,filename),os.path.join(phageListDir,success_folder,filename))
except:
	print "Unable to move file %s to success file folder." % filename


#Add the action data to the success output file, update tally of total script warnings and errors, then proceed
add_replace_output_list = [matchedData[0],\
						matchedData[1],\
						matchedData[2],\
						matchedData[3],\
						matchedData[8],\
						matchedData[4],\
						author_dictionary[matchedData[9]],\
						matchedData[5],\
						matchedData[7],\
						matchedData[10],\
						matchedData[6]]


success_action_file_writer.writerow(add_replace_output_list)
script_warnings += record_warnings
script_errors += record_errors
print "Processing of %s is complete." %filename

















###Above: error code from import script that needs to be implemented


















































write_out(output_file,"\nAll files have been iterated through.")
raw_input("\nPress ENTER to proceed to next import stage.")


#Final verifications
write_out(output_file,"\n\n\n\nFinal checks...")
write_out(output_file,"\n\nAll update actions have been implemented.")
write_out(output_file,"\n\nAll remove actions have been implemented.")

#If there were failed actions during the genome upload stage, add it back to the add_replace_data_dictionary.
#The add_replace_data_dictionary now contains a list of failed actions as well as unaddressed actions (actions that did not have matching genome files)
if len(failed_actions) > 0:
	for genome_data in failed_actions:
		add_replace_data_dict[genome_data[1]] = genome_data


#Verify all add/replace actions from the import csv file were addressed.
if len(add_replace_data_dict) > 0:

	failed_action_file = '%s_failed_import_table_actions.csv' % date
	failed_action_file_handle = open(os.path.join(phageListDir,failed_folder,failed_action_file),"w")
	failed_action_file_writer = csv.writer(failed_action_file_handle)
	write_out(output_file,"\n\nThe following add/replace action(s) in the import table were NOT successfully implemented:")

	for key in add_replace_data_dict:

		failed_output_list = [add_replace_data_dict[key][0],\
								add_replace_data_dict[key][1],\
								add_replace_data_dict[key][2],\
								add_replace_data_dict[key][3],\
								add_replace_data_dict[key][8],\
								add_replace_data_dict[key][4],\
								author_dictionary[add_replace_data_dict[key][9]],\
								add_replace_data_dict[key][5],\
								add_replace_data_dict[key][7],\
								add_replace_data_dict[key][10],\
								add_replace_data_dict[key][6]]

		failed_action_file_writer.writerow(failed_output_list)
		write_out(output_file,"\n" + str(failed_output_list))

	failed_action_file_handle.close()

else:
	write_out(output_file,"\nAll add/replace actions have been implemented.")


#Verify that none of the genbank files failed to be uploaded.
if len(failed_genome_files) > 0:

	write_out(output_file,"\n\nThe following genbank files were NOT successfully processed:")
	for filename in failed_genome_files:
		write_out(output_file,"\n" + filename)

		try:
			shutil.move(os.path.join(phageListDir,filename),os.path.join(phageListDir,failed_folder,filename))
		except:
			print "Unable to move file %s to failed file folder." % filename


else:
	write_out(output_file,"\n\nAll genbank files were successfully processed.")



#Output the total number of warnings and errors encountered during the import process
write_out(output_file,"\n\nTotal number of warnings encountered: %s" % script_warnings)
write_out(output_file,"\nTotal number of errors encountered: %s" % script_errors)


#Retrieve the total number of genomes now in the updated database
try:
	con = mdb.connect(mysqlhost, username, password, database)
	con.autocommit(False)
	cur = con.cursor()
except:
	print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password.\nImport script was not completed."
	output_file.close()
	success_action_file_handle.close()
	sys.exit(1)

try:
	cur.execute("START TRANSACTION")
	cur.execute("SELECT count(*) FROM phage")
	final_tally = cur.fetchall()
	cur.execute("COMMIT")
	cur.close()
	con.autocommit(True)

except:
	output_file.close()
	success_action_file_handle.close()
	mdb_exit("\nUnable to access the database to retrieve genome information.\nImport script was not completed.")

write_out(output_file,"\n\nTotal phages in database after changes: " + str(final_tally[0][0]))


#Close script.
success_action_file_handle.close()
write_out(output_file,"\n\n\n\nImport script completed.")
output_file.close()
