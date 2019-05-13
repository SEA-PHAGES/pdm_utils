






#TODO import appropriate classes

# Determine which files contain GenBank-formatted flat file genome information.
# Files in the directory may not have the correct extension,
# or they may contain 0, 1, or >1 parseable records.
# When SeqIO parses files, if there are 0 Genbank-formatted records,
# it does not throw an error, but simply moves on.
# This function iterates through each file in the user-indicated directory,
# and it checks how many actually are valid Genbank files.
# It's advantageous to do this first round of file iteration to identify
# any files with multiple records present, and that the file can be
# successfully read by Biopython SeqIO (if not, it could crash the script).
# Not sure how often there are multiple records present in non-SEA-PHAGES
# GenBank records, but it is a good idea to verify there is only one record
# per file before proceeding.
def validate_flat_file(filename):

	ADMISSIBLE_FILE_TYPES = set(["gb","gbf","gbk","txt"])
	valid = 0

	#If the file extension is not admissible, then skip. Otherwise, proceed
	if filename.split('.')[-1] not in ADMISSIBLE_FILE_TYPES:

		# TODO error handling
		write_out(output_file,"\nError: file %s does not have a valid file extension. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")
		continue

	#This try/except clause prevents the code from crashing if there
	#is a problem with a file that Biopython has trouble parsing.
	try:
		#Keep track of how many records Biopython parses
		parsed_records_tally = 0
		for seq_record in SeqIO.parse(os.path.join(phageListDir,filename), "genbank"):
			parsed_records_tally += 1

	except:

		# TODO error handling
		write_out(output_file,"\nError: Biopython is unable to parse file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")
		continue


	if parsed_records_tally == 0:

		# TODO error handling
		write_out(output_file,"\nError: Biopython was unable to parse any records from file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")

	elif parsed_records_tally > 1:

		# TODO error handling
		write_out(output_file,"\nError: Biopython found two records in file %s. This file will not be processed." % filename)
		script_errors += 1
		raw_input("\nPress ENTER to proceed to next file.")



	else:
		# The record has 1 parsed record
		valid = 1

	return(valid)








#TODO import appropriate classes
def parse_cds_feature(feature):

	gene_object = Cds()

    #Feature type
    gene_object.type_id = 'CDS'

    #Locus tag
    try:
        gene_object.set_locus_tag(feature.qualifiers['locus_tag'][0])
    except:
        gene_object.set_locus_tag("")


    #Orientation
    gene_object.set_strand_for_import(feature.strand)


    #Left and right boundary coordinates
    left_boundary,right_boundary,compound_parts = parse_coordinates(feature)
    gene_object.left_boundary = left_boundary
    gene_object.right_boundary = right_boundary
    gene_object.compound_parts = compound_parts

    #Translation
    try:
        gene_object.set_translation(feature.qualifiers['translation'][0])
    except:
        pass

    #TODO change this to nucleotide length?
    gene_object.gene_length = (len(translation) * 3) + 3  #Add 3 for the stop codon...




    #Translation table used
    try:
        gene_object.translation_table = feature.qualifiers["transl_table"][0]
    except:
        gene_object.translation_table = ""





	#HERE
    # Gene function, note, and product descriptions
    # The Cds object should be initialized with '' for each description
    # field.

    try:
        feature_product,feature_processed_product = \
        retrieve_description(feature.qualifiers['product'][0])
        gene_object.product_description = feature_product
        gene_object.processed_product_description = feature_processed_product

    except:
        pass


    try:
        feature_function,feature_processed_function = \
        retrieve_description(feature.qualifiers['function'][0])
        gene_object.function_description = feature_function
        gene_object.processed_function_description = feature_processed_function

    except:
        pass

    try:
        feature_note,feature_processed_note = \
        retrieve_description(feature.qualifiers['note'][0])
        gene_object.note_description = feature_note
        gene_object.processed_note_description = feature_processed_note

    except:
        pass



    #Gene number
    try:
        gene_object.gene_number = feature.qualifiers['gene'][0]
    except:
        pass


    return(gene_object)








# Parsing gene boundary coordinates can be tricky.
# There can be more than one set of coordinates if it is a
# compound location.

def parse_coordinates(feature):

    left_boundary = "0"
    right_boundary = "0"
    compound_parts = 0

    if str(feature.location)[:4] == 'join':

        compound_parts = len(feature.location.parts)

        # Skip this compound feature if it is comprised of more
        # than two features (too tricky to parse).
        if compound_parts <= 2:

            #Retrieve compound feature positions based on strand
            if feature.strand == 1:
                left_boundary = str(feature.location.parts[0].start)
                right_boundary = str(feature.location.parts[1].end)
            elif feature.strand == -1:
                left_boundary = str(feature.location.parts[1].start)
                right_boundary = str(feature.location.parts[0].end)

            #If strand is None, not sure how to parse it
            else:
                pass
        else:
            #TODO error handling - make note if unable to parse coordinates
            write_out(output_file,"\nWarning: gene %s is a compound feature that is unable to be parsed. This CDS will be skipped, but processing of the other genes will continue." % geneID)
            record_errors += question("\nError: unable to parse gene %s of phage %s." % (geneID,phageName))
            continue





    else:
        left_boundary = str(feature.location.start)
        right_boundary = str(feature.location.end)
        compound_parts = 1

    return(left_boundary,right_boundary,compound_parts)









def parse_organism_field(record_organism):

	phage_name = ""
	host_name = ""

	record_organism_list = record_organism.split(" ")


	# Get rid of the "Unclassified" element if present.
    if (len(record_organism_list) > 0 and \
	record_organism_list[-1] == "Unclassified."):

		record_organism_list.pop()


	#Now the last element is the phage name
	if len(record_organism_list) > 0:

	    phage_name = record_organism_list.pop()

	# Get rid of the "phage" element if present.
    if (len(record_organism_list) > 0 and \
	record_organism_list[-1] == "phage"):

		record_organism_list.pop()

	#Now the last element is the host name
	if len(record_organism_list) > 0:

		host_name = record_organism_list.pop()



	#TODO account for ", complete genome" text?

	#TODO account for merged names (e.g. "Mycobacteriophage")?

	return(phage_name,host_name)









# #TODO implement this function
#
# #Analayze tRNA features
# def parse_trna_feature(feature):
#     #TODO need to implement this function
#     #return(None)
#
#
# 	#Retrieve tRNA coordinates
# 	try:
#
# 		#Biopython converts coordinates to 0-index
# 		#Start(left) coordinates are 0-based inclusive (feature starts there)
# 		#Stop (right) coordinates are 0-based exclusive (feature stops 1bp prior to coordinate)
# 		tRNA_left = str(feature.location.start)
# 		tRNA_right = str(feature.location.end)
#
# 	except:
#
#         #TODO error handling
# 		write_out(output_file,"\nError: a tRNA has incorrect coordinates in phage %s."\
# 				% phageName)
# 		record_errors += 1
# 		continue
#
#
#
# 	#Retrieve top strand of tRNA feature. It is NOT necessarily
# 	#in the correct orientation
# 	tRNA_size = abs(tRNA_right - tRNA_left)
# 	tRNA_seq = phageSeq[tRNA_left:tRNA_right].upper()
#
#
#
#
#
# 	#Convert sequence to reverse complement if it is on bottom strand
# 	if feature.strand == 1:
# 		pass
# 	elif feature.strand == -1:
# 		tRNA_seq = tRNA_seq.reverse_complement()
# 	else:
#         #TODO error handling
# 		record_errors += 1
# 		write_out(output_file,"\Error: tRNA starting at %s does not have proper orientation in %s phage." \
# 				% (tRNA_left + 1,phageName))
# 		continue
#
#
#
#
#
# 	#Retrieve and check product
# 	try:
# 		tRNA_product = feature.qualifiers['product'][0].lower().strip()
#
# 	except:
#         #TODO error handling
# 		write_out(output_file,"\nError: tRNA starting at %s is missing product field in phage %s." \
# 			% (tRNA_left + 1,phageName))
# 		record_errors += 1
# 		tRNA_product = ''
#
#
# 	#Retrieve note
# 	#In the future, this field may need to be parsed in a similar
# 	#manner as the product field. For now, do nothing.
# 	try:
# 		tRNA_note = feature.qualifiers['note'][0].lower().strip()
# 	except:
# 		tRNA_note = ''
#
#
#







# #TODO implement this function
#
# def check_tRNA_product(product_field):
#
#
#
# 	#This is an initial attempt at checking the tRNA product description
# 	#Ultimately, a regular expression would be better to use
# 	#tRNA product example = 'tRNA-Ser (AGC)'
#
# 	#The code block below functions, but it does not fully account for
# 	#tRNA-OTHER descriptions, tRNA-Stop descriptions,
# 	#and it does not check the accuracy of
# 	#the amino acid and anticodon pairing.
# 	#The biggest problem is that the expected product and note descriptions
# 	#are expected to change after they reach NCBI, so it is not clear
# 	#how to best address that issue here, since nothing in the import
# 	#table reflects WHERE the annotated genome came from.
#
#
#
# 	product_error = 0
#
#
# 	#product starts off as lowercase 'trna-ser (agc)'
# 	#split1_list = 'trna' and 'ser (agc)'
# 	tRNA_product_split1_list = product_field.split('-')
#
# 	#If product is missing, an error will already have been thrown.
# 	#The product should have a hypthen, so only parse product if it can be
# 	#split into two elements.
# 	if len(tRNA_product_split1_list) == 2:
#
# 		tRNA_product_split1_prefix = tRNA_product_split1_list[0].strip() #'trna'
#
# 		#split2_list = 'ser' and 'agc)'
# 		tRNA_product_split2_list = tRNA_product_split1_list[1].split('(')
# 		tRNA_product_amino_acid_three = tRNA_product_split2_list[0].strip() #'ser'
#
# 		if tRNA_product_amino_acid_three != 'other' and \
# 			tRNA_product_amino_acid_three != 'stop' and \
# 			len(tRNA_product_amino_acid_three) != 3:
# 				product_error += 1
#
# 		#The code block below checks for the presence of an anticodon.
# 		#No need to use it currently, since there is so much variability
# 		#at the tRNA product field, but retain for future use.
# 		# if len(tRNA_product_split2_list) == 2:
# 		#
# 		#     #split3_list = 'agc' and ''
# 		#     tRNA_product_split3_list = tRNA_product_split2_list[1].split(')')
# 		#
# 		#     #Only check the anticodon if the amino acid is NOT 'other'
# 		#     if tRNA_product_amino_acid_three != 'other' and \
# 		#         len(tRNA_product_split3_list) == 2:
# 		#
# 		#         tRNA_product_anticodon = tRNA_product_split3_list[0].strip() #'agc'
# 		#         if len(tRNA_product_anticodon) != 3:
# 		#             product_error += 1
# 		#     else:
# 		#         product_error += 1
# 		#
# 		# else:
# 		#     product_error += 1
#
#
# 	else:
# 		product_error += 1
#
# 	return product_error







# #TODO the code below is pasted from import script.
# # It evaluates tRNA features, and needs to be implemented in the
# # check_tRNA function
#
# def check_trna_feature(feature):
#
#
# 	#Now that start and stop have been parsed, check if coordinates are fuzzy or not
# 	if (tRNA_left.isdigit() and tRNA_right.isdigit()):
# 		tRNA_left = int(tRNA_left)
# 		tRNA_right = int(tRNA_right)
# 	else:
# 		write_out(output_file,"\nError: tRNA starting at %s has fuzzy coordinates in phage %s."\
# 				% (tRNA_left,phageName))
# 		record_errors += 1
# 		continue
#
#
#
#
#
# 	if len(tRNA_seq) != tRNA_size:
# 		write_out(output_file,"\nError: unable to retrieve sequence for tRNA starting at %s in phage %s."\
# 				% (tRNA_left + 1,phageName))
# 		record_errors += 1
# 		continue
#
#
#
# 	#Check to see if forward strand terminal nucleotide is correct = A or C
# 	if tRNA_seq[-1] != 'A' and tRNA_seq[-1] != 'C':
# 		record_warnings += 1
# 		write_out(output_file,"\nWarning: tRNA starting at %s does not appear to have correct terminal nucleotide in %s phage." \
# 				% (tRNA_left + 1,phageName))
# 		record_errors += question("\nError: tRNA starting at %s has incorrect terminal nucleotide in %s phage." \
# 				% (tRNA_left + 1,phageName))
#
# 	if tRNA_size < 60 or tRNA_size > 100:
# 		record_warnings += 1
# 		write_out(output_file,"\nWarning: tRNA starting at %s does not appear to be the correct size in %s phage."  \
# 				% (tRNA_left + 1,phageName))
# 		record_errors += question("\nError: tRNA starting at %s is incorrect size in %s phage." \
# 				% (tRNA_left + 1,phageName))
#
#
#
#
#
# 		if len(tRNA_product) > 0:
# 			if check_tRNA_product(tRNA_product) > 0:
# 				write_out(output_file,"\nError: tRNA starting at %s has incorrect amino acid or anticodon in %s." \
# 					% (tRNA_left + 1, phageName))
# 				record_errors += 1
# 		else:
# 			write_out(output_file,"\nError: tRNA starting at %s has incorrect product in %s." \
# 				% (tRNA_left + 1, phageName))
# 			record_errors += 1

#TODO revamp this code into a function






#TODO work on this function
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















#TODO import appropriate classes
#Input = GenBank record object parsed from Biopython
def parse_flat_file_data(filepath):



    genome_object = Genome()

	retrieved_record in SeqIO.read(filepath, "genbank")

    #Keep track of the file from which the record is derived.
    genome_object.set_filename(filepath)

    try:
        genome_object.record_name = retrieved_record.name
    except:
        genome_object.record_name = ""


    try:
        record_organism = retrieved_record.annotations['organism']
        genome_object.record_organism = record_organism

        #Truncate organism name for the 'phage name' and 'host name' field

		value1, value2 = parse_organism_field(record_organism)
		genome_object.set_phage_name(value1)
		genome_object.set_phage_id()
		genome_object.host_name = value2


    except:
        genome_object.record_organism = ""


    try:
        genome_object.record_id = retrieved_record.id
    except:
        genome_object.record_id = ""

    try:
        # There may be a list of accessions associated with this file.
		# I think the first accession in the list is the most recent.
        record_accession = retrieved_record.annotations['accessions'][0]
    except:
        record_accession = ""

	genome_object.set_accession_empty(record_accession)


    try:
        genome_object.set_record_description(retrieved_record.description)
    except:
        genome_object.set_record_description("")

    try:
        genome_object.set_record_source(retrieved_record.annotations['source'])
    except:
        genome_object.set_record_source("")



    try:
        #The retrieved authors can be stored in multiple Reference elements
        record_references = retrieved_record.annotations['references']
        record_references_author_list = []
        for reference in record_references:
            record_references_author_list.append(reference.authors)
        record_author_string = ';'.join(record_references_author_list)
        genome_object.record_authors = record_author_string
    except:
        genome_object.record_authors = ""



    #Nucleotide sequence and errors
    genome_object.set_sequence(retrieved_record.seq)


    #Iterate through all features
    source_feature_list = []
    cds_feature_list = []
	trna_feature_list = []


    source_object_list = []
    cds_object_list = []
	trna_object_list = []


    #A good bit of the code for parsing features is copied from import_phage.py
    for feature in retrieved_record.features:

        #Retrieve the Source Feature info
        if feature.type == "source":
			source_feature_list.append(feature)

		elif feature.type == "tRNA":
			trna_feature_list.append(feature)

		elif feature.type == "CDS":
			cds_feature_list.append(feature)


        else:
			pass



    #Parse the features and assign to the genome object
    genome_object.source_features = source_feature_list

    for feature in cds_feature_list:
        cds_object = parse_cds_feature(feature)
        cds_object_list.append(cds_object)
    genome_object.cds_features = cds_object_list


    for feature in trna_feature_list:
        #TODO combine all tRNA analysis into one function
        trna_object = parse_trna_feature(feature)
        trna_object_list.append(trna_object)
    genome_object.trna_features = trna_object_list




    # If there is only one source feature present, parse it.
    if len(source_feature_list) == 1:

        source_feature = source_feature_list[0]

        try:
            genome_object.source_feature_organism = \
            str(source_feature.qualifiers['organism'][0])

        except:
            pass

        try:
            genome_object.source_feature_host = \
            str(source_feature.qualifiers['host'][0])

        except:
            pass

        try:
            genome_object.source_feature_lab_host = \
            str(source_feature.qualifiers['lab_host'][0])

        except:
            pass

    # TODO return object
    return pass


#Identify list of files that are Genbank-formatted.
def create_flat_file_list(all_files):


    #TODO refactor?
    #write_out(output_file,"\n\n\n\nAccessing genbank-formatted files for add/replace actions...")


    list_of_failed_genome_files = []
    list_of_valid_genome_files = []

    for filename in all_files:

        valid = validate_flat_file(filename)

        if valid == 1:
            list_of_valid_genome_files.append(filename)
        else:
            list_of_failed_genome_files.append(filename)

            # TODO error handling - file might not have contained a parsable
        # flat file. Test whether the genome object contains data.
        # If it is empty, add filename to list of failed files.

    return list_of_valid_genome_files,list_of_failed_genome_files


# Function iterates through list of files and returns
# a list of GenBank-formatted flat files and a list of file names
# they could not be parsed.


def parse_all_flat_files():
    list_of_flat_file_genomes = []

    for filename in list_of_valid_genome_files:

        filepath = os.path.join(path_to_folder,filename)
        flat_file_genome = parse_flat_file_data(filepath)
        list_of_flat_file_genomes.append(flat_file_genome)

    return list_of_flat_file_genomes
