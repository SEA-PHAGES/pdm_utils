"""Functions to interact with, use, and parse genomic data from
GenBank-formatted flat files."""





from Bio import SeqIO
from classes import Genome
from classes import Eval
from classes import Cds
from classes import Trna
from functions import basic





















# TODO unit test below.




# TODO create function to retrieve_flatfile_record
def retrieve_flatfile_record(filepath)


    #

    # return retrieve_record parse from SeqIO
    pass



# TODO unit test.
def retrieve_flat_file_features(retrieved_record):
    """Groups different types of features from a GenBank-formatted
    flat file into separate lists.
    This function returns a dictionary.
    Key = feature type (source, tRNA, CDS, other).
    Value = list of features.
    """

    # List of all Biopython feature objects parsed from record.
    biopython_feature_dict = {}

    sources = []
    cdss = []
    trnas = []
    others = []

    for feature in retrieved_record.features:
        if feature.type == "source":
            sources.append(feature)
        elif feature.type == "tRNA":
            trnas.append(feature)
        elif feature.type == "CDS":
            cdss.append(feature)
        else:
            others.append(feature)

    biopython_feature_dict["source"] = sources
    biopython_feature_dict["tRNA"] = trnas
    biopython_feature_dict["CDS"] = cdss
    biopython_feature_dict["other"] = others

    return biopython_feature_dict



# TODO unit test.
def parse_coordinates(feature):
    """Parse the boundary coordinates for a feature derived from a
    GenBank-formatted flat file. Parsing these coordinates can be tricky.
    There can be more than one set of coordinates if it is
    a compound location. Also, the boundaries may not be precise;
    instead they may be open or fuzzy.
    """

    if str(feature.location)[:4] == 'join':

        compound_parts = len(feature.location.parts)

        # Skip this compound feature if it is comprised of more
        # than two features (too tricky to parse).
        if compound_parts <= 2:

            # Retrieve compound feature positions based on strand.
            if feature.strand == 1:
                left_boundary = str(feature.location.parts[0].start)
                right_boundary = str(feature.location.parts[1].end)

            elif feature.strand == -1:
                left_boundary = str(feature.location.parts[1].start)
                right_boundary = str(feature.location.parts[0].end)

            else:
                left_boundary = 0
                right_boundary = 0
                message = \
                    "Unable to parse coordinates since strand is undefined."
                eval_result = Eval.construct_error(message)
        else:

            left_boundary = 0
            right_boundary = 0
            message = \
                "This is a compound feature that is unable to be parsed."
            eval_result = Eval.construct_warning(message, message)

    else:
        left_boundary = str(feature.location.start)
        right_boundary = str(feature.location.end)
        compound_parts = 1

    return (left_boundary, right_boundary, compound_parts)





# TODO unit test.
def parse_cds_feature(feature):
    """Parses a Biopython CDS Feature.
    """


    cds = Cds.CdsFeature()

    cds.type_id = "CDS"

    try:
        cds.locus_tag = feature.qualifiers["locus_tag"][0]
    except:
        cds.locus_tag = ""


    # Orientation
    # TODO confirm strand formatting is correct.
    cds.set_strand(feature.strand, "fr_long")


    cds.left_boundary, cds.right_boundary, cds.compound_parts = \
        parse_coordinates(feature)

    cds.coordinate_format = "" # TODO create name for this value.

    try:
        cds.set_translation(feature.qualifiers["translation"][0])
    except:
        cds.set_translation("")

    cds.set_nucleotide_length()


    try:
        cds.translation_table = feature.qualifiers["transl_table"][0]
    except:
        pass

    try:
        product, processed_product = \
            basic.reformat_description(feature.qualifiers["product"][0])
        cds.product_description = product
        cds.processed_product_description = processed_product

    except:
        pass

    try:
        function, processed_function = \
            basic.reformat_description(feature.qualifiers["function"][0])
        cds.function_description = function
        cds.processed_function_description = processed_function

    except:
        pass

    try:
        note, processed_note = \
            basic.reformat_description(feature.qualifiers["note"][0])
        cds.note_description = note
        cds.processed_note_description = processed_note

    except:
        pass

    try:
        cds.gene_number = feature.qualifiers["gene"][0]
    except:
        pass

    return cds







# TODO unit test.
def create_feature_dictionary(feature_list):
    """From a list of all features derived from a GenBank-formatted
    flat file, create a dictionary of features based on their type."""

    feature_type_set = set()
    feature_dict = {}

    for feature in feature_list:
        feature_type_set.add(feature.type)

    for type in feature_type_set:
        subfeature_list = []

        index = 0
        while index < len(feature_list):

            if feature.type == type:
                subfeature_list.append(feature)
            index += 1

        feature_dict[type] = subfeature_list

    return feature_dict




# TODO implement.
# TODO unit test.
def create_cds_objects(feature_list):
    """."""
    object_list = []
    for feature in feature_list:
        cds_object = parse_cds_feature(feature)
        object_list.append(cds_object)
    return object_list



# TODO implement.
# TODO unit test.
def create_trna_objects(feature_list):
    """."""
    object_list = []
    for feature in feature_list:
        cds_object = parse_trna_feature(feature)
        object_list.append(cds_object)
    return object_list

# TODO implement.
# TODO unit test.
def create_source_objects(feature_list):
    """."""
    object_list = []
    for feature in feature_list:
        cds_object = parse_source_feature(feature)
        object_list.append(cds_object)
    return object_list

















#TODO import appropriate classes
#Input = GenBank record object parsed from Biopython
def parse_flat_file_data(filepath):
    """Parses a GenBank-formatted flat file into a Genome object using
    data that has already been parsed by Bio.SeqIO.
    """

    genome_obj = Genome.Genome()

    # TODO this line probably doesn't work right.
    retrieved_record = SeqIO.read(filepath, "genbank")

    #Keep track of the file from which the record is derived.
    genome_obj.set_filename(filepath)

    try:
        genome_obj.record_name = retrieved_record.name
    except:
        genome_obj.record_name = ""




    # TODO I need to set these attributes:
    #Truncate organism name for the 'phage name' and 'host name' field
    # value1, value2 = parse_organism_field(record_organism)
    # genome_obj.set_phage_name(value1)
    # genome_obj.set_phage_id()
    # genome_obj.host_name = value2

    try:
        genome_obj.record_organism = retrieved_record.annotations['organism']


    except:
        genome_obj.record_organism = ""

    # Determines host and phage name from organism field.
    genome_obj.parse_record_organism()


    try:
        genome_obj.record_id = retrieved_record.id
    except:
        genome_obj.record_id = ""

    try:
        # Accessions are stored in a list.
        # There may be more than one accessions associated with this file.
        # The first accession in the list is assumed to be the most recent.
        accession = retrieved_record.annotations['accessions'][0]
    except:
        accession = ""

    genome_obj.set_accession(accession)


    try:
        genome_obj.record_description = retrieved_record.description
    except:
        genome_obj.record_description = ""

    # Determines host and phage name from description field.
    genome_obj.parse_record_description()

    try:
        genome_obj.record_source = retrieved_record.annotations['source']
    except:
        genome_obj.record_source = ""

    # Determines host and phage name from source field.
    genome_obj.parse_record_source()


    try:

        # TODO the authors could be returned as a list, instead of joining
        # into one long string.
        # The retrieved authors can be stored in multiple Reference elements.
        record_references = retrieved_record.annotations['references']
        record_references_author_list = []
        for reference in record_references:
            record_references_author_list.append(reference.authors)
        record_author_string = ';'.join(record_references_author_list)
        genome_obj.record_authors = record_author_string
    except:
        genome_obj.record_authors = ""



    # Nucleotide sequence, length, and % GC.
    genome_obj.set_sequence(retrieved_record.seq)





    ### TODO BELOW IN PROGRESS OF RE-FACTORING, PAUSED
    # TO DEVELOP FEATURE LIST RETRIEVAL


    # TODO create a dictionary of all features.

    # TODO create a list of CDS objects.

    # TODO create a list of tRNA objects.

    # TODO create a list of source objects.

    #Parse the features and assign to the genome object
    genome_obj.source_features = source_feature_list

    for feature in cds_feature_list:
    genome_obj.cds_features = cds_object_list


    for feature in trna_feature_list:
    genome_obj.trna_features = trna_object_list






    # If there is only one source feature present, parse it.
    if len(source_feature_list) == 1:

        source_feature = source_feature_list[0]

        try:
            genome_obj.source_feature_organism = \
            str(source_feature.qualifiers['organism'][0])

        except:
            pass

        try:
            genome_obj.source_feature_host = \
            str(source_feature.qualifiers['host'][0])

        except:
            pass

        try:
            genome_obj.source_feature_lab_host = \
            str(source_feature.qualifiers['lab_host'][0])

        except:
            pass

    # TODO return object

    return pass

















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
















#
# #TODO implement this function
#
# #Parse tRNA features
# def parse_trna_feature(feature):
#
#     #TODO need to implement this function
#     #return(None)
#
#     #Retrieve tRNA coordinates
#     try:
#         #Biopython converts coordinates to 0-index
#         #Start(left) coordinates are 0-based inclusive (feature starts there)
#         #Stop (right) coordinates are 0-based exclusive (feature stops 1bp prior to coordinate)
#         tRNA_left = str(feature.location.start)
#         tRNA_right = str(feature.location.end)
#
#     except:
#         #TODO error handling
#         write_out(output_file,"\nError: a tRNA has incorrect coordinates in phage %s."\
#                 % phageName)
#         record_errors += 1
#         continue
#
#     #Retrieve top strand of tRNA feature. It is NOT necessarily
#     #in the correct orientation
#     tRNA_size = abs(tRNA_right - tRNA_left)
#     tRNA_seq = phageSeq[tRNA_left:tRNA_right].upper()
#
#     #Convert sequence to reverse complement if it is on bottom strand
#     if feature.strand == 1:
#         pass
#     elif feature.strand == -1:
#         tRNA_seq = tRNA_seq.reverse_complement()
#     else:
#         #TODO error handling
#         record_errors += 1
#         write_out(output_file,"\Error: tRNA starting at %s does not have proper orientation in %s phage." \
#                 % (tRNA_left + 1,phageName))
#         continue
#
#     #Retrieve and check product
#     try:
#         tRNA_product = feature.qualifiers['product'][0].lower().strip()
#     except:
#         #TODO error handling
#         write_out(output_file,"\nError: tRNA starting at %s is missing product field in phage %s." \
#             % (tRNA_left + 1,phageName))
#         record_errors += 1
#         tRNA_product = ''
#
#     #Retrieve note
#     #In the future, this field may need to be parsed in a similar
#     #manner as the product field. For now, do nothing.
#     try:
#         tRNA_note = feature.qualifiers['note'][0].lower().strip()
#     except:
#         tRNA_note = ''
#
#     return pass





###
