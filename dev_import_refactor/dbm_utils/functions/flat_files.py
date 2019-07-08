"""Functions to interact with, use, and parse genomic data from
GenBank-formatted flat files."""





from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from Bio import Alphabet
from Bio.Seq import Seq
from classes import Genome, Eval, Cds, Trna, Source
from functions import basic
from constants import constants
from datetime import datetime






# TODO this function may need to be improved.
# The coordinates are
# returned as integers, but it is possible that coordinates from Biopython
# SeqFeature object are not integers, so this needs to be accounted for.
# Also, this function doesn't test to confirm that coordinates are
# ExactPosition objects (i.e. non-fuzzy coordinates).
def parse_coordinates(feature):
    """Parse the boundary coordinates for a biopython feature object
    derived from a GenBank-formatted flat file.
    Parsing these coordinates can be tricky.
    There can be more than one set of coordinates if it is
    a compound location. Also, the boundaries may not be precise;
    instead they may be open or fuzzy.
    """


    if (isinstance(feature.location, FeatureLocation) or \
        isinstance(feature.location, CompoundLocation)):

        if feature.strand is None:
            left_boundary = -1
            right_boundary = -1
            parts = 0

        elif isinstance(feature.location, FeatureLocation):
            left_boundary = int(feature.location.start)
            right_boundary = int(feature.location.end)
            parts = 1

        elif isinstance(feature.location, CompoundLocation):

            parts = len(feature.location.parts)

            # Skip this compound feature if it is comprised of more
            # than two features (too tricky to parse).
            if parts == 2:

                # Retrieve compound feature positions based on strand.
                if feature.strand == 1:
                    left_boundary = int(feature.location.parts[0].start)
                    right_boundary = int(feature.location.parts[1].end)

                elif feature.strand == -1:
                    left_boundary = int(feature.location.parts[1].start)
                    right_boundary = int(feature.location.parts[0].end)

                else:
                    pass
            else:
                left_boundary = -1
                right_boundary = -1

    else:
        left_boundary = -1
        right_boundary = -1
        parts = 0

    return (left_boundary, right_boundary, parts)


def parse_cds_feature(cds, feature, parent_translation_table = 11):
    """Parses a Biopython CDS SeqFeature.
    """

    cds.seqfeature = feature

    try:
        cds.locus_tag = feature.qualifiers["locus_tag"][0]
    except:
        cds.locus_tag = ""


    # Orientation
    cds.set_strand(feature.strand, "fr_short", case = True)

    cds.left_boundary, \
    cds.right_boundary, \
    cds.compound_parts = parse_coordinates(feature)


    # Coordinate format for GenBank flat file features parsed by Biopython
    # are 0-based half open intervals.
    cds.coordinate_format = "0_half_open"

    try:
        translation = feature.qualifiers["translation"][0]
    except:
        translation = ""

    translation = Seq(translation, Alphabet.IUPAC.protein)
    cds.set_translation(translation)

    cds.set_nucleotide_length()


    try:
        cds.translation_table = feature.qualifiers["transl_table"][0]
    except:
        cds.translation_table = ""

    try:
        cds.product_description, \
        cds.processed_product_description = \
            basic.reformat_description(feature.qualifiers["product"][0])

    except:
        cds.product_description = ""
        cds.processed_product_description = ""

    try:
        cds.function_description, \
        cds.processed_function_description = \
            basic.reformat_description(feature.qualifiers["function"][0])

    except:
        cds.function_description = ""
        cds.processed_function_description = ""

    try:
        cds.note_description, \
        cds.processed_note_description = \
            basic.reformat_description(feature.qualifiers["note"][0])

    except:
        cds.note_description = ""
        cds.processed_note_description = ""

    try:
        cds.gene_number = feature.qualifiers["gene"][0]
    except:
        cds.gene_number = ""



    cds.parent_translation_table = parent_translation_table





def create_cds_objects(biopython_feature_list):
    """Convert all Biopython CDS SeqFeatures to CdsFeature objects."""
    cds_object_list = []

    for feature in biopython_feature_list:
        cds = Cds.CdsFeature()
        parse_cds_feature(cds, feature)
        cds_object_list.append(cds)

    return cds_object_list




def parse_source_feature(source, feature):
    """Parses a Biopython Source Feature.
    """

    try:
        source.organism = str(feature.qualifiers["organism"][0])
    except:
        source.organism = ""
    try:
        source.host = str(feature.qualifiers["host"][0])
    except:
        source.host = ""
    try:
        source.lab_host = str(feature.qualifiers["lab_host"][0])
    except:
        source.lab_host = ""




def create_source_objects(biopython_feature_list):
    """Convert all Biopython Source SeqFeatures to SourceFeature objects."""
    source_object_list = []

    for feature in biopython_feature_list:
        source = Source.SourceFeature()
        parse_source_feature(source, feature)
        source_object_list.append(source)

    return source_object_list




def create_feature_dictionary(feature_list):
    """From a list of all Biopython SeqFeatures derived from a GenBank-formatted
    flat file, create a dictionary of SeqFeatures based on their type.
    Key = feature type (source, tRNA, CDS, other).
    Value = list of features."""

    feature_type_set = set()
    feature_dict = {}

    for feature in feature_list:
        feature_type_set.add(feature.type)

    for type in feature_type_set:
        feature_sublist = []

        index = 0
        while index < len(feature_list):

            feature = feature_list[index]

            if feature.type == type:
                feature_sublist.append(feature)
            index += 1

        feature_dict[type] = feature_sublist

    return feature_dict



def parse_flat_file_data(genome_obj, \
                            retrieved_record, \
                            filepath = "", \
                            translation_table = 11, \
                            phage_id_field = "record_organism_phage_name"):

    """Parses a GenBank-formatted flat file into a Genome object using
    data that has already been parsed by Bio.SeqIO.
    """

    # Keep track of the file from which the record is derived.
    genome_obj.set_filename(filepath)


    try:
        genome_obj.record_name = retrieved_record.name

        # It appears that if name is not present, Biopython auto-populates
        # this attribute as "<unknown name>"
        if genome_obj.record_name == "<unknown name>":
            genome_obj.record_name = ""

    except:
        genome_obj.record_name = ""



    try:
        genome_obj.record_organism = retrieved_record.annotations['organism']
    except:
        genome_obj.record_organism = ""

    # Identifies host and phage name from organism field.
    genome_obj.parse_record_organism()


    try:
        genome_obj.record_id = retrieved_record.id

        # It appears that if id is not present, Biopython auto-populates
        # this attribute as "<unknown id>"
        if genome_obj.record_id == "<unknown id>":
            genome_obj.record_id = ""

    except:
        genome_obj.record_id = ""

    try:
        # Since accessions are stored in a list, there may be more than
        # one accessions associated with this file.
        # The first accession in the list is assumed to be the most recent.
        accession = retrieved_record.annotations['accessions'][0]
    except:
        accession = ""

    genome_obj.set_accession(accession)


    try:
        genome_obj.record_description = retrieved_record.description

        # It appears that if description is not present, Biopython
        # auto-populates this attribute as "<unknown description>"
        if genome_obj.record_description == "<unknown description>":
            genome_obj.record_description = ""

    except:
        genome_obj.record_description = ""

    # Identifies host and phage name from description field.
    genome_obj.parse_record_description()

    try:
        genome_obj.record_source = retrieved_record.annotations['source']
    except:
        genome_obj.record_source = ""

    # Identifies host and phage name from record source field.
    genome_obj.parse_record_source()



    try:
        # The retrieved authors can be stored in multiple Reference elements.
        refs = retrieved_record.annotations['references']
        authors_list = []
        for ref in refs:

            # Note: Reference objects are instantiated with an empty
            # authors attribute. So if no authors are present in a Reference,
            # it will still concatenate an empty string, resulting in an
            # author_string = ";;;" etc. So only add the authors info if
            # it is not an empty string.
            if ref.authors != "":
                authors_list.append(ref.authors)

        authors_string = ';'.join(authors_list)
        genome_obj.record_authors = authors_string

    except:
        genome_obj.record_authors = ""



    # Biopython requires the parsed record contains a sequence, so
    # no need to test whether the seq attribute is present or not.
    # Nucleotide sequence, length, and % GC.
    genome_obj.set_sequence(retrieved_record.seq)



    try:
        date = retrieved_record.annotations["date"]
        genome_obj.record_date = datetime.strptime(date,'%d-%b-%Y')
    except:
        genome_obj.record_date = basic.convert_empty("", "empty_datetime_obj")
        # TODO throw an error if no date?


    # Now that record fields are parsed, set the phage_id.
    genome_obj.set_phage_id_from_field(phage_id_field)

    genome_obj.type = "flat_file"


    # Create lists of parsed features.
    # Note: Biopython instantiates the features attribute with
    # an empty list, so no need to test if features attribute is
    # present or not.
    feature_dict = create_feature_dictionary(retrieved_record.features)

    if "CDS" in feature_dict.keys():
        cds_object_list = create_cds_objects(feature_dict["CDS"])
    else:
        cds_object_list = []


    # TODO the parent_phage_id can't be set until the genome phage_id
    # is set. But the phage_id is not determined from within the
    # flat file. It is inputted externally, such as from a ticket.
    # Once the ticket is used to set the phage_id, the parent_phage_id
    # atttributes can be set for all features.
    # for cds in cds_object_list:
    #     cds.parent_phage_id = genome_obj.phage_id


    if "source" in feature_dict.keys():
        source_object_list = create_source_objects(feature_dict["source"])
    else:
        source_object_list = []

    if "tRNA" in feature_dict.keys():
        trna_object_list = create_trna_objects(feature_dict["tRNA"])
    else:
        trna_object_list = []

    if "tmrna" in feature_dict.keys():
        tmrna_object_list = create_tmrna_objects(feature_dict["tmrna"])
    else:
        tmrna_object_list = []

    genome_obj.translation_table = translation_table
    genome_obj.set_cds_features(cds_object_list)
    genome_obj.set_source_features(source_object_list)
    genome_obj.set_trna_features(trna_object_list)
    # genome_obj.set_tmrna_features(tmrna_object_list)



# TODO refactor to simply return true?
def check_flat_file_type(filepath):
    """Verify that the file contains a file extension common to
    GenBank-formatted flat files."""
    filename = filepath.split("/")[-1]

    if filename.split('.')[-1] not in constants.ADMISSIBLE_FILE_TYPES:
        result = "File does not have a valid file extension, " + \
            "so it will not be processed."
        status = "error"
    else:
        result = "Flat file extension is valid."
        status = "correct"

    # TODO add eval id?
    definition = "Check flat file type."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)

    return eval




def parse_flat_file(filepath):
    """Determine whether the file contains a single GenBank-formatted record
    that can be parsed by Biopython SeqIO.
    Files may contain 0, 1, or >1 parseable records.
    When SeqIO parses files, if there are 0 Genbank-formatted records,
    it does not throw an error, but simply moves on.
    Only when there is one record in the file is the parsed data returned.
    It is not clear how often there multiple records are present in
    non-SEA-PHAGES GenBank records, but it is a good idea to verify
    there is only one record per file before proceeding."""


    records = []
    record = None

    # If Biopython is unable to parse the file, an error is encountered.
    try:
        records = list(SeqIO.parse(filepath, "genbank"))
    except:
        records = None

    if records is None:

        result = "Biopython is unable to parse this file, " + \
            "so it will not be processed."
        status = "error"

    elif len(records) == 0:

        result = "Biopython did not find any GenBank-formatted records, " + \
            "so it will not be processed."
        status = "error"

    elif len(records) > 1:

        result = "Biopython found multiple records, " + \
            "so it will not be processed."
        status = "error"

    else:
        result = "Record was successfully parsed."
        status = "correct"
        record = records[0]


    # TODO add an eval id?
    definition = "Parse the flat file record."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)

    return (record, eval)


def create_parsed_flat_file_list(all_files,
                                phage_id_field = "record_organism_phage_name"):
    """Create a list of genome objects containing data parsed from
    flat files."""

    failed_files = []
    valid_files = []
    genomes = []
    all_results = []

    for filename in all_files:

        file_results = []
        genome_obj = Genome.Genome()
        result1 = check_flat_file_type(filename)
        record, result2 = parse_flat_file(filename)

        if result1.status != "correct":
            file_results.append(result1)

        if result2.status != "correct":
            file_results.append(result2)

        if len(file_results) == 0:
            parse_flat_file_data(genome_obj, record, filename, phage_id_field)
            valid_files.append(filename)
            genomes.append(genome_obj)

        else:
            failed_files.append(filename)
            all_results += file_results

    return (genomes, all_results, valid_files, failed_files)
















# TODO unit test below.











# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_trna_feature(trna, feature):
    """Parses a Biopython tRNA Feature.
    """

    pass


# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def create_trna_objects(biopython_feature_list):
    """Convert all Biopython tRNA SeqFeatures to tRNAFeature objects."""
    trna_object_list = []

    for feature in biopython_feature_list:
        trna = ""
        trna_object_list.append(trna)

    return trna_object_list



# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_tmrna_feature(tmrna, feature):
    """Parses a Biopython tRNA Feature.
    """

    pass


# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def create_tmrna_objects(biopython_feature_list):
    """Convert all Biopython tRNA SeqFeatures to tRNAFeature objects."""
    tmrna_object_list = []

    for feature in biopython_feature_list:
        tmrna = ""
        tmrna_object_list.append(tmrna)

    return tmrna_object_list



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
