"""Functions to interact with, use, and parse genomic data from
GenBank-formatted flat files."""


from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from Bio import Alphabet
from Bio.Seq import Seq
from pdm_utils.classes import genome, cds, trna, source
from pdm_utils.functions import basic
from pdm_utils.constants import constants
from datetime import datetime
from pdm_utils.classes import genomepair





# TODO this function may need to be improved.
# The coordinates are
# returned as integers, but it is possible that coordinates from Biopython
# SeqFeature object are not integers, so this needs to be accounted for.
# Also, this function doesn't test to confirm that coordinates are
# ExactPosition objects (i.e. non-fuzzy coordinates).
def parse_coordinates(seqfeature):
    """Parse the boundary coordinates from a GenBank-formatted flat file.

    The functions takes a Biopython SeqFeature object containing data
    that was parsed from the feature in the flat file. Parsing these
    coordinates can be tricky.
    There can be more than one set of coordinates if it is
    a compound location. Only features with 1 or 2 open reading frames
    (parts) are correctly parsed. Also, the boundaries may not be precise;
    instead they may be open or fuzzy. Non-precise coordinates are
    converted to '-1'. If the strand is undefined, the coordinates
    are converted to '-1' and parts is set to '0'. If an incorrect
    data type is provided, coorindates are set to '-1' and parts is
    set to '0'.

    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :returns:
        tuple (left, right, parts)
        WHERE
        left(int) is the first coordinate, regardless of strand.
        right(int) is the second coordinate, regardless of strand.
        parts(int) is the number of open reading frames that define
        the feature.
    """

    if (isinstance(seqfeature.location, FeatureLocation) or \
        isinstance(seqfeature.location, CompoundLocation)):

        if seqfeature.strand is None:
            left = -1
            right = -1
            parts = 0
        elif isinstance(seqfeature.location, FeatureLocation):
            left = int(seqfeature.location.start)
            right = int(seqfeature.location.end)
            parts = 1
        elif isinstance(seqfeature.location, CompoundLocation):
            parts = len(seqfeature.location.parts)

            # Skip this compound seqfeature if it is comprised of more
            # than two features (tricky to parse).
            if parts == 2:

                # Retrieve compound seqfeature positions based on strand.
                if seqfeature.strand == 1:
                    left = int(seqfeature.location.parts[0].start)
                    right = int(seqfeature.location.parts[1].end)
                elif seqfeature.strand == -1:
                    left = int(seqfeature.location.parts[1].start)
                    right = int(seqfeature.location.parts[0].end)
                else:
                    pass
            else:
                left = -1
                right = -1
        else:
            pass
    else:
        left = -1
        right = -1
        parts = 0
    return (left, right, parts)


# TODO implement set_wrap()?
def parse_cds_seqfeature(seqfeature, genome_id=""):
    """Parse data from a Biopython CDS SeqFeature object into a Cds object.

    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :param genome_id: An identifier for the genome in which the seqfeature
    is defined.
    :type genome_id: str
    :returns: A  pdm_utils Cds object
    :rtype: Cds
    """
    cds_ftr = cds.Cds()
    cds_ftr.genome_id = genome_id
    cds_ftr.seqfeature = seqfeature

    try:
        locus_tag = seqfeature.qualifiers["locus_tag"][0]
    except:
        locus_tag = ""
    cds_ftr.set_locus_tag(locus_tag)

    cds_ftr.set_strand(seqfeature.strand, "fr_short", case = True)
    cds_ftr.left, cds_ftr.right, cds_ftr.compound_parts = parse_coordinates(seqfeature)
    #cds_ftr.set_wrap()

    # Coordinate format for GenBank flat file features parsed by Biopython
    # are 0-based half open intervals.
    cds_ftr.coordinate_format = "0_half_open"


    # For translation, convert it to a Biopython Seq object.
    try:
        translation = seqfeature.qualifiers["translation"][0]
    except:
        translation = ""
    translation = Seq(translation, Alphabet.IUPAC.protein)
    cds_ftr.set_translation(translation)
    cds_ftr.set_nucleotide_length()

    try:
        translation_table = seqfeature.qualifiers["transl_table"][0]
    except:
        translation_table = -1
    cds_ftr.set_translation_table(translation_table)

    try:
        cds_ftr.product, cds_ftr.processed_product = \
            basic.reformat_description(seqfeature.qualifiers["product"][0])
    except:
        cds_ftr.product = ""
        cds_ftr.processed_product = ""

    try:
        cds_ftr.function, cds_ftr.processed_function = \
            basic.reformat_description(seqfeature.qualifiers["function"][0])
    except:
        cds_ftr.function = ""
        cds_ftr.processed_function = ""

    try:
        cds_ftr.note, cds_ftr.processed_note = \
            basic.reformat_description(seqfeature.qualifiers["note"][0])
    except:
        cds_ftr.note = ""
        cds_ftr.processed_note = ""

    try:
        cds_ftr.gene = seqfeature.qualifiers["gene"][0]
    except:
        cds_ftr.gene = ""

    cds_ftr.set_name()
    return cds_ftr


# TODO remove this function. No longer needed.
# def create_cds_objects(seqfeature_list):
#     """Convert all Biopython CDS SeqFeatures to Cds objects."""
#     cds_list = []
#     for seqfeature in seqfeature_list:
#         cds_ftr = parse_cds_seqfeature(seqfeature)
#         cds_list.append(cds_ftr)
#     return cds_list


def parse_source_seqfeature(seqfeature, genome_id=""):
    """Parses a Biopython Source SeqFeature.

    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :param genome_id: An identifier for the genome in which the seqfeature
    is defined.
    :type genome_id: str
    :returns: A pdm_utils Source object
    :rtype: Source
    """
    src_ftr = source.Source()
    src_ftr.genome_id = genome_id
    src_ftr.seqfeature = seqfeature

    try:
        src_ftr.organism = str(seqfeature.qualifiers["organism"][0])
    except:
        src_ftr.organism = ""
    try:
        src_ftr.host = str(seqfeature.qualifiers["host"][0])
    except:
        src_ftr.host = ""
    try:
        src_ftr.lab_host = str(seqfeature.qualifiers["lab_host"][0])
    except:
        src_ftr.lab_host = ""
    return src_ftr


# TODO remove this function.
# def create_source_objects(seqfeature_list):
#     """Convert all Biopython Source SeqFeatures to Source objects."""
#     source_list = []
#     for seqfeature in seqfeature_list:
#         src_ftr = parse_source_seqfeature(seqfeature)
#         source_list.append(src_ftr)
#     return source_list


def create_seqfeature_dictionary(seqfeature_list):
    """Create a dictionary of Biopython SeqFeature objects based on their type.

    From a list of all Biopython SeqFeatures derived from a GenBank-formatted
    flat file, create a dictionary of SeqFeatures based on their 'type'
    attribute.

    :param seqfeature_list: List of Biopython SeqFeatures
    :type seqfeature_list: list
    :param genome_id: An identifier for the genome in which the seqfeature
    is defined.
    :type genome_id: str
    :returns:
        A dictionary of Biopython SeqFeatures:
        Key: SeqFeature type (source, tRNA, CDS, other)
        Value: the SeqFeature
    :rtype: dict
    """

    seqfeature_type_set = set()
    seqfeature_dict = {}
    for seqfeature in seqfeature_list:
        seqfeature_type_set.add(seqfeature.type)
    for type in seqfeature_type_set:
        sublist = []
        index = 0
        while index < len(seqfeature_list):
            seqfeature = seqfeature_list[index]
            if seqfeature.type == type:
                sublist.append(seqfeature)
            index += 1
        seqfeature_dict[type] = sublist
    return seqfeature_dict


def parse_genome_data(seqrecord, filepath="",
        translation_table=11, phage_id_field="organism_name"):
    """Parse data from a Biopython SeqRecord object into a Genome object.

    All Source, CDS, tRNA, and tmRNA features are parsed into their
    associates Source, Cds, Trna, and Tmrna objects.

    :param seqrecord: A Biopython SeqRecord object.
    :type seqrecord: SeqRecord
    :param filepath: A filename associated with the returned Genome object.
    :type filepath: str
    :param translation_table:
        The applicable translation table for the genome's CDS features.
    :type translation_table: int
    :param phage_id_field:
        The SeqRecord attribute from which the unique genome
        identifier/name is stored.
    :type phage_id_field: str
    :returns: A pdm_utils Genome object.
    :rtype: Genome
    """

    # Keep track of the file from which the record is derived.
    gnm = genome.Genome()
    gnm.set_filename(filepath)
    gnm.type = "flat_file"

    try:
        gnm.name = seqrecord.name
        # It appears that if name is not present, Biopython auto-populates
        # this attribute as "<unknown name>"
        if gnm.name == "<unknown name>":
            gnm.name = ""
    except:
        gnm.name = ""

    try:
        gnm.organism = seqrecord.annotations['organism']
    except:
        gnm.organism = ""

    # Identifies host and phage name from organism field.
    gnm.parse_organism()

    try:
        gnm.id = seqrecord.id
        # It appears that if id is not present, Biopython auto-populates
        # this attribute as "<unknown id>"
        if gnm.id == "<unknown id>":
            gnm.id = ""
    except:
        gnm.id = ""

    try:
        # Since accessions are stored in a list, there may be more than
        # one accessions associated with this file.
        # The first accession in the list is assumed to be the most recent.
        accession = seqrecord.annotations['accessions'][0]
    except:
        accession = ""

    gnm.set_accession(accession)

    try:
        gnm.description = seqrecord.description
        # It appears that if description is not present, Biopython
        # auto-populates this attribute as "<unknown description>"
        if gnm.description == "<unknown description>":
            gnm.description = ""
    except:
        gnm.description = ""

    # Identifies host and phage name from description field.
    gnm.parse_description()

    try:
        gnm.source = seqrecord.annotations['source']
    except:
        gnm.source = ""

    # Identifies host and phage name from record source field.
    gnm.parse_source()

    try:
        # The retrieved authors can be stored in multiple Reference elements.
        refs = seqrecord.annotations['references']
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
        gnm.authors = authors_string
    except:
        gnm.authors = ""

    # Biopython requires the parsed record contains a sequence, so
    # no need to test whether the seq attribute is present or not.
    # Nucleotide sequence, length, and % GC.
    gnm.set_sequence(seqrecord.seq)

    try:
        date = seqrecord.annotations["date"]
        gnm.date = datetime.strptime(date,'%d-%b-%Y')
    except:
        gnm.date = basic.convert_empty("", "empty_datetime_obj")

    # Now that record fields are parsed, set the phage_id.
    gnm.set_id(attribute=phage_id_field)

    # Create lists of parsed features.
    # Note: Biopython instantiates the features attribute with
    # an empty list, so no need to test if features attribute is
    # present or not.
    seqfeature_dict = create_seqfeature_dictionary(seqrecord.features)

    cds_list = []
    if "CDS" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["CDS"]:
            cds_ftr = parse_cds_seqfeature(seqfeature, genome_id=gnm.id)
            cds_list.append(cds_ftr)

    source_list = []
    if "source" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["source"]:
            src_ftr = parse_source_seqfeature(seqfeature, genome_id=gnm.id)
            source_list.append(src_ftr)

    # TODO unit test after functions are constructed.
    trna_list = []
    if "tRNA" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tRNA"]:
            trna_ftr = parse_trna_seqfeature(seqfeature, genome_id=gnm.id)
            trna_list.append(trna_ftr)

    # TODO unit test after functions are constructed.
    tmrna_list = []
    if "tmrna" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tmrna"]:
            tmrna = parse_tmrna_seqfeature(seqfeature, genome_id=gnm.id)
            tmrna_list.append(tmrna)

    gnm.translation_table = translation_table
    gnm.set_cds_features(cds_list)
    gnm.set_source_features(source_list)


    gnm.set_trna_features(trna_list)
    # gnm.set_tmrna_features(tmrna_list)

    # The Cds.id is constructed from the Genome.id and the Cds order.
    gnm.set_feature_ids(use_type=True, use_cds=True)

    # TODO set tRNA feature ids.
    #gnm.set_feature_ids(use_type=True, use_trna=True)

    return gnm




def copy_data_to(bndl, type, flag="ticket"):
    """Copy data to a genome object derived from a 'flat_file'.

    The Bundle object is expected to contain at least two Genome objects
    in its 'genome_dict' dictionary. The first 'donor' genome is expected to be
    stored in the dictionary with its key equivalent to its 'type'.
    The second 'receiver' genome is expected to be derived
    from a GenBank-formatted flat file and stored in the Bundle
    object's 'genome_dict' dictionary with a 'flat_file' key.
    A GenomePair is created using these two genomes,
    and the data is copied from the donor genome to the
    'flat_file' genome.

    :param bndl:
        A Bundle object containing both Genome objects stored
        in the 'genome_dict' attribute.
    :type bndl: Bundle
    :param type:
        The value of the donor genome's 'type' attribute, which is
        used as the its key in the Bundle's 'genome_dict' dictionary.
    :type type: str
    :param flag:
        The value used to indicate which 'flat_file'
        attributes should be populated from the donor genome.
        Several 'flat_file' genome attributes are set to this value.
        Using a GenomePair method, the 'flat_file' attributes
        with this flag are re-populated from data of the corresponding
        attributed in the donor genome.
    :type flag: str
    """

    if "flat_file" in bndl.genome_dict.keys():
        genome1 = bndl.genome_dict["flat_file"]
        genome1.cluster = flag
        genome1.subcluster = flag
        genome1.name = flag # TODO should this attribute be copied?
        genome1.host_genus = flag
        genome1.accession = flag
        genome1.cluster_subcluster = flag
        genome1.annotation_status = flag
        genome1.annotation_author = flag
        genome1.annotation_qc = flag
        genome1.retrieve_record = flag
        genome1.set_value_flag(flag)

        if type in bndl.genome_dict.keys():

            genome2 = bndl.genome_dict[type]

            # Copy all data that is set to 'ticket' and
            # add to Bundle object.
            genome_pair = genomepair.GenomePair()
            genome_pair.genome1 = genome1
            genome_pair.genome2 = genome2
            genome_pair.copy_data("type", genome2.type, genome1.type, flag)
            bndl.set_genome_pair(genome_pair, genome1.type, genome2.type)

        # Now record an error if there are still fields
        # that need to be copied.
        genome1.set_value_flag(flag)
        genome1.check_value_flag()


def parse_files(file_list, id_field="organism_name"):
    """Parse data from a list of flat files.

    All GenBank-formatted flat files present in the list of files
    are first parsed into Biopython SeqRecord objects, and
    then parsed into pdm_utils Genome objects.

    :param file_list: A list of filenames.
    :type file_list: list
    :param id_field:
        The name of the attribute in the SeqRecord object
        from which the unique genome identifier/name is stored.
    :type id_field: str
    :returns:
        tuple (genomes, valid_files, failed_files)
        WHERE
        genomes(list) is a list of pdm_utils Genome objects parsed
        from the files.
        valid_files(list) is a list of filenames from which a
        Biopython SeqRecord object was successfully parsed.
        failed_files(list) is a list of filenames from which a
        Biopython SeqRecord object was not successfully parsed.
    """
    failed_files = []
    valid_files = []
    genomes = []
    for filename in file_list:
        try:
            seqrecords = list(SeqIO.parse(filename, "genbank"))
        except:
            seqrecords = []

        if len(seqrecords) == 1:
            gnm = parse_genome_data(seqrecords[0], filename, id_field)
            genomes.append(gnm)
            valid_files.append(filename)
        else:
            failed_files.append(filename)
            # # If there is no parseable record, a genome object is still
            # # created and populated with 'type and 'filename'.
            # gnm = genome.Genome()
            # gnm.type = "flat_file"
            # gnm.set_filename(filename)
    return (genomes, valid_files, failed_files)













# TODO unit test below.











# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_trna_seqfeature(seqfeature, genome_id=""):
    """Parses a Biopython tRNA SeqFeature.
    """
    trna_ftr = ""
    return trna_ftr



# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_tmrna_seqfeature(seqfeature, genome_id=""):
    """Parses a Biopython tRNA SeqFeature.
    """
    tmrna = ""
    return tmrna









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























# Functions that are no longer needed.



# # TODO this is probably no longer needed. There is no need to impose
# # restrictions on file extensions for flat files.
# def check_extension(filepath):
#     """Verify the file extension is common for GenBank-formatted flat files."""
#     valid = False
#     filename = filepath.split("/")[-1]
#     if filename.split('.')[-1] in constants.ADMISSIBLE_FILE_TYPES:
#         valid = True
#     return valid

# TODO the follow create_parsed_flat_file() and
# create_parsed_flat_file_list() functions may no longer be needed.
# def create_parsed_flat_file(filename, id_field="organism_name"):
#     """Create a genome object parsed from flat files."""
#
#     valid = check_extension(filename)
#     if valid:
#         try:
#             records = list(SeqIO.parse(filename, "genbank"))
#         except:
#             records = []
#
#         if len(records) == 1:
#             gnm = parse_genome_data(records[0], filename, id_field)
#         else:
#             gnm = genome.Genome()
#     else:
#         gnm = genome.Genome()
#
#     # TODO currently the parse_genome_data() sets filename and type,
#     # so this is redundant. But some attribute needs to be set
#     # even if there is a problem with parsing the record in the file.
#
#     # If there is no parseable record, a genome object is still
#     # created and populated with 'type and 'filename'.
#     gnm.type = "flat_file"
#     gnm.set_filename(filename)
#
#     return gnm
#
#
# def create_parsed_flat_file_list(all_files, id_field="organism_name"):
#     """Create a list of genome objects containing data parsed from
#     flat files."""
#
#     failed_files = []
#     valid_files = []
#     genomes = []
#     for filename in all_files:
#         gnm = create_parsed_flat_file(filename, id_field = id_field)
#         genomes.append(gnm)
#         if gnm.id == "":
#             # If the file was not parsed, the id will remain empty.
#             failed_files.append(filename)
#         else:
#             valid_files.append(filename)
#     return (genomes, valid_files, failed_files)
#
#
#
#




###
