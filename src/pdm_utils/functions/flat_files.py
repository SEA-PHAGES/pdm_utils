"""Functions to interact with, use, and parse genomic data from
GenBank-formatted flat files."""

import pathlib
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, ExactPosition
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pdm_utils.classes import genome, cds, trna, source
from pdm_utils.functions import basic
from pdm_utils.constants import constants
from datetime import datetime
from pdm_utils.classes import genomepair



def retrieve_genome_data(filepath):
    """Retrieve data from a GenBank-formatted flat file.

    :param filepath:
        Path to GenBank-formatted flat file that will be parsed
        using Biopython.
    :type filepath: Path
    :returns:
        If there is only one record, a Biopython SeqRecord of parsed data.
        If the file cannot be parsed, or if there are multiple records,
        None value is returned.
    :rtype: SeqRecord
    """
    try:
        seqrecords = list(SeqIO.parse(filepath, "genbank"))
    except:
        seqrecords = []
    # filename = filepath.split("/")[-1]
    if len(seqrecords) == 0:
        print(f"There are no records in {filepath.name}.")
        seqrecord = None
    elif len(seqrecords) > 1:
        print(f"There are multiple records in {filepath.name}." )
        seqrecord = None
    else:
        seqrecord = seqrecords[0]
        return seqrecord


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
        tuple (start, stop, parts)
        WHERE
        start(int) is the first coordinate, regardless of strand.
        stop(int) is the second coordinate, regardless of strand.
        parts(int) is the number of open reading frames that define
        the feature.
    """
    start_position = None
    stop_position = None
    start = -1
    stop = -1
    parts = 0

    if (isinstance(seqfeature.location, FeatureLocation) or \
        isinstance(seqfeature.location, CompoundLocation)):

        if seqfeature.strand is None:
            pass
        elif isinstance(seqfeature.location, FeatureLocation):
            parts = 1
            start_position = seqfeature.location.start
            stop_position = seqfeature.location.end
        elif isinstance(seqfeature.location, CompoundLocation):
            parts = len(seqfeature.location.parts)

            # Skip this compound seqfeature if it is comprised of more
            # than two features (tricky to parse).
            if parts == 2:

                # Retrieve compound seqfeature positions based on strand.
                if seqfeature.strand == 1:
                    start_position = seqfeature.location.parts[0].start
                    stop_position = seqfeature.location.parts[1].end
                elif seqfeature.strand == -1:
                    start_position = seqfeature.location.parts[1].start
                    stop_position = seqfeature.location.parts[0].end
                else:
                    pass
            else:
                pass
        else:
            pass
    else:
        pass
    if isinstance(start_position, ExactPosition):
        start = int(start_position)
    if isinstance(stop_position, ExactPosition):
        stop = int(stop_position)
    return (start, stop, parts)


def parse_cds_seqfeature(seqfeature):
    """Parse data from a Biopython CDS SeqFeature object into a Cds object.

    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :param genome_id:

        An identifier for the genome in which the seqfeature
        is defined.

    :type genome_id: str
    :returns: A  pdm_utils Cds object
    :rtype: Cds
    """
    cds_ftr = cds.Cds()
    cds_ftr.seqfeature = seqfeature

    try:
        locus_tag = seqfeature.qualifiers["locus_tag"][0]
    except:
        locus_tag = ""
    cds_ftr.set_locus_tag(locus_tag)

    cds_ftr.set_orientation(seqfeature.strand, "fr_short", case = True)
    cds_ftr.start, cds_ftr.stop, cds_ftr.parts = parse_coordinates(seqfeature)

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
        translation_table = 0
    cds_ftr.set_translation_table(translation_table)

    try:
        cds_ftr.raw_product, cds_ftr.product = \
            basic.reformat_description(seqfeature.qualifiers["product"][0])
    except:
        cds_ftr.raw_product = ""
        cds_ftr.product = ""

    try:
        cds_ftr.raw_function, cds_ftr.function = \
            basic.reformat_description(seqfeature.qualifiers["function"][0])
    except:
        cds_ftr.raw_function = ""
        cds_ftr.function = ""

    try:
        cds_ftr.raw_note, cds_ftr.note = \
            basic.reformat_description(seqfeature.qualifiers["note"][0])
    except:
        cds_ftr.raw_note = ""
        cds_ftr.note = ""

    try:
        cds_ftr.gene = seqfeature.qualifiers["gene"][0]
    except:
        cds_ftr.gene = ""

    cds_ftr.set_name()
    return cds_ftr



def parse_source_seqfeature(seqfeature):
    """Parses a Biopython Source SeqFeature.

    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :param genome_id:

        An identifier for the genome in which the seqfeature
        is defined.

    :type genome_id: str
    :returns: A pdm_utils Source object
    :rtype: Source
    """
    src_ftr = source.Source()
    src_ftr.seqfeature = seqfeature
    start, stop, parts = parse_coordinates(seqfeature)
    src_ftr.start = start
    src_ftr.stop = stop

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

    src_ftr.parse_organism()
    src_ftr.parse_host()
    src_ftr.parse_lab_host()

    return src_ftr


def create_seqfeature_dictionary(seqfeature_list):
    """Create a dictionary of Biopython SeqFeature objects based on their type.

    From a list of all Biopython SeqFeatures derived from a GenBank-formatted
    flat file, create a dictionary of SeqFeatures based on their 'type'
    attribute.

    :param seqfeature_list: List of Biopython SeqFeatures
    :type seqfeature_list: list
    :param genome_id:

        An identifier for the genome in which the seqfeature
        is defined.

    :type genome_id: str
    :returns:

        A dictionary of Biopython SeqFeatures:
        Key: SeqFeature type (source, tRNA, CDS, other)
        Value: SeqFeature

    :rtype: dict
    """

    seqfeature_type_set = set()
    seqfeature_dict = {}
    for seqfeature in seqfeature_list:
        seqfeature_type_set.add(seqfeature.type)
    for type in seqfeature_type_set:
        sublist = []
        for index in range(len(seqfeature_list)):
            seqfeature = seqfeature_list[index]
            if seqfeature.type == type:
                sublist.append(seqfeature)
        seqfeature_dict[type] = sublist
    return seqfeature_dict

# TODO should this function be improved? Maybe create a more basic
# parse_flat_file() function - the only parameter is the seqrecord, and
# there is only minimal parsing and data processing.
# Then the parse_genome_data() function calls parse_flat_file(),
# and processes some data in specific ways.
def parse_genome_data(seqrecord, filepath=pathlib.Path(),
        translation_table=11, genome_id_field="_organism_name", gnm_type="",
        host_genus_field="_organism_host_genus"):
    """Parse data from a Biopython SeqRecord object into a Genome object.

    All Source, CDS, tRNA, and tmRNA features are parsed into their
    associated Source, Cds, Trna, and Tmrna objects.

    :param seqrecord: A Biopython SeqRecord object.
    :type seqrecord: SeqRecord
    :param filepath: A filename associated with the returned Genome object.
    :type filepath: Path
    :param translation_table:
        The applicable translation table for the genome's CDS features.
    :type translation_table: int
    :param genome_id_field:
        The SeqRecord attribute from which the unique genome
        identifier/name is stored.
    :type genome_id_field: str
    :param host_genus_field:
        The SeqRecord attribute from which the unique host genus
        identifier/name is stored.
    :type host_genus_field: str
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :returns: A pdm_utils Genome object.
    :rtype: Genome
    """

    # Keep track of the file from which the record is derived.
    gnm = genome.Genome()
    gnm.set_filename(filepath)
    gnm.type = gnm_type

    try:
        gnm.organism = seqrecord.annotations["organism"]
    except:
        gnm.organism = ""

    # Identifies host and phage name from organism field.
    gnm.parse_organism()

    try:
        # Since accessions are stored in a list, there may be more than
        # one accessions associated with this file.
        # The first accession in the list is assumed to be the most recent.
        accession = seqrecord.annotations["accessions"][0]
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
        gnm.source = seqrecord.annotations["source"]
    except:
        gnm.source = ""

    # Identifies host and phage name from record source field.
    gnm.parse_source()

    try:
        # The retrieved authors can be stored in multiple Reference elements.
        refs = seqrecord.annotations["references"]
        authors_list = []
        for ref in refs:
            # Note: Reference objects are instantiated with an empty
            # authors attribute. So if no authors are present in a Reference,
            # it will still concatenate an empty string, resulting in an
            # author_string = ";;;" etc. So only add the authors info if
            # it is not an empty string.
            if ref.authors != "":
                authors_list.append(ref.authors)
        authors_string = ";".join(authors_list)
        gnm.authors = authors_string
    except:
        gnm.authors = ""

    # Biopython requires the parsed record contains a sequence, so
    # no need to test whether the seq attribute is present or not.
    # Nucleotide sequence, length, and % GC.
    gnm.set_sequence(seqrecord.seq)

    try:
        date = seqrecord.annotations["date"]
        gnm.date = datetime.strptime(date, "%d-%b-%Y")
    except:
        gnm.date = constants.EMPTY_DATE


    # # Now that record fields are parsed, set the genome name, id,
    # # and host_genus.
    if genome_id_field != "":
        gnm.name = getattr(gnm, genome_id_field)
        gnm.set_id(value=gnm.name)
    else:
        # The seqrecord name and id are used if genome_id_field is empty.
        try:
            gnm.name = seqrecord.name
            # It appears that if name is not present, Biopython auto-populates
            # this attribute as "<unknown name>"
            if gnm.name == "<unknown name>":
                gnm.name = ""
        except:
            gnm.name = ""

        try:
            gnm.id = seqrecord.id
            # It appears that if id is not present, Biopython auto-populates
            # this attribute as "<unknown id>"
            if gnm.id == "<unknown id>":
                gnm.id = ""
        except:
            gnm.id = ""

    gnm.set_host_genus(attribute=host_genus_field)

    # Create lists of parsed features.
    # Note: Biopython instantiates the features attribute with
    # an empty list, so no need to test if features attribute is
    # present or not.
    seqfeature_dict = create_seqfeature_dictionary(seqrecord.features)

    cds_list = []
    if "CDS" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["CDS"]:
            cds_ftr = parse_cds_seqfeature(seqfeature)
            cds_ftr.genome_id = gnm.id
            cds_ftr.genome_length = gnm.length
            cds_ftr.set_nucleotide_sequence(parent_genome_seq=gnm.seq)
            cds_list.append(cds_ftr)

    source_list = []
    if "source" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["source"]:
            src_ftr = parse_source_seqfeature(seqfeature)
            src_ftr.genome_id = gnm.id
            source_list.append(src_ftr)

    # TODO unit test after functions are constructed.
    # TODO implement for trnas
    trna_list = []
    if "tRNA" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tRNA"]:
            trna_ftr = parse_trna_seqfeature(seqfeature)
            trna_list.append(trna_ftr)

    # TODO unit test after functions are constructed.
    tmrna_list = []
    if "tmrna" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tmrna"]:
            tmrna = parse_tmrna_seqfeature(seqfeature)
            tmrna_list.append(tmrna)

    gnm.translation_table = translation_table
    gnm.set_cds_features(cds_list)
    gnm.set_source_features(source_list)
    gnm.set_trna_features(trna_list)
    # gnm.set_tmrna_features(tmrna_list)

    # The Cds.id is constructed from the Genome.id and the Cds order.
    gnm.set_feature_ids(use_type=True, use_cds=True)
    gnm.set_feature_ids(use_type=True, use_source=True)

    # TODO set tRNA feature ids.
    #gnm.set_feature_ids(use_type=True, use_trna=True)
    return gnm


def genome_to_seqrecord(phage_genome):
    """Creates a SeqRecord object from a pdm_utils Genome object.

    :param phage_genome: A pdm_utils Genome object.
    :type phage_genome: Genome
    :returns: A BioPython SeqRecord object
    :rtype: SeqRecord
    """

    assert phage_genome != None,\
    "Genome object passed is None and not initialized"
    try:
        record = SeqRecord(phage_genome.seq)
        record.seq.alphabet = IUPAC.IUPACAmbiguousDNA()
    except AttributeError:
        print("Genome object failed to be converted to SeqRecord.",
              "Genome valid attribute 'seq' is required to",
              "convert to SeqRecord object.")
        raise
    record.name = phage_genome.name
    if phage_genome.accession != "":
        record.id = phage_genome.accession
    record.features = get_seqrecord_features(phage_genome)
    record.description = get_seqrecord_description(phage_genome)
    record.annotations=\
            get_seqrecord_annotations(phage_genome)

    return record

#TODO Owen unittest.
def get_seqrecord_features(phage_genome):
    """Helper function that uses Genome data to populate
    the features SeqRecord atribute

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        features is a list of SeqFeature objects parsed
        from cds objects
    """

    features = []
    for phage_cds in phage_genome.cds_features:
        features.append(phage_cds.seqfeature)

    return features

# TODO Owen unittest.
def get_seqrecord_description(phage_genome):
    """Helper function to construct a description SeqRecord attribute.

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        description is a formatted string parsed
        from genome data
    """

    description = (f"{phage_genome.host_genus} phage {phage_genome.id}"
                    ", complete genome")
    return description

# TODO Owen unittest.
def get_seqrecord_annotations(phage_genome):
    """Helper function that uses Genome data to populate
    the annotations SeqRecord attribute

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        annotations(dictionary) is a dictionary with
        the formatting of BioPython's SeqRecord
        annotations attribute
    """

    annotations = {"molecule type": "DNA",\
            "topology" : "linear",\
            "data_file_division" : "PHG",\
            "date" : "",\
            "accessions" : [],\
            "sequence_version" : "1",\
            "keyword" : [],\
            "source" : "",\
            "organism" : "",\
            "taxonomy" : [],\
            "comment": ()}
    annotations["date"] = phage_genome.date
    annotations["source"] =\
            "{} phage {}".format\
            (phage_genome.host_genus, phage_genome.id)
    annotations["organism"] =\
            "{} phage {}".format\
            (phage_genome.host_genus, phage_genome.name)
    annotations["taxonomy"].append("Viruses")
    annotations["taxonomy"].append("dsDNA Viruses")
    annotations["taxonomy"].append("Caudovirales")
    annotations["comment"] =\
            get_seqrecord_annotations_comments(phage_genome)
    return annotations

# TODO Owen unittest.
def get_seqrecord_annotations_comments(phage_genome):
    """Helper function that uses Genome data to populate
    the comment annotation attribute

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        cluster_comment, auto_generated_comment
        annotation_status_comment, qc_and_retrieval values
        (tuple) is a tuple with
        the formatting of BioPython's SeqRecord
        annotations comment attribute
    """
    if phage_genome.subcluster == "":
        cluster_comment = "Cluster: {}; Subcluster: None".format\
                (phage_genome.cluster)
    else:
        cluster_comment = "Cluster: {}; Subcluster: {}".format\
                (phage_genome.cluster, phage_genome.subcluster)
    auto_generated_comment =\
            "Auto-generated genome record from the MySQL database"
    annotation_status_comment =\
            "Annotation Status: {}; Annotation Author: {}".format\
            (phage_genome.annotation_status,\
            phage_genome.annotation_author)
    retrieval_value = \
            "RetrieveRecord: {}".format(phage_genome.retrieve_record)

    return (cluster_comment, auto_generated_comment,\
            annotation_status_comment, retrieval_value)


# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_trna_seqfeature(seqfeature):
    """Parses a Biopython tRNA SeqFeature and returns a pdm_utils Trna object.
    """
    trna_ftr = ""
    return trna_ftr


# TODO need to implement. Christian is developing tRNA object.
# TODO unit test.
def parse_tmrna_seqfeature(seqfeature):
    """Parses a Biopython tRNA SeqFeature and returns a pdm_utils Tmrna object.
    """
    tmrna = ""
    return tmrna


def create_fasta_seqrecord(header, sequence_string):
    """Create a fasta-formatted Biopython SeqRecord object.

    :param header: Description of the sequence.
    :type header: str
    :param sequence_string: Nucleotide sequence.
    :type sequence_string: str
    :returns: Biopython SeqRecord containing the nucleotide sequence.
    :rtype: SeqRecord
    """
    seq = Seq(sequence_string, alphabet=IUPAC.unambiguous_dna)
    seqrecord = SeqRecord(seq, description=header)
    return seqrecord
