"""Functions to interact with, use, and parse genomic data from
GenBank-formatted flat files."""

import pathlib
from collections import OrderedDict
from datetime import datetime

from Bio import SeqIO
from Bio.SeqFeature import (
                SeqFeature, CompoundLocation, FeatureLocation, ExactPosition)
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pdm_utils.classes import genome, cds, trna, tmrna, source
from pdm_utils.constants import constants


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
        print(f"There are multiple records in {filepath.name}.")
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

    if (isinstance(seqfeature.location, FeatureLocation) or
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
    finally:
        cds_ftr.set_locus_tag(locus_tag, delimiter=None)

    cds_ftr.set_orientation(seqfeature.strand, "fr_short", case=True)
    cds_ftr.start, cds_ftr.stop, cds_ftr.parts = parse_coordinates(seqfeature)

    # Coordinate format for GenBank flat file features parsed by Biopython
    # are 0-based half open intervals.
    cds_ftr.coordinate_format = "0_half_open"

    # For translation, convert it to a Biopython Seq object.
    try:
        translation = seqfeature.qualifiers["translation"][0]
    except:
        translation = ""
    finally:
        translation = Seq(translation, Alphabet.IUPAC.protein)
        cds_ftr.set_translation(translation)

    cds_ftr.set_nucleotide_length(translation=True)

    try:
        translation_table = seqfeature.qualifiers["transl_table"][0]
    except:
        translation_table = 0
    finally:
        cds_ftr.set_translation_table(translation_table)

    try:
        product = seqfeature.qualifiers["product"][0]
    except:
        product = ""
    finally:
        cds_ftr.set_description_field("product", product)

    try:
        function = seqfeature.qualifiers["function"][0]
    except:
        function = ""
    finally:
        cds_ftr.set_description_field("function", function)

    try:
        note = seqfeature.qualifiers["note"][0]
    except:
        note = ""
    finally:
        cds_ftr.set_description_field("note", note)

    try:
        gene = seqfeature.qualifiers["gene"][0]
    except:
        gene = ""
    finally:
        cds_ftr.set_gene(gene)

    cds_ftr.set_name()
    return cds_ftr


# TODO: Christian unit test
def parse_trna_seqfeature(seqfeature):
    """
    Parse data from a Biopython tRNA SeqFeature object into a
    Trna object.
    :param seqfeature: Biopython SeqFeature
    :type seqfeature: SeqFeature
    :returns: a pdm_utils Trna object
    :rtype: Trna
    """
    trna_ftr = trna.Trna()
    trna_ftr.seqfeature = seqfeature

    try:
        locus_tag = seqfeature.qualifiers["locus_tag"][0]
    except (KeyError, IndexError):
        locus_tag = ""
    finally:
        trna_ftr.set_locus_tag(locus_tag, delimiter=None)

    trna_ftr.set_orientation(seqfeature.strand, "fr_short", True)
    trna_ftr.start, trna_ftr.stop, trna_ftr.parts = parse_coordinates(
                                                                    seqfeature)

    # Coordinate format for GenBank flat file features parsed by Biopython
    # are 0-based half open intervals.
    trna_ftr.coordinate_format = "0_half_open"

    trna_ftr.set_nucleotide_length(use_seq=True)

    try:
        product = seqfeature.qualifiers["product"][0]
    except (KeyError, IndexError):
        product = ""
    finally:
        trna_ftr.product = product

    try:
        note = seqfeature.qualifiers["note"][0]
    except (KeyError, IndexError):
        note = ""
    finally:
        trna_ftr.note = note

    try:
        gene = seqfeature.qualifiers["gene"][0]
    except (KeyError, IndexError):
        gene = ""
    finally:
        trna_ftr.gene = gene

    trna_ftr.set_name()
    return trna_ftr


# TODO: Christian unit test
def parse_tmrna_seqfeature(seqfeature):
    """
    Parses data from a BioPython tmRNA SeqFeature object into a
    Tmrna object.
    :param seqfeature: BioPython SeqFeature
    :type seqfeature: SeqFeature
    :return: pdm_utils Tmrna object
    :rtype: Tmrna
    """
    tmrna_ftr = tmrna.Tmrna()
    tmrna_ftr.seqfeature = seqfeature

    try:
        locus_tag = seqfeature.qualifiers["locus_tag"][0]
    except (KeyError, IndexError):
        locus_tag = ""
    finally:
        tmrna_ftr.set_locus_tag(locus_tag, delimiter=None)

    tmrna_ftr.set_orientation(seqfeature.strand, "fr_short", True)
    tmrna_ftr.start, tmrna_ftr.stop, tmrna_ftr.parts = parse_coordinates(
                                                                    seqfeature)

    # Coordinate format for GenBank flat file features parsed by Biopython
    # are 0-based half open intervals.
    tmrna_ftr.coordinate_format = "0_half_open"

    tmrna_ftr.set_nucleotide_length(use_seq=True)

    try:
        note = seqfeature.qualifiers["note"][0]
    except (KeyError, IndexError):
        note = ""
    finally:
        tmrna_ftr.note = note

    try:
        gene = seqfeature.qualifiers["gene"][0]
    except (KeyError, IndexError):
        gene = ""
    finally:
        tmrna_ftr.gene = gene

    tmrna_ftr.set_name()
    return tmrna_ftr


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
                      translation_table=11, genome_id_field="_organism_name",
                      gnm_type="", host_genus_field="_organism_host_genus"):
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
    finally:
        # Identifies host and phage name from organism field.
        gnm.parse_organism()

    try:
        # Since accessions are stored in a list, there may be more than
        # one accessions associated with this file.
        # The first accession in the list is assumed to be the most recent.
        accession = seqrecord.annotations["accessions"][0]
    except:
        accession = ""
    finally:
        gnm.set_accession(accession)

    try:
        gnm.description = seqrecord.description
        # It appears that if description is not present, Biopython
        # auto-populates this attribute as "<unknown description>"
        if gnm.description == "<unknown description>":
            gnm.description = ""
    except:
        gnm.description = ""
    finally:
        # Identifies host and phage name from description field.
        gnm.parse_description()

    try:
        gnm.source = seqrecord.annotations["source"]
    except:
        gnm.source = ""
    finally:
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

    trna_list = []
    if "tRNA" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tRNA"]:
            trna_ftr = parse_trna_seqfeature(seqfeature)
            trna_ftr.genome_id = gnm.id
            trna_ftr.genome_length = gnm.length
            trna_ftr.set_nucleotide_sequence(parent_genome_seq=gnm.seq)
            trna_ftr.set_nucleotide_length(use_seq=True)
            trna_ftr.parse_amino_acid()
            trna_ftr.parse_anticodon()
            trna_list.append(trna_ftr)

    tmrna_list = []
    if "tmRNA" in seqfeature_dict.keys():
        for seqfeature in seqfeature_dict["tmRNA"]:
            tmrna_ftr = parse_tmrna_seqfeature(seqfeature)
            tmrna_ftr.genome_id = gnm.id
            tmrna_ftr.genome_length = gnm.length
            tmrna_ftr.set_nucleotide_sequence(parent_genome_seq=gnm.seq)
            tmrna_ftr.set_nucleotide_length(use_seq=True)
            tmrna_ftr.parse_peptide_tag()
            tmrna_ftr.run_aragorn()
            tmrna_list.append(tmrna_ftr)

    gnm.translation_table = translation_table
    gnm.set_cds_features(cds_list)
    gnm.set_source_features(source_list)
    gnm.set_trna_features(trna_list)
    gnm.set_tmrna_features(tmrna_list)

    # The feature.id is constructed from the Genome.id and the feature order.
    gnm.set_feature_ids(use_type=True, use_cds=True)
    gnm.set_feature_ids(use_type=True, use_source=True)
    gnm.set_feature_ids(use_type=True, use_trna=True)
    gnm.set_feature_ids(use_type=True, use_tmrna=True)
    return gnm


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


# Needs unittests, however:
# Seqfeature retrieval and generation is clunky probably requires some
# over arching seqfeature generation.
# May delay unittests until structure is revamped
def genome_to_seqrecord(phage_genome):
    """Creates a SeqRecord object from a pdm_utils Genome object.

    :param phage_genome: A pdm_utils Genome object.
    :type phage_genome: Genome
    :returns: A BioPython SeqRecord object
    :rtype: SeqRecord
    """
    assert phage_genome is not None, \
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
    if phage_genome.accession == "" or phage_genome.accession is None:
        record.id = "".join(["DRAFT_", phage_genome.name])
    else:
        record.id = phage_genome.accession
    record.features = get_genome_seqrecord_features(phage_genome)
    record.description = get_genome_seqrecord_description(phage_genome)
    record.annotations = get_genome_seqrecord_annotations(phage_genome)

    return record


def cds_to_seqrecord(cds, parent_genome, gene_domains=[], desc_type="gb"):
    """Creates a SeqRecord object from a Cds and its parent Genome.

    :param cds: A populated Cds object.
    :type cds: Cds
    :param phage_genome: Populated parent Genome object of the Cds object.
    :param domains: List of domain objects populated with column attributes
    :type domains: list
    :param desc_type: Inteneded format of the CDS SeqRecord description.
    :type desc_type: str
    :returns: Filled Biopython SeqRecord object.
    :rtype: SeqRecord
    """
    record = SeqRecord(cds.translation)
    record.seq.alphabet = IUPAC.IUPACProtein()
    record.name = cds.id
    record.id = cds.id

    cds.set_seqfeature()

    source = f"{parent_genome.host_genus} phage {cds.genome_id}"
    source_feature = cds.create_seqfeature("source", 0,
                                           cds.translation_length, 1)
    source_feature.qualifiers["organism"] = [source]

    record.features = [source_feature]
    record.features.append(cds.create_seqfeature("Protein", 0,
                                                 cds.translation_length, 1))

    cds_feature = cds.create_seqfeature("CDS", 0, cds.translation_length, 1)
    format_cds_seqrecord_CDS_feature(cds_feature, cds, parent_genome)
    record.features.append(cds_feature)

    region_features = get_cds_seqrecord_regions(gene_domains, cds)
    for region_feature in region_features:
        record.features.append(region_feature)

    if desc_type == "fasta":
        description = ""

        product = cds.seqfeature.qualifiers["product"][0]
        if product > "":
            description = "".join([description, f"[product={product}]"])

        cluster = parent_genome.cluster
        if cluster > "":
            description = " ".join([description, f"[cluster={cluster}]"])

        host_genus = parent_genome.host_genus
        if host_genus > "":
            description = " ".join([description, f"[host_genus={host_genus}]"])

        if cds.locus_tag > "":
            description = " ".join([description, f"[locus={cds.locus_tag}]"])

        record.description = description
    elif desc_type == "gb":
        description = f"{source} gp {cds.name}"
        record.description = description

    record.annotations = get_cds_seqrecord_annotations(cds, parent_genome)

    return record


def get_genome_seqrecord_features(phage_genome):
    """Helper function that uses Genome data to populate
    the features SeqRecord atribute

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        features is a list of SeqFeature objects parsed
        from cds objects
    """

    source_feature = SeqFeature(FeatureLocation(0, phage_genome.length),
                                strand=1, type="source")
    source_feature.qualifiers = OrderedDict()
    source_feature.qualifiers["source"] = (f"{phage_genome.host_genus} phage "
                                           f"{phage_genome.name}")

    features = [source_feature]

    for phage_cds in phage_genome.cds_features:
        phage_cds.set_seqfeature(type="gene")
        features.append(phage_cds.seqfeature)
        phage_cds.set_seqfeature(type="CDS")
        features.append(phage_cds.seqfeature)

    for phage_trna in phage_genome.trna_features:
        phage_trna.set_seqfeature(type="gene")
        features.append(phage_trna.seqfeature)
        phage_trna.set_seqfeature()
        features.append(phage_trna.seqfeature)

    return features


def get_genome_seqrecord_description(phage_genome):
    """Helper function to construct a description SeqRecord attribute.

    :param phage_genome:
        Input a Genome object.
    :type phage_genome: genome
    :returns:
        description is a formatted string parsed
        from genome data
    """
    description = (f"{phage_genome.host_genus} phage {phage_genome.name}"
                   ", complete genome")
    return description


def get_genome_seqrecord_annotations(phage_genome):
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

    annotations = {"molecule type": "DNA",
                   "topology": "linear",
                   "data_file_division": "PHG",
                   "date": "",
                   "accessions": [],
                   "sequence_version": "1",
                   "keywords": [],
                   "source": "",
                   "organism": "",
                   "taxonomy": [],
                   "comment": ()}
    annotations["date"] = phage_genome.date
    annotations["keywords"] = ["complete_genome"]
    annotations["source"] = (
            f"{phage_genome.host_genus} phage {phage_genome.id}")
    annotations["organism"] = (
            f"{phage_genome.host_genus} phage {phage_genome.name}")
    annotations["taxonomy"].append("Viruses")
    annotations["taxonomy"].append("dsDNA Viruses")
    annotations["taxonomy"].append("Caudovirales")
    annotations["comment"] = get_genome_seqrecord_annotations_comments(
                                                                phage_genome)
    return annotations


def get_genome_seqrecord_annotations_comments(phage_genome):
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
        cluster_comment = "Cluster: {}; Subcluster: None".format(
                                                        phage_genome.cluster)
    else:
        cluster_comment = "Cluster: {}; Subcluster: {}".format(
                                phage_genome.cluster, phage_genome.subcluster)
    auto_generated_comment = (
            "Auto-generated genome record from the MySQL database")
    annotation_status_comment = (
            "Annotation Status: {}; Annotation Author: {}".format(
                                             phage_genome.annotation_status,
                                             phage_genome.annotation_author))
    retrieval_value = (
            "RetrieveRecord: {}".format(phage_genome.retrieve_record))

    return (cluster_comment, auto_generated_comment,
            annotation_status_comment, retrieval_value)


def get_cds_seqrecord_regions(gene_domains, cds):
    region_features = []
    for gene_domain in gene_domains:
        region_feature = cds.create_seqfeature("Region",
                                               gene_domain["QueryStart"],
                                               gene_domain["QueryEnd"], 1)
        region_feature.qualifiers["region_name"] = [gene_domain["Name"]]

        description = gene_domain["Description"]
        if description is None:
            description = ""
        else:
            description = description.decode("utf-8")
        region_feature.qualifiers["note"] = [description]

        region_feature.qualifiers["db_xref"] = ["CDD:"
                                                f"{gene_domain['DomainID']}"]
        region_features.append(region_feature)

    return region_features


def format_cds_seqrecord_CDS_feature(cds_feature, cds, parent_genome):
    cds_feature.qualifiers.pop("codon_start")
    cds_feature.qualifiers.pop("product")
    cds_feature.qualifiers.pop("translation")

    if cds.coordinate_format == "0_half_open":
        start = cds.seqfeature.location.start + 1
    elif cds.coordinate_format == "1_closed":
        start = cds.seqfeature.location.start
    else:
        start = cds.seqfeature.location.start

    coded_by = (f"{parent_genome.accession}:"
                f"{start}..""{cds.seqfeature.location.end}")
    if cds.seqfeature.strand == -1:
        coded_by = f"complement({coded_by})"

    cds_feature.qualifiers["coded by"] = [coded_by]


def get_cds_seqrecord_annotations(cds, parent_genome):
    """Function that creates a Cds SeqRecord annotations attribute dict.
    :param cds: A populated Cds object.
    :type cds: Cds
    :param phage_genome: Populated parent Genome object of the Cds object.
    :type phage_genome: Genome
    :returns: Formatted SeqRecord annotations dictionary.
    :rtype: dict{str}
    """
    annotations = {"topology": "linear",
                   "data_file_division": "PHG",
                   "date": "",
                   "accessions": [],
                   "sequence_version": "",
                   "keywords": [],
                   "source": "",
                   "organism": "",
                   "taxonomy": [],
                   "comment": ()}

    annotations["date"] = parent_genome.date
    annotations["organism"] = (f"{parent_genome.host_genus} phage "
                               f"{cds.genome_id}")
    annotations["source"] = f"Accession {parent_genome.accession}"

    annotations["taxonomy"].append("Viruses")
    annotations["taxonomy"].append("dsDNA Viruses")
    annotations["taxonomy"].append("Caudovirales")

    annotations["comment"] = get_cds_seqrecord_annotations_comments(cds)

    return annotations


def get_cds_seqrecord_annotations_comments(cds):
    """Function that creates a Cds SeqRecord comments attribute tuple.

    :param cds:
    :type cds:
    """
    pham_comment = f"Pham: {cds.pham_id}"
    auto_generated_comment = "Auto-generated CDS record from a MySQL database"

    return (pham_comment, auto_generated_comment)


def sort_seqrecord_features(seqrecord):
    """Function that sorts and processes the seqfeature objects of a seqrecord.

    :param seqrecord: Phage genome Biopython seqrecord object
    :type seqrecord: SeqRecord
    """

    try:
        def _sorting_key(seqfeature): return seqfeature.location.start
        seqrecord.features.sort(key=_sorting_key)
    except:
        if seqrecord is None:
            raise TypeError
        print("Genome seqrecord features unable to be sorted")
        pass
