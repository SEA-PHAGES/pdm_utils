"""Import evaluations descriptions."""
from pdm_utils.constants import constants


def get_string(value_set):
    value_list = list(value_set)
    value_list.sort()
    value_list = [str(x) for x in value_list]
    string = ", ".join(value_list)
    return string

suffix = constants.NAME_SUFFIX
statuses = get_string(constants.ANNOTATION_STATUS_SET)
authors = get_string(constants.ANNOTATION_AUTHOR_SET)
records = get_string(constants.RETRIEVE_RECORD_SET)


pfx1 = "A genome that needs to be added"
sfx1 = "that is already in the database"

pfx2 = "A genome that needs to be replaced"
pfx3 = "The genome is expected to have"
pfx4 = "The genome is not expected to have"
pfx5 = "The new and old genomes are expected to have"

pfx6 = "The bundled data is expected to contain"

EVAL_DESCRIPTIONS = {

    # Target genome to import
    "GNM_001": (f"{pfx1} cannot have a PhageID {sfx1} or is ''"),
    "GNM_002": (f"{pfx1} cannot have a Name {sfx1} or is ''."),
    "GNM_003": (f"{pfx1} cannot have a nucleotide sequence {sfx1} or is ''."),
    "GNM_004": (f"{pfx1} is not expected to have a 'final' annotation status."),
    "GNM_005": (f"{pfx1} cannot have an Accession {sfx1}."),
    "GNM_006": (f"{pfx2} must have a PhageID {sfx1}."),
    "GNM_007": (f"{pfx2} must have a nucleotide sequence {sfx1}."),
    "GNM_008": (f"{pfx2} is not expected to have a 'draft' annotation status."),
    "GNM_009": (f"A genome with 'draft' annotation status is expected to contain the '{suffix}' suffix in the Name."),
    "GNM_010": ("A genome with 'draft' annotation status is expected to have 0 CDS descriptions."),
    "GNM_011": ("A genome with 'draft' annotation status is expected to have a '' Accession."),
    "GNM_012": (f"A genome with 'final' annotation status is not expected to contain the '{suffix}' suffix in the Name."),
    "GNM_013": ("A genome with 'final' annotation status is expected to have > 0 CDS descriptions."),
    "GNM_014": (f"The genome is not expected to contain the '{suffix}' suffix in the PhageID."),
    "GNM_015": (f"{pfx3} a valid Annotation Status: {statuses}."),
    "GNM_016": (f"{pfx3} a valid Annotation Author: {authors}."),
    "GNM_017": (f"{pfx3} a valid Retrieve Record: {records}."),
    "GNM_018": (f"{pfx3} a Cluster that has been previously defined or be a Singleton."),
    "GNM_019": (f"{pfx3} no Subcluster or have a Subcluster that has been previously defined."),
    "GNM_020": (f"{pfx3} the following Translation Table: 11."),
    "GNM_021": (f"{pfx3} a Host Genus that has been previously defined."),
    "GNM_022": (f"{pfx3} a Cluster comprised strictly of alphabetic characters (e.g. 'A')."),
    "GNM_023": (f"{pfx3} a Subcluster comprised of alphabetic characters followed by numeric characters (e.g. 'A15')."),
    "GNM_024": (f"{pfx3} a taxonomically-consistent Cluster and Subcluster. If there is a Subcluster, it should be a member of the Cluster. If the genome is a Singleton, there should be no Subcluster.  "),
    "GNM_025": (f"{pfx3} a date recorded."),
    "GNM_026": (f"{pfx3} GC content equal to or greater than 0."),
    "GNM_027": (f"{pfx3} GC content equal to or less than 100."),
    "GNM_028": (f"{pfx3} a length greater than 0."),
    "GNM_029": (f"{pfx3} 1 or more annotated CDS features."),
    "GNM_030": (f"{pfx3} 0 features that are nested within another feature, that contain the same start and stop coordinates as another feature, or that contain the same stop coordinate as another feature."),
    "GNM_031": (f"{pfx4} any ambiguous nucleotides or gaps."),
    "GNM_032": (f"{pfx4} any PhageID typos in the Description field."),
    "GNM_033": (f"{pfx4} any PhageID typos in the Source field."),
    "GNM_034": (f"{pfx4} any PhageID typos in the Organism field."),
    "GNM_035": (f"{pfx4} any host genus typos in the Description field."),
    "GNM_036": (f"{pfx4} any host genus typos in the Source field."),
    "GNM_037": (f"{pfx4} any host genus typos in the Organism field."),
    "GNM_038": (f"Since annotation author is True, {pfx3.lower()} one or more specific author(s) listed in the References fields."),
    "GNM_039": (f"Since annotation author is True, {pfx4.lower()} to have generic authors (e.g. 'Lastname, Firstname') listed in the References fields."),
    "GNM_040": (f"Since annotation author is False, {pfx4.lower()} one or more specific author(s) listed in the References fields."),

    # Current genome in database
    "GNM2_001": ("The genome to be replaced is only expected to have a 'draft' annotation status."),

    # Bundled data
    "BNDL_001": (f"{pfx6} an import ticket with instructions regarding how to import the genome data."),
    "BNDL_002": (f"{pfx6} genome data parsed from a file that will be imported."),
    "BNDL_003": (f"{pfx6} genome data parsed from an import ticket that is not present within the file."),
    "BNDL_004": (f"{pfx6} genome data parsed from PhagesDB that is not present within the file."),
    "BNDL_005": (f"{pfx6} genome data parsed from the MySQL database."),
    "BNDL_006": (f"{pfx6} paired genome data from a file and from a MySQL database."),
    "BNDL_007": (f"MySQL statements generated from the bundled data are expected to be successfully executed in a MySQL database."),



    # Import ticket
    "TKT_001": ("The import ticket is expected to have a unique TicketID."),
    "TKT_002": ("The import ticket is expected to have a unique PhageID."),
    "TKT_003": ("The import ticket is expected to have a valid type (e.g. Add or Replace)."),
    "TKT_004": ("The import ticket is expected to have a valid description field (e.g. Product)."),
    "TKT_005": ("The import ticket is expected to have a valid Eval Mode (e.g. 'final')."),
    "TKT_006": ("The import ticket is expected to have a set of evaluation flags."),
    "TKT_007": ("The import ticket is not expected to have an empty PhageID."),
    "TKT_008": ("The import ticket is expected to have a compatible type and retain settings."),
    "TKT_009": ("The import ticket is expected to have valid add data."),
    "TKT_010": ("The import ticket is expected to have valid retain data."),
    "TKT_011": ("The import ticket is expected to have valid retrieve data."),
    "TKT_012": ("The import ticket is expected to have valid parse data."),

    # Source feature
    "SRC_001": (f"{pfx4} any PhageID typos in the Organism qualifier of the Source feature."),
    "SRC_002": (f"{pfx4} any host genus typos in the Organism qualifier of the Source feature."),
    "SRC_003": (f"{pfx4} any host genus typos in the Host qualifier of the Source feature."),
    "SRC_004": (f"{pfx4} any host genus typos in the Lab_Host qualifier of the Source feature."),

    # CDS feature
    "CDS_001": ("The CDS feature's translation is not expected to contain any ambiguous amino acids or gaps."),
    "CDS_002": ("The CDS feature's nucleotide and amino acid sequence are expected to be consistent."),
    "CDS_003": ("The CDS feature is expected to contain a Translation qualifier."),
    "CDS_004": ("The CDS feature's Translation_Table qualifier is expected to be 11."),
    "CDS_005": ("The CDS feature is expected to have an integer start coordinate > 0."),
    "CDS_006": ("The CDS feature's orientation is expected to be either F or R."),
    "CDS_007": ("The CDS feature is expected to have a Locus_Tag qualifier."),
    "CDS_008": ("The CDS feature is not expected to have any PhageID typos in the Locus_Tag qualifier."),
    "CDS_009": ("The CDS feature is expected to have a Gene qualifier."),
    "CDS_010": ("The CDS feature's Gene qualifier is expected to be an integer."),
    "CDS_011": ("The CDS feature is expected to have the same integer in the Gene and Locus_Tag qualifiers."),
    "CDS_012": ("The CDS feature is not expected to have a description in any qualifiers other than the selected qualifier, if the selected qualifier has no description."),
    "CDS_013": ("The CDS feature is expected to have an integer stop coordinate > 0."),
    "CDS_014": ("The CDS feature is expected to be composed of 1 or more parts."),


    # Paired genomes
    "GP_001": (f"{pfx5} the same PhageID."),
    "GP_002": (f"{pfx5} the same Nucleotide Sequence."),
    "GP_003": (f"{pfx5} the same Length."),
    "GP_004": (f"{pfx5} the same Cluster."),
    "GP_005": (f"{pfx5} the same Subcluster."),
    "GP_006": (f"{pfx5} the same Host Genus."),
    "GP_007": (f"{pfx5} the same Annotation Author."),
    "GP_008": (f"{pfx5} the same Translation Table."),
    "GP_009": (f"{pfx5} the same Retrieve Record."),
    "GP_010": (f"The new genome is expected to have a more recent date than the old genome."),
    "GP_011": ("The new genome Phage Name is not expected to contain the '_Draft' suffix."),
    "GP_012": ("The new genome Annotation Status is not expected to be 'draft'."),
    "GP_013": (f"{pfx5} the same Name."),
    "GP_014": (f"{pfx5} the same Annotation Status."),
    "GP_015": (f"{pfx5} the same Accession."),


    # TODO add Trna feature

    # TODO add Tmrna feature


}
