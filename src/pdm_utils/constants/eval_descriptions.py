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
pfx2 = "A genome that needs to be replaced"
pfx3 = "The genome is expected to have"
pfx4 = "The genome is not expected to have"
pfx5 = "The new and old genomes are expected to have"
pfx6 = "The bundled data is expected to contain"

sfx1 = "that is already in the database"

EVAL_DESCRIPTIONS = {
    # Target genome to import
    "GNM-EVAL-001": (f"{pfx1} cannot have a PhageID {sfx1} or is ''"),
    "GNM-EVAL-002": (f"{pfx1} cannot have a Name {sfx1} or is ''."),
    "GNM-EVAL-003": (f"{pfx1} cannot have a nucleotide sequence {sfx1} or is ''."),
    "GNM-EVAL-004": (f"{pfx1} is not expected to have a 'final' annotation status."),
    "GNM-EVAL-005": (f"{pfx1} cannot have an Accession {sfx1}."),
    "GNM-EVAL-006": (f"{pfx2} must have a PhageID {sfx1}."),
    "GNM-EVAL-007": (f"{pfx2} must have a nucleotide sequence {sfx1}."),
    "GNM-EVAL-008": (f"{pfx2} is not expected to have a 'draft' annotation status."),
    "GNM-EVAL-009": (f"A genome with 'draft' annotation status is expected to contain the '{suffix}' suffix in the Name."),
    "GNM-EVAL-010": ("A genome with 'draft' annotation status is expected to have 0 CDS descriptions."),
    "GNM-EVAL-011": ("A genome with 'draft' annotation status is expected to have a '' Accession."),
    "GNM-EVAL-012": (f"A genome with 'final' annotation status is not expected to contain the '{suffix}' suffix in the Name."),
    "GNM-EVAL-013": ("A genome with 'final' annotation status is expected to have > 0 CDS descriptions."),
    "GNM-EVAL-014": (f"The genome is not expected to contain the '{suffix}' suffix in the PhageID."),
    "GNM-EVAL-015": (f"{pfx3} a valid Annotation Status: {statuses}."),
    "GNM-EVAL-016": (f"{pfx3} a valid Annotation Author: {authors}."),
    "GNM-EVAL-017": (f"{pfx3} a valid Retrieve Record: {records}."),
    "GNM-EVAL-018": (f"{pfx3} a Cluster that has been previously defined or be a Singleton."),
    "GNM-EVAL-019": (f"{pfx3} no Subcluster or have a Subcluster that has been previously defined."),
    "GNM-EVAL-020": (f"{pfx3} the following Translation Table: 11."),
    "GNM-EVAL-021": (f"{pfx3} a Host Genus that has been previously defined."),
    "GNM-EVAL-022": (f"{pfx3} a Cluster comprised strictly of alphabetic characters (e.g. 'A')."),
    "GNM-EVAL-023": (f"{pfx3} a Subcluster comprised of alphabetic characters followed by numeric characters (e.g. 'A15')."),
    "GNM-EVAL-024": (f"{pfx3} a taxonomically-consistent Cluster and Subcluster. If there is a Subcluster, it should be a member of the Cluster. If the genome is a Singleton, there should be no Subcluster.  "),
    "GNM-EVAL-025": (f"{pfx3} a date recorded."),
    "GNM-EVAL-026": (f"{pfx3} GC content equal to or greater than 0."),
    "GNM-EVAL-027": (f"{pfx3} GC content equal to or less than 100."),
    "GNM-EVAL-028": (f"{pfx3} a length greater than 0."),
    "GNM-EVAL-029": (f"{pfx3} 1 or more annotated CDS features."),
    "GNM-EVAL-030": (f"{pfx3} 0 features that are nested within another feature, that contain the same start and stop coordinates as another feature, or that contain the same stop coordinate as another feature."),
    "GNM-EVAL-031": (f"{pfx4} any ambiguous nucleotides or gaps."),
    "GNM-EVAL-032": (f"{pfx4} any PhageID typos in the Description field."),
    "GNM-EVAL-033": (f"{pfx4} any PhageID typos in the Source field."),
    "GNM-EVAL-034": (f"{pfx4} any PhageID typos in the Organism field."),
    "GNM-EVAL-035": (f"{pfx4} any host genus typos in the Description field."),
    "GNM-EVAL-036": (f"{pfx4} any host genus typos in the Source field."),
    "GNM-EVAL-037": (f"{pfx4} any host genus typos in the Organism field."),
    "GNM-EVAL-038": (f"Since annotation author is True, {pfx3.lower()} one or more specific author(s) listed in the References fields."),
    "GNM-EVAL-039": (f"Since annotation author is True, {pfx4.lower()} to have generic authors (e.g. 'Lastname, Firstname') listed in the References fields."),
    "GNM-EVAL-040": (f"Since annotation author is False, {pfx4.lower()} one or more specific author(s) listed in the References fields."),

    # Current genome in database
    "GNM2-EVAL-001": ("The genome to be replaced is only expected to have a 'draft' annotation status."),

    # Bundled data
    "BNDL-EVAL-001": (f"{pfx6} an import ticket with instructions regarding how to import the genome data."),
    "BNDL-EVAL-002": (f"{pfx6} genome data parsed from a file that will be imported."),
    "BNDL-EVAL-003": (f"{pfx6} genome data parsed from an import ticket that is not present within the file."),
    "BNDL-EVAL-004": (f"{pfx6} genome data parsed from PhagesDB that is not present within the file."),
    "BNDL-EVAL-005": (f"{pfx6} genome data parsed from the MySQL database."),
    "BNDL-EVAL-006": (f"{pfx6} paired genome data from a file and from a MySQL database."),
    "BNDL-EVAL-007": (f"MySQL statements generated from the bundled data are expected to be successfully executed in a MySQL database."),

    # Import ticket
    "TKT-EVAL-001": ("The import ticket is expected to have a unique TicketID."),
    "TKT-EVAL-002": ("The import ticket is expected to have a unique PhageID."),
    "TKT-EVAL-003": ("The import ticket is expected to have a valid type (e.g. Add or Replace)."),
    "TKT-EVAL-004": ("The import ticket is expected to have a valid description field (e.g. Product)."),
    "TKT-EVAL-005": ("The import ticket is expected to have a valid Evaluation Mode (e.g. 'final')."),
    "TKT-EVAL-006": ("The import ticket is expected to have a set of evaluation flags."),
    "TKT-EVAL-007": ("The import ticket is not expected to have an empty PhageID."),
    "TKT-EVAL-008": ("The import ticket is expected to have a compatible type and retain settings."),
    "TKT-EVAL-009": ("The import ticket is expected to have valid add data."),
    "TKT-EVAL-010": ("The import ticket is expected to have valid retain data."),
    "TKT-EVAL-011": ("The import ticket is expected to have valid retrieve data."),
    "TKT-EVAL-012": ("The import ticket is expected to have valid parse data."),

    # Source feature
    "SRC-EVAL-001": (f"{pfx4} any PhageID typos in the Organism qualifier of the Source feature."),
    "SRC-EVAL-002": (f"{pfx4} any host genus typos in the Organism qualifier of the Source feature."),
    "SRC-EVAL-003": (f"{pfx4} any host genus typos in the Host qualifier of the Source feature."),
    "SRC-EVAL-004": (f"{pfx4} any host genus typos in the Lab_Host qualifier of the Source feature."),

    # CDS feature
    "CDS-EVAL-001": ("The CDS feature's translation is not expected to contain any ambiguous amino acids or gaps."),
    "CDS-EVAL-002": ("The CDS feature's nucleotide and amino acid sequence are expected to be consistent."),
    "CDS-EVAL-003": ("The CDS feature is expected to contain a Translation qualifier."),
    "CDS-EVAL-004": ("The CDS feature's Translation_Table qualifier is expected to be 11."),
    "CDS-EVAL-005": ("The CDS feature is expected to have an integer start coordinate > 0."),
    "CDS-EVAL-006": ("The CDS feature's orientation is expected to be either F or R."),
    "CDS-EVAL-007": ("The CDS feature is expected to have a Locus_Tag qualifier."),
    "CDS-EVAL-008": ("The CDS feature is not expected to have any PhageID typos in the Locus_Tag qualifier."),
    "CDS-EVAL-009": ("The CDS feature is expected to have a Gene qualifier."),
    "CDS-EVAL-010": ("The CDS feature's Gene qualifier is expected to be an integer."),
    "CDS-EVAL-011": ("The CDS feature is expected to have the same integer in the Gene and Locus_Tag qualifiers."),
    "CDS-EVAL-012": ("The CDS feature is not expected to have a description in any qualifiers other than the selected qualifier, if the selected qualifier has no description."),
    "CDS-EVAL-013": ("The CDS feature is expected to have an integer stop coordinate > 0."),
    "CDS-EVAL-014": ("The CDS feature is expected to be composed of 1 or more parts."),

    # Paired genomes
    "GP-EVAL-001": (f"{pfx5} the same PhageID."),
    "GP-EVAL-002": (f"{pfx5} the same Nucleotide Sequence."),
    "GP-EVAL-003": (f"{pfx5} the same Length."),
    "GP-EVAL-004": (f"{pfx5} the same Cluster."),
    "GP-EVAL-005": (f"{pfx5} the same Subcluster."),
    "GP-EVAL-006": (f"{pfx5} the same Host Genus."),
    "GP-EVAL-007": (f"{pfx5} the same Annotation Author."),
    "GP-EVAL-008": (f"{pfx5} the same Translation Table."),
    "GP-EVAL-009": (f"{pfx5} the same Retrieve Record."),
    "GP-EVAL-010": (f"The new genome is expected to have a more recent date than the old genome."),
    "GP-EVAL-011": ("The new genome Phage Name is not expected to contain the '_Draft' suffix."),
    "GP-EVAL-012": ("The new genome Annotation Status is not expected to be 'draft'."),
    "GP-EVAL-013": (f"{pfx5} the same Name."),
    "GP-EVAL-014": (f"{pfx5} the same Annotation Status."),
    "GP-EVAL-015": (f"{pfx5} the same Accession."),

    # tRNA feature
    "TRNA-EVAL-001": ("The tRNA feature is expected to have an integer start coordinate > 0."),
    "TRNA-EVAL-002": ("The tRNA feature is expected to have an integer stop coordinate > 0."),
    "TRNA-EVAL-003": ("The tRNA feature is expected to be composed of just 1 part."),
    "TRNA-EVAL-004": ("The tRNA feature's orientation is expected to be either F or R."),
    "TRNA-EVAL-005": ("The tRNA feature is expected to have a Locus Tag qualifier."),
    "TRNA-EVAL-006": ("The tRNA feature is not expected to have any PhageID typos in the Locus Tag qualifier."),
    "TRNA-EVAL-007": ("The tRNA feature is expected to have a Gene qualifier."),
    "TRNA-EVAL-008": ("The tRNA feature's Gene qualifier is expected to be an integer."),
    "TRNA-EVAL-009": ("The tRNA feature is expected to have the same integer in the Gene and Locus Tag qualifiers."),
    "TRNA-EVAL-010": ("The tRNA feature is expected to have a valid isotype (amino acid) prediction."),
    "TRNA-EVAL-011": ("The tRNA feature is expected to have a correct isotype (amino acid) prediction."),
    "TRNA-EVAL-012": ("The tRNA feature is expected to have a valid anticodon prediction."),
    "TRNA-EVAL-013": ("The tRNA feature is expected to have a correct anticodon prediction."),
    "TRNA-EVAL-014": ("The tRNA feature is expected to be ~75bp long (60-100bp allowed)."),
    "TRNA-EVAL-015": ("The tRNA feature is expected to end in CCA, CC, or C."),
    "TRNA-EVAL-016": ("The tRNA feature is expected to have a correct orientation."),
    "TRNA-EVAL-017": ("The tRNA feature is expected to be found by at least one tRNA prediction tool."),
    "TRNA-EVAL-018": ("The tRNA feature is expected to have coordinates that match Aragorn or tRNAscan-SE."),

    # tmRNA feature
    "TMRNA-EVAL-001": ("The tmRNA feature is expected to have an integer start coordinate > 0."),
    "TMRNA-EVAL-002": ("The tmRNA feature is expected to have an integer stop coordinate > 0."),
    "TMRNA-EVAL-003": ("The tmRNA feature is expected to be composed of just 1 part."),
    "TMRNA-EVAL-004": ("The tmRNA feature's orientation is expected to be either F or R."),
    "TMRNA-EVAL-005": ("The tmRNA feature is expected to have a Locus Tag qualifier."),
    "TMRNA-EVAL-006": ("The tmRNA feature is not expected to have any PhageID typos in the Locus Tag qualifier."),
    "TMRNA-EVAL-007": ("The tmRNA feature is expected to have a Gene qualifier."),
    "TMRNA-EVAL-008": ("The tmRNA feature's Gene qualifier is expected to be an integer."),
    "TMRNA-EVAL-009": ("The tmRNA feature is expected to have the same integer in the Gene and Locus Tag qualifiers."),
    "TMRNA-EVAL-010": ("The tmRNA feature is expected to have a peptide tag prediction."),
    "TMRNA-EVAL-011": ("The tmRNA feature is expected to have a correct peptide tag prediction."),
    "TMRNA-EVAL-012": ("The tmRNA feature is expected to have a correct orientation."),
}
