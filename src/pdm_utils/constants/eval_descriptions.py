"""Import evaluations descriptions."""
from pdm_utils.constants import constants


def get_string(value_set):
    string = ""
    for value in value_set:
        string = string + ", " + value
    return value

suffix = constants.NAME_SUFFIX
statuses = get_string(constants.ANNOTATION_STATUS_SET)
authors = get_string(constants.ANNOTATION_AUTHOR_SET)
records = get_string(constants.RETRIEVE_RECORD_SET)

eval_descriptions = {

    "GNM_001": ("A genome that needs to be added cannot have a PhageID that is already in the database or is ''"),
    "GNM_002": ("A genome that needs to be added cannot have a Name that is already in the database or is ''."),
    "GNM_003": ("A genome that needs to be added cannot have a nucleotide sequence that is already in the database or is ''."),
    "GNM_004": ("A genome that needs to be added is not expected to have a 'final' annotation status."),
    "GNM_005": ("A genome that needs to be added cannot have an Accession that is already in the database."),
    "GNM_006": ("A genome that needs to be replaced must have a PhageID that is already in the database."),
    "GNM_007": ("A genome that needs to be replaced must have a nucleotide sequence that is already in the database."),
    "GNM_008": ("A genome that needs to be replaced is not expected to have a 'draft' annotation status."),
    "GNM_009": (f"A genome with 'draft' annotation status is expected to contain the '{suffix}' suffix in the Name."),
    "GNM_010": ("A genome with 'draft' annotation status is expected to have 0 CDS descriptions."),
    "GNM_011": ("A genome with 'draft' annotation status is expected to have a '' Accession."),
    "GNM_012": (f"A genome with 'final' annotation status is not expected to contain the '{suffix}' suffix in the Name."),
    "GNM_013": ("A genome with 'draft' annotation status is expected to have > 0 CDS descriptions."),
    "GNM_014": (f"The genome is not expected to contain the '{suffix}' suffix in the PhageID."),
    "GNM_015": (f"The genome is expected to contain a valid Annotation Status: {statuses}."),
    "GNM_016": (f"The genome is expected to contain a valid Annotation Author: {authors}."),
    "GNM_017": (f"The genome is expected to contain a valid Retrieve Record: {records}."),
    "GNM_018": ("The genome is expected to be a Singleton or have a Cluster that has been previously defined."),
    "GNM_019": ("The genome is expected to have no Subcluster or have a Subcluster that has been previously defined."),
    "GNM_020": ("The genome is expected to have the following Translation Table: 11."),
    "GNM_021": ("The genome is expected to have a Host Genus that has been previously defined."),
    "GNM_022": ("The genome is expected to have a Cluster comprised strictly of alphabetic characters (e.g. 'A')."),
    "GNM_023": ("The genome is expected to have a Subcluster comprised of alphabetic characters followed by numeric characters (e.g. 'A15')."),
    "GNM_024": ("The genome is expected to have a taxonomically-consistent Cluster and Subcluster. If there is a Subcluster, it should be a member of the Cluster. If the genome is a Singleton, there should be no Subcluster.  "),
    "": (""),
    "": (""),
    "": (""),
    "": (""),
    "": (""),
    "": (""),
    "": (""),
    "": (""),
    "": (""),






}
