"""Collection of constants and dictionaries used for maintaining the phage
database."""

from Bio.Alphabet import IUPAC
from datetime import datetime

IMPORT_TABLE_SIZE = 12
IMPORT_TABLE_DICT = {
    "id":"",
    "type":"",
    "description_field":"",
    "run_mode":"",
    "phage_id":"",
    "host_genus":"",
    "cluster":"",
    "subcluster":"",
    "annotation_status":"",
    "annotation_author":"",
    "accession":"",
    "retrieve_record":""
    }


NAME_SUFFIX = "_Draft"
ANNOTATION_STATUS_SET = set(["draft", "final", "gbk"])
ANNOTATION_AUTHOR_SET = set([0,1])
ANNOTATION_QC_SET = set([0,1])
RETRIEVE_RECORD_SET = set([0,1])
EMPTY_DATE = datetime.strptime('1/1/0001', '%m/%d/%Y')
LOCUS_TAG_PREFIX_SET = set(["SEA", "PBI", "PHIRE"])

# Common list of values that represent empty or null values.
EMPTY_SET = set(["",
                "none",
                "null",
                None,
                "not applicable",
                "na",
                "n/a",
                "0",
                0,
                EMPTY_DATE])


# Set up dna and protein alphabets to verify sequence integrity
DNA_ALPHABET = set(IUPAC.IUPACUnambiguousDNA.letters)
PROTEIN_ALPHABET = set(IUPAC.protein.letters)


# Set of all types of allowable ticket actions using this script.
# Add = add a new genome without removing another.
# Remove = delete a genome without adding another.
# Replace = delete a genome and replace it with another.
# Update = make changes to one or more fields related to a genome
# already present in the database (e.g. HostStrain, Cluster, Subcluster, etc.)
IMPORT_TICKET_TYPE_SET = set(["add", "replace",])

# Create set of most common gene description genbank qualifiers.
DESCRIPTION_FIELD_SET = set(["product", "note", "function"])


# TODO this is probably no longer needed.
# # Create list of potential host names to ignore.
# # This is primarily for databases that contain phages of all host phyla
# # and not just Actinobacteria.
# HOST_IGNORE = ['enterobacteria','phage','bacteriophage','cyanophage']




# List of names that represent authors that have control over
# a genome record annotations. This constant is stored as a list
# so that multiple names can be stored, if needed.
AUTHOR_SET = set(["hatfull"])


# TODO this is probably no longer needed.
# # Dictionary for storing authorship info.
# # 1 = list of authors that should be listed on a genome record.
# # 0 = 'gbk', representing a genome record that a group does not have
# # control over.
# AUTHOR_DICTIONARY = {0:set(['gbk']),1: AUTHOR_SET}


# Settings to access data through the PhagesDB API.
API_PREFIX = "https://phagesdb.org/api/phages/"
API_SUFFIX = "/?format=json"
API_HOST_GENERA = "https://phagesdb.org/api/host_genera/"
API_CLUSTERS = "https://phagesdb.org/api/clusters/"

# TODO this is probably no longer needed.
# # Set of valid file extensions for flat files to be evaluated.
# ADMISSIBLE_FILE_TYPES = set(["gb","gbf","gbk","txt"])

# Phage name typo correction dictionary.
# Key = Phage name as it is spelled in the GenBank-formatted record.
# Value = Phage name as it is spelled in PhagesDB and/or Phamerator, and thus
# how it should be spelled in the import ticket.
# The phage name parsed from the GenBank-formatted record gets
# reassigned this corrected name.
# Reasons for the exceptions are indicated.
PHAGE_NAME_DICT = {
    # PhagesDB is unable to handle underscores.
    'ATCC29399B_C':'ATCC29399BC',
    'ATCC29399B_T':'ATCC29399BT',

    # ELB20 was reported as 'ELB20' in the original publication, but spelled
    # 'phiELB20' in the GenBank record.
    'phiELB20':'ELB20',

    # Names are spelled differently in GenBank record compared to
    # the original publication.
    'P100_1':'P100.1',
    'P100_A':'P100A',

    # 'LeBron' was changed to 'Bron' by ICTV. They won't change it back.
    'Bron':'LeBron',

    # Inadvertent typos introduced at some point that GenBank won't
    # correct since ICTV uses the typos.
    'BBPiebs31':'BPBiebs31',
    'CaptnMurica':'CapnMurica',
    'Fionnbarth':'Fionnbharth',

    # Case changes implemented by ICTV (and now GenBank).
    'Baka':'BAKA',
    'CJW1':'Cjw1',
    'Dlane':'DLane',
    'Kssjeb':'KSSJEB',
    'Littlee':'LittleE',
    'Billknuckles':'BillKnuckles',
    'Packman':'PackMan',
    'Mrgordo':'MrGordo',
    'Ericb':'EricB',
    'lockley':'Lockley',
    'Heldan':'HelDan',
    'Ta17a':'TA17a'
    }

# Host genus typo dictionary.
# Key = Host genus as it is spelled in the GenBank-formatted record.
# Value = Host genus as it is spelled in PhagesDB and/or Phamerator, and thus
# how it should be spelled in the import ticket.
# The host genus parsed from the GenBank-formatted record gets
# reassigned this corrected name.
HOST_GENUS_DICT = {
    "Mycolicibacterium":"Mycobacterium"
    }



# Define run modes:

# TODO implement the 'import_locus_tag' option.

# Options that impact how data is processed but are not utilized
# for evaluation of specific parts of a flat file:

# 'check_replace':           Should unexpected genome replacements be reported?
# 'import_locus_tag':        Should locus_tags be imported?


# Options that are utilized during the evaluation stage:

# 'check_locus_tag':         Should the structure of locus_tags be checked?
# 'check_description_field': Should CDS descriptions in unexpected
#                            fields be reported?
# 'check_description':       Should unexpected CDS descriptions be reported?
# 'check_trna':              Should tRNA features be evaluated?
# 'check_id_typo':           Should genome ID typos be reported?
# 'check_host_typo':         Should host typos be reported?
# 'check_author':            Should unexpected authors be reported?
# 'check_gene':              Should the CDS 'gene' qualifier be evaluated?
# 'check_seq':               Should the nucleotide sequence be evaluated?

def _get_run_mode_base():
    dict = {
        "check_locus_tag":True,
        "check_description_field":True,
        "check_replace":True,
        "check_trna":True,
        "import_locus_tag":True,
        "check_id_typo":True,
        "check_host_typo":True,
        "check_author":True,
        "check_description":True,
        "check_gene":True
        }
    return dict

RUN_MODE_BASE = _get_run_mode_base()


# Auto-annotations.
def _get_run_mode_pecaan():
    new_dict = _get_run_mode_base()
    new_dict["check_locus_tag"] = False
    new_dict["check_trna"] = False
    new_dict["import_locus_tag"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_host_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    return new_dict
RUN_MODE_PECAAN = _get_run_mode_pecaan()

# Manual annotations.
def _get_run_mode_phagesdb():
    new_dict = _get_run_mode_base()
    new_dict["import_locus_tag"] = False
    return new_dict
RUN_MODE_PHAGESDB = _get_run_mode_phagesdb()

# SEA-PHAGES GenBank records.
def _get_run_mode_sea_auto():
    new_dict = _get_run_mode_base()
    new_dict["check_locus_tag"] = False
    new_dict["check_description_field"] = False
    new_dict["check_replace"] = False
    new_dict["check_trna"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    new_dict["check_gene"] = False
    return new_dict
RUN_MODE_SEA_AUTO = _get_run_mode_sea_auto()

# Non-SEA-PHAGES GenBank records.
def _get_run_mode_misc():
    new_dict = _get_run_mode_base()
    new_dict["check_locus_tag"] = False
    new_dict["check_replace"] = False
    new_dict["check_trna"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_host_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    new_dict["check_gene"] = False
    return new_dict
RUN_MODE_MISC = _get_run_mode_misc()

# Custom QC settings. User can select the settings, so it is initialized as
# a copy of the base run_mode. At the command line, the user can
# provide the customized combination of options.
RUN_MODE_CUSTOM = _get_run_mode_base()


# A dictionary that holds all the other run_mode dictionaries.
# Import tables will use the keys to retrieve the right combination
# of parameters. If new options needed to be created, they need to
# be added to this dictionary.
# The 'custom' dictionary enables the user to specify a unique
# combination of options at the command line.
RUN_MODES = {
    "pecaan":RUN_MODE_PECAAN,
    "phagesdb":RUN_MODE_PHAGESDB,
    "sea_auto":RUN_MODE_SEA_AUTO,
    "misc":RUN_MODE_MISC,
    "custom":RUN_MODE_CUSTOM
    }
