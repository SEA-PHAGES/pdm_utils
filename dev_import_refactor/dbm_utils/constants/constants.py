"""Collection of constants and dictionaries used for maintaining the phage
database."""


from Bio.Alphabet import IUPAC
from datetime import datetime







IMPORT_TABLE_SIZE = 12

NAME_SUFFIX = "_Draft"

ANNOTATION_STATUS_SET = set(["draft", "final", "gbk"])
ANNOTATION_AUTHOR_SET = set([0,1])
ANNOTATION_QC_SET = set([0,1])
RETRIEVE_RECORD_SET = set([0,1])


EMPTY_DATE = datetime.strptime('1/1/0001', '%m/%d/%Y')

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
TICKET_TYPE_SET = set(["add", "remove", "replace", "update"])
IMPORT_TICKET_TYPE_SET = set(["add", "replace",])


# TODO is this constant still needed?
# Create set of most common gene description genbank qualifiers.
DESCRIPTION_FIELD_SET = set(["product", "note", "function"])


# Create list of potential host names to ignore.
# This is primarily for databases that contain phages of all host phyla
# and not just Actinobacteria.
HOST_IGNORE = ['enterobacteria','phage','bacteriophage','cyanophage']


LOCUS_TAG_PREFIX_SET = set(["SEA", "PBI", "PHIRE"])


# List of names that represent authors that have control over
# a genome record annotations. This constant is stored as a list
# so that multiple names can be stored, if needed.
AUTHOR_SET = set(["hatfull"])

# Dictionary for storing authorship info.
# 1 = list of authors that should be listed on a genome record.
# 0 = 'gbk', representing a genome record that a group does not have
# control over.
AUTHOR_DICTIONARY = {0:set(['gbk']),1: AUTHOR_SET}




# PhagesDB API to retrieve specific genome information.
API_PREFIX = "https://phagesdb.org/api/phages/"
API_SUFFIX = "/?format=json"


API_HOST_GENERA = "https://phagesdb.org/api/host_genera/"
API_CLUSTERS = "https://phagesdb.org/api/clusters/"

# Set of valid file extensions for flat files to be evaluated.
ADMISSIBLE_FILE_TYPES = set(["gb","gbf","gbk","txt"])



# Set of possible import ticket run modes.
RUN_MODE_SET = set(["phagesdb", "pecaan", "ncbi_auto", "ncbi_misc", "custom"])




# Phage name typo correction dictionary.
# Key = Phage name as it is spelled in the GenBank-formatted record.
# Value = Phage name as it is spelled in PhagesDB and/or Phamerator, and thus
# how it should be spelled in the import ticket.
# The phage name parsed from the GenBank-formatted record gets
# reassigned this corrected name.
# Reasons for the exceptions are indicated.

PHAGE_NAME_TYPO_DICT = {


    # PhagesDB is unable to handle underscores.
    'ATCC29399B_C':'ATCC29399BC',\
    'ATCC29399B_T':'ATCC29399BT',\


    # ELB20 was reported as 'ELB20' in the original publication, but spelled
    # 'phiELB20' in the GenBank record.
    'phiELB20':'ELB20',\


    # Names are spelled differently in GenBank record compared to
    # the original publication.
    'P100_1':'P100.1',\
    'P100_A':'P100A',\

    # 'LeBron' was changed to 'Bron' by ICTV. They won't change it back.
    'Bron':'LeBron',\


    # Inadvertent typos introduced at some point that GenBank won't
    # correct since ICTV uses the typos.
    'BBPiebs31':'BPBiebs31',\
    'CaptnMurica':'CapnMurica',\
    'Fionnbarth':'Fionnbharth',\

    # Case changes implemented by ICTV (and now GenBank).
    'Baka':'BAKA',\
    'CJW1':'Cjw1',\
    'Dlane':'DLane',\
    'Kssjeb':'KSSJEB',\
    'Littlee':'LittleE',\
    'Billknuckles':'BillKnuckles',\
    'Packman':'PackMan',\
    'Mrgordo':'MrGordo',\
    'Ericb':'EricB',\
    'lockley':'Lockley',\
    'Heldan':'HelDan',\
    'Ta17a':'TA17a'\
    }































# TODO need to clean up code below.




#Definitions for different run mode types

# Description of different run mode options:


# Options unrelated to flat file content:
# 'check_replace': Should unexpected genome replacements be prevented?
# 'import_locus_tag': Should locus_tags be imported?


# Options relating how to evaluate a flat file:
# 'filename': Should the filename be used to assign the Genome ID?
# 'check_locus_tag': Should the structure of locus_tags be checked?
# 'check_description_field': Does it matter if CDS descriptions are present in multiple fields?
# 'check_description': Should unexpected CDS descriptions within the flat file be reported?
# 'check_trna': Should tRNA features be evaluated?
# 'check_id_typo': Should genome ID typos within the flat file be reported?
# 'check_host_typo': Should host typos within the flat file be reported?
# 'check_author': Should unexpected authors within the flat file be reported?


RUN_MODE_BASE = {
    "use_filename":True,
    "check_locus_tag":False,
    "check_description_field":True,
    "check_replace":True,
    "check_trna":True,
    "import_locus_tag":True,
    "check_id_typo":True,
    "check_host_typo":True,
    "check_author":True,
    "check_description":True
    }

#Auto-annotations
def __get_run_mode_pecaan(dict=RUN_MODE_BASE):
    new_dict = dict.copy()
    new_dict["check_trna"] = False
    new_dict["import_locus_tag"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_host_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    return new_dict
RUN_MODE_PECAAN = __get_run_mode_pecaan()

#Manual annotations
def __get_run_mode_phagesdb(dict=RUN_MODE_BASE):
    new_dict = dict.copy()
    new_dict["import_locus_tag"] = False
RUN_MODE_PHAGESDB = __get_run_mode_phagesdb()

#SEA-PHAGES GenBank records
def __get_run_mode_sea_auto(dict=RUN_MODE_BASE):
    new_dict = dict.copy()
    new_dict["check_locus_tag"] = False
    new_dict["check_description_field"] = False
    new_dict["check_replace"] = False
    new_dict["check_trna"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    return new_dict
RUN_MODE_SEA_AUTO = __get_run_mode_sea_auto()

#Non-SEA-PHAGES GenBank records
def __get_run_mode_misc(dict=RUN_MODE_BASE):
    new_dict = dict.copy()
    new_dict["use_filename"] = False
    new_dict["check_locus_tag"] = False
    new_dict["check_replace"] = False
    new_dict["check_trna"] = False
    new_dict["check_id_typo"] = False
    new_dict["check_host_typo"] = False
    new_dict["check_author"] = False
    new_dict["check_description"] = False
    return new_dict
RUN_MODE_MISC = __get_run_mode_misc()

#Custom QC settings. User can select the settings, so it is initialized as
#an empty dictionary that only gets filled if there is a ticket indicating
#a custom set of parameters is needed.
RUN_MODE_CUSTOM = RUN_MODE_BASE.copy()


#A dictionary that holds all the other dictionaries.
#Import tables will use the keys to retrieve the right combination of parameters.
#If new options needed to be created, they need to be added to this dictionary.
#'none': reserved for import tickets that do not need a run mode specified (such as UPDATE tickets)
#'other': reserved for when users manually create import tickets and do not know
#which is the best option for their needs. Currently, this defaults to the 'phagesdb'
#run mode, since that is the most stringest criteria.
#'custom': reserved for when the user wants to specify a unique combination of options
#that are not reflected in the other run modes.
RUN_MODES = {
    "pecaan":RUN_MODE_PECAAN,
    "phagesdb":RUN_MODE_PHAGESDB,
    "sea_auto":RUN_MODE_SEA_AUTO,
    "misc":RUN_MODE_MISC,
    "custom":RUN_MODE_CUSTOM
    }
