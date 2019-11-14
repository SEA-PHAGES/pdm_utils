"""Collection of constants and dictionaries used for maintaining the phage
database."""

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from datetime import datetime
from pathlib import Path




IMPORT_TABLE_STRUCTURE = {
    "order": [
                        "type",
                        "phage_id",
                        "description_field",
                        "run_mode",
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "accession",
                        "annotation_author",
                        "retrieve_record",
                        "annotation_status"],
    "required": set([
                        "type",
                        "phage_id"]),
    "optional": set([
                        "description_field",
                        "run_mode",
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "annotation_status",
                        "annotation_author",
                        "accession",
                        "retrieve_record"]),
    "valid_ticket": set([
                        "type",
                        "phage_id",
                        "description_field",
                        "run_mode"]),
    "valid_retain": set([
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "annotation_author",
                        "accession",
                        "retrieve_record"]),
    "valid_retrieve": set([
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "accession"]),
    "valid_add": set([
                        "host_genus",
                        "cluster",
                        "subcluster",
                        "annotation_author",
                        "annotation_status",
                        "accession",
                        "retrieve_record"]),
    "keywords": set(["retrieve", "retain", "none"])
}

NAME_SUFFIX = "_Draft"
ANNOTATION_STATUS_SET = set(["draft", "final", "unknown"])
ANNOTATION_AUTHOR_SET = set([0,1])
RETRIEVE_RECORD_SET = set([0,1])
EMPTY_DATE = datetime.strptime("1/1/0001 00:00:00", "%m/%d/%Y %H:%M:%S")

# Some locus tags have "PHIRE", but this should not be the case.
LOCUS_TAG_PREFIX_SET = set(["SEA", "PBI"])
EMPTY_GENOME_SEQ = Seq("", IUPAC.ambiguous_dna)

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

# Path to blastclust binary
BLASTCLUST_PATH = Path("~/bin/blast-2.2.14/bin").expanduser()


# Set up dna and protein alphabets to verify sequence integrity
DNA_ALPHABET = set(IUPAC.IUPACUnambiguousDNA.letters)
PROTEIN_ALPHABET = set(IUPAC.protein.letters)


# Set of all types of allowable ticket actions using this script.
# Add = add a new genome without removing another.
# Remove = delete a genome without adding another.
# Replace = delete a genome and replace it with another.
# Update = make changes to one or more fields related to a genome
# already present in the database (e.g. HostStrain, Cluster, Subcluster, etc.)
IMPORT_TICKET_TYPE_SET = set(["add", "replace"])

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


# Settings to access data through the PhagesDB API.
API_PREFIX = "https://phagesdb.org/api/phages/"
API_SUFFIX = "/?format=json"
API_HOST_GENERA = "https://phagesdb.org/api/host_genera/"
API_CLUSTERS = "https://phagesdb.org/api/clusters/"

# API at phagesdb has you specify how many results to return.
# 1 page at length 100000 will return everything.
SEQUENCED_PAGE = 1
SEQUENCED_SIZE = 100000
API_SEQUENCED = ("https://phagesdb.org/api/sequenced_phages/"
                 f"?page={SEQUENCED_PAGE}"
                 f"&page_size={SEQUENCED_SIZE}")



UNPHAMERATED_PHAGE_LIST = "https://phagesdb.org/data/unphameratedlist"


# Auto-annotations from PECAAN
PECAAN_PREFIX = "https://discoverdev.kbrinsgd.org/phameratoroutput/phage/"


# Phage name typo correction dictionary.
# Key = Phage name as it is spelled in the GenBank-formatted record.
# Value = Phage name as it is spelled in PhagesDB and/or Phamerator, and thus
# how it should be spelled in the import ticket.
# The phage name parsed from the GenBank-formatted record gets
# reassigned this corrected name.
# Reasons for the exceptions are indicated.
PHAGE_ID_DICT = {
    # PhagesDB is unable to handle underscores.
    "ATCC29399B_C": "ATCC29399BC",
    "ATCC29399B_T": "ATCC29399BT",

    # ELB20 was reported as 'ELB20' in the original publication, but spelled
    # 'phiELB20' in the GenBank record.
    "phiELB20": "ELB20",

    # Names are spelled differently in GenBank record compared to
    # the original publication.
    "P100_1": "P100.1",
    "P100_A": "P100A",

    # 'LeBron' was changed to 'Bron' by ICTV. They won't change it back.
    "Bron": "LeBron",

    # Inadvertent typos introduced at some point that GenBank won't
    # correct since ICTV uses the typos.
    "BBPiebs31": "BPBiebs31",
    "CaptnMurica": "CapnMurica",
    "Fionnbarth": "Fionnbharth",

    # Case changes implemented by ICTV (and now GenBank).
    "Baka": "BAKA",
    "Billknuckles": "BillKnuckles",
    "Crimd": "CrimD",
    "CJW1": "Cjw1",
    "Deadp": "DeadP",
    "Dlane": "DLane",
    "Dotproduct": "DotProduct",
    "Ericb": "EricB",
    "Gumbie": "GUmbie",
    "Heldan": "HelDan",
    "Jaws": "JAWS",
    "Joedirt": "JoeDirt",
    "Kssjeb": "KSSJEB",
    "Littlee": "LittleE",
    "Macncheese": "MacnCheese",
    "Mrgordo": "MrGordo",
    "Packman": "PackMan",
    "lockley": "Lockley",
    "Rockyhorror": "RockyHorror",
    "Ta17a": "TA17a"
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


# List of host_genus synonyms.
HOST_GENUS_SYNONYMS = [{"Mycobacterium", "Mycobacterio", "Mycolicibacterium"}]







# Phamerator server info
DB_HOST = "phamerator.webfactional.com"
DB_HOST_DIR = "/home/phamerator/webapps/htdocs/databases_Hatfull/"
DB_WEBSITE = "http://phamerator.webfactional.com/databases_Hatfull/"
