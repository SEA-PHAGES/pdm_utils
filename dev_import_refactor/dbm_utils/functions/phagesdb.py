"""Functions to interact with PhagesDB"""


from classes import Eval
from functions import basic
from constants import constants
import urllib.request
import json


def parse_phagesdb_phage_name(data_dict):
    """Retrieve Phage Name from PhagesDB."""
    try:
        phage_name = data_dict['phage_name']
        result = "Phage name retrieved."
        status = "correct"

    except:
        phage_name = ""
        result = "Unable to retrieve Phage Name data from PhagesDB."
        status = "error"

    # TODO add an eval id?
    definition = "Parse the PhagesDB phage name."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)

    return (phage_name, eval)



def parse_phagesdb_cluster(data_dict):
    """Retrieve Cluster from PhagesDB.
    If the phage is clustered, 'pcluster' is a dictionary, and one key is
    the Cluster data (Cluster or 'Singleton').
    If for some reason no Cluster info is added at the time
    the genome is added to PhagesDB, 'pcluster' may automatically be
    set to NULL, which gets converted to "Unclustered" during retrieval.
    In Phamerator NULL means Singleton, and the long form
    "Unclustered" is invalid due to its character length,
    so this value is converted to 'UNK' ('Unknown').
    """

    try:
        if data_dict["pcluster"] is None:
            cluster = 'UNK'

        else:
            cluster = data_dict["pcluster"]["cluster"]

        result = "Cluster retrieved."
        status = "correct"

    except:
        cluster = ""
        result = "Unable to retrieve Cluster data from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Parse the PhagesDB cluster."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (cluster, eval)


def parse_phagesdb_subcluster(data_dict):
    """Retrieve Subcluster from PhagesDB.
    If for some reason no cluster info is added at the time
    the genome is added to PhagesDB, 'psubcluster' may automatically be
    set to NULL, which gets returned as None.
    If the phage is a Singleton, 'psubcluster' is None.
    If the phage is clustered but not subclustered, 'psubcluster' is None.
    If the phage is clustered and subclustered, 'psubcluster'
    is a dictionary, and one key is the Subcluster data."""

    try:
        if data_dict["psubcluster"] is None:
            subcluster = "none"
        else:
            subcluster = data_dict["psubcluster"]["subcluster"]
        result = "Subcluster retrieved."
        status = "correct"

    except:
        subcluster = ""
        result = "Unable to retrieve Subcluster data from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Parse the PhagesDB subcluster."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)

    return (subcluster, eval)


def parse_phagesdb_host(data_dict):
    """Retrieve Host from PhagesDB."""
    try:
        host = data_dict["isolation_host"]["genus"]
        result = "Host retrieved."
        status = "correct"

    except:
        host = ""
        result = "Unable to retrieve Host data from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Parse the PhagesDB host."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (host, eval)


def parse_phagesdb_accession(data_dict):
    """Retrieve Accession from PhagesDB."""
    try:
        accession = data_dict["genbank_accession"]
        result = "Accession retrieved."
        status = "correct"

    except:
        accession = ""
        result = "Unable to retrieve Accession data from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Parse the PhagesDB accession."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (accession, eval)



def parse_phagesdb_filename(data_dict):
    """Retrieve fasta filename from PhagesDB."""

    try:
        fastafile_url = data_dict["fasta_file"]
        result = "Fasta filename retrieved."
        status = "correct"

    except:
        fastafile_url = ""
        result = "Unable to retrieve fasta filename from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Parse the PhagesDB fasta filename."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (fastafile_url, eval)


def retrieve_phagesdb_fasta(fastafile_url):
    """Retrieve fasta file from PhagesDB."""

    try:
        request = urllib.request.Request(fastafile_url)
        with urllib.request.urlopen(request) as response:
            fasta_data = response.read()
            fasta_data = fasta_data.decode("utf-8")
        result = "Fasta data retrieved."
        status = "correct"

    except:
        fasta_data = ""
        result = "Unable to retrieve fasta data from PhagesDB."
        status = "error"


    # TODO add an eval id?
    definition = "Retrieve the PhagesDB fasta data."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (fasta_data, eval)




# TODO obsolete, once a biopython-based parser has been built?
def parse_fasta_file(fasta_file):
    """Parses sequence data from a fasta-formatted file.
    """
    # All sequence rows in the fasta file may not have equal widths,
    # so some processing of the data is required. If you split by newline,
    # the header is retained in the first list element.
    split_fasta_data = fasta_file.split('\n')

    header = ""
    sequence = ""

    if len(split_fasta_data) > 1:

        header = split_fasta_data[0]
        if header[0] == ">":
            header = header[1:] # remove '>' symbol.
            result = "Fasta record successfully parsed."
            status = "correct"

        else:
            result = "Record is not fasta-formatted."
            status = "error"

        header = header.strip() # remove any whitespace
        index = 1
        while index < len(split_fasta_data):

            # Strip off potential whitespace before appending, such as '\r'.
            sequence = sequence + split_fasta_data[index].strip()
            index += 1

    else:
        result = "Record is not fasta-formatted."
        status = "error"


    result = [header, sequence]


    # TODO add an eval id?
    definition = "Parse the PhagesDB fasta data."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (result, eval)






# TODO implement
# TODO unit test.
def parse_fasta_file2(fasta_file):
    """Parses sequence data from a fasta-formatted file using Biopython.
    """

    # TODO need to work out whether SeqIO can retrieve the data
    # directly from the server or whether the data needs to be first
    #retrieved and then parsed.

    pass









def parse_phagesdb_data(genome_obj,data_dict):
    """Parses a dictionary of genome data retrieved from PhagesDB into a
    Genome object.
    """

    list_of_results = []


    # Phage Name, PhageID and SearchID
    phage_name, result1 = parse_phagesdb_phage_name(data_dict)
    genome_obj.phage_name = phage_name
    genome_obj.set_phage_id(phage_name)
    list_of_results.append(result1)

    # Host
    host, result2 = parse_phagesdb_host(data_dict)
    genome_obj.set_host(host, "none_string")
    list_of_results.append(result2)

    # Accession
    accession, result3 = parse_phagesdb_accession(data_dict)
    genome_obj.set_accession(accession, "none_string")
    list_of_results.append(result3)

    # Cluster
    cluster, result4 = parse_phagesdb_cluster(data_dict)
    genome_obj.set_cluster(cluster)
    list_of_results.append(result4)

    #Subcluster
    subcluster, result5 = parse_phagesdb_subcluster(data_dict)
    genome_obj.set_subcluster(subcluster, "none_string")
    list_of_results.append(result5)

    # Fasta file URL
    fastafile_url, result6 = parse_phagesdb_filename(data_dict)
    genome_obj.filename = fastafile_url
    list_of_results.append(result6)


    # Fasta file record
    if genome_obj.filename != "":
        fasta_file, result7 = retrieve_phagesdb_fasta(genome_obj.filename)
        list_of_results.append(result7)

        # TODO unit test - not sure how to test this, since this function
        # retrieves and parses files from PhagesDB.
        # Genome sequence and parsed record
        if fasta_file != "":
            fasta_record, result8 = parse_fasta_file(fasta_file)
            genome_obj.seqrecord = fasta_record
            genome_obj.sequence = fasta_record[1]
            list_of_results.append(result8)


    # TODO this filtering step may not be ideal.
    list_of_evals = []
    for result in list_of_results:
        if result.status != "correct":
            list_of_evals.append(result)

    return list_of_evals



def retrieve_phagesdb_data(phage_url):
    """Retrieve all data from PhagesDB for a specific phage."""

    try:
        data_json = urllib.request.urlopen(phage_url)
        data_dict = json.loads(data_json.read())
        result = "PhagesDB data successfully retrieved."
        status = "correct"
    except:
        data_dict = {}
        result = "Unable to retrieve data from PhagesDB."
        status = "error"

    # TODO add an eval id?
    definition = "Retrieve PhagesDB data."
    eval = Eval.Eval(id = "", \
                    definition = definition, \
                    result = result, \
                    status = status)
    return (data_dict, eval)




def construct_phage_url(phage_name):
    """Create URL to retrieve phage-specific data from PhagesDB."""
    phage_url = constants.API_PREFIX + phage_name + constants.API_SUFFIX
    return phage_url


















###
