"""Functions to interact with PhagesDB"""


from classes import Eval
from functions import basic
import constants
import urllib.request
import json


def parse_phagesdb_phage_name(data_dict):
    """Retrieve Phage Name from PhagesDB."""
    try:
        phage_name = data_dict['phage_name']
        eval_object = None

    except:
        phage_name = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Phage Name data from PhagesDB.")

    return (phage_name, eval_object)



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
        eval_object = None

    except:
        cluster = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Cluster data from PhagesDB.")

    return (cluster, eval_object)


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
        eval_object = None

    except:
        subcluster = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Subcluster data from PhagesDB.")

    return (subcluster, eval_object)


def parse_phagesdb_host(data_dict):
    """Retrieve Host from PhagesDB."""
    try:
        host = data_dict["isolation_host"]["genus"]
        eval_object = None

    except:
        host = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Host data from PhagesDB.")

    return (host, eval_object)


def parse_phagesdb_accession(data_dict):
    """Retrieve Accession from PhagesDB."""
    try:
        accession = data_dict["genbank_accession"]
        eval_object = None

    except:
        accession = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Accession data from PhagesDB.")

    return (accession, eval_object)


def retrieve_phagesdb_fasta(fastafile_url):
    """Retrieve fasta file from PhagesDB."""

    try:
        request = urllib.request.Request(fastafile_url)
        with urllib.request.urlopen(request) as response:
            fasta_data = response.read()
            fasta_data = fasta_data.decode("utf-8")
        eval_object = None

    except:
        fasta_data = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve fasta data from PhagesDB.")

    return (fasta_data, eval_object)




















# TODO unit test below




























#TODO unit test
def parse_phagesdb_sequence(data_dict):
    """Retrieve genome sequence from PhagesDB."""

    try:
        fastafile_url = data_dict["fasta_file"]
        retrieved_fasta_file = retrieve_phagesdb_fasta(fastafile_url)
        sequence = functions_general.parse_fasta_file(retrieved_fasta_file)
        eval_object = None

    except:
        sequence = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve genome sequence data from PhagesDB.")

    return (sequence, eval_object)










#TODO unit test
def parse_phagesdb_data(genome_obj,data_dict):
    """Parses a dictionary of genome data retrieved from PhagesDB into a
    Genome object.
    """

    list_of_results = []


    # Phage Name and Search Name
    phage_name, result1 = parse_phagesdb_phage_name(data_dict)
    genome_obj.set_phage_name(phage_name)
    list_of_results.append(result1)

    # Host
    host, result2 = parse_phagesdb_host(data_dict)
    genome_obj.set_host(host)
    list_of_results.append(result2)

    # Accession
    accession, result3 = parse_phagesdb_accession(data_dict)
    genome_obj.set_accession(accession)
    list_of_results.append(result3)

    # Cluster
    cluster, result4 = parse_phagesdb_cluster(data_dict)
    genome_obj.set_cluster(cluster)
    list_of_results.append(result4)


    #Subcluster
    subcluster, result5 = parse_phagesdb_subcluster(data_dict)
    genome_obj.set_subcluster(subcluster)
    list_of_results.append(result5)


    # Genome sequence
    pdb_sequence, result6 = parse_phagesdb_sequence(data_dict)
    genome_obj.sequence = pdb_sequence
    list_of_results.append(result6)


    list_of_evals = []
    for result in list_of_results:
        if result is not None:
            list_of_evals.append(result)

    #TODO do I need to return this object?
    return (genome_obj, list_of_evals)








# TODO complete function
# TODO unit test.
def retrieve_phagesdb_genome(phage_id):
    """Retrieve all data from PhagesDB for a specific phage."""


    phage_url = constants.API_PREFIX + \
                phage_id + \
                constants.API_SUFFIX


    try:
        online_data_json = request.urlopen(phage_url)
        online_data_dict = json.loads(online_data_json.read())

        #Returns a genome object
        phagesdb_genome = \
            phagesdb.parse_phagesdb_data(phagesdb_genome, online_data_dict)

        if ticket.host == "retrieve":
            ticket.host = phagesdb_genome.host
        if ticket.cluster == "retrieve":
            ticket.cluster = phagesdb_genome.cluster
        if ticket.subcluster == "retrieve":
            ticket.subcluster = phagesdb_genome.subcluster
        if ticket.accession == "retrieve":
            ticket.accession = phagesdb_genome.accession

    except:
        pass
































###
