"""Functions to interact with PhagesDB"""


from pdm_utils.pipelines.db_import import evaluate
from pdm_utils.classes import genome
from pdm_utils.classes import genomepair
from pdm_utils.functions import basic
from pdm_utils.constants import constants
import urllib.request
import json
from pdm_utils.functions import misc


def parse_phage_name(data_dict):
    """Retrieve Phage Name from PhagesDB."""
    try:
        phage_name = data_dict['phage_name']
    except:
        phage_name = ""
    return phage_name


def parse_cluster(data_dict):
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
    except:
        cluster = ""
    return cluster


def parse_subcluster(data_dict):
    """Retrieve Subcluster from PhagesDB.

    If for some reason no cluster info is added at the time
    the genome is added to PhagesDB, 'psubcluster' may automatically be
    set to NULL, which gets returned as None.
    If the phage is a Singleton, 'psubcluster' is None.
    If the phage is clustered but not subclustered, 'psubcluster' is None.
    If the phage is clustered and subclustered, 'psubcluster'
    is a dictionary, and one key is the Subcluster data.
    """

    try:
        if data_dict["psubcluster"] is None:
            subcluster = "none"
        else:
            subcluster = data_dict["psubcluster"]["subcluster"]
    except:
        subcluster = ""
    return subcluster


def parse_host_genus(data_dict):
    """Retrieve host_genus from PhagesDB."""
    try:
        host_genus = data_dict["isolation_host"]["genus"]
    except:
        host_genus = ""
    return host_genus


def parse_accession(data_dict):
    """Retrieve Accession from PhagesDB."""
    try:
        accession = data_dict["genbank_accession"]
    except:
        accession = ""
    return accession


def parse_fasta_filename(data_dict):
    """Retrieve fasta filename from PhagesDB."""

    try:
        fastafile_url = data_dict["fasta_file"]
    except:
        fastafile_url = ""
    return fastafile_url


def retrieve_fasta_data(fastafile_url):
    """Retrieve fasta file from PhagesDB."""

    try:
        request = urllib.request.Request(fastafile_url)
        with urllib.request.urlopen(request) as response:
            fasta_data = response.read()
            fasta_data = fasta_data.decode("utf-8")
    except:
        fasta_data = ""
    return fasta_data


def parse_fasta_data(fasta_string):
    """Parses data returned from a fasta-formatted file."""
    # All sequence rows in the fasta file may not have equal widths,
    # so some processing of the data is required. If you split by newline,
    # the header is retained in the first list element.
    split_fasta_data = fasta_string.split('\n')
    header = ""
    sequence = ""

    if len(split_fasta_data) > 1:
        header = split_fasta_data[0]
        if header[0] == ">":
            header = header[1:] # Remove '>' symbol.

        header = header.strip() # Remove any whitespace
        index = 1
        while index < len(split_fasta_data):
            # Strip off potential whitespace before appending, such as '\r'.
            sequence = sequence + split_fasta_data[index].strip()
            index += 1
    result = (header, sequence)
    return result


def parse_genome_data(data_dict):
    """Parses a dictionary of genome data retrieved from PhagesDB into a
    Genome object.
    """

    gnm = genome.Genome()
    gnm.type = "phagesdb"

    # Phage Name, PhageID
    phage_name = parse_phage_name(data_dict)
    gnm.name = phage_name
    gnm.set_id(value=phage_name)

    # Host
    host_genus = parse_host_genus(data_dict)
    gnm.set_host_genus(host_genus, "empty_string")

    # Accession
    accession = parse_accession(data_dict)
    gnm.set_accession(accession, "empty_string")

    # Cluster
    cluster = parse_cluster(data_dict)
    gnm.set_cluster(cluster)

    #Subcluster
    subcluster = parse_subcluster(data_dict)
    gnm.set_subcluster(subcluster, "empty_string")

    # Fasta file URL
    fastafile_url = parse_fasta_filename(data_dict)
    gnm.filename = fastafile_url

    # Fasta file record
    if gnm.filename != "":
        fasta_file = retrieve_fasta_data(gnm.filename)

        # TODO unit test - not sure how to test this, since this function
        # retrieves and parses files from PhagesDB.
        # Genome sequence and parsed record
        if fasta_file != "":
            header, seq = parse_fasta_data(fasta_file)
            gnm.set_sequence(seq)
            gnm.description = header
            gnm.parse_description()

    # TODO not sure if these evaluations should be in this function or not.
    evaluate.check_phagesdb_genome(gnm, set([""]))
    return gnm


def retrieve_genome_data(phage_url):
    """Retrieve all data from PhagesDB for a specific phage."""

    try:
        data_json = urllib.request.urlopen(phage_url)
        data_dict = json.loads(data_json.read())
    except:
        data_dict = {}
    return data_dict


def construct_phage_url(phage_name):
    """Create URL to retrieve phage-specific data from PhagesDB."""
    phage_url = constants.API_PREFIX + phage_name + constants.API_SUFFIX
    return phage_url


def copy_data_from(bndl, type, flag="retrieve"):
    """Copy data from a 'phagesdb' genome object.

    If a genome object stored in the Bundle object has
    attributes that are set to be 'retrieved' and auto-filled,
    retrieve the data from PhagesDB to complete the genome.
    The 'type' parameter indicates the type of genome that may need
    to be populated from PhagesDB.
    """

    if type in bndl.genome_dict.keys():
        genome1 = bndl.genome_dict[type]
        genome1.set_value_flag(flag)

        if genome1._value_flag:
            phage_url = construct_phage_url(genome1.id)
            data_dict = retrieve_genome_data(phage_url)

            # If there was an error with retrieving data from PhagesDB,
            # an empty dictionary is returned.
            if len(data_dict.keys()) != 0:
                genome2 = parse_genome_data(data_dict)
                bndl.genome_dict[genome2.type] = genome2

                # Copy all retrieved data and add to Bundle object.
                genome_pair = genomepair.GenomePair()
                genome_pair.genome1 = genome1
                genome_pair.genome2 = genome2
                genome_pair.copy_data("type", genome2.type, genome1.type, flag)
                bndl.set_genome_pair(genome_pair, genome1.type, genome2.type)

        # Now record an error if there are still fields
        # that need to be retrieved.
        genome1.set_value_flag(flag)
        genome1.check_value_flag()


def retrieve_data_list(url):
    """Retrieve list of data from PhagesDB."""
    try:
        data_json = urllib.request.urlopen(url)
        data_list = json.loads(data_json.read())
    except:
        data_list = []
    return data_list


def create_host_genus_set(url=constants.API_HOST_GENERA):
    """Create a set of host genera currently in PhagesDB.

    The parameter is a list, and each element is a dictionary of data
    pertaining to a different host genus.
    """

    try:
        output = retrieve_data_list(url)
    except:
        output = []
    host_genera_set = set()
    for genus_dict in output:
        try:
            host_genera_set.add(genus_dict["genus_name"])
        except:
            pass
    return host_genera_set


def create_cluster_subcluster_sets(url=constants.API_CLUSTERS):
    """Create sets of clusters and subclusters currently in PhagesDB.

    The parameter is a list, and each element is a dictionary of data
    pertaining to a different cluster.
    """

    try:
        output = retrieve_data_list(url)
    except:
        output = []

    cluster_set = set()
    subcluster_set = set()
    for data in output:
        try:
            cluster_set.add(data["cluster"]) # This set contains 'Singleton'.
            try:
                subclusters_list = data["subclusters_set"]
                subcluster_set = subcluster_set | set(subclusters_list)
            except:
                pass
        except:
            pass
    return (cluster_set, subcluster_set)

















###
