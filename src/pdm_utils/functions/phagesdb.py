"""Functions to interact with PhagesDB"""

import json
import pathlib
import urllib.request

from pdm_utils.classes import genome
from pdm_utils.constants import constants

def parse_phage_name(data_dict):
    """Retrieve Phage Name from PhagesDB.

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Name of the phage.
    :rtype: str
    """
    try:
        phage_name = data_dict["phage_name"]
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
    In the MySQL database NULL means Singleton, and the long form
    "Unclustered" is invalid due to its character length,
    so this value is converted to 'UNK' ('Unknown').

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Cluster of the phage.
    :rtype: str
    """
    try:
        if data_dict["pcluster"] is None:
            cluster = "UNK"
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

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Subcluster of the phage.
    :rtype: str
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
    """Retrieve host_genus from PhagesDB.

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Host genus of the phage.
    :rtype: str
    """
    try:
        host_genus = data_dict["isolation_host"]["genus"]
    except:
        host_genus = ""
    return host_genus


def parse_accession(data_dict):
    """Retrieve Accession from PhagesDB.

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Accession of the phage.
    :rtype: str
    """
    try:
        accession = data_dict["genbank_accession"]
    except:
        accession = ""
    return accession


def parse_fasta_filename(data_dict):
    """Retrieve fasta filename from PhagesDB.

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :returns: Name of the fasta file for the phage.
    :rtype: str
    """
    try:
        fastafile_url = data_dict["fasta_file"]
    except:
        fastafile_url = ""
    return fastafile_url


def retrieve_url_data(url):
    """Retrieve fasta file from PhagesDB.

    :param url: URL for data to be retrieved.
    :type url: str
    :returns: Data from the URL.
    :rtype: str
    """
    try:
        # gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1) ==> required for
        # urllib2.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED]
        # certificate verify failed (_ssl.c:590)> ==> creating new TLS
        # context tells urllib2 to ignore certificate chain
        # NOTE this is BAD SECURITY, prone to man-in-the-middle attacks
        request = urllib.request.Request(url)
        with urllib.request.urlopen(request) as response:
            data = response.read()
            data = data.decode("utf-8")
    except:
        print(f"Unable to retrieve data from {url}")
        data = ""
    return data


def parse_fasta_data(fasta_data):
    """Parses data returned from a fasta-formatted file.

    :param fasta_data: Data from a fasta file.
    :type fasta_data: str
    :returns:
        tuple (header, sequence)
        WHERE
        header(str) is the first line parsed from the parsed file.
        sequence(str) is the nucleotide sequence parsed from the file.
    :rtype: tuple
    """
    # All sequence rows in the fasta file may not have equal widths,
    # so some processing of the data is required. If you split by newline,
    # the header is retained in the first list element.
    split_fasta_data = fasta_data.split('\n')
    header = ""
    sequence = ""
    if len(split_fasta_data) > 1:
        header = split_fasta_data[0]
        if header[0] == ">":
            header = header[1:] # Remove '>' symbol.
        header = header.strip() # Remove any whitespace
        for index in range(1, len(split_fasta_data)):
            # Strip off potential whitespace before appending, such as '\r'.
            sequence = sequence + split_fasta_data[index].strip()
    return (header, sequence)


def parse_genome_data(data_dict, gnm_type="", seq=False):
    """Parses a dictionary of PhagesDB genome data into a pdm_utils Genome object.

    :param data_dict: Dictionary of data retrieved from PhagesDB.
    :type data_dict: dict
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :param seq: Indicates whether the genome sequence should be retrieved.
    :type seq: bool
    :returns: A pdm_utils Genome object with the parsed data.
    :rtype: Genome
    """
    gnm = genome.Genome()
    gnm.type = gnm_type

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
    gnm.set_subcluster(subcluster)

    # Fasta file URL
    fastafile_url = parse_fasta_filename(data_dict)
    fastafile_url_path = pathlib.Path(fastafile_url)
    gnm.set_filename(fastafile_url_path)

    # Fasta file record
    # if fastafile_url != "":
    if (fastafile_url != "" and seq == True):
        fasta_file = retrieve_url_data(fastafile_url)

        # TODO unit test - not sure how to test this, since this function
        # retrieves and parses files from PhagesDB.
        # Genome sequence and parsed record
        if fasta_file != "":
            header, seq = parse_fasta_data(fasta_file)
            gnm.set_sequence(seq)
            gnm.description = header
            gnm.parse_description()

    gnm.misc = data_dict
    return gnm

def retrieve_genome_data(phage_url):
    """Retrieve all data from PhagesDB for a specific phage.

    :param phage_url: URL for data pertaining to a specific phage.
    :type phage_url: str
    :returns: Dictionary of data parsed from the URL.
    :rtype: dict
    """
    try:
        data_json = urllib.request.urlopen(phage_url)
        data_dict = json.loads(data_json.read())
        # TODO should data_dict.close() be called after retrieving data?
    except:
        data_dict = {}
    return data_dict


def construct_phage_url(phage_name):
    """Create URL to retrieve phage-specific data from PhagesDB.

    :param phage_name: Name of the phage of interest.
    :type phage_name: str
    :returns: URL pertaining to the phage.
    :rtype: str
    """
    phage_url = constants.API_PREFIX + phage_name + constants.API_SUFFIX
    return phage_url


def get_genome(phage_id, gnm_type="", seq=False):
    """Get genome data from PhagesDB.

    :param phage_id: The name of the phage to be retrieved from PhagesDB.
    :type phage_id: str
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :param seq: Indicates whether the genome sequence should be retrieved.
    :type seq: bool
    :returns:
        A pdm_utils Genome object with the parsed data.
        If not genome is retrieved, None is returned.
    :rtype: Genome
    """
    phage_url = construct_phage_url(phage_id)
    data_dict = retrieve_genome_data(phage_url)
    if len(data_dict.keys()) != 0:
        gnm = parse_genome_data(data_dict, gnm_type=gnm_type, seq=seq)
    else:
        gnm = None
    return gnm


# TODO unittest.
def get_phagesdb_data(url):
    """Retrieve all sequenced genome data from PhagesDB.

    :param url: URL to connect to PhagesDB API.
    :type url: str
    :returns:
        List of dictionaries, where each dictionary contains
        data for each phage. If a problem is encountered during retrieval,
        an empty list is returned.
    :rtype: list
    """
    data_json = urllib.request.urlopen(url)
    # Response is a bytes object that json.loads can't read without first
    # being decoded to a UTF-8 string.
    data_dict = json.loads(data_json.read().decode("utf-8"))
    data_json.close()

    # Returned dict:
    # Keys:
    # count = # of phages
    # results = list of sequenced phage data
    data_list = data_dict["results"]
    diff = len(data_list) - data_dict["count"]
    if diff != 0:
        print("\nUnable to retrieve all phage data from PhagesDB due to "
              "default parameters.")
        data_list = []
    return data_list


# TODO unittest.
def parse_genomes_dict(data_dict, gnm_type="", seq=False):
    """Returns a dictionary of pdm_utils Genome objects

    :param data_dict:
        Dictionary of dictionaries.
        Key = PhageID.
        Value = Dictionary of genome data retrieved from PhagesDB.
    :type data_dict: dict
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :param seq: Indicates whether the genome sequence should be retrieved.
    :type seq: bool
    :returns:
        Dictionary of pdm_utils Genome object.
        Key = PhageID.
        Value = Genome object.
    :rtype: dict
    """
    genome_dict = {}
    for key in data_dict.keys():
        gnm = parse_genome_data(data_dict[key], gnm_type=gnm_type, seq=seq)
        genome_dict[gnm.id] = gnm
    return genome_dict




def retrieve_data_list(url):
    """Retrieve list of data from PhagesDB.

    :param url: A URL from which to retrieve data.
    :type url: str
    :returns: A list of data retrieved from the URL.
    :rtype: list
    """
    try:
        data_json = urllib.request.urlopen(url)
        data_list = json.loads(data_json.read())
    except:
        data_list = []
    return data_list


def create_host_genus_set(url=constants.API_HOST_GENERA):
    """Create a set of host genera currently in PhagesDB.

    :param url: A URL from which to retrieve host genus data.
    :type url: str
    :returns: All unique host genera listed on PhagesDB.
    :rtype: set
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

    :param url: A URL from which to retrieve cluster and subcluster data.
    :type url: str
    :returns:
        tuple (cluster_set, subcluster_set)
        WHERE
        cluster_set(set) is a set of all unique clusters on PhagesDB.
        subcluster_set(set) is a set of all unique subclusters on PhagesDB.
    :rtype: tuple
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


# TODO unittest.
def get_unphamerated_phage_list(url):
    """Retreive list of unphamerated phages from PhagesDB.

    :param url:
        A URL from which to retrieve a list of PhagesDB genomes that are not
        in the most up-to-date instance of the Actinobacteriophage
        MySQL database.
    :type url: str
    :returns: List of PhageIDs.
    :rtype: list
    """
    # Retrieved file is a tab-delimited text file.
    # Each row is a newly-sequenced phage.
    response = urllib.request.urlopen(url)
    processed_list = []
    for new_phage in response:
        new_phage = new_phage.strip()  # Remove \t at end of each row
        new_phage = new_phage.decode("utf-8")  # convert bytes object to str
        processed_list.append(new_phage)
    return processed_list
