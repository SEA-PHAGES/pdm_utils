"""Functions to interact with PhagesDB"""


import Eval
import functions_general



#TODO complete function

def retrieve_phagesdb_data(phage_url):
    """Retrieve all data from PhagesDB for a specific phage."""
    #TODO retrieve json data, then convert to dictionary
    pass





#TODO unit test
def retrieve_phagesdb_fasta(fastafile_url):
    """Retrieve fasta file from PhagesDB."""

    response = urllib2.urlopen(fastafile_url)
    retrieved_fasta_file = response.read()
    response.close()

    return retrieved_fasta_file



#TODO unit test
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



#TODO unit test
def parse_phagesdb_cluster(data_dict):
    """Retrieve Cluster from PhagesDB.
    On PhagesDB, phages may have a Cluster and no Subcluster info
    (which is set to None). If the phage has a Subcluster,
    it should also have a Cluster. If by accident no Cluster or
    Subcluster info is added at the time the genome is added to
    PhagesDB, the Cluster may automatically be set to NULL,
    which gets converted to "Unclustered" during retrieval.
    This is problematic because in Phamerator NULL means Singleton,
    and the long form "Unclustered" will be filtered out later in the script
    due to its character length, so it needs to be abbreviated.
     """

    try:
        if data_dict['pcluster'] is None:
            cluster = 'UNK'
        else:
            cluster = data_dict['pcluster']['cluster']
        eval_object = None

    except:
        cluster = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Cluster data from PhagesDB.")

    return (cluster, eval_object)



#TODO unit test
def parse_phagesdb_subcluster(data_dict):
    """Retrieve Subcluster from PhagesDB. Subcluster could be empty
    if by error no Cluster or Subcluster data has yet been entered
    on PhagesDB. But it may be empty because there is no subcluster
    designation yet for members of the Cluster."""
    try:
        if data_dict['psubcluster'] is None:
            subcluster = "none"
        else:
            subcluster = data_dict['psubcluster']['subcluster']
        eval_object = None

    except:
        subcluster = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Subcluster data from PhagesDB.")

    return (subcluster, eval_object)





#TODO unit test
def parse_phagesdb_host(data_dict):
    """Retrieve Host from PhagesDB."""
    try:
        host = data_dict['isolation_host']['genus']
        eval_object = None

    except:
        host = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Host data from PhagesDB.")

    return (host, eval_object)




#TODO unit test
def parse_phagesdb_accession(data_dict):
    """Retrieve Accession from PhagesDB."""
    try:
        accession = data_dict['genbank_accession']
        eval_object = None

    except:
        accession = ""
        eval_object = Eval.construct_error( \
            "Unable to retrieve Accession data from PhagesDB.")

    return (accession, eval_object)





#TODO unit test
def parse_phagesdb_sequence(data_dict):
    """Retrieve genome sequence from PhagesDB."""

    # Note: old code first tested  "online_data_dict['fasta_file'] is not None".
    # May need to re-implement this logic.
    try:
        fastafile_url = data_dict['fasta_file']
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
    phage_name, result1 = parse_phagesdb_phage_name(data_dict):
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





































###
