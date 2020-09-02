"""Misc. functions to interact with NCBI databases."""

from Bio import Entrez, SeqIO

from pdm_utils.functions import basic

# GLOBAL VARIABLES
# ----------------------------------------------------------------------------
RETTYPE_MAPPINGS = {"gb": "gb", "tbl": "ft"}


# TODO unittest.
def set_entrez_credentials(tool=None, email=None, api_key=None):
    """Set BioPython Entrez credentials to improve speed and reliability.

    :param tool: Name of the software/tool being used.
    :type tool: str
    :param email: Email contact information for NCBI.
    :type email: str
    :param api_key: Unique NCBI-issued identifier to enhance retrieval speed.
    :type api_key: str
    """
    if tool is not None:
        Entrez.tool = tool
    if email is not None:
        Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key


# TODO unittest.
def run_esearch(db="", term="", usehistory=""):
    """Search for valid records in NCBI.

    Uses NCBI esearch implemented through BioPython Entrez.

    :param db: Name of the database to search.
    :type db: str
    :param term: Search term.
    :type term: str
    :param usehistory: Indicates if prior searches should be used.
    :type usehistory: str
    :return: Results of the search for each valid record.
    :rtype: dict
    """
    search_handle = Entrez.esearch(db=db, term=term, usehistory=usehistory)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    return search_record


# TODO unittest.
def get_summaries(db="", query_key="", webenv=""):
    """Retrieve record summaries from NCBI for a list of accessions.

    Uses NCBI esummary implemented through BioPython Entrez.

    :param db: Name of the database to get summaries from.
    :type db: str
    :param query_key:
        Identifier for the search.
        This can be directly generated from run_esearch().
    :type query_key: str
    :param webenv: Identifier that can be directly generated from run_esearch()
    :type webenv: str
    :return: List of dictionaries, where each dictionary is a record summary.
    :rtype: list
    """
    summary_handle = Entrez.esummary(db=db, query_key=query_key, webenv=webenv)
    summary_records = Entrez.read(summary_handle)
    summary_handle.close()
    return summary_records


# TODO test.
def get_accessions_to_retrieve(summary_records):
    """Extract accessions from summary records.

    :param summary_records:
        List of dictionaries, where each dictionary is a record summary.
    :type summary_records: list
    :return: List of accessions.
    :rtype: list
    """
    accessions = []
    for doc_sum in summary_records:
        doc_sum_accession = doc_sum["Caption"]
        accessions.append(doc_sum_accession)
    return accessions


def get_data_handle(accession_list, db="nucleotide", rettype="gb",
                    retmode="text"):
    fetch_query = ",".join(accession_list)
    fetch_handle = Entrez.efetch(db=db, id=fetch_query, rettype=rettype,
                                 retmode=retmode)

    return fetch_handle


# TODO Owen test
def get_verified_data_handle(acc_id_dict, ncbi_cred_dict={}, batch_size=200,
                             file_type="gb"):
    """Retrieve genomes from GenBank.

    output_folder = Path to where files will be saved.
    acc_id_dict = Dictionary where key = Accession and value = List[PhageIDs]
    """

    # More setup variables if NCBI updates are desired. NCBI Bookshelf resource
    # "The E-utilities In-Depth: Parameters, Syntax and More", by Dr. Eric
    # Sayers, recommends that a single request not contain more than about 200
    # UIDS so we will use that as our batch size, and all Entrez requests must
    # include the user's email address and tool name.
    set_entrez_credentials(
        tool=ncbi_cred_dict.get("tool"),
        email=ncbi_cred_dict.get("email"),
        api_key=ncbi_cred_dict.get("api_key"))

    # Use esearch to verify the accessions are valid and efetch to retrieve
    # the record
    # Create batches of accessions
    unique_accession_list = list(acc_id_dict.keys())

    # Add [ACCN] field to each accession number
    appended_accessions = \
        [accession + "[ACCN]" for accession in unique_accession_list]

    # When retrieving in batch sizes, first create the list of values
    # indicating which indices of the unique_accession_list should be used
    # to create each batch.
    # For instace, if there are five accessions, batch size of two produces
    # indices = 0,2,4

    chunked_accessions = basic.partition_list(appended_accessions, batch_size)
    for chunk in chunked_accessions:
        delimiter = " | "
        esearch_term = delimiter.join(chunk)

        # Use esearch for each accession
        search_record = run_esearch(db="nucleotide", term=esearch_term,
                                    usehistory="y")
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]
        summary_records = get_summaries(db="nucleotide",
                                        query_key=search_query_key,
                                        webenv=search_webenv)

        accessions_to_retrieve = get_accessions_to_retrieve(summary_records)
        if len(accessions_to_retrieve) > 0:
            fetch_handle = get_data_handle(accessions_to_retrieve,
                                           rettype=RETTYPE_MAPPINGS[file_type])
            return fetch_handle
        else:
            return None


# TODO unittest.
def get_records(accession_list, db="nucleotide", rettype="gb", retmode="text"):
    """Retrieve records from NCBI from a list of active accessions.

    Uses NCBI efetch implemented through BioPython Entrez.

    :param accession_list: List of NCBI accessions.
    :type accession_list: list
    :param db: Name of the database to get summaries from (e.g. 'nucleotide').
    :type db: str
    :param rettype: Type of record to retrieve (e.g. 'gb').
    :type rettype: str
    :param retmode: Format of data to retrieve (e.g. 'text').
    :type retmode: str
    :return: List of BioPython SeqRecords generated from GenBank records.
    :rtype: list
    """
    retrieved_records = []
    fetch_query = ",".join(accession_list)
    fetch_handle = Entrez.efetch(db=db, id=fetch_query, rettype=rettype,
                                 retmode=retmode)
    fetch_records = SeqIO.parse(fetch_handle, "genbank")
    for record in fetch_records:
        retrieved_records.append(record)
    fetch_handle.close()
    return retrieved_records
