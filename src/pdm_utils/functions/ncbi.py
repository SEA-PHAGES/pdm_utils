"""Misc. functions to interact with NCBI databases."""
from Bio import Entrez, SeqIO

# TODO unittest.
def set_entrez_credentials(tool=None, email=None, api_key=None):
    """Set Entrez credentials to improve speed and reliability."""
    if tool is not None:
        Entrez.tool = tool
    if email is not None:
        Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key


# TODO unittest.
def run_esearch(db="", term="", usehistory=""):
    """Run esearch."""
    search_handle = Entrez.esearch(db=db, term=term, usehistory=usehistory)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    return search_record


# TODO unittest.
def get_summaries(db="", query_key="", webenv=""):
    """Retrieve summaries from NCBI for a list of accessions using esummary."""
    summary_handle = Entrez.esummary(db=db, query_key=query_key, webenv=webenv)
    summary_records = Entrez.read(summary_handle)
    summary_handle.close()
    return summary_records


# TODO unittest.
def get_records(accession_list, db="nucleotide", rettype="gb", retmode="text"):
    """Retrieve records from NCBI from a list of active accessions."""
    retrieved_records = []
    fetch_query = ",".join(accession_list)
    fetch_handle = Entrez.efetch(db=db, id=fetch_query, rettype=rettype, retmode=retmode)
    fetch_records = SeqIO.parse(fetch_handle, "genbank")
    for record in fetch_records:
        retrieved_records.append(record)
    fetch_handle.close()
    return retrieved_records









###
