"""Misc. functions to interact with NCBI databases."""
from Bio import Entrez, SeqIO

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
def get_records(accessions, db="nucleotide", rettype="gb", retmode="text"):
    """Retrieve records from NCBI from a list of active accessions."""
    fetch_query = ",".join(accessions)
    fetch_handle = Entrez.efetch(db=db, id=fetch_query, rettype=rettype, retmode=retmode)
    fetch_records = SeqIO.parse(fetch_handle, "genbank")
    for record in fetch_records:
        retrieved_records.append(record)
    fetch_handle.close()
    return retrieved_records









###
