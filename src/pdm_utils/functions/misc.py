"""Misc. functions useful for comparing and processing genomes."""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import csv
from pdm_utils.constants import constants

def match_genome_by_id(bndl, genome_dict, key1, key2=None):
    """Match genome object to another genome object using id.

    :param bndl:
        A bundle object containing at least two pdm_utils genome objects.
    :type bndl: pdm_utils bundle
    :param key1:
        Indicates the type of genome stored in the bundle
        genome dictionary to base the match from.
    :type key1: str
    :param key2:
        Indicates the type of genome to be stored in
        the bundle genome dictionary.
    :type key2: str
    """
    try:
        ref_genome = bndl.genome_dict[key1]
        if ref_genome.id in genome_dict.keys():
            matched_genome = genome_dict[ref_genome.id]
            if key2 is None:
                bndl.genome_dict[matched_genome.type] = matched_genome
            else:
                bndl.genome_dict[key2] = matched_genome
    except:
        pass


# TODO this is probably not needed anymore.
def match_genomes(bndl_list, genome_dict, key1, key2=None):
    """Match genome object to another genome object using id.

    :param bndl_list: List of bundle objects.
    :type bndl_list: list
    :param genome_dict: A dictionary of pdm_utils genome objects.
    :type genome_dict: dict
    :param key1:
        Indicates the type of genome stored in the Bundle
        genome dictionary to base the match from.
    :type key1: str
    :param key2:
        Indicates the type of genome to be stored in
        the bundle genome dictionary.
    :type key2: str
    """
    index = 0
    while index < len(bndl_list):
        bndl = bndl_list[index]
        try:
            ref_genome = bndl.genome_dict[key1]
            if ref_genome.id in genome_dict.keys():
                matched_genome = genome_dict[ref_genome.id]
                if key2 is None:
                    bndl.genome_dict[matched_genome.type] = matched_genome
                else:
                    bndl.genome_dict[key2] = matched_genome
        except:
            pass
        index += 1


def create_fasta_seqrecord(header, sequence_string):
    """Create a fasta-formatted Biopython SeqRecord object.

    :param header: Description of the sequence.
    :type header: str
    :param sequence_string: Nucleotide sequence.
    :type sequence_string: str
    :returns: Biopython SeqRecord containing the nucleotide sequence.
    :rtype: SeqRecord
    """
    seq = Seq(sequence_string, alphabet = IUPAC.unambiguous_dna)
    seqrecord = SeqRecord(seq, description = header)
    return seqrecord
