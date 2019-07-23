"""Misc. functions useful for comparing and processing genomes."""


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC





def match_genome_by_id(bundle_obj, genome_dict, key1, key2 = None):
    """Match genome object to another genome object using id.
    The 'key1' parameter provides the type of genome stored in the Bundle
    genome dictionary to base the match from.
    The 'key2' parameter provides the type of genome to be stored in
    the Bundle genome dictionary."""

    try:
        ref_genome = bundle_obj.genome_dict[key1]
        if ref_genome.id in genome_dict.keys():
            matched_genome = genome_dict[ref_genome.id]

            if key2 is None:
                bundle_obj.genome_dict[matched_genome.type] = matched_genome
            else:
                bundle_obj.genome_dict[key2] = matched_genome
    except:
        pass




# TODO this is probably not needed anymore.
def match_genomes(list_of_bundle_objects, genome_dict, key1, key2 = None):
    """Match genome object to another genome object using id.
    The 'key1' parameter provides the type of genome stored in the Bundle
    genome dictionary to base the match from.
    The 'key2' parameter provides the type of genome to be stored in
    the Bundle genome dictionary."""

    index = 0
    while index < len(list_of_bundle_objects):

        bundle_obj = list_of_bundle_objects[index]

        try:
            ref_genome = bundle_obj.genome_dict[key1]
            if ref_genome.id in genome_dict.keys():
                matched_genome = genome_dict[ref_genome.id]

                if key2 is None:
                    bundle_obj.genome_dict[matched_genome.type] = matched_genome
                else:
                    bundle_obj.genome_dict[key2] = matched_genome

        except:
            pass

        index += 1




def create_fasta_seqrecord(header, sequence_string):
    """Create a fasta-formatted Biopython SeqRecord object."""

    seq = Seq(sequence_string, alphabet = IUPAC.unambiguous_dna)
    seqrecord = SeqRecord(seq, description = header)

    return seqrecord





















###
