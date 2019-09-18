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












# # TODO move this to another module and construct this function.
# def custom_runmode():
#
#     if run_mode_custom_total > 0:
#
#     	print "\n\nOne or more import tickets indicate custom QC parameters."
#     	print "The following options will be applied to all tickets requiring a custom QC:"
#
#     	#1. Use the file's basename as the PhageID instead of the phage name in the file
#     	print "\n\n\n\nNormally, the PhageID is determined from the Organism field in the record."
#     	run_mode_custom_dict['use_basename'] = select_option(
#     		"\nInstead, do you want to use the file name as the PhageID? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#2. Create GeneIDs from locus tags or PhageID_GeneNumber concatenation
#     	#Many times non-Hatfull authored genomes do not have consistent or specific locus tags, which complicates
#     	#assigning GeneIDs. This option provides a locus tag override and will assign GeneIDs
#     	#by joining PhageID and CDS number.
#     	print "\n\n\n\nNormally, the GeneIDs are assigned by using the locus tags."
#     	run_mode_custom_dict['custom_gene_id'] = select_option(
#     		"\nInstead, do you want to create GeneIDs by combining the PhageID and gene number? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#3. Ensure GeneIDs have phage name spelled correctly?
#     	#New SEA-PHAGES annotated genomes should be check for spelling, but maybe not
#     	#for other types of genomes.
#     	print "\n\n\n\nNormally, the GeneIDs are required to contain the phage name without typos."
#     	run_mode_custom_dict['ignore_gene_id_typo'] = select_option(
#     		"\nInstead, do you want to allow missing or mispelled phage names in the GeneID? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#4. Gene descriptions are set to import table field regardless
#     	#This should be run for new SEA-PHAGES annotated genomes, but may want to
#     	#be skipped if importing many genomes from NCBI
#     	print "\n\n\n\nNormally, the gene descriptions are verified to be present in the import table qualifier."
#     	run_mode_custom_dict['ignore_description_field_check'] = select_option(
#     		"\nInstead, do you want to use the import table qualifier without verification? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#5. Replacing final with final warning
#     	#Once a genome gets in NCBI, it is expected that a Final status genome is
#     	#replaced with another Final status genome, so it could get annoying to
#     	#keep getting the warning, so it can be turned off.
#     	print "\n\n\n\nNormally, a warning is indicated if a Final status genome is being replaced."
#     	run_mode_custom_dict['ignore_replace_warning'] = select_option(
#     		"\nInstead, do you want to silence the status warnings? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#6. tRNA QC
#     	#Many genomes from NCBI, including SEA-PHAGES, may not have consistently
#     	#annotated tRNAs. So the tRNA QC can be skipped.
#     	print "\n\n\n\nNormally, tRNAs are checked only in new manually annotated genomes."
#     	run_mode_custom_dict['ignore_trna_check'] = select_option(
#     		"\nInstead, do you want to ignore the tRNA quality checks? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#7. Retain locus tags
#     	#The locus tag field in the database should only reflect 'official' locus tags from
#     	#bona fide Genbank records, and not simply Genbank-formatted files.
#     	#Locus tags from Pecaan auto-annotated and SMART team manually annotated
#     	#genomes should not be retained.
#     	print "\n\n\n\nNormally, CDS locus tags are retained only for bona fide Genbank records."
#     	run_mode_custom_dict['ignore_locus_tag_import'] = select_option(
#     		"\nInstead, do you want to ignore all locus tags? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#8. Ignore phage name typos in the header
#     	#The phage name can be found in several header fields. This should be
#     	#checked for new manual annotations. Since NCBI doesn't like to change these
#     	#types of typos, parsing NCBI records should skip this step.
#     	print "\n\n\n\nNormally, the phage name is verified only in new manual annotations."
#     	run_mode_custom_dict['ignore_phage_name_typos'] = select_option(
#     		"\nInstead, do you want to ignore any phage name typos in the header? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#9. Ignore host name typos in the header
#     	#The host name can be found in several header fields. This should be
#     	#checked for new manual annotations and in SEA-PHAGES NCBI records.
#     	print "\n\n\n\nNormally, the host name is verified in new manual annotations or SEA-PHAGES NCBI records."
#     	run_mode_custom_dict['ignore_host_typos'] = select_option(
#     		"\nInstead, do you want to ignore any host name typos in the header? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#10. Ignore generic author in the author list
#     	#Sometimes the generic author 'Lastname' or 'Firstname' gets added in DNA Master
#     	#This should be checked for new manual annotations only.
#     	print "\n\n\n\nNormally, the author list is checked to ensure the generic author Lastname,Firstname is absent in new manual annotations."
#     	run_mode_custom_dict['ignore_generic_author'] = select_option(
#     		"\nDo you want to ignore any generic authors in the author list? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#     	#11. Ignore gene description checks
#     	#The gene description may contain errors.
#     	#This should be checked for new manual annotations only.
#     	print "\n\n\n\nNormally, gene descriptions are checked in new manual annotations."
#     	run_mode_custom_dict['ignore_description_check'] = select_option(
#     		"\nDo you want to ignore gene description checks? (yes or no) ", \
#     		set(['yes','y','no','n']))
#
#
#
#
#
#     	#Now add the customized dictionary of options to the dictionary of all run modes
#     	run_mode_options_dict['custom'] = run_mode_custom_dict














###
