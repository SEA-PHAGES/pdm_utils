"""Pipeline to compare data between MySQL, PhagesDB, and GenBank databases."""

# TODO this object-oriented pipeline is not fully integrated into
# the pdm_utils package. The majority of the code is redundant and can be
# dramatically reduced once integrated. Steps that are redundant:
# 1. Many class definitions (CDS features, Genomes)
# 2. ORM functions to map data from MySQL, PhagesDB, and GenBank to the classes
# 3. GenBank data retrieval

# Note this script compares and matches data from GenBank data and MySQL data.
# As a result, there are many similarly named variables.
# Variables are prefixed to indicate database:
# GenBank =  "gbk", "g"
# MySQL = "mysql", "m"
# PhagesDB = "pdb", "p"

import argparse
import csv
from datetime import date
import os
import pathlib
import re
import sys
import time

from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import phagesdb


DEFAULT_OUTPUT_FOLDER = os.getcwd()

#Create output directories
CURRENT_DATE = date.today().strftime("%Y%m%d")
WORKING_FOLDER = f"{CURRENT_DATE}_compare"

GENBANK_OUTPUT = "gbk_records"
MYSQL_OUTPUT = "mysql_records"
PHAGESDB_OUTPUT = "phagesdb_records"

ERROR_FOLDER = "errors"
RECORD_FOLDER = "records"
RESULTS_FOLDER = "results"

COMPARE_SETTINGS = "compare_settings.csv"
DUPLICATE_MYSQL_NAMES = "mysql_name_duplicates.csv"
DUPLICATE_MYSQL_ACC = "mysql_accession_duplicates.csv"
DUPLICATE_PDB_NAMES = "phagesdb_name_duplicates.csv"
FAILED_ACC_RETRIEVE = "mysql_genomes_unmatched_to_genbank.csv"
UNMATCHED_GENOMES = "mysql_genomes_unmatched_to_pdb.csv"

GNM_SUMMARY_OUTPUT = "summary_genome.csv"
CDS_SUMMARY_OUTPUT = "summary_cds.csv"

GENOME_OUTPUT = "genome.csv"
GENE_OUTPUT = "cds.csv"

PHAGE_QUERY = ("SELECT PhageID, Name, HostGenus, Sequence, Length, "
               "Status, Cluster, Accession, RetrieveRecord, "
               "DateLastModified, AnnotationAuthor FROM phage")
GENE_QUERY = ("SELECT PhageID, GeneID, Name, Start, Stop, Orientation, "
              "Translation, Notes from gene")
VERSION_QUERY = "SELECT Version FROM version"


# Define classes




# Base genome class.
# Merged UnannotatedGenome, AnnotatedGenome, PhagesdbGenome,
# MysqlGenome, and GbkGenome classes.
class Genome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields
        self.type = "" # Identifier to describes source of this genome
                       # (e.g. mysql, phagesdb, genbank, etc.)
        self.id = ""
        self.name = ""
        self.host_genus = ""
        self.length = 0
        self.seq = Seq("", IUPAC.ambiguous_dna)
        self.accession = ""
        self.annotation_status = "" # Final, Draft, Unknown version of genome data
        self.cluster = ""
        self.subcluster = ""
        self.retrieve_record = ""
        self.date = ""
        self.annotation_author = "" # 1 (Hatfull), 0 (GenBank)

        # GenBank data
        self.description = ""
        self.source = ""
        self.organism = ""
        self.authors = ""

        # Computed data fields
        self.cds_features = []
        self._cds_features_tally = 0
        self._cds_descriptions_tally = 0
        self._cds_functions_tally = 0
        self._cds_products_tally = 0
        self._cds_notes_tally = 0


    # Define all attribute setters:
    def set_accession(self, value):
        if value is None or value.strip() == "":
            self.accession = ""
        else:
            value = value.strip()
            self.accession = value.split(".")[0]

    def set_cds_features(self, value):
        self.cds_features = value # Should be a list
        self._cds_features_tally = len(self.cds_features)

    def set_cluster(self, value):
        if value is None:
            self.cluster = "Singleton"
        elif value == "UNK":
            self.cluster = ""
        else:
            self.cluster = value





# Genome attributes
setattr(Genome, "sequence", "") # TODO string but seq is Seq
setattr(Genome, "record_name", "")
setattr(Genome, "record_id", "")
setattr(Genome, "source_feature_organism", "")
setattr(Genome, "source_feature_host", "")
setattr(Genome, "source_feature_lab_host", "")
setattr(Genome, "_missing_locus_tags_tally", 0)
setattr(Genome, "_locus_tag_typos_tally", 0)
setattr(Genome, "_description_field_error_tally", 0)
setattr(Genome, "_status_accession_error", False)
setattr(Genome, "_status_description_error", False)
setattr(Genome, "_genes_with_errors_tally", 0)
setattr(Genome, "_nucleotide_errors", False)
setattr(Genome, "_cds_features_with_translation_error_tally", 0)
setattr(Genome, "_cds_features_boundary_error_tally", 0)


# Genome setters
def set_sequence(self, value):
    self.sequence = value.upper()
    self.length = len(self.sequence)
setattr(Genome, "set_sequence", set_sequence)

# Genome error checks

def compute_nucleotide_errors(self, dna_alphabet):
    nucleotide_set = set(self.sequence)
    nucleotide_error_set = nucleotide_set - dna_alphabet
    if len(nucleotide_error_set) > 0:
        self._nucleotide_errors = True
setattr(Genome, "compute_nucleotide_errors", compute_nucleotide_errors)

def compute_cds_feature_errors(self):
    for cds_feature in self.cds_features:
        if cds_feature._amino_acid_errors:
            self._cds_features_with_translation_error_tally += 1
        if cds_feature._boundary_error:
            self._cds_features_boundary_error_tally += 1
setattr(Genome, "compute_cds_feature_errors", compute_cds_feature_errors)

def check_status_accession(self):

    # Be sure to first set the accession attribute before the
    # annotation_status attribute, else this will throw an error.
    # Now that the AnnotationAuthor field contains authorship data, the
    # 'unknown' annotation annotation_status now reflects
    # an 'unknown' annotation (in regards to if it was auto-annotated
    # or manually annotated).
    # So for the annotation_status-accession error, if the
    # annotation_status is 'unknown', there is no reason to assume
    # whether there should be an accession or not. Only for 'final'
    # (manually annotated) genomes should there be an accession.
    if self.annotation_status == "final" and self.accession == "":
        self._status_accession_error = True
setattr(Genome, "check_status_accession", check_status_accession)

def compute_status_description_error(self):
    # Iterate through all CDS features, see if they have descriptions,
    # then compare to the annotation_status.
    for feature in self.cds_features:
        if feature.raw_description != "":
            self._cds_descriptions_tally += 1
    if self.annotation_status == "draft" and self._cds_descriptions_tally > 0:
        self._status_description_error = True
    elif self.annotation_status == "final" and self._cds_descriptions_tally == 0:
        self._status_description_error = True
    else:
        pass
setattr(Genome, "compute_status_description_error", compute_status_description_error)

#Even though this method iterates through the CDS features
# like the compute_status_description_error does,
#it has to be kept separate, since you need to wait to run
# this method after all genome and gene matching is completed.
def compute_genes_with_errors_tally(self):
    for feature in self.cds_features:
        # Need to first compute the number of errors per gene
        feature.compute_total_cds_errors()
        if feature._total_errors > 0:
            self._genes_with_errors_tally += 1
setattr(Genome, "compute_genes_with_errors_tally", compute_genes_with_errors_tally)

def compute_gbk_cds_feature_errors(self):
    for cds_feature in self.cds_features:

        # counting descriptions should skip if it is blank
        # or "hypothetical protein".
        if cds_feature.product != "":
            self._cds_products_tally += 1

        if cds_feature.function != "":
            self._cds_functions_tally += 1

        if cds_feature.note != "":
            self._cds_notes_tally += 1

        if cds_feature._locus_tag_missing:
            self._missing_locus_tags_tally += 1
        else:

            search_name = basic.edit_suffix(self.name, "remove").lower()
            pattern4 = re.compile(search_name)

            search_result = pattern4.search(cds_feature.locus_tag.lower())

            if search_result == None:
                self._locus_tag_typos_tally += 1
                cds_feature.set_locus_tag_typo() # Sets to True

        if cds_feature._description_field_error:
            self._description_field_error_tally += 1
setattr(Genome, "compute_gbk_cds_feature_errors", compute_gbk_cds_feature_errors)








class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        # Datafields from MySQL database:
        self.type = "" # Feature type: CDS, GenomeBoundary,or tRNA
        self.id = "" # Gene ID comprised of PhageID and Gene name
        self.name = ""
        self.start = "" # Position of left boundary, 0-indexed
        self.stop = "" # Position of right boundary, 0-indexed
        self.orientation = "" # 'forward', 'reverse', or 'NA'
        self.translation = ""
        self.translation_length = ""
        self.genome_id = ""
        self.raw_description = ""
        self.description = "" # non-generic gene descriptions
        self.locus_tag = ""


        # GenBank data
        self.gene = ""
        self.raw_product = ""
        self.raw_function = ""
        self.raw_note = ""
        self.product = ""
        self.function = ""
        self.note = ""

    # Define all attribute setters:
    def set_strand(self, value):
        self.orientation = basic.reformat_strand(value, "fr_long")

    def set_translation(self, value):
        self.translation = value.upper()
        self.translation_length = len(self.translation)

    def set_genome_id(self, value):
        self.genome_id = value

    def set_notes(self, value1, value2):
        self.raw_description = value1
        self.description = value2

    def set_product_description(self, value1, value2):
        self.raw_product = value1
        self.product = value2

    def set_function_description(self, value1, value2):
        self.raw_function = value1
        self.function = value2

    def set_note_description(self, value1, value2):
        self.raw_note = value1
        self.note = value2









###Cds
setattr(CdsFeature, "_search_genome_id", "")
setattr(CdsFeature, "_start_end_strand_id", "")
setattr(CdsFeature, "_end_strand_id", "")
setattr(CdsFeature, "_locus_tag_missing", False)
setattr(CdsFeature, "_locus_tag_typo", False)
setattr(CdsFeature, "_description_field_error", False)
setattr(CdsFeature, "_amino_acid_errors", False)
setattr(CdsFeature, "_boundary_error", False)
setattr(CdsFeature, "_unmatched_error", False)
setattr(CdsFeature, "_total_errors", 0)

def set_start_end_strand_id(self):
    # Create a tuple of feature location data.
    # For start and end of feature, it doesn't matter whether
    # the feature is complex with a translational frameshift or not.
    # Retrieving the "start" and "end" attributes return the
    # very beginning and end of the feature,
    # disregarding the inner "join" coordinates.
    self._start_end_strand_id = (str(self.start),str(self.stop),self.orientation)

    # Since this id matched genes with different start sites,
    # the orientation impacts whether the left or right boundary is used
    if self.orientation == "forward":
        self._end_strand_id = (str(self.stop),self.orientation)
    elif self.orientation == "reverse":
        self._end_strand_id = (str(self.start),self.orientation)
    else:
        pass
setattr(CdsFeature, "set_start_end_strand_id", set_start_end_strand_id)

def set_search_genome_id(self):
    self._search_genome_id = basic.edit_suffix(self.genome_id, "remove").lower()
setattr(CdsFeature, "set_search_genome_id", set_search_genome_id)


def check_locus_tag(self):
    if self.locus_tag == "":
        self._locus_tag_missing = True
setattr(CdsFeature, "check_locus_tag", check_locus_tag)

def compute_amino_acid_errors(self, protein_alphabet):
    amino_acid_set = set(self.translation)
    amino_acid_error_set = amino_acid_set - protein_alphabet
    if len(amino_acid_error_set) > 0:
        self._amino_acid_errors = True
setattr(CdsFeature, "compute_amino_acid_errors", compute_amino_acid_errors)

def compute_boundary_error(self):
    # Check if start and end coordinates are fuzzy
    if not (str(self.start).isdigit() and str(self.stop).isdigit()):
        self._boundary_error = True
setattr(CdsFeature, "compute_boundary_error", compute_boundary_error)

def set_locus_tag_typo(self):
    self._locus_tag_typo = True
setattr(CdsFeature, "set_locus_tag_typo", set_locus_tag_typo)

def compute_description_error(self):
    # If the product description is empty or generic,
    # and the function or note descriptions are not, there is an error.
    if (self.product == "" and self.function != "" or self.note != ""):
        self._description_field_error = True
setattr(CdsFeature, "compute_description_error", compute_description_error)

def compute_total_cds_errors(self):
    if self._amino_acid_errors:
        self._total_errors += 1
    if self._boundary_error:
        self._total_errors += 1
    if self._description_field_error:
        self._total_errors += 1
    if self._locus_tag_missing:
        self._total_errors += 1
    if self._locus_tag_typo:
        self._total_errors += 1
    if self._unmatched_error:
        self._total_errors += 1
setattr(CdsFeature, "compute_total_cds_errors", compute_total_cds_errors)









class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.m_genome = ""
        self.p_genome = ""
        self.g_genome = ""

        # MySQL and GenBank matched data comparison results
        self._m_g_seq_mismatch = False
        self._m_g_seq_length_mismatch = False
        self._g_header_fields_name_mismatch = False
        self._g_host_mismatch = False

        # List of MatchedCdsFeature objects
        self._m_g_perfect_matched_ftrs = []
        self._m_g_imperfect_matched_ftrs = []

        # List of CdsFeature objects
        self._m_ftrs_unmatched_in_g = []
        self._g_ftrs_unmatched_in_m = []

        self._m_g_perfect_matched_ftrs_tally = 0
        self._m_g_imperfect_matched_ftrs_tally = 0
        self._m_ftrs_unmatched_in_g_tally = 0
        self._g_ftrs_unmatched_in_m_tally = 0
        self._m_g_different_descriptions_tally = 0
        self._m_g_different_translations_tally = 0
        self._m_g_author_error = False

        # MySQL and PhagesDB matched data comparison results
        self._m_p_seq_mismatch = False
        self._m_p_seq_length_mismatch = False
        self._m_p_host_mismatch = False
        self._m_p_accession_mismatch = False
        self._m_p_cluster_mismatch = False


        # PhagesDB and GenBank matched data comparison results
        self._p_g_seq_mismatch = False
        self._p_g_seq_length_mismatch = False

        # Total errors summary
        self._contains_errors = False
        self._total_number_genes_with_errors = 0


    # Define all attribute setters:

    def compare_mysql_gbk_genomes(self):

        # verify that there is a MySQL and GenBank genome in
        # the matched genome object.
        m_gnm = self.m_genome
        g_gnm = self.g_genome
        m_cds_lst = m_gnm.cds_features

        if isinstance(m_gnm, Genome) and isinstance(g_gnm, Genome):

            if m_gnm.sequence != g_gnm.sequence:
                self._m_g_seq_mismatch = True
            if m_gnm.length != g_gnm.length:
                self._m_g_seq_length_mismatch = True

            # Compare phage names
            pattern1 = re.compile("^" + m_gnm.name + "$")
            pattern2 = re.compile("^" + m_gnm.name)

            if (basic.find_expression(pattern2,g_gnm.description.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.source.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.organism.split(" ")) == 0 or
                    basic.find_expression(pattern1,g_gnm.source_feature_organism.split(" ")) == 0):

                self._g_header_fields_name_mismatch = True

            # Compare host_genus data
            search_host = m_gnm.host_genus
            if search_host == "Mycobacterium":
                search_host = search_host[:-3]
            pattern3 = re.compile("^" + search_host)

            if ((basic.find_expression(pattern3,g_gnm.description.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.source.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.organism.split(" ")) == 0 or
                    basic.find_expression(pattern3,g_gnm.source_feature_organism.split(" ")) == 0) or
                    (g_gnm.source_feature_host != "" and basic.find_expression(pattern3,g_gnm.source_feature_host.split(" ")) == 0) or
                    (g_gnm.source_feature_lab_host != "" and basic.find_expression(pattern3,g_gnm.source_feature_lab_host.split(" ")) == 0)):

                self._g_host_mismatch = True

            # Check author list for errors
            # For genomes with AnnotationAuthor = 1 (Hatfull), Graham is expected
            # to be an author.
            # For genomes with AnnotationAuthor = 0 (non-Hatfull/GenBank), Graham
            # is NOT expected to be an author.
            pattern5 = re.compile("hatfull")
            search_result = pattern5.search(g_gnm.authors.lower())
            if m_gnm.annotation_author == 1 and search_result == None:
                self._m_g_author_error = True
            elif m_gnm.annotation_author == 0 and search_result != None:
                self._m_g_author_error = True
            else:
                # Any other combination of MySQL and GenBank author can be skipped
                pass


            # Compare CDS features

            # First find all unique start-end-orientation CDS identifiers
            # for MySQL and GenBank genomes.
            m_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            m_start_end_strand_duplicate_id_set = set()

            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id not in m_start_end_strand_id_set:
                    m_start_end_strand_id_set.add(cds_ftr._start_end_strand_id)
                else:
                    m_start_end_strand_duplicate_id_set.add(cds_ftr._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            m_start_end_strand_id_set = m_start_end_strand_id_set - m_start_end_strand_duplicate_id_set


            g_cds_list = g_gnm.cds_features
            g_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            g_start_end_strand_duplicate_id_set = set()
            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id not in g_start_end_strand_id_set:
                    g_start_end_strand_id_set.add(cds_ftr._start_end_strand_id)

                else:
                    g_start_end_strand_duplicate_id_set.add(cds_ftr._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            g_start_end_strand_id_set = g_start_end_strand_id_set - g_start_end_strand_duplicate_id_set

            # Create the perfect matched and unmatched sets
            m_unmatched_start_end_strand_id_set = m_start_end_strand_id_set - g_start_end_strand_id_set
            g_unmatched_start_end_strand_id_set = g_start_end_strand_id_set - m_start_end_strand_id_set
            perfect_matched_cds_id_set = m_start_end_strand_id_set & g_start_end_strand_id_set


            # From the unmatched sets, created second round of end-orientation id sets
            # MySQL end_strand data
            m_end_strand_id_set = set()

            # All end_strand ids that are not unique
            m_end_strand_duplicate_id_set = set()
            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id in m_unmatched_start_end_strand_id_set:
                    if cds_ftr._end_strand_id not in m_end_strand_id_set:
                        m_end_strand_id_set.add(cds_ftr._end_strand_id)
                    else:
                        m_end_strand_duplicate_id_set.add(cds_ftr._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            m_end_strand_id_set = m_end_strand_id_set - m_end_strand_duplicate_id_set


            g_end_strand_id_set = set()

            # All end_strand ids that are not unique
            g_end_strand_duplicate_id_set = set()
            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id in g_unmatched_start_end_strand_id_set:
                    if cds_ftr._end_strand_id not in g_end_strand_id_set:
                        g_end_strand_id_set.add(cds_ftr._end_strand_id)
                    else:
                        g_end_strand_duplicate_id_set.add(cds_ftr._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            g_end_strand_id_set = g_end_strand_id_set - g_end_strand_duplicate_id_set


            # Create the imperfect matched set
            imperfect_matched_cds_id_set = m_end_strand_id_set & g_end_strand_id_set


            # Now go back through all CDS features and assign
            # to the appropriate dictionary or list.
            m_perfect_matched_cds_dict = {}
            m_imperfect_matched_cds_dict = {}
            m_unmatched_cds_list = []


            for cds_ftr in m_cds_lst:
                if cds_ftr._start_end_strand_id in perfect_matched_cds_id_set:
                    m_perfect_matched_cds_dict[cds_ftr._start_end_strand_id] = cds_ftr
                elif cds_ftr._end_strand_id in imperfect_matched_cds_id_set:
                    m_imperfect_matched_cds_dict[cds_ftr._end_strand_id] = cds_ftr
                else:
                    m_unmatched_cds_list.append(cds_ftr)


            g_perfect_matched_cds_dict = {}
            g_imperfect_matched_cds_dict = {}
            g_unmatched_cds_list = []

            for cds_ftr in g_cds_list:
                if cds_ftr._start_end_strand_id in perfect_matched_cds_id_set:
                    g_perfect_matched_cds_dict[cds_ftr._start_end_strand_id] = cds_ftr
                elif cds_ftr._end_strand_id in imperfect_matched_cds_id_set:
                    g_imperfect_matched_cds_dict[cds_ftr._end_strand_id] = cds_ftr
                else:
                    g_unmatched_cds_list.append(cds_ftr)


            # Create MatchedCdsFeatures objects
            # Compute matched gene errors and tallies
            # Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object._m_feature = m_perfect_matched_cds_dict[start_end_strand_tup]
                matched_cds_object._g_feature = g_perfect_matched_cds_dict[start_end_strand_tup]
                matched_cds_object.compare_mysql_gbk_cds_ftrs()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1



                if matched_cds_object._m_g_different_translations:
                    self._m_g_different_translations_tally += 1
                if matched_cds_object._m_g_different_descriptions:
                    self._m_g_different_descriptions_tally += 1
                self._m_g_perfect_matched_ftrs.append(matched_cds_object)


            # Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object._m_feature = m_imperfect_matched_cds_dict[end_strand_tup]
                matched_cds_object._g_feature = g_imperfect_matched_cds_dict[end_strand_tup]
                matched_cds_object.compare_mysql_gbk_cds_ftrs()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1

                if matched_cds_object._m_g_different_descriptions:
                    self._m_g_different_descriptions_tally += 1
                self._m_g_imperfect_matched_ftrs.append(matched_cds_object)


            # Compute unmatched error and gene total errors for
            # all unmatched features.
            for cds_ftr in m_unmatched_cds_list:
                cds_ftr._unmatched_error = True
                cds_ftr.compute_total_cds_errors()
                if cds_ftr._total_errors > 0:
                    self._total_number_genes_with_errors += 1

            for cds_ftr in g_unmatched_cds_list:
                cds_ftr._unmatched_error = True
                cds_ftr.compute_total_cds_errors()
                if cds_ftr._total_errors > 0:
                    self._total_number_genes_with_errors += 1


            # Set unmatched CDS lists
            self._m_ftrs_unmatched_in_g = m_unmatched_cds_list
            self._g_ftrs_unmatched_in_m = g_unmatched_cds_list

            # Now compute the number of features in each category
            self._m_g_perfect_matched_ftrs_tally = len(self._m_g_perfect_matched_ftrs)
            self._m_g_imperfect_matched_ftrs_tally = len(self._m_g_imperfect_matched_ftrs)
            self._m_ftrs_unmatched_in_g_tally = len(self._m_ftrs_unmatched_in_g)
            self._g_ftrs_unmatched_in_m_tally = len(self._g_ftrs_unmatched_in_m)


        # If there is no matching GenBank genome,
        # assign all MySQL genes to Unmatched.
        else:

            # Set unmatched CDS lists, but do NOT count them in unmatched tally.
            # Unmatched tally should reflect unmatched genes if there is
            # actually a matching GenBank genome.
            self._m_ftrs_unmatched_in_g = m_cds_lst

            # Now that all errors have been computed for each gene,
            # compute how many genes have errors
            # If there is no matching GenBank genome,
            # gene error tallies are only computed for the MySQL genome.
            m_gnm.compute_genes_with_errors_tally()
            self._total_number_genes_with_errors = m_gnm._genes_with_errors_tally


    def compare_mysql_phagesdb_genomes(self):

        # verify that there is a MySQL and PhagesDB genome
        # in the matched genome object.
        m_gnm = self.m_genome
        p_gnm = self.p_genome

        if isinstance(m_gnm, Genome) and isinstance(p_gnm, Genome):

            if m_gnm.sequence != p_gnm.sequence:
                self._m_p_seq_mismatch = True
            if m_gnm.length != p_gnm.length:
                self._m_p_seq_length_mismatch = True
            if m_gnm.accession != p_gnm.accession:
                self._m_p_accession_mismatch = True
            if m_gnm.host_genus != p_gnm.host_genus:
                self._m_p_host_mismatch = True
            if m_gnm.cluster != p_gnm.cluster:
                self._m_p_cluster_mismatch = True


    def compare_phagesdb_gbk_genomes(self):

        # verify that there is a PhagesDB and GenBank genome
        # in the matched genome object.
        p_gnm = self.p_genome
        g_gnm = self.g_genome

        if isinstance(p_gnm, Genome) and isinstance(g_gnm, Genome):
            if p_gnm.sequence != g_gnm.sequence:
                self._p_g_seq_mismatch = True
            if p_gnm.length != g_gnm.length:
                self._p_g_seq_length_mismatch = True



    def compute_total_genome_errors(self):
        m_gnm = self.m_genome
        g_gnm = self.g_genome
        p_gnm = self.p_genome

        if self._total_number_genes_with_errors > 0:
            self._contains_errors = True

        # MySQL genome
        if m_gnm._nucleotide_errors:
            self._contains_errors = True
        if m_gnm._status_description_error:
            self._contains_errors = True
        if m_gnm._status_accession_error:
            self._contains_errors = True

        # GenBank genome
        if isinstance(g_gnm, Genome):
            if g_gnm._nucleotide_errors:
                self._contains_errors = True

        # PhagesDB genome
        if isinstance(p_gnm, Genome):
            if p_gnm._nucleotide_errors:
                self._contains_errors = True

        # MySQL-GenBank
        if self._m_g_seq_mismatch:
            self._contains_errors = True
        if self._m_g_seq_length_mismatch:
            self._contains_errors = True
        if self._g_header_fields_name_mismatch:
            self._contains_errors = True
        if self._g_host_mismatch:
            self._contains_errors = True
        if self._m_g_imperfect_matched_ftrs_tally > 0:
            self._contains_errors = True
        if self._m_ftrs_unmatched_in_g_tally > 0:
            self._contains_errors = True
        if self._g_ftrs_unmatched_in_m_tally > 0:
            self._contains_errors = True
        if self._m_g_different_descriptions_tally > 0:
            self._contains_errors = True
        if self._m_g_different_translations_tally > 0:
            self._contains_errors = True
        if self._m_g_author_error:
            self._contains_errors = True

        # MySQL-PhagesDB
        if self._m_p_seq_mismatch:
            self._contains_errors = True
        if self._m_p_seq_length_mismatch:
            self._contains_errors = True
        if self._m_p_host_mismatch:
            self._contains_errors = True
        if self._m_p_accession_mismatch:
            self._contains_errors = True
        if self._m_p_cluster_mismatch:
            self._contains_errors = True

        # PhagesDB-GenBank
        if self._p_g_seq_mismatch:
            self._contains_errors = True
        if self._p_g_seq_length_mismatch:
            self._contains_errors = True



class MatchedCdsFeatures:

    # Initialize all attributes:
    def __init__(self):


        # Initialize all non-calculated attributes:
        self._m_feature = ""
        self._g_feature = ""

        # Matched data comparison results
        self._m_g_different_translations = False
        self._m_g_different_start_sites = False
        self._m_g_different_descriptions = False

        # Total errors summary
        self._total_errors = 0


    # Define all attribute setters:

    def compare_mysql_gbk_cds_ftrs(self):

        if self._m_feature.orientation == "forward":
            if str(self._m_feature.start) != str(self._g_feature.start):
                self._m_g_different_start_sites = True
        elif self._m_feature.orientation == "reverse":
            if str(self._m_feature.stop) != str(self._g_feature.stop):
                self._m_g_different_start_sites = True
        else:
            pass


        product_description_set = set()
        product_description_set.add(self._m_feature.description)
        product_description_set.add(self._g_feature.product)


        if len(product_description_set) != 1:
            self._m_g_different_descriptions = True

        if self._m_feature.translation != self._g_feature.translation:
            self._m_g_different_translations = True



        # Compute total errors
        # First add all matched feature errors
        if self._m_g_different_translations:
            self._total_errors += 1

        if self._m_g_different_start_sites:
            self._total_errors += 1

        if self._m_g_different_descriptions:
            self._total_errors += 1

        # Now add all errors from each individual feature
        # You first compute errors for each individual feature.
        # This step is performed here instead of in the mainline code
        # because you need to wait for the feature matching step
        # after the genome matching step.
        self._m_feature.compute_total_cds_errors()
        self._g_feature.compute_total_cds_errors()
        self._total_errors += self._m_feature._total_errors
        self._total_errors += self._g_feature._total_errors



class DatabaseSummary:

    # Initialize all attributes:
    def __init__(self,matched_genomes_list):

        # Initialize all non-calculated attributes:
        self._matched_genomes_list = matched_genomes_list

        # Initialize all calculated attributes:

        # MySQL data
        # General genome data
        self._m_g_update_flag_tally = 0

        # Genome data checks
        self._m_gnms_with_nucleotide_errors_tally = 0
        self._m_gnms_with_translation_errors_tally = 0
        self._m_gnms_with_boundary_errors_tally = 0
        self._m_gnms_with_status_accession_error_tally = 0
        self._m_gnms_with_status_description_error_tally = 0

        # PhagesDB data
        # Genome data checks
        self._p_gnms_with_nucleotide_errors_tally = 0

        # GenBank data
        # Genome data checks
        self._g_gnms_with_nucleotide_errors_tally = 0
        self._g_gnms_with_translation_errors_tally = 0
        self._g_gnms_with_boundary_errors_tally = 0
        self._g_gnms_with_missing_locus_tags_tally = 0
        self._g_gnms_with_locus_tag_typos_tally = 0
        self._g_gnms_with_description_field_errors_tally = 0

        # MySQL-PhagesDB checks
        self._m_p_seq_mismatch_tally = 0
        self._m_p_seq_length_mismatch_tally = 0
        self._m_p_cluster_mismatch_tally = 0
        self._m_p_accession_mismatch_tally = 0
        self._m_p_host_mismatch_tally = 0

        # MySQL-GenBank checks
        self._m_g_seq_mismatch_tally = 0
        self._m_g_seq_length_mismatch_tally = 0
        self._m_g_header_phage_mismatch_tally = 0
        self._m_g_header_host_mismatch_tally = 0
        self._m_g_gnms_with_imperfectly_matched_ftrs_tally = 0
        self._m_g_gnms_with_unmatched_m_ftrs_tally = 0
        self._m_g_gnms_with_unmatched_g_ftrs_tally = 0
        self._m_g_gnms_with_different_descriptions_tally = 0
        self._m_g_gnms_with_different_translations_tally = 0
        self._m_g_gnms_with_author_errors_tally = 0

        # PhagesDB-GenBank checks
        self._p_g_seq_mismatch_tally = 0
        self._p_g_seq_length_mismatch_tally = 0


        # MySQL feature
        # Gene data checks
        self._m_translation_errors_tally = 0
        self._m_boundary_errors_tally = 0


        # GenBank feature
        # Gene data checks
        self._g_translation_errors_tally = 0
        self._g_boundary_errors_tally = 0
        self._g_missing_locus_tags_tally = 0
        self._g_locus_tag_typos_tally = 0
        self._g_description_field_errors_tally = 0

        # MySQL-GenBank checks
        self._m_g_different_descriptions_tally = 0
        self._m_g_different_start_sites_tally = 0
        self._m_g_different_translation_tally = 0
        self._m_g_unmatched_m_ftrs_tally = 0
        self._m_g_unmatched_g_ftrs_tally = 0


        # Calculate summary metrics
        self._m_total_gnms_analyzed = 0
        self._m_gnms_unmatched_to_p_tally = 0
        self._m_gnms_unmatched_to_g_tally = 0
        self._total_genomes_with_errors = 0


    # Define setter functions
    def compute_summary(self):
        for matched_genomes in self._matched_genomes_list:

            self._m_total_gnms_analyzed += 1
            m_gnm = matched_genomes.m_genome
            p_gnm = matched_genomes.p_genome
            g_gnm = matched_genomes.g_genome

            if matched_genomes._contains_errors:
                self._total_genomes_with_errors += 1


            # MySQL data
            if isinstance(m_gnm, Genome):

                self._m_g_update_flag_tally += m_gnm.retrieve_record

                # Genome data checks
                if m_gnm._nucleotide_errors:
                    self._m_gnms_with_nucleotide_errors_tally += 1
                if m_gnm._cds_features_with_translation_error_tally > 0:
                    self._m_gnms_with_translation_errors_tally += 1
                    self._m_translation_errors_tally += m_gnm._cds_features_with_translation_error_tally
                if m_gnm._cds_features_boundary_error_tally > 0:
                    self._m_gnms_with_boundary_errors_tally += 1
                    self._m_boundary_errors_tally += m_gnm._cds_features_boundary_error_tally
                if m_gnm._status_accession_error:
                    self._m_gnms_with_status_accession_error_tally += 1
                if m_gnm._status_description_error:
                    self._m_gnms_with_status_description_error_tally += 1

            # PhagesDB data
            if isinstance(p_gnm, Genome):

                # Genome data checks
                if p_gnm._nucleotide_errors:
                    self._p_gnms_with_nucleotide_errors_tally += 1
            else:
                self._m_gnms_unmatched_to_p_tally += 1

            # GenBank data
            if isinstance(g_gnm, Genome):

                # Genome data checks
                if g_gnm._nucleotide_errors:
                    self._g_gnms_with_nucleotide_errors_tally += 1
                if g_gnm._cds_features_with_translation_error_tally > 0:
                    self._g_gnms_with_translation_errors_tally += 1
                    self._g_translation_errors_tally += g_gnm._cds_features_with_translation_error_tally
                if g_gnm._cds_features_boundary_error_tally > 0:
                    self._g_gnms_with_boundary_errors_tally += 1
                    self._g_boundary_errors_tally += g_gnm._cds_features_boundary_error_tally
                if g_gnm._missing_locus_tags_tally > 0:
                    self._g_gnms_with_missing_locus_tags_tally += 1
                    self._g_missing_locus_tags_tally += g_gnm._missing_locus_tags_tally
                if g_gnm._locus_tag_typos_tally > 0:
                    self._g_gnms_with_locus_tag_typos_tally += 1
                    self._g_locus_tag_typos_tally += g_gnm._locus_tag_typos_tally
                if g_gnm._description_field_error_tally > 0:
                    self._g_gnms_with_description_field_errors_tally += 1
                    self._g_description_field_errors_tally += g_gnm._description_field_error_tally

            else:
                self._m_gnms_unmatched_to_g_tally += 1

            # MySQL-PhagesDB checks
            if matched_genomes._m_p_seq_mismatch:
                self._m_p_seq_mismatch_tally += 1
            if matched_genomes._m_p_seq_length_mismatch:
                self._m_p_seq_length_mismatch_tally += 1
            if matched_genomes._m_p_cluster_mismatch:
                self._m_p_cluster_mismatch_tally += 1
            if matched_genomes._m_p_accession_mismatch:
                self._m_p_accession_mismatch_tally += 1
            if matched_genomes._m_p_host_mismatch:
                self._m_p_host_mismatch_tally += 1

            # MySQL-GenBank checks
            if matched_genomes._m_g_seq_mismatch:
                self._m_g_seq_mismatch_tally += 1
            if matched_genomes._m_g_seq_length_mismatch:
                self._m_g_seq_length_mismatch_tally += 1
            if matched_genomes._g_header_fields_name_mismatch:
                self._m_g_header_phage_mismatch_tally += 1
            if matched_genomes._g_host_mismatch:
                self._m_g_header_host_mismatch_tally += 1
            if matched_genomes._m_g_imperfect_matched_ftrs_tally > 0:
                self._m_g_gnms_with_imperfectly_matched_ftrs_tally += 1
                self._m_g_different_start_sites_tally += matched_genomes._m_g_imperfect_matched_ftrs_tally
            if matched_genomes._m_ftrs_unmatched_in_g_tally > 0:
                self._m_g_gnms_with_unmatched_m_ftrs_tally += 1
                self._m_g_unmatched_m_ftrs_tally += matched_genomes._m_ftrs_unmatched_in_g_tally
            if matched_genomes._g_ftrs_unmatched_in_m_tally > 0:
                self._m_g_gnms_with_unmatched_g_ftrs_tally += 1
                self._m_g_unmatched_g_ftrs_tally += matched_genomes._g_ftrs_unmatched_in_m_tally
            if matched_genomes._m_g_different_descriptions_tally > 0:
                self._m_g_gnms_with_different_descriptions_tally += 1
                self._m_g_different_descriptions_tally += matched_genomes._m_g_different_descriptions_tally
            if matched_genomes._m_g_different_translations_tally > 0:
                self._m_g_gnms_with_different_translations_tally += 1
                self._m_g_different_translation_tally += matched_genomes._m_g_different_translations_tally
            if matched_genomes._m_g_author_error:
                self._m_g_gnms_with_author_errors_tally += 1


            # PhagesDB-GenBank checks
            if matched_genomes._p_g_seq_mismatch:
                self._p_g_seq_mismatch_tally += 1
            if matched_genomes._p_g_seq_length_mismatch:
                self._p_g_seq_length_mismatch_tally += 1

# End of class definitions




def output_to_file(data_list, folder, filename):
    """Output list data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        for element in data_list:
            writer.writerow(element)

def prepare_unmatched_to_gbk_output(gnms):
    """Prepare list of MySQL unmatched to GenBank data to be saved to file."""
    l = []
    if len(gnms) > 0:
        l.append(["PhageID", "PhageName", "Author", "Status", "Accession"])
        for gnm in gnms:
            l.append([gnm.id,
                      gnm.name,
                      gnm.annotation_author,
                      gnm.annotation_status,
                      gnm.accession])
    return l

def prepare_unmatched_to_pdb_output(gnms):
    """Prepare list of MySQL unmatched to PhagesDB data to be saved to file."""
    l = []
    if len(gnms) > 0:
        l.append(["PhageID", "PhageName", "Author", "Status"])
        for gnm in gnms:
            l.append([gnm.id,
                      gnm.name,
                      gnm.annotation_author,
                      gnm.annotation_status])
    return l

def save_seqrecord(seqrecord, output_path, file_prefix, ext, seqrecord_ext, interactive):
    """Save record to file."""
    file_path = basic.make_new_file(output_path, file_prefix, ext, attempt=100)
    if file_path is not None:
        SeqIO.write(seqrecord, file_path, seqrecord_ext)
    else:
        print(f"Duplicated filenames for {file_prefix}. Unable to output data.")
        if interactive:
            input('Press ENTER to proceed')



def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for comparing databases."""

    COMPARE_HELP = ("Pipeline to compare MySQL, PhagesDB, and "
                    "GenBank databases for inconsistencies.")
    DATABASE_HELP = "Name of the MySQL database from which to compare data."
    OUTPUT_FOLDER_HELP = ("Path to the folder to store results.")
    NCBI_CRED_FILE_HELP = ("Path to the file containing NCBI credentials.")

    PHAGESDB_HELP = "Indicates that PhagesDB data should be compared."
    GENBANK_HELP = "Indicates that GenBank data should be compared."
    SAVE_RECORDS_HELP = \
        ("Indicates that records retrieved from external "
         "databases will be saved.")
    INTERNAL_AUTHORS_HELP = \
        "Indicates that genomes with internal authorship will be evaluated."
    EXTERNAL_AUTHORS_HELP = \
        "Indicates that genomes with external authorship will be evaluated."

    DRAFT_HELP = \
        "Indicates that genomes with 'draft' annotation_status will be evaluated."
    FINAL_HELP = \
        "Indicates that genomes with 'final' annotation_status will be evaluated."
    UNKNOWN_HELP = \
        "Indicates that genomes with 'unknown' annotation_status will be evaluated."
    INTERACTIVE_HELP = \
        "Indicates whether evaluation is paused when errors are encountered."

    parser = argparse.ArgumentParser(description=COMPARE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
                        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                        help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-p", "--phagesdb", action="store_true",
        default=False, help=PHAGESDB_HELP)
    parser.add_argument("-g", "--genbank", action="store_true",
        default=False, help=GENBANK_HELP)
    parser.add_argument("-c", "--ncbi_credentials_file", type=pathlib.Path,
        help=NCBI_CRED_FILE_HELP)
    parser.add_argument("-s", "--save_records", action="store_true",
        default=False, help=SAVE_RECORDS_HELP)
    parser.add_argument("-ia", "--internal_authors", action="store_true",
        default=False, help=INTERNAL_AUTHORS_HELP)
    parser.add_argument("-ea", "--external_authors", action="store_true",
        default=False, help=EXTERNAL_AUTHORS_HELP)
    parser.add_argument("-d", "--draft", action="store_true",
        default=False, help=DRAFT_HELP)
    parser.add_argument("-f", "--final", action="store_true",
        default=False, help=FINAL_HELP)
    parser.add_argument("-u", "--unknown", action="store_true",
        default=False, help=UNKNOWN_HELP)
    parser.add_argument("-i", "--interactive", action="store_true",
        default=False, help=INTERACTIVE_HELP)



    # Assumed command line arg structure:
    # python3 -m pdm_utils <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


def get_dbs(pdb, gbk):
    """Create set of databases to compare to MySQL."""
    dbs = set()
    if pdb == True:
        dbs.add("phagesdb")
    if gbk == True:
        dbs.add("genbank")
    return dbs


def get_authors(internal, external):
    """Create set of authorship to compare."""
    authorship = set()
    if internal == True:
        authorship.add(1)
    if external == True:
        authorship.add(0)
    return authorship

def get_status(draft, final, unknown):
    """Create set of annotation_status to compare."""
    status = set()
    if draft == True:
        status.add("draft")
    if final == True:
        status.add("final")
    if unknown == True:
        status.add("unknown")
    return status




def main(unparsed_args_list):
    """Run compare pipeline."""
    start_time = time.strftime("%x %X")
    args = parse_args(unparsed_args_list)
    database = args.database
    save_records = args.save_records
    ncbi_credentials_file = args.ncbi_credentials_file
    interactive = args.interactive

    # Setup output
    output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    working_dir = pathlib.Path(WORKING_FOLDER)
    working_path = basic.make_new_dir(output_folder, working_dir,
                                      attempt=50)
    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    results_path = pathlib.Path(working_path, RESULTS_FOLDER)
    error_path = pathlib.Path(working_path, ERROR_FOLDER)
    records_path = pathlib.Path(working_path, RECORD_FOLDER)
    results_path.mkdir()
    error_path.mkdir()
    records_path.mkdir()

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    alchemist = AlchemyHandler(database=database)
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the compare pipeline")

    # Gather sets of settings to determine which types of comparison to do.
    valid_dbs = get_dbs(args.phagesdb, args.genbank)
    valid_authors = get_authors(args.internal_authors, args.external_authors)
    valid_status = get_status(args.draft, args.final, args.unknown)

    # Record the summary of settings for the comparison.
    record_compare_settings(working_path, database, engine, valid_status,
                            valid_authors, valid_dbs)

    pmd_tup = process_mysql_data(working_path, engine, valid_status,
                                 valid_authors, interactive, save_records)
    mysql_genome_dict = pmd_tup[0]
    mysql_accessions = pmd_tup[1]
    mysql_acc_duplicates = pmd_tup[2]

    set_mysql_gnm_attr(mysql_genome_dict)

    #Now retrieve all PhagesDB data
    if "phagesdb" in valid_dbs:
        ppd_tup = process_phagesdb_data(working_path, interactive, save_records)
        pdb_genome_dict = ppd_tup[0]
        pdb_name_duplicates = ppd_tup[1]

    else:
        pdb_genome_dict = {}
        pdb_name_duplicates = set()

    #Retrieve and parse GenBank records if selected by user
    if "genbank" in valid_dbs:
        gbk_genome_dict = process_gbk_data(working_path, ncbi_credentials_file,
                                mysql_accessions, interactive, save_records)
    else:
        gbk_genome_dict = {} #Key = accession; #Value = genome data

    # Now that all GenBank and PhagesDB data is retrieved,
    # match up to MySQL genome data.
    match_tup = match_all_genomes(mysql_genome_dict, pdb_genome_dict,
                      gbk_genome_dict, pdb_name_duplicates,
                      mysql_acc_duplicates)
    matched_genomes_list = match_tup[0]

    #Output unmatched data to file
    record_unmatched_pdb_data(match_tup[1], working_path)
    record_unmatched_gbk_data(match_tup[2], working_path)

    # Run checks on each genome and all matched data.
    check_mysql_gnms(mysql_genome_dict)
    check_pdb_gnms(pdb_genome_dict)
    check_gbk_gnms(gbk_genome_dict)
    check_matched_gnms(matched_genomes_list)

    # Now that all individual matched_genome_objects have all
    # computed attributes, compute database summary.
    summarize_data(matched_genomes_list, working_path)

    end_time = time.strftime("%x %X")
    print(f"Start time: {start_time}")
    print(f"Stop time: {end_time}")

    return


def record_compare_settings(working_path, database, engine, valid_status,
                            valid_authors, valid_dbs):
    """Save user-selected settings."""

    # Get data from the MySQL database.
    version = str(mysqldb_basic.scalar(engine, VERSION_QUERY))

    # Format setting strings
    lst = [
        ["Comparison date:", CURRENT_DATE],
        ["MySQL database:", database],
        ["MySQL database version:", version],
        selected_authors_lst(valid_authors),
        ["Databases compared to MySQL:", ", ".join(valid_dbs)],
        ["Genomes with the following Annotation Status will be compared:",
             ", ".join(valid_status)]
         ]

    output_to_file(lst, pathlib.Path(working_path), COMPARE_SETTINGS)

def selected_authors_lst(lst1):
    lst1 = [str(i) for i in list(lst1)]
    lst2 =  ["Genomes with the following AnnotationAuthor will be compared:",
             ", ".join(lst1)]
    return lst2

def process_mysql_data(working_path, engine, valid_status, valid_authors,
                       interactive, save):
    """Retrieve and process MySQL data."""

    gnm_data_dict = mysqldb_basic.query_dict_list(engine, PHAGE_QUERY)
    gene_data_dict = mysqldb_basic.query_dict_list(engine, GENE_QUERY)

    tup = create_mysql_gnms(gnm_data_dict, valid_status, valid_authors)
    gnm_dict = tup[0]
    name_duplicates = tup[2]
    accessions = tup[3]
    accession_dupes = tup[4]

    if save == True:
        print("Saving MySQL genomes to file...")
        save_gnms_to_fasta(gnm_dict, pathlib.Path(working_path, RECORD_FOLDER),
                           MYSQL_OUTPUT, interactive)

    # PhagesDB relies on the phageName, and not the phageID.
    # But the MySQL database does not require phageName values to be unique.
    # Check if there are any phageName duplications.
    # If there are, they will not be able to be compared to PhagesDB data.
    # Output duplicate phage search names to file
    if len(name_duplicates) > 0:
        print("Warning: There are duplicate phage names in the MySQL database.")
        print("Some MySQL genomes will not be able to be matched to PhagesDB.")
        lst1 = [[i] for i in list(name_duplicates)]
        output_to_file(lst1, pathlib.Path(working_path, ERROR_FOLDER),
                       DUPLICATE_MYSQL_NAMES)
        if interactive:
            input('Press ENTER to proceed')

    # Accessions aren't required to be unique in the MySQL
    # database (but they should be), so there could be duplicates
    # Output duplicate names to file
    if len(accession_dupes) > 0:
        print("Warning: There are duplicate accessions in the MySQL database.")
        print("Some MySQL genomes will not be able to be "
              "matched to GenBank records.")
        lst2 = [[i] for i in list(accession_dupes)]
        output_to_file(lst2, pathlib.Path(working_path, ERROR_FOLDER),
                       DUPLICATE_MYSQL_ACC)
        if interactive:
            input('Press ENTER to proceed')

    mysql_genes = get_gene_obs_list(gene_data_dict)
    match_gnm_genes(gnm_dict, mysql_genes)

    return gnm_dict, accessions, accession_dupes


def create_mysql_gnms(list_of_dicts, valid_status, valid_authors):
    """Create genome objects from MySQL data."""

    print('\n\nPreparing genome data sets from the MySQL database...')
    gnm_dict = {}
    names = set()
    duplicate_names = set()
    accessions = set()
    duplicate_accessions = set()

    # Iterate through each MySQL genome and create a genome object
    for dict in list_of_dicts:

        # Check to see if the genome has a user-selected
        # annotation_status and authorship.
        if (dict["Status"] not in valid_status or
                dict["AnnotationAuthor"] not in valid_authors):
            continue
        else:
            gnm = create_mysql_gnm(dict)
            gnm_dict[dict["PhageID"]] = gnm

            # This keeps track of whether there are duplicate phage
            # names that will be used to match up to PhagesDB data.
            if gnm.id in names:
                duplicate_names.add(gnm.id)
            else:
                names.add(gnm.id)

            # This keeps track of whether there are duplicate
            # accession numbers that will be used to match up to GenBank data.
            if gnm.accession != "":
                if gnm.accession in accessions:
                    duplicate_accessions.add(gnm.accession)
                else:
                    accessions.add(gnm.accession)

    return (gnm_dict, names, duplicate_names, accessions, duplicate_accessions)


def create_mysql_gnm(dict):
    """Create genome object from MySQL data."""
    gnm = Genome()
    gnm.type = "mysql"
    gnm.id = dict["PhageID"]
    gnm.name = dict["Name"]
    gnm.host_genus = dict["HostGenus"]
    gnm.set_sequence(dict["Sequence"].decode("utf-8"))
    gnm.set_accession(dict["Accession"])
    gnm.annotation_status = dict["Status"]
    gnm.set_cluster(dict["Cluster"])
    gnm.retrieve_record = dict["RetrieveRecord"]
    gnm.date = dict["DateLastModified"]
    gnm.annotation_author = dict["AnnotationAuthor"]
    return gnm


def set_mysql_gnm_attr(gnm_dict):
    """Set compare-specific attributes not in pdm_utils Genome class."""
    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        for cds_ftr in gnm.cds_features:
            set_mysql_cds_attr(cds_ftr)



def save_gnms_to_fasta(gnm_dict, main_path, new_dir, interactive):
    """Save genome data to fasta file."""

    output_path = pathlib.Path(main_path, new_dir)
    output_path.mkdir()

    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        fasta = SeqRecord(Seq(gnm.sequence), id=gnm.id, description="")
        save_seqrecord(fasta, output_path, gnm.id, "fasta", "fasta", interactive)


def get_gene_obs_list(list_of_dicts):
    """Convert list of tuples of data to list of Cds objects."""

    print(f"\n\nPreparing {len(list_of_dicts)} gene data set(s) "
          "from the MySQL database...")
    lst = []
    for dict in list_of_dicts:
        cds_ftr = create_mysql_cds(dict)
        lst.append(cds_ftr)
    return lst


def create_mysql_cds(dict):
    """Create MySQL CDS feature object."""
    cds_ftr = CdsFeature()
    cds_ftr.set_genome_id(dict["PhageID"])
    cds_ftr.id = dict["GeneID"]
    cds_ftr.name = dict["Name"]
    cds_ftr.type = "CDS"
    cds_ftr.start = dict["Start"]
    cds_ftr.stop = dict["Stop"]
    cds_ftr.set_strand(dict["Orientation"])
    cds_ftr.set_translation(dict["Translation"])
    tup2 = basic.reformat_description(dict["Notes"].decode("utf-8"))
    cds_ftr.set_notes(tup2[0], tup2[1])
    return cds_ftr

def set_mysql_cds_attr(cds_ftr):
    """Set compare-specific MySQL Cds attributes not in pdm_utils Cds class."""
    cds_ftr.set_search_genome_id()
    cds_ftr.set_start_end_strand_id()


def check_mysql_cds(cds_ftr):
    """Check for errors in MySQL CDS feature."""
    cds_ftr.compute_amino_acid_errors(constants.PROTEIN_ALPHABET)
    cds_ftr.compute_boundary_error()


def match_gnm_genes(gnm_dict, gene_list):
    """Match MySQL gene and genome data."""
    print("Matching gene and genome data sets from the MySQL database...")
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]
        genes = []
        for gene in gene_list:
            if gene.genome_id == id:
                genes.append(gene)
        gnm.set_cds_features(genes)

def check_mysql_gnms(gnm_dict):
    """Check for errors in MySQL matched genome and CDS features."""
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]

        for cds_ftr in gnm.cds_features:
            check_mysql_cds(cds_ftr)

        gnm.check_status_accession()
        gnm.compute_nucleotide_errors(constants.DNA_ALPHABET)
        gnm.compute_cds_feature_errors()
        gnm.compute_status_description_error()


def process_phagesdb_data(working_path, interactive, save):
    """Retrieve data from PhagesDB and process results."""

    gbd_tup = get_phagesdb_data(interactive)
    gnm_dict = gbd_tup[0]
    names = gbd_tup[1]
    name_dupes = gbd_tup[2]

    if save == True:
        print("Saving PhagesDB genomes to file...")
        save_gnms_to_fasta(gnm_dict, pathlib.Path(working_path, RECORD_FOLDER),
                           PHAGESDB_OUTPUT, interactive)

    # PhagesDB phage names should be unique, but just make sure after
    # they are converted to a search name.
    if len(name_dupes) > 0:

        print("Warning: There are duplicate phage search names in PhagesDB.")
        print("Some PhagesDB genomes will not be able to be "
              "matched to genomes in MySQL.")
        lst = [[i] for i in list(name_dupes)]
        output_to_file(lst, pathlib.Path(working_path, ERROR_FOLDER),
                       DUPLICATE_PDB_NAMES)
        if interactive:
            input('Press ENTER to proceed')

    return gnm_dict, name_dupes


def get_phagesdb_data(interactive):
    """Retrieve data from PhagesDB."""
    print('\n\nRetrieving data from PhagesDB...')
    gnm_dict = {}
    names = set()
    dupe_names = set()
    data_list = phagesdb.get_phagesdb_data(constants.API_SEQUENCED)

    if len(data_list) > 0:
        for i in range(len(data_list)):
            gnm = create_pdb_gnm(data_list[i])
            pdb_id = gnm.id
            if pdb_id in names:
                dupe_names.add(pdb_id)
            else:
                names.add(pdb_id)
                gnm_dict[pdb_id] = gnm
    else:
        print("Error retrieving PhagesDB data.")
        if interactive:
            input('Press ENTER to proceed')
    return gnm_dict, names, dupe_names



def create_pdb_gnm(dict):
    """Parse PhagesDB data to genome object."""

    gnm = Genome()
    gnm.type = "phagesdb"

    #Name, Host, Accession
    gnm.id = dict["phage_name"]
    gnm.name = dict["phage_name"]
    gnm.host_genus = dict["isolation_host"]["genus"]
    gnm.set_accession(dict["genbank_accession"])

    #Cluster
    if dict["pcluster"] is not None:
        # Sometimes cluster information is not present.
        # In the PhagesDB database, it is recorded as NULL.
        # When phages data is downloaded from PhagesDB,
        # NULL cluster data is converted to "Unclustered".
        # In these cases, leave cluster as ''
        gnm.cluster = dict["pcluster"]["cluster"]

    #Subcluster
    if dict["psubcluster"] is not None:
        #A phage may be clustered but not subclustered.
        #In these cases, leave subcluster as ''
        gnm.subcluster = dict["psubcluster"]["subcluster"]

    #Check to see if there is a fasta file stored on PhagesDB for this phage
    if dict["fasta_file"] is not None:
        fasta_data = phagesdb.retrieve_url_data(dict["fasta_file"])
        header, seq = phagesdb.parse_fasta_data(fasta_data)
        gnm.set_sequence(seq)

    return gnm


def check_pdb_gnms(gnm_dict):
    """Check for errors in PhagesDB genome."""
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]
        if gnm.sequence != "":
            gnm.compute_nucleotide_errors(constants.DNA_ALPHABET)


def process_gbk_data(working_path, creds_file, accessions, interactive, save):
    """Retrieve and process GenBank data."""

    if save == True:
        output_path = pathlib.Path(working_path, RECORD_FOLDER, GENBANK_OUTPUT)
        output_path.mkdir()

    creds = ncbi.get_ncbi_creds(creds_file)
    records, retrieval_errors = get_genbank_data(creds, accessions)

    #Report the accessions that could not be retrieved.
    output_to_file(retrieval_errors, pathlib.Path(working_path, ERROR_FOLDER),
                   FAILED_ACC_RETRIEVE)

    gnm_dict = {}
    for record in records:
        gnm = create_gbk_gnm(record)
        set_gbk_gnm_attr(gnm, record)
        gnm_dict[gnm.accession] = gnm
        if save == True:
            save_gbk_genome(gnm, record, output_path, interactive)

    return gnm_dict


def get_genbank_data(ncbi_cred_dict, accession_set, batch_size=200):
    """Retrieve genomes from GenBank."""

    print("\n\nRetrieving records from GenBank")

    # Use esearch to verify the accessions are valid and
    # efetch to retrieve the record.
    ncbi.set_entrez_credentials(
        tool=ncbi_cred_dict["ncbi_tool"],
        email=ncbi_cred_dict["ncbi_email"],
        api_key=ncbi_cred_dict["ncbi_api_key"])


    #Create batches of accessions
    acc_list = list(accession_set)

    #Add [ACCN] field to each accession number
    for index in range(len(acc_list)):
        acc_list[index] = acc_list[index] + "[ACCN]"

    #Keep track of specific records
    retrieved_records = []
    retrieval_errors = []

    batch_indices = basic.create_indices(acc_list, batch_size)
    for indices in batch_indices:
        start = indices[0]
        stop = indices[1]

        current_batch_size = stop - start
        delimiter = " | "
        esearch_term = delimiter.join(acc_list[start:stop])

        #Use esearch for each accession
        search_handle = Entrez.esearch(db="nucleotide", term=esearch_term,
                                       usehistory="y")
        search_record = Entrez.read(search_handle)
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]

        # Keep track of the accessions that failed to be located in GenBank
        if search_count < current_batch_size:
            search_acc_fail = search_record["ErrorList"]["PhraseNotFound"]

            # Each element in this list is formatted "accession[ACCN]"
            for element in search_acc_fail:
                retrieval_errors.append(element[:-6])

        # Now retrieve all records using efetch
        fetch_handle = Entrez.efetch(db="nucleotide",
                                    rettype="gb",
                                    retmode="text",
                                    retstart=0,
                                    retmax=search_count,
                                    webenv=search_webenv,
                                    query_key=search_query_key)
        fetch_records = SeqIO.parse(fetch_handle,"genbank")

        for record in fetch_records:
            retrieved_records.append(record)

        search_handle.close()
        fetch_handle.close()
    return retrieved_records, retrieval_errors


def create_gbk_gnm(retrieved_record):
    """Parse GenBank data to genome object."""

    gnm = Genome()
    gnm.type = "genbank"

    try:
        # There may be a list of accessions associated with this file.
        # Discard the version suffix if it is present in the
        # Accession field (it might not be present).
        accession = retrieved_record.annotations["accessions"][0]
        accession = accession.split(".")[0]
        gnm.accession = accession
    except:
        gnm.accession = ""

    try:
        gnm.description = retrieved_record.description
    except:
        gnm.description = ""

    try:
        gnm.source = retrieved_record.annotations["source"]
    except:
        gnm.source = ""

    try:
        organism = retrieved_record.annotations["organism"]
        gnm.organism = organism

        # Only truncate organism name for the 'phage name' field
        if organism.split(" ")[-1].lower() == "unclassified.":
            gnm.name = organism.split(" ")[-2]
            gnm.id = organism.split(" ")[-2]
        else:
            gnm.name = organism.split(" ")[-1]
            gnm.id = organism.split(" ")[-1]
    except:
        gnm.organism = ""

    try:
        # The retrieved authors can be stored in multiple Reference elements
        record_references = retrieved_record.annotations["references"]
        record_references_author_list = []
        for reference in record_references:
            record_references_author_list.append(reference.authors)
        record_author_string = ";".join(record_references_author_list)
        gnm.authors = record_author_string
    except:
        gnm.authors = ""

    # Nucleotide sequence
    gnm.set_sequence(retrieved_record.seq)

    # Iterate through all features
    gbk_cds_ftrs = []
    for feature in retrieved_record.features:
        if feature.type == "CDS":
            gene_object = create_gbk_cds(feature)
            gbk_cds_ftrs.append(gene_object)

    gnm.set_cds_features(gbk_cds_ftrs)
    return gnm



def set_gbk_gnm_attr(gnm, retrieved_record):
    """Set compare-specific attributes not in pdm_utils Genome class."""
    gnm.record_name = ""
    gnm.record_id = ""

    source_feature_list = []
    for feature in retrieved_record.features:
        # Retrieve the Source Feature info
        if feature.type == "source":
            source_feature_list.append(feature)

    if len(source_feature_list) == 1:
        try:
            src_org = str(source_feature_list[0].qualifiers["organism"][0])
            gnm.source_feature_organism = src_org
        except:
            pass
        try:
            src_host = str(source_feature_list[0].qualifiers["host"][0])
            gnm.source_feature_host = src_host
        except:
            pass
        try:
            src_lab_host = str(source_feature_list[0].qualifiers["lab_host"][0])
            gnm.source_feature_lab_host = src_lab_host
        except:
            pass

    for cds_ftr in gnm.cds_features:
        set_gbk_cds_attr(cds_ftr)



def check_gbk_gnms(gnm_dict):
    """Check for errors in GenBank genome and CDS features."""
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]

        for cds_ftr in gnm.cds_features:
            check_gbk_cds(cds_ftr)
        gnm.compute_nucleotide_errors(constants.DNA_ALPHABET)
        gnm.compute_cds_feature_errors()
        gnm.compute_gbk_cds_feature_errors()


def create_gbk_cds(feature):
    """Parse data from GenBank CDS feature."""

    gene_object = CdsFeature()

    # Feature type
    gene_object.type = "CDS"

    # Locus tag
    try:
        gene_object.locus_tag = feature.qualifiers["locus_tag"][0]
    except:
        gene_object.locus_tag = ""


    # Orientation
    if feature.strand == 1:
        gene_object.set_strand("forward")
    elif feature.strand == -1:
        gene_object.set_strand("reverse")
    # ssRNA phages
    elif feature.strand is None:
        gene_object.set_strand("NA")

    # Gene boundary coordinates
    # Compound features are tricky to parse.
    if str(feature.location)[:4] == "join":

        # Skip this compound feature if it is comprised of more
        # than two features (too tricky to parse).
        if len(feature.location.parts) <= 2:
            # Retrieve compound feature positions based on orientation
            if feature.strand == 1:
                gene_object.start = str(feature.location.parts[0].start)
                gene_object.stop = str(feature.location.parts[1].end)
            elif feature.strand == -1:
                gene_object.start = str(feature.location.parts[1].start)
                gene_object.stop = str(feature.location.parts[0].end)
            # If strand is None, not sure how to parse it
            else:
                pass
    else:
        gene_object.start = str(feature.location.start)
        gene_object.stop = str(feature.location.end)

    # Translation
    try:
        gene_object.set_translation(feature.qualifiers["translation"][0])
    except:
        pass

    # Gene function, note, and product descriptions
    try:
        tup1 = basic.reformat_description(feature.qualifiers["product"][0])
        gene_object.set_product_description(tup1[0],tup1[1])
    except:
        pass
    try:
        tup2 = basic.reformat_description(feature.qualifiers["function"][0])
        gene_object.set_function_description(tup2[0],tup2[1])
    except:
        pass

    try:
        tup3 = basic.reformat_description(feature.qualifiers["note"][0])
        gene_object.set_note_description(tup3[0],tup3[1])
    except:
        pass

    # Gene number
    try:
        gene_object.gene = feature.qualifiers["gene"][0]
    except:
        pass

    return gene_object


def set_gbk_cds_attr(cds_ftr):
    """Set compare-specific Cds attributes not in pdm_utils Cds class."""
    cds_ftr.set_start_end_strand_id()


def check_gbk_cds(cds_ftr):
    """Check for errors in GenBank CDS feature."""
    cds_ftr.check_locus_tag()
    cds_ftr.compute_amino_acid_errors(constants.PROTEIN_ALPHABET)
    cds_ftr.compute_boundary_error()
    cds_ftr.compute_description_error()



def save_gbk_genome(gnm, record, output_path, interactive):
    """Save GenBank record to file."""
    file_prefix = gnm.id + "__" + gnm.accession
    save_seqrecord(record, output_path, file_prefix, "gb", "genbank", interactive)



def match_all_genomes(mysql_gnms, pdb_gnms, gbk_gnms, pdb_name_duplicates,
                      mysql_acc_duplicates):
    """Match MySQL, PhagesDB, and GenBank genomes."""

    print("Matching PhagesDB and GenBank genomes to MySQL genomes...")
    matched_gnms = []
    mysql_unmatched_to_pdb_gnms = []
    mysql_unmatched_to_gbk_gnms = []
    match_count = 1
    total_count = len(mysql_gnms.keys())
    # Iterate through each phage in the MySQL database
    for id in mysql_gnms.keys():

        print("Matching genome %s of %s" %(match_count,total_count))
        mysql_gnm = mysql_gnms[id]
        matched_objects = MatchedGenomes()
        matched_objects.m_genome = mysql_gnm

        # Match up PhagesDB genome
        # First try to match up the phageID, and if that doesn't work,
        # try to match up the phageName
        if mysql_gnm.id in pdb_gnms.keys():
            pdb_genome = pdb_gnms[mysql_gnm.id]
            # Make sure the pdb_genome doesn't have a search name
            # that was duplicated
            if pdb_genome.id in pdb_name_duplicates:
                pdb_genome = ""

        else:
            pdb_genome = ""
            mysql_unmatched_to_pdb_gnms.append(mysql_gnm)

        matched_objects.p_genome = pdb_genome

        # Now match up GenBank genome
        acc = mysql_gnm.accession
        if acc != "" and acc not in mysql_acc_duplicates:
                # Retrieval of record may have failed
                try:
                    g_gnm = gbk_gnms[acc]
                except:
                    g_gnm = ""
                    mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        else:
            g_gnm = ""
            mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        matched_objects.g_genome = g_gnm
        matched_gnms.append(matched_objects)
        match_count += 1

    return matched_gnms, mysql_unmatched_to_pdb_gnms, mysql_unmatched_to_gbk_gnms

def record_unmatched_gbk_data(gnms, working_path):
    """Save unmatched MySQL data not matched to GenBank."""
    lst1 = prepare_unmatched_to_gbk_output(gnms)
    output_to_file(lst1, pathlib.Path(working_path, ERROR_FOLDER),
                   FAILED_ACC_RETRIEVE)

def record_unmatched_pdb_data(gnms, working_path):
    """Save unmatched MySQL data not matched to PhagesDB."""
    lst2 = prepare_unmatched_to_pdb_output(gnms)
    output_to_file(lst2, pathlib.Path(working_path, ERROR_FOLDER),
                   UNMATCHED_GENOMES)

def check_matched_gnms(matched_gnms):
    """Compare all matched data."""

    print("Comparing matched genomes and identifying inconsistencies...")
    count = 1
    total = len(matched_gnms)
    for matched_gnm in matched_gnms:
        print(f"Comparing matched genome set {count} of {total}")
        matched_gnm.compare_mysql_gbk_genomes()
        matched_gnm.compare_mysql_phagesdb_genomes()
        matched_gnm.compare_phagesdb_gbk_genomes()
        matched_gnm.compute_total_genome_errors()
        count += 1





def summarize_data(matched_genomes_list, working_path):
    """Create summary of data and save."""
    summary_object = DatabaseSummary(matched_genomes_list)
    summary_object.compute_summary()
    output_all_data(working_path, summary_object)



def output_all_data(output_path, summary_object):
    """Output all analysis results."""
    print("Outputting results to file...")

    cds_file_path = pathlib.Path(output_path, RESULTS_FOLDER, GENE_OUTPUT)
    cds_fh = cds_file_path.open("w")
    cds_writer = csv.writer(cds_fh)
    cds_headers = create_gene_headers()
    cds_writer.writerow(cds_headers)

    #Create genome summary output file
    lst1 = create_genome_summary_data(summary_object)
    output_to_file(lst1, pathlib.Path(output_path, RESULTS_FOLDER),
                   GNM_SUMMARY_OUTPUT)

    #Create CDS summary output file
    lst3 = create_cds_summary_data(summary_object)
    output_to_file(lst3, pathlib.Path(output_path, RESULTS_FOLDER),
                   CDS_SUMMARY_OUTPUT)

    # Now iterate through matched objects.
    # All MySQL genomes are stored in a MatchedGenomes object,
    # even if there are no PhagesDB or GenBank matches.
    # All but a few MySQL phages should be matched to PhagesDB
    # Many MySQL phages should be matched to GenBank
    genomes_data = []
    for matched_genomes in summary_object._matched_genomes_list:

        gnm_data = create_genome_data(matched_genomes)
        genomes_data.append(gnm_data)

        #Once all matched genome data has been outputted,
        # iterate through all matched gene data
        m_gnm = matched_genomes.m_genome
        all_ftrs = get_all_features(matched_genomes)
        for mixed_ftr in all_ftrs:
            ftr_data = create_feature_data(m_gnm, mixed_ftr)
            cds_writer.writerow(ftr_data)

    lst2 = [create_genome_headers()]
    lst2.extend(genomes_data)
    output_to_file(lst2, pathlib.Path(output_path, RESULTS_FOLDER),
                   GENOME_OUTPUT)
    cds_fh.close()
    return


def get_all_features(matched_genomes):
    """Create list of all features."""

    perfect_match = matched_genomes._m_g_perfect_matched_ftrs
    imperfectly_match = matched_genomes._m_g_imperfect_matched_ftrs
    mysql_unmatch = matched_genomes._m_ftrs_unmatched_in_g
    gbk_unmatch = matched_genomes._g_ftrs_unmatched_in_m

    lst = []
    lst.extend(perfect_match)
    lst.extend(imperfectly_match)
    lst.extend(mysql_unmatch)
    lst.extend(gbk_unmatch)

    return lst


def create_gene_headers():
    """Create list of column headers."""
    headers = [
        # Gene summary
        "mysql_phage_name",
        "total_errors",

        # MySQL
        # General gene data
        "mysql_phage_id",
        "mysql_phage_id",
        "mysql_type_id",
        "mysql_gene_id",
        "mysql_gene_name",
        "mysql_left_boundary",
        "mysql_right_boundary",
        "mysql_strand",
        "mysql_translation",
        "mysql_translation_length",
        "mysql_gene_notes",

        # Gene data checks
        "mysql_translation_error",
        "mysql_gene_coords_error",

        # GenBank
        # General gene data
        "gbk_locus_tag",
        "gbk_gene_number",
        "gbk_type_id",
        "gbk_left_boundary",
        "gbk_right_boundary",
        "gbk_strand",
        "gbk_translation",
        "gbk_translation_length",
        "gbk_product_description",
        "gbk_function_description",
        "gbk_note_description",

        # Gene data checks
        "gbk_translation_error",
        "gbk_gene_coords_error",
        "gbk_missing_locus_tag",
        "gbk_locus_tag_typo",
        "gbk_description_field_error",

        # MySQL-GenBank checks
        "mysql_gbk_unmatched_error",
        "mysql_gbk_description_error",
        "mysql_gbk_start_coordinate_error",
        "mysql_gbk_translation_error"
        ]
    return headers


def create_genome_headers():
    """Create list of column headers."""
    headers = [
        "mysql_phage_id",
        "contains_errors",

        # MySQL
        # General genome data
        "mysql_phage_name",
        "mysql_phage_id",
        "mysql_phage_name",
        "mysql_status",
        "mysql_cluster",
        "mysql_host",
        "mysql_accession",
        "mysql_dna_seq_length",
        "mysql_gene_tally",
        "mysql_description_tally",
        "mysql_gbk_update_flag",
        "mysql_date_last_modified",
        "mysql_annotation_author",

        # Genome data checks
        "mysql_dna_seq_error",
        "mysql_gene_translation_error_tally",
        "mysql_gene_coords_error_tally",
        "mysql_status_description_error",
        "mysql_status_accession_error",

        # PhagesDB
        # General genome data
        "pdb_phage_name",
        "pdb_phage_name",
        "pdb_cluster",
        "pdb_subcluster",
        "pdb_host",
        "pdb_accession",
        "pdb_dna_seq_length",

        # Genome data checks
        "pdb_dna_seq_error",

        # GenBank
        # General genome data
        "gbk_phage_name",
        "gbk_phage_name",
        "gbk_record_id",
        "gbk_record_name",
        "gbk_record_accession",
        "gbk_record_definition",
        "gbk_record_source",
        "gbk_record_organism",
        "gbk_source_feature_organism",
        "gbk_source_feature_host",
        "gbk_source_feature_lab_host",
        "gbk_authors",
        "gbk_dna_seq_length",
        "gbk_gene_tally",

        # Genome data checks
        "gbk_dna_seq_error",
        "gbk_gene_translation_error_tally",
        "gbk_gene_coords_error_tally",
        "gbk_gene_product_tally",
        "gbk_gene_function_tally",
        "gbk_gene_note_tally",
        "gbk_missing_locus_tag_tally",
        "gbk_locus_tag_typo_tally",
        "gbk_description_field_error_tally",

        # MySQL-PhagesDB
        "mysql_pdb_dna_seq_error",
        "mysql_pdb_dna_seq_length_error",
        "mysql_pdb_cluster_error",
        "mysql_pdb_accession_error",
        "mysql_pdb_host_error",

        # MySQL-GenBank
        "mysql_gbk_dna_seq_error",
        "mysql_gbk_dna_seq_length_error",
        "mysql_gbk_record_header_name_error",
        "mysql_gbk_record_header_host_error",

        # Author error is dependent on MySQL genome annotation author and
        # GenBank list of authors, so this metric should be reported with
        # the other mysql_gbk error tallies.
        "mysql_gbk_author_error",
        "mysql_gbk_perfectly_matched_gene_tally",
        "mysql_gbk_imperfectly_matched_gene_tally",
        "mysql_gbk_unmatched_mysql_gene_tally",
        "mysql_gbk_unmatched_gbk_gene_tally",
        "mysql_gbk_gene_description_error_tally",
        "mysql_gbk_perfectly_matched_gene_translation_error_tally",

        # Number of genes with errors is computed slightly differently
        # depending on whethere there are matching MySQL and GenBank genomes.
        # Therefore,this metric should be reported with the other mysql_gbk
        #  error tallies even if there is no matching GenBank genome.
        "mysql_gbk_genes_with_errors_tally",

        # PhagesDB-GenBank
        "pdb_gbk_dna_seq_error",
        "pdb_gbk_dna_seq_length_error"
        ]
    return headers


def create_genome_summary_fields():
    """Create genome summary row ids."""
    lst = [
        # Column header
        "Database comparison metric",

        # Database summaries
        "mysql_total_genomes_analyzed",
        "mysql_genomes_unmatched_to_pdb_tally",
        "mysql_genomes_unmatched_to_gbk_tally",
        "total_genomes_with_errors",


        # MySQL data
        # General genome data
        "mysql_gbk_update_flag_tally",

        # Genome data checks
        "mysql_genomes_with_nucleotide_errors_tally",
        "mysql_genomes_with_translation_errors_tally",
        "mysql_genomes_with_boundary_errors_tally",
        "mysql_genomes_with_status_accession_error_tally",
        "mysql_genomes_with_status_description_error_tally",

        # PhagesDB data
        # Genome data checks
        "pdb_genomes_with_nucleotide_errors_tally",

        # GenBank data
        # Genome data checks
        "gbk_genomes_with_description_field_errors_tally",
        "gbk_genomes_with_nucleotide_errors_tally",
        "gbk_genomes_with_translation_errors_tally",
        "gbk_genomes_with_boundary_errors_tally",
        "gbk_genomes_with_missing_locus_tags_tally",
        "gbk_genomes_with_locus_tag_typos_tally",

        # MySQL-PhagesDB checks
        "mysql_pdb_sequence_mismatch_tally",
        "mysql_pdb_sequence_length_mismatch_tally",
        "mysql_pdb_cluster_mismatch_tally",
        "mysql_pdb_accession_mismatch_tally",
        "mysql_pdb_host_mismatch_tally",

        # MySQL-GenBank checks
        "mysql_gbk_sequence_mismatch_tally",
        "mysql_gbk_sequence_length_mismatch_tally",
        "mysql_gbk_record_header_phage_mismatch_tally",
        "mysql_gbk_record_header_host_mismatch_tally",
        "mysql_gbk_genomes_with_author_errors_tally",
        "mysql_gbk_genomes_with_imperfectly_matched_features_tally",
        "mysql_gbk_genomes_with_unmatched_mysql_features_tally",
        "mysql_gbk_genomes_with_unmatched_gbk_features_tally",
        "mysql_gbk_genomes_with_different_descriptions_tally",
        "mysql_gbk_genomes_with_different_translations_tally",

        # PhagesDB-GenBank checks
        "pdb_gbk_sequence_mismatch_tally",
        "pdb_gbk_sequence_length_mismatch_tally"
        ]
    return lst


def create_cds_summary_fields():
    """Create CDS summary row ids."""
    lst = [
        # Column header
        "Database comparison metric",

        # MySQL feature
        # Gene data checks
        "mysql__translation_errors_tally",
        "mysql__boundary_errors_tally",

        # GenBank feature
        # Gene data checks
        "gbk_translation_errors_tally",
        "gbk_boundary_errors_tally",
        "gbk_missing_locus_tags_tally",
        "gbk_locus_tag_typos_tally",
        "gbk_description_field_errors_tally",

        # MySQL-GenBank checks
        "mysql_gbk_different_descriptions_tally",
        "mysql_gbk_different_start_sites_tally",
        "mysql_gbk_different_translation_tally",
        "mysql_gbk_unmatched_mysql_features_tally",
        "mysql_gbk_unmatched_gbk_features_tally"
        ]
    return lst


def create_genome_summary_data(summary_object):
    """Create summary of all genome results."""

    lst1 = [
        # Column header
        "tally",

        # First output database summary data
        summary_object._m_total_gnms_analyzed,
        summary_object._m_gnms_unmatched_to_p_tally,
        summary_object._m_gnms_unmatched_to_g_tally,
        summary_object._total_genomes_with_errors,

        # MySQL data
        # General genome data
        summary_object._m_g_update_flag_tally,

        # Genome data checks
        summary_object._m_gnms_with_nucleotide_errors_tally,
        summary_object._m_gnms_with_translation_errors_tally,
        summary_object._m_gnms_with_boundary_errors_tally,
        summary_object._m_gnms_with_status_accession_error_tally,
        summary_object._m_gnms_with_status_description_error_tally,

        # PhagesDB data
        # Genome data checks
        summary_object._p_gnms_with_nucleotide_errors_tally,

        # GenBank data
        # Genome data checks
        summary_object._g_gnms_with_description_field_errors_tally,
        summary_object._g_gnms_with_nucleotide_errors_tally,
        summary_object._g_gnms_with_translation_errors_tally,
        summary_object._g_gnms_with_boundary_errors_tally,
        summary_object._g_gnms_with_missing_locus_tags_tally,
        summary_object._g_gnms_with_locus_tag_typos_tally,


        # MySQL-PhagesDB checks
        summary_object._m_p_seq_mismatch_tally,
        summary_object._m_p_seq_length_mismatch_tally,
        summary_object._m_p_cluster_mismatch_tally,
        summary_object._m_p_accession_mismatch_tally,
        summary_object._m_p_host_mismatch_tally,

        # MySQL-GenBank checks
        summary_object._m_g_seq_mismatch_tally,
        summary_object._m_g_seq_length_mismatch_tally,
        summary_object._m_g_header_phage_mismatch_tally,
        summary_object._m_g_header_host_mismatch_tally,
        summary_object._m_g_gnms_with_author_errors_tally,
        summary_object._m_g_gnms_with_imperfectly_matched_ftrs_tally,
        summary_object._m_g_gnms_with_unmatched_m_ftrs_tally,
        summary_object._m_g_gnms_with_unmatched_g_ftrs_tally,
        summary_object._m_g_gnms_with_different_descriptions_tally,
        summary_object._m_g_gnms_with_different_translations_tally,

        # PhagesDB-GenBank checks
        summary_object._p_g_seq_mismatch_tally,
        summary_object._p_g_seq_length_mismatch_tally
        ]

    lst2 = create_genome_summary_fields()
    lst3 = []
    for i in range(len(lst2)):
        lst3.append([lst2[i],lst1[i]])
    return lst3


def create_cds_summary_data(summary_object):
    """Create summary of all CDS results."""
    lst1 = [
        # Column header
        "tally",

        # MySQL feature
        # Gene data checks
        summary_object._m_translation_errors_tally,
        summary_object._m_boundary_errors_tally,

        # GenBank feature
        # Gene data checks
        summary_object._g_translation_errors_tally,
        summary_object._g_boundary_errors_tally,
        summary_object._g_missing_locus_tags_tally,
        summary_object._g_locus_tag_typos_tally,
        summary_object._g_description_field_errors_tally,

        # MySQL-GenBank checks
        summary_object._m_g_different_descriptions_tally,
        summary_object._m_g_different_start_sites_tally,
        summary_object._m_g_different_translation_tally,
        summary_object._m_g_unmatched_m_ftrs_tally,
        summary_object._m_g_unmatched_g_ftrs_tally
        ]

    lst2 = create_cds_summary_fields()
    lst3 = []
    for i in range(len(lst2)):
        lst3.append([lst2[i],lst1[i]])
    return lst3


def create_feature_data(mysql_genome, mixed_ftr):
    """Create feature data to output."""

    lst = []

    # Gene summaries
    # Add MySQL genome search name to each gene row regardless of
    # the type of CDS data (matched or unmatched).
    lst.append(mysql_genome.id) # matched MySQL search name
    lst.append(mixed_ftr._total_errors) # total # of errors for this gene

    if isinstance(mixed_ftr,MatchedCdsFeatures):
        mysql_ftr = mixed_ftr._m_feature
        gbk_ftr = mixed_ftr._g_feature
    else:
        if isinstance(mixed_ftr, CdsFeature):
            mysql_ftr = mixed_ftr
            gbk_ftr = ""
        elif isinstance(mixed_ftr, CdsFeature):
            mysql_ftr = ""
            gbk_ftr = mixed_ftr
        else:
            mysql_ftr = ""
            gbk_ftr = ""

    # MySQL feature
    if isinstance(mysql_ftr, CdsFeature):

        # General gene data
        lst.append(mysql_ftr.genome_id)
        lst.append(mysql_ftr._search_genome_id)
        lst.append(mysql_ftr.type)
        lst.append(mysql_ftr.id)
        lst.append(mysql_ftr.name)
        lst.append(mysql_ftr.start)
        lst.append(mysql_ftr.stop)
        lst.append(mysql_ftr.orientation)
        lst.append(mysql_ftr.translation)
        lst.append(mysql_ftr.translation_length)
        lst.append(mysql_ftr.raw_description)

        # Gene data checks
        lst.append(mysql_ftr._amino_acid_errors)
        lst.append(mysql_ftr._boundary_error)

    else:
        lst.extend(["","","","","","","","","","","","",""])


    # GenBank feature
    if isinstance(gbk_ftr, CdsFeature):

        # General gene data
        lst.append(gbk_ftr.locus_tag)
        lst.append(gbk_ftr.gene)
        lst.append(gbk_ftr.type)
        lst.append(gbk_ftr.start)
        lst.append(gbk_ftr.stop)
        lst.append(gbk_ftr.orientation)
        lst.append(gbk_ftr.translation)
        lst.append(gbk_ftr.translation_length)
        lst.append(gbk_ftr.raw_product)
        lst.append(gbk_ftr.raw_function)
        lst.append(gbk_ftr.raw_note)

        # Gene data checks
        lst.append(gbk_ftr._amino_acid_errors)
        lst.append(gbk_ftr._boundary_error)
        lst.append(gbk_ftr._locus_tag_missing)
        lst.append(gbk_ftr._locus_tag_typo)
        lst.append(gbk_ftr._description_field_error)
    else:
        lst.extend(["","","","","","","","","","","","","","","",""])


    # MySQL-GenBank checks
    if isinstance(mixed_ftr,MatchedCdsFeatures):

        # If this is a matched CDS feature, both MySQL and
        # GenBank features should have identical unmatched_error value.
        lst.append(mixed_ftr._m_feature._unmatched_error)

        lst.append(mixed_ftr._m_g_different_descriptions)
        lst.append(mixed_ftr._m_g_different_start_sites)
        lst.append(mixed_ftr._m_g_different_translations)
    else:
        lst.append(mixed_ftr._unmatched_error)
        lst.extend(["","",""])
    return lst



def create_genome_data(matched_genomes):
    """Create genome data to output."""

    mysql_gnm = matched_genomes.m_genome
    pdb_genome = matched_genomes.p_genome
    gbk_gnm = matched_genomes.g_genome

    lst = []

    # Genome summary data
    lst.append(mysql_gnm.id)
    lst.append(matched_genomes._contains_errors)

    # MySQL data
    # General genome data
    lst.append(mysql_gnm.name)
    lst.append(mysql_gnm.id)
    lst.append(mysql_gnm.name)
    lst.append(mysql_gnm.annotation_status)
    lst.append(mysql_gnm.cluster)
    lst.append(mysql_gnm.host_genus)
    lst.append(mysql_gnm.accession)
    lst.append(mysql_gnm.length)
    lst.append(mysql_gnm._cds_features_tally)
    lst.append(mysql_gnm._cds_descriptions_tally)
    lst.append(mysql_gnm.retrieve_record)
    lst.append(mysql_gnm.date)
    lst.append(mysql_gnm.annotation_author)


    # Genome data checks
    lst.append(mysql_gnm._nucleotide_errors)
    lst.append(mysql_gnm._cds_features_with_translation_error_tally)
    lst.append(mysql_gnm._cds_features_boundary_error_tally)
    lst.append(mysql_gnm._status_description_error)
    lst.append(mysql_gnm._status_accession_error)


    # PhagesDB data
    if isinstance(pdb_genome, Genome):

        # General genome data
        lst.append(pdb_genome.name)
        lst.append(pdb_genome.name)
        lst.append(pdb_genome.cluster)
        lst.append(pdb_genome.subcluster)
        lst.append(pdb_genome.host_genus)
        lst.append(pdb_genome.accession)
        lst.append(pdb_genome.length)

        # Genome data checks
        lst.append(pdb_genome._nucleotide_errors)
    else:
        lst.extend(["","","","","","","",""])



    # GenBank data
    if isinstance(gbk_gnm, Genome):

        # General genome data
        lst.append(gbk_gnm.name)
        lst.append(gbk_gnm.name)
        lst.append(gbk_gnm.record_id)
        lst.append(gbk_gnm.record_name)
        lst.append(gbk_gnm.accession)
        lst.append(gbk_gnm.description)
        lst.append(gbk_gnm.source)
        lst.append(gbk_gnm.organism)
        lst.append(gbk_gnm.source_feature_organism)
        lst.append(gbk_gnm.source_feature_host)
        lst.append(gbk_gnm.source_feature_lab_host)
        lst.append(gbk_gnm.authors)
        lst.append(gbk_gnm.length)
        lst.append(gbk_gnm._cds_features_tally)


        # Genome data checks
        lst.append(gbk_gnm._nucleotide_errors)
        lst.append(gbk_gnm._cds_features_with_translation_error_tally)
        lst.append(gbk_gnm._cds_features_boundary_error_tally)
        lst.append(gbk_gnm._cds_products_tally)
        lst.append(gbk_gnm._cds_functions_tally)
        lst.append(gbk_gnm._cds_notes_tally)
        lst.append(gbk_gnm._missing_locus_tags_tally)
        lst.append(gbk_gnm._locus_tag_typos_tally)
        lst.append(gbk_gnm._description_field_error_tally)

    else:
        lst.extend(["","","","","","","","","","",
                    "","","","","","","","","","",
                    "","",""])

    # MySQL-PhagesDB checks
    if isinstance(pdb_genome, Genome):
        lst.append(matched_genomes._m_p_seq_mismatch)
        lst.append(matched_genomes._m_p_seq_length_mismatch)
        lst.append(matched_genomes._m_p_cluster_mismatch)
        lst.append(matched_genomes._m_p_accession_mismatch)
        lst.append(matched_genomes._m_p_host_mismatch)
    else:
        lst.extend(["","","","",""])


    # MySQL-GenBank checks
    if isinstance(gbk_gnm, Genome):
        lst.append(matched_genomes._m_g_seq_mismatch)
        lst.append(matched_genomes._m_g_seq_length_mismatch)
        lst.append(matched_genomes._g_header_fields_name_mismatch)
        lst.append(matched_genomes._g_host_mismatch)
        lst.append(matched_genomes._m_g_author_error)
        lst.append(matched_genomes._m_g_perfect_matched_ftrs_tally)
        lst.append(matched_genomes._m_g_imperfect_matched_ftrs_tally)
        lst.append(matched_genomes._m_ftrs_unmatched_in_g_tally)
        lst.append(matched_genomes._g_ftrs_unmatched_in_m_tally)
        lst.append(matched_genomes._m_g_different_descriptions_tally)
        lst.append(matched_genomes._m_g_different_translations_tally)
    else:
        lst.extend(["","","","","","","","","","",""])

    # Number of genes with errors
    lst.append(matched_genomes._total_number_genes_with_errors)

    # Output PhagesDB-GenBank checks
    if isinstance(pdb_genome, Genome) and isinstance(gbk_gnm, Genome):
        lst.append(matched_genomes._p_g_seq_mismatch)
        lst.append(matched_genomes._p_g_seq_length_mismatch)
    else:
        lst.extend(["",""])
    return lst
