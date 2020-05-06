"""Pipeline to compare data between MySQL, PhagesDB, and GenBank databases."""

# TODO this object-oriented pipeline is not fully integrated into
# the pdm_utils package.

# Note this script compares and matches data from GenBank data and MySQL data.
# As a result, there are many similarly named variables.
# Variables are prefixed to indicate database:
# GenBank =  "gbk".
# MySQL = "ph" or "mysql"
# PhagesDB = "pdb"


import argparse
import csv
from datetime import date
import json
import os
import pathlib
import re
import sys
import time
import urllib.request

from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import ncbi
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb, mysqldb_basic

DEFAULT_OUTPUT_FOLDER = os.getcwd()

#Set up dna and protein alphabets to verify sequence integrity
DNA_ALPHABET = set(IUPAC.IUPACUnambiguousDNA.letters)
PROTEIN_ALPHABET = set(IUPAC.ExtendedIUPACProtein.letters)


#Create output directories
CURRENT_DATE = date.today().strftime("%Y%m%d")
RESULTS_FOLDER = f"{CURRENT_DATE}_compare"

GENBANK_OUTPUT = "gbk_records"
MYSQL_OUTPUT = "mysql_records"
PHAGESDB_OUTPUT = "phagesdb_records"

FIRST_HEADER = CURRENT_DATE + " Database comparison"

DUPLICATE_MYSQL_NAMES = "error_duplicate_mysql_phage_names.csv"
DUPLICATE_MYSQL_ACC = "error_duplicate_mysql_phage_accessions.csv"
DUPLICATE_PDB_NAMES = "error_duplicate_phagesdb_phage_names.csv"
FAILED_ACC_RETRIEVE = "error_accession_retrieval.csv"
UNMATCHED_GENOMES = "error_unmatched_genomes.csv"
SUMMARY_OUTPUT = "results_summary.csv"
GENOME_OUTPUT = "results_genome.csv"
GENE_OUTPUT = "results_gene.csv"

PHAGE_QUERY = ("SELECT PhageID, Name, HostGenus, Sequence, Length, "
               "Status, Cluster, Accession, RetrieveRecord, "
               "DateLastModified, AnnotationAuthor FROM phage")
GENE_QUERY = ("SELECT PhageID, GeneID, Name, Start, Stop, Orientation, "
              "Translation, Notes from gene")
VERSION_QUERY = "SELECT Version FROM version"


def database_header(database, version):
    header = f"{database}_v{version}"
    return header

#Define several functions

#Make sure there is no "_Draft" suffix
def remove_draft_suffix(value):
    # Is the word "_Draft" appended to the end of the name?
    value_truncated = value.lower()
    if value_truncated[-6:] == "_draft":
        value_truncated = value_truncated[:-6]
    return value_truncated

def parse_strand(value):
    value = value.lower()
    if value == "f" or value == "forward":
        value = "forward"
    elif value == "r" or value == "reverse":
        value = "reverse"
    else:
        value = "NA"
    return value


#Function to split gene description field
def retrieve_description(description):
    if description is None:
        description = ''
    else:
        description = description.lower().strip()

    split_description = description.split(' ')
    if description == 'hypothetical protein':
        search_description = ''

    elif description == 'phage protein':
        search_description = ''

    elif description == 'unknown':
        search_description = ''

    elif description == '\\n':
        search_description = ''

    elif description.isdigit():
        search_description = ''

    elif len(split_description) == 1:

        if split_description[0][:2] == 'gp' and split_description[0][2:].isdigit():
            search_description = ''

        elif split_description[0][:3] == 'orf' and split_description[0][3:].isdigit():
            search_description = ''

        else:
            search_description = description

    elif len(split_description) == 2:

        if split_description[0] == 'gp' and split_description[1].isdigit():
            search_description = ''

        elif split_description[0] == 'orf' and split_description[1].isdigit():
            search_description = ''

        else:
            search_description = description

    else:
        search_description = description
    return description,search_description



#Function to search through a list of elements using a regular expression
def find_name(expression,list_of_items):
    search_tally = 0
    for element in list_of_items:
        search_result = expression.search(element)
        if search_result:
            search_tally += 1
    return search_tally


def output_to_file(data_list, folder, filename, selected_status,
                   database, version, selected_authors):
    """Output list data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow([FIRST_HEADER])
        writer.writerow([database_header(database, version)])
        writer.writerow([selected_authors])
        writer.writerow([selected_status])
        for element in data_list:
            writer.writerow(element)


def output_failed_accession(data_list, folder, filename, selected_status,
                            selected_dbs, database, selected_authors, version):
    """Output failed accession data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow([FIRST_HEADER])
        writer.writerow([database_header(database, version)])
        writer.writerow([selected_authors])
        writer.writerow([selected_status])
        writer.writerow([selected_dbs])
        writer.writerow(["Accessions unable to be retrieved from GenBank:"])
        for element in data_list:
            writer.writerow(element)



def output_summary_data(data_list, folder, filename, selected_status,
                        selected_dbs, database, selected_authors, version):
    """Output summary data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow([FIRST_HEADER])
        writer.writerow([database_header(database, version)])
        writer.writerow([selected_authors])
        writer.writerow([selected_status])
        writer.writerow([selected_dbs])
        for element in data_list:
            writer.writerow(element)




def output_genome_results(data_list, folder, filename, selected_status,
                          selected_dbs, database, selected_authors, version):
    """Output genome data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow([FIRST_HEADER])
        writer.writerow([database_header(database, version)])
        writer.writerow([selected_authors])
        writer.writerow([selected_status])
        writer.writerow([selected_dbs])
        for element in data_list:
            writer.writerow(element)






def output_unmatched(data_list, folder, filename, selected_status,
                     selected_dbs, database, selected_authors, version):
    """Output unmatched genomes to file."""

    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow([FIRST_HEADER])
        writer.writerow([database_header(database, version)])
        writer.writerow([selected_authors])
        writer.writerow([selected_status])
        writer.writerow([selected_dbs])
        writer.writerow([""])
        for element in data_list:
            writer.writerow(element)


def prepare_unmatched_output(unmatched_pdb, unmatched_gbk):
    """Prepare list of unmatched data to be outputted to file."""
    unmatched_rows = []
    if len(unmatched_pdb) > 0:
        unmatched_rows.append(["The following MySQL genomes were not matched to PhagesDB:"])
        unmatched_rows.append(["PhageID", "PhageName", "Author", "Status"])
        for element in unmatched_pdb:
            unmatched_rows.append([element.id,
                                   element.name,
                                   element.annotation_author,
                                   element.annotation_status])
    if len(unmatched_gbk) > 0:
        unmatched_rows.append([""])
        unmatched_rows.append(["The following MySQL genomes were not matched to GenBank:"])
        unmatched_rows.append(["PhageID","PhageName","Author","Status","Accession"])
        for element in unmatched_gbk:
            unmatched_rows.append([element.id,
                                   element.name,
                                   element.annotation_author,
                                   element.annotation_status,
                                   element.accession])
    return unmatched_rows


# TODO move to basic module
# TODO test
def make_new_file(output_dir, new_file, ext, attempt=1):
    """Make a new file.

    Checks to verify the new file name is valid and does not
    already exist. If it already exists, it attempts to extend
    the name with an integer suffix.

    :param output_dir:
        Full path to the directory where the new directory will be created.
    :type output_dir: Path
    :param new_file: Name of the new file to be created.
    :type new_file: Path
    :param attempt: Number of attempts to create the file.
    :type attempt: int
    :returns:
        If successful, the full path of the created file.
        If unsuccessful, None.
    :rtype: Path, None
    """
    valid = False
    count = 0
    while (not valid and count < attempt):
        if count > 0:
            new_file_mod = new_file.stem + "_" + str(count)
            new_file_mod = pathlib.Path(new_file_mod)
        else:
            new_file_mod = new_file

        new_file_mod = new_file_mod + "." + ext
        new_path = pathlib.Path(output_dir, new_file_mod)
        if new_path.is_file() == False:
            valid = True
        count += 1
    if not valid:
        return None
    else:
        return new_path


def save_seqrecord(seqrecord, output_path, file_prefix, ext, seqrecord_ext):
    """Save record to file."""
    file_path = make_new_file(output_path, file_prefix, ext, attempt=100)
    if file_path is not None:
        SeqIO.write(seqrecord, file_path, seqrecord_ext)
    else:
        print(f"Duplicated filenames for {file_prefix}. Unable to output data.")
        input("Press ENTER to continue")





# Define classes



# Base genome class
class UnannotatedGenome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields
        self.type = "" # Identifier to describes source of this genome
                       # (e.g. mysql, phagesdb, genbank, etc.)


        self.name = ""
        self.host_genus = ""
        self.sequence = "" # TODO string but seq is Seq
        self.seq = Seq("", IUPAC.ambiguous_dna)
        self.accession = ""


        # Computed datafields
        self._search_name = "" # No "_Draft" and converted to lowercase
        self.length = 0
        self._nucleotide_errors = False


    # Define all attribute setters:
    def set_phage_name(self,value):
        self.name = value
        self._search_name = remove_draft_suffix(self.name)
    def set_sequence(self,value):
        self.sequence = value.upper()
        self.length = len(self.sequence)
    def set_accession(self,value):
        if value is None or value.strip() == "":
            self.accession = ""
        else:
            value = value.strip()
            self.accession = value.split(".")[0]
    def compute_nucleotide_errors(self,DNA_ALPHABET):
        nucleotide_set = set(self.sequence)
        nucleotide_error_set = nucleotide_set - DNA_ALPHABET
        if len(nucleotide_error_set) > 0:
            self._nucleotide_errors = True






class AnnotatedGenome(UnannotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        UnannotatedGenome.__init__(self)

        # Non-computed datafields

        # Computed datafields
        self.cds_features = []
        self._cds_features_tally = 0
        self._cds_features_with_translation_error_tally = 0
        self._cds_features_boundary_error_tally = 0

    # Define all attribute setters:
    def set_cds_features(self,value):
        self.cds_features = value # Should be a list
        self._cds_features_tally = len(self.cds_features)
    def compute_cds_feature_errors(self):
        for cds_feature in self.cds_features:
            if cds_feature._amino_acid_errors:
                self._cds_features_with_translation_error_tally += 1
            if cds_feature._boundary_error:
                self._cds_features_boundary_error_tally += 1






class PhameratorGenome(AnnotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        AnnotatedGenome.__init__(self)
        # Non-computed datafields
        self.id = ""
        self.annotation_status = "" # Final, Draft, Unknown version of genome data
        self.cluster_subcluster = "" # Combined cluster_subcluster data.
        self.retrieve_record = ""
        self.date = ""
        self.annotation_author = "" # 1 (Hatfull), 0 (GenBank)

        # Computed datafields
        self._search_id = "" # No "_Draft" and converted to lowercase
        self._status_accession_error = False
        self._status_description_error = False
        self._cds_descriptions_tally = 0
        self._genes_with_errors_tally = 0 # Number of genes with >1 error


    # Define all attribute setters:
    def set_phage_id(self,value):
        self.id = value
        self._search_id = remove_draft_suffix(self.id)
    def set_status(self,value):
        self.annotation_status = value

        # Be sure to first set the accession attribute before the annotation_status attribute,
        # else this will throw an error.
        # Now that the AnnotationAuthor field contains authorship data, the
        # 'unknown' annotation annotation_status now reflects an 'unknown' annotation (in
        # regards to if it was auto-annotated or manually annotated).
        # So for the annotation_status-accession error, if the annotation_status is 'unknown',
        # there is no reason to assume whether there should be an accession or not.
        # Only for 'final' (manually annotated) genomes should there be an accession.
        if self.annotation_status == "final" and self.accession == "":
            self._status_accession_error = True



    def set_cluster_subcluster(self,value):
        if value is None:
            self.cluster_subcluster = "Singleton"
        elif value == "UNK":
            self.cluster_subcluster = ""
        else:
            self.cluster_subcluster = value
    def set_date_last_modified(self,value):
        if value is None:
            self.date = ""
        else:
            self.date = value

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


    #Even though this method iterates through the cds features
    # like the compute_status_description_error does,
    #it has to be kept separate, since you need to wait to run
    # this method after all genome and gene matching is completed.
    def compute_genes_with_errors_tally(self):
        for feature in self.cds_features:
            # Need to first compute the number of errors per gene
            feature.compute_total_cds_errors()
            if feature._total_errors > 0:
                self._genes_with_errors_tally += 1




class PhagesdbGenome(UnannotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        UnannotatedGenome.__init__(self)

        # Non-computed datafields
        self.cluster = ""
        self.subcluster = ""




class NcbiGenome(AnnotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        AnnotatedGenome.__init__(self)

        # Non-computed data fields
        self.record_name = "" # TODO no equivalent in ORM
        self.record_id = "" # TODO no equivalent in ORM
        self.accession = ""
        self.description = ""
        self.source = ""
        self.organism = ""
        self.source_feature_organism = "" # TODO move to source feature
        self.source_feature_host = "" # TODO move to source feature
        self.source_feature_lab_host = "" # TODO move to source feature
        self.authors = ""


        # Computed data fields
        self._cds_functions_tally = 0
        self._cds_products_tally = 0
        self._cds_notes_tally = 0
        self._missing_locus_tags_tally = 0 # TODO no equivalent in ORM
        self._locus_tag_typos_tally = 0 # TODO no equivalent in ORM
        self._description_field_error_tally = 0 # TODO no equivalent in ORM



    # Define setter functions
    def compute_ncbi_cds_feature_errors(self):
        for cds_feature in self.cds_features:

            # counting descriptions should skip if it is blank
            # or "hypothetical protein".
            if cds_feature._search_product_description != "":
                self._cds_products_tally += 1

            if cds_feature._search_function_description != "":
                self._cds_functions_tally += 1

            if cds_feature._search_note_description != "":
                self._cds_notes_tally += 1



            if cds_feature._locus_tag_missing:
                self._missing_locus_tags_tally += 1
            else:
                pattern4 = re.compile(self._search_name)
                search_result = pattern4.search(cds_feature.locus_tag.lower())

                if search_result == None:
                    self._locus_tag_typos_tally += 1
                    cds_feature.set_locus_tag_typo() # Sets to True

            if cds_feature._description_field_error:
                self._description_field_error_tally += 1





class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        # Datafields from MySQL database:
        self.type = "" # Feature type: CDS, GenomeBoundary,or tRNA
        self.start = "" # Position of left boundary, 0-indexed
        self.stop = "" # Position of right boundary, 0-indexed
        self.orientation = "" # 'forward', 'reverse', or 'NA'
        self.translation = ""
        self.translation_length = ""

        # Computed datafields
        self._amino_acid_errors = False # TODO no equivalent in ORM
        self._start_end_strand_id = "" # TODO no equivalent in ORM
        self._end_strand_id = "" # TODO no equivalent in ORM
        self._boundary_error = False # TODO no equivalent in ORM
        self._unmatched_error = False # Tracks if it contains a match or not # TODO no equivalent in ORM

    # Define all attribute setters:
    def set_strand(self,value):
        self.orientation = parse_strand(value)
    def set_translation(self,value):
        self.translation = value.upper()
        self.translation_length = len(self.translation)
    def compute_amino_acid_errors(self,PROTEIN_ALPHABET):
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - PROTEIN_ALPHABET
        if len(amino_acid_error_set) > 0:
            self._amino_acid_errors = True
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
    def set_unmatched_error(self):
        self._unmatched_error = True
    def compute_boundary_error(self):
        # Check if start and end coordinates are fuzzy
        if not (str(self.start).isdigit() and str(self.stop).isdigit()):
            self._boundary_error = True





class MysqlCdsFeature(CdsFeature):

    # Initialize all attributes:
    def __init__(self):
        CdsFeature.__init__(self)

        # Initialize all non-calculated attributes:

        # Datafields from MySQL database:
        self.genome_id = ""
        self.id = "" # Gene ID comprised of PhageID and Gene name
        self.name = ""
        self.raw_description = ""
        self.description = "" # non-generic gene descriptions

        # Computed datafields
        self._search_genome_id = "" # TODO no equivalent in ORM
        self._total_errors = 0 # TODO no equivalent in ORM

    # Define all attribute setters:
    def set_phage_id(self,value):
        self.genome_id = value
        self._search_genome_id = remove_draft_suffix(self.genome_id)
    def set_notes(self,value1,value2):
        self.raw_description = value1
        self.description = value2
    def compute_total_cds_errors(self):
        if self._amino_acid_errors:
            self._total_errors += 1
        if self._boundary_error:
            self._total_errors += 1
        if self._unmatched_error:
            self._total_errors += 1






class NcbiCdsFeature(CdsFeature):

    # Initialize all attributes:
    def __init__(self):
        CdsFeature.__init__(self)

        # Initialize all non-calculated attributes:
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self.gene_number = ""
        self.product_description = ""
        self.function_description = ""
        self.note_description = ""
        self._search_product_description = ""
        self._search_function_description = ""
        self._search_note_description = ""

        # Inititalize all calculated attributes:
        self._locus_tag_missing = False

        # Computed at genome level since it uses phage name:
        self._locus_tag_typo = False
        self._description_field_error = False
        self._total_errors = 0


    # Define all attribute setters:
    def set_locus_tag(self,value):
        self.locus_tag = value
        if self.locus_tag == "":
            self._locus_tag_missing = True
    def set_gene_number(self,value):
        self.gene_number = value
    def set_product_description(self,value1,value2):
        self.product_description = value1
        self._search_product_description = value2
    def set_function_description(self,value1,value2):
        self.function_description = value1
        self._search_function_description = value2
    def set_note_description(self,value1,value2):
        self.note_description = value1
        self._search_note_description = value2
    def set_locus_tag_typo(self):
        self._locus_tag_typo = True

    def compute_description_error(self):

        # If the product description is empty or generic,
        # and the function or note descriptions are not, there is an error.
        if self._search_product_description == "" and \
            (self._search_function_description != "" or \
            self._search_note_description != ""):

            self._description_field_error = True

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





class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.m_genome = ""
        self.p_genome = ""
        self.g_genome = ""

        # MySQL and GenBank matched data comparison results
        self._phamerator_ncbi_sequence_mismatch = False
        self._phamerator_ncbi_sequence_length_mismatch = False
        self._ncbi_record_header_fields_phage_name_mismatch = False
        self._ncbi_host_mismatch = False

        # List of MatchedCdsFeature objects
        self._phamerator_ncbi_perfect_matched_features = []
        self._phamerator_ncbi_imperfect_matched_features = []

        # List of CdsFeature objects
        self._phamerator_features_unmatched_in_ncbi = []
        self._ncbi_features_unmatched_in_phamerator = []

        self._phamerator_ncbi_perfect_matched_features_tally = 0
        self._phamerator_ncbi_imperfect_matched_features_tally = 0
        self._phamerator_features_unmatched_in_ncbi_tally = 0
        self._ncbi_features_unmatched_in_phamerator_tally = 0
        self._phamerator_ncbi_different_descriptions_tally = 0
        self._phamerator_ncbi_different_translations_tally = 0
        self._ph_ncbi_author_error = False





        # MySQL and PhagesDB matched data comparison results
        self._phamerator_phagesdb_sequence_mismatch = False
        self._phamerator_phagesdb_sequence_length_mismatch = False
        self._phamerator_phagesdb_host_mismatch = False
        self._phamerator_phagesdb_accession_mismatch = False
        self._phamerator_phagesdb_cluster_subcluster_mismatch = False


        # PhagesDB and GenBank matched data comparison results
        self._phagesdb_ncbi_sequence_mismatch = False
        self._phagesdb_ncbi_sequence_length_mismatch = False

        # Total errors summary
        self._contains_errors = False
        self._total_number_genes_with_errors = 0



    # Define all attribute setters:
    def set_phamerator_genome(self,value):
        self.m_genome = value
    def set_phagesdb_genome(self,value):
        self.p_genome = value
    def set_ncbi_genome(self,value):
        self.g_genome = value

    def compare_phamerator_ncbi_genomes(self):

        # verify that there is a MySQL and GenBank genome in
        # the matched genome object.
        ph_genome = self.m_genome
        ncbi_genome = self.g_genome
        ph_cds_list = ph_genome.cds_features

        if isinstance(ph_genome,PhameratorGenome) and isinstance(ncbi_genome,NcbiGenome):

            if ph_genome.sequence != ncbi_genome.sequence:
                self._phamerator_ncbi_sequence_mismatch = True
            if ph_genome.length != ncbi_genome.length:
                self._phamerator_ncbi_sequence_length_mismatch = True

            # Compare phage names
            pattern1 = re.compile("^" + ph_genome.name + "$")
            pattern2 = re.compile("^" + ph_genome.name)

            if find_name(pattern2,ncbi_genome.description.split(" ")) == 0 or \
                find_name(pattern1,ncbi_genome.source.split(" ")) == 0 or \
                find_name(pattern1,ncbi_genome.organism.split(" ")) == 0 or \
                find_name(pattern1,ncbi_genome.source_feature_organism.split(" ")) == 0:

                self._ncbi_record_header_fields_phage_name_mismatch = True

            # Compare host_genus data
            search_host = ph_genome.host_genus
            if search_host == "Mycobacterium":
                search_host = search_host[:-3]
            pattern3 = re.compile("^" + search_host)

            if (find_name(pattern3,ncbi_genome.description.split(" ")) == 0 or \
                find_name(pattern3,ncbi_genome.source.split(" ")) == 0 or \
                find_name(pattern3,ncbi_genome.organism.split(" ")) == 0 or \
                find_name(pattern3,ncbi_genome.source_feature_organism.split(" ")) == 0) or \
                (ncbi_genome.source_feature_host != "" and find_name(pattern3,ncbi_genome.source_feature_host.split(" ")) == 0) or \
                (ncbi_genome.source_feature_lab_host != "" and find_name(pattern3,ncbi_genome.source_feature_lab_host.split(" ")) == 0):

                self._ncbi_host_mismatch = True

            # Check author list for errors
            # For genomes with AnnotationAuthor = 1 (Hatfull), Graham is expected
            # to be an author.
            # For genomes with AnnotationAuthor = 0 (non-Hatfull/GenBank), Graham
            # is NOT expected to be an author.
            pattern5 = re.compile("hatfull")
            search_result = pattern5.search(ncbi_genome.authors.lower())
            if ph_genome.annotation_author == 1 and search_result == None:
                self._ph_ncbi_author_error = True
            elif ph_genome.annotation_author == 0 and search_result != None:
                self._ph_ncbi_author_error = True
            else:
                # Any other combination of MySQL and GenBank author can be skipped
                pass


            # Compare CDS features

            # First find all unique start-end-orientation cds identifiers
            # for MySQL and GenBank genomes.
            ph_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            ph_start_end_strand_duplicate_id_set = set()

            for cds in ph_cds_list:
                if cds._start_end_strand_id not in ph_start_end_strand_id_set:
                    ph_start_end_strand_id_set.add(cds._start_end_strand_id)
                else:
                    ph_start_end_strand_duplicate_id_set.add(cds._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            ph_start_end_strand_id_set = ph_start_end_strand_id_set - ph_start_end_strand_duplicate_id_set


            ncbi_cds_list = ncbi_genome.cds_features
            ncbi_start_end_strand_id_set = set()

            # All end_strand ids that are not unique
            ncbi_start_end_strand_duplicate_id_set = set()
            for cds in ncbi_cds_list:
                if cds._start_end_strand_id not in ncbi_start_end_strand_id_set:
                    ncbi_start_end_strand_id_set.add(cds._start_end_strand_id)

                else:
                    ncbi_start_end_strand_duplicate_id_set.add(cds._start_end_strand_id)
            # Remove the duplicate end_strand ids from the main id_set
            ncbi_start_end_strand_id_set = ncbi_start_end_strand_id_set - ncbi_start_end_strand_duplicate_id_set

            # Create the perfect matched and unmatched sets
            ph_unmatched_start_end_strand_id_set = ph_start_end_strand_id_set - ncbi_start_end_strand_id_set
            ncbi_unmatched_start_end_strand_id_set = ncbi_start_end_strand_id_set - ph_start_end_strand_id_set
            perfect_matched_cds_id_set = ph_start_end_strand_id_set & ncbi_start_end_strand_id_set


            # From the unmatched sets, created second round of end-orientation id sets
            # MySQL end_strand data
            ph_end_strand_id_set = set()

            # All end_strand ids that are not unique
            ph_end_strand_duplicate_id_set = set()
            for cds in ph_cds_list:
                if cds._start_end_strand_id in ph_unmatched_start_end_strand_id_set:
                    if cds._end_strand_id not in ph_end_strand_id_set:
                        ph_end_strand_id_set.add(cds._end_strand_id)
                    else:
                        ph_end_strand_duplicate_id_set.add(cds._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            ph_end_strand_id_set = ph_end_strand_id_set - ph_end_strand_duplicate_id_set


            ncbi_end_strand_id_set = set()

            # All end_strand ids that are not unique
            ncbi_end_strand_duplicate_id_set = set()
            for cds in ncbi_cds_list:
                if cds._start_end_strand_id in ncbi_unmatched_start_end_strand_id_set:
                    if cds._end_strand_id not in ncbi_end_strand_id_set:
                        ncbi_end_strand_id_set.add(cds._end_strand_id)
                    else:
                        ncbi_end_strand_duplicate_id_set.add(cds._end_strand_id)

            # Remove the duplicate end_strand ids from the main id_set
            ncbi_end_strand_id_set = ncbi_end_strand_id_set - ncbi_end_strand_duplicate_id_set


            # Create the imperfect matched set
            imperfect_matched_cds_id_set = ph_end_strand_id_set & ncbi_end_strand_id_set


            # Now go back through all cds features and assign
            # to the appropriate dictionary or list.
            ph_perfect_matched_cds_dict = {}
            ph_imperfect_matched_cds_dict = {}
            ph_unmatched_cds_list = []


            for cds in ph_cds_list:
                if cds._start_end_strand_id in perfect_matched_cds_id_set:
                    ph_perfect_matched_cds_dict[cds._start_end_strand_id] = cds
                elif cds._end_strand_id in imperfect_matched_cds_id_set:
                    ph_imperfect_matched_cds_dict[cds._end_strand_id] = cds
                else:
                    ph_unmatched_cds_list.append(cds)


            ncbi_perfect_matched_cds_dict = {}
            ncbi_imperfect_matched_cds_dict = {}
            ncbi_unmatched_cds_list = []

            for cds in ncbi_cds_list:
                if cds._start_end_strand_id in perfect_matched_cds_id_set:
                    ncbi_perfect_matched_cds_dict[cds._start_end_strand_id] = cds
                elif cds._end_strand_id in imperfect_matched_cds_id_set:
                    ncbi_imperfect_matched_cds_dict[cds._end_strand_id] = cds
                else:
                    ncbi_unmatched_cds_list.append(cds)


            # Create MatchedCdsFeatures objects
            # Compute matched gene errors and tallies
            # Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1



                if matched_cds_object._phamerator_ncbi_different_translations:
                    self._phamerator_ncbi_different_translations_tally += 1
                if matched_cds_object._phamerator_ncbi_different_descriptions:
                    self._phamerator_ncbi_different_descriptions_tally += 1
                self._phamerator_ncbi_perfect_matched_features.append(matched_cds_object)


            # Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object._total_errors > 0:
                    self._total_number_genes_with_errors += 1

                if matched_cds_object._phamerator_ncbi_different_descriptions:
                    self._phamerator_ncbi_different_descriptions_tally += 1
                self._phamerator_ncbi_imperfect_matched_features.append(matched_cds_object)


            # Compute unmatched error and gene total errors for
            # all unmatched features.
            for cds in ph_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds._total_errors > 0:
                    self._total_number_genes_with_errors += 1

            for cds in ncbi_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds._total_errors > 0:
                    self._total_number_genes_with_errors += 1




            # Set unmatched cds lists
            self._phamerator_features_unmatched_in_ncbi = ph_unmatched_cds_list
            self._ncbi_features_unmatched_in_phamerator = ncbi_unmatched_cds_list

            # Now compute the number of features in each category
            self._phamerator_ncbi_perfect_matched_features_tally = len(self._phamerator_ncbi_perfect_matched_features)
            self._phamerator_ncbi_imperfect_matched_features_tally = len(self._phamerator_ncbi_imperfect_matched_features)
            self._phamerator_features_unmatched_in_ncbi_tally = len(self._phamerator_features_unmatched_in_ncbi)
            self._ncbi_features_unmatched_in_phamerator_tally = len(self._ncbi_features_unmatched_in_phamerator)


        # If there is no matching GenBank genome,
        # assign all MySQL genes to Unmatched.
        else:

            # Set unmatched cds lists, but do NOT count them in unmatched tally.
            # Unmatched tally should reflect unmatched genes if there is
            # actually a matching GenBank genome.
            self._phamerator_features_unmatched_in_ncbi = ph_cds_list

            # Now that all errors have been computed for each gene,
            # compute how many genes have errors
            # If there is no matching GenBank genome,
            # gene error tallies are only computed for the MySQL genome.
            ph_genome.compute_genes_with_errors_tally()
            self._total_number_genes_with_errors = ph_genome._genes_with_errors_tally




    def compare_phamerator_phagesdb_genomes(self):

        # verify that there is a MySQL and PhagesDB genome
        # in the matched genome object.
        ph_genome = self.m_genome
        pdb_genome = self.p_genome

        if isinstance(ph_genome,PhameratorGenome) and isinstance(pdb_genome,PhagesdbGenome):

            if ph_genome.sequence != pdb_genome.sequence:
                self._phamerator_phagesdb_sequence_mismatch = True
            if ph_genome.length != pdb_genome.length:
                self._phamerator_phagesdb_sequence_length_mismatch = True
            if ph_genome.accession != pdb_genome.accession:
                self._phamerator_phagesdb_accession_mismatch = True
            if ph_genome.host_genus != pdb_genome.host_genus:
                self._phamerator_phagesdb_host_mismatch = True
            if ph_genome.cluster_subcluster != pdb_genome.cluster and \
                ph_genome.cluster_subcluster != pdb_genome.subcluster:

                self._phamerator_phagesdb_cluster_subcluster_mismatch = True


    def compare_phagesdb_ncbi_genomes(self):

        # verify that there is a PhagesDB and GenBank genome
        # in the matched genome object.
        pdb_genome = self.p_genome
        ncbi_genome = self.g_genome

        if isinstance(pdb_genome,PhagesdbGenome) and isinstance(ncbi_genome,NcbiGenome):
            if pdb_genome.sequence != ncbi_genome.sequence:
                self._phagesdb_ncbi_sequence_mismatch = True
            if pdb_genome.length != ncbi_genome.length:
                self._phagesdb_ncbi_sequence_length_mismatch = True



    def compute_total_genome_errors(self):

        ph_genome = self.m_genome
        ncbi_genome = self.g_genome
        pdb_genome = self.p_genome


        if self._total_number_genes_with_errors > 0:
            self._contains_errors = True


        # MySQL genome
        if ph_genome._nucleotide_errors:
            self._contains_errors = True
        if ph_genome._status_description_error:
            self._contains_errors = True
        if ph_genome._status_accession_error:
            self._contains_errors = True


        # GenBank genome
        if isinstance(ncbi_genome,NcbiGenome):
            if ncbi_genome._nucleotide_errors:
                self._contains_errors = True



        # PhagesDB genome
        if isinstance(pdb_genome,PhagesdbGenome):
            if pdb_genome._nucleotide_errors:
                self._contains_errors = True


        # MySQL-GenBank
        if self._phamerator_ncbi_sequence_mismatch:
            self._contains_errors = True
        if self._phamerator_ncbi_sequence_length_mismatch:
            self._contains_errors = True
        if self._ncbi_record_header_fields_phage_name_mismatch:
            self._contains_errors = True
        if self._ncbi_host_mismatch:
            self._contains_errors = True
        if self._phamerator_ncbi_imperfect_matched_features_tally > 0:
            self._contains_errors = True
        if self._phamerator_features_unmatched_in_ncbi_tally > 0:
            self._contains_errors = True
        if self._ncbi_features_unmatched_in_phamerator_tally > 0:
            self._contains_errors = True
        if self._phamerator_ncbi_different_descriptions_tally > 0:
            self._contains_errors = True
        if self._phamerator_ncbi_different_translations_tally > 0:
            self._contains_errors = True
        if self._ph_ncbi_author_error:
            self._contains_errors = True


        # MySQL-PhagesDB
        if self._phamerator_phagesdb_sequence_mismatch:
            self._contains_errors = True
        if self._phamerator_phagesdb_sequence_length_mismatch:
            self._contains_errors = True
        if self._phamerator_phagesdb_host_mismatch:
            self._contains_errors = True
        if self._phamerator_phagesdb_accession_mismatch:
            self._contains_errors = True
        if self._phamerator_phagesdb_cluster_subcluster_mismatch:
            self._contains_errors = True

        # PhagesDB-GenBank
        if self._phagesdb_ncbi_sequence_mismatch:
            self._contains_errors = True
        if self._phagesdb_ncbi_sequence_length_mismatch:
            self._contains_errors = True








class MatchedCdsFeatures:

    # Initialize all attributes:
    def __init__(self):


        # Initialize all non-calculated attributes:
        self._m_feature = ""
        self._g_ftr = ""

        # Matched data comparison results
        self._phamerator_ncbi_different_translations = False
        self._phamerator_ncbi_different_start_sites = False
        self._phamerator_ncbi_different_descriptions = False

        # Total errors summary
        self._total_errors = 0



    # Define all attribute setters:
    def set_phamerator_feature(self,value):
        self._m_feature = value
    def set_ncbi_feature(self,value):
        self._g_ftr = value


    def compare_phamerator_ncbi_cds_features(self):

        if self._m_feature.orientation == "forward":
            if str(self._m_feature.start) != str(self._g_ftr.start):
                self._phamerator_ncbi_different_start_sites = True
        elif self._m_feature.orientation == "reverse":
            if str(self._m_feature.stop) != str(self._g_ftr.stop):
                self._phamerator_ncbi_different_start_sites = True
        else:
            pass


        product_description_set = set()
        product_description_set.add(self._m_feature.description)
        product_description_set.add(self._g_ftr._search_product_description)


        if len(product_description_set) != 1:
            self._phamerator_ncbi_different_descriptions = True

        if self._m_feature.translation != self._g_ftr.translation:
            self._phamerator_ncbi_different_translations = True



        # Compute total errors
        # First add all matched feature errors
        if self._phamerator_ncbi_different_translations:
            self._total_errors += 1

        if self._phamerator_ncbi_different_start_sites:
            self._total_errors += 1

        if self._phamerator_ncbi_different_descriptions:
            self._total_errors += 1

        # Now add all errors from each individual feature
        # You first compute errors for each individual feature.
        # This step is performed here instead of in the mainline code
        # because you need to wait for the feature matching step
        # after the genome matching step.
        self._m_feature.compute_total_cds_errors()
        self._g_ftr.compute_total_cds_errors()
        self._total_errors += self._m_feature._total_errors
        self._total_errors += self._g_ftr._total_errors











class DatabaseSummary:

    # Initialize all attributes:
    def __init__(self,matched_genomes_list):

        # Initialize all non-calculated attributes:
        self._matched_genomes_list = matched_genomes_list

        # Initialize all calculated attributes:

        # MySQL data
        # General genome data
        self._ph_ncbi_update_flag_tally = 0

        # Genome data checks
        self._ph_genomes_with_nucleotide_errors_tally = 0
        self._ph_genomes_with_translation_errors_tally = 0
        self._ph_genomes_with_boundary_errors_tally = 0
        self._ph_genomes_with_status_accession_error_tally = 0
        self._ph_genomes_with_status_description_error_tally = 0

        # PhagesDB data
        # Genome data checks
        self._pdb_genomes_with_nucleotide_errors_tally = 0

        # GenBank data
        # Genome data checks
        self._ncbi_genomes_with_nucleotide_errors_tally = 0
        self._ncbi_genomes_with_translation_errors_tally = 0
        self._ncbi_genomes_with_boundary_errors_tally = 0
        self._ncbi_genomes_with_missing_locus_tags_tally = 0
        self._ncbi_genomes_with_locus_tag_typos_tally = 0
        self._ncbi_genomes_with_description_field_errors_tally = 0

        # MySQL-PhagesDB checks
        self._ph_pdb_sequence_mismatch_tally = 0
        self._ph_pdb_sequence_length_mismatch_tally = 0
        self._ph_pdb_cluster_subcluster_mismatch_tally = 0
        self._ph_pdb_accession_mismatch_tally = 0
        self._ph_pdb_host_mismatch_tally = 0

        # MySQL-GenBank checks
        self._ph_ncbi_sequence_mismatch_tally = 0
        self._ph_ncbi_sequence_length_mismatch_tally = 0
        self._ph_ncbi_record_header_phage_mismatch_tally = 0
        self._ph_ncbi_record_header_host_mismatch_tally = 0
        self._ph_ncbi_genomes_with_imperfectly_matched_features_tally = 0
        self._ph_ncbi_genomes_with_unmatched_phamerator_features_tally = 0
        self._ph_ncbi_genomes_with_unmatched_ncbi_features_tally = 0
        self._ph_ncbi_genomes_with_different_descriptions_tally = 0
        self._ph_ncbi_genomes_with_different_translations_tally = 0
        self._ph_ncbi_genomes_with_author_errors_tally = 0

        # PhagesDB-GenBank checks
        self._pdb_ncbi_sequence_mismatch_tally = 0
        self._pdb_ncbi_sequence_length_mismatch_tally = 0


        # MySQL feature
        # Gene data checks
        self._ph_translation_errors_tally = 0
        self._ph_boundary_errors_tally = 0


        # GenBank feature
        # Gene data checks
        self._ncbi_translation_errors_tally = 0
        self._ncbi_boundary_errors_tally = 0
        self._ncbi_missing_locus_tags_tally = 0
        self._ncbi_locus_tag_typos_tally = 0
        self._ncbi_description_field_errors_tally = 0

        # MySQL-GenBank checks
        self._ph_ncbi_different_descriptions_tally = 0
        self._ph_ncbi_different_start_sites_tally = 0
        self._ph_ncbi_different_translation_tally = 0
        self._ph_ncbi_unmatched_phamerator_features_tally = 0
        self._ph_ncbi_unmatched_ncbi_features_tally = 0


        # Calculate summary metrics
        self._ph_total_genomes_analyzed = 0
        self._ph_genomes_unmatched_to_pdb_tally = 0
        self._ph_genomes_unmatched_to_ncbi_tally = 0
        self._total_genomes_with_errors = 0


    # Define setter functions
    def compute_summary(self):
        for matched_genomes in self._matched_genomes_list:

            self._ph_total_genomes_analyzed += 1
            ph_genome = matched_genomes.m_genome
            pdb_genome = matched_genomes.p_genome
            ncbi_genome = matched_genomes.g_genome

            if matched_genomes._contains_errors:
                self._total_genomes_with_errors += 1


            # MySQL data
            if isinstance(ph_genome,PhameratorGenome):

                self._ph_ncbi_update_flag_tally += ph_genome.retrieve_record

                # Genome data checks
                if ph_genome._nucleotide_errors:
                    self._ph_genomes_with_nucleotide_errors_tally += 1
                if ph_genome._cds_features_with_translation_error_tally > 0:
                    self._ph_genomes_with_translation_errors_tally += 1
                    self._ph_translation_errors_tally += ph_genome._cds_features_with_translation_error_tally
                if ph_genome._cds_features_boundary_error_tally > 0:
                    self._ph_genomes_with_boundary_errors_tally += 1
                    self._ph_boundary_errors_tally += ph_genome._cds_features_boundary_error_tally
                if ph_genome._status_accession_error:
                    self._ph_genomes_with_status_accession_error_tally += 1
                if ph_genome._status_description_error:
                    self._ph_genomes_with_status_description_error_tally += 1

            # PhagesDB data
            if isinstance(pdb_genome,PhagesdbGenome):

                # Genome data checks
                if pdb_genome._nucleotide_errors:
                    self._pdb_genomes_with_nucleotide_errors_tally += 1
            else:
                self._ph_genomes_unmatched_to_pdb_tally += 1

            # GenBank data
            if isinstance(ncbi_genome,NcbiGenome):

                # Genome data checks
                if ncbi_genome._nucleotide_errors:
                    self._ncbi_genomes_with_nucleotide_errors_tally += 1
                if ncbi_genome._cds_features_with_translation_error_tally > 0:
                    self._ncbi_genomes_with_translation_errors_tally += 1
                    self._ncbi_translation_errors_tally += ncbi_genome._cds_features_with_translation_error_tally
                if ncbi_genome._cds_features_boundary_error_tally > 0:
                    self._ncbi_genomes_with_boundary_errors_tally += 1
                    self._ncbi_boundary_errors_tally += ncbi_genome._cds_features_boundary_error_tally
                if ncbi_genome._missing_locus_tags_tally > 0:
                    self._ncbi_genomes_with_missing_locus_tags_tally += 1
                    self._ncbi_missing_locus_tags_tally += ncbi_genome._missing_locus_tags_tally
                if ncbi_genome._locus_tag_typos_tally > 0:
                    self._ncbi_genomes_with_locus_tag_typos_tally += 1
                    self._ncbi_locus_tag_typos_tally += ncbi_genome._locus_tag_typos_tally
                if ncbi_genome._description_field_error_tally > 0:
                    self._ncbi_genomes_with_description_field_errors_tally += 1
                    self._ncbi_description_field_errors_tally += ncbi_genome._description_field_error_tally

            else:
                self._ph_genomes_unmatched_to_ncbi_tally += 1

            # MySQL-PhagesDB checks
            if matched_genomes._phamerator_phagesdb_sequence_mismatch:
                self._ph_pdb_sequence_mismatch_tally += 1
            if matched_genomes._phamerator_phagesdb_sequence_length_mismatch:
                self._ph_pdb_sequence_length_mismatch_tally += 1
            if matched_genomes._phamerator_phagesdb_cluster_subcluster_mismatch:
                self._ph_pdb_cluster_subcluster_mismatch_tally += 1
            if matched_genomes._phamerator_phagesdb_accession_mismatch:
                self._ph_pdb_accession_mismatch_tally += 1
            if matched_genomes._phamerator_phagesdb_host_mismatch:
                self._ph_pdb_host_mismatch_tally += 1

            # MySQL-GenBank checks
            if matched_genomes._phamerator_ncbi_sequence_mismatch:
                self._ph_ncbi_sequence_mismatch_tally += 1
            if matched_genomes._phamerator_ncbi_sequence_length_mismatch:
                self._ph_ncbi_sequence_length_mismatch_tally += 1
            if matched_genomes._ncbi_record_header_fields_phage_name_mismatch:
                self._ph_ncbi_record_header_phage_mismatch_tally += 1
            if matched_genomes._ncbi_host_mismatch:
                self._ph_ncbi_record_header_host_mismatch_tally += 1
            if matched_genomes._phamerator_ncbi_imperfect_matched_features_tally > 0:
                self._ph_ncbi_genomes_with_imperfectly_matched_features_tally += 1
                self._ph_ncbi_different_start_sites_tally += matched_genomes._phamerator_ncbi_imperfect_matched_features_tally
            if matched_genomes._phamerator_features_unmatched_in_ncbi_tally > 0:
                self._ph_ncbi_genomes_with_unmatched_phamerator_features_tally += 1
                self._ph_ncbi_unmatched_phamerator_features_tally += matched_genomes._phamerator_features_unmatched_in_ncbi_tally
            if matched_genomes._ncbi_features_unmatched_in_phamerator_tally > 0:
                self._ph_ncbi_genomes_with_unmatched_ncbi_features_tally += 1
                self._ph_ncbi_unmatched_ncbi_features_tally += matched_genomes._ncbi_features_unmatched_in_phamerator_tally
            if matched_genomes._phamerator_ncbi_different_descriptions_tally > 0:
                self._ph_ncbi_genomes_with_different_descriptions_tally += 1
                self._ph_ncbi_different_descriptions_tally += matched_genomes._phamerator_ncbi_different_descriptions_tally
            if matched_genomes._phamerator_ncbi_different_translations_tally > 0:
                self._ph_ncbi_genomes_with_different_translations_tally += 1
                self._ph_ncbi_different_translation_tally += matched_genomes._phamerator_ncbi_different_translations_tally
            if matched_genomes._ph_ncbi_author_error:
                self._ph_ncbi_genomes_with_author_errors_tally += 1


            # PhagesDB-GenBank checks
            if matched_genomes._phagesdb_ncbi_sequence_mismatch:
                self._pdb_ncbi_sequence_mismatch_tally += 1
            if matched_genomes._phagesdb_ncbi_sequence_length_mismatch:
                self._pdb_ncbi_sequence_length_mismatch_tally += 1















# End of class definitions





































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
    parser.add_argument("-i", "--internal_authors", action="store_true",
        default=False, help=INTERNAL_AUTHORS_HELP)
    parser.add_argument("-e", "--external_authors", action="store_true",
        default=False, help=EXTERNAL_AUTHORS_HELP)
    parser.add_argument("-d", "--draft", action="store_true",
        default=False, help=DRAFT_HELP)
    parser.add_argument("-f", "--final", action="store_true",
        default=False, help=FINAL_HELP)
    parser.add_argument("-u", "--unknown", action="store_true",
        default=False, help=UNKNOWN_HELP)



    # Assumed command line arg structure:
    # python3 -m pdm_utils <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


def get_dbs(phagesdb, genbank):
    """Create set of databases to compare to MySQL."""
    dbs = set()
    if phagesdb == True:
        dbs.add("phagesdb")
    if genbank == True:
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

    # Setup output
    output_folder = basic.set_path(args.output_folder, kind="dir",
                                        expect=True)
    working_dir = pathlib.Path(RESULTS_FOLDER)
    working_path = basic.make_new_dir(output_folder, working_dir,
                                      attempt=50)
    if working_path is None:
        print(f"Invalid working directory '{working_dir}'")
        sys.exit(1)

    # Verify database connection and schema compatibility.
    print("Connecting to the MySQL database...")
    alchemist = AlchemyHandler(database=database)
    alchemist.connect(pipeline=True)
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the compare pipeline")

    # Determine which database should be compared in addition to
    # MySQL: PhagesDB, GenBank.
    valid_dbs = get_dbs(args.phagesdb, args.genbank)
    selected_dbs = "Databases compared to MySQL: " + ", ".join(valid_dbs)

    # Determine which type of genomes should be checked based on
    # who annotated the genome: Hatfull or Gbk authors
    valid_authors = get_authors(args.internal_authors, args.external_authors)
    author_str = []
    for i in valid_authors:
        author_str.append(str(i))
    selected_authors =  ("Genomes with following AnnotationAuthor "
                         "will be compared: " + ", ".join(author_str))

    # Determine which types of genomes should be checked based on
    # the annotation_status of the annotations: draft, final, unknown
    # ncbi_credentials_file
    valid_status = get_status(args.draft, args.final, args.unknown)
    selected_status = ("Genomes with following Annotation Status will "
                       "be compared: " + ", ".join(valid_status))

    # Get data from the MySQL database.
    version = str(mysqldb_basic.scalar(engine, VERSION_QUERY))
    mysql_gnm_data_dict = mysqldb_basic.query_dict_list(engine, PHAGE_QUERY)
    mysql_gene_data_dict = mysqldb_basic.query_dict_list(engine, GENE_QUERY)

    mysql_tup = create_mysql_gnms(mysql_gnm_data_dict, valid_status, valid_authors)
    mysql_gnms_dict = mysql_tup[0]
    mysql_name_duplicates = mysql_tup[2]
    mysql_accessions = mysql_tup[3]
    mysql_acc_duplicates = mysql_tup[4]

    if save_records == True:
        save_mysql_gnms(mysql_gnms_dict, working_path, MYSQL_OUTPUT)

    # PhagesDB relies on the phageName, and not the phageID.
    # But the MySQL database does not require phageName values to be unique.
    # Check if there are any phageName duplications.
    # If there are, they will not be able to be compared to PhagesDB data.
    # Output duplicate phage search names to file
    if len(mysql_name_duplicates) > 0:
        print("Warning: There are duplicate phage names in the MySQL database.")
        print("Some MySQL genomes will not be able to be matched to PhagesDB.")
        output_to_file(list(mysql_name_duplicates),
                        working_path,
                        DUPLICATE_MYSQL_NAMES,
                        selected_status,
                        database, version,
                        selected_authors)
        input('Press ENTER to proceed')

    # Accessions aren't required to be unique in the MySQL
    # database (but they should be), so there could be duplicates
    # Output duplicate names to file
    if len(mysql_acc_duplicates) > 0:
        print("Warning: There are duplicate accessions in the MySQL database.")
        print("Some MySQL genomes will not be able to be matched to GenBank records.")
        output_to_file(list(mysql_acc_duplicates),
                        working_path,
                        DUPLICATE_MYSQL_ACC,
                        selected_status,
                        database, version,
                        selected_authors)
        input('Press ENTER to proceed')

    mysql_genes = get_gene_obs_list(mysql_gene_data_dict)
    match_gnm_genes(mysql_gnms_dict, mysql_genes)

    #Now retrieve all PhagesDB data
    if "phagesdb" in valid_dbs:
        pdb_tup = get_phagesdb_data()
        pdb_genome_dict = pdb_tup[0]
        pdb_search_name_set = pdb_tup[1]
        pdb_search_name_duplicate_set = pdb_tup[2]

        if save_records == True:
            save_phagesdb_genomes(pdb_genome_dict, working_path, PHAGESDB_OUTPUT)
    else:
        pdb_genome_dict = {}
        pdb_search_name_set = set()
        pdb_search_name_duplicate_set = set()


        # PhagesDB phage names should be unique, but just make sure after
        # they are converted to a search name.
        if len(pdb_search_name_duplicate_set) > 0:

            print("Warning: There are duplicate phage search names in phagesdb.")
            print("Some phagesdb genomes will not be able to be "
                  "matched to genomes in MySQL.")
            output_to_file(list(pdb_search_name_duplicate_set),
                            working_path,
                            DUPLICATE_PDB_NAMES,
                            selected_status,
                            database, version,
                            selected_authors)
            input('Press ENTER to proceed')


    #Retrieve and parse GenBank records if selected by user
    ncbi_genome_dict = {} #Key = accession; #Value = genome data
    if "genbank" in valid_dbs:

        if save_records == True:
            ncbi_output_path = pathlib.Path(working_path, GENBANK_OUTPUT)
            ncbi_output_path.mkdir()

        ncbi_cred_dict = ncbi.get_ncbi_creds(ncbi_credentials_file)
        gbk_tup = get_genbank_data(ncbi_cred_dict, mysql_accessions)
        retrieved_record_list = gbk_tup[0]
        retrieval_error_list = gbk_tup[1]

        #Report the accessions that could not be retrieved.
        output_failed_accession(retrieval_error_list, working_path,
                                FAILED_ACC_RETRIEVE, selected_status,
                                selected_dbs, database, selected_authors,
                                version)

        for retrieved_record in retrieved_record_list:
            genome_object = create_gbk_gnm(retrieved_record)
            ncbi_genome_dict[genome_object.accession] = genome_object
            if save_records == True:
                save_gbk_genome(genome_object, retrieved_record, ncbi_output_path)


    # Now that all GenBank and PhagesDB data is retrieved,
    # match up to MySQL genome data.

    match_tup = match_all_genomes(mysql_gnms_dict, pdb_genome_dict, ncbi_genome_dict,
                      pdb_search_name_duplicate_set, mysql_acc_duplicates)
    matched_genomes_list = match_tup[0]

    #Output unmatched data to file
    unmatched_rows = prepare_unmatched_output(match_tup[1], match_tup[2])
    output_unmatched(unmatched_rows, working_path,
                     FAILED_ACC_RETRIEVE, selected_status,
                     selected_dbs, database, selected_authors,
                     version)

    #Now that all genomes have been matched, iterate through each
    # matched objects and run methods to compare the genomes
    compare_data(matched_genomes_list)

    #Now that all individual matched_genome_objects have all computed attributes,
    #compute database summary
    summary_object = DatabaseSummary(matched_genomes_list)
    summary_object.compute_summary()

    #Now output results
    output_all_data(working_path, database, version, selected_authors,
                    selected_status, selected_dbs, summary_object)

    end_time = time.strftime("%x %X")
    print(f"Start time: {start_time}")
    print(f"Stop time: {end_time}")

    return


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

        # Check to see if the genome has a user-selected annotation_status and authorship
        if (dict["Status"] not in valid_status or
                dict["AnnotationAuthor"] not in valid_authors):
            continue
        else:
            gnm = create_mysql_gnm(dict)
            gnm_dict[dict["PhageID"]] = gnm

            # This keeps track of whether there are duplicate phage
            # names that will be used to match up to PhagesDB data.
            if gnm._search_name in names:
                duplicate_names.add(gnm._search_name)
            else:
                names.add(gnm._search_name)

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
    gnm = PhameratorGenome()
    gnm.set_phage_id(dict["PhageID"])
    gnm.set_phage_name(dict["Name"])
    gnm.host_genus = dict["HostGenus"]
    gnm.set_sequence(dict["Sequence"].decode("utf-8"))
    gnm.set_accession(dict["Accession"])
    gnm.set_status(dict["Status"])
    gnm.set_cluster_subcluster(dict["Cluster"])
    gnm.retrieve_record = dict["RetrieveRecord"]
    gnm.set_date_last_modified(dict["DateLastModified"])
    gnm.annotation_author = dict["AnnotationAuthor"]
    gnm.compute_nucleotide_errors(DNA_ALPHABET)
    return gnm


def save_mysql_gnms(gnm_dict, main_path, new_dir):
    """Save MySQL genome data to file."""

    output_path = pathlib.Path(main_path, new_dir)
    output_path.mkdir()

    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        name = gnm._search_name
        fasta = SeqRecord(Seq(gnm.sequence), id=name, description="")
        save_seqrecord(fasta, output_path, name, "fasta", "fasta")




def get_gene_obs_list(list_of_dicts):
    """Convert list of tuples of data to list of Cds objects."""

    print(f"\n\nPreparing {len(list_of_dicts)} gene data set(s) "
          "from the MySQL database...")
    lst = []
    for dict in list_of_dicts:
        gene = create_mysql_gene(dict)
        lst.append(gene)
    return lst


def create_mysql_gene(dict):
    """Create MySQL CDS feature object."""
    gene = MysqlCdsFeature()
    gene.set_phage_id(dict["PhageID"])
    gene.id = dict["GeneID"]
    gene.name = dict["Name"]
    gene.type = "CDS"
    gene.start = dict["Start"]
    gene.stop = dict["Stop"]
    gene.set_strand(dict["Orientation"])
    gene.set_translation(dict["Translation"])
    tup2 = retrieve_description(dict["Notes"].decode("utf-8"))
    gene.set_notes(tup2[0], tup2[1])
    gene.compute_amino_acid_errors(PROTEIN_ALPHABET)
    gene.set_start_end_strand_id()
    gene.compute_boundary_error()
    return gene



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
        gnm.compute_cds_feature_errors()
        gnm.compute_status_description_error()




def get_phagesdb_data():
    """Retrieve data from PhagesDB."""
    #Data for each phage is stored in a dictionary per phage,
    # and all dictionaries are stored in a list under "results"

    print('\n\nRetrieving data from PhagesDB...')
    gnm_dict = {}
    names = set()
    dupe_names = set()

    #Retrieve a list of all sequenced phages listed on PhagesDB
    #You have to specify how many results to return at once.
    # If you set it to 1 page long and 100,000 genomes/page,
    # then this will return everything.
    url = "http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000"
    data_dict = get_json_data(url)
    for element_dict in data_dict["results"]:

        gnm = create_pdb_gnm(element_dict)
        pdb_search_name = gnm._search_name
        if pdb_search_name in names:
            dupe_names.add(pdb_search_name)
        else:
            names.add(pdb_search_name)
            gnm_dict[pdb_search_name] = gnm

    # Make sure all sequenced phage data has been retrieved
    check_pdb_retrieved(data_dict, gnm_dict)

    return gnm_dict, names, dupe_names


def check_pdb_retrieved(data_dict, gnm_dict):
    """Check if all PhagesDB data has been retrieved."""

    if (len(data_dict["results"]) != data_dict["count"] or
            len(data_dict["results"]) != len(gnm_dict)):

        print("\nUnable to retrieve all phage data from PhagesDB "
              "due to default retrieval parameters.")
        print("Update parameters in script to enable these functions.")
        input("Press ENTER to proceed")



def save_phagesdb_genomes(gnm_dict, main_path, new_dir):
    """Save PhagesDB genomes."""

    output_path = pathlib.Path(main_path, new_dir)
    output_path.mkdir()

    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        name = gnm._search_name
        seq = gnm.sequence
        fasta = SeqRecord(Seq(seq), id=name, description="")
        save_seqrecord(fasta, output_path, name, "fasta", "fasta")



def get_json_data(url):
    """Get data from PhagesDB."""
    json_data = urllib.request.urlopen(url)
    data_dict = json.loads(json_data.read())
    json_data.close()
    return data_dict



def create_pdb_gnm(dict):
    """Parse PhagesDB data to genome object."""

    gnm = PhagesdbGenome()

    #Name, Host, Accession
    gnm.set_phage_name(dict["phage_name"])
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
        seq = parse_fasta_seq(dict["fasta_file"])
        gnm.set_sequence(seq)
        gnm.compute_nucleotide_errors(DNA_ALPHABET)

    return gnm




def parse_fasta_seq(url):
    """Parse sequence from fasta file url."""
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request) as response:
        file = response.read()
        file = file.decode("utf-8")

    # All sequence rows in the fasta file may not have equal widths,
    # so some processing of the data is required
    # If split by newline, the header is retained in the first list element.
    data = file.split("\n")
    seq = ""
    for index in range(len(data)):
        # Strip off potential whitespace before appending, such as '\r'
        seq = seq + data[index].strip()
    return seq



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

    # When retrieving in batch sizes, first create the list of
    # values indicating which indices of the accession_retrieval_list
    # should be used to create each batch.
    # For instace, if there are five accessions,
    # batch size of two produces indices = 0,2,4.
    for batch_index_start in range(0,len(acc_list),batch_size):

        if batch_index_start + batch_size > len(acc_list):
            batch_index_stop = len(acc_list)
        else:
            batch_index_stop = batch_index_start + batch_size

        current_batch_size = batch_index_stop - batch_index_start
        delimiter = " | "
        esearch_term = delimiter.join(acc_list[batch_index_start:batch_index_stop])

        #Use esearch for each accession
        search_handle = Entrez.esearch(db="nucleotide", term=esearch_term,
                                       usehistory="y")
        search_record = Entrez.read(search_handle)
        search_count = int(search_record["Count"])
        search_webenv = search_record["WebEnv"]
        search_query_key = search_record["QueryKey"]

        # Keep track of the accessions that failed to be located in GenBank
        if search_count < current_batch_size:
            search_accession_failure = search_record["ErrorList"]["PhraseNotFound"]

            # Each element in this list is formatted "accession[ACCN]"
            for element in search_accession_failure:
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

    gnm = NcbiGenome()

    try:
        gnm.record_name = retrieved_record.name
    except:
        gnm.record_name = ""

    try:
        gnm.record_id = retrieved_record.id
    except:
        gnm.record_id = ""

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
        if organism.split(" ")[-1] == "Unclassified.":
            gnm.set_phage_name(organism.split(" ")[-2])
        else:
            gnm.set_phage_name(organism.split(" ")[-1])

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

    # Nucleotide sequence and errors
    gnm.set_sequence(retrieved_record.seq)
    gnm.compute_nucleotide_errors(DNA_ALPHABET)

    # Iterate through all features
    source_feature_list = []
    ncbi_cds_features = []

    for feature in retrieved_record.features:

        if feature.type != "CDS":
            # Retrieve the Source Feature info
            if feature.type == "source":
                source_feature_list.append(feature)
        else:
            gene_object = create_gbk_gene(feature)
            ncbi_cds_features.append(gene_object)

    # Set the following variables after iterating through all features
    # If there was one and only one source feature present,
    # parse certain qualifiers.
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

    gnm.set_cds_features(ncbi_cds_features)
    gnm.compute_cds_feature_errors()
    gnm.compute_ncbi_cds_feature_errors()

    return gnm


def create_gbk_gene(feature):
    """Parse data from GenBank CDS feature."""

    gene_object = NcbiCdsFeature()

    # Feature type
    gene_object.type = "CDS"

    # Locus tag
    try:
        gene_object.set_locus_tag(feature.qualifiers["locus_tag"][0])
    except:
        gene_object.set_locus_tag("")

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
        tup1 = retrieve_description(feature.qualifiers["product"][0])
        gene_object.set_product_description(tup1[0],tup1[1])
    except:
        pass
    try:
        tup2 = retrieve_description(feature.qualifiers["function"][0])
        gene_object.set_function_description(tup2[0],tup2[1])
    except:
        pass

    try:
        tup3 = retrieve_description(feature.qualifiers["note"][0])
        gene_object.set_note_description(tup3[0],tup3[1])
    except:
        pass

    # Gene number
    try:
        gene_object.set_gene_number(feature.qualifiers["gene"][0])
    except:
        pass

    # Compute other fields
    gene_object.compute_amino_acid_errors(PROTEIN_ALPHABET)
    gene_object.set_start_end_strand_id()
    gene_object.compute_boundary_error()
    gene_object.compute_description_error()
    return gene_object


def save_gbk_genome(genome_object, record, output_path):
    """Save GenBank record to file."""

    name_prefix = genome_object.organism.lower()
    if name_prefix.split(" ")[-1] == "unclassified.":
        name_prefix = name_prefix.split(" ")[-2]
    else:
        name_prefix = name_prefix.split(" ")[-1]

    file_prefix = name_prefix + "__" + genome_object.accession
    save_seqrecord(record, output_path, file_prefix, "gb", "genbank")


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
        matched_objects.set_phamerator_genome(mysql_gnm)

        # Match up PhagesDB genome
        # First try to match up the phageID, and if that doesn't work,
        # try to match up the phageName
        if mysql_gnm._search_id in pdb_gnms.keys():
            pdb_genome = pdb_gnms[mysql_gnm._search_id]
            # Make sure the pdb_genome doesn't have a search name
            # that was duplicated
            if pdb_genome._search_name in pdb_name_duplicates:
                pdb_genome = ""

        elif mysql_gnm._search_name in pdb_gnms.keys() and not \
            mysql_gnm._search_name in mysql_name_duplicates:

            pdb_genome = pdb_gnms[mysql_gnm._search_name]
            # Make sure the pdb_genome doesn't have a search name
            # that was duplicated.
            if pdb_genome._search_name in pdb_name_duplicates:
                pdb_genome = ""

        else:
            pdb_genome = ""
            mysql_unmatched_to_pdb_gnms.append(mysql_gnm)

        matched_objects.set_phagesdb_genome(pdb_genome)

        # Now match up GenBank genome
        acc = mysql_gnm.accession
        if acc != "" and acc not in mysql_acc_duplicates:
                # Retrieval of record may have failed
                try:
                    ncbi_genome = gbk_gnms[acc]
                except:
                    ncbi_genome = ""
                    mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        else:
            ncbi_genome = ""
            mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        matched_objects.set_ncbi_genome(ncbi_genome)
        matched_gnms.append(matched_objects)
        match_count += 1

    return matched_gnms, mysql_unmatched_to_pdb_gnms, mysql_unmatched_to_gbk_gnms



def compare_data(matched_gnms):
    """Compare all matched data."""

    print("Comparing matched genomes and identifying inconsistencies...")
    count = 1
    total = len(matched_gnms)
    for matched_gnm in matched_gnms:
        print(f"Comparing matched genome set {count} of {total}")
        matched_gnm.compare_phamerator_ncbi_genomes()
        matched_gnm.compare_phamerator_phagesdb_genomes()
        matched_gnm.compare_phagesdb_ncbi_genomes()
        matched_gnm.compute_total_genome_errors()
        count += 1










def output_all_data(output_path, database, version, selected_authors,
                    selected_status, selected_dbs, summary_object):
    """Output all analysis results."""
    print("Outputting results to file...")

    gene_file_path = pathlib.Path(output_path, GENE_OUTPUT)
    gene_report_fh = gene_file_path.open("w")
    gene_report_writer = csv.writer(gene_report_fh)
    gene_report_writer.writerow([FIRST_HEADER])
    gene_report_writer.writerow([database_header(database, version)])
    gene_report_writer.writerow([selected_authors])
    gene_report_writer.writerow([selected_status])
    gene_report_writer.writerow([selected_dbs])

    gene_headers = create_gene_headers()
    gene_report_writer.writerow(gene_headers)


    #Create database summary output file
    summary_output = create_summary_data(summary_object)
    output_summary_data(summary_output, output_path, SUMMARY_OUTPUT, selected_status,
                        selected_dbs, database, selected_authors, version)

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
        ph_genome = matched_genomes.m_genome
        all_ftrs = get_all_features(matched_genomes)
        for mixed_ftr in all_ftrs:
            ftr_data = create_feature_data(ph_genome, mixed_ftr)
            gene_report_writer.writerow(ftr_data)

    genome_output_list = create_genome_headers()
    genome_output_list.extend(genomes_data)
    output_genome_results(genome_output_list, output_path, GENOME_OUTPUT, selected_status,
                          selected_dbs, database, selected_authors, version)
    gene_report_fh.close()
    return


def get_all_features(matched_genomes):
    """Create list of all features."""

    perfect_match = matched_genomes._phamerator_ncbi_perfect_matched_features
    imperfectly_match = matched_genomes._phamerator_ncbi_imperfect_matched_features
    mysql_unmatch = matched_genomes._phamerator_features_unmatched_in_ncbi
    gbk_unmatch = matched_genomes._ncbi_features_unmatched_in_phamerator

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
        "ph_search_name",
        "total_errors",

        # MySQL
        # General gene data
        "ph_phage_id",
        "ph_search_id",
        "ph_type_id",
        "ph_gene_id",
        "ph_gene_name",
        "ph_left_boundary",
        "ph_right_boundary",
        "ph_strand",
        "ph_translation",
        "ph_translation_length",
        "ph_gene_notes",

        # Gene data checks
        "ph_translation_error",
        "ph_gene_coords_error",

        # GenBank
        # General gene data
        "ncbi_locus_tag",
        "ncbi_gene_number",
        "ncbi_type_id",
        "ncbi_left_boundary",
        "ncbi_right_boundary",
        "ncbi_strand",
        "ncbi_translation",
        "ncbi_translation_length",
        "ncbi_product_description",
        "ncbi_function_description",
        "ncbi_note_description",

        # Gene data checks
        "ncbi_translation_error",
        "ncbi_gene_coords_error",
        "ncbi_missing_locus_tag",
        "ncbi_locus_tag_typo",
        "ncbi_description_field_error",

        # MySQL-GenBank checks
        "ph_ncbi_unmatched_error",
        "ph_ncbi_description_error",
        "ph_ncbi_start_coordinate_error",
        "ph_ncbi_translation_error"
        ]
    return headers


def create_genome_headers():
    """Create list of column headers."""
    headers = [
        "ph_phage_id",
        "contains_errors",

        # MySQL
        # General genome data
        "ph_phage_name",
        "ph_search_id",
        "ph_search_name",
        "ph_status",
        "ph_cluster_subcluster",
        "ph_host",
        "ph_accession",
        "ph_dna_seq_length",
        "ph_gene_tally",
        "ph_description_tally",
        "ph_ncbi_update_flag",
        "ph_date_last_modified",
        "ph_annotation_author",

        # Genome data checks
        "ph_dna_seq_error",
        "ph_gene_translation_error_tally",
        "ph_gene_coords_error_tally",
        "ph_status_description_error",
        "ph_status_accession_error",

        # PhagesDB
        # General genome data
        "pdb_phage_name",
        "pdb_search_name",
        "pdb_cluster",
        "pdb_subcluster",
        "pdb_host",
        "pdb_accession",
        "pdb_dna_seq_length",

        # Genome data checks
        "pdb_dna_seq_error",

        # GenBank
        # General genome data
        "ncbi_phage_name",
        "ncbi_search_name",
        "ncbi_record_id",
        "ncbi_record_name",
        "ncbi_record_accession",
        "ncbi_record_definition",
        "ncbi_record_source",
        "ncbi_record_organism",
        "ncbi_source_feature_organism",
        "ncbi_source_feature_host",
        "ncbi_source_feature_lab_host",
        "ncbi_authors",
        "ncbi_dna_seq_length",
        "ncbi_gene_tally",

        # Genome data checks
        "ncbi_dna_seq_error",
        "ncbi_gene_translation_error_tally",
        "ncbi_gene_coords_error_tally",
        "ncbi_gene_product_tally",
        "ncbi_gene_function_tally",
        "ncbi_gene_note_tally",
        "ncbi_missing_locus_tag_tally",
        "ncbi_locus_tag_typo_tally",
        "ncbi_description_field_error_tally",

        # MySQL-PhagesDB
        "ph_pdb_dna_seq_error",
        "ph_pdb_dna_seq_length_error",
        "ph_pdb_cluster_subcluster_error",
        "ph_pdb_accession_error",
        "ph_pdb_host_error",

        # MySQL-GenBank
        "ph_ncbi_dna_seq_error",
        "ph_ncbi_dna_seq_length_error",
        "ph_ncbi_record_header_name_error",
        "ph_ncbi_record_header_host_error",

        # Author error is dependent on MySQL genome annotation author and
        # GenBank list of authors, so this metric should be reported with
        # the other ph_ncbi error tallies.
        "ph_ncbi_author_error",
        "ph_ncbi_perfectly_matched_gene_tally",
        "ph_ncbi_imperfectly_matched_gene_tally",
        "ph_ncbi_unmatched_phamerator_gene_tally",
        "ph_ncbi_unmatched_ncbi_gene_tally",
        "ph_ncbi_gene_description_error_tally",
        "ph_ncbi_perfectly_matched_gene_translation_error_tally",

        # Number of genes with errors is computed slightly differently
        # depending on whethere there are matching MySQL and GenBank genomes.
        # Therefore,this metric should be reported with the
        # other ph_ncbi error tallies even if there is no matching GenBank genome.
        "ph_ncbi_genes_with_errors_tally",

        # PhagesDB-GenBank
        "pdb_ncbi_dna_seq_error",
        "pdb_ncbi_dna_seq_length_error"
        ]
    return [headers]






def create_summary_fields():
    """Create summary row ids."""
    fields = [
        "",
        "Genome summary:",

        # Column header
        "Database comparison metric",

        # Database summaries
        "ph_total_genomes_analyzed",
        "ph_genomes_unmatched_to_pdb_tally",
        "ph_genomes_unmatched_to_ncbi_tally",
        "total_genomes_with_errors",


        # MySQL data
        # General genome data
        "ph_ncbi_update_flag_tally",

        # Genome data checks
        "ph_genomes_with_nucleotide_errors_tally",
        "ph_genomes_with_translation_errors_tally",
        "ph_genomes_with_boundary_errors_tally",
        "ph_genomes_with_status_accession_error_tally",
        "ph_genomes_with_status_description_error_tally",

        # PhagesDB data
        # Genome data checks
        "pdb_genomes_with_nucleotide_errors_tally",

        # GenBank data
        # Genome data checks
        "ncbi_genomes_with_description_field_errors_tally",
        "ncbi_genomes_with_nucleotide_errors_tally",
        "ncbi_genomes_with_translation_errors_tally",
        "ncbi_genomes_with_boundary_errors_tally",
        "ncbi_genomes_with_missing_locus_tags_tally",
        "ncbi_genomes_with_locus_tag_typos_tally",

        # MySQL-PhagesDB checks
        "ph_pdb_sequence_mismatch_tally",
        "ph_pdb_sequence_length_mismatch_tally",
        "ph_pdb_cluster_subcluster_mismatch_tally",
        "ph_pdb_accession_mismatch_tally",
        "ph_pdb_host_mismatch_tally",

        # MySQL-GenBank checks
        "ph_ncbi_sequence_mismatch_tally",
        "ph_ncbi_sequence_length_mismatch_tally",
        "ph_ncbi_record_header_phage_mismatch_tally",
        "ph_ncbi_record_header_host_mismatch_tally",
        "ph_ncbi_genomes_with_author_errors_tally",
        "ph_ncbi_genomes_with_imperfectly_matched_features_tally",
        "ph_ncbi_genomes_with_unmatched_phamerator_features_tally",
        "ph_ncbi_genomes_with_unmatched_ncbi_features_tally",
        "ph_ncbi_genomes_with_different_descriptions_tally",
        "ph_ncbi_genomes_with_different_translations_tally",

        # PhagesDB-GenBank checks
        "pdb_ncbi_sequence_mismatch_tally",
        "pdb_ncbi_sequence_length_mismatch_tally",

        # Separate all checks that tally all genes
        "",
        "Gene summary:",

        # Column header
        "Database comparison metric",

        # MySQL feature
        # Gene data checks
        "ph_translation_errors_tally",
        "ph_boundary_errors_tally",

        # GenBank feature
        # Gene data checks
        "ncbi_translation_errors_tally",
        "ncbi_boundary_errors_tally",
        "ncbi_missing_locus_tags_tally",
        "ncbi_locus_tag_typos_tally",
        "ncbi_description_field_errors_tally",

        # MySQL-GenBank checks
        "ph_ncbi_different_descriptions_tally",
        "ph_ncbi_different_start_sites_tally",
        "ph_ncbi_different_translation_tally",
        "ph_ncbi_unmatched_phamerator_features_tally",
        "ph_ncbi_unmatched_ncbi_features_tally"
        ]

    return fields


def create_summary_data(summary_object):
    """Create summary of all results."""
    lst1 = []
    lst1.append("")
    lst1.append("")
    lst1.append("tally") # Column header

    # First output database summary data
    lst1.append(summary_object._ph_total_genomes_analyzed)
    lst1.append(summary_object._ph_genomes_unmatched_to_pdb_tally)
    lst1.append(summary_object._ph_genomes_unmatched_to_ncbi_tally)
    lst1.append(summary_object._total_genomes_with_errors)

    # MySQL data
    # General genome data
    lst1.append(summary_object._ph_ncbi_update_flag_tally)

    # Genome data checks
    lst1.append(summary_object._ph_genomes_with_nucleotide_errors_tally)
    lst1.append(summary_object._ph_genomes_with_translation_errors_tally)
    lst1.append(summary_object._ph_genomes_with_boundary_errors_tally)
    lst1.append(summary_object._ph_genomes_with_status_accession_error_tally)
    lst1.append(summary_object._ph_genomes_with_status_description_error_tally)

    # PhagesDB data
    # Genome data checks
    lst1.append(summary_object._pdb_genomes_with_nucleotide_errors_tally)

    # GenBank data
    # Genome data checks
    lst1.append(summary_object._ncbi_genomes_with_description_field_errors_tally)
    lst1.append(summary_object._ncbi_genomes_with_nucleotide_errors_tally)
    lst1.append(summary_object._ncbi_genomes_with_translation_errors_tally)
    lst1.append(summary_object._ncbi_genomes_with_boundary_errors_tally)
    lst1.append(summary_object._ncbi_genomes_with_missing_locus_tags_tally)
    lst1.append(summary_object._ncbi_genomes_with_locus_tag_typos_tally)


    # MySQL-PhagesDB checks
    lst1.append(summary_object._ph_pdb_sequence_mismatch_tally)
    lst1.append(summary_object._ph_pdb_sequence_length_mismatch_tally)
    lst1.append(summary_object._ph_pdb_cluster_subcluster_mismatch_tally)
    lst1.append(summary_object._ph_pdb_accession_mismatch_tally)
    lst1.append(summary_object._ph_pdb_host_mismatch_tally)

    # MySQL-GenBank checks
    lst1.append(summary_object._ph_ncbi_sequence_mismatch_tally)
    lst1.append(summary_object._ph_ncbi_sequence_length_mismatch_tally)
    lst1.append(summary_object._ph_ncbi_record_header_phage_mismatch_tally)
    lst1.append(summary_object._ph_ncbi_record_header_host_mismatch_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_author_errors_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_imperfectly_matched_features_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_unmatched_phamerator_features_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_unmatched_ncbi_features_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_different_descriptions_tally)
    lst1.append(summary_object._ph_ncbi_genomes_with_different_translations_tally)

    # PhagesDB-GenBank checks
    lst1.append(summary_object._pdb_ncbi_sequence_mismatch_tally)
    lst1.append(summary_object._pdb_ncbi_sequence_length_mismatch_tally)


    # Gene summaries
    lst1.append("")
    lst1.append("")
    lst1.append("tally") # Column header

    # MySQL feature
    # Gene data checks
    lst1.append(summary_object._ph_translation_errors_tally)
    lst1.append(summary_object._ph_boundary_errors_tally)

    # GenBank feature
    # Gene data checks
    lst1.append(summary_object._ncbi_translation_errors_tally)
    lst1.append(summary_object._ncbi_boundary_errors_tally)
    lst1.append(summary_object._ncbi_missing_locus_tags_tally)
    lst1.append(summary_object._ncbi_locus_tag_typos_tally)
    lst1.append(summary_object._ncbi_description_field_errors_tally)

    # MySQL-GenBank checks
    lst1.append(summary_object._ph_ncbi_different_descriptions_tally)
    lst1.append(summary_object._ph_ncbi_different_start_sites_tally)
    lst1.append(summary_object._ph_ncbi_different_translation_tally)
    lst1.append(summary_object._ph_ncbi_unmatched_phamerator_features_tally)
    lst1.append(summary_object._ph_ncbi_unmatched_ncbi_features_tally)

    fields = create_summary_fields()
    lst2 = []
    if len(fields) == len(lst1):
        for i in range(len(fields)):
            lst2.append([fields[i],lst1[i]])
    else:
        lst2.append([("Different number of database summary metrics "
                      "than database summary data fields.")])
        lst2.append(["Unable to output database summary."])
    return lst2




def create_feature_data(mysql_genome, mixed_ftr):
    """Create feature data to output."""

    lst = []

    # Gene summaries
    # Add MySQL genome search name to each gene row regardless of the type of CDS data (matched or unmatched)
    lst.append(mysql_genome._search_name) # matched MySQL search name
    lst.append(mixed_ftr._total_errors) # total # of errors for this gene

    if isinstance(mixed_ftr,MatchedCdsFeatures):
        mysql_ftr = mixed_ftr._m_feature
        gbk_ftr = mixed_ftr._g_ftr
    else:
        if isinstance(mixed_ftr,MysqlCdsFeature):
            mysql_ftr = mixed_ftr
            gbk_ftr = ""
        elif isinstance(mixed_ftr,NcbiCdsFeature):
            mysql_ftr = ""
            gbk_ftr = mixed_ftr
        else:
            mysql_ftr = ""
            gbk_ftr = ""

    # MySQL feature
    if isinstance(mysql_ftr,MysqlCdsFeature):

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
    if isinstance(gbk_ftr,NcbiCdsFeature):

        # General gene data
        lst.append(gbk_ftr.locus_tag)
        lst.append(gbk_ftr.gene_number)
        lst.append(gbk_ftr.type)
        lst.append(gbk_ftr.start)
        lst.append(gbk_ftr.stop)
        lst.append(gbk_ftr.orientation)
        lst.append(gbk_ftr.translation)
        lst.append(gbk_ftr.translation_length)
        lst.append(gbk_ftr.product_description)
        lst.append(gbk_ftr.function_description)
        lst.append(gbk_ftr.note_description)

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

        # If this is a matched cds feature, both MySQL and
        # GenBank features should have identical unmatched_error value.
        lst.append(mixed_ftr._m_feature._unmatched_error)

        lst.append(mixed_ftr._phamerator_ncbi_different_descriptions)
        lst.append(mixed_ftr._phamerator_ncbi_different_start_sites)
        lst.append(mixed_ftr._phamerator_ncbi_different_translations)
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
    lst.append(mysql_gnm._search_id)
    lst.append(mysql_gnm._search_name)
    lst.append(mysql_gnm.annotation_status)
    lst.append(mysql_gnm.cluster_subcluster)
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
    if isinstance(pdb_genome, PhagesdbGenome):

        # General genome data
        lst.append(pdb_genome.name)
        lst.append(pdb_genome._search_name)
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
    if isinstance(gbk_gnm, NcbiGenome):

        # General genome data
        lst.append(gbk_gnm.name)
        lst.append(gbk_gnm._search_name)
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
    if isinstance(pdb_genome, PhagesdbGenome):
        lst.append(matched_genomes._phamerator_phagesdb_sequence_mismatch)
        lst.append(matched_genomes._phamerator_phagesdb_sequence_length_mismatch)
        lst.append(matched_genomes._phamerator_phagesdb_cluster_subcluster_mismatch)
        lst.append(matched_genomes._phamerator_phagesdb_accession_mismatch)
        lst.append(matched_genomes._phamerator_phagesdb_host_mismatch)
    else:
        lst.extend(["","","","",""])


    # MySQL-GenBank checks
    if isinstance(gbk_gnm, NcbiGenome):
        lst.append(matched_genomes._phamerator_ncbi_sequence_mismatch)
        lst.append(matched_genomes._phamerator_ncbi_sequence_length_mismatch)
        lst.append(matched_genomes._ncbi_record_header_fields_phage_name_mismatch)
        lst.append(matched_genomes._ncbi_host_mismatch)
        lst.append(matched_genomes._ph_ncbi_author_error)
        lst.append(matched_genomes._phamerator_ncbi_perfect_matched_features_tally)
        lst.append(matched_genomes._phamerator_ncbi_imperfect_matched_features_tally)
        lst.append(matched_genomes._phamerator_features_unmatched_in_ncbi_tally)
        lst.append(matched_genomes._ncbi_features_unmatched_in_phamerator_tally)
        lst.append(matched_genomes._phamerator_ncbi_different_descriptions_tally)
        lst.append(matched_genomes._phamerator_ncbi_different_translations_tally)
    else:
        lst.extend(["","","","","","","","","","",""])

    # Number of genes with errors
    lst.append(matched_genomes._total_number_genes_with_errors)

    # Output PhagesDB-GenBank checks
    if isinstance(pdb_genome, PhagesdbGenome) and isinstance(gbk_gnm, NcbiGenome):
        lst.append(matched_genomes._phagesdb_ncbi_sequence_mismatch)
        lst.append(matched_genomes._phagesdb_ncbi_sequence_length_mismatch)
    else:
        lst.extend(["",""])
    return lst




if __name__ == "__main__":
    main(sys.argv.insert(0, "empty"))
