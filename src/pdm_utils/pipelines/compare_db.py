"""Pipeline to compare data between MySQL, PhagesDB, and GenBank databases."""

# TODO this pipeline is not fully integrated into the pdm_utils package, and
# several parts need to be improved.

# 1. It originally utilized its own Cds and Genome classes. It now relies on
# the base pdm_utils Cds and Genome classes, but then requires additional
# attributes and methods, which are added in this module. The pipeline keeps
# track of errors using booleans and tallies, in contrast to the import pipeline
# which keeps track of errors using Evaluation objects. Ultimately, this pipeline
# should probably realy on Evaluation objects.

# 2. It also utilizes GenomeTriad, CdsPair, and DbCompareSummary classes
# which need to be substantially refactored and generalized.

# 3. In order to compare PhagesDB to MySQL, it needs to retrieve genome
# sequences from PhagesDB. This requires retrieving and parsing thousands of
# fasta files, which is very slow. It does this for all genomes on PhagesDB,
# even if they are not matched to MySQL genomes, so it is a big waste of time.

# 4. The data consistency checks are not comprehensive. Many more checks
# should be added.

# 5. None of the functions in this pipeline are tested.

# Note this script compares and matches data from PhagesDB, GenBank, and MySQL.
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

from pdm_utils.classes import cds
from pdm_utils.classes import genome
from pdm_utils.classes import cdspair
from pdm_utils.classes import dbcomparesummary
from pdm_utils.classes import genometriad
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import flat_files
from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import ncbi
from pdm_utils.functions import parsing
from pdm_utils.functions import phagesdb
from pdm_utils.functions import querying

DEFAULT_OUTPUT_FOLDER = os.getcwd()

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
              "CONVERT(Translation USING utf8) as Translation, "
              "Notes from gene")
VERSION_QUERY = "SELECT Version FROM version"

TARGET_TABLE = "phage"

# Identifiers for object types based on their database source.
GNM_MYSQL = "mysql"
GNM_PDB = "phagesdb"
GNM_GBK = "genbank"
CDS_MYSQL = "cds_mysql"
CDS_GBK = "cds_genbank"
CDSPAIR_MYSQL_GBK = "mysql_genbank"



# TODO these Genome and Cds class modifications are a temporary
# solution to be able to integrate the compare pipeline into pdm_utils.
# These modifications should ultimately be
# integrated directly into the base pdm_utils classes or
# the compare pipeline should be revamped to not need these extra
# attributes and methods.

### Genome class modifications

# TODO refactor and test.
def modify_genome_class(GenomeClass):
    """Add new attributes and methods to Genome class."""

    # New attributes
    setattr(GenomeClass, "record_name", "")
    setattr(GenomeClass, "record_id", "")
    setattr(GenomeClass, "source_feature_organism", "")
    setattr(GenomeClass, "source_feature_host", "")
    setattr(GenomeClass, "source_feature_lab_host", "")
    setattr(GenomeClass, "_missing_locus_tags_tally", 0)
    setattr(GenomeClass, "_locus_tag_typos_tally", 0)
    setattr(GenomeClass, "_description_field_error_tally", 0)
    setattr(GenomeClass, "_status_accession_error", False)
    setattr(GenomeClass, "_status_description_error", False)
    setattr(GenomeClass, "_genes_with_errors_tally", 0)
    setattr(GenomeClass, "_nucleotide_errors", False)
    setattr(GenomeClass, "_cds_features_with_translation_error_tally", 0)
    setattr(GenomeClass, "_cds_features_boundary_error_tally", 0)

    # New methods
    setattr(GenomeClass, "compute_nucleotide_errors",
            compute_nucleotide_errors)
    setattr(GenomeClass, "compute_cds_feature_errors",
            compute_cds_feature_errors)
    setattr(GenomeClass, "check_status_accession",
            check_status_accession)
    setattr(GenomeClass, "compute_status_description_error",
            compute_status_description_error)
    setattr(GenomeClass, "compute_genes_with_errors_tally",
            compute_genes_with_errors_tally)
    setattr(GenomeClass, "compute_gbk_cds_feature_errors",
            compute_gbk_cds_feature_errors)


# Genome error checks

# TODO refactor and test.
def compute_nucleotide_errors(self, dna_alphabet):
    nucleotide_set = set(self.seq)
    nucleotide_error_set = nucleotide_set - dna_alphabet
    if len(nucleotide_error_set) > 0:
        self._nucleotide_errors = True

# TODO refactor and test.
def compute_cds_feature_errors(self):
    for cds_feature in self.cds_features:
        if cds_feature._amino_acid_errors:
            self._cds_features_with_translation_error_tally += 1
        if cds_feature._boundary_error:
            self._cds_features_boundary_error_tally += 1

# TODO refactor and test.
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

# TODO refactor and test.
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


# TODO refactor and test.
def compute_genes_with_errors_tally(self):
    #Even though this method iterates through the CDS features
    # like the compute_status_description_error does,
    #it has to be kept separate, since you need to wait to run
    # this method after all genome and gene matching is completed.
    for feature in self.cds_features:
        # Need to first compute the number of errors per gene
        feature.check_for_errors()
        if feature._total_errors > 0:
            self._genes_with_errors_tally += 1

# TODO refactor and test.
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


### Cds class modification

# TODO refactor and test.
def modify_cds_class(CdsClass):
    """Add new attributes and methods to Cds class."""

    setattr(CdsClass, "_search_genome_id", "")
    setattr(CdsClass, "_start_end_strand_id", "")
    setattr(CdsClass, "_end_strand_id", "")
    setattr(CdsClass, "_locus_tag_missing", False)
    setattr(CdsClass, "_locus_tag_typo", False)
    setattr(CdsClass, "_description_field_error", False)
    setattr(CdsClass, "_amino_acid_errors", False)
    setattr(CdsClass, "_boundary_error", False)
    setattr(CdsClass, "_unmatched_error", False)
    setattr(CdsClass, "_total_errors", 0)

    # New methods
    setattr(CdsClass, "set_start_end_strand_id", set_start_end_strand_id)
    setattr(CdsClass, "set_search_genome_id", set_search_genome_id)
    setattr(CdsClass, "check_locus_tag", check_locus_tag)
    setattr(CdsClass, "compute_amino_acid_errors", compute_amino_acid_errors)
    setattr(CdsClass, "compute_boundary_error", compute_boundary_error)
    setattr(CdsClass, "set_locus_tag_typo", set_locus_tag_typo)
    setattr(CdsClass, "compute_description_error", compute_description_error)
    setattr(CdsClass, "check_for_errors", check_for_errors)

# TODO refactor and test.
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

# TODO refactor and test.
def set_search_genome_id(self):
    self._search_genome_id = basic.edit_suffix(self.genome_id, "remove").lower()


# TODO refactor and test.
def check_locus_tag(self):
    if self.locus_tag == "":
        self._locus_tag_missing = True

# TODO refactor and test.
def compute_amino_acid_errors(self, protein_alphabet):
    amino_acid_set = set(self.translation)
    amino_acid_error_set = amino_acid_set - protein_alphabet
    if len(amino_acid_error_set) > 0:
        self._amino_acid_errors = True

# TODO refactor and test.
def compute_boundary_error(self):
    # Check if start and end coordinates are fuzzy
    if not (str(self.start).isdigit() and str(self.stop).isdigit()):
        self._boundary_error = True

# TODO refactor and test.
def set_locus_tag_typo(self):
    self._locus_tag_typo = True

# TODO refactor and test.
def compute_description_error(self):
    # If the product description is empty or generic,
    # and the function or note descriptions are not, there is an error.
    if (self.product == "" and self.function != "" or self.note != ""):
        self._description_field_error = True

# TODO refactor and test.
def check_for_errors(self):
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




### Pipeline begins below

# TODO refactor and test.
def main(unparsed_args_list):
    """Run compare pipeline."""
    modify_genome_class(genome.Genome)
    modify_cds_class(cds.Cds)

    start_time = time.strftime("%x %X")
    args = parse_args(unparsed_args_list)
    database = args.database
    save_records = args.save_records
    interactive = args.interactive

    # Create config object with data obtained from file and/or defaults.
    config = configfile.build_complete_config(args.config_file)
    mysql_creds = config["mysql"]
    ncbi_creds = config["ncbi"]

    # Filters input: phage.Status=draft AND phage.HostGenus=Mycobacterium
    # Args structure: [['phage.Status=draft'], ['phage.HostGenus=Mycobacterium']]
    filters = args.filters

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
    alchemist = AlchemyHandler(database=database,
                               username=mysql_creds["user"],
                               password=mysql_creds["password"])
    alchemist.connect(pipeline=True)
    alchemist.build_metadata()
    engine = alchemist.engine
    mysqldb.check_schema_compatibility(engine, "the compare pipeline")

    # Get PhageIDs, which also validates filters.
    phage_ids = get_ids(alchemist, filters, TARGET_TABLE)

    # Gather user-selected databases.
    valid_dbs = get_dbs(args.phagesdb, args.genbank)

    # Get version from the MySQL database.
    version = str(mysqldb_basic.scalar(engine, VERSION_QUERY))

    # Record the summary of settings for the comparison.
    record_compare_settings(working_path, database, version, filters, valid_dbs)

    # Now proceed with getting all genome data.
    pmd_tup = process_mysql_data(working_path, engine, phage_ids,
                                 interactive, save_records)
    mysql_genome_dict = pmd_tup[0]
    mysql_accessions = pmd_tup[1]
    mysql_acc_duplicates = pmd_tup[2]

    # Close connections since no need to interact with MySQL anymore.
    engine.dispose()

    #Now retrieve all PhagesDB data
    if "phagesdb" in valid_dbs:
        ppd_tup = process_phagesdb_data(working_path, interactive, save_records)
        pdb_genome_dict = ppd_tup[0]
        pdb_name_duplicates = ppd_tup[1]

    else:
        pdb_genome_dict = {}
        pdb_name_duplicates = set()

    # Retrieve and parse GenBank records if selected by user
    if "genbank" in valid_dbs:
        gbk_genome_dict = process_gbk_data(working_path, ncbi_creds,
                                mysql_accessions, interactive, save_records)
        # gbk_genome_dict: Key = accession; #Value = genome data
    else:
        gbk_genome_dict = {}

    # Now that all GenBank and PhagesDB data is retrieved,
    # match up to MySQL genome data.
    match_tup = match_all_genomes(mysql_genome_dict, pdb_genome_dict,
                      gbk_genome_dict, pdb_name_duplicates,
                      mysql_acc_duplicates)
    matched_genomes_list = match_tup[0]

    # Output unmatched data to file
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


# TODO refactor and test.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for comparing databases."""

    compare_help = (
        "Pipeline to compare MySQL, PhagesDB, and "
        "GenBank databases for inconsistencies.")
    database_help = "Name of the MySQL database from which to compare data."
    output_folder_help = "Path to the folder to store results."
    phagesdb_help = "Indicates that PhagesDB data should be compared."
    genbank_help = "Indicates that GenBank data should be compared."
    save_records_help = (
        "Indicates that all genomes compared will be saved to file.")
    filters_help = (
        "Indicates which genomes to include/exclude from the comparison, "
        "with each conditional formatted as 'table.Field=value'.")
    interactive_help = (
        "Indicates whether evaluation is paused when errors are encountered.")
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=compare_help)
    parser.add_argument("database", type=str, help=database_help)
    parser.add_argument("-o", "--output_folder", type=pathlib.Path,
                        default=pathlib.Path(DEFAULT_OUTPUT_FOLDER),
                        help=output_folder_help)
    parser.add_argument("-p", "--phagesdb", action="store_true",
        default=False, help=phagesdb_help)
    parser.add_argument("-g", "--genbank", action="store_true",
        default=False, help=genbank_help)
    parser.add_argument("-f", "--filters", nargs="?",
                        type=parsing.parse_cmd_string, help=filters_help,
                        default=[])
    parser.add_argument("-s", "--save_records", action="store_true",
        default=False, help=save_records_help)
    parser.add_argument("-i", "--interactive", action="store_true",
        default=False, help=interactive_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)


    # Assumed command line arg structure:
    # python3 -m pdm_utils <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args



# TODO test.
def get_ids(alchemist, filters, target_table):
    """Get list of unique MySQL target table primary key values to evaluate."""

    primary_key = get_primary_key(alchemist.metadata, target_table)

    # Create filter object and then add command line filter strings
    db_filter = Filter(alchemist=alchemist, key=primary_key)
    db_filter.values = []

    # Attempt to add filters and exit if needed.
    add_filters(db_filter, filters)

    # Performs the query
    db_filter.update()

    # db_filter.values now contains list of PhageIDs that pass the filters.
    return db_filter.values


# TODO test.
def get_primary_key(metadata, target_table):
    """Get the primary key to the target table."""

    # Get SQLAlchemy metadata Table object
    # table_obj.primary_key.columns is a
    # SQLAlchemy ColumnCollection iterable object
    # Set primary key = 'phage.PhageID'
    table = querying.get_table(metadata, target_table)
    for column in table.primary_key.columns:
        primary_key = column
    return primary_key


# TODO test.
def add_filters(filter_obj, filters):
    """Add filters from command line to filter object."""
    errors = 0
    for or_filters in filters:
        for filter in or_filters:
            # Catch the error if it is an invalid table.column
            try:
                filter_obj.add(filter)
            except:
                print(f"Invalid filter: {filter}")
                errors += 1
    if errors > 0:
        print("Unable to run compare pipeline.")
        sys.exit(1)


# TODO refactor and test.
def get_dbs(pdb, gbk):
    """Create set of databases to compare to MySQL."""
    dbs = set()
    if pdb == True:
        dbs.add("phagesdb")
    if gbk == True:
        dbs.add("genbank")
    return dbs


# TODO refactor and test.
def output_to_file(data_list, folder, filename):
    """Output list data to file."""
    file_path = pathlib.Path(folder, filename)
    with file_path.open("w") as handle:
        writer = csv.writer(handle)
        for element in data_list:
            writer.writerow(element)


# TODO refactor and test.
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


# TODO refactor and test.
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


# TODO refactor and test.
def save_seqrecord(seqrecord, output_path, file_prefix, ext,
                   seqrecord_ext, interactive):
    """Save record to file."""
    file_path = basic.make_new_file(output_path, file_prefix, ext, attempt=100)
    if file_path is not None:
        SeqIO.write(seqrecord, file_path, seqrecord_ext)
    else:
        print(f"Duplicated filenames for {file_prefix}. Unable to output data.")
        if interactive:
            input('Press ENTER to proceed')


# TODO refactor and test.
def record_compare_settings(working_path, database, version, filters, valid_dbs):
    """Save user-selected settings."""

    # Format setting strings
    lst = [
        ["Comparison date:", CURRENT_DATE],
        ["MySQL database:", database],
        ["MySQL database version:", version],
        ["Databases compared to MySQL:", ", ".join(valid_dbs)]
        ]

    # Add the filters
    count = 1
    for or_filters in filters:
        for filter in or_filters:
            lst.append([f"Filter {count}:", filter])
            count += 1

    output_to_file(lst, pathlib.Path(working_path), COMPARE_SETTINGS)

# TODO refactor and test.
def selected_authors_lst(lst1):
    lst1 = [str(i) for i in list(lst1)]
    lst2 =  ["Genomes with the following AnnotationAuthor will be compared:",
             ", ".join(lst1)]
    return lst2

# TODO refactor and test.
def process_mysql_data(working_path, engine, phage_ids, interactive, save):
    """Retrieve and process MySQL data."""

    print('\n\nPreparing genome data sets from the MySQL database...')
    genome_list = mysqldb.parse_genome_data(engine,
                                            phage_id_list=phage_ids,
                                            phage_query=PHAGE_QUERY,
                                            gene_query=GENE_QUERY,
                                            gnm_type=GNM_MYSQL)

    tup = filter_mysql_genomes(genome_list)
    gnm_dict = tup[0]
    name_duplicates = tup[2]
    accessions = tup[3]
    accession_dupes = tup[4]

    set_mysql_gnm_attr(gnm_dict)

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

    return gnm_dict, accessions, accession_dupes


# TODO refactor and test.
def filter_mysql_genomes(genome_list):
    """Only keep selected subset of genomes with no value duplications."""

    gnm_dict = {}
    names = set()
    duplicate_names = set()
    accessions = set()
    duplicate_accessions = set()

    for gnm in genome_list:

        gnm_dict[gnm.id] = gnm

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


# TODO refactor and test.
def set_mysql_gnm_attr(gnm_dict):
    """Set compare-specific attributes not in pdm_utils Genome class."""
    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        for cds_ftr in gnm.cds_features:
            set_mysql_cds_attr(cds_ftr)

# TODO refactor and test.
def save_gnms_to_fasta(gnm_dict, main_path, new_dir, interactive):
    """Save genome data to fasta file."""

    output_path = pathlib.Path(main_path, new_dir)
    output_path.mkdir()

    for key in gnm_dict.keys():
        gnm = gnm_dict[key]
        fasta = SeqRecord(gnm.seq, id=gnm.id, description="")
        save_seqrecord(fasta, output_path, gnm.id, "fasta",
                       "fasta", interactive)


# TODO refactor and test.
def set_mysql_cds_attr(cds_ftr):
    """Set compare-specific MySQL Cds attributes not in pdm_utils Cds class."""
    cds_ftr.set_search_genome_id()
    cds_ftr.set_start_end_strand_id()
    cds_ftr.type = CDS_MYSQL

# TODO refactor and test.
def check_mysql_cds(cds_ftr):
    """Check for errors in MySQL CDS feature."""
    cds_ftr.compute_amino_acid_errors(constants.PROTEIN_ALPHABET)
    cds_ftr.compute_boundary_error()


# TODO refactor and test.
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


# TODO refactor and test.
def process_phagesdb_data(working_path, interactive, save):
    """Retrieve data from PhagesDB and process results."""

    gbd_tup = get_pdb_data(interactive)
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


# TODO refactor and test.
def get_pdb_data(interactive):
    """Retrieve data from PhagesDB."""
    print('\n\nRetrieving data from PhagesDB...')
    gnm_dict = {}
    names = set()
    dupe_names = set()
    data_list = phagesdb.get_phagesdb_data(constants.API_SEQUENCED)

    if len(data_list) > 0:
        for i in range(len(data_list)):
            # Note: this step take a long time because it needs to retrieve
            # and parse fasta files to get the genome sequence.
            gnm = phagesdb.parse_genome_data(data_list[i], gnm_type=GNM_PDB,
                                             seq=True)
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


# TODO refactor and test.
def check_pdb_gnms(gnm_dict):
    """Check for errors in PhagesDB genome."""
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]
        if str(gnm.seq) != "":
            gnm.compute_nucleotide_errors(constants.DNA_ALPHABET)

# TODO refactor and test.
def process_gbk_data(working_path, ncbi_creds, accessions, interactive, save):
    """Retrieve and process GenBank data."""

    if save == True:
        output_path = pathlib.Path(working_path, RECORD_FOLDER, GENBANK_OUTPUT)
        output_path.mkdir()

    records, retrieval_errors = get_genbank_data(ncbi_creds, accessions)

    #Report the accessions that could not be retrieved.
    output_to_file(retrieval_errors, pathlib.Path(working_path, ERROR_FOLDER),
                   FAILED_ACC_RETRIEVE)

    gnm_dict = {}
    for record in records:
        gnm = flat_files.parse_genome_data(record,gnm_type=GNM_GBK)
        set_gbk_gnm_attr(gnm)
        gnm_dict[gnm.accession] = gnm
        if save == True:
            save_gbk_genome(gnm, record, output_path, interactive)

    return gnm_dict


# TODO refactor and test.
def get_genbank_data(ncbi_cred_dict, accession_set, batch_size=200):
    """Retrieve genomes from GenBank."""

    print("\n\nRetrieving records from GenBank")

    # Use esearch to verify the accessions are valid and
    # efetch to retrieve the record.
    ncbi.set_entrez_credentials(
        tool=ncbi_cred_dict["tool"],
        email=ncbi_cred_dict["email"],
        api_key=ncbi_cred_dict["api_key"])


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


# TODO refactor and test.
def set_gbk_gnm_attr(gnm):
    """Set compare-specific attributes not in pdm_utils Genome class."""
    gnm.record_name = ""
    gnm.record_id = ""

    if len(gnm.source_features) == 1:
        src_ftr = gnm.source_features[0]
        gnm.source_feature_organism = src_ftr.organism
        gnm.source_feature_host = src_ftr.host
        gnm.source_feature_lab_host = src_ftr.lab_host

    for cds_ftr in gnm.cds_features:
        set_gbk_cds_attr(cds_ftr)


# TODO refactor and test.
def check_gbk_gnms(gnm_dict):
    """Check for errors in GenBank genome and CDS features."""
    for id in gnm_dict.keys():
        gnm = gnm_dict[id]

        for cds_ftr in gnm.cds_features:
            check_gbk_cds(cds_ftr)
        gnm.compute_nucleotide_errors(constants.DNA_ALPHABET)
        gnm.compute_cds_feature_errors()
        gnm.compute_gbk_cds_feature_errors()


# TODO refactor and test.
def set_gbk_cds_attr(cds_ftr):
    """Set compare-specific Cds attributes not in pdm_utils Cds class."""
    cds_ftr.set_start_end_strand_id()
    cds_ftr.type = CDS_GBK


# TODO refactor and test.
def check_gbk_cds(cds_ftr):
    """Check for errors in GenBank CDS feature."""
    cds_ftr.check_locus_tag()
    cds_ftr.compute_amino_acid_errors(constants.PROTEIN_ALPHABET)
    cds_ftr.compute_boundary_error()
    cds_ftr.compute_description_error()


# TODO refactor and test.
def save_gbk_genome(gnm, record, output_path, interactive):
    """Save GenBank record to file."""
    file_prefix = gnm.id + "__" + gnm.accession
    save_seqrecord(record, output_path, file_prefix,
                   "gb", "genbank", interactive)


# TODO refactor and test.
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
        gnm_triad = genometriad.GenomeTriad()
        gnm_triad.m_genome = mysql_gnm

        # Match up PhagesDB genome
        # First try to match up the phageID, and if that doesn't work,
        # try to match up the phageName
        if mysql_gnm.id in pdb_gnms.keys():
            pdb_genome = pdb_gnms[mysql_gnm.id]
            # Make sure the pdb_genome doesn't have a search name
            # that was duplicated
            if pdb_genome.id in pdb_name_duplicates:
                # Do NOT set genome.type, so that it remains ""
                pdb_genome = genome.Genome()

        else:
            # Do NOT set genome.type, so that it remains ""
            pdb_genome = genome.Genome()
            mysql_unmatched_to_pdb_gnms.append(mysql_gnm)

        gnm_triad.p_genome = pdb_genome

        # Now match up GenBank genome
        acc = mysql_gnm.accession
        if acc != "" and acc not in mysql_acc_duplicates:
                # Retrieval of record may have failed
                try:
                    g_gnm = gbk_gnms[acc]
                except:
                    # Do NOT set genome.type, so that it remains ""
                    g_gnm = genome.Genome()
                    mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        else:
            # Do NOT set genome.type, so that it remains ""
            g_gnm = genome.Genome()
            mysql_unmatched_to_gbk_gnms.append(mysql_gnm)

        gnm_triad.g_genome = g_gnm
        matched_gnms.append(gnm_triad)
        match_count += 1

    return (matched_gnms, mysql_unmatched_to_pdb_gnms, mysql_unmatched_to_gbk_gnms)

# TODO refactor and test.
def record_unmatched_gbk_data(gnms, working_path):
    """Save unmatched MySQL data not matched to GenBank."""
    lst1 = prepare_unmatched_to_gbk_output(gnms)
    output_to_file(lst1, pathlib.Path(working_path, ERROR_FOLDER),
                   FAILED_ACC_RETRIEVE)

# TODO refactor and test.
def record_unmatched_pdb_data(gnms, working_path):
    """Save unmatched MySQL data not matched to PhagesDB."""
    lst2 = prepare_unmatched_to_pdb_output(gnms)
    output_to_file(lst2, pathlib.Path(working_path, ERROR_FOLDER),
                   UNMATCHED_GENOMES)

# TODO refactor and test.
def check_matched_gnms(gnm_triads):
    """Compare all matched data."""

    print("Comparing matched genomes and identifying inconsistencies...")
    count = 1
    total = len(gnm_triads)
    for gnm_triad in gnm_triads:
        print(f"Comparing matched genome set {count} of {total}")

        gnm_triad.compare_mysql_gbk_genomes(GNM_MYSQL, GNM_GBK,
                                            CDSPAIR_MYSQL_GBK)
        gnm_triad.compare_mysql_phagesdb_genomes(GNM_MYSQL, GNM_PDB)
        gnm_triad.compare_phagesdb_gbk_genomes(GNM_PDB, GNM_GBK)
        gnm_triad.compute_total_genome_errors(GNM_GBK, GNM_PDB)
        count += 1


# TODO refactor and test.
def summarize_data(matched_genomes_list, working_path):
    """Create summary of data and save."""
    summary = dbcomparesummary.DbCompareSummary(matched_genomes_list)
    summary.compute_summary(GNM_MYSQL, GNM_PDB, GNM_GBK)
    output_all_data(working_path, summary)


# TODO refactor and test.
def output_all_data(output_path, summary):
    """Output all analysis results."""
    print("Outputting results to file...")

    cds_file_path = pathlib.Path(output_path, RESULTS_FOLDER, GENE_OUTPUT)
    cds_fh = cds_file_path.open("w")
    cds_writer = csv.writer(cds_fh)
    cds_headers = create_gene_headers()
    cds_writer.writerow(cds_headers)

    #Create genome summary output file
    lst1 = create_genome_summary_data(summary)
    output_to_file(lst1, pathlib.Path(output_path, RESULTS_FOLDER),
                   GNM_SUMMARY_OUTPUT)

    #Create CDS summary output file
    lst3 = create_cds_summary_data(summary)
    output_to_file(lst3, pathlib.Path(output_path, RESULTS_FOLDER),
                   CDS_SUMMARY_OUTPUT)

    # Now iterate through matched objects.
    # All MySQL genomes are stored in a GenomeTriad object,
    # even if there are no PhagesDB or GenBank matches.
    # All but a few MySQL phages should be matched to PhagesDB
    # Many MySQL phages should be matched to GenBank
    genomes_data = []
    for gnm_triad in summary._matched_genomes_list:

        gnm_data = create_genome_data(gnm_triad)
        genomes_data.append(gnm_data)

        #Once all matched genome data has been outputted,
        # iterate through all matched gene data
        m_gnm = gnm_triad.m_genome
        all_ftrs = get_all_features(gnm_triad)
        for mixed_ftr in all_ftrs:
            ftr_data = create_feature_data(m_gnm, mixed_ftr)
            cds_writer.writerow(ftr_data)

    lst2 = [create_genome_headers()]
    lst2.extend(genomes_data)
    output_to_file(lst2, pathlib.Path(output_path, RESULTS_FOLDER),
                   GENOME_OUTPUT)
    cds_fh.close()
    return


# TODO refactor and test.
def get_all_features(gnm_triad):
    """Create list of all features."""

    perfect_match = gnm_triad._m_g_perfect_matched_ftrs
    imperfectly_match = gnm_triad._m_g_imperfect_matched_ftrs
    mysql_unmatch = gnm_triad._m_ftrs_unmatched_in_g
    gbk_unmatch = gnm_triad._g_ftrs_unmatched_in_m

    lst = []
    lst.extend(perfect_match)
    lst.extend(imperfectly_match)
    lst.extend(mysql_unmatch)
    lst.extend(gbk_unmatch)

    return lst


# TODO refactor and test.
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


# TODO refactor and test.
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


# TODO refactor and test.
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


# TODO refactor and test.
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


# TODO refactor and test.
def create_genome_summary_data(summary):
    """Create summary of all genome results."""

    lst1 = [
        # Column header
        "tally",

        # First output database summary data
        summary._m_total_gnms_analyzed,
        summary._m_gnms_unmatched_to_p_tally,
        summary._m_gnms_unmatched_to_g_tally,
        summary._total_genomes_with_errors,

        # MySQL data
        # General genome data
        summary._m_g_update_flag_tally,

        # Genome data checks
        summary._m_gnms_with_nucleotide_errors_tally,
        summary._m_gnms_with_translation_errors_tally,
        summary._m_gnms_with_boundary_errors_tally,
        summary._m_gnms_with_status_accession_error_tally,
        summary._m_gnms_with_status_description_error_tally,

        # PhagesDB data
        # Genome data checks
        summary._p_gnms_with_nucleotide_errors_tally,

        # GenBank data
        # Genome data checks
        summary._g_gnms_with_description_field_errors_tally,
        summary._g_gnms_with_nucleotide_errors_tally,
        summary._g_gnms_with_translation_errors_tally,
        summary._g_gnms_with_boundary_errors_tally,
        summary._g_gnms_with_missing_locus_tags_tally,
        summary._g_gnms_with_locus_tag_typos_tally,


        # MySQL-PhagesDB checks
        summary._m_p_seq_mismatch_tally,
        summary._m_p_seq_length_mismatch_tally,
        summary._m_p_cluster_mismatch_tally,
        summary._m_p_accession_mismatch_tally,
        summary._m_p_host_mismatch_tally,

        # MySQL-GenBank checks
        summary._m_g_seq_mismatch_tally,
        summary._m_g_seq_length_mismatch_tally,
        summary._m_g_header_phage_mismatch_tally,
        summary._m_g_header_host_mismatch_tally,
        summary._m_g_gnms_with_author_errors_tally,
        summary._m_g_gnms_with_imperfectly_matched_ftrs_tally,
        summary._m_g_gnms_with_unmatched_m_ftrs_tally,
        summary._m_g_gnms_with_unmatched_g_ftrs_tally,
        summary._m_g_gnms_with_different_descriptions_tally,
        summary._m_g_gnms_with_different_translations_tally,

        # PhagesDB-GenBank checks
        summary._p_g_seq_mismatch_tally,
        summary._p_g_seq_length_mismatch_tally
        ]

    lst2 = create_genome_summary_fields()
    lst3 = []
    for i in range(len(lst2)):
        lst3.append([lst2[i],lst1[i]])
    return lst3


# TODO refactor and test.
def create_cds_summary_data(summary):
    """Create summary of all CDS results."""
    lst1 = [
        # Column header
        "tally",

        # MySQL feature
        # Gene data checks
        summary._m_translation_errors_tally,
        summary._m_boundary_errors_tally,

        # GenBank feature
        # Gene data checks
        summary._g_translation_errors_tally,
        summary._g_boundary_errors_tally,
        summary._g_missing_locus_tags_tally,
        summary._g_locus_tag_typos_tally,
        summary._g_description_field_errors_tally,

        # MySQL-GenBank checks
        summary._m_g_different_descriptions_tally,
        summary._m_g_different_start_sites_tally,
        summary._m_g_different_translation_tally,
        summary._m_g_unmatched_m_ftrs_tally,
        summary._m_g_unmatched_g_ftrs_tally
        ]

    lst2 = create_cds_summary_fields()
    lst3 = []
    for i in range(len(lst2)):
        lst3.append([lst2[i],lst1[i]])
    return lst3

# TODO refactor and test.
def create_feature_data(mysql_gnm, mixed_ftr):
    """Create feature data to output."""

    lst = []

    # Gene summaries
    # Add MySQL genome search name to each gene row regardless of
    # the type of CDS data (matched or unmatched).
    lst.append(mysql_gnm.id) # matched MySQL search name
    lst.append(mixed_ftr._total_errors) # total # of errors for this gene

    # Default type attribute is empty ""
    cds_pair = cdspair.CdsPair()
    mysql_ftr = cds.Cds()
    gbk_ftr = cds.Cds()

    # Determine what exactly the mixed feature object is and assign to
    # the appropriate variable.
    if isinstance(mixed_ftr, cdspair.CdsPair):
        cds_pair = mixed_ftr
        mysql_ftr = mixed_ftr.cds1
        gbk_ftr = mixed_ftr.cds2
    elif isinstance(mixed_ftr, cds.Cds):
        if mixed_ftr.type == CDS_MYSQL:
            mysql_ftr = mixed_ftr
        elif mixed_ftr.type == CDS_GBK:
            gbk_ftr = mixed_ftr
        else:
            pass
    else:
        pass


    # MySQL feature
    if mysql_ftr.type == CDS_MYSQL:

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
    if gbk_ftr.type == CDS_GBK:

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
    if cds_pair.type == CDSPAIR_MYSQL_GBK:

        # If this is a matched CDS feature, both MySQL and
        # GenBank features should have identical unmatched_error value.
        lst.append(cds_pair.cds1._unmatched_error)

        lst.append(cds_pair.different_description)
        lst.append(cds_pair.different_start_site)
        lst.append(cds_pair.different_translation)
    else:
        lst.append(mixed_ftr._unmatched_error)
        lst.extend(["","",""])
    return lst



# TODO refactor and test.
def create_genome_data(gnm_triad):
    """Create genome data to output."""

    mysql_gnm = gnm_triad.m_genome
    pdb_gnm = gnm_triad.p_genome
    gbk_gnm = gnm_triad.g_genome

    lst = []

    # Genome summary data
    lst.append(mysql_gnm.id)
    lst.append(gnm_triad._contains_errors)

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
    if pdb_gnm.type == GNM_PDB:

        # General genome data
        lst.append(pdb_gnm.name)
        lst.append(pdb_gnm.name)
        lst.append(pdb_gnm.cluster)
        lst.append(pdb_gnm.subcluster)
        lst.append(pdb_gnm.host_genus)
        lst.append(pdb_gnm.accession)
        lst.append(pdb_gnm.length)

        # Genome data checks
        lst.append(pdb_gnm._nucleotide_errors)
    else:
        lst.extend(["","","","","","","",""])



    # GenBank data
    if gbk_gnm.type == GNM_GBK:

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
    if pdb_gnm.type == GNM_PDB:

        lst.append(gnm_triad._m_p_seq_mismatch)
        lst.append(gnm_triad._m_p_seq_length_mismatch)
        lst.append(gnm_triad._m_p_cluster_mismatch)
        lst.append(gnm_triad._m_p_accession_mismatch)
        lst.append(gnm_triad._m_p_host_mismatch)
    else:
        lst.extend(["","","","",""])


    # MySQL-GenBank checks
    if gbk_gnm.type == GNM_GBK:
        lst.append(gnm_triad._m_g_seq_mismatch)
        lst.append(gnm_triad._m_g_seq_length_mismatch)
        lst.append(gnm_triad._g_header_fields_name_mismatch)
        lst.append(gnm_triad._g_host_mismatch)
        lst.append(gnm_triad._m_g_author_error)
        lst.append(gnm_triad._m_g_perfect_matched_ftrs_tally)
        lst.append(gnm_triad._m_g_imperfect_matched_ftrs_tally)
        lst.append(gnm_triad._m_ftrs_unmatched_in_g_tally)
        lst.append(gnm_triad._g_ftrs_unmatched_in_m_tally)
        lst.append(gnm_triad._m_g_different_descriptions_tally)
        lst.append(gnm_triad._m_g_different_translations_tally)
    else:
        lst.extend(["","","","","","","","","","",""])

    # Number of genes with errors
    lst.append(gnm_triad._total_number_genes_with_errors)

    # Output PhagesDB-GenBank checks
    if pdb_gnm.type == GNM_PDB and gbk_gnm.type == GNM_GBK:
        lst.append(gnm_triad._p_g_seq_mismatch)
        lst.append(gnm_triad._p_g_seq_length_mismatch)
    else:
        lst.extend(["",""])
    return lst
