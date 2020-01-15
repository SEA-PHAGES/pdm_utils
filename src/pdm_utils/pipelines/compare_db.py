#!/usr/bin/env python
#Database comparison script
#University of Pittsburgh
#Travis Mavrich
#20170203
#The purpose of this script is to compare the MySQL, phagesdb, and NCBI databases for inconsistencies and report what needs to be updated.

# Note this script compares and matches data from Genbank data and MySQL data. As a result, there are many similarly
# named variables. Variables are prefixed to indicate database:
#NCBI =  "ncbi".
#MySQL = "ph"
#phagesdb = "pdb"


#Built-in libraries
import time, sys, os, getpass, csv, re
import json
import urllib.request
import argparse
import pathlib
from pdm_utils.functions import ncbi
from pdm_utils.functions import basic



#Third-party libraries
try:
    import pymysql as pms
    # import MySQLdb as mdb
    from Bio import SeqIO, Entrez
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
except:
    print("\nUnable to import one or more of the following third-party modules: \n\n\
        pymysql\n\
        Biopython SeqIO\n\
        Biopython Seq,\n\
        Biopython SeqRecord,\n\
        Biopython Entrez\n\
        Biopython IUPAC\n")
    print("\n\nInstall modules and try again.\n\n")
    sys.exit(1)



#Define several functions


#Exits MySQL
def mdb_exit(message):
    print("\nError")
    # print("\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
    print("\nThe import script did not complete.")
    print("\nExiting MySQL.")
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    print("\nExiting import script.")
    sys.exit(1)

#Closes all file handles currently open
def close_all_files(file_list):
    for file_handle in file_list:
        file_handle.close()

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


#Allows user to select specific options
def select_option(message,valid_response_set):

    response_valid = False
    while response_valid == False:
        response = input(message)
        if response.isdigit():
            response = int(response)
        else:
            response = response.lower()

        if response in valid_response_set:
            response_valid = True
            if response == 'y':
                response  = 'yes'
            elif response == 'n':
                response  = 'no'
        else:
            print('Invalid response.')
    return response



#Output list to file
def output_to_file(data_list,filename,genome_status_selected,database_string,genome_author_selected):
    filename_fh = open(os.path.join(main_output_path,date + "_" + filename), 'w')
    filename_writer = csv.writer(filename_fh)
    filename_writer.writerow([date + ' Database comparison'])
    filename_writer.writerow([database_string])
    filename_writer.writerow([genome_author_selected])
    filename_writer.writerow([genome_status_selected])
    for element in data_list:
        filename_writer.writerow([element])
    filename_fh.close()





#Ensure the output filename is unique
def create_unique_filename(filename_directory,filename_base,filename_ext):

    file_exists = True
    rename_counter = 0
    while file_exists == True:

        if rename_counter == 0:
            unique_filename = filename_base + filename_ext
        else:
            unique_filename = filename_base + '_' + str(rename_counter) + filename_ext

        unique_filename_path = os.path.join(filename_directory,unique_filename)
        file_exists = os.path.isfile(unique_filename_path)

        if file_exists == True:
            print('Warning: duplicate output file:')
            print(unique_filename_path)
            print('Filename will be modified.')
            input('Press ENTER to proceed')
            rename_counter += 1

    return unique_filename_path

















#Define data classes



#Base genome class
class UnannotatedGenome:

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields
        self.__phage_name = ''
        self.__host = ''
        self.__sequence = ''
        self.__accession = ''


        # Computed datafields
        self.__search_name = '' # The phage name void of "_Draft" and converted to lowercase
        self.__length = 0
        self.__nucleotide_errors = False


    # Define all attribute setters:
    def set_phage_name(self,value):
        self.__phage_name = value
        self.__search_name = remove_draft_suffix(self.__phage_name)
    def set_host(self,value):
        self.__host = value
    def set_sequence(self,value):
        self.__sequence = value.upper()
        self.__length = len(self.__sequence)
    def set_accession(self,value):
        if value is None or value.strip() == '':
            self.__accession = ''
        else:
            value = value.strip()
            self.__accession = value.split('.')[0]
    def compute_nucleotide_errors(self,dna_alphabet_set):
        nucleotide_set = set(self.__sequence)
        nucleotide_error_set = nucleotide_set - dna_alphabet_set
        if len(nucleotide_error_set) > 0:
            self.__nucleotide_errors = True


    # Define all attribute getters:
    def get_phage_name(self):
        return self.__phage_name
    def get_host(self):
        return self.__host
    def get_sequence(self):
        return self.__sequence
    def get_length(self):
        return self.__length
    def get_accession(self):
        return self.__accession
    def get_search_name(self):
        return self.__search_name
    def get_nucleotide_errors(self):
        return self.__nucleotide_errors




class AnnotatedGenome(UnannotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        UnannotatedGenome.__init__(self)

        # Non-computed datafields

        # Computed datafields
        self.__cds_features = []
        self.__cds_features_tally = 0
        self.__cds_features_with_translation_error_tally = 0
        self.__cds_features_boundary_error_tally = 0

    # Define all attribute setters:
    def set_cds_features(self,value):
        self.__cds_features = value #Should be a list
        self.__cds_features_tally = len(self.__cds_features)
    def compute_cds_feature_errors(self):
        for cds_feature in self.__cds_features:
            if cds_feature.get_amino_acid_errors():
                self.__cds_features_with_translation_error_tally += 1
            if cds_feature.get_boundary_error():
                self.__cds_features_boundary_error_tally += 1

    # Define all attribute getters:
    def get_cds_features(self):
        return self.__cds_features
    def get_cds_features_tally(self):
        return self.__cds_features_tally
    def get_cds_features_with_translation_error_tally(self):
        return self.__cds_features_with_translation_error_tally
    def get_cds_features_boundary_error_tally(self):
        return self.__cds_features_boundary_error_tally





class PhameratorGenome(AnnotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        AnnotatedGenome.__init__(self)
        # Non-computed datafields
        self.__phage_id = ''
        self.__status = '' #Final, Draft, Unknown version of genome data
        self.__cluster_subcluster = '' #Combined cluster_subcluster data.
        self.__ncbi_update_flag = ''
        self.__date_last_modified = ''
        self.__annotation_author = '' #1 (Hatfull), 0 (Genbank)

        # Computed datafields
        self.__search_id = '' # The phage ID void of "_Draft" and converted to lowercase
        self.__status_accession_error = False
        self.__status_description_error = False
        self.__description_tally = 0
        self.__genes_with_errors_tally = 0 #How many genes in this genome have at least one error?


    # Define all attribute setters:
    def set_phage_id(self,value):
        self.__phage_id = value
        self.__search_id = remove_draft_suffix(self.__phage_id)
    def set_status(self,value):
        self.__status = value

        #Be sure to first set the accession attribute before the status attribute,
        #else this will throw an error.
        #Now that the AnnotationAuthor field contains authorship data, the
        #'unknown' annotation status now reflects an 'unknown' annotation (in
        #regards to if it was auto-annotated or manually annotated).
        #So for the status-accession error, if the status is 'unknown',
        #there is no reason to assume whether there should be an accession or not.
        #Only for 'final' (manually annotated) genomes should there be an accession.
        #Old code, using 'unknown':
        # if (self.__status == 'final' or self.__status == 'unknown') and self.get_accession() == '':
        #     self.__status_accession_error = True
        if self.__status == 'final' and self.get_accession() == '':
            self.__status_accession_error = True



    def set_cluster_subcluster(self,value):
        if value is None:
            self.__cluster_subcluster = 'Singleton'
        elif value == 'UNK':
            self.__cluster_subcluster = ''
        else:
            self.__cluster_subcluster = value
    def set_ncbi_update_flag(self,value):
        self.__ncbi_update_flag = value
    def set_date_last_modified(self,value):
        if value is None:
            self.__date_last_modified = ''
        else:
            self.__date_last_modified = value
    def set_annotation_author(self,value):
        self.__annotation_author = value

    def compute_status_description_error(self):
        #Iterate through all CDS features, see if they have descriptions, then compare to the status
        for feature in self.get_cds_features():
            if feature.get_notes() != '':
                self.__description_tally += 1
        if self.__status == 'draft' and self.__description_tally > 0:
            self.__status_description_error = True
        elif self.__status == 'final' and self.__description_tally == 0:
            self.__status_description_error = True
        else:
            pass


    #Even though this method iterates through the cds features like the compute_status_description_error does,
    #it has to be kept separate, since you need to wait to run this method after all genome and gene matching
    #is completed.
    def compute_genes_with_errors_tally(self):
        for feature in self.get_cds_features():
            #Need to first compute the number of errors per gene
            feature.compute_total_cds_errors()
            if feature.get_total_errors() > 0:
                self.__genes_with_errors_tally += 1





    # Define all attribute getters:
    def get_phage_id(self):
        return self.__phage_id
    def get_cluster_subcluster(self):
        return self.__cluster_subcluster
    def get_status(self):
        return self.__status
    def get_search_id(self):
        return self.__search_id
    def get_ncbi_update_flag(self):
        return self.__ncbi_update_flag
    def get_date_last_modified(self):
        return self.__date_last_modified
    def get_status_description_error(self):
        return self.__status_description_error
    def get_status_accession_error(self):
        return self.__status_accession_error
    def get_description_tally(self):
        return self.__description_tally
    def get_genes_with_errors_tally(self):
        return self.__genes_with_errors_tally
    def get_annotation_author(self):
        return self.__annotation_author




class PhagesdbGenome(UnannotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        UnannotatedGenome.__init__(self)

        # Non-computed datafields
        self.__cluster = ''
        self.__subcluster = ''

        # Computed datafields

    # Define all attribute setters:
    def set_cluster(self,value):
        self.__cluster = value
    def set_subcluster(self,value):
        self.__subcluster = value

    # Define all attribute getters:
    def get_cluster(self):
        return self.__cluster
    def get_subcluster(self):
        return self.__subcluster



class NcbiGenome(AnnotatedGenome):

    # Initialize all attributes:
    def __init__(self):
        AnnotatedGenome.__init__(self)

        #Non-computed data fields
        self.__record_name = ''
        self.__record_id = ''
        self.__record_accession = ''
        self.__record_description = ''
        self.__record_source = ''
        self.__record_organism = ''
        self.__source_feature_organism = ''
        self.__source_feature_host = ''
        self.__source_feature_lab_host = ''
        self.__record_authors = ''


        #Computed data fields
        self.__function_descriptions_tally = 0
        self.__product_descriptions_tally = 0
        self.__note_descriptions_tally = 0
        self.__missing_locus_tags_tally = 0
        self.__locus_tag_typos_tally = 0
        self.__description_field_error_tally = 0



    #Define setter functions
    def set_record_name(self,value):
        self.__record_name = value
    def set_record_id(self,value):
        self.__record_id = value
    def set_record_accession(self,value):
        self.__record_accession = value
    def set_record_description(self,value):
        self.__record_description = value
    def set_record_source(self,value):
        self.__record_source = value
    def set_record_organism(self,value):
        self.__record_organism = value
    def set_source_feature_organism(self,value):
        self.__source_feature_organism = value
    def set_source_feature_host(self,value):
        self.__source_feature_host = value
    def set_source_feature_lab_host(self,value):
        self.__source_feature_lab_host = value
    def compute_ncbi_cds_feature_errors(self):
        for cds_feature in self.get_cds_features():

            #counting descriptions should skip if it is blank or "hypothetical protein"
            if cds_feature.get_search_product_description() != '':
                self.__product_descriptions_tally += 1

            if cds_feature.get_search_function_description() != '':
                self.__function_descriptions_tally += 1

            if cds_feature.get_search_note_description() != '':
                self.__note_descriptions_tally += 1



            if cds_feature.get_locus_tag_missing():
                self.__missing_locus_tags_tally += 1
            else:
                pattern4 = re.compile(self.get_search_name())
                search_result = pattern4.search(cds_feature.get_locus_tag().lower())

                if search_result == None:
                    self.__locus_tag_typos_tally += 1
                    cds_feature.set_locus_tag_typo() #Sets this attribute to True

            if cds_feature.get_description_field_error():
                self.__description_field_error_tally += 1

    def set_record_authors(self,value):
        self.__record_authors = value



    #Define getter functions
    def get_record_name(self):
        return self.__record_name
    def get_record_id(self):
        return self.__record_id
    def get_record_accession(self):
        return self.__record_accession
    def get_record_description(self):
        return self.__record_description
    def get_record_source(self):
        return self.__record_source
    def get_record_organism(self):
        return self.__record_organism
    def get_source_feature_organism(self):
        return self.__source_feature_organism
    def get_source_feature_host(self):
        return self.__source_feature_host
    def get_source_feature_lab_host(self):
        return self.__source_feature_lab_host
    def get_function_descriptions_tally(self):
        return self.__function_descriptions_tally
    def get_product_descriptions_tally(self):
        return self.__product_descriptions_tally
    def get_note_descriptions_tally(self):
        return self.__note_descriptions_tally
    def get_missing_locus_tags_tally(self):
        return self.__missing_locus_tags_tally
    def get_locus_tag_typos_tally(self):
        return self.__locus_tag_typos_tally
    def get_description_field_error_tally(self):
        return self.__description_field_error_tally
    def get_record_authors(self):
        return self.__record_authors



class CdsFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from MySQL database:
        self.__type_id = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.__left_boundary = '' #Position of left boundary of feature, 0-indexed
        self.__right_boundary = '' #Position of right boundary of feature, 0-indexed
        self.__strand = '' #'forward', 'reverse', or 'NA'
        self.__translation = ''
        self.__translation_length = ''

        # Computed datafields
        self.__amino_acid_errors = False
        self.__start_end_strand_id = ''
        self.__end_strand_id = ''
        self.__boundary_error = False
        self.__unmatched_error = False #keeps track if it is contains a match or not

    # Define all attribute setters:
    def set_left_boundary(self,value):
        self.__left_boundary = value
    def set_right_boundary(self,value):
        self.__right_boundary = value
    def set_strand(self,value):
        self.__strand = parse_strand(value)
    def set_translation(self,value):
        self.__translation = value.upper()
        self.__translation_length = len(self.__translation)
    def set_type_id(self,value):
        self.__type_id = value
    def compute_amino_acid_errors(self,protein_alphabet_set):
        amino_acid_set = set(self.__translation)
        amino_acid_error_set = amino_acid_set - protein_alphabet_set
        if len(amino_acid_error_set) > 0:
            self.__amino_acid_errors = True
    def set_start_end_strand_id(self):
        #Create a tuple of feature location data.
        #For start and end of feature, it doesn't matter whether the feature is complex with a translational
        #frameshift or not. Retrieving the "start" and "end" attributes return the very beginning and end of
        #the feature, disregarding the inner "join" coordinates.
        self.__start_end_strand_id = (str(self.__left_boundary),str(self.__right_boundary),self.__strand)

        #Since this id matched genes with different start sites,
        #the strand impacts whether the left or right boundary is used
        if self.__strand == 'forward':
            self.__end_strand_id = (str(self.__right_boundary),self.__strand)
        elif self.__strand == 'reverse':
            self.__end_strand_id = (str(self.__left_boundary),self.__strand)
        else:
            pass
    def set_unmatched_error(self):
        self.__unmatched_error = True
    def compute_boundary_error(self):
        #Check if start and end coordinates are fuzzy
        if not (str(self.__left_boundary).isdigit() and str(self.__right_boundary).isdigit()):
            self.__boundary_error = True



    # Define all attribute getters:
    def get_left_boundary(self):
        return self.__left_boundary
    def get_right_boundary(self):
        return self.__right_boundary
    def get_type_id(self):
        return self.__type_id
    def get_strand(self):
        return self.__strand
    def get_amino_acid_errors(self):
        return self.__amino_acid_errors
    def get_translation(self):
        return self.__translation
    def get_translation_length(self):
        return self.__translation_length
    def get_start_end_strand_id(self):
        return self.__start_end_strand_id
    def get_end_strand_id(self):
        return self.__end_strand_id
    def get_unmatched_error(self):
        return self.__unmatched_error
    def get_boundary_error(self):
        return self.__boundary_error


class PhameratorCdsFeature(CdsFeature):

    # Initialize all attributes:
    def __init__(self):
        CdsFeature.__init__(self)

        # Initialize all non-calculated attributes:

        #Datafields from MySQL database:
        self.__phage_id = ''
        self.__gene_id = '' #Gene ID comprised of PhageID and Gene name
        self.__gene_name = ''
        self.__notes = ''
        self.__search_notes = '' #non-generic gene descriptions

        # Computed datafields
        self.__search_id = ''
        self.__total_errors = 0

    # Define all attribute setters:
    def set_phage_id(self,value):
        self.__phage_id = value
        self.__search_id = remove_draft_suffix(self.__phage_id)
    def set_gene_id(self,value):
        self.__gene_id = value
    def set_gene_name(self,name):
        self.__gene_name = name
    def set_notes(self,value1,value2):
        self.__notes = value1
        self.__search_notes = value2
    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self.__total_errors += 1
        if self.get_boundary_error():
            self.__total_errors += 1
        if self.get_unmatched_error():
            self.__total_errors += 1



    # Define all attribute getters:
    def get_gene_id(self):
        return self.__gene_id
    def get_gene_name(self):
        return self.__gene_name
    def get_notes(self):
        return self.__notes
    def get_search_notes(self):
        return self.__search_notes
    def get_phage_id(self):
        return self.__phage_id
    def get_search_id(self):
        return self.__search_id
    def get_total_errors(self):
        return self.__total_errors



class NcbiCdsFeature(CdsFeature):

    # Initialize all attributes:
    def __init__(self):
        CdsFeature.__init__(self)

        # Initialize all non-calculated attributes:
        self.__locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.__gene_number = ''
        self.__product_description = ''
        self.__function_description = ''
        self.__note_description = ''
        self.__search_product_description = ''
        self.__search_function_description = ''
        self.__search_note_description = ''

        # Inititalize all calculated attributes:
        self.__locus_tag_missing = False
        self.__locus_tag_typo = False #This can only be computed at the genome level since it uses the phage name
        self.__description_field_error = False
        self.__total_errors = 0


    # Define all attribute setters:
    def set_locus_tag(self,value):
        self.__locus_tag = value
        if self.__locus_tag == '':
            self.__locus_tag_missing = True
    def set_gene_number(self,value):
        self.__gene_number = value
    def set_product_description(self,value1,value2):
        self.__product_description = value1
        self.__search_product_description = value2
    def set_function_description(self,value1,value2):
        self.__function_description = value1
        self.__search_function_description = value2
    def set_note_description(self,value1,value2):
        self.__note_description = value1
        self.__search_note_description = value2
    def set_locus_tag_typo(self):
        self.__locus_tag_typo = True

    def compute_description_error(self):

        #If the product description is empty or generic, and the function or note descriptions are not, there is an error
        if self.__search_product_description == '' and \
            (self.__search_function_description != '' or \
            self.__search_note_description != ''):

            self.__description_field_error = True

    def compute_total_cds_errors(self):
        if self.get_amino_acid_errors():
            self.__total_errors += 1
        if self.get_boundary_error():
            self.__total_errors += 1
        if self.__description_field_error:
            self.__total_errors += 1
        if self.__locus_tag_missing:
            self.__total_errors += 1
        if self.__locus_tag_typo:
            self.__total_errors += 1
        if self.get_unmatched_error():
            self.__total_errors += 1



    # Define all attribute getters:
    def get_locus_tag(self):
        return self.__locus_tag
    def get_gene_number(self):
        return self.__gene_number
    def get_product_description(self):
        return self.__product_description
    def get_function_description(self):
        return self.__function_description
    def get_note_description(self):
        return self.__note_description
    def get_search_product_description(self):
        return self.__search_product_description
    def get_search_function_description(self):
        return self.__search_function_description
    def get_search_note_description(self):
        return self.__search_note_description
    def get_locus_tag_missing(self):
        return self.__locus_tag_missing
    def get_locus_tag_typo(self):
        return self.__locus_tag_typo
    def get_description_field_error(self):
        return self.__description_field_error
    def get_total_errors(self):
        return self.__total_errors




class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.__phamerator_genome = ''
        self.__phagesdb_genome = ''
        self.__ncbi_genome = ''

        #MySQL and NCBI matched data comparison results
        self.__phamerator_ncbi_sequence_mismatch = False
        self.__phamerator_ncbi_sequence_length_mismatch = False
        self.__ncbi_record_header_fields_phage_name_mismatch = False
        self.__ncbi_host_mismatch = False
        self.__phamerator_ncbi_perfect_matched_features = [] #List of MatchedCdsFeature objects
        self.__phamerator_ncbi_imperfect_matched_features = [] #List of MatchedCdsFeature objects
        self.__phamerator_features_unmatched_in_ncbi = [] #List of CdsFeature objects
        self.__ncbi_features_unmatched_in_phamerator = [] #List of CdsFeature objects
        self.__phamerator_ncbi_perfect_matched_features_tally = 0
        self.__phamerator_ncbi_imperfect_matched_features_tally = 0
        self.__phamerator_features_unmatched_in_ncbi_tally = 0
        self.__ncbi_features_unmatched_in_phamerator_tally = 0
        self.__phamerator_ncbi_different_descriptions_tally = 0
        self.__phamerator_ncbi_different_translations_tally = 0
        self.__ph_ncbi_author_error = False





        #MySQL and phagesdb matched data comparison results
        self.__phamerator_phagesdb_sequence_mismatch = False
        self.__phamerator_phagesdb_sequence_length_mismatch = False
        self.__phamerator_phagesdb_host_mismatch = False
        self.__phamerator_phagesdb_accession_mismatch = False
        self.__phamerator_phagesdb_cluster_subcluster_mismatch = False


        #phagesdb and NCBI matched data comparison results
        self.__phagesdb_ncbi_sequence_mismatch = False
        self.__phagesdb_ncbi_sequence_length_mismatch = False

        #Total errors summary
        self.__contains_errors = False
        self.__total_number_genes_with_errors = 0



    # Define all attribute setters:
    def set_phamerator_genome(self,value):
        self.__phamerator_genome = value
    def set_phagesdb_genome(self,value):
        self.__phagesdb_genome = value
    def set_ncbi_genome(self,value):
        self.__ncbi_genome = value

    def compare_phamerator_ncbi_genomes(self):

        #verify that there is a MySQL and NCBI genome in the matched genome object
        ph_genome = self.__phamerator_genome
        ncbi_genome = self.__ncbi_genome
        ph_cds_list = ph_genome.get_cds_features()

        if isinstance(ph_genome,PhameratorGenome) and isinstance(ncbi_genome,NcbiGenome):

            if ph_genome.get_sequence() != ncbi_genome.get_sequence():
                self.__phamerator_ncbi_sequence_mismatch = True
            if ph_genome.get_length() != ncbi_genome.get_length():
                self.__phamerator_ncbi_sequence_length_mismatch = True

            #Compare phage names
            pattern1 = re.compile('^' + ph_genome.get_phage_name() + '$')
            pattern2 = re.compile('^' + ph_genome.get_phage_name())

            if find_name(pattern2,ncbi_genome.get_record_description().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_record_source().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_record_organism().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_source_feature_organism().split(' ')) == 0:

                self.__ncbi_record_header_fields_phage_name_mismatch = True

            #Compare host data
            search_host = ph_genome.get_host()
            if search_host == 'Mycobacterium':
                search_host = search_host[:-3]
            pattern3 = re.compile('^' + search_host)

            if (find_name(pattern3,ncbi_genome.get_record_description().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_record_source().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_record_organism().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_source_feature_organism().split(' ')) == 0) or \
                (ncbi_genome.get_source_feature_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_host().split(' ')) == 0) or \
                (ncbi_genome.get_source_feature_lab_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_lab_host().split(' ')) == 0):

                self.__ncbi_host_mismatch = True

            #Check author list for errors
            #For genomes with AnnotationAuthor = 1 (Hatfull), Graham is expected
            #to be an author.
            #For genomes with AnnotationAuthor = 0 (non-Hatfull/Genbank), Graham
            #is NOT expected to be an author.
            pattern5 = re.compile('hatfull')
            search_result = pattern5.search(ncbi_genome.get_record_authors().lower())
            if ph_genome.get_annotation_author() == 1 and search_result == None:
                self.__ph_ncbi_author_error = True
            elif ph_genome.get_annotation_author() == 0 and search_result != None:
                self.__ph_ncbi_author_error = True
            else:
                #Any other combination of MySQL and ncbi author can be skipped
                pass


            #Compare CDS features

            #First find all unique start-end-strand cds identifiers for MySQL and ncbi genomes
            ph_start_end_strand_id_set = set()
            ph_start_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() not in ph_start_end_strand_id_set:
                    ph_start_end_strand_id_set.add(cds.get_start_end_strand_id())
                else:
                    ph_start_end_strand_duplicate_id_set.add(cds.get_start_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ph_start_end_strand_id_set = ph_start_end_strand_id_set - ph_start_end_strand_duplicate_id_set


            ncbi_cds_list = ncbi_genome.get_cds_features()
            ncbi_start_end_strand_id_set = set()
            ncbi_start_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() not in ncbi_start_end_strand_id_set:
                    ncbi_start_end_strand_id_set.add(cds.get_start_end_strand_id())

                else:
                    ncbi_start_end_strand_duplicate_id_set.add(cds.get_start_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ncbi_start_end_strand_id_set = ncbi_start_end_strand_id_set - ncbi_start_end_strand_duplicate_id_set

            #Create the perfect matched and unmatched sets
            ph_unmatched_start_end_strand_id_set = ph_start_end_strand_id_set - ncbi_start_end_strand_id_set
            ncbi_unmatched_start_end_strand_id_set = ncbi_start_end_strand_id_set - ph_start_end_strand_id_set
            perfect_matched_cds_id_set = ph_start_end_strand_id_set & ncbi_start_end_strand_id_set


            #From the unmatched sets, created second round of end-strand id sets
            #MySQL end_strand data
            ph_end_strand_id_set = set()
            ph_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() in ph_unmatched_start_end_strand_id_set:
                    if cds.get_end_strand_id() not in ph_end_strand_id_set:
                        ph_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ph_end_strand_duplicate_id_set.add(cds.get_end_strand_id())

            #Remove the duplicate end_strand ids from the main id_set
            ph_end_strand_id_set = ph_end_strand_id_set - ph_end_strand_duplicate_id_set


            ncbi_end_strand_id_set = set()
            ncbi_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() in ncbi_unmatched_start_end_strand_id_set:
                    if cds.get_end_strand_id() not in ncbi_end_strand_id_set:
                        ncbi_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ncbi_end_strand_duplicate_id_set.add(cds.get_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ncbi_end_strand_id_set = ncbi_end_strand_id_set - ncbi_end_strand_duplicate_id_set


            #Create the imperfect matched set
            imperfect_matched_cds_id_set = ph_end_strand_id_set & ncbi_end_strand_id_set


            #Now go back through all cds features and assign to the appropriate dictionary or list
            ph_perfect_matched_cds_dict = {}
            ph_imperfect_matched_cds_dict = {}
            ph_unmatched_cds_list = []


            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ph_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ph_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ph_unmatched_cds_list.append(cds)


            ncbi_perfect_matched_cds_dict = {}
            ncbi_imperfect_matched_cds_dict = {}
            ncbi_unmatched_cds_list = []

            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ncbi_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ncbi_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ncbi_unmatched_cds_list.append(cds)


            #Create MatchedCdsFeatures objects
            #Compute matched gene errors and tallies
            #Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1



                if matched_cds_object.get_phamerator_ncbi_different_translations():
                    self.__phamerator_ncbi_different_translations_tally += 1
                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_perfect_matched_features.append(matched_cds_object)


            #Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_imperfect_matched_features.append(matched_cds_object)


            #Compute unmatched error and gene total errors for all unmatched features
            for cds in ph_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

            for cds in ncbi_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1




            #Set unmatched cds lists
            self.__phamerator_features_unmatched_in_ncbi = ph_unmatched_cds_list
            self.__ncbi_features_unmatched_in_phamerator = ncbi_unmatched_cds_list

            #Now compute the number of features in each category
            self.__phamerator_ncbi_perfect_matched_features_tally = len(self.__phamerator_ncbi_perfect_matched_features)
            self.__phamerator_ncbi_imperfect_matched_features_tally = len(self.__phamerator_ncbi_imperfect_matched_features)
            self.__phamerator_features_unmatched_in_ncbi_tally = len(self.__phamerator_features_unmatched_in_ncbi)
            self.__ncbi_features_unmatched_in_phamerator_tally = len(self.__ncbi_features_unmatched_in_phamerator)


        #If there is no matching NCBI genome, assign all MySQL genes to Unmatched
        else:

            #Set unmatched cds lists, but do NOT count them in the unmatched tally.
            #The unmatched tally should reflect unmatched genes if there is actually a metching NCBI genome.
            self.__phamerator_features_unmatched_in_ncbi = ph_cds_list

            #Now that all errors have been computed for each gene, compute how many genes have errors
            #If there is no matching NCBI genome, gene error tallies are only computed for the MySQL genome
            ph_genome.compute_genes_with_errors_tally()
            self.__total_number_genes_with_errors = ph_genome.get_genes_with_errors_tally()




    def compare_phamerator_phagesdb_genomes(self):

        #verify that there is a MySQL and phagesdb genome in the matched genome object
        ph_genome = self.__phamerator_genome
        pdb_genome = self.__phagesdb_genome

        if isinstance(ph_genome,PhameratorGenome) and isinstance(pdb_genome,PhagesdbGenome):

            if ph_genome.get_sequence() != pdb_genome.get_sequence():
                self.__phamerator_phagesdb_sequence_mismatch = True
            if ph_genome.get_length() != pdb_genome.get_length():
                self.__phamerator_phagesdb_sequence_length_mismatch = True
            if ph_genome.get_accession() != pdb_genome.get_accession():
                self.__phamerator_phagesdb_accession_mismatch = True
            if ph_genome.get_host() != pdb_genome.get_host():
                self.__phamerator_phagesdb_host_mismatch = True
            if ph_genome.get_cluster_subcluster() != pdb_genome.get_cluster() and \
                ph_genome.get_cluster_subcluster() != pdb_genome.get_subcluster():

                self.__phamerator_phagesdb_cluster_subcluster_mismatch = True


    def compare_phagesdb_ncbi_genomes(self):

        #verify that there is a phagesdb and NCBI genome in the matched genome object
        pdb_genome = self.__phagesdb_genome
        ncbi_genome = self.__ncbi_genome

        if isinstance(pdb_genome,PhagesdbGenome) and isinstance(ncbi_genome,NcbiGenome):
            if pdb_genome.get_sequence() != ncbi_genome.get_sequence():
                self.__phagesdb_ncbi_sequence_mismatch = True
            if pdb_genome.get_length() != ncbi_genome.get_length():
                self.__phagesdb_ncbi_sequence_length_mismatch = True



    def compute_total_genome_errors(self):

        ph_genome = self.__phamerator_genome
        ncbi_genome = self.__ncbi_genome
        pdb_genome = self.__phagesdb_genome


        if self.__total_number_genes_with_errors > 0:
            self.__contains_errors = True


        #MySQL genome
        if ph_genome.get_nucleotide_errors():
            self.__contains_errors = True
        if ph_genome.get_status_description_error():
            self.__contains_errors = True
        if ph_genome.get_status_accession_error():
            self.__contains_errors = True


        #NCBI genome
        if isinstance(ncbi_genome,NcbiGenome):
            if ncbi_genome.get_nucleotide_errors():
                self.__contains_errors = True



        #phagesdb genome
        if isinstance(pdb_genome,PhagesdbGenome):
            if pdb_genome.get_nucleotide_errors():
                self.__contains_errors = True


        #MySQL-NCBI
        if self.__phamerator_ncbi_sequence_mismatch:
            self.__contains_errors = True
        if self.__phamerator_ncbi_sequence_length_mismatch:
            self.__contains_errors = True
        if self.__ncbi_record_header_fields_phage_name_mismatch:
            self.__contains_errors = True
        if self.__ncbi_host_mismatch:
            self.__contains_errors = True
        if self.__phamerator_ncbi_imperfect_matched_features_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_features_unmatched_in_ncbi_tally > 0:
            self.__contains_errors = True
        if self.__ncbi_features_unmatched_in_phamerator_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_ncbi_different_descriptions_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_ncbi_different_translations_tally > 0:
            self.__contains_errors = True
        if self.__ph_ncbi_author_error:
            self.__contains_errors = True


        #MySQL-phagesdb
        if self.__phamerator_phagesdb_sequence_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_sequence_length_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_host_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_accession_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_cluster_subcluster_mismatch:
            self.__contains_errors = True

        #phagesdb-NCBI
        if self.__phagesdb_ncbi_sequence_mismatch:
            self.__contains_errors = True
        if self.__phagesdb_ncbi_sequence_length_mismatch:
            self.__contains_errors = True


    # Define all attribute getters:
    def get_phamerator_genome(self):
        return self.__phamerator_genome
    def get_phagesdb_genome(self):
        return self.__phagesdb_genome
    def get_ncbi_genome(self):
        return self.__ncbi_genome
    def get_phamerator_ncbi_sequence_mismatch(self):
        return self.__phamerator_ncbi_sequence_mismatch
    def get_phamerator_ncbi_perfect_matched_features(self):
        return self.__phamerator_ncbi_perfect_matched_features
    def get_phamerator_ncbi_imperfect_matched_features(self):
        return self.__phamerator_ncbi_imperfect_matched_features
    def get_phamerator_features_unmatched_in_ncbi(self):
        return self.__phamerator_features_unmatched_in_ncbi
    def get_ncbi_features_unmatched_in_phamerator(self):
        return self.__ncbi_features_unmatched_in_phamerator
    def get_phamerator_ncbi_perfect_matched_features_tally(self):
        return self.__phamerator_ncbi_perfect_matched_features_tally
    def get_phamerator_ncbi_imperfect_matched_features_tally(self):
        return self.__phamerator_ncbi_imperfect_matched_features_tally
    def get_phamerator_features_unmatched_in_ncbi_tally(self):
        return self.__phamerator_features_unmatched_in_ncbi_tally
    def get_ncbi_features_unmatched_in_phamerator_tally(self):
        return self.__ncbi_features_unmatched_in_phamerator_tally
    def get_phamerator_ncbi_sequence_length_mismatch(self):
        return self.__phamerator_ncbi_sequence_length_mismatch
    def get_ncbi_record_header_fields_phage_name_mismatch(self):
        return self.__ncbi_record_header_fields_phage_name_mismatch
    def get_ncbi_host_mismatch(self):
        return self.__ncbi_host_mismatch
    def get_phamerator_ncbi_different_descriptions_tally(self):
        return self.__phamerator_ncbi_different_descriptions_tally
    def get_phamerator_ncbi_different_translation_tally(self):
        return self.__phamerator_ncbi_different_translations_tally

    def get_phagesdb_ncbi_sequence_mismatch(self):
        return self.__phagesdb_ncbi_sequence_mismatch
    def get_phagesdb_ncbi_sequence_length_mismatch(self):
        return self.__phagesdb_ncbi_sequence_length_mismatch

    def get_phamerator_phagesdb_sequence_mismatch(self):
        return self.__phamerator_phagesdb_sequence_mismatch
    def get_phamerator_phagesdb_sequence_length_mismatch(self):
        return self.__phamerator_phagesdb_sequence_length_mismatch
    def get_phamerator_phagesdb_host_mismatch(self):
        return self.__phamerator_phagesdb_host_mismatch
    def get_phamerator_phagesdb_accession_mismatch(self):
        return self.__phamerator_phagesdb_accession_mismatch
    def get_phamerator_phagesdb_cluster_subcluster_mismatch(self):
        return self.__phamerator_phagesdb_cluster_subcluster_mismatch
    def get_total_number_genes_with_errors(self):
        return self.__total_number_genes_with_errors
    def get_contains_errors(self):
        return self.__contains_errors
    def get_ph_ncbi_author_error(self):
        return self.__ph_ncbi_author_error


class MatchedCdsFeatures:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.__phamerator_feature = ''
        self.__ncbi_feature = ''

        #Matched data comparison results
        self.__phamerator_ncbi_different_translations = False #True = there are different translations
        self.__phamerator_ncbi_different_start_sites = False #True = there are different start sites
        self.__phamerator_ncbi_different_descriptions = False #True = there are different gene descriptions

        #Total errors summary
        self.__total_errors = 0



    # Define all attribute setters:
    def set_phamerator_feature(self,value):
        self.__phamerator_feature = value
    def set_ncbi_feature(self,value):
        self.__ncbi_feature = value


    def compare_phamerator_ncbi_cds_features(self):

        if self.__phamerator_feature.get_strand() == 'forward':
            if str(self.__phamerator_feature.get_left_boundary()) != str(self.__ncbi_feature.get_left_boundary()):
                self.__phamerator_ncbi_different_start_sites = True
        elif self.__phamerator_feature.get_strand() == 'reverse':
            if str(self.__phamerator_feature.get_right_boundary()) != str(self.__ncbi_feature.get_right_boundary()):
                self.__phamerator_ncbi_different_start_sites = True
        else:
            pass


        product_description_set = set()
        product_description_set.add(self.__phamerator_feature.get_search_notes())
        product_description_set.add(self.__ncbi_feature.get_search_product_description())


        if len(product_description_set) != 1:
            self.__phamerator_ncbi_different_descriptions = True

        if self.__phamerator_feature.get_translation() != self.__ncbi_feature.get_translation():
            self.__phamerator_ncbi_different_translations = True



        #Compute total errors
        #First add all matched feature errors
        if self.__phamerator_ncbi_different_translations:
            self.__total_errors += 1

        if self.__phamerator_ncbi_different_start_sites:
            self.__total_errors += 1

        if self.__phamerator_ncbi_different_descriptions:
            self.__total_errors += 1

        #Now add all errors from each individual feature
        #You first compute errors for each individual feature.
        #This step is performed here instead of in the mainline code
        #because you need to wait for the feature matching step after the genome matching step
        self.__phamerator_feature.compute_total_cds_errors()
        self.__ncbi_feature.compute_total_cds_errors()
        self.__total_errors += self.__phamerator_feature.get_total_errors()
        self.__total_errors += self.__ncbi_feature.get_total_errors()





    # Define all attribute getters:
    def get_phamerator_feature(self):
        return self.__phamerator_feature
    def get_ncbi_feature(self):
        return self.__ncbi_feature
    def get_phamerator_ncbi_different_start_sites(self):
        return self.__phamerator_ncbi_different_start_sites
    def get_phamerator_ncbi_different_descriptions(self):
        return self.__phamerator_ncbi_different_descriptions
    def get_phamerator_ncbi_different_translations(self):
        return self.__phamerator_ncbi_different_translations
    def get_total_errors(self):
        return self.__total_errors


class DatabaseSummary:

    # Initialize all attributes:
    def __init__(self,matched_genomes_list):

        # Initialize all non-calculated attributes:
        self.__matched_genomes_list = matched_genomes_list

        # Initialize all calculated attributes:

        #MySQL data
        #General genome data
        self.__ph_ncbi_update_flag_tally = 0

        #Genome data checks
        self.__ph_genomes_with_nucleotide_errors_tally = 0
        self.__ph_genomes_with_translation_errors_tally = 0
        self.__ph_genomes_with_boundary_errors_tally = 0
        self.__ph_genomes_with_status_accession_error_tally = 0
        self.__ph_genomes_with_status_description_error_tally = 0

        #Phagesdb data
        #Genome data checks
        self.__pdb_genomes_with_nucleotide_errors_tally = 0

        #NCBI data
        #Genome data checks
        self.__ncbi_genomes_with_nucleotide_errors_tally = 0
        self.__ncbi_genomes_with_translation_errors_tally = 0
        self.__ncbi_genomes_with_boundary_errors_tally = 0
        self.__ncbi_genomes_with_missing_locus_tags_tally = 0
        self.__ncbi_genomes_with_locus_tag_typos_tally = 0
        self.__ncbi_genomes_with_description_field_errors_tally = 0

        #MySQL-phagesdb checks
        self.__ph_pdb_sequence_mismatch_tally = 0
        self.__ph_pdb_sequence_length_mismatch_tally = 0
        self.__ph_pdb_cluster_subcluster_mismatch_tally = 0
        self.__ph_pdb_accession_mismatch_tally = 0
        self.__ph_pdb_host_mismatch_tally = 0

        #MySQL-NCBI checks
        self.__ph_ncbi_sequence_mismatch_tally = 0
        self.__ph_ncbi_sequence_length_mismatch_tally = 0
        self.__ph_ncbi_record_header_phage_mismatch_tally = 0
        self.__ph_ncbi_record_header_host_mismatch_tally = 0
        self.__ph_ncbi_genomes_with_imperfectly_matched_features_tally = 0
        self.__ph_ncbi_genomes_with_unmatched_phamerator_features_tally = 0
        self.__ph_ncbi_genomes_with_unmatched_ncbi_features_tally = 0
        self.__ph_ncbi_genomes_with_different_descriptions_tally = 0
        self.__ph_ncbi_genomes_with_different_translations_tally = 0
        self.__ph_ncbi_genomes_with_author_errors_tally = 0

        #phagesdb-NCBI checks
        self.__pdb_ncbi_sequence_mismatch_tally = 0
        self.__pdb_ncbi_sequence_length_mismatch_tally = 0


        #MySQL feature
        #Gene data checks
        self.__ph_translation_errors_tally = 0
        self.__ph_boundary_errors_tally = 0


        #NCBI feature
        #Gene data checks
        self.__ncbi_translation_errors_tally = 0
        self.__ncbi_boundary_errors_tally = 0
        self.__ncbi_missing_locus_tags_tally = 0
        self.__ncbi_locus_tag_typos_tally = 0
        self.__ncbi_description_field_errors_tally = 0

        #MySQL-NCBI checks
        self.__ph_ncbi_different_descriptions_tally = 0
        self.__ph_ncbi_different_start_sites_tally = 0
        self.__ph_ncbi_different_translation_tally = 0
        self.__ph_ncbi_unmatched_phamerator_features_tally = 0
        self.__ph_ncbi_unmatched_ncbi_features_tally = 0


        #Calculate summary metrics
        self.__ph_total_genomes_analyzed = 0
        self.__ph_genomes_unmatched_to_pdb_tally = 0
        self.__ph_genomes_unmatched_to_ncbi_tally = 0
        self.__total_genomes_with_errors = 0


    #Define setter functions
    def compute_summary(self):
        for matched_genomes in self.__matched_genomes_list:

            self.__ph_total_genomes_analyzed += 1
            ph_genome = matched_genomes.get_phamerator_genome()
            pdb_genome = matched_genomes.get_phagesdb_genome()
            ncbi_genome = matched_genomes.get_ncbi_genome()

            if matched_genomes.get_contains_errors():
                self.__total_genomes_with_errors += 1


            #MySQL data
            if isinstance(ph_genome,PhameratorGenome):

                self.__ph_ncbi_update_flag_tally += ph_genome.get_ncbi_update_flag()

                #Genome data checks
                if ph_genome.get_nucleotide_errors():
                    self.__ph_genomes_with_nucleotide_errors_tally += 1
                if ph_genome.get_cds_features_with_translation_error_tally() > 0:
                    self.__ph_genomes_with_translation_errors_tally += 1
                    self.__ph_translation_errors_tally += ph_genome.get_cds_features_with_translation_error_tally()
                if ph_genome.get_cds_features_boundary_error_tally() > 0:
                    self.__ph_genomes_with_boundary_errors_tally += 1
                    self.__ph_boundary_errors_tally += ph_genome.get_cds_features_boundary_error_tally()
                if ph_genome.get_status_accession_error():
                    self.__ph_genomes_with_status_accession_error_tally += 1
                if ph_genome.get_status_description_error():
                    self.__ph_genomes_with_status_description_error_tally += 1

            #Phagesdb data
            if isinstance(pdb_genome,PhagesdbGenome):

                #Genome data checks
                if pdb_genome.get_nucleotide_errors():
                    self.__pdb_genomes_with_nucleotide_errors_tally += 1
            else:
                self.__ph_genomes_unmatched_to_pdb_tally += 1

            #NCBI data
            if isinstance(ncbi_genome,NcbiGenome):

                #Genome data checks
                if ncbi_genome.get_nucleotide_errors():
                    self.__ncbi_genomes_with_nucleotide_errors_tally += 1
                if ncbi_genome.get_cds_features_with_translation_error_tally() > 0:
                    self.__ncbi_genomes_with_translation_errors_tally += 1
                    self.__ncbi_translation_errors_tally += ncbi_genome.get_cds_features_with_translation_error_tally()
                if ncbi_genome.get_cds_features_boundary_error_tally() > 0:
                    self.__ncbi_genomes_with_boundary_errors_tally += 1
                    self.__ncbi_boundary_errors_tally += ncbi_genome.get_cds_features_boundary_error_tally()
                if ncbi_genome.get_missing_locus_tags_tally() > 0:
                    self.__ncbi_genomes_with_missing_locus_tags_tally += 1
                    self.__ncbi_missing_locus_tags_tally += ncbi_genome.get_missing_locus_tags_tally()
                if ncbi_genome.get_locus_tag_typos_tally() > 0:
                    self.__ncbi_genomes_with_locus_tag_typos_tally += 1
                    self.__ncbi_locus_tag_typos_tally += ncbi_genome.get_locus_tag_typos_tally()
                if ncbi_genome.get_description_field_error_tally() > 0:
                    self.__ncbi_genomes_with_description_field_errors_tally += 1
                    self.__ncbi_description_field_errors_tally += ncbi_genome.get_description_field_error_tally()

            else:
                self.__ph_genomes_unmatched_to_ncbi_tally += 1

            #MySQL-phagesdb checks
            if matched_genomes.get_phamerator_phagesdb_sequence_mismatch():
                self.__ph_pdb_sequence_mismatch_tally += 1
            if matched_genomes.get_phamerator_phagesdb_sequence_length_mismatch():
                self.__ph_pdb_sequence_length_mismatch_tally += 1
            if matched_genomes.get_phamerator_phagesdb_cluster_subcluster_mismatch():
                self.__ph_pdb_cluster_subcluster_mismatch_tally += 1
            if matched_genomes.get_phamerator_phagesdb_accession_mismatch():
                self.__ph_pdb_accession_mismatch_tally += 1
            if matched_genomes.get_phamerator_phagesdb_host_mismatch():
                self.__ph_pdb_host_mismatch_tally += 1

            #MySQL-NCBI checks
            if matched_genomes.get_phamerator_ncbi_sequence_mismatch():
                self.__ph_ncbi_sequence_mismatch_tally += 1
            if matched_genomes.get_phamerator_ncbi_sequence_length_mismatch():
                self.__ph_ncbi_sequence_length_mismatch_tally += 1
            if matched_genomes.get_ncbi_record_header_fields_phage_name_mismatch():
                self.__ph_ncbi_record_header_phage_mismatch_tally += 1
            if matched_genomes.get_ncbi_host_mismatch():
                self.__ph_ncbi_record_header_host_mismatch_tally += 1
            if matched_genomes.get_phamerator_ncbi_imperfect_matched_features_tally() > 0:
                self.__ph_ncbi_genomes_with_imperfectly_matched_features_tally += 1
                self.__ph_ncbi_different_start_sites_tally += matched_genomes.get_phamerator_ncbi_imperfect_matched_features_tally()
            if matched_genomes.get_phamerator_features_unmatched_in_ncbi_tally() > 0:
                self.__ph_ncbi_genomes_with_unmatched_phamerator_features_tally += 1
                self.__ph_ncbi_unmatched_phamerator_features_tally += matched_genomes.get_phamerator_features_unmatched_in_ncbi_tally()
            if matched_genomes.get_ncbi_features_unmatched_in_phamerator_tally() > 0:
                self.__ph_ncbi_genomes_with_unmatched_ncbi_features_tally += 1
                self.__ph_ncbi_unmatched_ncbi_features_tally += matched_genomes.get_ncbi_features_unmatched_in_phamerator_tally()
            if matched_genomes.get_phamerator_ncbi_different_descriptions_tally() > 0:
                self.__ph_ncbi_genomes_with_different_descriptions_tally += 1
                self.__ph_ncbi_different_descriptions_tally += matched_genomes.get_phamerator_ncbi_different_descriptions_tally()
            if matched_genomes.get_phamerator_ncbi_different_translation_tally() > 0:
                self.__ph_ncbi_genomes_with_different_translations_tally += 1
                self.__ph_ncbi_different_translation_tally += matched_genomes.get_phamerator_ncbi_different_translation_tally()
            if matched_genomes.get_ph_ncbi_author_error():
                self.__ph_ncbi_genomes_with_author_errors_tally += 1


            #phagesdb-NCBI checks
            if matched_genomes.get_phagesdb_ncbi_sequence_mismatch():
                self.__pdb_ncbi_sequence_mismatch_tally += 1
            if matched_genomes.get_phagesdb_ncbi_sequence_length_mismatch():
                self.__pdb_ncbi_sequence_length_mismatch_tally += 1



    #Define getter functions
    def get_matched_genomes_list(self):
        return self.__matched_genomes_list

    #MySQL data
    #General genome data
    def get_ph_ncbi_update_flag_tally(self):
        return self.__ph_ncbi_update_flag_tally

    #Genome data checks
    def get_ph_genomes_with_nucleotide_errors_tally(self):
        return self.__ph_genomes_with_nucleotide_errors_tally
    def get_ph_genomes_with_translation_errors_tally(self):
        return self.__ph_genomes_with_translation_errors_tally
    def get_ph_genomes_with_boundary_errors_tally(self):
        return self.__ph_genomes_with_boundary_errors_tally
    def get_ph_genomes_with_status_accession_error_tally(self):
        return self.__ph_genomes_with_status_accession_error_tally
    def get_ph_genomes_with_status_description_error_tally(self):
        return self.__ph_genomes_with_status_description_error_tally

    #Phagesdb data
    #Genome data checks
    def get_pdb_genomes_with_nucleotide_errors_tally(self):
        return self.__pdb_genomes_with_nucleotide_errors_tally

    #NCBI data
    #Genome data checks
    def get_ncbi_genomes_with_description_field_errors_tally(self):
        return self.__ncbi_genomes_with_description_field_errors_tally

    def get_ncbi_genomes_with_nucleotide_errors_tally(self):
        return self.__ncbi_genomes_with_nucleotide_errors_tally
    def get_ncbi_genomes_with_translation_errors_tally(self):
        return self.__ncbi_genomes_with_translation_errors_tally
    def get_ncbi_genomes_with_boundary_errors_tally(self):
        return self.__ncbi_genomes_with_boundary_errors_tally
    def get_ncbi_genomes_with_missing_locus_tags_tally(self):
        return self.__ncbi_genomes_with_missing_locus_tags_tally
    def get_ncbi_genomes_with_locus_tag_typos_tally(self):
        return self.__ncbi_genomes_with_locus_tag_typos_tally

    #MySQL-phagesdb checks
    def get_ph_pdb_sequence_mismatch_tally(self):
        return self.__ph_pdb_sequence_mismatch_tally
    def get_ph_pdb_sequence_length_mismatch_tally(self):
        return self.__ph_pdb_sequence_length_mismatch_tally
    def get_ph_pdb_cluster_subcluster_mismatch_tally(self):
        return self.__ph_pdb_cluster_subcluster_mismatch_tally
    def get_ph_pdb_accession_mismatch_tally(self):
        return self.__ph_pdb_accession_mismatch_tally
    def get_ph_pdb_host_mismatch_tally(self):
        return self.__ph_pdb_host_mismatch_tally

    #MySQL-NCBI checks
    def get_ph_ncbi_sequence_mismatch_tally(self):
        return self.__ph_ncbi_sequence_mismatch_tally
    def get_ph_ncbi_sequence_length_mismatch_tally(self):
        return self.__ph_ncbi_sequence_length_mismatch_tally
    def get_ph_ncbi_record_header_phage_mismatch_tally(self):
        return self.__ph_ncbi_record_header_phage_mismatch_tally
    def get_ph_ncbi_record_header_host_mismatch_tally(self):
        return self.__ph_ncbi_record_header_host_mismatch_tally
    def get_ph_ncbi_genomes_with_imperfectly_matched_features_tally(self):
        return self.__ph_ncbi_genomes_with_imperfectly_matched_features_tally
    def get_ph_ncbi_genomes_with_unmatched_phamerator_features_tally(self):
        return self.__ph_ncbi_genomes_with_unmatched_phamerator_features_tally
    def get_ph_ncbi_genomes_with_unmatched_ncbi_features_tally(self):
        return self.__ph_ncbi_genomes_with_unmatched_ncbi_features_tally
    def get_ph_ncbi_genomes_with_different_descriptions_tally(self):
        return self.__ph_ncbi_genomes_with_different_descriptions_tally
    def get_ph_ncbi_genomes_with_different_translations_tally(self):
        return self.__ph_ncbi_genomes_with_different_translations_tally
    def get_ph_ncbi_genomes_with_author_errors_tally(self):
        return self.__ph_ncbi_genomes_with_author_errors_tally



    #phagesdb-NCBI checks
    def get_pdb_ncbi_sequence_mismatch_tally(self):
        return self.__pdb_ncbi_sequence_mismatch_tally
    def get_pdb_ncbi_sequence_length_mismatch_tally(self):
        return self.__pdb_ncbi_sequence_length_mismatch_tally

    #MySQL feature
    #Gene data checks
    def get_ph_translation_errors_tally(self):
        return self.__ph_translation_errors_tally
    def get_ph_boundary_errors_tally(self):
        return self.__ph_boundary_errors_tally

    #NCBI feature
    #Gene data checks
    def get_ncbi_translation_errors_tally(self):
        return self.__ncbi_translation_errors_tally
    def get_ncbi_boundary_errors_tally(self):
        return self.__ncbi_boundary_errors_tally
    def get_ncbi_missing_locus_tags_tally(self):
        return self.__ncbi_missing_locus_tags_tally
    def get_ncbi_locus_tag_typos_tally(self):
        return self.__ncbi_locus_tag_typos_tally
    def get_ncbi_description_field_errors_tally(self):
        return self.__ncbi_description_field_errors_tally

    #MySQL-NCBI checks
    def get_ph_ncbi_different_descriptions_tally(self):
        return self.__ph_ncbi_different_descriptions_tally
    def get_ph_ncbi_different_start_sites_tally(self):
        return self.__ph_ncbi_different_start_sites_tally
    def get_ph_ncbi_different_translation_tally(self):
        return self.__ph_ncbi_different_translation_tally
    def get_ph_ncbi_unmatched_phamerator_features_tally(self):
        return self.__ph_ncbi_unmatched_phamerator_features_tally
    def get_ph_ncbi_unmatched_ncbi_features_tally(self):
        return self.__ph_ncbi_unmatched_ncbi_features_tally

    #Database summaries
    def get_ph_total_genomes_analyzed(self):
        return self.__ph_total_genomes_analyzed
    def get_ph_genomes_unmatched_to_pdb_tally(self):
        return self.__ph_genomes_unmatched_to_pdb_tally
    def get_ph_genomes_unmatched_to_ncbi_tally(self):
        return self.__ph_genomes_unmatched_to_ncbi_tally

    def get_total_genomes_with_errors(self):
        return self.__total_genomes_with_errors

#End of class definitions





































#here
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for comparing databases."""

    COMPARE_HELP = ("Pipeline to compare MySQL, PhagesDB, and "
                    "NCBI databases for inconsistencies.")
    DATABASE_HELP = "Name of the MySQL database from which to compare data."
    OUTPUT_FOLDER_HELP = ("Path to the folder to store results.")
    NCBI_CRED_FILE_HELP = ("Path to the file containing NCBI credentials.")

    parser = argparse.ArgumentParser(description=COMPARE_HELP)
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("output_folder", type=pathlib.Path, help=OUTPUT_FOLDER_HELP)
    parser.add_argument("-c", "--ncbi_credentials_file", type=pathlib.Path,
        help=NCBI_CRED_FILE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args







def main(unparsed_args_list):

    ##Mainline code
    start_time = time.strftime("%x %X")

    #Get the command line parameters
    args = parse_args(unparsed_args_list)
    database = args.database
    output_dir = basic.set_path(args.output_folder, kind="dir", expect=True)
    output_dir = str(output_dir)

    #Get email info for NCBI
    ncbi_cred_dict = ncbi.get_ncbi_creds(args.ncbi_credentials_file)


    #Set up MySQL parameters
    mysqlhost = 'localhost'
    print("\n\n")
    username = getpass.getpass(prompt='mySQL username:')
    print("\n\n")
    password = getpass.getpass(prompt='mySQL password:')
    print("\n\n")




    #Set up dna and protein alphabets to verify sequence integrity
    dna_alphabet_set = set(IUPAC.IUPACUnambiguousDNA.letters)
    protein_alphabet_set = set(IUPAC.ExtendedIUPACProtein.letters)


    #Create output directories
    date = time.strftime("%Y%m%d")

    main_output_folder = '%s_compare' % date
    main_output_path = os.path.join(output_dir,main_output_folder)


    try:
        os.mkdir(main_output_path)
    except:
        print("\nUnable to create output folder: %s" % main_output_path)
        sys.exit(1)

    os.chdir(main_output_path)

















    #Determine which database should be compared: MySQL, phagesdb, NCBI
    analyze_database_options = [\
        'MySQL database',\
        'MySQL database and phagesdb',\
        'MySQL database and NCBI',\
        'MySQL database, phagesdb, and NCBI']
    print('\n\nThe following databases can be compared:')
    print('0: ' + analyze_database_options[0])
    print('1: ' + analyze_database_options[1])
    print('2: ' + analyze_database_options[2])
    print('3: ' + analyze_database_options[3])
    analyze_database = select_option(\
        "\nWhich databases do you want to compare? ", \
        set([0,1,2,3]))

    analyze_database_output = \
        'Databases compared: ' + \
        analyze_database_options[analyze_database]



    valid_database_set = set()
    if analyze_database == 1:
        valid_database_set.add('phagesdb')
    elif analyze_database == 2:
        valid_database_set.add('ncbi')
    elif analyze_database == 3:
        valid_database_set.add('phagesdb')
        valid_database_set.add('ncbi')
    else:
        #By default, all MySQL genomes are checked, \
        #so 'else' results in only MySQL genomes checked
        pass





    #Determine if fasta files should be saved.
    save_phamerator_records = select_option(\
        "\n\nDo you want to save retrieved MySQL records to disk? (yes or no) ", \
        set(['yes','y','no','n']))


    #Create a folder to store MySQL records
    phamerator_output_folder = '%s_mysql_records' % date
    phamerator_output_path = os.path.join(main_output_path,phamerator_output_folder)
    os.mkdir(phamerator_output_path)






    #Set up phagesdb fasta file folder if selected by user
    if 'phagesdb' in valid_database_set:


        #Determine if fasta files should be saved.
        save_phagesdb_records = select_option(\
            "\n\nDo you want to save retrieved phagesdb records to disk? (yes or no) ", \
            set(['yes','y','no','n']))


        #Create a folder to store phagesdb records
        phagesdb_output_folder = '%s_phagesdb_records' % date
        phagesdb_output_path = os.path.join(main_output_path,phagesdb_output_folder)
        os.mkdir(phagesdb_output_path)




    #Set up NCBI parameters if selected by user
    if 'ncbi' in valid_database_set:


        batch_size = ''
        batch_size_valid = False
        while batch_size_valid == False:
            batch_size = input('\nRecord retrieval batch size (must be greater than 0 and recommended is 100-200): ')
            print("\n\n")
            if batch_size.isdigit():
                batch_size = int(batch_size)
                if batch_size > 0:
                    batch_size_valid = True
                else:
                    print('Invalid choice.')
                    print("\n\n")
            else:
                print('Invalid choice.')
                print("\n\n")

        #Determine if NCBI records should be saved
        save_ncbi_records = select_option(\
            "\n\nDo you want to save retrieved NCBI records to disk? (yes or no) ", \
            set(['yes','y','no','n']))


        #Create a folder to store NCBI records
        ncbi_output_folder = '%s_ncbi_records' % date
        ncbi_output_path = os.path.join(main_output_path,ncbi_output_folder)
        os.mkdir(ncbi_output_path)




    #Determine which type of genomes should be checked based on
    #who annotated the genome: Hatfull or Gbk authors
    analyze_genome_author_options = [\
        'none',\
        'Hatfull only',\
        'Genbank only',\
        'Hatfull and Genbank']
    print('\n\nThe following types of genomes based on authorship can be analyzed:')
    #print('0: ' + analyze_genome_author_options[0])
    print('1: ' + analyze_genome_author_options[1])
    print('2: ' + analyze_genome_author_options[2])
    print('3: ' + analyze_genome_author_options[3])
    analyze_genome_author = select_option(\
        "\nWhich type of analysis do you want? ", \
        set([1,2,3]))

    analyze_genome_author_output = \
        'Annotation author analyzed: ' + \
        analyze_genome_author_options[analyze_genome_author]


    valid_genome_author_set = set()
    if analyze_genome_author == 1:
        valid_genome_author_set.add(1)
    elif analyze_genome_author == 2:
        valid_genome_author_set.add(0)
    elif analyze_genome_author == 3:
        valid_genome_author_set.add(1)
        valid_genome_author_set.add(0)
    else:
        pass





    #Determine which types of genomes should be checked based on
    #the status of the annotations: draft, final, unknown
    analyze_genome_status_options = [\
        'none',\
        'draft only',\
        'final only',\
        'unknown only',\
        'draft and final',\
        'draft and unknown',\
        'final and unknown',\
        'draft, final, and unknown']
    print('\n\nThe following types of genomes based on annotation status can be analyzed:')
    #print('0: ' + analyze_genome_status_options[0])
    print('1: ' + analyze_genome_status_options[1])
    print('2: ' + analyze_genome_status_options[2])
    print('3: ' + analyze_genome_status_options[3])
    print('4: ' + analyze_genome_status_options[4])
    print('5: ' + analyze_genome_status_options[5])
    print('6: ' + analyze_genome_status_options[6])
    print('7: ' + analyze_genome_status_options[7])
    analyze_genome_status = select_option(\
        "\nWhich type of analysis do you want? ", \
        set([1,2,3,4,5,6,7]))

    analyze_genome_status_output = \
        'Annotation status analyzed: ' + \
        analyze_genome_status_options[analyze_genome_status]


    valid_genome_status_set = set()
    if analyze_genome_status == 1:
        valid_genome_status_set.add('draft')
    elif analyze_genome_status == 2:
        valid_genome_status_set.add('final')
    elif analyze_genome_status == 3:
        valid_genome_status_set.add('unknown')
    elif analyze_genome_status == 4:
        valid_genome_status_set.add('draft')
        valid_genome_status_set.add('final')
    elif analyze_genome_status == 5:
        valid_genome_status_set.add('draft')
        valid_genome_status_set.add('unknown')
    elif analyze_genome_status == 6:
        valid_genome_status_set.add('final')
        valid_genome_status_set.add('unknown')
    elif analyze_genome_status == 7:
        valid_genome_status_set.add('draft')
        valid_genome_status_set.add('final')
        valid_genome_status_set.add('unknown')
    else:
        pass







    #Create list to hold all open file handles
    file_handle_list = []


    #Retrieve database version
    #Retrieve current genome data in database
    #0 = PhageID
    #1 = Name
    #2 = HostGenus
    #3 = Sequence
    #4 = Length
    #5 = annotation status
    #6 = Cluster
    #7 = Accession
    #8 = auto-update NCBI record flag
    #9 = DateLastModified
    #10 = annotation authorship

    #Retrieve current gene data in database
    #0 = PhageID
    #1 = GeneID
    #2 = Name
    #3 = Start
    #4 = Stop
    #5 = Orientation
    #6 = Translation
    #7 = Notes

    try:
        con = pms.connect("localhost", username, password, database)
        # con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()
    except pms.err.Error as err:
        print("Error connecting to MySQL database")
        print("Error {}: {}".format(err.args[0], err.args[1]))
        sys.exit(1)

    try:
        cur.execute("START TRANSACTION")
        cur.execute("SELECT Version FROM version")
        ph_version = str(cur.fetchone()[0])
        cur.execute("SELECT PhageID,Name,HostGenus,Sequence,Length,\
                            Status,Cluster,Accession,RetrieveRecord,\
                            DateLastModified,AnnotationAuthor FROM phage")
        ph_genome_data_tuples = cur.fetchall()
        cur.execute("SELECT PhageID,GeneID,Name,Start,Stop,Orientation,Translation,Notes from gene")
        ph_gene_data_tuples = cur.fetchall()
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)

    except:
        mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")

    con.close()



    print('\n\nPreparing genome data sets from the MySQL database...')
    ph_genome_count = 1
    processed_ph_genome_count = 0 #May not be equal to ph_genome_count
    ph_genome_object_dict = {} #Key = PhageID; #Value = genome_object
    ph_search_name_set = set()
    ph_search_name_duplicate_set = set()
    ph_accession_set = set()
    ph_accession_duplicate_set = set()
    ph_total_genome_count = len(ph_genome_data_tuples)

    #Iterate through each MySQL genome and create a genome object
    for genome_tuple in ph_genome_data_tuples:
        print("Processing MySQL genome %s out of %s" \
                %(ph_genome_count,ph_total_genome_count))

        #Check to see if the genome has a user-selected status and authorship
        if genome_tuple[5] not in valid_genome_status_set or \
            genome_tuple[10] not in valid_genome_author_set:
            continue
        else:
            genome_object = PhameratorGenome()
            genome_object.set_phage_id(genome_tuple[0])
            genome_object.set_phage_name(genome_tuple[1])
            genome_object.set_host(genome_tuple[2])
            genome_object.set_sequence(genome_tuple[3].decode("utf-8"))
            genome_object.set_accession(genome_tuple[7]) #Set accession before status
            genome_object.set_status(genome_tuple[5])
            genome_object.set_cluster_subcluster(genome_tuple[6])
            genome_object.set_ncbi_update_flag(genome_tuple[8])
            genome_object.set_date_last_modified(genome_tuple[9])
            genome_object.set_annotation_author(genome_tuple[10])
            genome_object.compute_nucleotide_errors(dna_alphabet_set)
            ph_genome_object_dict[genome_tuple[0]] = genome_object

            #This keeps track of whether there are duplicate phage names that will be used
            #to match up to phagesdb data.
            if genome_object.get_search_name() in ph_search_name_set:
                ph_search_name_duplicate_set.add(genome_object.get_search_name())
            else:
                ph_search_name_set.add(genome_object.get_search_name())

            #This keeps track of whether there are duplicate accession numbers that will be
            #used to match up to NCBI data
            if genome_object.get_accession() != '':
                if genome_object.get_accession() in ph_accession_set:
                    ph_accession_duplicate_set.add(genome_object.get_accession())
                else:
                    ph_accession_set.add(genome_object.get_accession())
            ph_genome_count += 1
            processed_ph_genome_count += 1
            #Currently, the processed_ph_genome_count variable does not appear to be used.


            #If selected by user, save retrieved record to file
            if save_phamerator_records == 'yes':

                #To output a fasta file, a Biopython SeqRecord must be created first
                phamerator_fasta_seqrecord = SeqRecord(Seq(genome_object.get_sequence()),\
                                            id=genome_object.get_search_name(),\
                                            description='')

                #Create unique filename
                phamerator_fasta_seqrecord_path = create_unique_filename(\
                                                    phamerator_output_path,
                                                    genome_object.get_search_name(),
                                                    '.fasta')
                #Output the fasta file
                SeqIO.write(phamerator_fasta_seqrecord,\
                                    phamerator_fasta_seqrecord_path,\
                                    'fasta')







    #phagesdb relies on the phageName, and not the phageID. But the MySQL database does not require phageName values to be unique.
    #Check if there are any phageName duplications. If there are, they will not be able to be compared to phagesdb data.
    #Output duplicate phage search names to file
    if len(ph_search_name_duplicate_set) > 0:
        print('Warning: There are duplicate phage search names in the MySQL database.')
        print('Some MySQL genomes will not be able to be matched to phagesdb.')
        output_to_file(list(ph_search_name_duplicate_set),\
                        'duplicate_mysql_phage_names.csv',\
                        analyze_genome_status_output,\
                        database + '_v' + ph_version,\
                        analyze_genome_author_output)
        input('Press ENTER to proceed')

    #Accessions aren't required to be unique in the MySQL database (but they should be), so there could be duplicates
    #Output duplicate names to file
    if len(ph_accession_duplicate_set) > 0:
        print('Warning: There are duplicate accessions in the MySQL database.')
        print('Some MySQL genomes will not be able to be matched to NCBI records.')
        output_to_file(list(ph_accession_duplicate_set),\
                        'duplicate_mysql_phage_accessions.csv',\
                        analyze_genome_status_output,\
                        database + '_v' + ph_version,\
                        analyze_genome_author_output)
        input('Press ENTER to proceed')




    print('\n\nPreparing gene data sets from the MySQL database...')
    ph_gene_count = 1
    ph_gene_objects_list = []
    ph_gene_data_phage_id_set = set()
    ph_total_gene_count = len(ph_gene_data_tuples)

    for gene_tuple in ph_gene_data_tuples:
        print("Processing MySQL gene %s of %s" %(ph_gene_count,ph_total_gene_count))
        ph_gene_data_phage_id_set.add(gene_tuple[0])

        gene_object = PhameratorCdsFeature()
        gene_object.set_phage_id(gene_tuple[0])
        gene_object.set_gene_id(gene_tuple[1])
        gene_object.set_gene_name(gene_tuple[2])
        gene_object.set_type_id('CDS')
        gene_object.set_left_boundary(gene_tuple[3])
        gene_object.set_right_boundary(gene_tuple[4])
        gene_object.set_strand(gene_tuple[5])
        gene_object.set_translation(gene_tuple[6])

        gene_object_description, gene_object_search_description = retrieve_description(gene_tuple[7].decode("utf-8"))
        gene_object.set_notes(gene_object_description, gene_object_search_description)


        gene_object.compute_amino_acid_errors(protein_alphabet_set)
        gene_object.set_start_end_strand_id()
        gene_object.compute_boundary_error()

        ph_gene_objects_list.append(gene_object)
        ph_gene_count += 1


    #ph_gene_objects_dict = {} #Key = PhageID; #Value = list of gene objects
    print('Matching gene and genome data sets from the MySQL database...')
    ph_match_count = 1
    ph_total_match_count = len(ph_genome_object_dict.keys())
    for phage_id in ph_genome_object_dict.keys():
        print("Matching gene data for genome %s of %s" %(ph_match_count,ph_total_match_count))
        genome_object = ph_genome_object_dict[phage_id]
        new_gene_object_list = []
        for gene_object in ph_gene_objects_list:
            if gene_object.get_phage_id() == phage_id:
                new_gene_object_list.append(gene_object)

        genome_object.set_cds_features(new_gene_object_list)
        genome_object.compute_cds_feature_errors()
        genome_object.compute_status_description_error()
        ph_match_count += 1











    #Now retrieve all phagesdb data

    #Data for each phage is stored in a dictionary per phage, and all dictionaries are stored in a list under "results"
    pdb_genome_count = 1
    pdb_genome_dict = {}
    pdb_search_name_set = set()
    pdb_search_name_duplicate_set = set()



    if 'phagesdb' in valid_database_set:
        #Retrieve a list of all sequenced phages listed on phagesdb
        #You have to specify how many results to return at once. If you set it to 1 page long and 100,000 genomes/page, then this will return everything
        print('\n\nRetrieving data from phagesdb...')


        pdb_sequenced_phages_url = "http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000"
        # pdb_sequenced_phages_json = urllib.urlopen(pdb_sequenced_phages_url)
        pdb_sequenced_phages_json = urllib.request.urlopen(pdb_sequenced_phages_url)
        pdb_sequenced_phages_dict = json.loads(pdb_sequenced_phages_json.read())
        pdb_sequenced_phages_json.close()



        pdb_total_genome_count = len(pdb_sequenced_phages_dict['results'])

        for element_dict in pdb_sequenced_phages_dict['results']:
            print("Processing phagesdb genome %s of %s" %(pdb_genome_count,pdb_total_genome_count))

            genome_object = PhagesdbGenome()

            #Name, Host, Accession
            genome_object.set_phage_name(element_dict['phage_name'])
            genome_object.set_host(element_dict['isolation_host']['genus'])
            genome_object.set_accession(element_dict['genbank_accession'])

            #Cluster
            if element_dict['pcluster'] is not None:
                #Sometimes cluster information is not present. In the phagesdb database, it is recorded as NULL.
                #When phages data is downloaded from phagesdb, NULL cluster data is converted to "Unclustered".
                #In these cases, leave cluster as ''
                genome_object.set_cluster(element_dict['pcluster']['cluster'])

            #Subcluster
            if element_dict['psubcluster'] is not None:
                #A phage may be clustered but not subclustered.
                #In these cases, leave subcluster as ''
                genome_object.set_subcluster(element_dict['psubcluster']['subcluster'])

            #Check to see if there is a fasta file stored on phagesdb for this phage
            if element_dict['fasta_file'] is not None:
                fastafile_url = element_dict['fasta_file']

                # response = urllib2.urlopen(fastafile_url)
                # retrieved_fasta_file = response.read()
                # response.close()

                request = urllib.request.Request(fastafile_url)
                with urllib.request.urlopen(request) as response:
                    retrieved_fasta_file = response.read()
                    retrieved_fasta_file = retrieved_fasta_file.decode("utf-8")

                #All sequence rows in the fasta file may not have equal widths, so some processing of the data is required
                #If you split by newline, the header is retained in the first list element
                split_fasta_data = retrieved_fasta_file.split('\n')
                pdb_sequence = ''
                index = 1
                while index < len(split_fasta_data):
                    pdb_sequence = pdb_sequence + split_fasta_data[index].strip() #Strip off potential whitespace before appending, such as '\r'
                    index += 1
                genome_object.set_sequence(pdb_sequence)
                genome_object.compute_nucleotide_errors(dna_alphabet_set)

            pdb_search_name = genome_object.get_search_name()
            if pdb_search_name in pdb_search_name_set:
                pdb_search_name_duplicate_set.add(pdb_search_name)
            else:
                pdb_search_name_set.add(pdb_search_name)
                pdb_genome_dict[pdb_search_name] = genome_object
            pdb_genome_count += 1



            #If selected by user, save retrieved record to file
            if save_phagesdb_records == 'yes':

                #To output a fasta file, a Biopython SeqRecord must be created first
                phagesdb_fasta_seqrecord = SeqRecord(Seq(genome_object.get_sequence()),\
                                            id=genome_object.get_search_name(),\
                                            description='')

                #Create unique filename
                phagesdb_fasta_seqrecord_path = create_unique_filename(\
                                                    phagesdb_output_path,
                                                    genome_object.get_search_name(),
                                                    '.fasta')

                #Output the fasta file
                SeqIO.write(phagesdb_fasta_seqrecord,\
                                    phagesdb_fasta_seqrecord_path,\
                                    'fasta')



        #phagesdb phage names are unique, but just make sure after they are converted to a search name
        if len(pdb_search_name_duplicate_set) > 0:

            print('Warning: There are duplicate phage search names in phagesdb.')
            print('Some phagesdb genomes will not be able to be matched to genomes in MySQL.')
            output_to_file(list(pdb_search_name_duplicate_set),\
                            'duplicate_phagesdb_phage_names.csv',\
                            analyze_genome_status_output,\
                            database + '_v' + ph_version,\
                            analyze_genome_author_output)
            input('Press ENTER to proceed')

        #Make sure all sequenced phage data has been retrieved
        if (len(pdb_sequenced_phages_dict['results']) != pdb_sequenced_phages_dict['count'] or \
            len(pdb_sequenced_phages_dict['results']) != len(pdb_genome_dict)):

            print("\nUnable to retrieve all phage data from phagesdb due to default retrieval parameters.")
            print('Update parameters in script to enable these functions.')
            input('Press ENTER to proceed')










    #Retrieve and parse NCBI records if selected by user
    ncbi_genome_dict = {} #Key = accession; #Value = genome data
    ncbi_genome_count = 1


    if 'ncbi' in valid_database_set:

        print("\n\nRetrieving records from NCBI")

        #Use esearch to verify the accessions are valid and efetch to retrieve the record
        ncbi.set_entrez_credentials(
            tool=ncbi_cred_dict["ncbi_tool"],
            email=ncbi_cred_dict["ncbi_email"],
            api_key=ncbi_cred_dict["ncbi_api_key"])



        #Create batches of accessions
        accession_retrieval_list = list(ph_accession_set)

        #Add [ACCN] field to each accession number
        index = 0
        while index < len(accession_retrieval_list):
            accession_retrieval_list[index] = accession_retrieval_list[index] + "[ACCN]"
            index += 1

        #Keep track of specific records
        retrieved_record_list = []
        retrieval_error_list = []



        #When retrieving in batch sizes, first create the list of values indicating which indices of the accession_retrieval_list should be used to create each batch
        #For instace, if there are five accessions, batch size of two produces indices = 0,2,4
        for batch_index_start in range(0,len(accession_retrieval_list),batch_size):


            if batch_index_start + batch_size > len(accession_retrieval_list):
                batch_index_stop = len(accession_retrieval_list)
            else:
                batch_index_stop = batch_index_start + batch_size

            current_batch_size = batch_index_stop - batch_index_start
            delimiter = " | "
            esearch_term = delimiter.join(accession_retrieval_list[batch_index_start:batch_index_stop])


            #Use esearch for each accession
            search_handle = Entrez.esearch(db = 'nucleotide', term = esearch_term,usehistory = 'y')
            search_record = Entrez.read(search_handle)
            search_count = int(search_record['Count'])
            search_webenv = search_record['WebEnv']
            search_query_key = search_record['QueryKey']



            #Keep track of the accessions that failed to be located in NCBI
            if search_count < current_batch_size:
                search_accession_failure = search_record['ErrorList']['PhraseNotFound']

                #Each element in this list is formatted "accession[ACCN]"
                for element in search_accession_failure:
                    retrieval_error_list.append(element[:-6])



            #Now retrieve all records using efetch
            fetch_handle = Entrez.efetch(db = 'nucleotide', \
                                        rettype = 'gb', \
                                        retmode = 'text', \
                                        retstart = 0, \
                                        retmax = search_count, \
                                        webenv = search_webenv, \
                                        query_key = search_query_key)
            fetch_records = SeqIO.parse(fetch_handle,'genbank')

            for record in fetch_records:
                retrieved_record_list.append(record)

            search_handle.close()
            fetch_handle.close()



        #Report the accessions that could not be retrieved.
        failed_accession_report_fh = open(os.path.join(main_output_path,date + '_failed_accession_retrieval.csv'), 'w')
        failed_accession_report_writer = csv.writer(failed_accession_report_fh)
        failed_accession_report_writer.writerow([date + ' Database comparison'])
        failed_accession_report_writer.writerow([database + '_v' + ph_version])
        failed_accession_report_writer.writerow([analyze_genome_author_output])
        failed_accession_report_writer.writerow([analyze_genome_status_output])
        failed_accession_report_writer.writerow([analyze_database_output])
        failed_accession_report_writer.writerow(['Accessions unable to be retrieved from NCBI:'])
        for retrieval_error_accession in retrieval_error_list:
            failed_accession_report_writer.writerow([retrieval_error_accession])
        failed_accession_report_fh.close()




        ncbi_total_genome_count = len(retrieved_record_list)
        for retrieved_record in retrieved_record_list:
            print("Processing NCBI genome %s of %s" %(ncbi_genome_count,ncbi_total_genome_count))

            genome_object = NcbiGenome()

            try:
                genome_object.set_record_name(retrieved_record.name)
            except:
                genome_object.set_record_name('')

            try:
                genome_object.set_record_id(retrieved_record.id)
            except:
                genome_object.set_record_id('')

            try:
                #There may be a list of accessions associated with this file. I think the first accession in the list is the most recent.
                #Discard the version suffix if it is present in the Accession field (it might not be present)
                record_accession = retrieved_record.annotations['accessions'][0]
                record_accession = record_accession.split('.')[0]
                genome_object.set_record_accession(record_accession)
            except:
                genome_object.set_record_accession('')

            try:
                genome_object.set_record_description(retrieved_record.description)
            except:
                genome_object.set_record_description('')

            try:
                genome_object.set_record_source(retrieved_record.annotations['source'])
            except:
                genome_object.set_record_source('')

            try:
                record_organism = retrieved_record.annotations['organism']
                genome_object.set_record_organism(record_organism)

                #Only truncate organism name for the 'phage name' field
                if record_organism.split(' ')[-1] == 'Unclassified.':
                    genome_object.set_phage_name(record_organism.split(' ')[-2])
                else:
                    genome_object.set_phage_name(record_organism.split(' ')[-1])

            except:
                genome_object.set_record_organism('')


            try:
                #The retrieved authors can be stored in multiple Reference elements
                record_references = retrieved_record.annotations['references']
                record_references_author_list = []
                for reference in record_references:
                    record_references_author_list.append(reference.authors)
                record_author_string = ';'.join(record_references_author_list)
                genome_object.set_record_authors(record_author_string)
            except:
                genome_object.set_record_authors('')


            #Nucleotide sequence and errors
            genome_object.set_sequence(retrieved_record.seq)
            genome_object.compute_nucleotide_errors(dna_alphabet_set)


            #Iterate through all features
            source_feature_list = []
            ncbi_cds_features = []

            #A good bit of the code for parsing features is copied from import_phage.py
            for feature in retrieved_record.features:

                gene_object = NcbiCdsFeature()


                if feature.type != 'CDS':
                    #Retrieve the Source Feature info
                    if feature.type == 'source':
                        source_feature_list.append(feature)

                else:

                    #Feature type
                    gene_object.set_type_id('CDS')

                    #Locus tag
                    try:
                        gene_object.set_locus_tag(feature.qualifiers['locus_tag'][0])
                    except:
                        gene_object.set_locus_tag('')




                    #Orientation
                    if feature.strand == 1:
                        gene_object.set_strand('forward')
                    elif feature.strand == -1:
                        gene_object.set_strand('reverse')
                    #ssRNA phages
                    elif feature.strand is None:
                        gene_object.set_strand('NA')


                    #Gene boundary coordinates
                    #Compound features are tricky to parse.
                    if str(feature.location)[:4] == 'join':

                        #Skip this compound feature if it is comprised of more than two features (too tricky to parse).
                        if len(feature.location.parts) <= 2:
                            #Retrieve compound feature positions based on strand
                            if feature.strand == 1:
                                gene_object.set_left_boundary(str(feature.location.parts[0].start))
                                gene_object.set_right_boundary(str(feature.location.parts[1].end))
                            elif feature.strand == -1:
                                gene_object.set_left_boundary(str(feature.location.parts[1].start))
                                gene_object.set_right_boundary(str(feature.location.parts[0].end))
                            #If strand is None, not sure how to parse it
                            else:
                                pass
                    else:
                        gene_object.set_left_boundary(str(feature.location.start))
                        gene_object.set_right_boundary(str(feature.location.end))

                    #Translation
                    try:
                        gene_object.set_translation(feature.qualifiers['translation'][0])
                    except:
                        pass

                    #Gene function, note, and product descriptions
                    try:
                        feature_product,feature_search_product = retrieve_description(feature.qualifiers['product'][0])
                        gene_object.set_product_description(feature_product,feature_search_product)
                    except:
                        pass
                    try:
                        feature_function,feature_search_function = retrieve_description(feature.qualifiers['function'][0])
                        gene_object.set_function_description(feature_function,feature_search_function)
                    except:
                        pass

                    try:
                        feature_note,feature_search_note = retrieve_description(feature.qualifiers['note'][0])
                        gene_object.set_note_description(feature_note,feature_search_note)
                    except:
                        pass

                    #Gene number
                    try:
                        gene_object.set_gene_number(feature.qualifiers['gene'][0])
                    except:
                        pass

                    #Compute other fields
                    gene_object.compute_amino_acid_errors(protein_alphabet_set)
                    gene_object.set_start_end_strand_id()
                    gene_object.compute_boundary_error()
                    gene_object.compute_description_error()

                    #Now add to full list of gene objects
                    ncbi_cds_features.append(gene_object)




            #Set the following variables after iterating through all features

            #If there was one and only one source feature present, parse certain qualifiers
            if len(source_feature_list) == 1:
                try:
                    genome_object.set_source_feature_organism(str(source_feature_list[0].qualifiers['organism'][0]))
                except:
                    pass
                try:
                    genome_object.set_source_feature_host(str(source_feature_list[0].qualifiers['host'][0]))
                except:
                    pass
                try:
                    genome_object.set_source_feature_lab_host(str(source_feature_list[0].qualifiers['lab_host'][0]))
                except:
                    pass

            genome_object.set_cds_features(ncbi_cds_features)
            genome_object.compute_cds_feature_errors()
            genome_object.compute_ncbi_cds_feature_errors()


            #After parsing all data, add to the ncbi dictionary
            ncbi_genome_dict[genome_object.get_record_accession()] = genome_object

            #If selected by user, save retrieved record to file
            if save_ncbi_records == 'yes':

                #Save the file, but use a truncated version of the full organism name
                name_prefix = genome_object.get_record_organism().lower()
                if name_prefix.split(' ')[-1] == 'unclassified.':
                    name_prefix = name_prefix.split(' ')[-2]
                else:
                    name_prefix = name_prefix.split(' ')[-1]


                #Create unique filename
                ncbi_gb_seqrecord_path = create_unique_filename(\
                        ncbi_output_path,\
                        name_prefix + '__' + genome_object.get_record_accession(),\
                        '.gb')

                #Output the fasta file
                SeqIO.write(retrieved_record,\
                                    ncbi_gb_seqrecord_path,\
                                    'genbank')



            ncbi_genome_count += 1





    #Now that all NCBI and phagesdb data is retrieved, match up to MySQL genome data
    print("Matching phagesdb and NCBI genomes to MySQL genomes...")
    matched_genomes_list = [] #Will store a list of MatchedGenomes objects
    ph_unmatched_to_pdb_genomes = [] #List of the MySQL genome objects with no phagesdb matches
    ph_unmatched_to_ncbi_genomes = [] #List of the MySQL genome objects with no NCBI matches
    matched_genome_count = 1
    matched_total_genome_count = len(ph_genome_object_dict.keys())
    #Iterate through each phage in the MySQL database
    for phage_id in ph_genome_object_dict.keys():

        print("Matching genome %s of %s" %(matched_genome_count,matched_total_genome_count))

        phamerator_genome = ph_genome_object_dict[phage_id]
        matched_objects = MatchedGenomes()
        matched_objects.set_phamerator_genome(phamerator_genome)

        #Match up phagesdb genome
        #First try to match up the phageID, and if that doesn't work, try to match up the phageName
        if phamerator_genome.get_search_id() in pdb_genome_dict.keys():
            pdb_genome = pdb_genome_dict[phamerator_genome.get_search_id()]
            #Make sure the pdb_genome doesn't have a search name that was duplicated
            if pdb_genome.get_search_name() in pdb_search_name_duplicate_set:
                pdb_genome = ''

        elif phamerator_genome.get_search_name() in pdb_genome_dict.keys() and not \
            phamerator_genome.get_search_name() in ph_search_name_duplicate_set:

            pdb_genome = pdb_genome_dict[phamerator_genome.get_search_name()]
            #Make sure the pdb_genome doesn't have a search name that was duplicated
            if pdb_genome.get_search_name in pdb_search_name_duplicate_set:
                pdb_genome = ''

        else:
            pdb_genome = ''
            ph_unmatched_to_pdb_genomes.append(phamerator_genome)

        matched_objects.set_phagesdb_genome(pdb_genome)

        #Now match up NCBI genome
        if phamerator_genome.get_accession() != '' and \
            phamerator_genome.get_accession() not in ph_accession_duplicate_set:
                #Retrieval of record may have failed
                try:
                    ncbi_genome = ncbi_genome_dict[phamerator_genome.get_accession()]
                except:
                    ncbi_genome = ''
                    ph_unmatched_to_ncbi_genomes.append(phamerator_genome)

        else:
            ncbi_genome = ''
            ph_unmatched_to_ncbi_genomes.append(phamerator_genome)

        matched_objects.set_ncbi_genome(ncbi_genome)
        matched_genomes_list.append(matched_objects)
        matched_genome_count += 1


    #Output unmatched data to file
    unmatched_genome_output_fh = open(os.path.join(main_output_path,date + '_database_comparison_unmatched_genomes.csv'), 'w')
    unmatched_genome_output_writer = csv.writer(unmatched_genome_output_fh)
    unmatched_genome_output_writer.writerow([date + ' Database comparison'])
    unmatched_genome_output_writer.writerow([database + '_v' + ph_version])
    unmatched_genome_output_writer.writerow([analyze_genome_author_output])
    unmatched_genome_output_writer.writerow([analyze_genome_status_output])
    unmatched_genome_output_writer.writerow([analyze_database_output])
    unmatched_genome_output_writer.writerow([''])

    if 'phagesdb' in valid_database_set:

        unmatched_genome_output_writer.writerow(['The following MySQL genomes were not matched to phagesdb:'])
        unmatched_genome_output_writer.writerow(['PhageID','PhageName','Author','Status'])
        for element in ph_unmatched_to_pdb_genomes:
            unmatched_genome_output_writer.writerow([element.get_phage_id(),\
                                                        element.get_phage_name(),\
                                                        element.get_annotation_author(),\
                                                        element.get_status()])

    if 'ncbi' in valid_database_set:
        unmatched_genome_output_writer.writerow([''])
        unmatched_genome_output_writer.writerow(['\nThe following MySQL genomes were not matched to NCBI:'])
        unmatched_genome_output_writer.writerow(['PhageID','PhageName','Author','Status','Accession'])
        for element in ph_unmatched_to_ncbi_genomes:
            unmatched_genome_output_writer.writerow([element.get_phage_id(),\
                                                        element.get_phage_name(),\
                                                        element.get_annotation_author(),\
                                                        element.get_status(),\
                                                        element.get_accession()])
    unmatched_genome_output_fh.close()






    #Now that all genomes have been matched, iterate through each matched objects
    #and run methods to compare the genomes
    print("Comparing matched genomes and identifying inconsistencies...")
    comparison_count = 1
    comparison_total_count = len(matched_genomes_list)
    for matched_genome_object in matched_genomes_list:
        print("Comparison matched genome set %s of %s" %(comparison_count,comparison_total_count))
        matched_genome_object.compare_phamerator_ncbi_genomes() #This method automatically calls method to match and compare cds features
        matched_genome_object.compare_phamerator_phagesdb_genomes()
        matched_genome_object.compare_phagesdb_ncbi_genomes()
        matched_genome_object.compute_total_genome_errors()
        comparison_count += 1

    #Now that all individual matched_genome_objects have all computed attributes,
    #compute database summary
    summary_object = DatabaseSummary(matched_genomes_list)
    summary_object.compute_summary()















    #Now output results
    print("Outputting results to file...")
    #Open files to record update information



    #Create database summary output file
    database_summary_report_fh = open(os.path.join(main_output_path,date + '_database_comparison_summary_output.csv'), 'w')
    file_handle_list.append(database_summary_report_fh)
    database_summary_report_writer = csv.writer(database_summary_report_fh)
    database_summary_report_writer.writerow([date + ' Database comparison'])
    database_summary_report_writer.writerow([database + '_v' + ph_version])
    database_summary_report_writer.writerow([analyze_genome_author_output])
    database_summary_report_writer.writerow([analyze_genome_status_output])
    database_summary_report_writer.writerow([analyze_database_output])

    #Create vector of column headers
    summary_report_fields = [\

        '',\
        'Genome summary:',\

        #Column header
        'Database comparison metric',\


        #Database summaries
        'ph_total_genomes_analyzed',\
        'ph_genomes_unmatched_to_pdb_tally',\
        'ph_genomes_unmatched_to_ncbi_tally',\
        'total_genomes_with_errors',\


        #MySQL data
        #General genome data
        'ph_ncbi_update_flag_tally',\

        #Genome data checks
        'ph_genomes_with_nucleotide_errors_tally',\
        'ph_genomes_with_translation_errors_tally',\
        'ph_genomes_with_boundary_errors_tally',\
        'ph_genomes_with_status_accession_error_tally',\
        'ph_genomes_with_status_description_error_tally',\

        #Phagesdb data
        #Genome data checks
        'pdb_genomes_with_nucleotide_errors_tally',\

        #NCBI data
        #Genome data checks
        'ncbi_genomes_with_description_field_errors_tally',\
        'ncbi_genomes_with_nucleotide_errors_tally',\
        'ncbi_genomes_with_translation_errors_tally',\
        'ncbi_genomes_with_boundary_errors_tally',\
        'ncbi_genomes_with_missing_locus_tags_tally',\
        'ncbi_genomes_with_locus_tag_typos_tally',\

        #MySQL-phagesdb checks
        'ph_pdb_sequence_mismatch_tally',\
        'ph_pdb_sequence_length_mismatch_tally',\
        'ph_pdb_cluster_subcluster_mismatch_tally',\
        'ph_pdb_accession_mismatch_tally',\
        'ph_pdb_host_mismatch_tally',\

        #MySQL-NCBI checks
        'ph_ncbi_sequence_mismatch_tally',\
        'ph_ncbi_sequence_length_mismatch_tally',\
        'ph_ncbi_record_header_phage_mismatch_tally',\
        'ph_ncbi_record_header_host_mismatch_tally',\
        'ph_ncbi_genomes_with_author_errors_tally',\
        'ph_ncbi_genomes_with_imperfectly_matched_features_tally',\
        'ph_ncbi_genomes_with_unmatched_phamerator_features_tally',\
        'ph_ncbi_genomes_with_unmatched_ncbi_features_tally',\
        'ph_ncbi_genomes_with_different_descriptions_tally',\
        'ph_ncbi_genomes_with_different_translations_tally',\

        #phagesdb-NCBI checks
        'pdb_ncbi_sequence_mismatch_tally',\
        'pdb_ncbi_sequence_length_mismatch_tally',\



        #Separate all checks that tally all genes
        '',\
        'Gene summary:',\

        #Column header
        'Database comparison metric',\

        #MySQL feature
        #Gene data checks
        'ph_translation_errors_tally',\
        'ph_boundary_errors_tally',\

        #NCBI feature
        #Gene data checks
        'ncbi_translation_errors_tally',\
        'ncbi_boundary_errors_tally',\
        'ncbi_missing_locus_tags_tally',\
        'ncbi_locus_tag_typos_tally',\
        'ncbi_description_field_errors_tally',\

        #MySQL-NCBI checks
        'ph_ncbi_different_descriptions_tally',\
        'ph_ncbi_different_start_sites_tally',\
        'ph_ncbi_different_translation_tally',\
        'ph_ncbi_unmatched_phamerator_features_tally',\
        'ph_ncbi_unmatched_ncbi_features_tally']






    genome_report_fh = open(os.path.join(main_output_path,date + '_database_comparison_genome_output.csv'), 'w')
    file_handle_list.append(genome_report_fh)
    genome_report_writer = csv.writer(genome_report_fh)
    genome_report_writer.writerow([date + ' Database comparison'])
    genome_report_writer.writerow([database + '_v' + ph_version])
    genome_report_writer.writerow([analyze_genome_author_output])
    genome_report_writer.writerow([analyze_genome_status_output])
    genome_report_writer.writerow([analyze_database_output])



    #Create vector of column headers
    genome_report_column_headers = [\

        'ph_phage_id',\
        'contains_errors',\

        #MySQL
        #General genome data
        'ph_phage_name',\
        'ph_search_id',\
        'ph_search_name',\
        'ph_status',\
        'ph_cluster_subcluster',\
        'ph_host',\
        'ph_accession',\
        'ph_dna_seq_length',\
        'ph_gene_tally',\
        'ph_description_tally',\
        'ph_ncbi_update_flag',\
        'ph_date_last_modified',\
        'ph_annotation_author',\


        #Genome data checks
        'ph_dna_seq_error',\
        'ph_gene_translation_error_tally',\
        'ph_gene_coords_error_tally',\
        'ph_status_description_error',\
        'ph_status_accession_error',\


        #phagesdb
        #General genome data
        'pdb_phage_name',\
        'pdb_search_name',\
        'pdb_cluster',\
        'pdb_subcluster',\
        'pdb_host',\
        'pdb_accession',\
        'pdb_dna_seq_length',\

        #Genome data checks
        'pdb_dna_seq_error',\


        #ncbi
        #General genome data
        'ncbi_phage_name',\
        'ncbi_search_name',\
        'ncbi_record_id',\
        'ncbi_record_name',\
        'ncbi_record_accession',\
        'ncbi_record_definition',\
        'ncbi_record_source',\
        'ncbi_record_organism',\
        'ncbi_source_feature_organism',\
        'ncbi_source_feature_host',\
        'ncbi_source_feature_lab_host',\
        'ncbi_authors',\
        'ncbi_dna_seq_length',\
        'ncbi_gene_tally',\

        #Genome data checks
        'ncbi_dna_seq_error',\
        'ncbi_gene_translation_error_tally',\
        'ncbi_gene_coords_error_tally',\
        'ncbi_gene_product_tally',\
        'ncbi_gene_function_tally',\
        'ncbi_gene_note_tally',\
        'ncbi_missing_locus_tag_tally',\
        'ncbi_locus_tag_typo_tally',\
        'ncbi_description_field_error_tally',\

        #MySQL-phagesdb
        'ph_pdb_dna_seq_error',\
        'ph_pdb_dna_seq_length_error',\
        'ph_pdb_cluster_subcluster_error',\
        'ph_pdb_accession_error',\
        'ph_pdb_host_error',\

        #MySQL-ncbi
        'ph_ncbi_dna_seq_error',\
        'ph_ncbi_dna_seq_length_error',\
        'ph_ncbi_record_header_name_error',\
        'ph_ncbi_record_header_host_error',\


        #Author error is dependent on MySQL genome annotation author and
        #NCBI list of authors, so this metric should be reported with
        #the other ph_ncbi error tallies.
        'ph_ncbi_author_error',\
        'ph_ncbi_perfectly_matched_gene_tally',\
        'ph_ncbi_imperfectly_matched_gene_tally',\
        'ph_ncbi_unmatched_phamerator_gene_tally',\
        'ph_ncbi_unmatched_ncbi_gene_tally',\
        'ph_ncbi_gene_description_error_tally',\
        'ph_ncbi_perfectly_matched_gene_translation_error_tally',\

        #Number of genes with errors is computed slightly differently
        #depending on whethere there are matching MySQL and NCBI genomes.
        #Therefore,this metric should be reported with the other ph_ncbi error tallies
        #even if there is no matching NCBI genome.
        'ph_ncbi_genes_with_errors_tally',\



        #phagesdb-ncbi
        'pdb_ncbi_dna_seq_error',\
        'pdb_ncbi_dna_seq_length_error']
    genome_report_writer.writerow(genome_report_column_headers)





    gene_report_fh = open(os.path.join(main_output_path,date + '_database_comparison_gene_output.csv'), 'w')
    file_handle_list.append(gene_report_fh)
    gene_report_writer = csv.writer(gene_report_fh)
    gene_report_writer.writerow([date + ' Database comparison'])
    gene_report_writer.writerow([database + '_v' + ph_version])
    gene_report_writer.writerow([analyze_genome_author_output])
    gene_report_writer.writerow([analyze_genome_status_output])
    gene_report_writer.writerow([analyze_database_output])

    #Create vector of column headers
    gene_report_column_headers = [\

        #Gene summary
        'ph_search_name',\
        'total_errors',\

        #MySQL
        #General gene data
        'ph_phage_id',\
        'ph_search_id',\
        'ph_type_id',\
        'ph_gene_id',\
        'ph_gene_name',\
        'ph_left_boundary',\
        'ph_right_boundary',\
        'ph_strand',\
        'ph_translation',\
        'ph_translation_length',\
        'ph_gene_notes',\

        #Gene data checks
        'ph_translation_error',\
        'ph_gene_coords_error',\



        #NCBI
        #General gene data
        'ncbi_locus_tag',\
        'ncbi_gene_number',\
        'ncbi_type_id',\
        'ncbi_left_boundary',\
        'ncbi_right_boundary',\
        'ncbi_strand',\
        'ncbi_translation',\
        'ncbi_translation_length',\
        'ncbi_product_description',\
        'ncbi_function_description',\
        'ncbi_note_description',\

        #Gene data checks
        'ncbi_translation_error',\
        'ncbi_gene_coords_error',\
        'ncbi_missing_locus_tag',\
        'ncbi_locus_tag_typo',\
        'ncbi_description_field_error',\

        #MySQL-NCBI checks
        'ph_ncbi_unmatched_error',\
        'ph_ncbi_description_error',\
        'ph_ncbi_start_coordinate_error',\
        'ph_ncbi_translation_error']
    gene_report_writer.writerow(gene_report_column_headers)







    #Output all data to file

    summary_data_output = []
    summary_data_output.append('')
    summary_data_output.append('')
    summary_data_output.append('tally') #Column header

    #First output database summary data
    summary_data_output.append(summary_object.get_ph_total_genomes_analyzed())
    summary_data_output.append(summary_object.get_ph_genomes_unmatched_to_pdb_tally())
    summary_data_output.append(summary_object.get_ph_genomes_unmatched_to_ncbi_tally())
    summary_data_output.append(summary_object.get_total_genomes_with_errors())



    #MySQL data
    #General genome data
    summary_data_output.append(summary_object.get_ph_ncbi_update_flag_tally())

    #Genome data checks
    summary_data_output.append(summary_object.get_ph_genomes_with_nucleotide_errors_tally())
    summary_data_output.append(summary_object.get_ph_genomes_with_translation_errors_tally())
    summary_data_output.append(summary_object.get_ph_genomes_with_boundary_errors_tally())
    summary_data_output.append(summary_object.get_ph_genomes_with_status_accession_error_tally())
    summary_data_output.append(summary_object.get_ph_genomes_with_status_description_error_tally())

    #Phagesdb data
    #Genome data checks
    summary_data_output.append(summary_object.get_pdb_genomes_with_nucleotide_errors_tally())

    #NCBI data
    #Genome data checks
    summary_data_output.append(summary_object.get_ncbi_genomes_with_description_field_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_genomes_with_nucleotide_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_genomes_with_translation_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_genomes_with_boundary_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_genomes_with_missing_locus_tags_tally())
    summary_data_output.append(summary_object.get_ncbi_genomes_with_locus_tag_typos_tally())





    #MySQL-phagesdb checks
    summary_data_output.append(summary_object.get_ph_pdb_sequence_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_pdb_sequence_length_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_pdb_cluster_subcluster_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_pdb_accession_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_pdb_host_mismatch_tally())

    #MySQL-NCBI checks
    summary_data_output.append(summary_object.get_ph_ncbi_sequence_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_sequence_length_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_record_header_phage_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_record_header_host_mismatch_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_author_errors_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_imperfectly_matched_features_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_unmatched_phamerator_features_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_unmatched_ncbi_features_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_different_descriptions_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_genomes_with_different_translations_tally())

    #phagesdb-NCBI checks
    summary_data_output.append(summary_object.get_pdb_ncbi_sequence_mismatch_tally())
    summary_data_output.append(summary_object.get_pdb_ncbi_sequence_length_mismatch_tally())




    #Gene summaries
    summary_data_output.append('')
    summary_data_output.append('')
    summary_data_output.append('tally') #Column header


    #MySQL feature
    #Gene data checks
    summary_data_output.append(summary_object.get_ph_translation_errors_tally())
    summary_data_output.append(summary_object.get_ph_boundary_errors_tally())

    #NCBI feature
    #Gene data checks
    summary_data_output.append(summary_object.get_ncbi_translation_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_boundary_errors_tally())
    summary_data_output.append(summary_object.get_ncbi_missing_locus_tags_tally())
    summary_data_output.append(summary_object.get_ncbi_locus_tag_typos_tally())
    summary_data_output.append(summary_object.get_ncbi_description_field_errors_tally())

    #MySQL-NCBI checks
    summary_data_output.append(summary_object.get_ph_ncbi_different_descriptions_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_different_start_sites_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_different_translation_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_unmatched_phamerator_features_tally())
    summary_data_output.append(summary_object.get_ph_ncbi_unmatched_ncbi_features_tally())




    if len(summary_report_fields) == len(summary_data_output):
        summary_output_index = 0
        #database_summary_report_writer.writerow(['database_comparison_metric','tally']) #Column headers
        while summary_output_index < len(summary_report_fields):
            database_summary_report_writer.writerow([\
                summary_report_fields[summary_output_index],\
                summary_data_output[summary_output_index]])
            summary_output_index += 1
    else:
        database_summary_report_writer.writerow(['Different number of database summary metrics than database summary data fields.'])
        database_summary_report_writer.writerow(['Unable to output database summary.'])



    #Now iterate through matched objects.
    #All MySQL genomes are stored in a MatchedGenomes object, even if there are no phagesdb or NCBI matches.
    #All but a few MySQL phages should be matched to phagesdb
    #Only ~half of MySQL phages should be matched to NCBI
    for matched_genomes in summary_object.get_matched_genomes_list():

        genome_data_output = []
        ph_genome = matched_genomes.get_phamerator_genome()
        pdb_genome = matched_genomes.get_phagesdb_genome()
        ncbi_genome = matched_genomes.get_ncbi_genome()




        #Genome summary data
        genome_data_output.append(ph_genome.get_phage_id())# PhageID
        genome_data_output.append(matched_genomes.get_contains_errors())# any errors detected in this genome?

        #MySQL data
        #General genome data
        genome_data_output.append(ph_genome.get_phage_name())# Name
        genome_data_output.append(ph_genome.get_search_id())# search_id
        genome_data_output.append(ph_genome.get_search_name())# search name
        genome_data_output.append(ph_genome.get_status())# status
        genome_data_output.append(ph_genome.get_cluster_subcluster())# cluster_subcluster
        genome_data_output.append(ph_genome.get_host())# Host
        genome_data_output.append(ph_genome.get_accession())# Accession
        genome_data_output.append(ph_genome.get_length())# sequence_length
        genome_data_output.append(ph_genome.get_cds_features_tally())# # genes
        genome_data_output.append(ph_genome.get_description_tally()) # # genes with descriptions
        genome_data_output.append(ph_genome.get_ncbi_update_flag())# ncbi_update_flag
        genome_data_output.append(ph_genome.get_date_last_modified()) # DateLastModified
        genome_data_output.append(ph_genome.get_annotation_author()) # AnnotationAuthor


        #Genome data checks
        genome_data_output.append(ph_genome.get_nucleotide_errors())# sequence contains std nucleotides?
        genome_data_output.append(ph_genome.get_cds_features_with_translation_error_tally())# # translations with non-std amino acids
        genome_data_output.append(ph_genome.get_cds_features_boundary_error_tally())# # genes with non-standard start-stops
        genome_data_output.append(ph_genome.get_status_description_error())# status-description error
        genome_data_output.append(ph_genome.get_status_accession_error())# status-accession error


        #Phagesdb data
        if isinstance(pdb_genome,PhagesdbGenome):

            #General genome data
            genome_data_output.append(pdb_genome.get_phage_name())# Name
            genome_data_output.append(pdb_genome.get_search_name())# search name
            genome_data_output.append(pdb_genome.get_cluster())# cluster
            genome_data_output.append(pdb_genome.get_subcluster())# subcluster
            genome_data_output.append(pdb_genome.get_host())# Host
            genome_data_output.append(pdb_genome.get_accession())# Accession
            genome_data_output.append(pdb_genome.get_length())# sequence_length

            #Genome data checks
            genome_data_output.append(pdb_genome.get_nucleotide_errors())# sequence contains std nucleotides?
        else:
            genome_data_output.extend(['','','','','','','',''])



        #NCBI data
        if isinstance(ncbi_genome,NcbiGenome):

            #General genome data
            genome_data_output.append(ncbi_genome.get_phage_name())# the assigned phage name for this record
            genome_data_output.append(ncbi_genome.get_search_name())# the assigned search phage name for this record
            genome_data_output.append(ncbi_genome.get_record_id())# record_id
            genome_data_output.append(ncbi_genome.get_record_name())# record_name
            genome_data_output.append(ncbi_genome.get_record_accession())# record_accession
            genome_data_output.append(ncbi_genome.get_record_description())# record_description
            genome_data_output.append(ncbi_genome.get_record_source())# record_source
            genome_data_output.append(ncbi_genome.get_record_organism())# record_organism
            genome_data_output.append(ncbi_genome.get_source_feature_organism())# source_feature_organism
            genome_data_output.append(ncbi_genome.get_source_feature_host())# source_feature_host
            genome_data_output.append(ncbi_genome.get_source_feature_lab_host())# source_feature_lab_host
            genome_data_output.append(ncbi_genome.get_record_authors())# author list
            genome_data_output.append(ncbi_genome.get_length())# sequence_length
            genome_data_output.append(ncbi_genome.get_cds_features_tally())# # genes


            #Genome data checks
            genome_data_output.append(ncbi_genome.get_nucleotide_errors())# sequence contains std nucleotides?
            genome_data_output.append(ncbi_genome.get_cds_features_with_translation_error_tally())# # translations with non-std amino acids
            genome_data_output.append(ncbi_genome.get_cds_features_boundary_error_tally())# # genes with non-standard start-stops
            genome_data_output.append(ncbi_genome.get_product_descriptions_tally())# # genes with product descriptions
            genome_data_output.append(ncbi_genome.get_function_descriptions_tally())# # genes with function descriptions
            genome_data_output.append(ncbi_genome.get_note_descriptions_tally())# # genes with notes descriptions
            genome_data_output.append(ncbi_genome.get_missing_locus_tags_tally())# # genes with missing locus tags
            genome_data_output.append(ncbi_genome.get_locus_tag_typos_tally())# # genes with locus tag typos
            genome_data_output.append(ncbi_genome.get_description_field_error_tally())# # genes with descriptions in wrong field

        else:
            genome_data_output.extend(['','','','','','','','','','',\
                                        '','','','','','','','','','',\
                                        '','',''])

        #MySQL-phagesdb checks
        if isinstance(pdb_genome,PhagesdbGenome):
            genome_data_output.append(matched_genomes.get_phamerator_phagesdb_sequence_mismatch())# sequence
            genome_data_output.append(matched_genomes.get_phamerator_phagesdb_sequence_length_mismatch())# sequence length
            genome_data_output.append(matched_genomes.get_phamerator_phagesdb_cluster_subcluster_mismatch())# cluster_subcluster
            genome_data_output.append(matched_genomes.get_phamerator_phagesdb_accession_mismatch())# accession
            genome_data_output.append(matched_genomes.get_phamerator_phagesdb_host_mismatch())# host
        else:
            genome_data_output.extend(['','','','',''])


        #MySQL-NCBI checks
        if isinstance(ncbi_genome,NcbiGenome):
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_sequence_mismatch())# sequence
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_sequence_length_mismatch())# sequence length
            genome_data_output.append(matched_genomes.get_ncbi_record_header_fields_phage_name_mismatch())# PhageID or PhageName in record header fields mismatch
            genome_data_output.append(matched_genomes.get_ncbi_host_mismatch())# Host in record header or source feature mismatch
            genome_data_output.append(matched_genomes.get_ph_ncbi_author_error())# Author list is missing 'Hatfull'
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_perfect_matched_features_tally())# # genes perfectly matched
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_imperfect_matched_features_tally())# # genes imperfectly matched (different start sites)
            genome_data_output.append(matched_genomes.get_phamerator_features_unmatched_in_ncbi_tally())# # MySQL genes not matched
            genome_data_output.append(matched_genomes.get_ncbi_features_unmatched_in_phamerator_tally())# # NCBI genes not matched
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_different_descriptions_tally())# # genes with MySQL descriptions not in NCBI description fields
            genome_data_output.append(matched_genomes.get_phamerator_ncbi_different_translation_tally())# # genes perfectly matched with different translations
        else:
            genome_data_output.extend(['','','','','',\
                                        '','','','','',''])

        #Number of genes with errors
        genome_data_output.append(matched_genomes.get_total_number_genes_with_errors())# # genes with at least one error



        #Output phagesdb-NCBI checks
        if isinstance(pdb_genome,PhagesdbGenome) and isinstance(ncbi_genome,NcbiGenome):
            genome_data_output.append(matched_genomes.get_phagesdb_ncbi_sequence_mismatch())# sequence
            genome_data_output.append(matched_genomes.get_phagesdb_ncbi_sequence_length_mismatch())# sequence length
        else:
            genome_data_output.extend(['',''])

        genome_report_writer.writerow(genome_data_output)


        #Once all matched genome data has been outputted, iterate through all matched gene data


        perfectly_matched_features = matched_genomes.get_phamerator_ncbi_perfect_matched_features()
        imperfectly_matched_features = matched_genomes.get_phamerator_ncbi_imperfect_matched_features()
        ph_unmatched_features = matched_genomes.get_phamerator_features_unmatched_in_ncbi()
        ncbi_unmatched_features = matched_genomes.get_ncbi_features_unmatched_in_phamerator()

        all_features_list = []
        all_features_list.extend(perfectly_matched_features)
        all_features_list.extend(imperfectly_matched_features)
        all_features_list.extend(ph_unmatched_features)
        all_features_list.extend(ncbi_unmatched_features)



        #Iterate through the list of mixed feature objects
        for mixed_feature_object in all_features_list:

            feature_data_output = [] #Will hold all data for each gene

            #Gene summaries
            #Add MySQL genome search name to each gene row regardless of the type of CDS data (matched or unmatched)
            feature_data_output.append(ph_genome.get_search_name())# matched MySQL search name
            feature_data_output.append(mixed_feature_object.get_total_errors())# total # of errors for this gene

            #Now retrieve specific data
            if isinstance(mixed_feature_object,MatchedCdsFeatures):
                phamerator_feature = mixed_feature_object.get_phamerator_feature()
                ncbi_feature = mixed_feature_object.get_ncbi_feature()
            else:

                if isinstance(mixed_feature_object,PhameratorCdsFeature):
                    phamerator_feature = mixed_feature_object
                    ncbi_feature = ''
                elif isinstance(mixed_feature_object,NcbiCdsFeature):
                    phamerator_feature = ''
                    ncbi_feature = mixed_feature_object
                else:
                    phamerator_feature = ''
                    ncbi_feature = ''

            #MySQL feature
            if isinstance(phamerator_feature,PhameratorCdsFeature):

                #General gene data
                feature_data_output.append(phamerator_feature.get_phage_id())# phage_id
                feature_data_output.append(phamerator_feature.get_search_id())# search_id
                feature_data_output.append(phamerator_feature.get_type_id())# type_id
                feature_data_output.append(phamerator_feature.get_gene_id())# gene_id
                feature_data_output.append(phamerator_feature.get_gene_name())# gene_name
                feature_data_output.append(phamerator_feature.get_left_boundary())# left boundary
                feature_data_output.append(phamerator_feature.get_right_boundary())# right boundary
                feature_data_output.append(phamerator_feature.get_strand())# strand
                feature_data_output.append(phamerator_feature.get_translation())# translation
                feature_data_output.append(phamerator_feature.get_translation_length())# translation_length
                feature_data_output.append(phamerator_feature.get_notes())# notes

                #Gene data checks
                feature_data_output.append(phamerator_feature.get_amino_acid_errors())# translation contains std amino acids
                feature_data_output.append(phamerator_feature.get_boundary_error())# contains std start and stop coordinates

            else:

                feature_data_output.extend(['','','','','',\
                                            '','','','','',\
                                            '','',''])


            #NCBI feature
            if isinstance(ncbi_feature,NcbiCdsFeature):

                #General gene data
                feature_data_output.append(ncbi_feature.get_locus_tag())# locus_tag
                feature_data_output.append(ncbi_feature.get_gene_number())# gene_number
                feature_data_output.append(ncbi_feature.get_type_id())# type_id
                feature_data_output.append(ncbi_feature.get_left_boundary())# left boundary
                feature_data_output.append(ncbi_feature.get_right_boundary())# right boundary
                feature_data_output.append(ncbi_feature.get_strand())# strand
                feature_data_output.append(ncbi_feature.get_translation())# translation
                feature_data_output.append(ncbi_feature.get_translation_length())# translation_length
                feature_data_output.append(ncbi_feature.get_product_description())# product description
                feature_data_output.append(ncbi_feature.get_function_description())# function description
                feature_data_output.append(ncbi_feature.get_note_description())# note description

                #Gene data checks
                feature_data_output.append(ncbi_feature.get_amino_acid_errors())# translation contains std amino acids
                feature_data_output.append(ncbi_feature.get_boundary_error())# contains std start and stop coordinates
                feature_data_output.append(ncbi_feature.get_locus_tag_missing())# missing locus tag
                feature_data_output.append(ncbi_feature.get_locus_tag_typo())# locus tag typo
                feature_data_output.append(ncbi_feature.get_description_field_error())# description in function or note but not product


            else:
                feature_data_output.extend(['','','','','',\
                                            '','','','','',\
                                            '','','','','',''])


            #MySQL-NCBI checks
            if isinstance(mixed_feature_object,MatchedCdsFeatures):

                #If this is a matched cds feature, both MySQL and ncbi features should have identical unmatched_error value.
                feature_data_output.append(mixed_feature_object.get_phamerator_feature().get_unmatched_error())

                feature_data_output.append(mixed_feature_object.get_phamerator_ncbi_different_descriptions())# MySQL description in product, function, or note description
                feature_data_output.append(mixed_feature_object.get_phamerator_ncbi_different_start_sites())# same start site
                feature_data_output.append(mixed_feature_object.get_phamerator_ncbi_different_translations())# same translation
            else:
                feature_data_output.append(mixed_feature_object.get_unmatched_error())
                feature_data_output.extend(['','',''])


            gene_report_writer.writerow(feature_data_output)










    #Close all file handles
    close_all_files(file_handle_list)



    #Close script.
    end_time = time.strftime("%x %X")

    print('Start time: %s' %start_time)
    print('Stop time: %s' %end_time)

    print("\n\n\nDatabase comparison script completed.")


if __name__ == "__main__":
    main(sys.argv.insert(0, "empty"))
