#!/usr/bin/env python
#Database comparison script
#University of Pittsburgh
#Travis Mavrich
#20170203
#The purpose of this script is to compare the Phamerator, phagesdb, and NCBI databases for inconsistencies and report what needs to be updated.

# Note this script compares and matches data from Genbank data and Phamerator data. As a result, there are many similarly
# named variables. Variables are prefixed to indicate database:
#NCBI =  "ncbi".
#Phamerator = "ph"
#phagesdb = "pdb"


#Built-in libraries
import time, sys, os, getpass, csv, re, shutil
import json, urllib, urllib2

#Third-party libraries
import MySQLdb as mdb





#Define data classes
def remove_draft_suffix(value):
    # Is the word "_Draft" appended to the end of the name?
    value_truncated = value.lower()
    if value[-6:] == "_draft":
        value_truncated = value[:-6]
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

#Genome:

class PhameratorGenome:


    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields
        self.__phage_id = 'empty'
        self.__status = 'empty' #Final, Draft, Gbk version of genome data
        self.__accession = 'empty'
        self.__cluster_subcluster = 'empty' #Combined cluster_subcluster data.
        self.__retrieve_record_flag = 'empty'
        self.__phage_name = 'empty'
        self.__host = 'empty'
        self.__sequence = 'empty'
        self.__length = 'empty'
        self.__ncbi_update_flag = 'empty'
        self.__gene_data_list = []


        # Computed datafields
        self.__search_id = 'empty' # The phage ID void of "_Draft" and converted to lowercase
        self.__search_name = 'empty' # The phage name void of "_Draft" and converted to lowercase



    # Define all attribute setters:
    def set_phage_id(self,value):
        self.__phage_id = value
    def set_phage_name(self,value):
        self.__phage_name = value
    def set_host(self,value):
        self.__host = value
    def set_sequence(self,value):
        self.__sequence = value
    def set_length(self,value):
        self.__length = value
    def set_status(self,value):
        self.__status = value
    def set_accession(self,value):
        if value is None or value.strip() == '':
            self.__accession = 'none'
        else:
            self.__accession = value.split('.')[0]
    def set_cluster_subcluster(self,value):
        if value is None:
            self.__cluster_subcluster = 'Singleton'
        else:
            self.__cluster_subcluster = value
    def set_ncbi_update_flag(self,value):
        self.__ncbi_update_flag = value

    def set_search_id(self):
        self.__search_id = remove_draft_suffix(self.__phage_id)
    def set_search_name(self):
        self.__search_name = remove_draft_suffix(self.__phage_name)
    def set_gene_data_list(self,value):
        self.__gene_data_list = value


    # Define all attribute getters:
    def get_phage_id(self):
        return self.__phage_id
    def get_phage_name(self):
        return self.__phage_name
    def get_host(self):
        return self.__host
    def get_sequence(self):
        return self.__sequence
    def get_length(self):
        return self.__length
    def get_cluster(self):
        return self.__cluster
    def get_status(self):
        return self.__status
    def get_accession(self):
        return self.__accession
    def get_search_id(self):
        return self.__search_id
    def get_search_name(self):
        return self.__search_name
    def get_gene_data_list(self):
        return self.__gene_data_list


#Inherits Genome class
def PhagesdbGenome

    # Initialize all attributes:
    def __init__(self):

        # Non-computed datafields
        self.__phage_name = ''
        self.__host = ''
        self.__sequence = ''


        # Computed datafields
        self.__search_name = '' # The phage name void of "_Draft" and converted to lowercase
        self.__length = 0


    # Define all attribute setters:
    def set_phage_name(self,value):
        self.__phage_name = value
    def set_host(self,value):
        self.__host = value
    def set_sequence(self,value):
        self.__sequence = value
    def set_search_name(self):
        self.__search_name = remove_draft_suffix(self.__phage_name)
    def set_length(self):
        self.__length = len(self.__sequence)

    # Define all attribute getters:
    def get_phage_id(self):
        return self.__phage_id
    def get_name(self):
        return self.__name
    def get_host(self):
        return self.__host
    def get_sequence(self):
        return self.__sequence
    def get_length(self):
        return self.__length
    def get_cluster(self):
        return self.__cluster
    def get_status(self):
        return self.__status
    def get_accession(self):
        return self.__accession
    def get_search_name(self):
        return self.__search_name
    def get_length(self):
        return self.__length

        ###Add setter and getter functions. Pull from Genome class if applicable


#Inherits Genome class
def NcbiGenome

###Initialize and inherit

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
        self.__gene_data_list = []

        #Computed data fields
        self.__tally_function_descriptions = 0 #specific to NCBI records
        self.__tally_product_descriptions = 0 #specific to NCBI records
        self.__tally_note_descriptions = 0 #specific to NCBI records


        #Define setter functions
        def set_record_name(self,value):
            self.__record_name = value
        def set_record_id(self,value):
            self.__record_id = ''
        def set_record_accession(self,value):
            ###This part still needs worked on. Trim first accession in the list?
            ###self.__record_accession = ''
        def set_record_description(self,value):
            self.__record_description = ''
        def set_record_source(self,value):
            self.__record_source = ''
        def set_record_organism(self,value):
            self.__record_organism = ''
        def set_source_feature_organism(self,value):
            self.__source_feature_organism = ''
        def set_source_feature_host(self,value):
            self.__source_feature_host = ''
        def set_source_feature_lab_host(self,value):
            self.__source_feature_lab_host = ''
        def set_gene_data_list(self,value):
            ###Still working on this
            ###self.__gene_data_list = []





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
        def get_gene_data_list(self):






        ###For setting source field data, iterate through each gene



class PhameratorFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        #self.__phage_id = 'empty' #Phamerator PhageID, so that phage genome data can be accessed from PhageGenome objects
        self.__gene_id = 'empty' #Gene ID comprised of PhageID and Gene name
        self.__gene_name = 'empty'
        self.__type_id = 'empty' #Feature type: CDS, GenomeBoundary,or tRNA
        self.__left_boundary = 'empty' #Position of left boundary of feature, 0-indexed
        self.__right_boundary = 'empty' #Position of right boundary of feature, 0-indexed
        self.__strand = 'empty' #'forward', 'reverse', or 'NA'
        self.__translation = 'empty'
        self.__length = 'empty'
        self.__notes = 'empty'

        # Computed datafields
        self.__search_id = 'empty' # The phage ID void of "_Draft" and converted to lowercase
        self.__search_name = 'empty' # The phage name void of "_Draft" and converted to lowercase



    # Define all attribute setters:
    def set_phage_id(self,value):
        self.__phage_id = value
    def set_gene_id(self,value):
        self.__gene_id = value
    def set_gene_name(self,name):
        self.__gene_name = name
    def set_left_boundary(self,value):
        self.__left_boundary = value
    def set_right_boundary(self,value):
        self.__right_boundary = value
    def set_strand(self,value):
        self.__strand = parse_strand(value)
    def set_length(self,value):
        self.__length = value
    def set_translation(self,value):
        self.__translation = value
    def set_type_id(self,value):
        self.__type_id = value
    def set_notes(self,value):
        self.__notes = value
    def set_search_id(self):
        self.__search_id = remove_draft_suffix(self.__phage_id)
    def set_search_name(self):
        self.__search_name = remove_draft_suffix(self.__phage_name)

    def set_genome_boundary_straddle(self):
        # Determine if the gene straddles the genome boundary and set the attribute appropriately.
        if self.__left_boundary > self.__right_boundary:
            self.__genome_boundary_straddle = "yes"
        else:
            self.__genome_boundary_straddle = "no"


    # Define all attribute getters:
    def get_phage_id(self):
        return self.__phage_id
    def get_gene_id(self):
        return self.__gene_id
    def get_left_boundary(self):
        return self.__left_boundary
    def get_right_boundary(self):
        return self.__right_boundary
    def get_length(self):
        return self.__length
    def get_gene_name(self):
        return self.__gene_name
    def get_type_id(self):
        return self.__type_id
    def get_strand(self):
        return self.__strand
    def get_notes(self):
        return self.__notes
    def get_search_id(self):
        return self.__search_id
    def get_search_name(self):
        return self.__search_name
    def get_genome_boundary_straddle(self):
        return self.__genome_boundary_straddle



class NcbiFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        #self.__phage_id = 'empty' #Phamerator PhageID, so that phage genome data can be accessed from PhageGenome objects
        self.__locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.__gene_number = ''
        self.__type_id = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.__left_boundary = '' #Position of left boundary of feature, 0-indexed
        self.__right_boundary = '' #Position of right boundary of feature, 0-indexed
        self.__strand = '' #'forward', 'reverse', or 'NA'
        self.__translation = ''
        self.__product_description = ''
        self.__function_description = ''
        self.__note_description = ''



    # Define all attribute setters:
    def set_phage_id(self,value):
        self.__phage_id = value
    def set_gene_id(self,value):
        self.__gene_id = value
    def set_gene_name(self,name):
        self.__gene_name = name
    def set_left_boundary(self,value):
        self.__left_boundary = value
    def set_right_boundary(self,value):
        self.__right_boundary = value
    def set_strand(self,value):
        self.__strand = parse_strand(value)
    def set_length(self,value):
        self.__length = value
    def set_translation(self,value):
        self.__translation = value
    def set_type_id(self,value):
        self.__type_id = value
    def set_notes(self,value):
        self.__notes = value
    def set_search_name(self):
        self.__search_name = remove_draft_suffix(self.__phage_name)

    def set_genome_boundary_straddle(self):
        # Determine if the gene straddles the genome boundary and set the attribute appropriately.
        if self.__left_boundary > self.__right_boundary:
            self.__genome_boundary_straddle = "yes"
        else:
            self.__genome_boundary_straddle = "no"


    # Define all attribute getters:
    def get_phage_id(self):
        return self.__phage_id
    def get_gene_id(self):
        return self.__gene_id
    def get_left_boundary(self):
        return self.__left_boundary
    def get_right_boundary(self):
        return self.__right_boundary
    def get_length(self):
        return self.__length
    def get_gene_name(self):
        return self.__gene_name
    def get_type_id(self):
        return self.__type_id
    def get_strand(self):
        return self.__strand
    def get_notes(self):
        return self.__notes
    def get_genome_boundary_straddle(self):
        return self.__genome_boundary_straddle




####
class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        self.__search_name = 'empty'
        self.__phamerator_genome = 'empty'
        self.__phagesdb_genome = 'empty'
        self.__ncbi_genome = 'empty'


        #Phamerator-specific data
        self.__phamerator_phage_id = 'empty'
        self.__phamerator_name = 'empty'
        self.__phamerator_host = 'empty'
        self.__phamerator_accession = 'empty'

        #Phagesdb-specific data
        self.__phagesdb_name = 'empty'
        self.__phagesdb_host = 'empty'

        #NCBI-specific data
        self.__ncbi_record_name = 'empty'
        self.__ncbi_record_id = 'empty'
        self.__ncbi_record_description = 'empty'
        self.__ncbi_record_source = 'empty'
        self.__ncbi_record_organism = 'empty'
        self.__ncbi_source_feature_organism = 'empty'
        self.__ncbi_source_feature_host = 'empty'
        self.__ncbi_source_feature_lab_host = 'empty'

        #Add other header fields
        self.__ncbi_tally_product_descriptions = 0
        self.__ncbi_tally_function_descriptions = 0
        self.__ncbi_tally_note_descriptions = 0

        #Matched data comparison results
        self.__sequence_match_phamerator_phagesdb = True
        self.__sequence_match_phamerator_ncbi = True
        self.__sequence_match_phagesdb_ncbi = True
        self.__phamerator_ncbi_tally_mismatched_genes = 0
        self.__ncbi_tally_phage_name_typos = 0
        self.__ncbi_tally_host_typos = 0
        self.__ncbi_tally_geneID_typos = 0
        self.__tally_total_errors = 0


        ###Check for unexpected nucleotides






        # Define all attribute setters:
        def set_phamerator_genome(self,value):
            self.__phamerator_genome = value
        def set_phagesdb_genome(self,value):
            self.__phagesdb_genome = value
        def set_phamerator_genome(self,value):
            self.__ncbi_genome = value



    #def set_phage_id(self,value):
    #    self.__phage_id = value


        # Define all attribute getters:
        def get_phamerator_genome(self):
            return self.__phamerator_genome
        def get_phagesdb_genome(self):
            return self.__phagesdb_genome
        def get_phamerator_genome(self):
            return self.__ncbi_genome




    #def get_phage_id(self):
    #    return self.__phage_id
########










#Get the command line parameters
try:
    database = sys.argv[1] #What Phamerator database should be compared to phagesdb?
    updateFileDir = sys.argv[2] #What is the directory into which the report should go
except:
    print "\n\n\
            This is a python script to compare the Phamerator, phagesdb, and NCBI databases for inconsistencies.\n\
            It requires two arguments:\n\
            First argument: name of MySQL database that will be checked (e.g. 'Actino_Draft').\n\
            Second argument: directory path to where the consistency report should be made (csv-formatted).\n\"
    sys.exit(1)




#Expand home directory
home_dir = os.path.expanduser('~')


#Verify the folder for the consistency report exists

#Add '/' at the end if it's not there
if updateFileDir[-1] != "/":
    updateFileDir = updateFileDir + "/"


#Expand the path if it references the home directory
if updateFileDir[0] == "~":
    updateFileDir = home_dir + updateFileDir[1:]

#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
updateFileDir = os.path.abspath(updateFileDir)


if os.path.isdir(updateFileDir) == False:
    print "\n\nInvalid input for output folder.\n\n"
    sys.exit(1)







#Set up MySQL parameters
mysqlhost = 'localhost'
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"

#Get email infor for NCBI
contact_email = raw_input("Provide email for NCBI: ")

batch_size = ""
batch_size_valid = False
while batch_size_valid == False:
    batch_size = raw_input("Record retrieval batch size (must be greater than 0 and recommended is 100-200): ")
    print "\n\n"
    if batch_size.isdigit():
        batch_size = int(batch_size)
        if batch_size > 0:
            batch_size_valid = True
        else:
            print "Invalid choice."
            print "\n\n"
    else:
        print "Invalid choice."
        print "\n\n"


#Define several functions

#Print out statements to both the terminal and to the output file
def write_out(filename,statement):
    print statement
    filename.write(statement)


#For questionable data, user is requested to clarify if the data is correct or not
def question(message):
    number = -1
    while number < 0:
        value = raw_input("Is this correct? (yes or no): ")
        if (value.lower() == "yes" or value.lower() == "y"):
            number = 0
        elif (value.lower() == "no" or value.lower() == "n"):
            write_out(report_file,message)
            number = 1
        else:
            print "Invalid response."
    #This value will be added to the current error total. If 0, no error was encountered. If 1, an error was encountered.
    return number

#Exits MySQL
def mdb_exit(message):
    write_out(report_file,"\nError: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]`)
    write_out(report_file,message)
    write_out(report_file,"\nThe import script did not complete.")
    write_out(report_file,"\nExiting MySQL.")
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    cur.close()
    con.close()
    write_out(report_file,"\nExiting import script.")
    report_file.close()
    sys.exit(1)

#Closes all file handles currently open
def close_all_files(file_list):
    for file_handle in file_list:
        file_handle.close()


#You have to specify how many results to return at once. If you set it to 1 page long and 100,000 genomes/page, then this will return everything
pdb_sequenced_phages_url = "http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000"





#Create output directories
date = time.strftime("%Y%m%d")

output_folder = '%s_database_comparison' % date
output_path = os.path.join(updateFileDir,output_folder)


try:
    os.mkdir(os.path.join(updateFileDir,output_folder))
except:
    print "\nUnable to create output folder: %s" % os.path.join(updateFileDir,output_folder)
    sys.exit(1)

os.chdir(output_path)



#Create a folder names "genomes" to store NCBI records

genomes_folder = "genomes"
os.mkdir(genomes_folder)




#Open file to record update information
report_file = open(os.path.join(updateFileDir,output_folder,date + "_database_comparison.csv"), "w")
write_out(report_file,date + " Database comparison:\n\n\n")




#Retrieve database version
#Retrieve current genome data in database
#0 = PhageID
#1 = Name
#2 = HostStrain
#3 = Sequence
#4 = Length
#5 = status
#6 = Cluster
#7 = Accession
#8 = auto-update NCBI record flag
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
    con = mdb.connect(mysqlhost, username, password, database)
    con.autocommit(False)
    cur = con.cursor()
except:
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
    report_file.close()
    sys.exit(1)

try:
    cur.execute("START TRANSACTION")
    cur.execute("SELECT version FROM version")
    ph_version = str(cur.fetchone()[0])
    cur.execute("SELECT PhageID,Name,HostStrain,Sequence,SequenceLength,status,Cluster,Accession,RetrieveRecord FROM phage")
    ph_genome_data_tuples = cur.fetchall()
    cur.execute("SELECT PhageID,GeneID,Name,Start,Stop,Orientation,Translation,Notes from gene")
    ph_gene_data_tuples = cur.fetchall()
    cur.execute("COMMIT")
    cur.close()
    con.autocommit(True)

except:
    mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")

con.close()

#write_out(report_file,"\nPhamerator database: " + database)
#write_out(report_file,"\nPhamerator database version: " + db_version)








#Variable to track number of warnings and total_errors encountered
warnings = 0
total_errors = 0

#Create data sets
phageId_set = set()
phageName_set = set()
phageHost_set = set()
phageStatus_set = set()
phageCluster_set = set()
print "Preparing genome data sets from the phamerator database..."
ph_genome_object_dict = {} #Key = PhageID; #Value = genome_object
ph_search_name_set = set()
ph_search_name_duplicate_set = set()
ph_accession_set = set()
ph_accession_duplicate_set = set()
for genome_tuple in ph_genome_data_tuples:


    genome_object = PhameratorGenome()
    genome_object.set_phage_id(genome_tuple[0])
    genome_object.set_phage_name(genome_tuple[1])
    genome_object.set_host(genome_tuple[2])
    genome_object.set_sequence(genome_tuple[3])
    genome_object.set_length(genome_tuple[4])
    genome_object.set_status(genome_tuple[5])
    genome_object.set_cluster_subcluster(genome_tuple[6])
    genome_object.set_accession(genome_tuple[7])
    genome_object.set_ncbi_update_flag(genome_tuple[8])
    genome_object.set_search_id()
    genome_object.set_search_name()
    ph_genome_objects[genome_tuple[0]] = genome_object

    #This keeps track of whether there are duplicate phage names that will be used
    #to match up to phagesdb data.
    if genome_object.get_search_name() in ph_search_name_set:
        ph_search_name_duplicate_set.add(genome_object.get_search_name())
    else:
        ph_search_name_set.add(genome_object.get_search_name())

    #This keeps track of whether there are duplicate accession numbers that will be
    #used to match up to NCBI data
    if genome_object.get_accession() != 'none':
        if genome_object.get_accession() in ph_accession_set:
            ph_accession_duplicate_set.add(genome_object.get_accession())
        else:
            ph_accession_set.add(genome_object.get_accession())



#phagesdb relies on the phageName, and not the phageID. But Phamerator does not require phageName values to be unique.
#Check if there are any phageName duplications. If there are, they will not be able to be compared to phagesdb data.
if len(ph_search_name_duplicate_set) > 0:
    print "Warning: Data is not able to be matched to phagesdb because of the following non-unique phage Names in phamerator:"
    for element in ph_search_name_duplicate_set:
        print element
    raw_input("Press ENTER to proceed")
    ###Add a sys.exit line???

#Phamerator accessions aren't unique, so there could be duplicates
if len(ph_accession_duplicate_set) > 0:
    print "Warning: There are duplicate accessions in Phamerator. Unable to proceed with NCBI record retrieval."
    for accession in ph_accession_duplicate_set:
        print accession
        ph_accession_set.pop(accession)
    raw_input("Press ENTER to proceed")
    ###Add a sys.exit line???






ph_gene_objects_list = []
ph_gene_data_phage_id_set = set()
for gene_tuple in ph_gene_data_tuples:

    ph_gene_data_phage_id_set.add(gene_tuple[0])

    gene_object = PhameratorFeature()
    gene_object.set_phage_id(gene_tuple[0])
    gene_object.set_gene_id(gene_tuple[1])
    gene_object.set_gene_name(gene_tuple[2])
    gene_object.set_type_id('CDS')
    gene_object.set_left_boundary(gene_tuple[3])
    gene_object.set_right_boundary(gene_tuple[4])
    gene_object.set_strand(gene_tuple[5])
    gene_object.set_translation(gene_tuple[6])
    gene_object.set_notes(gene_tuple[7])
    gene_object.set_search_id()
    gene_object.set_search_name()
    ph_gene_objects_list.append(gene_object)


#ph_gene_objects_dict = {} #Key = PhageID; #Value = list of gene objects
for phage_id in ph_genome_object_dict.keys():
    genome_object = ph_genome_object_dict[phage_id]
    new_gene_object_list = []
    for gene_object in ph_gene_objects_list:
        if gene_object.get_phage_id == phage_id:
            new_gene_object_list.append(gene_object)

    genome_object.set_gene_data_list(new_gene_object_list)


#




















#Now retrieve all phagesdb data


#Retrieve a list of all sequenced phages listed on phagesdb
#You have to specify how many results to return at once. If you set it to 1 page long and 100,000 genomes/page, then this will return everything
print "Retrieving data from phagesdb..."
pdb_sequenced_phages_url = "http://phagesdb.org/api/sequenced_phages/?page=1&page_size=100000"
pdb_sequenced_phages_json = urllib.urlopen(pdb_sequenced_phages_url)
pdb_sequenced_phages_dict = json.loads(pdb_sequenced_phages_json.read())
pdb_sequenced_phages_json.close()

#Data for each phage is stored in a dictionary per phage, and all dictionaries are stored in a list under "results"
pdb_genome_dict = {}
pdb_search_name_set = set()
pdb_search_name_duplicate_set = set()
for element_dict in sequenced_phages_dict['results']:

    genome_object = PhagesdbGenome()
    genome_object.set_phage_name(element_dict['phage_name'])
    genome_object.set_host(element_dict['isolation_host']['genus'])
    genome_object.set_search_name('phage_name')


    #Check to see if there is a fasta file stored on phagesdb for this phage
    if element_dict['fasta_file'] is not None:
        fastafile_url = element_dict['fasta_file']

        response = urllib2.urlopen(fastafile_url)
        retrieved_fasta_file = response.read()
        response.close()

        #All sequence rows in the fasta file may not have equal widths, so some processing of the data is required
        #If you split by newline, the header is retained in the first list element
        split_fasta_data = retrieved_fasta_file.split('\n')
        pdb_sequence = ""
        index = 1
        while index < len(split_fasta_data):
            pdb_sequence = pdb_sequence + split_fasta_data[index].strip() #Strip off potential whitespace before appending, such as '\r'
            index += 1
        genome_object.set_sequence(pdb_sequence)
        genome_object.set_length()

    pdb_search_name = genome_object.get_search_name()
    if pdb_search_name in pdb_search_name_set:
        pdb_search_name_duplicate_set.add(pdb_search_name)
    else:
        pdb_search_name_set.add(pdb_search_name)
        pdb_genome_dict[pdb_search_name] = element_dict



#phagesdb phage names are unique, but just make sure after they are converted to a search name
if len(pdb_search_name_duplicate_set) > 0:
    print "Warning: phagesdb data is not able to be matched because of the following non-unique phage names in phagesdb:"
    for element in pdb_search_name_duplicate_set:
        print element
    raw_input("Press ENTER to proceed")
    ###Add a sys.exit line???


#Make sure all sequenced phage data has been retrieved
if (len(pdb_sequenced_phages_dict['results']) != pdb_sequenced_phages_dict['count'] or \
    len(pdb_sequenced_phages_dict['results']) != len(pdb_genome_dict)):

    print "\nUnable to retrieve all phage data from phagesdb due to default retrieval parameters."
    print "Update parameters in script to enable these functions."
    raw_input("Press ENTER to proceed")
















#Create list of Phamerator accession numbers and retrieve NCBI records
    #save records to folder if selected by user

#For each NCBI record:
    #parse genome features and store in objects



#For each phamerator genome:
    #retrieve matched phagesdb and NCBI record
    #make genome comparisons






####
print "\n\nRetrieving updated records from NCBI"


ph_accession_set

#Use esearch to verify the accessions are valid and efetch to retrieve the record
Entrez.email = contact_email
Entrez.tool = "NCBIRecordRetrievalScript"



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



#When retrieving in batch sizes, first create the list of values indicating which indices of the unique_accession_list should be used to create each batch
#For instace, if there are five accessions, batch size of two produces indices = 0,2,4
for batch_index_start in range(0,len(accession_retrieval_list),batch_size):


    if batch_index_start + batch_size > len(accession_retrieval_list):
        batch_index_stop = len(accession_retrieval_list)
    else:
        batch_index_stop = batch_index_start + batch_size

    current_batch_size = batch_index_stop - batch_index_start
    delimiter = " | "
    esearch_term = delimiter.join(unique_accession_list[batch_index_start:batch_index_stop])


    #Use esearch for each accession
    search_handle = Entrez.esearch(db = "nucleotide", term = esearch_term,usehistory="y")
    search_record = Entrez.read(search_handle)
    search_count = int(search_record["Count"])
    search_webenv = search_record["WebEnv"]
    search_query_key = search_record["QueryKey"]



    #Keep track of the accessions that failed to be located in NCBI
    if search_count < current_batch_size:
        search_accession_failure = search_record["ErrorList"]["PhraseNotFound"]

        #Each element in this list is formatted "accession[ACCN]"
        for element in search_accession_failure:
            retrieval_error_list.append(element[:-6])



    #Now retrieve all records using efetch
    fetch_handle = Entrez.efetch(db = "nucleotide", \
                                rettype = "gb", \
                                retmode = "text", \
                                retstart = 0, \
                                retmax = search_count, \
                                webenv = search_webenv, \
                                query_key = search_query_key)
    fetch_records = SeqIO.parse(fetch_handle,"genbank")

    for record in fetch_records:

        retrieved_record_list.append(record)

    search_handle.close()
    fetch_handle.close()



###Create retrieved record dictionary







#Report the genomes that could not be retrieved.
tally_retrieval_failure = len(retrieval_error_list)
for retrieval_error_accession in retrieval_error_list:

    genome_data = unique_accession_dict[retrieval_error_accession]
    phamerator_id = genome_data[0]
    phamerator_name = genome_data[1]
    phamerator_host = genome_data[2]
    phamerator_status = genome_data[3]
    phamerator_cluster = genome_data[4]
    phamerator_date = genome_data[5]
    phamerator_accession = genome_data[6]
    phamerator_retrieve = genome_data[7]
    phamerator_program = genome_data[8]




#Now that all records have been retrieved, check which records are newer than the upload date of the current version in phamerator.

ncbi_genome_dict = {} #Key = accession; #Value = genome data
for retrieved_record in retrieved_record_list:

    genome_object = NcbiGenome()

    try:
        genome_object.set_record_name(retrieved_record.name)
    except:
        genome_object.set_record_name('none')

    try:
        genome_object.set_record_id(retrieved_record.id)
    except:
        genome_object.set_record_id('none')

    try:
        #There may be a list of accessions associated with this file. I think the first accession in the list is the most recent.
        #Discard the version suffix if it is present in the Accession field (it might not be present)
        record_accession = retrieved_record.annotations["accessions"][0]
        record_accession = record_accession.split('.')[0]
        genome_object.set_record_accession(record_accession)
    except:
        genome_object.set_record_accession('none')

    try:
        genome_object.set_record_description(retrieved_record.description)
    except:
        genome_object.set_record_description('none')

    try:
        genome_object.set_record_source(retrieved_record.annotations['source'])
    except:
        genome_object.set_record_source('none')

    try:
        record_organism = retrieved_record.annotations['organism']
        if record_organism.split(' ')[-1] == 'Unclassified.'':
            genome_object.set_record_organism(record_organism.split(' ')[-2])
        else:
            genome_object.set_record_organism(record_organism.split(' ')[-1])
    except:
        genome_object.set_record_organism('none')

    ###Iterate through all features

    ###Set the following variables AFTER iterating through all features
    try:
        genome_object.set_source_feature_organism()
    except:
        genome_object.set_source_feature_organism('none')
    try:
        genome_object.set_source_feature_host()
    except:
        genome_object.set_source_feature_host('none')
    try:
        genome_object.set_source_feature_lab_host()
    except:
        genome_object.set_source_feature_lab_host('none')


    #After parsing all data, add to the ncbi dictionary
    ncbi_genome_dict[genome_object.get_accession()] = genome_object

    #If selected by user, save retrieved record to file
    ###Still need to add this functionality
    if save_records == 'yes':
        ncbi_filename = phamerator_name.lower() + "__" + retrieved_record_accession + ".gb"
        SeqIO.write(retrieved_record,os.path.join(ncbi_output_path,genomes_folder,ncbi_filename),"genbank")




#Now that all NCBI and phagesdb data is retrieved, match up to Phamerator genome data




ph_matched_to_pdb_count = 0
ph_unmatched_to_pdb_count = 0
ph_unmatched_to_pdb_genomes = [] #List of the phamerator genome objects with no phagesdb matches
ph_matched_to_ncbi_count = 0
ph_unmatched_to_ncbi_count = 0
ph_unmatched_to_ncbi_genomes = [] #List of the phamerator genome objects with no NCBI matches

#Iterate through each phage in Phamerator
matched_genomes_list = [] #Will store a list of MatchedGenomes objects

for phage_id in ph_genome_object_dict.keys():


    phamerator_genome = ph_genome_object_dict[phage_id]
    matched_objects = MatchedGenomes()
    matched_objects.set_phamerator_genome(phamerator_genome)

    #Match up phagesdb genome
    #First try to match up the phageID, and if that doesn't work, try to match up the phageName
    if phamerator_genome.get_search_id() in pdb_genome_dict.keys():
        pdb_genome = pdb_genome_dict[phamerator_genome.get_search_id()]
        ph_matched_to_pdb_count += 1

    elif phamerator_genome.get_search_name() in pdb_genome_dict.keys():
        pdb_genome = pdb_genome_dict[phamerator_genome.get_search_name()]
        ph_matched_to_pdb_count += 1

    else:
        ###Need to create output file to store the unmatched data
        ###write_out(report_file,"\nError: unable to find phageID %s or phageName %s from phagesdb." %(phameratorId,phameratorName))
        pdb_genome = 'none'
        ph_unmatched_to_pdb_count += 1
        ph_unmatched_to_pdb_genomes.append(phamerator_genome)

    matched_objects.set_phagesdb_genome(pdb_genome)

    #Now match up NCBI genome
    if phamerator_genome.get_accession() != 'none':
        try:
            ncbi_genome = ncbi_genome_dict[phamerator_genome.get_accession()]


    else:
        ncbi_genome = 'none'


#####












###OLD compare_databases



    #Compare genome sequence (if this option was selected)
    #phagesdb stores the 'official' nucleotide sequence, so make sure the sequence in Phamerator matches the sequence in phagesdb
    if check_sequence == True:




            #Compare phagesdb recorded genome size with genome size based on the stored fasta sequence
            if phagesdbSize != phagesdbSequence_size:
                write_out(report_file,"\nError: phagesdb genome size %s does not match fasta sequence genome size %s for phagesdb phageName %s." %(phagesdbSize,phagesdbSequence_size,phagesdbName))
                total_errors += 1

            #Compare phamerator recorded genome size with genome size based on the stored fasta sequence
            if phameratorSize != phagesdbSequence_size:
                write_out(report_file,"\nError: Phamerator genome size %s does not match fasta sequence genome size %s for phageID %s." %(phameratorSize,phagesdbSequence_size,phameratorId))
                total_errors += 1


            #Compare genome sequences stored in Phamerator and phagesdb fasta file
            if phameratorSequence.lower() != phagesdbSequence.lower():
                write_out(report_file,"\nError: Genome sequences stored in Phamerator and phagesdb do not match for phageID %s." %phameratorId)

                print phagesdbSequence[:10]
                print phagesdbSequence[-10:]
                print len(phagesdbSequence)
                print phameratorSequence[:10]
                print phameratorSequence[-10:]
                print len(phameratorSequence)


                total_errors += 1


        else:
            write_out(report_file,"\nError: no fasta file found on phagesdb for phageID %s." %phameratorID)
            total_errors += 1



write_out(report_file,"\nMatched phage tally: %s." %matched_count)
write_out(report_file,"\nUnmatched phage tally: %s." %unmatched_count)
write_out(report_file,"\nUnmatched phages:")
for element in unmatched_phageId_list:
    write_out(report_file,"\n%s" %element)



















#Close script.
write_out(report_file,"\n\n\n\Database comparison script completed.")

#close all file handles
