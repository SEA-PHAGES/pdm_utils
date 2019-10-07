"""Functions for converting a local SQL database query for a selction of phages to formatted files"""
"""Pipeline for converting a database, filtered for some phages, and writing appropriate output files"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.classes import genome, cds, mysqlconnectionhandler, filter
from pdm_utils.functions import flat_files, phamerator, basic
from functools import singledispatch
from typing import List, Dict
from pathlib import Path
import cmd, readline, os, sys, typing, argparse, csv

file_format_choices = ["gb", "fasta", "csv"]

def run_file_export(unparsed_args_list):
    """Uses parsed args to call file export functions 
    and run the entirety of the file export pipeline
    """ 

    args = parse_file_export_args(unparsed_args_list)
       
    if args.import_table != None: 
        phage_name_filter_list = \
                parse_phage_list_input(Path(args.import_table))
    elif args.single_genomes != None:
        phage_name_filter_list = \
                parse_phage_list_input(args.single_genomes)
    elif args.all:
        phage_name_filter_list = []

    if not args.interactive:
        if args.verbose:
            print("Establishing database connection to {}".\
                    format(args.database))
        sql_handle = establish_database_connection(args.database)
    
        seqfeature_file_output\
                (retrieve_seqrecord_from_database\
                (sql_handle, phage_name_filter_list,\
                verbose = args.verbose),\
                file_format = args.file_format,\
                export_path = Path(args.export_directory),\
                export_dir_name = args.folder_name,\
                verbose = args.verbose)

    else:
        interactive_file_export = Interactive_File_Export\
                (file_format = args.file_format, database = args.database,\
                phage_filter_list = phage_name_filter_list,\
                export_directory_name = args.folder_name,\
                export_directory_path = Path(args.export_directory))

        interactive_file_export.cmdloop()
        export_dict = interactive_file_export.data

        seqfeature_file_output\
                (retrieve_seqrecord_from_database\
                (export_dict["sql_handle"], export_dict["phage_filter_list"],\
                verbose = True),\
                file_format = export_dict["file_format"],\
                export_path = export_dict["directory_path"],\
                export_dir_name = export_dict["directory_name"],
                verbose = True)
          
def parse_file_export_args(unparsed_args_list):
    """Verifies the correct arguments are selected 
    for database to file
    :param unparsed_args_list:
        Input a series of unparsed args
    :type unparsed_args_list: list[str]: 
    :return parsed_args:
        Returns an Argument Parser object
        containing attributes with the
        parsed arguments
    :type parsed_args: ArgumentParser:
    """

    DATABASE_TO_FILE_HELP =(" Pipeline to export SQL database\
            data to formatted text files")
    DATABASE_HELP = ("Name of the MySQL database to import the genomes.")
    FILE_FORMAT = ("""
        Type of file to be exported into a directory
        The following are file export format options:
            -gb is the standard GenBank flat file format
            -fasta is a generic file containing a sequence                  and an information header
            -csv is a standard table 
                (comma seperated values) format 
                for information
        """) 
    
    IMPORT_TABLE_HELP = """
        Path to the CSV-formatted table containing
        instructions to process each genome.
        Structure of genome data table:
            1. --IGNORED-- Unused column 
            2. PhageID
            3. Host genus
            4. Cluster
            5. Subcluster
            6. Annotation status
            7. Gene description field (product, note, function)
            9. Accession
        """
    SINGLE_GENOMES_HELP = "Input the name of a single genome or multiple genomes to be exported"
    ALL_HELP = "Automatically selects all genomes from a database to be exported"

    VERBOSE_INTERACTIVITY = "Runs file export with minimal interactivity:\
            Complete print statements but no filtering options"
    FULL_INTERACTIVITY = "Runs file export with complete interactiviy:\
            Complete print statements and filtering options"
    SILENT_INTERACTIVITY = "Runs file export with no interactivity:\
            No print statements or filtering options"
    
    EXPORT_DIRECTORY_HELP =\
            "Input the path of the directory to store export files"
    FOLDER_NAME_HELP =\
            "Input the name of the folder which will contain the export files"  

    
    parser = argparse.ArgumentParser(description = DATABASE_TO_FILE_HELP)
    parser.add_argument("database", type=str, help = DATABASE_HELP)
    parser.add_argument("file_format", type=str, help = FILE_FORMAT,\
            default = "gb",\
            choices = file_format_choices)


    phage_list_input_args = parser.add_mutually_exclusive_group(required = True)
    phage_list_input_args.add_argument("-csv", "--import_table",\
            nargs = 1, type=str, help = IMPORT_TABLE_HELP)
    phage_list_input_args.add_argument("-sgs", "--single_genomes",\
            nargs = '+', type=str, help = SINGLE_GENOMES_HELP)
    phage_list_input_args.add_argument("-a", "--all", help = ALL_HELP, action = 'store_true')

    verbose_options = parser.add_mutually_exclusive_group()
    verbose_options.add_argument\
            ("-v", "--verbose", help = VERBOSE_INTERACTIVITY,\
            action = 'store_true')
    verbose_options.add_argument\
            ("-i", "--interactive", help = FULL_INTERACTIVITY,\
            action = 'store_true')
    verbose_options.add_argument\
            ("-s", "--silent", help = SILENT_INTERACTIVITY,\
            action = 'store_true')

    parser.add_argument("-dir", "--export_directory",\
            default = os.getcwd(), type=str,\
            help = EXPORT_DIRECTORY_HELP) 
    parser.add_argument("-name", "--folder_name",\
            default = "file_export", type = str, \
            help = FOLDER_NAME_HELP)

    

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    
    return(parsed_args) 
            
@singledispatch
def parse_phage_list_input(phage_list_input): 
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a list of phage names.
    :type phage_list_input: list[str]
    """

    print("Phage list input for database query is not a supported type")
    raise TypeError

@parse_phage_list_input.register(list)
def _(phage_list_input):
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a list of single genome name strings.
    :type phage_list_input: list[str]
    :return phage_list_input:
        Returns a list of single genome name strings.
    :type phage_list_input: list[str]
    """

    return phage_list_input

#Owen: Reundant, to be reimplemented 
###@parse_phage_list_input.register(str)
#def _(phage_list_input):
#    """Helper function to populate the filter list for a SQL query
#    :param phage_list_input:
#        Input a single genome name string.
#    :type phage_list_input: str:
#    :return phage_list:
#        Returns a list containing a singular phage name.
#    """
#    phage_list = []
#    phage_list.append(phage_list_input)
#    return phage_list

@parse_phage_list_input.register(Path)
def _(phage_list_input):
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a csv file path.
    :type phage_list_input: Path
    :return phage_list:
        Returns a list of phage names
    :type phage_list: list[str]
    """


    if not q.exists:
        raise ValueError("File {} is not found".\
                format(phage_list_input))

    phage_list = []
    with open(phage_list_input) as csv:
        csv_reader = csv.reader(csv, delimiter = ",")
        for name in csv_reader[1:]:
            phage_list.append(name[0])

    return phage_list

def establish_database_connection(database_name: str):
    """Creates a mysqlconnectionhandler object 
    and populates its credentials

    :param tag database_name:
        Input SQL database name.
    "type database_name: str
    """

    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = database_name
    sql_handle.get_credentials()
    try:
        sql_handle.validate_credentials
    except:
        print("SQL connection to database {}\
                with username and password failed".format(database_name))

    return sql_handle

def retrieve_seqrecord_from_database\
        (sql_database_handle: mysqlconnectionhandler.MySQLConnectionHandler\
        , phage_name_filter_list: List[str] = [],
        verbose = False):
    """Reads a local SQL database and converts it to a SeqRecord list

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :param phage_name_filter_list:
        Input a list of phage names within the selected
        SQL database.
    :type phage_name_filter_list: str[]
    """
    
    genome_query = "SELECT * FROM phage"
    cds_query = "SELECT * FROM gene"
    if verbose:
        print("Retrieving genomic data from database...")
    genome_list = phamerator.parse_genome_data\
            (sql_database_handle,\
                    phage_id_list = phage_name_filter_list\
                    ,phage_query = genome_query,\
                    gene_query = cds_query)
    if verbose:
        print("Retrieving database version...")
    database_versions = retrieve_database_version\
            (sql_database_handle)
    seq_record_list = []
    if verbose:
        print("Converting genomic data...")
    for genome in genome_list:

        if verbose:
            print("Converting {}".format(genome.name))
        set_cds_seqfeatures(genome)
        seqrecord = flat_files.genome_to_seqrecord(genome)
        append_database_version(seqrecord, database_versions)
        seq_record_list.append(seqrecord)


    return seq_record_list

def set_cds_seqfeatures(phage_genome: genome.Genome):
    """Helper function that queries for and returns cds data from a SQL database for a specific phage

    :param phage_genome:
        Input a genome object to query cds data for.
    :type phage_genome: genome
    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    """

    try:
        def _sorting_key(cds): return cds.left
        phage_genome.cds_features.sort(key=_sorting_key)
    except:
        if phage_genome == None:
            raise TypeError
        print("Genome cds features unable to be sorted")
        pass
    for cds in phage_genome.cds_features:
        cds.set_seqfeature()


def retrieve_database_version\
        (sql_handle: mysqlconnectionhandler.MySQLConnectionHandler):
    """Helper function that queries a SQL database for the database version and schema version

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :returns:
        database_versions_list(dictionary) is a dictionary
        of size 2 that contains values tied to keys
        "version" and "schema_version"
    """

    database_versions_list = phamerator.retrieve_data\
            (sql_handle, query='SELECT * FROM version')
    return database_versions_list[0]

def append_database_version(genome_seqrecord: SeqRecord,\
        version_data: Dict):
    """Helper function that appends the working database version in a comment within a SeqFeature annotation

    :param genome_seqfeature:
        Input a SeqRecord object generated from the working
        SQL database.
    :type genome_seqfeature: SeqRecord
    :param version_data:
        Input a version data dictionary parsed from a SQL database.
    :type version_data: dictionary
    """

    if len(version_data) < 2:
        print("Version data dictionary\
        containing SQL database version\
        data does not contain enough values")
        raise ValueError
    try:
        genome_seqrecord.annotations["comment"] =\
                genome_seqrecord.annotations["comment"] +\
                ("Database Version: {}; Schema Version: {}".format\
                (version_data["version"],\
                version_data["schema_version"]),)
    except:
        if genome_seqrecord == None:
            raise TypeError
        raise

def seqfeature_file_output(seqrecord_list: List[SeqRecord], file_format: str,\
        export_path: Path, export_dir_name: str = "file_export",\
        verbose = False):
    """Outputs files with a particuar format from a SeqRecord list

    :param seq_record_list:
        Input a list of SeqRecords.
    :type seq_record_list: SeqRecord[]
    :param file_format:
        Input SeqIO file output format.
    :type file_format: str
    :param export_path:
        Input the path for the placement of the directory
        of exported files.
    :type input_path: Path
    :param verbose:
        Input a boolean to represent the verbocity of
        the file export script.
    :type verbose: Boolean
    """
    if verbose:
        print("Resolving export path...")
    export_path = export_path.resolve()
    if not export_path.exists():
        print("Path parameter passed to seqfeature_file_output\
            is not a valid path")
        raise ValueError

    try: 
        export_path = Path(os.path.join(export_path, export_dir_name))
        if verbose:
            print("Resolving current export directory status...")
        if not export_path.is_dir():
            export_path.mkdir()
    except:
        print("Mkdir function failed to \
                create database_export_output\
                directory in {}".format(export_path))
        raise ValueError 

    if verbose:
        print("Writing selected data to files...")
    for record in seqrecord_list:
        print("Writing {}".format(record.name))
        output_dir="{}.{}".format\
                (record.name, file_format)
        output_path=export_path.joinpath(output_dir)
        output_handle=open(output_path, "w+")
        if file_format == "csv":
            pass 
        else:
            SeqIO.write(record, output_handle, file_format)
        output_handle.close()

def main(args):
    """Function to initialize file export"""

    if len(args) == 2 and (args[1] == "--interactive" or args[1] == "-i"):
        interactive_file_export = Interactive_File_Export()
        interactive_file_export.cmdloop()

    else: 
        args.insert(0, "blank_argument")
        run_file_export(args)

class Interactive_File_Export(cmd.Cmd):

    def __init__(self, file_format = "gb", database = None,\
        phage_filter_list = [], sql_handle = None, \
        export_directory_name = "file_export",\
        export_directory_path = Path(os.getcwd())):

        super(Interactive_File_Export, self).__init__()

        self.file_format = file_format
        self.database = database
        self.phage_filter_list = phage_filter_list
        self.sql_handle = sql_handle
        self.directory_name = export_directory_name
        self.directory_path = export_directory_path

        self.intro =\
        """---------------Hatfull Helper's File Export---------------
        Type help or ? to list commands.\n"""
        self.prompt = "(database) (export)user@localhost: "
        self.data = None

    def preloop(self):
        
        if self.database == None:
            print("---------------------Database Login ---------------------")
            self.database = input("MySQL database: ")

        if self.sql_handle == None or \
                self.sql_handle.database != self.database:
             
            self.sql_handle = establish_database_connection(self.database)
        self.prompt = "({}) (export){}@localhost: ".\
                format(self.database, self.sql_handle._username)

    def do_filter(self, *args):
        """Filters and queries database for genomes.
        USAGE: filter
        """
        db_filter = filter.Filter(self.database)
        interactive_filter = filter.Interactive_Filter\
                (db_filter = db_filter, sql_handle = self.sql_handle)
        interactive_filter.cmdloop()
        self.phage_filter_list = interactive_filter.data
         
    def do_format(self, *args):
        """Sets the current file format for genome export
        USAGE: format
        """

        format = input("File Format: ")
        if format in file_export_choices:
            self.file_format = format
            print("\
                    Changed format to {}.\n".format(self.file_format))
        else:
            print("File format not supported.\n")
        
    def do_directory_path(self, *args):
        """Sets the export directory name for genome export
        USAGE: format
        """
        
        path = Path(input("Export Directory Path: "))
        if path.resolve():
            self.directory_path = path
        else:
            print("\
                    Path not found.")

    def do_directory_name(self, *args):
        """Sets the export directory name for genome export
        USAGE: format
        """

        self.directory_name = input("Export Directory Name: ") 
    def do_clear(self, *args):        
        """Clears display terminal
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')
       
    def do_export(self, *args):
        """Exit interface and finish exporting files
        USAGE: export
        """
        print("\
                Initiating Export...\n")

        return True

    def do_exit(self, *args):
        """Exits program entirely without returning values
        USAGE: exit
        """
        print("\
                Exiting...\n")

        sys.exit(1)

    def postloop(self):

        self.data = {"sql_handle" : self.sql_handle,\
                "phage_filter_list" : self.phage_filter_list,\
                "file_format" : self.file_format,\
                "directory_path" : self.directory_path,\
                "directory_name" : self.directory_name}

if __name__ == "__main__":
    main(sys.argv)
