"""Functions for converting a local SQL database query for a selction of phages
to formatted files"""
"""Pipeline for converting a database, filtered for some phages, and writing
appropriate output files"""

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


#Global file constants
file_format_choices = ["gb", "fasta", "clustal", "embl",
                       "fasta-2line", "fastq", "fastq-solexa",
                       "fastq-illumina","ig", "igmt", "nexus",
                       "phd", "phylip", "pir", "seqxml","sff",
                       "stockholm", "tab", "qual"]

def run_file_export(unparsed_args_list):
    """Uses parsed args to run the entirety of the file export pipeline
    """

    args = parse_file_export_args(unparsed_args_list)

    if args.import_table:
        phage_filter_list = \
                    parse_phage_list_input(args.import_table[0])
    elif args.single_genomes:
        phage_filter_list = \
                    parse_phage_list_input(args.single_genomes)
    else:
        phage_filter_list = []

    if not args.interactive:
        if args.database == None:
            args.database = input("MySQL database: ")

        if args.verbose:
            print("Establishing connection to {}...".\
                         format(args.database))

        sql_handle = establish_database_connection(args.database)

        execute_file_export(args.file_format, sql_handle, phage_filter_list,
                            args.export_directory, args.folder_name,
                            verbose=args.verbose, csv_log=args.csv_log)
    else:
        interactive_filter = Cmd_Export(file_format=args.file_format,
                                        database=args.database,
                                        phage_filter_list=phage_filter_list,
                                        sql_handle=None,
                                        export_directory_name=args.folder_name,
                                        export_directory_path=\
                                                        args.export_directory)

        interactive_filter.cmdloop()

def execute_file_export(file_format, sql_handle, phage_filter_list,
                        export_path, folder_name,
                        verbose=False, csv_log=False):
    """Executes the entirety of the file export pipeline by calling its
       various functions

       :param file_format:
            Input a recognized SeqIO file format.
       :type file_format: str:
       :param database:
            Input the name of the local phameratory database.
       :type database: str:
       :param phage_filter_list:
            Input a list of names of phages.
       :type phage_filter_list: List[str]
       :param export_directory_path:
            Input the path for the created directory.
       :type export_path: Path:
       :param folder_name:
            Input the name for the created directory.
       :type export_directory_name: str:
       :param verbose:
            Input a toggle for optional printed statements.
       :type verbose: Boolean:
       :param csv_log:
            Input a toggle for an optional csv log.
       :type csv_log: Boolean:
       """

    if verbose:
        print("Retrieving genomic data from {}...".\
                format(sql_handle.database))
    genomes = phamerator.parse_genome_data(
                      sql_handle,
                      phage_id_list=phage_filter_list,
                      phage_query="SELECT * FROM phage",
                      gene_query="SELECT * FROM gene")

    if verbose:
        print("Converting genomic data to SeqRecord format...")
    seqrecords = []
    for gnm in genomes:
        set_cds_seqfeatures(gnm)
        if verbose:
            print(f"Converting {gnm.name}")
        seqrecords.append(flat_files.genome_to_seqrecord(gnm))

    if verbose:
        print("Retrieving database version...")
    db_version = retrieve_database_version(sql_handle)

    if verbose:
        print("Appending database version...")
    for record in seqrecords:
        append_database_version(record, db_version)

    seqrecord_to_file(seqrecords,
                      file_format,
                      export_path,
                      export_dir_name=folder_name,
                      verbose=verbose)

    if csv_log:
        if verbose:
            print("Writing csv log...")
        write_csv_log(genomes, export_path, export_dir_name=folder_name)

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
        The following are standard file export format options:
            -gb is the standard GenBank flat file format
            -fasta is a generic file containing a sequence
                and an information header
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
    SINGLE_GENOMES_HELP = "Input the name of a single genome or"
    "multiple genomes to be exported"
    ALL_HELP = "Automatically selects all genomes"
    "from a database to be exported"

    VERBOSE_HELP = "Runs file export with minimal interactivity:" \
            "Complete print statements but no filtering options"
    INTERACTIVE_HELP = "Runs file export with complete interactiviy:" \
           " Loads interactive menu to perform file_export."
    EXPORT_DIRECTORY_HELP =\
            "Input the path of the directory to store export files"
    FOLDER_NAME_HELP =\
            "Input the name of the folder which will contain the export files"

    CSV_LOG_HELP =\
            "Toggle to export a log file with data from the exported genomes"


    parser = argparse.ArgumentParser(description = DATABASE_TO_FILE_HELP)
    parser.add_argument("-db", "--database", type=str,
                        help=DATABASE_HELP, default=None)
    parser.add_argument("-ff", "--file_format", type=str, help=FILE_FORMAT,
                            default = "gb", choices=file_format_choices)


    phage_list_input_args = parser.add_mutually_exclusive_group()
    phage_list_input_args.add_argument("-tin", "--import_table",
                                       nargs=1, type=convert_file_path,
                                       help=IMPORT_TABLE_HELP)
    phage_list_input_args.add_argument("-sgin", "--single_genomes",
                                       nargs='+', type=str,
                                       help=SINGLE_GENOMES_HELP)
    phage_list_input_args.add_argument("-a", "--all",
                                       help = ALL_HELP, action ='store_true')

    parser.add_argument("-v", "--verbose",
                        default=False,
                        help=VERBOSE_HELP,
                        action='store_true')
    parser.add_argument("-i", "--interactive",
                        default=False,
                        help=INTERACTIVE_HELP,
                        action='store_true')
    parser.add_argument("-dir", "--export_directory",
                        default=Path.cwd(), type=convert_dir_path,
                        help=EXPORT_DIRECTORY_HELP)
    parser.add_argument("-name", "--folder_name",
                        default="file_export", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-log", "--csv_log",
                        help=CSV_LOG_HELP, action='store_true')
    parsed_args = parser.parse_args(unparsed_args_list[2:])

    return(parsed_args)

def convert_path(path: str):
    """
    Function to convert a string to a working Path object
    :param path:
        Input a string to be converted into a Path object.
    :type path: str
    :return path_object:
    Returns a path object from the inputted
    :type path_object: Path
    """
    path_object = Path(path)
    if "~" in path:
        path_object = path_object.expanduser()

    if path_object.exists():
        return path_object
    elif path_object.resolve().exists():
        path_object = path_object.resolve()

    print("String input failed to be converted to a working Path object " \
          "Path does not exist")

    raise ValueError

def convert_dir_path(path: str):
    """
    Helper function to convert a string to a working Path object
    :param path:
        Input a string to be converted into a Path object.
    :type path: str
    :return path_object:
        Returns a path object directing to a directory.
    :type path_object: Path
    """

    path_object = convert_path(path)

    if path_object.is_dir():
        return path_object
    else:
        print("Path input does not direct to a folder")
        raise ValueError

def convert_file_path(path: str):
    """
    Helper function to convert a string to a working Path object
    :param path:
        Input a string to be converted into a Path object.
    :type path: str
    :return path_object:
        Returns a path object directing to a file.
    :type path_object: Path
    """
    path_object = convert_path(path)

    if path_object.is_file():
        return path_object
    else:
        print("Path input does not direct to a file")
        raise ValueError

@singledispatch
def parse_phage_list_input(phage_list_input):
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a list of phage names.
    :type phage_list_input: list[str]
    """

    print("Phage list input for database query is not a supported type")
    raise TypeError

@parse_phage_list_input.register(Path)
def _(phage_list_input):
    phage_list = []
    with open(phage_list_input, newline = '') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ",", quotechar = '|')
        for name in csv_reader:
            phage_list.append(name[0])
    return phage_list

@parse_phage_list_input.register(list)
def _(phage_list_input):
    return phage_list_input

def establish_database_connection(database_name: str):
    """Creates a mysqlconnectionhandler object
    and populates its credentials

    :param tag database_name:
        Input SQL database name.
    "type database_name: str
    """

    if not isinstance(database_name, str):
        print("establish_database_connection requires string input")
        raise TypeError
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = database_name
    sql_handle.get_credentials()
    try:
        sql_handle.open_connection()
    except:
        print(f"SQL connection to database {database_name}"
            "with username and password failed")
        raise RuntimeError

    return sql_handle

def set_cds_seqfeatures(phage_genome: genome.Genome):
    """Helper function that queries for and returns
    cds data from a SQL database for a specific phage

    :param phage_genome:
        Input a genome object to query cds data for.
    :type phage_genome: genome
    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    """

    try:
        def _sorting_key(cds_feature): return cds_feature.left
        phage_genome.cds_features.sort(key=_sorting_key)
    except:
        if phage_genome == None:
            raise TypeError
        print("Genome cds features unable to be sorted")
        pass
    for cds_feature in phage_genome.cds_features:
        cds_feature.set_seqfeature()

def retrieve_database_version(sql_handle):
    """Helper function that queries a SQL database
    for the database version and schema version

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :returns:
        database_versions_list(dictionary) is a dictionary
        of size 2 that contains values tied to keys
        "Version" and "SchemaVersion"
    """

    database_versions_list = phamerator.retrieve_data(
            sql_handle, query='SELECT * FROM version')
    return database_versions_list[0]

def append_database_version(genome_seqrecord: SeqRecord, version_data: Dict):
    """Helper function that appends the working database version
    in a comment within a SeqFeature annotation

    :param genome_seqfeature:
        Input a SeqRecord object generated from the working
        SQL database.
    :type genome_seqfeature: SeqRecord
    :param version_data:
        Input a version data dictionary parsed from a SQL database.
    :type version_data: dictionary
    """

    if len(version_data) < 2:
        print("Version data dictionary "
        "containing SQL database version "
        "data does not contain enough values")
        raise ValueError
    try:
        genome_seqrecord.annotations["comment"] =\
                genome_seqrecord.annotations["comment"] + (
                    "Database Version: {}; Schema Version: {}".format(\
                            version_data["Version"], version_data["SchemaVersion"]),)
    except:
        if isinstance(genome_seqrecord, SeqRecord):
            raise ValueError

        elif genome_seqrecord == None:
            raise TypeError
        raise

def seqrecord_to_file(seqrecord_list: List[SeqRecord],
                           file_format: str,
                           export_path: Path,
                           export_dir_name="file_export",
                           verbose=False):
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
        export_path = export_path.joinpath(export_dir_name)
        if verbose:
            print("Resolving current export directory status...")
        if not export_path.is_dir():
            export_path.mkdir()
    except:
        print("Mkdir function failed to"
              f" create database_export_output directory in {export_path}")
        raise ValueError

    if verbose:
        print("Writing selected data to files...")
    for record in seqrecord_list:
        if verbose:
            print(f"Writing {record.name}")
        output_dir = f"{record.name}.{file_format}"
        output_path = export_path.joinpath(output_dir)
        output_handle = output_path.open(mode='w')
        SeqIO.write(record, output_handle, file_format)
        output_handle.close()

def write_csv_log(genome_list, export_path, export_dir_name="file_export"):
    """Writes a formatted csv file from genome objects"""


    export_path = export_path.joinpath(export_dir_name)

    if not export_path.exists():
        export_path.mkdir()

    log_path = Path(os.path.join(export_path, "log.csv"))
    logversion = 1

    while(log_path.exists()):
        logversion += 1
        log_path = Path(os.path.join(export_path, f"log{logversion}.csv"))

    csv_format = []
    csv_format.append(["PhageID",
                       "Accession",
                       "Name",
                       "HostStrain",
                       "SequenceLength",
                       "DateLastModified",
                       "Notes",
                       "GC",
                       "Cluster",
                       "Subcluster",
                       "Status",
                       "RetrieveRecord",
                       "AnnotationQC",])
    for gnm in genome_list:
        csv_format.append([gnm.id,
                           gnm.accession,
                           gnm.name,
                           gnm.host_genus,
                           gnm.length,
                           gnm.date,
                           gnm.description,
                           gnm.gc,
                           gnm.cluster,
                           gnm.subcluster,
                           gnm.annotation_status,
                           gnm.retrieve_record,
                           gnm.annotation_author])
    log_path.touch()
    with open(log_path, 'w', newline="") as csv_file:
        logwriter = csv.writer(csv_file, delimiter=",",
                               quotechar = "|",
                               quoting = csv.QUOTE_MINIMAL)
        for row in csv_format:
            logwriter.writerow(row)

def main(args):
    """Function to initialize file export"""

    args.insert(0, "blank_argument")
    run_file_export(args)

class Cmd_Export(cmd.Cmd):

    def __init__(self, file_format="gb",database=None,
                 phage_filter_list=[], sql_handle=None,
                 export_directory_name="file_export",
                 export_directory_path = Path.cwd()):

        super(Cmd_Export, self).__init__()

        self.file_format = file_format
        self.database = database
        self.phage_filter_list = phage_filter_list
        self.sql_handle = sql_handle
        self.directory_name = export_directory_name
        self.directory_path = export_directory_path
        self.csv_toggle = False

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

    def do_search(self, *args):
        """Filters and queries database for genomes.
        """
        db_filter = filter.Filter(self.database, self.sql_handle)
        interactive_filter = filter.Cmd_Filter(
                db_filter=db_filter, sql_handle=self.sql_handle)
        interactive_filter.cmdloop()
        self.phage_filter_list = interactive_filter.data

    def do_folder(self, *args):
        """Selects options for current folder
        FOLDER OPTIONS: Format, Path, Name, Export, Log
        """

        options = ["format", "path", "name", "export", "log"]
        option = args[0].lower()

        if option in options:
            if option == "format":
                self.folder_format()
            elif option == "path":
                self.folder_directory_path()
            elif option == "name":
                self.folder_directory_name()
            elif option =="export":
                self.folder_export()
            elif option == "log":
                if self.csv_toggle:
                    print("Csv logging off. \n")
                    self.csv_toggle = False
                else:
                    print("Csv logging on. \n")
                    self.csv_toggle = True
        else:
            print("""Folder command option not supported
            FOLDER OPTIONS: Format, Path, Name, Export, Log
            """)

    def folder_format(self):
        """Sets the current file format for genome export
        """

        format = input("File Format: ")
        if format in file_export_choices:
            self.file_format = format
            print("\
                    Changed format to {}.\n".format(self.file_format))
        else:
            print("File format not supported.\n")

    def folder_directory_path(self):
        """Sets the export directory name for genome export
        USAGE: format
        """

        path = Path(input("Export Directory Path: "))
        if path.resolve():
            self.directory_path = path
        else:
            print("\
                    Path not found.")

    def folder_directory_name(self):
        """Sets the export directory name for genome export
        USAGE: format
        """

        self.directory_name = input("Export Directory Name: ")

    def folder_export(self, *args):
        """Exit interface and finish exporting files
        USAGE: export
        """
        print("\
                Initiating Export...\n")
        execute_file_export(self.file_format, self.sql_handle,
                                self.phage_filter_list, self.directory_path,
                                self.directory_name,
                                verbose=False, csv_log=False)

    def do_clear(self, *args):
        """Clears display terminal
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')
        print(self.intro)

    def do_exit(self, *args):
        """Exits program entirely without returning values
        USAGE: exit
        """
        print("       Exiting...\n")

        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv)
