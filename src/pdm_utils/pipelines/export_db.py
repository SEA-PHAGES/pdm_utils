"""Pipeline for converting a database, filtered for some phages, and writing
appropriate output files"""

#for (to be moved) cds_to_seqrecord dependancy
from Bio.Alphabet import IUPAC

import argparse
import cmd
import csv
from functools import singledispatch
import os
from pathlib import Path
import re
import sys
import time
import copy
from typing import List, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from pdm_utils.classes import genome, cds, filter
from pdm_utils.functions import flat_files, mysqldb

# Valid file formats using Biopython
BIOPYTHON_CHOICES = ["gb", "fasta", "clustal", "fasta-2line", "nexus",
                       "phylip", "pir", "stockholm", "tab"]
# Valid Biopython formats that crash the script due to specific values in
# some genomes that can probably be fixed relatively easily and implemented.
# BIOPYTHON_CHOICES_FIXABLE = ["embl", "imgt", "seqxml"]
# Biopython formats that are not applicable.
# BIOPYTHON_CHOICES_INVALID = ["fastq", "fastq-solexa", "fastq-illumina",
#                              "phd", "sff", "qual"]
# Biopython formats that are not writable
# BIOPYTHON_CHOICES_NOT_WRITABLE = ["ig"]

def run_export(unparsed_args_list):
    """Uses parsed args to run the entirety of the file export pipeline.

    :param unparsed_args_list:
        Input a list of command line args.
    :type unparsed_args_list: List[str]
    """

    args = parse_export(unparsed_args_list)
    if args.verbose:
        print("Please input database credentials:")
    engine = mysqldb.connect_to_db(args.database)


    csvx = False
    ffx = None
    dbx = False
    ix = False

    values_list = None
    filters = None
    group = None

    if args.pipeline in BIOPYTHON_CHOICES+["csv"]:
        values_list = parse_value_list_input(args.input)
        filters = filter.parse_filters(args.filters)
        groups = filter.parse_groups(args.groups)

        if args.pipeline == "csv":
            csvx = True
        else:
            ffx = args.pipeline
    elif args.pipeline == "sql":
        dbx = True
        groups = []
    elif args.pipeline == "I":
        ix = True
    else:
        print("ERROR: Export pipeline option discrepency.\n"
              "Pipeline parsed from command line args is not supported")
        raise ValueError

    if not ix:
        execute_export(engine, args.output_path, args.output_name,
                            values_list=values_list, verbose=args.verbose,
                            csv_export=csvx, ffile_export=ffx, db_export=dbx,
                            table=args.table,
                            filters=filters, groups=groups)
    else:
        pass

def parse_export(unparsed_args_list):
    """Parses export_db arguments and stores them with an argparse object.

    :param unparsed_args_list:
        Input a list of command line args.
    :type unparsed_args_list: List[str]
    :returns tuple(export.pipeline, parsed_args):
        Return a tuple of the export pipeline and a parsed args object.
    :type tuple(export.pipeline, parsed_args): (str, Namespace)
    """
    EXPORT_SELECT_HELP = """
        Select a export pipeline option to export genomic data:
            -csv export (csv)
            -formatted file export (ffx)
            -database export (sql)
            -interactive export (I)
        """
    CSV_EXPORT_HELP = """
        Export option to export a csv file
        containing information about selected genomes.
            """
    DATABASE_HELP = "Name of the MySQL database to export from."
    DATABASE_EXPORT_HELP = """
        Export option to dump the current database
        into a .sql file and its version into a .version file.
            """
    FORMATTED_FILE_EXPORT_HELP = """
        Export option to export formatted files containing
        information about individual genomes.
            """
    INTERACTIVE_HELP = """
        Export option that enables full interactive walkthrough export.
        """
    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    SINGLE_GENOMES_HELP = """
        Selection input option that imports values from cmd line input.
            Follow selection argument with space separated
            names of genomes in the database.
        """
    TABLE_HELP = """
        Selection option that changes the table in the database selected from.
            Follow selection argument with a valid table from the database.
        """
    FILTERS_HELP = """
        Genome selection option that constructs primary attributes list from
        inputted filters.
            Follow selection argument with formatted filter request:
                Table.Field=Value
        """
    GROUPS_HELP = """
        Genome selection option that groups primary attributes list list from
        inputted groups.
            Follow selection argument with supported grouping option.
        """
    VERBOSE_HELP = """
        Export option that enables progress print statements.
        """
    OUTPUT_PATH_HELP = """
        Export option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired export directory.
        """
    OUTPUT_NAME_HELP = """
        Export option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
        """
    FILE_FORMAT_HELP = """
        Positional argument specifying the format of the file to export
        """
    export_options = BIOPYTHON_CHOICES + ["csv", "sql"]

    selection_parser = argparse.ArgumentParser()
    selection_parser.add_argument("pipeline", type=str,
                                  choices=export_options,
                                  nargs="?", default="I",
                                  help=EXPORT_SELECT_HELP)
    export = selection_parser.parse_args([unparsed_args_list[2]])

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str, help=DATABASE_HELP)

    parser.add_argument("-o", "--output_name", type=str, help=OUTPUT_NAME_HELP)
    parser.add_argument("-p", "--output_path", type=convert_dir_path,
                        help=OUTPUT_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)


    if export.pipeline in (BIOPYTHON_CHOICES + ["csv"]):
        table_choices = dict.fromkeys(BIOPYTHON_CHOICES, ["phage"])
        table_choices.update({"csv": ["domain", "gene", "gene_domain",
                                       "phage", "pham", "pham_color",
                                       "tmrna", "trna", "trna_structures"]})
        parser.add_argument("-t", "--table", help=TABLE_HELP,
                            choices=table_choices[export.pipeline])

        parser.add_argument("-if", "--import_file", type=convert_file_path,
                                help=IMPORT_FILE_HELP, dest="input",
                                default=[])
        parser.add_argument("-in", "--import_names", nargs="*",
                                help=SINGLE_GENOMES_HELP, dest="input",
                                default=[])
        parser.add_argument("-f", "--filter", nargs="*",
                                help=FILTERS_HELP,
                                dest="filters")

        if export.pipeline != "csv":
            parser.add_argument("-g", "--groups", type=str.lower,
                            help=GROUPS_HELP, nargs="*",
                            dest="groups")

    date = time.strftime("%Y%m%d")
    default_folder_name = f"{date}_export"
    default_folder_path = Path.cwd()

    parser.set_defaults(output_name=default_folder_name,
                        output_path=default_folder_path,
                        verbose=False, input=[],
                        pipeline=export.pipeline,
                        table="phage", filters=[], groups=[],
                        ix=False, csvx=False, dbx=False, ffx=False)

    parsed_args = parser.parse_args(unparsed_args_list[3:])
    return parsed_args

def execute_export(engine, output_path, output_name,
                        values_list=[], verbose=False,
                        csv_export=False, ffile_export=None, db_export=False,
                        table="phage", filters=[], groups=[]):
    """Executes the entirety of the file export pipeline.

    :param engine:
        Input a valid SQLAlchemy Engine object.
    :type engine: Engine:
    :param export_path:
        Input a valid path to place export folder.
    :type export_path: Path
    :param folder_name:
        Input a name for the export folder.
    :type folder_name: str
    :param phage_filter_list:
        Input a list of phageIDs.
    :type phage_filter_list: List[str]
    :param verbose:
        Input a boolean value for verbose option.
    :type verbose: boolean
    :param csv_export:
        Input a boolean value to toggle csv_export.
    :type csv_export: boolean
    :param ffile_export:
        Input a SeqIO supported file format to toggle ffile_export.
    :type ffile_export: str
    :param db_export:
        Input a boolean value to toggle db_export.
    :type db_export: boolean
    :param filters:
        Input a list of lists with filter values
    :type filters: List[List[str]]
    :param groups:
        Input a list of supported group values.
    :type groups: List[str]
    """

    if verbose:
        print("Retrieving database version...")
    db_version = mysqldb.get_version_table_data(engine)

    if verbose:
        print("Creating export folder...")
    export_path = output_path.joinpath(output_name)

    if(export_path.is_dir()):
        export_version = 1
        while(export_path.is_dir()):
            export_version += 1
            export_name = f"{output_name}_{export_version}"
            export_path = export_path.with_name(export_name)
    else:
        export_name = output_name

    if db_export:
        if verbose:
            print("Writing SQL database file...")
        write_database(engine, db_version["Version"],
                        output_path, output_name=export_name)

    if csv_export or ffile_export:
        db_filter = build_filter(engine, table, values_list,
                                 filters, verbose=verbose)
        if csv_export:
            file_name = f"{engine.url.database}_{table}"
            execute_csv_export(db_filter, engine,
                               output_path, export_name,
                               csv_name=file_name, table=table,
                               verbose=verbose)

        if ffile_export != None:
            if groups:
                folder_path = output_path.joinpath(export_name)
                folder_path.mkdir(exist_ok=True)
                ffx_grouping(engine, folder_path, groups, db_filter,
                             db_version, ffile_export,
                             table=table, verbose=verbose)
            else:
                execute_ffx_export(engine,
                                   db_filter.results(verbose=verbose),
                                   ffile_export,output_path, export_name,
                                   db_version, verbose=verbose,
                                   data_name=f"{engine.url.database}.{table}",
                                   table=table)

def execute_ffx_export(engine, values, file_format,
                       output_path, output_name, db_version,
                       verbose=False, data_name="database", table="phage"):
    """
    Executes the ffx export  pipeline by calling its
    various functions

        :param engine:
            Input a valid SQLAlchemy Engine object.
        :type engine: Engine
        :param db_filter:
            Input a db_filter with a loaded list of phageIDs.
        :type db_filter: Filter
        :param file_format:
            Input a SeqIO supported file format.
        :type file_format: str
        :param export_path:
            Input a valid path to place export folder.
        :type export_path: Path
        :param folder_name:
            Input a name for the export folder.
        :type folder_name: str
        :param db_version:
            Input a db_version dictionary.
        :type db_version: dict
        :param verbose:
            Input a boolean value for verbose option.
        :type verbose: boolean
        :param data_name:
            Input a name for the file export name option.
        :type data_name: str
    """

    if verbose:
        print(
          f"Retrieving {data_name} data from {engine.url.database}...")

    if table == "phage":
        genomes = mysqldb.parse_genome_data(
                                engine,
                                phage_id_list=values,
                                phage_query="SELECT * FROM phage",
                                gene_query="SELECT * FROM gene")
    else:
        raise ValueError

    if verbose:
        print(f"Converting {data_name} data to SeqRecord format...")
    seqrecords = []

    if table == "phage":
        for gnm in genomes:
            set_cds_seqfeatures(gnm)
            if verbose:
                print(f"Converting {gnm.name}...")
            seqrecords.append(flat_files.genome_to_seqrecord(gnm))
        if verbose:
            print("Appending database version...")
        for record in seqrecords:
            append_database_version(record, db_version)

    else:
        raise ValueError

    write_seqrecord(seqrecords,
                    file_format,
                    output_path,
                    export_dir_name=output_name,
                    verbose=verbose)

def execute_csv_export(db_filter, engine,
                       output_path, output_name,
                       csv_name="database", table="phage", verbose=False):
    remove_fields = {"phage"           : ["Sequence"],
                     "gene"            : ["Translation"],
                     "domain"          : [],
                     "gene_domain"     : [],
                     "pham"            : [],
                     "pham_color"      : [],
                     "trna"            : ["Sequence"],
                     "tmrna"           : [],
                     "trna_structures" : []}

    valid_fields = db_filter.db_graph.get_table(table).show_columns()

    for unwanted_field in remove_fields[table]:
        valid_fields.remove(unwanted_field)

    if db_filter.values:
        csv_request = ("SELECT " + ",".join(valid_fields) +\
                     f" FROM {table} WHERE {db_filter.key} IN ('" + \
                      "','".join(db_filter.results(verbose=verbose)) + "')")
    else:
        csv_request = ("SELECT " + ",".join(valid_fields) + f" FROM {table}")

    if verbose:
        print(
          f"Retrieving {csv_name} data from {engine.url.database}...")
    csv_data_dicts = mysqldb.query_dict_list(engine, csv_request)

    csv_data = []
    csv_data.append(csv_data_dicts[0].keys())
    row_data = []
    for dict in csv_data_dicts:
        for key in dict.keys():
            row_data.append(dict[key])
        csv_data.append(row_data)
        row_data = []
    write_csv(csv_data, output_path,
              output_name=output_name, csv_name=csv_name,
              verbose=verbose)

def ffx_grouping(engine, group_path, group_list, db_filter,
                 db_version, file_format, table="phage", verbose=False):
    """
    Recursive helper function that handles grouping
    for ffx
    """
    current_group_list = group_list.copy()
    current_group = current_group_list.pop(0)

    groups = db_filter.group(current_group[0], current_group[1],
                             verbose=verbose)

    for group in groups.keys():
        if verbose:
            print(f"For group {current_group[1]}='{group}' "
                  f"in {current_group[0]}")
        db_filter.set_values(groups[group])

        grouped_path = group_path.joinpath(group)
        grouped_path.mkdir(exist_ok=True)

        if current_group_list:
            curr_db_filter = db_filter.copy()
            curr_db_filter.add_filter(current_group[0], current_group[1],
                                      "=", group)

            ffx_grouping(engine, grouped_path, current_group_list,
                         curr_db_filter, db_version, file_format,
                         table=table, verbose=verbose)

        else:
            execute_ffx_export(engine, db_filter.results(), file_format,
                               group_path, group, db_version, verbose=verbose,
                               data_name=f"{current_group[1]}='{group}'",
                               table=table)

def build_filter(engine, table, values_list, filters, verbose=False):
    if verbose:
        print("Building SQL data handlers...")
    db_filter = filter.Filter(engine, table=table)
    if values_list:
        db_filter.set_values(values_list)

    for filter_list in filters:
        db_filter.add_filter(filter_list[0], filter_list[1],
                             filter_list[3], filter_list[2],
                             verbose=verbose)

    if not db_filter.updated:
        db_filter.update(verbose=verbose)
        if db_filter.hits(verbose=verbose) == 0:
            print("Database returned no results.")
            exit(1)
        if verbose:
            print("")
    db_filter.sort(db_filter.key)

    return db_filter

def convert_path(path: str):
    """Function to convert a string to a working Path object.

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
    """Helper function to convert a string to a working Path object.

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
        exit(1)

def convert_file_path(path: str):
    """Helper function to convert a string to a working Path object.

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
def parse_value_list_input(value_list_input):
    """Helper function to populate the filter list for a SQL query.

    :param phage_list_input:
        Input a list of phage names.
    :type value_list_input: list[str]
    """

    print("Phage list input for database query is not a supported type")
    raise TypeError

@parse_value_list_input.register(Path)
def _(value_list_input):
    value_list = []
    with open(value_list_input, newline = '') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ",", quotechar = '|')
        for name in csv_reader:
            value_list.append(name[0])
    return value_list

@parse_value_list_input.register(list)
def _(value_list_input):
    return value_list_input

def set_cds_seqfeatures(phage_genome: genome.Genome):
    """Helper function that queries for and returns
    cds data from a SQL database for a specific phage

    :param phage_genome:
        Input a genome object to query cds data for.
    :type phage_genome: genome
    """

    try:
        def _sorting_key(cds_feature): return cds_feature.start
        phage_genome.cds_features.sort(key=_sorting_key)
    except:
        if phage_genome == None:
            raise TypeError
        print("Genome cds features unable to be sorted")
        pass
    for cds_feature in phage_genome.cds_features:
        cds_feature.set_seqfeature()

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
    version_keys = version_data.keys()
    if "Version" not in version_keys or "SchemaVersion" not in version_keys:
        print("Version of selected database "
        "is outdated.\nVersion data is incompatable")
        raise ValueError
    try:
        genome_seqrecord.annotations["comment"] =\
                genome_seqrecord.annotations["comment"] + (
                    "Database Version: {}; Schema Version: {}".format(
                                        version_data["Version"],
                                        version_data["SchemaVersion"]),)
    except:
        if isinstance(genome_seqrecord, SeqRecord):
            raise ValueError

        elif genome_seqrecord == None:
            raise TypeError
        raise

def write_seqrecord(seqrecord_list: List[SeqRecord],
                           file_format: str,
                           export_path: Path,
                           export_dir_name="export",
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
            print(f"...Writing {record.name}...")
        output_dir = f"{record.name}.{file_format}"
        output_path = export_path.joinpath(output_dir)
        output_handle = output_path.open(mode='w')
        SeqIO.write(record, output_handle, file_format)
        output_handle.close()

def write_csv(csv_data, output_path, output_name="export",csv_name="database",
              verbose=False):
    """Writes a formatted csv file from genome objects"""

    export_path = output_path.joinpath(output_name)

    if not export_path.exists():
        export_path.mkdir()

    csv_path = export_path.joinpath(f"{csv_name}.csv")
    csv_version = 1

    while(csv_path.exists()):
        csv_version += 1
        csv_path = export_path.joinpath(f"{csv_name}{csv_version}.csv")

    if verbose:
            print(f"...Writing {csv_name}.csv...")

    csv_path.touch()
    with open(csv_path, 'w', newline="") as csv_file:
        csvwriter=csv.writer(csv_file, delimiter=",",
                             quotechar="\"",
                             quoting=csv.QUOTE_MINIMAL)
        for row in csv_data:
            csvwriter.writerow(row)

def write_database(engine, version, output_path,
                    output_name="export"):

    export_path = output_path.joinpath(output_name)

    if not export_path.exists():
        export_path.mkdir()

    # TODO this is a current temporary fix.
    sql_path = export_path.joinpath(f"{engine.url.database}.sql")
    os.system(f"mysqldump -u {engine.url.username} -p{engine.url.password} "
              f"--skip-comments {engine.url.database} > {str(sql_path)}")
    version_path = sql_path.with_name(f"{engine.url.database}.version")
    version_path.touch()
    version_path.write_text(f"{version}")

def main(args):
    """Function to initialize file export"""
    run_export(args)

class Cmd_Export(cmd.Cmd):
    def __init__(self, file_format="gb",database=None,
                 phage_filter_list=[], engine=None,
                 export_directory_name="export_db",
                 export_directory_path = Path.cwd()):

        super(Cmd_Export, self).__init__()

        self.file_format = file_format
        self.database = database
        self.phage_filter_list = phage_filter_list
        self.engine = engine
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

        if self.engine == None or \
           self.engine.url.database != self.database:
            self.engine = mysqldb.connect_to_db(self.database)

        self.prompt = "({}) (export){}@localhost: ".\
                format(self.database, self.engine.url.username)

    def do_search(self, *args):
        """Filters and queries database for genomes.
        """
        db_filter = filter.Filter(self.engine)
        interactive_filter = filter.Cmd_Filter(
                db_filter=db_filter, engine=self.engine)
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
        execute_export(self.file_format, self.engine,
                                self.values_list, self.directory_path,
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

#PROTOTYPE FUNCTIONS
def cds_to_seqrecord(cds):
    try:
        record = SeqRecord(cds.seq)
        record.seq.alphabet = IUPAC.IUPACAmbiguousDNA()
    except AttributeError:
        print("Genome object failed to be converted to SeqRecord\n."
              "Genome valid attribute 'seq' is required to "
              "convert to SeqRecord object.")

    record.name = cds.id
    if cds.locus_tag != "":
        record.id = cds.locus_tag
    cds.set_seqfeature()
    record.features = [cds.seqfeature]
    record.description = f"Single gene {cds.id}"
    record.annotations = get_cds_seqrecord_annotations(cds)

    return record

def get_cds_seqrecord_annotations(cds):
    annotations = {"molecule type": "DNA",
                   "topology" : "linear",
                   "data_file_division" : "PHG",
                   "date" : "",
                   "accessions" : [],
                   "sequence_version" : "1",
                   "keyword" : [],
                   "source" : "",
                   "organism" : "",
                   "taxonomy" : [],
                   "comment" : ()}
    return annotations

def parse_cds_data_from_geneid(engine, geneid_list):
    if not geneid_list:
        return []

    query = (f"SELECT * FROM gene WHERE GeneID IN ('" + \
              "','".join(geneid_list) + "')")
    result_list = mysqldb.query_dict_list(engine, query)
    cds_list = []
    for data_dict in result_list:
        cds_list.append(mysqldb.parse_gene_table_data(data_dict))

    return cds_list

if __name__ == "__main__":

    args = sys.argv
    args.insert(0, "blank_argument")
    main(args)
