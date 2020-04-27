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
import shutil
from typing import List, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import basic
from pdm_utils.functions import flat_files
from pdm_utils.functions import mysqldb
from pdm_utils.functions import parsing
from pdm_utils.functions import querying

# Valid file formats using Biopython
BIOPYTHON_CHOICES = ["gb", "fasta", "clustal", "fasta-2line", "nexus",
                       "phylip", "pir", "stockholm", "tab"]
FLAT_FILE_TABLES = ["phage"]
UNDESIRED_COLUMNS = {"phage"           : ["Sequence"],
                     "gene"            : ["Translation"],
                     "domain"          : [],
                     "gene_domain"     : [],
                     "pham"            : [],
                     "pham_color"      : [],
                     "trna"            : ["Sequence"],
                     "tmrna"           : [],
                     "trna_structures" : [],
                     "version"         : []}

# Valid Biopython formats that crash the script due to specific values in
# some genomes that can probably be fixed relatively easily and implemented.
# BIOPYTHON_CHOICES_FIXABLE = ["embl", "imgt", "seqxml"]
# Biopython formats that are not applicable.
# BIOPYTHON_CHOICES_INVALID = ["fastq", "fastq-solexa", "fastq-illumina",
#                              "phd", "sff", "qual"]
# Biopython formats that are not writable
# BIOPYTHON_CHOICES_NOT_WRITABLE = ["ig"]

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the file export pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    #Returns after printing appropriate error message from parsing/connecting.
    try:
        args = parse_export(unparsed_args_list)
    except ValueError:
        print("ABORTING EXPORT: Export failed during startup.")
        sys.exit(1)

    # TODO temp fix to allow sql export.
    if args.pipeline != "sql":
        mysqldb.check_schema_compatibility(args.alchemist.engine, "export")

    csvx = False
    ffx = None
    dbx = False
    ix = False

    values = []
    filters = []
    groups = []

    if args.pipeline in BIOPYTHON_CHOICES+["csv"]:
        values = parse_value_list_input(args.input)
        filters = args.filters
        groups = args.groups

        if args.pipeline == "csv":
            csvx = True
        else:
            ffx = args.pipeline
    elif args.pipeline == "sql":
        dbx = True
    elif args.pipeline == "I":
        ix = True
    else:
        print("ABORTED EXPORT: Unknown pipeline option discrepency.\n"
              "Pipeline parsed from command line args is not supported")
        sys.exit()

    if not ix:
        execute_export(args.alchemist, args.output_path, args.output_name,
                            values=values, verbose=args.verbose,
                            csv_export=csvx, ffile_export=ffx, db_export=dbx,
                            table=args.table,
                            filters=filters, groups=groups)
    else:
        pass

def parse_export(unparsed_args_list):
    """Parses export_db arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
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

    initial_parser = argparse.ArgumentParser()
    initial_parser.add_argument("database", type=str, help=DATABASE_HELP)

    initial_parser.add_argument("pipeline", type=str,
                               choices=export_options,
                               help=EXPORT_SELECT_HELP)

    initial = initial_parser.parse_args(unparsed_args_list[2:4])

    alchemist = establish_connection(initial.database)

    optional_parser = argparse.ArgumentParser()

    optional_parser.add_argument("-o", "--output_name",
                               type=str, help=OUTPUT_NAME_HELP)
    optional_parser.add_argument("-p", "--output_path", type=convert_dir_path,
                               help=OUTPUT_PATH_HELP)
    optional_parser.add_argument("-v", "--verbose", action="store_true",
                               help=VERBOSE_HELP)


    if initial.pipeline in (BIOPYTHON_CHOICES + ["csv"]):
        db_tables = mysqldb.get_db_tables(alchemist.engine, initial.database)
        table_choices = dict.fromkeys(BIOPYTHON_CHOICES, FLAT_FILE_TABLES)
        table_choices.update({"csv": db_tables})
        optional_parser.add_argument("-t", "--table", help=TABLE_HELP,
                                choices=table_choices[initial.pipeline])

        optional_parser.add_argument("-if", "--import_file",
                                type=convert_file_path,
                                help=IMPORT_FILE_HELP, dest="input",
                                default=[])
        optional_parser.add_argument("-in", "--import_names", nargs="*",
                                help=SINGLE_GENOMES_HELP, dest="input")
        optional_parser.add_argument("-f", "--filter", nargs="?",
                                type=parsing.parse_cmd_string,
                                help=FILTERS_HELP,
                                dest="filters")
        optional_parser.add_argument("-g", "--group", nargs="*",
                                help=GROUPS_HELP,
                                dest="groups")

    date = time.strftime("%Y%m%d")
    default_folder_name = f"{date}_export"
    default_folder_path = Path.cwd()

    optional_parser.set_defaults(alchemist=alchemist,
                        output_name=default_folder_name,
                        output_path=default_folder_path,
                        verbose=False, input=[],
                        pipeline=initial.pipeline,
                        table="phage", filters=[], groups=[],
                        ix=False, csvx=False, dbx=False, ffx=False)

    parsed_args = optional_parser.parse_args(unparsed_args_list[4:])

    return parsed_args

def execute_export(alchemist, output_path, output_name,
                        values=[], verbose=False,
                        csv_export=False, ffile_export=None, db_export=False,
                        table="phage", filters=[], groups=[]):
    """Executes the entirety of the file export pipeline.

    :param alchemist: A connected and fully built AlchemyHandler.
    :type alchemist: AlchemyHandler
    :param output_path: Path to a valid dir for new dir creation.
    :type export_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param phage_filter_list: A list of phageIDs.
    :type phage_filter_list: list[str]
    :param verbose: A boolean value for verbose option.
    :type verbose: boolean
    :param csv_export: A boolean value to toggle csv_export.
    :type csv_export: boolean
    :param ffile_export: SeqIO supported file format.
    :type ffile_export: str
    :param db_export: A boolean value to toggle db_export.
    :type db_export: boolean
    :param filters: A list of lists with filter values, grouped by ORs then AND.
    :type filters: list[list[str]]
    :param groups: A list of supported MySQL column names.
    :type groups: List[str]
    """

    if verbose:
        print("Retrieving database version...")
    db_version = mysqldb.get_version_table_data(alchemist.engine)

    if verbose:
        print("Creating export folder...")

    export_path = output_path.joinpath(output_name)
    export_path = basic.make_new_dir(output_path, export_path, attempt=50)

    if db_export:
        if verbose:
            print("Writing SQL database file...")
        write_database(alchemist, db_version["Version"], export_path)

    elif csv_export or ffile_export != None:
        try:
            db_filter = build_filtered_values(alchemist, table, filters,
                                                                values=values,
                                                                verbose=verbose)

            conditionals_map = {}
            if groups:
                if verbose:
                    print("Starting grouping process...")
                build_groups_map(db_filter, export_path, groups=groups,
                                            conditionals_map=conditionals_map,
                                            verbose=verbose)
            else:
                conditional_maps.update({export_path : \
                                         db_filter.build_where_clauses()})

        except:
            print("ABORTING EXPORT: Export failed during query structuring.")
            print("...Removing created export directory...")
            shutil.rmtree(str(export_path))
            sys.exit(1)

        try:
            if verbose:
                print("Prepared path structure, beginning export...")
            for mapped_path in conditionals_map.keys():
                conditionals = conditionals_map[mapped_path]

                if csv_export:
                    execute_csv_export(alchemist,
                                            mapped_path, export_path,
                                            table=table,
                                            conditionals=conditionals,
                                            values=values,
                                            verbose=verbose)

                elif ffile_export != None:
                    execute_ffx_export(alchemist,
                                            mapped_path, export_path,
                                            ffile_export, db_version,
                                            table=table,
                                            conditionals=conditionals,
                                            values=values,
                                            verbose=verbose)
        except:
            print("ABORTING EXPORT: Export failed during data export.")
            print("...Removing created export directory...")
            shutil.rmtree(str(export_path))
            sys.exit(1)

def establish_connection(database):
    alchemist = AlchemyHandler()

    alchemist.database = database
    try:
        alchemist.connect(login_attempts=3)
    except ValueError:
        print("Credentials invalid. "
              "Please check your MySQL credentials and try again.")
        raise ValueError
    except:
        print("Unknown error occured when logging into MySQL.")
        raise

    try:
        alchemist.connect(ask_database=True, login_attempts=0)
    except ValueError:
        print("Please check your MySQL database access, "
              "and/or your database availability.\n"
             f"Unable to connect to database '{database}' "
              "with valid credentials.")
        raise ValueError
    except:
        print("Unknown error occured when connecting to MySQL database.")
        raise
    return alchemist

def build_filtered_values(alchemist, table, filters, values=None,
                                                     verbose=False):
    table_obj = alchemist.get_table(table)
    primary_key = list(table_obj.primary_key.columns)[0]

    db_filter = Filter(alchemist=alchemist, key=primary_key)
    db_filter.values = values

    if verbose:
        print("Processing and building filters...")
    for or_filters in filters:
        for filter in or_filters:
            try:
                db_filter.add(filter)
            except:
                print(f"Filter '{filter}' is not a valid filter.")
                raise ValueError

    query = querying.build_count(alchemist.graph, primary_key,
                                       where=db_filter.build_where_clauses())
    hits = alchemist.scalar(query)

    if hits == 0:
        print("Applied filters resulted in no database hits.")
        raise ValueError

    if verbose:
        print(f"Total Database hits: {hits}")

    return db_filter

def build_groups_map(db_filter, export_path, groups=[], conditionals_map={},
                                                       verbose=False,
                                                       previous=None,
                                                       depth=0):

    groups = groups.copy()
    current_group = groups.pop(0)
    if verbose:
        if depth > 0:
            dots = ".." * depth
            print(f"{dots}Grouping by {current_group} in {previous}...")
        else:
            print(f"Grouping by {current_group}...")

    try:
        group_column = db_filter.convert_column_input(current_group)
    except:
        print(f"Group '{current_group}' is not a valid group")
        raise ValueError

    query = querying.build_distinct(db_filter.graph, group_column,
                                    where=db_filter.build_where_clauses())
    transposed_values = querying.first_column(db_filter.engine, query)

    if not transposed_values:
        export_path.rmdir()

    for group in transposed_values:
        group_path = export_path.joinpath(str(group))
        group_path.mkdir()
        db_filter_copy = db_filter.copy()
        db_filter_copy.add(f"{current_group}={group}")

        if groups:
            previous = f"{current_group} {group}"
            build_groups_map(db_filter_copy, group_path, groups=groups,
                                                    conditionals_map=\
                                                            conditionals_map,
                                                    verbose=verbose,
                                                    previous=previous,
                                                    depth=depth+1)
        else:
            previous = f"{current_group} {group}"
            conditionals_map.update({group_path : \
                                        db_filter_copy.build_where_clauses()})

def execute_csv_export(alchemist, export_path, output_path,
                                        table="phage", conditionals=None,
                                        values=[],
                                        verbose=False):
    table_obj = querying.get_table(alchemist.metadata, table)

    select_columns = []
    headers = []

    if verbose:
        relative_path = str(export_path.relative_to(output_path))
        print(f"Preparing {table} export for '{relative_path}'...")

    for column in table_obj.columns:
        if column.name not in UNDESIRED_COLUMNS[table]:
            select_columns.append(column)
            headers.append(column.name)

    primary_key = list(table_obj.primary_key.columns)[0]

    query = querying.build_select(alchemist.graph, select_columns,
                                                        where=conditionals)

    if values:
        results = querying.execute_value_subqueries(alchemist.engine, query,
                                                    primary_key, values)
    else:
        results = alchemist.execute(query)

    if verbose:
        print(f"...Writing csv '{export_path.name}.csv'...")

    file_path = export_path.joinpath(f"{export_path.name}.csv")
    basic.export_data_dict(results, file_path, headers,
                                               include_headers=True)

def execute_ffx_export(alchemist, export_path, output_path,
                       file_format, db_version, table="phage",
                       conditionals=None, verbose=False):

    if verbose:
        print(
          f"Retrieving {data_name} data from {sql_handle.database}...")

    table_obj = querying.get_table(alchemist.metadata, table)

    primary_key = list(table_obj.primary_key.columns)[0]

    if table == "phage":
        query = querying.build_distinct(alchemist.graph, primary_key,
                                                        where=conditionals,
                                                        order_by=column)
        if values:
            values = querying.first_column_value_subqueries(alchemist.engine,
                                                        query, primary_key, values)
        else:
            values = alchemist.first_column(query)

        genomes = mysqldb.parse_genome_data(
                                alchemist.engine,
                                phage_id_list=values,
                                phage_query="SELECT * FROM phage",
                                gene_query="SELECT * FROM gene")
    else:
        print(f"Unknown error occured, table '{table}' is not recognized.")
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
        print("Unknown error occured when compiling genome information.")
        raise ValueError

    write_seqrecord(seqrecords, file_format,
                    export_path, verbose=verbose)

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

    print("String input failed to be converted to a working Path object. \n"
          "Path may not exist.")
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
        print("Path input required to be a directory "
              "does not direct to a valid directory.")
        raise ValueError

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
        print("Path input required to be a file "
              "does not direct to a valid file.")
        raise ValueError

@singledispatch
def parse_value_list_input(value_list_input):
    """Helper function to populate the filter list for a SQL query.

    :param phage_list_input:
        Input a list of phage names.
    :type value_list_input: list[str]
    """

    print("Value list input for export is not a supported type.")
    raise TypeError

@parse_value_list_input.register(Path)
def _(value_list_input):
    value_list = []
    with open(value_list_input, newline = '') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ",")
        for name in csv_reader:
            value_list.append(name[0])
    return value_list

@parse_value_list_input.register(list)
def _(value_list_input):
    return value_list_input

@parse_value_list_input.register(str)
def _(value_list_input):
    return [value_list_input]

def set_cds_seqfeatures(phage_genome):
    """Helper function that queries for and returns
    cds data from a SQL database for a specific phage

    :param phage_genome:
        Input a genome object to query cds data for.
    :type phage_genome: genome
    :param sql_database_handle:
        Input a SQLAlchemy Engine object.
    :type sql_database_handle: Engine
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
        raise ValueError("Version of selected database is outdated. "
                         "Version data is incompatable.")
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
        print("Writing selected data to files...")

    for record in seqrecord_list:
        if verbose:
            print(f"...Writing {record.name}...")
        file_name = f"{record.name}.{file_format}"
        output_path = export_path.joinpath(file_name)
        output_handle = output_path.open(mode='w')
        SeqIO.write(record, output_handle, file_format)
        output_handle.close()

def write_database(alchemist, version, export_path):
    sql_path = export_path.joinpath(f"{alchemist.database}.sql")
    os.system(f"mysqldump -u {alchemist.username} -p{alchemist.password} "
              f"--skip-comments {alchemist.database} > {str(sql_path)}")
    version_path = sql_path.with_name(f"{alchemist.database}.version")
    version_path.touch()
    version_path.write_text(f"{version}")


#--------------------------------------------------------
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

def parse_cds_data_from_geneid(sql_handle, geneid_list):
    if not geneid_list:
        return []

    query = (f"SELECT * FROM gene WHERE GeneID IN ('" + \
              "','".join(geneid_list) + "')")
    result_list = sql_handle.execute_query(query)

    cds_list = []
    for data_dict in result_list:
        cds_list.append(mysqldb.parse_gene_table_data(data_dict))

    return cds_list


    args = sys.argv
    args.insert(0, "blank_argument")
    main(args)
