"""Pipeline for exporting database information into files."""
import argparse
import shutil
import sys
import time
from pathlib import Path

from pdm_utils.classes.filter import Filter
from pdm_utils.functions import (configfile, fileio, flat_files, mysqldb,
                                 mysqldb_basic, pham_alignment,
                                 pipelines_basic, querying)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_export"
DEFAULT_TABLE = "phage"

TEMP_DIR = "/tmp/pdm_utils_export_temp"

PHAGE_QUERY = "SELECT * FROM phage"
GENE_QUERY = "SELECT * FROM gene"
TRNA_QUERY = "SELECT * FROM trna"
TMRNA_QUERY = "SELECT * FROM tmrna"

# Valid Biopython formats that crash the script due to specific values in
# some genomes that can probably be fixed relatively easily and implemented.
# BIOPYTHON_PIPELINES_FIXABLE = ["embl", "imgt", "seqxml"]
# Biopython formats that are not applicable.
# BIOPYTHON_PIPELINES_INVALID = ["fastq", "fastq-solexa", "fastq-illumina",
#                              "phd", "sff", "qual"]
# Biopython formats that are not writable
# BIOPYTHON_PIPELINES_NOT_WRITABLE = ["ig"]

# Valid file formats using Biopython
BIOPYTHON_PIPELINES = ["gb", "fasta", "clustal", "fasta-2line", "nexus",
                       "phylip", "pir", "stockholm", "tab"]
FILTERABLE_PIPELINES = BIOPYTHON_PIPELINES + ["csv", "tbl"]
PIPELINES = FILTERABLE_PIPELINES + ["sql"]
FLAT_FILE_TABLES = ["phage", "gene"]
FIVE_COLUMN_TABLES = ["phage"]

CDD_DATA_COLUMNS = ["gene_domain.QueryStart", "gene_domain.QueryEnd",
                    "domain.DomainID", "domain.Name", "domain.Description"]

# Once trna has data, these tables can be reintroduced.
TABLES = ["phage", "gene", "domain", "gene_domain", "pham",
          # "trna", "tmrna", "trna_structures",
          "version"]
SEQUENCE_COLUMNS = {"phage": ["Sequence"],
                    "gene": ["Translation"],
                    "domain": [],
                    "gene_domain": [],
                    "pham": [],
                    "trna": ["Sequence"],
                    "tmrna": [],
                    "trna_structures": [],
                    "version": []}


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the file export pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    # Returns after printing appropriate error message from parsing/connecting.
    args = parse_export(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    # Exporting as a SQL file is not constricted by schema version.
    if args.pipeline != "sql":
        mysqldb.check_schema_compatibility(alchemist.engine, "export")

    values = None
    if args.pipeline in FILTERABLE_PIPELINES:
        values = pipelines_basic.parse_value_input(args.input)
        if not values:
            values = None

    if args.pipeline not in PIPELINES:
        print("ABORTED EXPORT: Unknown pipeline option discrepency.\n"
              "Pipeline parsed from command line args is not supported")
        sys.exit(1)

    if args.pipeline != "I":
        execute_export(alchemist, args.pipeline, folder_path=args.folder_path,
                       folder_name=args.folder_name, table=args.table,
                       values=values, filters=args.filters, groups=args.groups,
                       sort=args.sort, include_columns=args.include_columns,
                       exclude_columns=args.exclude_columns,
                       sequence_columns=args.sequence_columns,
                       raw_bytes=args.raw_bytes,
                       concatenate=args.concatenate, db_name=args.db_name,
                       verbose=args.verbose, dump=args.dump, force=args.force,
                       threads=args.number_processes, phams_out=args.phams_out)
    else:
        pass


def parse_export(unparsed_args_list):
    """Parses export_db arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    EXPORT_SELECT_HELP = """
        Select a export pipeline option to export database data.
            Select csv to export data from a table into a .csv file.
            Select sql to dump the current database into a .sql file.
            Select a formatted file option to export individual entries.
        """
    DATABASE_HELP = "Name of the MySQL database to export from."

    CONFIG_FILE_HELP = """
        Export option that enables use of a config file for sourcing
        credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
        """
    DUMP_HELP = """
        Export option that dumps exported files directly to the desired
        working directory.
        """
    FORCE_HELP = """
        Export option that aggresively creates and overwrites directories.
        """
    VERBOSE_HELP = """
        Export option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Export option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired export directory.
        """
    FOLDER_NAME_HELP = """
        Export option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
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
    WHERE_HELP = """
        Data filtering option that filters data by the inputted expressions.
            Follow selection argument with formatted filter expression:
                {Table}.{Column}={Value}
        """
    GROUP_BY_HELP = """
        Data selection option that groups data by the inputted columns.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """
    ORDER_BY_HELP = """
        Data selection option that sorts data by the inputted columns.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """

    DB_NAME_HELP = """
        MySQL export option to allow renaming of the exported database.
            Follow selection argument with the name of the desired database.
        """

    CONCATENATE_HELP = """
        SeqRecord export option to toggle concatenation of files.
            Toggle to enable concatenation of files.
        """

    SEQUENCE_COLUMNS_HELP = """
        Csv export option to toggle removal of sequence-based data.
            Toggle to include sequence based data.
        """
    INCLUDE_COLUMNS_HELP = """
        Csv export option to add additional columns of a database
        to an exported csv file.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """
    EXCLUDE_COLUMNS_HELP = """
        Csv export option to exclude select columns of a database
        from an exported csv file.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """
    RAW_BYTES_HELP = """
        Csv export option to conserve blob and encoded data from the database
        when exporting to a csv file.
            Toggle to leave conserve format of bytes-type data.
        """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str, help=DATABASE_HELP)

    subparsers = parser.add_subparsers(dest="pipeline", required=True,
                                       help=EXPORT_SELECT_HELP)

    filterable_parsers = []
    biopython_parsers = []
    subparser_list = []
    for pipeline in BIOPYTHON_PIPELINES:
        bio_parser = (subparsers.add_parser(pipeline))
        biopython_parsers.append(bio_parser)
        filterable_parsers.append(bio_parser)
        subparser_list.append(bio_parser)

    tbl_parser = subparsers.add_parser("tbl")
    csv_parser = subparsers.add_parser("csv")
    sql_parser = subparsers.add_parser("sql")

    subparser_list.append(tbl_parser)
    subparser_list.append(csv_parser)
    subparser_list.append(sql_parser)

    filterable_parsers.append(tbl_parser)
    filterable_parsers.append(csv_parser)

    for subparser in subparser_list:
        subparser.add_argument("-c", "--config_file",
                               type=pipelines_basic.convert_file_path,
                               help=CONFIG_FILE_HELP)
        subparser.add_argument("-m", "--folder_name", type=str,
                               help=FOLDER_NAME_HELP)
        subparser.add_argument("-o", "--folder_path", type=Path,
                               help=FOLDER_PATH_HELP)
        subparser.add_argument("-v", "--verbose", action="store_true",
                               help=VERBOSE_HELP)
        subparser.add_argument("-d", "--dump", action="store_true",
                               help=DUMP_HELP)
        subparser.add_argument("-f", "--force", action="store_true",
                               help=FORCE_HELP)
        subparser.add_argument("-np", "--number_processes", type=int)

    for subparser in filterable_parsers:
        table_choices = dict.fromkeys(BIOPYTHON_PIPELINES, FLAT_FILE_TABLES)
        table_choices["csv"] = TABLES
        table_choices["tbl"] = FIVE_COLUMN_TABLES

        subparser.add_argument("-if", "--import_file",
                               type=pipelines_basic.convert_file_path,
                               help=IMPORT_FILE_HELP, dest="input",
                               default=[])
        subparser.add_argument("-in", "--import_names", nargs="*",
                               help=SINGLE_GENOMES_HELP, dest="input")
        subparser.add_argument("-w", "--where", nargs="?",
                               help=WHERE_HELP, dest="filters")
        subparser.add_argument("-g", "--group_by", nargs="*",
                               help=GROUP_BY_HELP, dest="groups")
        subparser.add_argument("-s", "--order_by", nargs="*",
                               help=ORDER_BY_HELP)

    for subparser in biopython_parsers:
        subparser.add_argument("-cc", "--concatenate", help=CONCATENATE_HELP,
                               action="store_true")
        subparser.add_argument("-t", "--table", help=TABLE_HELP,
                               choices=FLAT_FILE_TABLES)

    tbl_parser.add_argument("-t", "--table", help=TABLE_HELP,
                            choices=FIVE_COLUMN_TABLES)

    csv_parser.add_argument("-t", "--table", help=TABLE_HELP,
                            choices=TABLES)
    csv_parser.add_argument("-sc", "--sequence_columns",
                            help=SEQUENCE_COLUMNS_HELP, action="store_true")
    csv_parser.add_argument("-ic", "--include_columns", nargs="*",
                            help=INCLUDE_COLUMNS_HELP)
    csv_parser.add_argument("-ec", "--exclude_columns", nargs="*",
                            help=EXCLUDE_COLUMNS_HELP)
    csv_parser.add_argument("-rb", "--raw_bytes", help=RAW_BYTES_HELP,
                            action="store_true")

    sql_parser.add_argument("-n", "--db_name", type=str, help=DB_NAME_HELP)
    sql_parser.add_argument("-pho", "--phams_out", action="store_true")

    for subparser in subparser_list:
        subparser.set_defaults(
                        folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[],
                        table=DEFAULT_TABLE, filters="", groups=[], sort=[],
                        include_columns=[], exclude_columns=[],
                        sequence_columns=False, concatenate=False,
                        raw_bytes=False, db_name=None, phams_out=False,
                        number_processes=1)

    parsed_args = parser.parse_args(unparsed_args_list[2:])

    return parsed_args


def execute_export(alchemist, pipeline, folder_path=None,
                   folder_name=DEFAULT_FOLDER_NAME, values=None, verbose=False,
                   dump=False, force=False, table=DEFAULT_TABLE, filters="",
                   groups=[], sort=[], include_columns=[], exclude_columns=[],
                   sequence_columns=False, raw_bytes=False, concatenate=False,
                   db_name=None, phams_out=False, threads=1):
    """Executes the entirety of the file export pipeline.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param pipeline: File type that dictates data processing.
    :type pipeline: str
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param force: A boolean to toggle aggresive building of directories.
    :type force: bool
    :param values: List of values to filter database results.
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param dump: A boolean value to toggle dump in current working dir.
    :type dump: bool
    :param table: MySQL table name.
    :type table: str
    :param filters: A list of lists with filter values, grouped by ORs.
    :type filters: str
    :param groups: A list of supported MySQL column names to group by.
    :type groups: list[str]
    :param sort: A list of supported MySQL column names to sort by.
    :type sort: list[str]
    :param include_columns: A csv export column selection parameter.
    :type include_columns: list[str]
    :param exclude_columns: A csv export column selection parameter.
    :type exclude_columns: list[str]
    :param sequence_columns: A boolean to toggle inclusion of sequence data.
    :type sequence_columns: bool
    :param concatenate: A boolean to toggle concaternation for SeqRecords.
    :type concaternate: bool
    :param threads: Number of processes/threads to spawn during the pipeline
    :type threads: int
    """
    if verbose:
        print("Retrieving database version...")
    db_version = mysqldb_basic.get_first_row_data(alchemist.engine, "version")

    if pipeline == "csv":
        if verbose:
            print("Processing columns for csv export...")
        csv_columns = filter_csv_columns(alchemist, table,
                                         include_columns=include_columns,
                                         exclude_columns=exclude_columns,
                                         sequence_columns=sequence_columns)

    if pipeline in FILTERABLE_PIPELINES:
        db_filter = pipelines_basic.build_filter(alchemist, table, filters,
                                                 values=values,
                                                 verbose=verbose)
        if sort:
            pipelines_basic.add_sort_columns(db_filter, sort, verbose=verbose)

    if verbose:
        print("Creating export folder...")
    export_path = pipelines_basic.create_working_path(folder_path, folder_name,
                                                      dump=dump, force=force)

    data_cache = {}
    if pipeline == "sql":
        execute_sql_export(alchemist, export_path, folder_path, db_version,
                           db_name=db_name, dump=dump, force=force,
                           phams_out=phams_out, threads=threads,
                           verbose=verbose)
    elif pipeline in FILTERABLE_PIPELINES:
        conditionals_map = pipelines_basic.build_groups_map(
                                                db_filter, export_path,
                                                groups=groups,
                                                verbose=verbose, force=force)

        if verbose:
            print("Prepared query and path structure, beginning export...")

        values = db_filter.values
        for mapped_path in conditionals_map.keys():
            db_filter.reset()
            db_filter.values = values

            conditionals = conditionals_map[mapped_path]
            db_filter.values = db_filter.build_values(where=conditionals)

            if db_filter.hits() == 0:
                print(f"No database entries received from {table} "
                      f"for '{mapped_path}'.")
                continue

            if sort:
                sort_columns = get_sort_columns(alchemist, sort)
                db_filter.sort(sort_columns)

            export_name = None
            if dump:
                if mapped_path == export_path:
                    export_name = folder_name

            pipelines_basic.create_working_dir(mapped_path, dump=dump,
                                               force=force)

            if pipeline in BIOPYTHON_PIPELINES + ["tbl"]:
                execute_ffx_export(alchemist, mapped_path, export_path,
                                   db_filter.values, pipeline, db_version,
                                   table, concatenate=concatenate,
                                   data_cache=data_cache,
                                   export_name=export_name, threads=threads,
                                   verbose=verbose, dump=dump)
            elif pipeline == "csv":
                execute_csv_export(db_filter, mapped_path, export_path,
                                   csv_columns, table, raw_bytes=raw_bytes,
                                   data_cache=data_cache,
                                   verbose=verbose, dump=dump)
    else:
        print("Unrecognized export pipeline, aborting export")
        sys.exit(1)


def execute_csv_export(db_filter, export_path, folder_path, columns, csv_name,
                       data_cache=None, sort=[], raw_bytes=False,
                       verbose=False, dump=False):
    """Executes csv export of a MySQL database table with select columns.

    :param db_filter: A connected and fully built Filter object.
    :type db_filter: Filter
    :param export_path: Path to a dir for file creation.
    :type export_path: Path
    :param folder_path: Path to a top-level dir.
    :type folder_path: Path
    :param table: MySQL table name.
    :type table: str
    :param conditionals: MySQL WHERE clause-related SQLAlchemy objects.
    :type conditionals: list[BinaryExpression]
    :param sort: A list of SQLAlchemy Columns to sort by.
    :type sort: list[Column]
    :param values: List of values to fitler database results.
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param dump: A boolean value to toggle dump in current working dir.
    :type dump: bool
    """
    if data_cache is None:
        data_cache = {}

    if verbose:
        relative_path = str(export_path.relative_to(folder_path))
        print(f"Preparing {csv_name} export for '{relative_path}'...")

    headers = [db_filter._key.name]
    for column in columns:
        if column.name != db_filter._key.name:
            headers.append(column.name)

    results = db_filter.select(columns)

    if not raw_bytes:
        decode_results(results, columns, verbose=verbose)

    if len(results) == 0:
        print(f"No database entries received for {csv_name}.")
        if not dump:
            export_path.rmdir()
    else:
        if dump:
            if export_path != folder_path:
                export_path.rmdir()
                export_path = export_path.parent
        if verbose:
            print(f"...Writing csv {csv_name}.csv in '{export_path.name}'...")

        file_path = export_path.joinpath(f"{csv_name}.csv")
        fileio.export_data_dict(results, file_path, headers,
                                include_headers=True)


def execute_ffx_export(alchemist, export_path, folder_path, values,
                       file_format, db_version, table, concatenate=False,
                       data_cache=None, verbose=False, dump=False,
                       threads=1, export_name=None):
    """Executes SeqRecord export of the compilation of data from a MySQL entry.

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param export_path: Path to a dir for file creation.
    :type export_path: Path
    :param folder_path: Path to a top-level dir.
    :type folder_path: Path
    :param file_format: Biopython supported file type.
    :type file_format: str
    :param db_version: Dictionary containing database version information.
    :type db_version: dict
    :param table: MySQL table name.
    :type table: str
    :param values: List of values to fitler database results.
    :type values: list[str]
    :param conditionals: MySQL WHERE clause-related SQLAlchemy objects.
    :type conditionals: list[BinaryExpression]
    :param sort: A list of SQLAlchemy Columns to sort by.
    :type sort: list[Column]
    :param concatenate: A boolean to toggle concatenation of SeqRecords.
    :type concaternate: bool
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    if data_cache is None:
        data_cache = {}

    if export_name is None:
        export_name = export_path.name

    if verbose:
        print(f"Retrieving {export_name} data...")

    if table == "phage":
        seqrecords = get_genome_seqrecords(alchemist, values,
                                           data_cache=data_cache,
                                           verbose=verbose)
    elif table == "gene":
        seqrecords = get_cds_seqrecords(alchemist, values,
                                        data_cache=data_cache, verbose=verbose,
                                        file_format=file_format)
    else:
        print(f"Unknown error occured, table '{table}' is not recognized "
              "for SeqRecord export pipelines.")
        sys.exit(1)

    if file_format == "tbl":
        fileio.write_feature_table(seqrecords, export_path, verbose=verbose)
    else:
        if verbose:
            print("Appending database version...")
        for record in seqrecords:
            append_database_version(record, db_version)
        fileio.write_seqrecords(seqrecords, file_format, export_path,
                                export_name=export_name, verbose=verbose,
                                concatenate=concatenate, threads=threads)


def execute_sql_export(alchemist, export_path, folder_path, db_version,
                       db_name=None, dump=False, force=False, phams_out=False,
                       threads=1, verbose=False):
    pipelines_basic.create_working_dir(export_path, dump=dump, force=force)

    if phams_out:
        temp_dir = Path(TEMP_DIR)
        if temp_dir.is_dir():
            shutil.rmtree(temp_dir)
        temp_dir.mkdir()

        phams_out_fasta_dir = temp_dir.joinpath("fastas")
        pipelines_basic.create_working_dir(phams_out_fasta_dir,
                                           dump=dump, force=force)
        phams_out_aln_dir = temp_dir.joinpath("alns")
        pipelines_basic.create_working_dir(phams_out_aln_dir,
                                           dump=dump, force=force)

        phams_dict = pham_alignment.get_all_pham_gene_translations(alchemist)

        if verbose:
            print("...Writing and aligning pham fasta files...")
        pham_alignment.write_phams(phams_out_fasta_dir, phams_out_aln_dir,
                                   phams_dict, cores=threads, verbose=verbose)

        pham_fastas_zip = export_path.joinpath("fastas.zip")
        pham_alns_zip = export_path.joinpath("alns.zip")

        shutil.make_archive(pham_fastas_zip.with_suffix(""), "zip",
                            temp_dir, phams_out_fasta_dir.name)
        shutil.make_archive(pham_alns_zip.with_suffix(""), "zip",
                            temp_dir, phams_out_aln_dir.name)

    if verbose:
        print("Writing SQL database file...")

    fileio.write_database(alchemist, db_version["Version"], export_path,
                          db_name=db_name)


# EXPORT-SPECIFIC HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# TODO Document and Unittest
def get_genome_seqrecords(alchemist, values, data_cache=None, verbose=False):
    if data_cache is None:
        data_cache = {}

    seqrecords = []
    for genome_id in values:
        genome = data_cache.get(genome_id)
        if genome is None:
            genome = get_single_genome(alchemist, genome_id, get_features=True,
                                       data_cache=data_cache)

        seqrecord = flat_files.genome_to_seqrecord(genome)
        flat_files.sort_seqrecord_features(seqrecord)
        seqrecords.append(seqrecord)

    return seqrecords


# TODO Document and Unittest
def get_cds_seqrecords(alchemist, values, data_cache=None, nucleotide=False,
                       verbose=False, file_format=None):
    if data_cache is None:
        data_cache = {}

    cds_list = parse_feature_data(alchemist, values=values)

    db_filter = Filter(alchemist)
    db_filter.key = 'gene.GeneID'

    if verbose:
        print("...Converting SQL data...")

    seqrecords = []
    for cds in cds_list:
        parent_genome = data_cache.get(cds.genome_id)

        if parent_genome is None:
            parent_genome = get_single_genome(alchemist, cds.genome_id,
                                              data_cache=data_cache)

        cds.genome_length = parent_genome.length
        cds.set_seqfeature()

        db_filter.values = [cds.id]
        gene_domains = db_filter.select(CDD_DATA_COLUMNS)

        record = flat_files.cds_to_seqrecord(cds, parent_genome,
                                             gene_domains=gene_domains,
                                             desc_type=file_format)
        seqrecords.append(record)

    return seqrecords


def get_sort_columns(alchemist, sort_inputs):
    """Function that converts input for sorting to SQLAlchemy Columns.

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param sort_inputs: A list of supported MySQL column names.
    :type sort_inputs: list[str]
    :returns: A list of SQLAlchemy Column objects.
    :rtype: list[Column]
    """
    sort_columns = []
    for sort_input in sort_inputs:
        try:
            sort_column = querying.get_column(alchemist.metadata, sort_input)
        except ValueError:
            print("Error occured while selecting sort columns.")
            print(f"Column inputted, '{sort_input}', is invalid.")
            sys.exit(1)
        finally:
            sort_columns.append(sort_column)

    return sort_columns


# CSV-EXPORT HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def filter_csv_columns(alchemist, table, include_columns=[],
                       exclude_columns=[], sequence_columns=False):
    """Function that filters and constructs a list of Columns to select.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param table: MySQL table name.
    :type table: str
    :param include_columns: A list of supported MySQL column names.
    :type include_columns: list[str]
    :param exclude_columns: A list of supported MySQL column names.
    :type exclude_columns: list[str]
    :param sequence_columns: A boolean to toggle inclusion of sequence data.
    :type sequence_columns: bool
    :returns: A list of SQLAlchemy Column objects.
    :rtype: list[Column]
    """
    table_obj = alchemist.metadata.tables[table]
    starting_columns = list(table_obj.columns)
    primary_key = list(table_obj.primary_key.columns)[0]

    include_column_objs = starting_columns
    for column in include_columns:
        try:
            column_obj = querying.get_column(alchemist.metadata, column)
        except ValueError:
            print("Error occured while selecting csv columns.")
            print(f"Column inputted, '{column}', is invalid.")
            sys.exit(1)
        finally:
            if column_obj not in include_column_objs:
                include_column_objs.append(column_obj)

    sequence_column_objs = []
    if not sequence_columns:
        for sequence_column in SEQUENCE_COLUMNS[table]:
            sequence_column_obj = dict(table_obj.c)[sequence_column]
            sequence_column_objs.append(sequence_column_obj)

    exclude_column_objs = sequence_column_objs
    for column in exclude_columns:
        try:
            column_obj = querying.get_column(alchemist.metadata, column)
        except ValueError:
            print("Error occured while selecting csv columns.")
            print(f"Column inputted, '{column}', is invalid.")
            sys.exit(1)
        finally:
            exclude_column_objs.append(column_obj)
            if column_obj.compare(primary_key):
                print(f"Primary key to {table} cannot be excluded")
                sys.exit(1)

            if column_obj not in exclude_column_objs:
                exclude_column_objs.append(column_obj)

    columns = []
    for column_obj in include_column_objs:
        if column_obj not in exclude_column_objs:
            columns.append(column_obj)

    return columns


def decode_results(results, columns, verbose=False):
    """Function that decodes encoded results from SQLAlchemy generated data.

    :param results: List of data dictionaries from a SQLAlchemy results proxy.
    :type results: list[dict]
    :param columns: SQLAlchemy Column objects.
    :type columns: list[Column]
    """
    for column in columns:
        if column.type.python_type == bytes:
            if verbose:
                print(f"Decoding retrieved {column} data...")
            for result in results:
                if not result[column.name] is None:
                    result[column.name] = result[column.name].decode("utf-8")


# Functions to be evaluated for another module:
# -----------------------------------------------------------------------------

# Copy of mysqldb.parse_feature_data() with value subquerying.
# Value subquerying needed for +300000 GeneID entries.
# Unsure about function redundancy.
def parse_feature_data(alchemist, values=[], limit=8000):
    """Returns Cds objects containing data parsed from a MySQL database.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param values: List of GeneIDs upon which the query can be conditioned.
    :type values: list[str]
    """
    gene_table = querying.get_table(alchemist.metadata, "gene")
    primary_key = list(gene_table.primary_key.columns)[0]
    cds_data_columns = list(gene_table.c)

    cds_data_query = querying.build_select(alchemist.graph, cds_data_columns)

    cds_data = querying.execute(alchemist.engine, cds_data_query,
                                in_column=primary_key, values=values,
                                limit=limit)

    cds_list = []
    for data_dict in cds_data:
        cds_ftr = mysqldb.parse_gene_table_data(data_dict)
        cds_list.append(cds_ftr)

    return cds_list


def append_database_version(genome_seqrecord, version_data):
    """Function that appends the database version to the SeqRecord comments.

    :param genome_seqrecord: Filled SeqRecord object.
    :type genome_seqfeature: SeqRecord
    :param version_data: Dictionary containing database version information.
    :type version_data: dict
    """
    version_keys = version_data.keys()
    version = "NULL"
    schema_version = "NULL"
    if "Version" in version_keys or "SchemaVersion" not in version_keys:
        version = version_data["Version"]
    if "SchemaVersion" in version_keys:
        schema_version = version_data["SchemaVersion"]

    genome_seqrecord.annotations["comment"] =\
        genome_seqrecord.annotations["comment"] + (
            "Database Version: {}; Schema Version: {}".format(
                                                version, schema_version),)


# TODO Document
# TODO Evaluate for redundancy
def get_single_genome(alchemist, phageid, get_features=False, data_cache=None):
    gene_query = None
    trna_query = None
    tmrna_query = None
    if get_features:
        gene_query = GENE_QUERY
        trna_query = TRNA_QUERY
        tmrna_query = TMRNA_QUERY

    genome = mysqldb.parse_genome_data(
                            alchemist.engine, phage_id_list=[phageid],
                            phage_query=PHAGE_QUERY, gene_query=gene_query,
                            trna_query=trna_query, tmrna_query=tmrna_query)[0]

    if data_cache is not None:
        data_cache[phageid] = genome

    return genome
