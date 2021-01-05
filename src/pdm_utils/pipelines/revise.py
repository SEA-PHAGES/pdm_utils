"""Pipeline to automate product annotation resubmissions to GenBank. """
import argparse
import logging
import sys
import time
from datetime import date
from pathlib import Path

import pdm_utils         # to get the version number
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio
from pdm_utils.functions import flat_files
from pdm_utils.functions import ncbi
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb

# -----------------------------------------------------------------------------
# GLOBAL VARIABLES

DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_revise")

BASE_CONDITIONALS = ("phage.Status = final AND "
                     "phage.AnnotationAuthor = 1 AND"
                     "phage.RetrieveRecord = 1")

REVISE_PIPELINES = ["local", "remote"]

CURATION_NAME = "revise.csv"

CURATION_HEADER = ["Phage", "Accession Number", "Locus Tag",
                   "Start", "Stop", "Product"]
TICKET_HEADER = ["table", "field", "value", "key_name", "key_value"]

FIVE_COLUMN_TABLE_HEADER = []

REVISION_COLUMNS = ["phage.PhageID", "phage.Accession", "gene.LocusTag",
                    "gene.Start", "gene.Stop", "gene.Notes", "gene.GeneID"]

LOCAL_INPUT_FILE_TYPES = ["function_report", "csv"]
LOCAL_OUTPUT_FILE_TYPES = ["ticket", "p_curation"]

REMOTE_OUTPUT_FILE_TYPES = ["tbl", "p_curation"]

INPUT_FILE_KEYS = {"function_report": {"data_key": "Pham",
                                       "filter_key": "gene.PhamID"},
                   "csv": {"data_key": "GeneID",
                           "filter_key": "gene.GeneID"}
                   }

PHAGE_QUERY = "SELECT * FROM phage"
GENE_QUERY = "SELECT * FROM gene"
TRNA_QUERY = "SELECT * FROM trna"
TMRNA_QUERY = "SELECT * FROM tmrna"


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
MAIN_LOG_FILE = "REVISE_LOG.log"
VERSION = pdm_utils.__version__
CURRENT_DATE = date.today().strftime("%Y%m%d")
# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the revise pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_revise(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    if args.pipeline == "local":
        execute_local_revise(alchemist, args.revisions_file,
                             folder_path=args.folder_path,
                             folder_name=args.folder_name,
                             config=config, input_type=args.input_type,
                             output_type=args.output_type,
                             filters=args.filters, groups=args.groups,
                             verbose=args.verbose, force=args.force,
                             production=args.production)
    elif args.pipeline == "remote":
        values = pipelines_basic.parse_value_input(args.input)
        execute_remote_revise(alchemist, folder_path=args.folder_path,
                              folder_name=args.folder_name, config=config,
                              values=values, filters=args.filters,
                              groups=args.groups, verbose=args.verbose,
                              output_type=args.output_type, force=args.force)


def parse_revise(unparsed_args_list):
    """Parses revise arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    DATABASE_HELP = """
        Name of the MySQL database to export from.
        """
    PIPELINE_HELP = """Pipeline of how to retrieve information to submit to
        GenBank
        """

    REMOTE_HELP = """Revise pipeline that accepts changes to submit to GenBank
        based on differences between data in the database and data retrieved
        from GenBank.
        """
    LOCAL_HELP = """Revise pipeline that accepts changes to submit to GenBank
        based on differences between data received from a file and the MySQL
        database.
        """

    REVISIONS_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    INPUT_TYPE_HELP = """
        Revision option that selects the input file type.
            Follow selection argument with a supported file type.
        """
    OUTPUT_TYPE_HELP = """
        Revision option that selects the output_file_type.
            Follow selection argument with a supported file type.
        """
    PRODUCTION_HELP = """
        Revise option to toggle additional filters to support production-level
        revision.
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

    FORCE_HELP = """
        Revise option that aggresively creates and overwrites directories.
        """
    CONFIG_FILE_HELP = """
        Revise option that enables use of a config file for sourcing
        credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
        """
    VERBOSE_HELP = """
        Revise option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Revise option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired directory.
        """
    FOLDER_NAME_HELP = """
        Revise option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
        """

    FILTERS_HELP = """
        Data filtering option that filters data by the inputted expressions.
            Follow selection argument with formatted filter expression:
                {Table}.{Column}={Value}
        """
    GROUP_BY_HELP = """
        Data selection option that groups data by the inputted columns.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """

    parser = argparse.ArgumentParser()

    parser.add_argument("database", type=str,  help=DATABASE_HELP)

    subparsers = parser.add_subparsers(dest="pipeline", required=True,
                                       help=PIPELINE_HELP)

    local_parser = subparsers.add_parser("local", help=LOCAL_HELP)
    remote_parser = subparsers.add_parser("remote", help=REMOTE_HELP)

    local_parser.add_argument("revisions_file", help=REVISIONS_FILE_HELP,
                              type=pipelines_basic.convert_file_path)

    local_parser.add_argument("-it", "--input_type", help=INPUT_TYPE_HELP,
                              choices=LOCAL_INPUT_FILE_TYPES)
    local_parser.add_argument("-ft", "--output_type", help=OUTPUT_TYPE_HELP,
                              choices=LOCAL_OUTPUT_FILE_TYPES)
    local_parser.add_argument("-p", "--production", help=PRODUCTION_HELP,
                              action="store_true")

    remote_parser.add_argument("-if", "--import_file", dest="input",
                               type=pipelines_basic.convert_file_path,
                               help=IMPORT_FILE_HELP)
    remote_parser.add_argument("-in", "--import_names", nargs="*",
                               dest="input", help=SINGLE_GENOMES_HELP)

    remote_parser.add_argument("-ft", "--output_type", help=OUTPUT_TYPE_HELP,
                               choices=REMOTE_OUTPUT_FILE_TYPES)

    for subparser in [local_parser, remote_parser]:
        subparser.add_argument("-f", "--force", action="store_true",
                               help=FORCE_HELP)
        subparser.add_argument("-c", "--config_file", help=CONFIG_FILE_HELP,
                               type=pipelines_basic.convert_file_path)
        subparser.add_argument("-m", "--folder_name", type=str,
                               help=FOLDER_NAME_HELP)
        subparser.add_argument("-o", "--folder_path", type=Path,
                               help=FOLDER_PATH_HELP)
        subparser.add_argument("-v", "--verbose", action="store_true",
                               help=VERBOSE_HELP)

        subparser.add_argument("-w", "--where", nargs="?", dest="filters",
                               help=FILTERS_HELP)
        subparser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                               help=GROUP_BY_HELP)

        subparser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                               folder_path=None, config_file=None, input=[],
                               input_type="function_report",
                               output_type="p_curation",
                               filters="", groups=[], verbose=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_local_revise(alchemist, revisions_file_path, folder_path=None,
                         folder_name=DEFAULT_FOLDER_NAME, config=None,
                         input_type="function_report",
                         output_type="p_curation", production=False,
                         filters="", groups=[],
                         force=False, verbose=False):
    """Executes the entirety of the genbank local revise pipeline.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param revisions_data_dicts: Data dictionaries containing pham/notes data.
    :type revisions_data_dicts: list[dict]
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param input_type: Specifies the file format of the input file
    :type input_type: str
    :param output_type: Specifies the file format of the outputted file
    :type output_type: str
    :param production: Toggles additional filters for production-level revision
    :type production: bool
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param force: A boolean to toggle aggresive building of directories.
    :type force: bool
    """
    keys = INPUT_FILE_KEYS.get(input_type)
    if keys is None:
        raise ValueError(f"Revision input type {input_type} is not supported.")

    revisions_data_dicts = fileio.retrieve_data_dict(revisions_file_path)

    values = []
    for data_dict in revisions_data_dicts:
        values.append(data_dict[keys['data_key']])

    db_filter = pipelines_basic.build_filter(alchemist, keys['filter_key'],
                                             filters, values=values,
                                             verbose=verbose)

    if production:
        db_filter.add(BASE_CONDITIONALS)

    revise_columns = db_filter.get_columns(REVISION_COLUMNS)

    if verbose:
        print("Creating export folder...")
    export_path = pipelines_basic.create_working_path(folder_path, folder_name,
                                                      force=force)

    conditionals_map = pipelines_basic.build_groups_map(db_filter, export_path,
                                                        force=force,
                                                        groups=groups,
                                                        verbose=verbose)

    if verbose:
        print("Prepared query and path structure, beginning revise export...")

    for mapped_path in conditionals_map.keys():
        conditionals = conditionals_map[mapped_path]

        if input_type == "function_report":
            export_dicts = use_function_report_data(
                                        db_filter, revisions_data_dicts,
                                        revise_columns, conditionals,
                                        verbose=verbose)
        elif input_type == "csv":
            export_dicts = use_csv_data(db_filter, revisions_data_dicts,
                                        revise_columns, conditionals,
                                        verbose=verbose)

        if not export_dicts:
            if verbose:
                print(f"'{mapped_path.name}' data selected does not require "
                      "revision; no file exported...")

            continue

        pipelines_basic.create_working_dir(mapped_path, force=force)

        write_revise_file(export_dicts, mapped_path, file_format=output_type,
                          verbose=verbose)


# TODO Owen unittest
def execute_remote_revise(alchemist, folder_path=None,
                          folder_name=DEFAULT_FOLDER_NAME, config=None,
                          output_type="p_curation", values=None,
                          filters="", groups=[], verbose=False, force=False):
    ncbi_creds = {}
    if config is not None:
        ncbi_creds = config["ncbi"]

    db_filter = pipelines_basic.build_filter(alchemist, "phage", filters,
                                             values=values, verbose=verbose)
    db_filter.add(BASE_CONDITIONALS)

    revise_path = pipelines_basic.create_working_path(folder_path, folder_name,
                                                      force=force)

    conditionals_map = pipelines_basic.build_groups_map(
                                    db_filter, revise_path, groups=groups,
                                    verbose=verbose, force=force)

    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]
        db_filter.values = db_filter.build_values(where=conditionals)

        if db_filter.hits() == 0:
            print(f"No database entries received for '{mapped_path}'.")

        pipelines_basic.create_working_dir(mapped_path, force=force)
        build_revise_log_file(mapped_path)

        logger.info(f"pdm_utils version: {VERSION}")
        logger.info(f"Revise run date: {CURRENT_DATE}")
        logger.info(f"Connected to database: {alchemist.database}")

        accession_data = db_filter.select(["phage.PhageID", "phage.Accession"])

        acc_id_dict = {}
        for data_dict in accession_data:
            accession = data_dict["Accession"]
            if not (accession is None or accession == ""):
                acc_id_dict[accession] = data_dict["PhageID"]

        tbl_records = get_tbl_records(acc_id_dict, ncbi_cred_dict=ncbi_creds)

        validated_phages = []
        for tbl_record in tbl_records:
            validated_phages.append(tbl_record.name)

        id_record_map = build_id_record_map(alchemist, validated_phages)

        if output_type == "tbl":
            revised_records = revise_seqrecords(id_record_map, tbl_records,
                                                verbose=verbose)

            if not revised_records:
                print("No discrepancies detected between "
                      f"local data and GenBank data for '{mapped_path}'.")
                continue

        elif output_type == "p_curation":
            curation_data_dicts = find_product_discrepancies(id_record_map,
                                                             tbl_records,
                                                             verbose=verbose)

            if not curation_data_dicts:
                print("No discrepancies detected between "
                      f"local data and GenBank data for '{mapped_path}'.")
                continue

        if output_type == "tbl":
            fileio.write_feature_table(revised_records, mapped_path,
                                       verbose=verbose)
        elif output_type == "p_curation":
            file_path = mapped_path.joinpath("revise.csv")
            fileio.export_data_dict(curation_data_dicts, file_path,
                                    CURATION_HEADER, include_headers=True)

# LOCAL REVISE HELPER FUNCTIONS
# -----------------------------------------------------------------------------


def use_function_report_data(db_filter, data_dicts, columns, conditionals,
                             verbose=False):
    """Reads in FunctionReport data and pairs it with existing data.

    :param db_filter: A connected and fully built Filter object.
    :type db_filter: Filter
    :param data_dicts: List of data dictionaries from a FunctionReport file.
    :type data_dicts: list[dict]
    :param columns: List of SQLAlchemy Columns to retrieve data for.
    :type columns: list[Column]
    :param conditionals: List of SQLAlchemy BinaryExpressions to filter with.
    :type conditionals: List[BinaryExpression]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    if verbose:
        print("Retreiving feature data using pham function report...")

    export_dicts = []
    for data_dict in data_dicts:
        final_call = data_dict["Final Call"]
        if final_call.lower() == "hypothetical protein":
            final_call = ""
        conditionals.append(querying.build_where_clause(
                            db_filter.graph, f"gene.Notes!='{final_call}'"))

        query = querying.build_select(db_filter.graph, columns,
                                      where=conditionals)

        results = querying.execute(db_filter.engine, query,
                                   in_column=db_filter.key,
                                   values=[data_dict["Pham"]])

        for result in results:
            if (not result["Accession"]) or (not result["LocusTag"]):
                continue
            result["Notes"] = data_dict["Final Call"]
            result["Start"] = result["Start"] + 1
            export_dicts.append(result)

    return export_dicts


# TODO Unittest
def use_csv_data(db_filter, data_dicts, columns, conditionals, verbose=False):
    """Reads in gene table csv data and pairs it with existing data.

    :param db_filter: A connected and fully built Filter object.
    :type db_filter: Filter
    :param data_dicts: List of data dictionaries from a FunctionReport file.
    :type data_dicts: list[dict]
    :param columns: List of SQLAlchemy Columns to retrieve data for.
    :type columns: list[Column]
    :param conditionals: List of SQLAlchemy BinaryExpressions to filter with.
    :type conditionals: List[BinaryExpression]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    if verbose:
        print("Retrieving feauture data using gene table csv...")

    query = querying.build_select(db_filter.graph, columns, where=conditionals)
    results = querying.execute(db_filter.engine, query,
                               in_column=db_filter.key,
                               values=db_filter.values)

    results_dict = {}
    for result in results:
        results_dict['GeneID'] = result

    export_dicts = []
    for data_dict in data_dicts:
        result_dict = results_dict.get(data_dict['GeneID'])
        if result_dict is None:
            continue
        elif result_dict["Notes"].decode("utf-8") != data_dict["Notes"]:
            result_dict["Notes"] = data_dict["Notes"]
            export_dicts.append(result_dict)

    return export_dicts


def write_revise_file(data_dicts, output_path, file_format="p_curation",
                      file_name=CURATION_NAME, verbose=False):
    """Writes a revision csv in the desired file format with necessary changes.

    :param data_dicts: List of data dictionaries to convert to curation format.
    :type data_dicts: list[dict]
    :param output_path: Path to a dir for file creation.
    :type output_path: Path
    :param file_format: Format of the csv to be written.
    :type file_format: str
    :param file_name: Name of the file to write curation data to.
    :type file_name: str
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    data_dicts = sorted(data_dicts, key=lambda data_dict: data_dict["PhageID"])

    if file_format == "p_curation":
        format_data = format_curation_data
    elif file_format == "ticket":
        format_data = format_update_ticket_data
    else:
        raise ValueError(f"File format '{file_format}' is not "
                         "recognized for revise output.")

    for d in data_dicts:
        format_data(d)

    if verbose:
        print(f"Writing {file_name} in {output_path.name}...")
    file_path = output_path.joinpath(file_name)

    if file_format == "p_curation":
        header = CURATION_HEADER
    elif file_format == "ticket":
        header = TICKET_HEADER

    fileio.export_data_dict(data_dicts, file_path, header,
                            include_headers=True)


def format_curation_data(row_dict):
    """Function to format revise dictionary keys to gb curation format.

    :param row_dict: Data dictionary for a revise file.
    :type row_dict: dict
    :param product: Gene product to append to the revise data dictionary.
    :type product: str
    """
    row_dict.pop("GeneID")

    row_dict["Phage"] = row_dict.pop("PhageID")
    row_dict["Accession Number"] = row_dict.pop("Accession")
    row_dict["Locus Tag"] = row_dict.pop("LocusTag")
    row_dict["Product"] = row_dict.pop("Notes")


def format_update_ticket_data(row_dict):
    """Function to format revise dictionary keys to update ticket format.

    :param row_dict: Data dictionary for a revise file.
    :type row_dict: dict
    :param product: Gene product to append to the revise data dictionary.
    :type product: str
    """
    row_dict.pop("PhageID")
    row_dict.pop("Accession")
    row_dict.pop("LocusTag")
    row_dict.pop("Start")
    row_dict.pop("Stop")

    product = row_dict.pop("Notes")
    if product.lower() == "hypothetical protein":
        product = ""

    row_dict["table"] = "gene"
    row_dict["field"] = "Notes"
    row_dict["value"] = product
    row_dict["key_name"] = "GeneID"
    row_dict["key_value"] = row_dict.pop("GeneID")


# REMOTE REVISE HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def build_revise_log_file(working_dir):
    log_file = working_dir.joinpath(MAIN_LOG_FILE)
    logging.basicConfig(filename=log_file, filemode="w",
                        level=logging.DEBUG,
                        format="pdm_utils revise: %(levelname)s: %(message)s")


def get_tbl_records(acc_id_dict, ncbi_cred_dict={}):
    filehandle = ncbi.get_verified_data_handle(acc_id_dict,
                                               ncbi_cred_dict=ncbi_cred_dict,
                                               file_type="tbl")
    tbl_records = []
    if filehandle is None:
        return tbl_records

    record_gen = fileio.parse_feature_table(filehandle)

    for tbl_record in record_gen:
        phage_name = acc_id_dict[tbl_record.id]
        tbl_record.name = phage_name
        flat_files.sort_seqrecord_features(tbl_record)
        tbl_records.append(tbl_record)

    return tbl_records


def build_id_record_map(alchemist, phageids):
    id_record_map = {}
    if not phageids:
        return id_record_map

    genomes = mysqldb.parse_genome_data(alchemist.engine,
                                        phage_id_list=phageids,
                                        phage_query=PHAGE_QUERY,
                                        gene_query=GENE_QUERY,
                                        trna_query=TRNA_QUERY,
                                        tmrna_query=TMRNA_QUERY)

    for genome in genomes:
        record = flat_files.genome_to_seqrecord(genome)
        id_record_map[record.id] = record

    return id_record_map


def revise_seqrecords(source_id_record_map, target_records, verbose=False):
    revised_records = []

    for target_record in target_records:
        source_record = source_id_record_map.get(target_record.id)

        if source_record is None:
            raise Exception

        edits = revise_seqrecord(target_record, source_record, verbose=verbose)
        if len(edits) > 1:
            flat_files.sort_seqrecord_features(target_record)
            revised_records.append(target_record)

        for edit_msg in edits:
            logger.info(edit_msg)

    return revised_records


def find_product_discrepancies(source_id_record_map, target_records,
                               verbose=False):
    total_product_discrepancies = []

    for target_record in target_records:
        source_record = source_id_record_map.get(target_record.id)

        if source_record is None:
            raise Exception

        record_discrepancies = curate_record_product_discrepancies(
                                target_record, source_record, verbose=verbose)
        total_product_discrepancies += record_discrepancies

    return total_product_discrepancies


def curate_record_product_discrepancies(target_record, template_record,
                                        verbose=False):
    record_discrepancies = []
    target_feature_map = create_feature_map(target_record)
    template_feature_map = create_feature_map(template_record)

    for key, target_feature_dict in target_feature_map.items():
        target_cds_feature = target_feature_dict.get("CDS")
        if target_cds_feature is None:
            continue

        temp_feature_dict = template_feature_map.get(key)
        if target_cds_feature is None:
            continue

        temp_cds_feature = temp_feature_dict.get("CDS")
        if temp_cds_feature is None:
            continue

        target_products = target_cds_feature.qualifiers["product"]
        template_products = temp_cds_feature.qualifiers["product"]
        if target_products[0] != template_products[0]:
            discrepancy = {"Phage": template_record.name,
                           "Accession Number": template_record.id,
                           "Locus Tag": key,
                           "Start": target_cds_feature.location.start+1,
                           "Stop": target_cds_feature.location.end,
                           "Product": template_products[0]}

            record_discrepancies.append(discrepancy)

    return record_discrepancies


def revise_seqrecord(target_record, template_record, verbose=False):
    """Function to edit a target record based on data from a template record.

    :param target_record: SeqRecord object to be changed based on the template.
    :type target_record: SeqRecord
    :param template_record: SeqRecord object to be used as a source of data.
    :type template_record: SeqRecord
    """
    target_feature_map = create_feature_map(target_record)
    template_feature_map = create_feature_map(template_record)

    if target_record.id != template_record.id:
        abort_msg = ("Invalid accession detected in the data for "
                     f"{template_record.name}, aborting evaluation...")
        return [abort_msg]

    record_msg = f"Evaluating {target_record.name} from GenBank..."
    if verbose:
        print(record_msg)

    edits = [record_msg]

    for key, temp_feature_dict in template_feature_map.items():
        if key == "" or key is None:
            abort_msg = ("Invalid Locus Tag detected in the local data for "
                         f"{template_record.name}, aborting evaluation...")
            return [abort_msg]

        feature_dict = target_feature_map.get(key)
        if feature_dict is None:
            add_msg = f"\tAdding '{key}' feature(s)..."
            if verbose:
                print(add_msg)
            for feature in temp_feature_dict.values():
                clean_feature(feature)
                target_record.append(feature)
            edits.append(add_msg)
        else:
            if "CDS" in temp_feature_dict.keys():
                temp_cds = feature_dict["CDS"]
                cds = temp_feature_dict.get("CDS")
                if cds is None:
                    add_msg = f"\tAdding '{key}' CDS feature..."
                    if verbose:
                        print(add_msg)

                    clean_feature(temp_cds)
                    target_record.append(temp_cds)
                    edits.append(add_msg)
            if "tRNA" in temp_feature_dict.keys():
                temp_trna = feature_dict["tRNA"]
                trna = temp_feature_dict.get("tRNA")
                if trna is None:
                    add_msg = f"\tAdding '{key}' tRNA feature..."
                    if verbose:
                        print(add_msg)

                    clean_feature(temp_cds)
                    target_record.append(temp_trna)
                    edits.append(add_msg)
            if "tmRNA" in feature_dict.keys():
                pass

    for key, feature_dict in target_feature_map.items():
        if key == "" or key is None:
            abort_msg = ("Invalid Locus Tag detected in the GenBank data for "
                         f"{target_record.name}, aborting evaluation...")
            return [abort_msg]

        temp_feature_dict = template_feature_map.get(key)
        if temp_feature_dict is None:

            remove_msg = f"\tRemoving '{key}' feature(s)..."
            if verbose:
                print(remove_msg)
            edits.append(remove_msg)
            for feature in feature_dict.values():
                target_record.features.remove(feature)
        else:
            gene = feature_dict["gene"]
            temp_gene = temp_feature_dict["gene"]
            gene_edits = revise_gene_feature(gene, temp_gene, key,
                                             verbose=verbose)
            if gene_edits:
                edits = edits + gene_edits

            if "CDS" in feature_dict.keys():
                cds = feature_dict["CDS"]
                temp_cds = temp_feature_dict.get("CDS")
                if temp_cds is None:
                    remove_msg = f"\tRemoving '{key}' CDS feature..."
                    if verbose:
                        print(remove_msg)
                    edits.append(remove_msg)
                    target_record.features.remove(cds)
                else:
                    cds_edits = revise_cds_feature(cds, temp_cds, key,
                                                   verbose=verbose)
                    if cds_edits:
                        edits = edits + cds_edits

            if "tRNA" in feature_dict.keys():
                trna = feature_dict["tRNA"]
                temp_trna = temp_feature_dict.get("tRNA")
                if temp_trna is None:
                    remove_msg = f"\tRemoving '{key}' CDS feature..."
                    if verbose:
                        print(remove_msg)
                    edits.append(remove_msg)
                    target_record.features.remove(trna)
                else:
                    trna_edits = revise_cds_feature(trna, temp_trna, key,
                                                    verbose=verbose)
                    if trna_edits:
                        edits = edits + trna_edits
            if "tmRNA" in feature_dict.keys():
                pass

    return edits


def revise_gene_feature(target_gene, template_gene, locus_tag, verbose=False):
    edits = []

    target_start = target_gene.location.start
    template_start = template_gene.location.start
    target_stop = target_gene.location.end
    template_stop = template_gene.location.end
    if target_start != template_start and target_stop == template_stop:
        start_edit_msg = (f"\tEditing '{locus_tag}' gene feature start from "
                          f"{target_start} to {template_start}...")
        if verbose:
            print(start_edit_msg)

        edits.append(start_edit_msg)
        target_gene.location = template_gene.location

    return edits


def revise_cds_feature(target_cds, template_cds, locus_tag, verbose=False):
    edits = []

    target_start = target_cds.location.start
    template_start = template_cds.location.start
    target_stop = target_cds.location.end
    template_stop = template_cds.location.end
    if target_start != template_start and target_stop == template_stop:
        start_edit_msg = (f"\tEditing '{locus_tag}' CDS feature start from "
                          f"{target_start} to {template_start}...")

        if verbose:
            print(start_edit_msg)

        edits.append(start_edit_msg)
        target_cds.location = template_cds.location

    target_products = target_cds.qualifiers["product"]
    template_products = template_cds.qualifiers["product"]
    if target_products != template_products:
        product_edit_msg = (f"\tEditing '{locus_tag}' CDS feature "
                            f"product from {';'.join(target_products)} "
                            f"to {';'.join(template_products)}")

        if verbose:
            print(product_edit_msg)

        edits.append(product_edit_msg)
        target_cds.qualifiers["product"] = template_products

    return edits


def revise_trna_feature(target_trna, template_trna, locus_tag, verbose=False):
    edits = []

    target_start = target_trna.location.start
    template_start = template_trna.location.start
    target_stop = target_trna.location.end
    template_stop = template_trna.location.end
    if target_start != template_start and target_stop == template_stop:
        start_edit_msg = (f"Editing '{locus_tag}' tRNA feature start...")

        if verbose:
            print(start_edit_msg)

        edits.append(start_edit_msg)
        target_trna.location = template_trna.location

    target_product = target_trna.qualifiers["product"][0]
    template_product = template_trna.qualifiers["product"][0]
    if target_product != template_product:
        product_edit_msg = (f"\tEditing '{locus_tag}' tRNA feature product "
                            f"from {target_product} to {template_product}...")

        if verbose:
            print(product_edit_msg)

        edits.append(product_edit_msg)
        target_trna.qualifiers["product"] = [template_product]

    return edits


def create_feature_map(record):
    """Revise helper function to map all the qualities of one locus tag.

    :param record: SeqRecord object to map gene features for.
    :type record: SeqRecord
    :returns: Returns a dictionary mapping locus_tags to features.
    :rtype: dict
    """
    cds_map = {}
    for feature in record.features:
        if feature.type == "gene":
            locus_tag = feature.qualifiers["locus_tag"][0]
            cds_map[locus_tag] = {"gene": feature}

    for locus, feature_dict in cds_map.items():
        gene_feature = feature_dict["gene"]
        for feature in record.features:
            if feature.location.end == gene_feature.location.end:
                if feature.type == "gene":
                    continue

                feature_dict[feature.type] = feature

    return cds_map


def clean_feature(feature):
    """Revise helper function to format a SeqFeature for Feature Table export.

    :param feature: Biopython SeqFeature object populated with MySQL data.
    :type feature: Biopython
    """
    qualifiers = feature.qualifiers
    temp_qualifiers = {}

    if feature.type == "gene":
        temp_qualifiers["gene"] = [qualifiers.get("gene")]
        temp_qualifiers["locus_tag"] = [qualifiers.get("locus_tag")]
    elif feature.type == "CDS":
        temp_qualifiers["product"] = [qualifiers.get("product")]
        temp_qualifiers["transl_table"] = [qualifiers.get("transl_table")]
        temp_qualifiers["note"] = [qualifiers.get("note")]
    elif feature.type == "tRNA":
        temp_qualifiers["product"] = [qualifiers.get("product")]
        temp_qualifiers["note"] = [qualifiers.get("note")]
    elif feature.type == "tmRNA":
        pass

    feature.qualifiers = temp_qualifiers


if __name__ == "__main__":
    args = sys.argv
    args.insert(0, "")
    main(args)
