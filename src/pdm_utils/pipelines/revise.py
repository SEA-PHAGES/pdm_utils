"""Pipeline to automate product annotation resubmissions to GenBank. """
import argparse
import os
import sys
import time
from pathlib import Path

from sqlalchemy import select

from pdm_utils.functions import basic
from pdm_utils.functions import fileio
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb_basic

#-----------------------------------------------------------------------------
#GLOBAL VARIABLES

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_revise"
DEFAULT_FOLDER_PATH = Path.cwd()
CURATION_NAME = "revise.csv"

CURATION_HEADER = ["Phage", "Accession Number", "Locus Tag", 
                   "Start", "Stop", "Product"]
FIVE_COLUMN_TABLE_HEADER = []


REVISION_COLUMNS = ["phage.PhageID", "phage.Accession", "gene.LocusTag",
                    "gene.Start", "gene.Stop", "gene.Notes"]

INPUT_FILE_TYPES = ["function_report", "csv"]
OUTPUT_FILE_TYPES = ["curation"]

INPUT_FILE_KEYS = {"function_report"   :\
                                       {"data_key"   : "Pham",
                                        "filter_key" : "gene.PhamID"},
                   "csv"               :\
                                       {"data_key"   : "GeneID",
                                        "filter_key" : "gene.GeneID"}
                  }


BASE_CONDITIONALS = ("phage.Status = final AND "
                     "phage.AnnotationAuthor = 1 AND"
                     "phage.RetrieveRecord = 1")

#-----------------------------------------------------------------------------
#MAIN FUNCTIONS

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the revise pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_revise(unparsed_args_list)
    
    alchemist = pipelines_basic.build_alchemist(args.database)

    execute_revise(alchemist, args.revisions_file, args.folder_path, 
                                                   args.folder_name,
                                                   input_type=args.input_type,
                                                   output_type=args.output_type,
                                                   filters=args.filters,
                                                   groups=args.groups,
                                                   verbose=args.verbose)

def parse_revise(unparsed_args_list):
    """Parses revise arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    DATABASE_HELP = """
        Name of the MySQL database to export from.
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

    VERBOSE_HELP = """
        Export option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Export option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired directory.
        """
    FOLDER_NAME_HELP = """
        Export option to change the name
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
    parser.add_argument("revisions_file", help=REVISIONS_FILE_HELP,
                                    type=pipelines_basic.convert_file_path)

    parser.add_argument("-m", "--folder_name", 
                                    type=str,  help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path", 
                                    type=pipelines_basic.convert_dir_path,
                                               help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true", 
                                               help=VERBOSE_HELP)

    parser.add_argument("-it", "--input_type", choices=INPUT_FILE_TYPES,
                                               help=INPUT_TYPE_HELP)
    parser.add_argument("-ot", "--output_type", choices=OUTPUT_FILE_TYPES,
                                               help=OUTPUT_TYPE_HELP)

    parser.add_argument("-f", "--filter", nargs="?", dest="filters",
                                               help=FILTERS_HELP)
    parser.add_argument("-g", "--group_by", nargs="*",
                                help=GROUP_BY_HELP,
                                dest="groups")
    
    date = time.strftime("%Y%m%d")
    default_folder_name = f"{date}_pham_revise"
    default_folder_path = Path.cwd()

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH, 
                        input_type="function_report",
                        output_type="curation",
                        filters="", groups=[], verbose=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args

def execute_revise(alchemist, revisions_file_path, folder_path, folder_name,
                                                   input_type="function_report",
                                                   output_type="curation",
                                                   filters="", groups=[],
                                                   verbose=False):
    """Executes the entirety of the genbank revise pipeline.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param revisions_data_dicts: Data dictionaries containing pham/notes data.
    :type revisions_data_dicts: list[dict]
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    keys = INPUT_FILE_KEYS.get(input_type)
    if keys is None: 
        raise ValueError(f"Revision input type {input_type} is not supported.")

    revisions_data_dicts = fileio.retrieve_data_dict(revisions_file_path)
    
    values = []
    for data_dict in revisions_data_dicts:
        values.append(data_dict[keys['data_key']])

    db_filter = pipelines_basic.build_filter(alchemist, keys['filter_key'], 
                                                            filters, 
                                                            values=values,
                                                            verbose=verbose)
    db_filter.add(BASE_CONDITIONALS)
    
    revise_columns = db_filter.get_columns(REVISION_COLUMNS)
     
    if verbose:
        print("Creating export folder...")
    export_path = folder_path.joinpath(folder_name)
    export_path = basic.make_new_dir(folder_path, export_path, attempt=50)

    conditionals_map = {}
    pipelines_basic.build_groups_map(db_filter, export_path, conditionals_map,
                                                         groups=groups,
                                                         verbose=verbose)

    if verbose:
        print("Prepared query and path structure, beginning review export...")

    for mapped_path in conditionals_map.keys():
        conditionals = conditionals_map[mapped_path]

        if input_type == "function_report":
            export_dicts = use_function_report_data(
                                            db_filter, revisions_data_dicts, 
                                            revise_columns, conditionals, 
                                            verbose=verbose)
        elif input_type == "csv":
            export_dicts = use_csv_data(    db_filter, revisions_data_dicts,
                                            revise_columns, conditionals,
                                            verbose=verbose)
        
        if not export_dicts:
            if verbose:
                print("'{mapped_path.name}' data selected does not require "
                      "revision; no file exported...")

            mapped_path.rmdir()
            continue

        if output_type == "curation":
            write_curation_data(export_dicts, mapped_path)
        elif output_type == "five_column":
            write_five_column_table(export_dicts, mapped_path)
        
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
        if verbose:
            print(f"...Retrieving data for pham {data_dict['Pham']}...")

        final_call = data_dict["Final Call"]
        if final_call == "Hypothetical Protein":
            final_call = ""
        conditionals.append(querying.build_where_clause(db_filter.graph,
                                f"gene.Notes!={final_call}"))

        query = querying.build_select(db_filter.graph, columns, 
                                                       where=conditionals)

        results = querying.execute(db_filter.engine, query, 
                                                in_column=db_filter.key,
                                                values=[data_dict["Pham"]])

        for result in results:
            if (not result["Accession"]) or (not result["LocusTag"]):
                continue
            result["Notes"] = data_dict["Final Call"]
            export_dicts.append(result)

    return export_dicts

#TODO Unittest
def use_csv_data(db_filter, data_dicts, columns, conditionals,
                                                               verbose=False):
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

    query = querying.build_select(alchemist.graph, columns,
                                                   where=conditionals)
    results = querying.execute(alchemist.engine, query,
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
            export_dicts.append(data_dict)

    return export_dicts

def write_curation_data(data_dicts, export_path, file_name=CURATION_NAME,
                                                 verbose=False):
    """Writes a curation submission csv.

    :param data_dicts: List of data dictionaries to convert to curation format.
    :type data_dicts: list[dict]
    :param export_path: Path to a dir for file creation.
    :type export_path: Path
    :param file_name: Name of the file to write curation data to.
    :type file_name: str
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    data_dicts = sorted(data_dicts, 
                              key=lambda data_dict: data_dict["PhageID"])

    for d in data_dicts:
        format_curation_data(d)

    if verbose:
        print(f"Writing {file_name} in {export_path.name}...")
    file_path = export_path.joinpath(file_name)
    fileio.export_data_dict(data_dicts, file_path, CURATION_HEADER, 
                                                    include_headers=True)

#-----------------------------------------------------------------------------
#CURATION-SPECIFIC HELPER FUNCTIONS

def format_curation_data(row_dict): 
    """Function to format revise dictionary keys.

    :param row_dict: Data dictionary for a revise file.
    :type row_dict: dict
    :param product: Gene product to append to the revise data dictionary.
    :type product: str
    """
    row_dict["Phage"] = row_dict.pop("PhageID")
    row_dict["Accession Number"] = row_dict.pop("Accession")
    row_dict["Locus Tag"] = row_dict.pop("LocusTag")
    row_dict["Product"] = row_dict.pop("Notes")

if __name__ == "__main__":
    args = sys.argv
    args.insert(0, "")
    main(args)
