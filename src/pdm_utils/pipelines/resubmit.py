import argparse
import os
import sys
import time
from pathlib import Path

from sqlalchemy import select

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import basic
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb_basic
from pdm_utils.pipelines import export_db

#-----------------------------------------------------------------------------
#GLOBAL VARIABLES

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_resubmit"
DEFAULT_FOLDER_PATH = Path.cwd()
CSV_NAME = "resubmit.csv"

RESUBMIT_HEADER = ["Phage", "Accession Number", "Locus Tag", 
                   "Start", "Stop", "Product"]

RESUBMIT_COLUMNS = ["phage.PhageID", "phage.Accession", "gene.LocusTag",
                    "gene.Start", "gene.Stop"]

BASE_CONDITIONALS = ("phage.Status = final AND "
                     "phage.AnnotationAuthor = 1 AND"
                     "phage.RetrieveRecord = 1")

#-----------------------------------------------------------------------------
#MAIN FUNCTIONS

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the resubmit pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_resubmit(unparsed_args_list)
   
    alchemist = AlchemyHandler(database=args.database)
    alchemist.connect(ask_database=True, pipeline=True)

    revisions_data_dicts = basic.retrieve_data_dict(args.revisions_file)

    execute_resubmit(alchemist, revisions_data_dicts, args.folder_path, 
                                                      args.folder_name,
                                                      filters=args.filters,
                                                      groups=args.groups,
                                                      verbose=args.verbose)

def parse_resubmit(unparsed_args_list):
    """Parses resubmit arguments and stores them with an argparse object.

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
                                    type=export_db.convert_file_path)

    parser.add_argument("-o", "--folder_name", 
                                    type=str,  help=FOLDER_NAME_HELP)
    parser.add_argument("-p", "--folder_path", 
                                    type=export_db.convert_dir_path,
                                               help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true", 
                                               help=VERBOSE_HELP)

    parser.add_argument("-f", "--filter", nargs="?", dest="filters",
                                               help=FILTERS_HELP)
    parser.add_argument("-g", "--group_by", nargs="*",
                                help=GROUP_BY_HELP,
                                dest="groups")
    
    date = time.strftime("%Y%m%d")
    default_folder_name = f"{date}_pham_resubmit"
    default_folder_path = Path.cwd()

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH, 
                        filters="", groups=[], verbose=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args

def execute_resubmit(alchemist, revisions_data_dicts, folder_path, folder_name,
                                                     filters="", groups=[],
                                                     verbose=False):
    """Executes the entirety of the genbank resubmit pipeline.

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
    db_filter = Filter(alchemist=alchemist)
    db_filter.key = "gene.PhamID"
    db_filter.add(BASE_CONDITIONALS)

    if filters != "":
        try:
            db_filter.add(filters)
        except:
            print("Please check your syntax for the conditional string:\n"
                 f"{filters}")
    
    resubmit_columns = db_filter.get_columns(RESUBMIT_COLUMNS)
    
    phams = []
    for data_dict in revisions_data_dicts:
        phams.append(data_dict["Pham"])

    db_filter.values = phams
    
    if verbose:
        print("Creating export folder...")
    export_path = folder_path.joinpath(folder_name)
    export_path = basic.make_new_dir(folder_path, export_path, attempt=50)

    conditionals_map = {}
    export_db.build_groups_map(db_filter, export_path, conditionals_map,
                                                         groups=groups,
                                                         verbose=verbose)

    if verbose:
        print("Prepared query and path structure, beginning review export...")

    for mapped_path in conditionals_map.keys():
        if verbose:
            print("Retreiving phage data for pham revisions...")
        export_dicts = []
        for data_dict in revisions_data_dicts:
            if verbose:
                print(f"...Retrieving data for pham {data_dict['Pham']}...")

            conditionals = conditionals_map[mapped_path]

            final_call = data_dict["Final Call"]
            if final_call == "Hypothetical Protein":
                final_call = ""
            conditionals.append(querying.build_where_clause(alchemist.graph,
                                    f"gene.Notes!={final_call}"))

            query = querying.build_select(alchemist.graph, resubmit_columns, 
                                                           where=conditionals)

            results = querying.execute(alchemist.engine, query, 
                                                    in_column=db_filter.key,
                                                    values=[data_dict["Pham"]])

            for result in results:
                format_resubmit_data(result, data_dict["Final Call"]) 
                export_dicts.append(result)

        if not export_dicts:
            if verbose:
                print("'{mapped_path.name}' data selected for resubmision "
                      "matches selected call; no resubmision exported...")

            mapped_path.rmdir()
            continue

        export_dicts = sorted(export_dicts, 
                              key=lambda export_dict: export_dict["Phage"])


        if verbose:
            print(f"Writing {CSV_NAME} in {mapped_path.name}...")
        file_path = mapped_path.joinpath(CSV_NAME)
        basic.export_data_dict(export_dicts, file_path, RESUBMIT_HEADER, 
                                                        include_headers=True)

#-----------------------------------------------------------------------------
#RESUBMIT-SPECIFIC HELPER FUNCTIONS

def format_resubmit_data(row_dict, product): 
    """Function to format resubmit dictionary keys.

    :param row_dict: Data dictionary for a resubmit file.
    :type row_dict: dict
    :param product: Gene product to append to the resubmit data dictionary.
    :type product: str
    """
    row_dict["Phage"] = row_dict.pop("PhageID")
    row_dict["Accession Number"] = row_dict.pop("Accession")
    row_dict["Locus Tag"] = row_dict.pop("LocusTag")
    row_dict["Start"] = row_dict.pop("Start")
    row_dict["Stop"] = row_dict.pop("Stop")
    row_dict["Product"] = product

def _sort_data(data_dict):
    return data_dict["Phage"]

if __name__ == "__main__":
    args = sys.argv
    args.insert(0, "")
    main(args)
