import argparse
import sys
import time
from pathlib import Path

from sqlalchemy import Column
from sqlalchemy import and_
from sqlalchemy.sql import func

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import parsing
from pdm_utils.functions import querying
from pdm_utils.functions import basic
from pdm_utils.pipelines import export_db

#-----------------------------------------------------------------------------
#GLOBAL VARIABLES

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_review"
DEFAULT_FOLDER_PATH = Path.cwd()


PF_HEADER = ["Pham", 
             "Final Call",
             "#Members", 
             "Clusters", 
             "#Functions", 
             "Functional Calls"]
PG_HEADER = ["Gene",
             "Phage",
             "Gene#",
             "Cluster",
             "Subcluster",
             "Functional Call",
             "Translation"]

PF_DATA_COLUMNS = ["gene.GeneID", "phage.Cluster", "gene.Notes"]
PG_DATA_COLUMNS = ["phage.PhageID", "gene.Name", "phage.Cluster", 
                   "phage.Subcluster", "gene.Notes", "gene.Translation"]

BASE_CONDITIONALS = ("phage.Status = final AND "
                     "phage.AnnotationAuthor = 1 AND "
                     "phage.RetrieveRecord = 1")

#-----------------------------------------------------------------------------
#MAIN FUNCTIONS

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the review pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_review(unparsed_args_list)

    alchemist = AlchemyHandler(database=args.database)
    alchemist.connect(ask_database=True, pipeline=True)

    values = export_db.parse_value_input(args.input)
    
    execute_review(alchemist, args.folder_path, args.folder_name,
                   review=args.review, values=values,
                   filters=args.filters, groups=args.groups, sort=args.sort,
                   g_reports=args.gene_reports, s_report=args.summary_report,
                   verbose=args.verbose)

def parse_review(unparsed_args_list):
    """Parses review arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    DATABASE_HELP = """
        Name of the MySQL database to export from.
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

    GENE_REPORTS_HELP = """
        Export option to toggle export of supplemental information about 
        the genes in the phams selected.
        """
    SUMMARY_REPORT_HELP = """
        Export option to toggle export of supplemental information about
        the profile of the phams selected.
        """

    REVIEW_HELP = """
        Review option to toggle review of phams.  If enabled,
        phams inputted or auto-generated will not be reviewed
        for inconsistencies.
        """
    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    IMPORT_NAMES_HELP = """
        Selection input option that imports values from cmd line input.
            Follow selection argument with space separated
            names of genomes in the database.
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
    
    parser = argparse.ArgumentParser()

    parser.add_argument("database", type=str,  help=DATABASE_HELP)
    parser.add_argument("-o", "--folder_name", 
                                    type=str,  help=FOLDER_NAME_HELP)
    parser.add_argument("-p", "--folder_path", 
                                    type=export_db.convert_dir_path,
                                               help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true", 
                                               help=VERBOSE_HELP)

    parser.add_argument("-gr", "--gene_reports", action="store_true",
                                               help=GENE_REPORTS_HELP)
    parser.add_argument("-sr", "--summary_report", action="store_true",
                                               help=SUMMARY_REPORT_HELP)

    parser.add_argument("-r", "--review", action="store_false",
                                               help=REVIEW_HELP)
    parser.add_argument("-if", "--import_files", dest="input",
                                    type=export_db.convert_file_path,
                                               help=IMPORT_FILE_HELP)
    parser.add_argument("-in", "--import_names", nargs="*", dest="input",
                                               help=IMPORT_NAMES_HELP)

    parser.add_argument("-f", "--where", nargs="?", dest="filters",
                                               help=WHERE_HELP)
    parser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                                               help=GROUP_BY_HELP)
    parser.add_argument("-s", "--order_by", nargs="*", dest="sort",
                                               help=ORDER_BY_HELP)
    
    
    default_folder_name = DEFAULT_FOLDER_NAME
    default_folder_path = DEFAULT_FOLDER_PATH

    parser.set_defaults(folder_name=default_folder_name,
                        folder_path=default_folder_path,
                        input=[], filters="", groups=[], sort=[],
                        review=True, gene_report=False, summary_report=False,
                        verbose=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args

def execute_review(alchemist, folder_path, folder_name, 
                              review=True, values=[],
                              filters="", groups=[], sort=[],
                              g_reports=False, s_report=False,
                              verbose=False):
    """Executes the entirety of the pham review pipeline.
    
    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param csv_title: Title for an appended csv file prefix.
    :type csv_title: str
    :param review: A boolean to toggle filtering of phams by pham discrepancies.
    :type review: bool
    :param values: List of values to filter database results.
    :type values: list[str]
    :param filters: A list of lists with filter values, grouped by ORs.
    :type filters: list[list[str]]
    :param groups: A list of supported MySQL column names to group by.
    :type groups: list[str]
    :param sort: A list of supported MySQL column names to sort by. 
    :param g_reports: A boolean to toggle export of additional pham information.
    :type g_reports: bool
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    db_filter = Filter(alchemist=alchemist)
    db_filter.key = ("gene.PhamID")
 
    if values:
        db_filter.values = values

    if verbose:
        print(f"Identified {len(values)} phams to review...")
           
    if filters != "":
        try:
            db_filter.add(filters)
        except:
            print("Please check your syntax for the conditional string:\n"
                 f"{filters}")
            sys.exit(1)
        finally:
            db_filter.update() 

    db_filter._filters = []
    db_filter._updated = False 
    db_filter._or_index = -1

    db_filter.add(BASE_CONDITIONALS)
    db_filter.update()

    if not db_filter.values:
        print("Current settings produced no database hits.")
        sys.exit(1)

    if review: 
        review_phams(db_filter, verbose=verbose)

    if sort:
        db_filter.sort(sort)

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
    original_phams = db_filter.values
    total_g_data = {}
    for mapped_path in conditionals_map.keys():
        conditionals = conditionals_map[mapped_path]
        db_filter.values = original_phams
        db_filter.values = db_filter.build_values(where=conditionals)

        pf_data = get_pf_data(alchemist, db_filter, verbose=verbose) 
        write_report(pf_data, mapped_path, PF_HEADER,
                     csv_name=f"FunctionReport",
                     verbose=verbose)

        if g_reports:
            execute_g_report_export(alchemist, db_filter, mapped_path, 
                                                    total_g_data=total_g_data,
                                                    verbose=verbose)

        if s_report:
            execute_s_report_export(alchemist, db_filter, conditionals, 
                                                    mapped_path,
                                                    verbose=verbose)
               
def execute_g_report_export(alchemist, db_filter, export_path, total_g_data={},
                                                         verbose=False):
    """Executes export of gene data for a reviewed pham.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param export_path: Path to a valid dir for new file creation.
    :type export_path: Path
    :param total_g_data: Total data extracted for gene reports.
    :type total_g_data: dict
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    phams = db_filter.values

    pg_path = export_path.joinpath("GeneReports")
    pg_path.mkdir()

    for pham in phams:
        db_filter.values = [pham]

        try:
            g_data = total_g_data[pham]
        except:
            g_data = get_g_data(alchemist, db_filter, pham, 
                                                        verbose=verbose)
            total_g_data[pham] = g_data

        write_report(g_data, pg_path, PG_HEADER,
                     csv_name=f"{pham}_GeneReport",
                     verbose=verbose)

    db_filter.values = phams

def execute_s_report_export(alchemist, db_filter, conditionals, export_path, 
                                                    verbose=False):
    """Executes export of summary data for all reviewed phams.

    :param alchemist: A connected and fully build
    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param export_path: A path to a valid dir for new file creation.
    :type export_path: Path
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    phams = db_filter.values

    if verbose:
        print(f"Retrieving SummaryReport data...")
    phages_histogram = {}
 
    for pham in phams:
        db_filter.values = [pham]
        phages = db_filter.build_values(column="phage.PhageID",
                                        where=conditionals)

        increment_histogram(phages, phages_histogram)

    db_filter.values = phams
    db_filter.transpose("phage.PhageID", set_values=True)
    db_filter.sort("phage.DateLastModified")

    submitted_sorted_phages = db_filter.values
    submitted_sorted_phages.reverse()
    last_submitted_phages = chunk_list(submitted_sorted_phages, 5)[0]

    occurance_sorted_phages = sort_histogram_keys(phages_histogram)
    top_occuring_phages = chunk_list(occurance_sorted_phages, 5)[0]

    version_data = mysqldb_basic.get_first_row_data(alchemist.engine, "version")

    if verbose:
        print(f"Writing SummaryReport.txt in {export_path.name}")

    s_path = export_path.joinpath("SummaryReport.txt")
    s_file = open(s_path, "w")

    s_file.write(f"Phams reviewed on: {time.strftime('%d-%m-%Y')}\n")
    s_file.write(f"Database reviewed: {alchemist.database}\n")
    s_file.write(f"Schema version: {version_data['SchemaVersion']} "
                 f"Database version: {version_data['Version']}\n\n") 
    s_file.write(f"Phams reviewed using the following base conditionals:\n")
    s_file.write(f"    {BASE_CONDITIONALS}\n")

    s_file.write(f"\n\n")
    s_file.write(f"Most occuring phages: {', '.join(top_occuring_phages)}\n")
    s_file.write(f"Phages recently submitted: {', '.join(last_submitted_phages)}\n")
    s_file.close()

    db_filter.values = phams
    db_filter.key = "pham.PhamID"

def review_phams(db_filter, verbose=False):
    """Finds and stores phams with discrepant function calls in a Filter.

    """
    notes = db_filter.get_column("gene.Notes")
    
    if verbose:
        print("Reviewing phams...")
    
    reviewed_phams = []
    for index in range(len(db_filter.values)):
        pham = db_filter.values[index]
        if verbose:
            print(f"...Analyzing Pham {pham}...")


        query = querying.build_count(db_filter.graph, notes.distinct(),
                                                    where=(db_filter.key==pham))
        func_count = mysqldb_basic.scalar(db_filter.engine, query)

        if func_count <= 1:
           continue              
        
        if verbose:
            print(f"......Detected discrepencies in Pham {pham}") 
        reviewed_phams.append(pham)

    if verbose:
        print(f"Detected {len(reviewed_phams)} disrepent phams...")

    db_filter.values = reviewed_phams

def write_report(data, export_path, header, csv_name="PhamReport",
                                            verbose=False):
    """Outputs a csv file

    """
    if not export_path.is_dir():
        print("Passed in path is not a directory.")
        sys.exit(1)

    file_path = export_path.joinpath(f"{csv_name}.csv")
    if verbose:
        print(f"Writing {file_path.name} in {export_path.name}...")

    basic.export_data_dict(data, file_path, header, include_headers=True)

#-----------------------------------------------------------------------------
#REVIEW-SPECIFIC HELPER FUNCTIONS

def get_pf_data(alchemist, db_filter, verbose=False):
    """
    """
    if verbose:
        print("Retrieving data for phams...")
     
    pf_columns = get_pf_data_columns(alchemist) 
    row_dicts = db_filter.retrieve(pf_columns)

    pf_data = []
    for pham in row_dicts.keys():
        if verbose:
            print(f"...Processing data for pham {pham}...")
        row_dict = row_dicts[pham]

        format_pf_data(row_dict, pham)
        pf_data.append(row_dict)

    return pf_data

def get_g_data(alchemist, db_filter, pham, verbose=False):  
    if verbose:
        print(f"Retrieving genes in pham {pham}...")
    db_filter.transpose("gene.GeneID", set_values=True) 
   
    pg_columns = get_g_data_columns(alchemist) 
    row_dicts = db_filter.retrieve(pg_columns)

    g_data = []
    for gene in row_dicts.keys():
        if verbose:
            print(f"...Processing data for gene {gene}...")

        row_dict = row_dicts[gene]

        format_g_data(row_dict, gene)
        g_data.append(row_dict)

    db_filter.key = "gene.PhamID"
    return g_data 

def format_pf_data(row_dict, pham):
    """Function to format function report dictionary keys.
    
    :param row_dict: Data dictionary for a function report.
    :type row_dict: dict
    :param pham: PhamID to append to the function report data dictionary.
    :type pham: int
    """
    row_dict["Pham"] = pham
    row_dict["Final Call"] = ""
    row_dict["#Members"] = len(row_dict.pop("GeneID"))
    row_dict["#Functions"] = len(row_dict["Notes"])
 
    row_dict["Clusters"] = ";".join([str(cluster) \
                                    for cluster in row_dict.pop("Cluster")])

    notes = row_dict.pop("Notes")
    for index in range(len(notes)):
        if notes[index] is None or notes[index] == "": 
            notes[index] = "Hypothetical Protein"

    row_dict["Functional Calls"] = ";".join(notes)

def format_g_data(row_dict, gene):
    """Function to format gene report dictionary keys.

    :param row_dict: Data dictionary for a gene report.
    :type row_dict: dict
    :param gene: GeneID to append to the gene report data dictionary.
    :type gene: str
    """
    row_dict["Gene"] = gene
    row_dict["Phage"] = row_dict.pop("PhageID")[0]
    row_dict["Gene#"] = row_dict.pop("Name")[0]
    row_dict["Cluster"] = row_dict["Cluster"][0]
    row_dict["Subcluster"] = row_dict["Subcluster"][0]
    row_dict["Translation"] = row_dict["Translation"][0]

    note = row_dict.pop("Notes")[0]
    if note is None or note == "":
        note = "Hypothetical Protein"

    row_dict["Functional Call"] = note

def get_pf_data_columns(alchemist):
    """Gets labelled columns for pham function data retrieval.
    
    :returns: List of labelled columns for function data retrieval.
    :rtype: list[Column]
    """
    pf_columns = []

    for column_name in PF_DATA_COLUMNS:
        pf_columns.append(querying.get_column(alchemist.metadata, column_name))

    return pf_columns

def get_g_data_columns(alchemist):
    """Gets labelled columns for pham gene data retrieval.

    :returns: List of labelled columns for gene data retrieval.
    :rtype: list[Column]
    """
    pg_columns = []

    for column_name in PG_DATA_COLUMNS:
        pg_columns.append(querying.get_column(alchemist.metadata, column_name))

    return pg_columns 

def increment_histogram(data, histogram):
    """Increments a dictionary histogram based on given data.

    :param data: Data to be used to index or create new keys in the histogram.
    :type data: list
    :param histogram: Dictionary containing keys whose values contain counts.
    :type histogram: dict
    """
    for item in data:
        try:
            histogram[item] += 1
        except:
            histogram[item] = 1

def sort_histogram_keys(histogram):
    """Sorts a dictionary histogram by its values and returns the sorted keys.
    
    :param histogram: Dictionary containing keys whose values contain counts.
    :type histogram: dict
    :returns: A list containing the keys from the histogram sorted by value.
    :rtype: list
    """
    sorted_keys = [key for key, value in sorted(histogram.items(), 
                                                key=lambda item:item[1], 
                                                reverse=True)] 

    return sorted_keys

def chunk_list(data_list, size):
    """Chunks list into a list of lists with the given size.

    :param data_list: List to be split into equal-sized lists.
    :type data_list: list
    :param size: Length of the resulting list chunks.
    :param size: int
    :returns: Returns list of lists with length of the given size.
    :rtype: list[list]
    """
    chunked_list = [data_list[i*size:(i+1)*size]\
            for i in range((len(data_list) + size - 1) // size)]

    return chunked_list


if __name__ == "__main__":
    args = sys.argv
    args.insert(0, "")
    main(args)
