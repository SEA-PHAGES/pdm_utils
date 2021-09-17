"""Pipeline for mapping the differences of PhamIDs between dataabse"""
import argparse
import shutil
import sys
import time
from pathlib import Path

from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio
from pdm_utils.functions import pipelines_basic

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pham_finder"
DEFAULT_FOLDER_PATH = Path.cwd()

PHAM_FINDER_HEADER = ["Reference Pham", "Corresponding Phams"]


def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the pham_finder pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_pham_finder(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(None, ask_database=False,
                                                config=config)
    
    values = None
    if args.input:
        values = pipelines_basic.parse_value_input(args.input)

    execute_pham_finder(alchemist, args.folder_path, args.folder_name,
                        args.adatabase, args.bdatabase, values=values,
                        filters=args.filters, groups=args.groups, 
                        sort=args.sort, show_per=args.show_percentages,
                        use_locus=args.use_locus, verbose=args.verbose)


def parse_pham_finder(unparsed_args_list):
    """Parses pham_finder arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    A_DATABASE_HELP = """
        Name of the MySQL database to retrieve reference phams from.
        """
    B_DATABASE_HELP = """
        Name of the MySQL database to find phams for.
        """

    CONFIG_FILE_HELP = """
        Find option that enables use of a config file for sourcing credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
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
    SINGLE_PHAMIDS_HELP = """
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

    SHOW_PERCENTAGES_HELP = """
        Pham finder option that enables the display of the percent of 
        the referenced pham that the found pham covers.
        """
    USE_LOCUS_HELP = """
        Pham finder option that converts between phams using LocusTags.
        """

    parser = argparse.ArgumentParser()
    parser.add_argument("adatabase", type=str, 
                        help=A_DATABASE_HELP)
    parser.add_argument("bdatabase", type=str, 
                        help=B_DATABASE_HELP)

    parser.add_argument("-c", "--config_file", 
                        type=pipelines_basic.convert_file_path,
                        help=CONFIG_FILE_HELP)
    parser.add_argument("-m", "--folder_name", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path", 
                        type=pipelines_basic.convert_dir_path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)

    parser.add_argument("-if", "--import_file", dest="input",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP)
    parser.add_argument("-in", "--import_names", dest="input", nargs="*",
                        help=SINGLE_PHAMIDS_HELP)
    parser.add_argument("-f", "--where", nargs="?", dest="filters",
                        help=WHERE_HELP)
    parser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                        help=GROUP_BY_HELP)
    parser.add_argument("-s", "--order_by", nargs="*", dest="sort",
                        help=ORDER_BY_HELP)

    parser.add_argument("-per", "--show_percentages", action="store_true", 
                        help=SHOW_PERCENTAGES_HELP)
    parser.add_argument("-ul", "--use_locus", action="store_true",
                        help=USE_LOCUS_HELP)

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH,
                        config_file=None, verbose=False, input=[],
                        filters="", groups=[], sort=[],
                        show_percentages=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


# TODO Owen Needs unittests
def execute_pham_finder(alchemist, folder_path, folder_name,
                        adatabase, bdatabase, values=None,
                        filters="", groups=[], sort=[],
                        show_per=False, use_locus=False, verbose=False):
    """Executes the entirety of the file export pipeline.

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param adatabase: Name of reference database to source phams from.
    :type adatabase: str
    :param bdatabase: Name of database to find corresponding phams for.
    :type bdatabase: str
    :param values: List of values to filter database results:
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param table: MySQL table name.
    :type table: str
    :param filters: A list of lists with filter values, grouped by ORs.
    :type filters: str
    :param groups: A list of supported MySQL column names to group by.
    :type groups: A list of supported MySQL column names to group by.
    :type groups: list[str]
    :param sort: A list of supported MySQL column names to sort by.
    :type sort: list[str]
    :param show_per: Enables display gene coverage of the corresponding phams.
    :type show_per: bool
    :param use_locus: Toggles conversion between phams using LocusTag instead
    :type use_locus: bool
    """
    if not (adatabase in alchemist.databases and
            bdatabase in alchemist.databases):
        print("User credentials does not have access to both "
              f"databases {adatabase} and {bdatabase}.\n"
              "Please check your database access and try again.")
        sys.exit(1)

    alchemist.database = adatabase
    alchemist.connect()
    a_filter = pipelines_basic.build_filter(alchemist, "gene.PhamID", filters,
                                            values=values, verbose=verbose)

    alchemist.database = bdatabase
    alchemist.connect()
    if use_locus:
        b_filter = pipelines_basic.build_filter(alchemist, "gene.LocusTag", "")
    else:
        b_filter = pipelines_basic.build_filter(alchemist, "gene", "")

    if sort:
        try:
            a_filter.sort(sort)
        except:
            print("Please check your syntax for sorting columns:\n"
                  f"{', '.join(sort)}")
            sys.exit(1)

    if verbose:
        print("Creating pham_finder folder...")
    export_path = folder_path.joinpath(folder_name)
    export_path = basic.make_new_dir(folder_path, export_path, attempt=50)

    conditionals_map = pipelines_basic.build_groups_map(a_filter, export_path,
                                                        groups=groups,
                                                        verbose=verbose)

    if verbose:
        print("Prepared query and path structure, beginning export...")

    values = a_filter.values
    for mapped_path in conditionals_map.keys():
        a_filter.reset()
        a_filter.values = values

        conditionals = conditionals_map[mapped_path]
        a_filter.values = a_filter.build_values(where=conditionals)

        if a_filter.hits() == 0:
            print("No database entries received from gene.PhamID "
                  f"for '{mapped_path}'.")
            shutil.rmtree(mapped_path)
            continue

        if sort:
            a_filter.sort(sort)

        mapped_phams = find_phams(a_filter, b_filter, show_per=show_per)
        if not mapped_phams:
            print("Phams are consistent between the two databases "
                  f"for '{mapped_path}'.")
            shutil.rmtree(mapped_path)
            continue

        out_data_dicts = []
        for ref_pham, corr_phams in mapped_phams.items():
            data_dict = {}
            data_dict[PHAM_FINDER_HEADER[0]] = ref_pham
            data_dict[PHAM_FINDER_HEADER[1]] = corr_phams
            out_data_dicts.append(data_dict)

        file_path = mapped_path.joinpath("PhamMap.csv")
        fileio.export_data_dict(out_data_dicts, file_path, PHAM_FINDER_HEADER,
                                include_headers=True)


# TODO Owen Needs unittests
def find_phams(a_filter, b_filter, show_per=False, use_locus=False):
    """Find phams helper function that finds phams via GeneID intermediates.

    :param a_filter: Fully built Filter connected to the reference database.
    :type a_filter: Filter
    :param b_filter: Fully build Filter connected to a database.
    :type b_filter: Filter
    :param show_per: Enables display gene coverage of the corresponding phams.
    :type show_per: bool
    :returns: Returns a dictionary mapping original phams to corresponding phams
    :rtype: dict{int:str}
    """
    if use_locus:
        genes = a_filter.retrieve("gene.GeneID")
    else:
        genes = a_filter.retrieve("gene.GeneID")

    mapped_phams = {}
    for pham, data_dict in genes.items():
        if use_locus:
            values = data_dict["LocusTag"]
            values.remove(None)
            b_filter.values = values
        else:
            b_filter.values = data_dict["GeneID"]        

        if show_per:
            pham_groups = b_filter.group("gene.PhamID")
            phams_list = list(pham_groups.keys())
            total_genes = 0
            
            for grouped_genes in pham_groups.values():
                total_genes += len(grouped_genes)
        else:
            phams_list = b_filter.transpose("gene.PhamID")

        if len(phams_list) == 1:
            if phams_list[0] == pham:
                print(f"{phams_list[0]} == {pham}")
                continue
            else:
                corr_phams = str(phams_list[0])
        elif len(phams_list) == 0:
            corr_phams = "None"
        else:
            if show_per:
                for i in range(len(phams_list)):
                    join_pham = phams_list[i]
                    percent = (len(pham_groups[join_pham])\
                                  /total_genes) * 100
                    percent = round(percent, 1)
                    phams_list[i] = "".join([str(join_pham), 
                                        "(", str(percent), "%)"])
            corr_phams = ";".join([str(join_pham) for join_pham in phams_list])

        mapped_phams[pham] = corr_phams
        
    return mapped_phams


if __name__ == "__main__":
    main(sys.argv)
