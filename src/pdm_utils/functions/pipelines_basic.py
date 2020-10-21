import csv
import shutil
import sys
from functools import singledispatch
from pathlib import Path

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import basic


# PIPELINE OUTPUT HANDLING
# -----------------------------------------------------------------------------
def create_working_dir(working_path, dump=False, force=False):
    if not dump:
        if not force:
            if working_path.is_dir():
                print("COWARDLY ABORTING PIPELINE: "
                      f"Directory '{working_path}' already exists.")
                sys.exit(1)

    working_path.mkdir(parents=True, exist_ok=True)


def create_working_path(folder_path, folder_name, dump=False, force=False,
                        attempt=50):
    if folder_path is None:
        working_path = create_default_path(folder_name, force=force,
                                           attempt=attempt)
    else:
        working_path = folder_path.joinpath(folder_name)

    if dump:
        working_path = working_path.parent

    return working_path


def create_default_path(name, force=False, attempt=50):
    default_path = Path.cwd().joinpath(name)
    if not force:
        default_path = basic.make_new_dir(Path.cwd(), default_path,
                                          attempt=attempt, mkdir=False)

    return default_path


def convert_dir_path(path):
    """Function to convert argparse input to a working directory path.

    :param path: A string to be converted into a Path object.
    :type path: str
    :returns: A Path object converted from the inputed string.
    :rtype: Path
    """
    return basic.set_path(Path(path), kind="dir")


def convert_file_path(path):
    """Function to convert argparse input to a working file path.

    :param path: A string to be converted into a Path object.
    :type path: str
    :returns: A Path object converted from the inputed string.
    :rtype: Path
    """
    return basic.set_path(Path(path), kind="file")


# PIPELINE INPUT HANDLING
# -----------------------------------------------------------------------------
@singledispatch
def parse_value_input(value_list_input):
    """Function to convert values input to a recognized data types.

    :param value_list_input: Values stored in recognized data types.
    :type value_list_input: list[str]
    :type value_list_input: Path
    :returns: List of values to filter database results.
    :rtype: list[str]
    """

    print("Value list input is of an unexpected type.")
    sys.exit(1)


@parse_value_input.register(Path)
def _(value_list_input):
    value_list = []
    with open(value_list_input, newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for name in csv_reader:
            value_list.append(name[0])
    return value_list


@parse_value_input.register(list)
def _(value_list_input):
    return value_list_input


# MYSQL FILTERS AND HANDLER FUNCTIONS
def build_alchemist(database, ask_database=True, config=None, dialect="mysql"):
    if config is not None:
        username = config["mysql"].get("user")
        password = config["mysql"].get("password")
        if not (username is None or password is None):
            alchemist = AlchemyHandler(username=username, password=password,
                                       dialect=dialect)
            alchemist.connect(login_attempts=0, pipeline=True)

            alchemist.database = database
            alchemist.connect(ask_database=ask_database, pipeline=True)

            return alchemist

    alchemist = AlchemyHandler(database=database)
    alchemist.connect(ask_database=ask_database, pipeline=True)

    return alchemist


def build_filter(alchemist, key, filters, values=None, verbose=False):
    """Applies MySQL WHERE clause filters using a Filter.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param table: MySQL table name.
    :type table: str
    :param filters: A list of lists with filter values, grouped by ORs.
    :type filters: list[list[str]]
    :param groups: A list of supported MySQL column names.
    :type groups: list[str]
    :returns: filter-Loaded Filter object.
    :rtype: Filter
    """
    db_filter = Filter(alchemist=alchemist)
    db_filter.key = key
    db_filter.values = values

    if filters != "":
        try:
            db_filter.add(filters)
        except:
            print("Please check your syntax for the conditional string: "
                  f"{filters}")
            exit(1)

        db_filter.parenthesize()

    return db_filter


def add_sort_columns(db_filter, sort_columns, verbose=False):
    if verbose:
        print("Processing columns for sorting...")

    try:
        db_filter.sort(sort_columns)
    except:
        print("Please check your syntax for sorting columns:\n"
              f"{', '.join(sort_columns)}")
        sys.exit(1)


def build_groups_map(db_filter, export_path, groups=[], verbose=False,
                     force=False, dump=False):
    """Function that generates a map between conditionals and grouping paths.

    :param db_filter: A connected and fully loaded Filter object.
    :type db_filter: Filter
    :param export_path: Path to a dir for new dir creation.
    :type folder_path: Path
    :param groups: A list of supported MySQL column names.
    :type groups: list[str]
    :param conditionals_map: A mapping between group conditionals and Paths.
    :type conditionals_map: dict{Path:list}
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param previous: Value set by function to provide info for print statements
    :type previous: str
    :param depth: Value set by function to provide info for print statements.
    :type depth: int
    :returns conditionals_map: A mapping between group conditionals and Paths.
    :rtype: dict{Path:list}
    """
    conditionals_map = {}
    try:
        build_groups_tree(db_filter, export_path, conditionals_map,
                          groups=groups, verbose=False, force=False)
    except:        
        print("COWARDLY ABORTING PIPELINE: "
              "Found duplicate directories during path structuring.")
        sys.exit(1)

    return conditionals_map


def build_groups_tree(db_filter, export_path, conditionals_map, groups=[],
                      verbose=False, force=False, previous=None, depth=0):
    """Recursive function that generates directories based on groupings.

    :param db_filter: A connected and fully loaded Filter object.
    :type db_filter: Filter
    :param export_path: Path to a dir for new dir creation.
    :type folder_path: Path
    :param groups: A list of supported MySQL column names.
    :type groups: list[str]
    :param conditionals_map: A mapping between group conditionals and Paths.
    :type conditionals_map: dict{Path:list}\
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param previous: Value set by function to provide info for print statements
    :type previous: str
    :param depth: Value set by function to provide info for print statements.
    :type depth: int
    :returns conditionals_map: A mapping between group conditionals and Paths.
    :rtype: dict{Path:list}
    """
    groups = groups.copy()
    conditionals = db_filter.build_where_clauses()
    if not groups:
        conditionals_map.update({export_path: conditionals})
        return

    current_group = groups.pop(0)
    if verbose:
        if depth > 0:
            dots = ".." * depth
            print(f"{dots}Grouping by {current_group} in {previous}...")
        else:
            print(f"Grouping by {current_group}...")

    try:
        group_column = db_filter.get_column(current_group)
    except:
        print(f"Group '{current_group}' is not a valid group.")
        sys.exit(1)

    transposed_values = db_filter.build_values(column=group_column,
                                               where=conditionals)

    for group in transposed_values:
        group_path = export_path.joinpath(str(group))
        if group_path.is_dir():
            if force:
                shutil.rmtree(group_path)
            else:
                raise OSError("COWARDLY ABORTING PIPELINE: "
                              "Found pre-existing directories during "
                              "group path stucturing.")

        db_filter_copy = db_filter.copy()
        db_filter_copy.add(f"{current_group}={group}")

        previous = f"{current_group} {group}"
        build_groups_tree(db_filter_copy, group_path, conditionals_map,
                          groups=groups, verbose=verbose,
                          force=force, previous=previous, depth=depth+1)
