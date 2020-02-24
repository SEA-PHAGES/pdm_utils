from datetime import datetime
from decimal import Decimal
import re

#Global file constants
COMPARATIVE_OPERATORS = [">", ">=", "<", "<="]
OPERATORS             = ["=", "!="] + COMPARATIVE_OPERATORS
COMPARABLE_TYPES      = [int, Decimal, float, datetime]
TYPES                 = [str, bytes] + COMPARABLE_TYPES
GROUP_OPTIONS = ["limited_set", "num_set", "str_set"]

def parse_column(unparsed_column):
    """Helper function to return a two-dimensional array of group parameters.

    :param unparsed_groups:
        Input a list of group expressions to parse and split.
    :type unparsed_groups: List[str]
    :return groups:
        Returns a two-dimensional array of group parameters.
    :type groups: List[List[str]]
    """
    column_format = re.compile("\w+\.\w+", re.IGNORECASE)

    if re.match(column_format, unparsed_column) != None:
        column = re.split("\W+", unparsed_column)
    else:
        raise ValueError(f"Unsupported table/column format: "
                         f"'{unparsed_column}'")

    return column

def parse_filter(unparsed_filter):
    """Helper function to return a two-dimensional array of filter parameters.

    :param unparsed_filters:
        Input a list of filter expressions to parse and split.
    :type unparsed_filters: List[str]
    :return filters:
        Returns a two-dimensional array of filter parameters.
    :type filters: List[List[str]]
    """
    filter_format = re.compile("\w+\.\w+ *([=<>!]+| +LIKE +| +IS NOT +) *\w+")

    if re.match(filter_format, unparsed_filter) != None:
        operators = re.compile("([=<>!]+| +LIKE +| +IS NOT +)")

        operator_split = re.split(operators, unparsed_filter)
        column_split = re.split("\W+", operator_split[0])
        
        filter = column_split + operator_split[1:]
    else:
        raise ValueError(f"Unsupported filtering format: '{unparsed_filter}'")
                
    return filter

def parse_cmd_filters(unparsed_cmd_string):
    cmd_filters_format = re.compile(
                "(((\w+\.\w+ *[=<>!( +LIKE +)( +IS NOT +)]+ *\w+)( +AND +))|"
                "((\w+\.\w+ *[=<>!(LIKE)(IS NOT)]+ *\w+)( +OR +))|"
                "( *\w+\.\w+ *[=<>!(LIKE)(IS NOT)]+ *\w+))+")
    
    and_splits = []
    or_splits = []
    if re.match(cmd_filters_format, unparsed_cmd_string) != None: 
        and_splits = re.split(" +AND +", unparsed_cmd_string)
        for split in and_splits:
            or_splits.append(re.split(" +OR +", split))


    return or_splits

def check_operator(operator, column_object):
        """Parses a operator string to match a MySQL query operators.

        :param operator:
            Input a raw operator string for an accepted MySQL operator.
        :type operator: str
        :param table:
            Input a case-sensitive string for a TableNode id.
        :type table: str
        :param field:
            Input a case-sensitive string for a ColumnNode id.
        :type field: str
        :param verbose:
            Set a boolean to control the terminal output.
        :type verbose: Boolean
        """
        if operator not in OPERATORS:
            raise ValueError(f"Operator {operator} is not supported.")

        column_type = column_object.type.python_type

        if column_type not in TYPES:
            raise ValueError(f"Column '{column_object.name}' "
                             f"has an unsupported type, {column_type}.")
        if operator in COMPARATIVE_OPERATORS and \
           column_type not in COMPARABLE_TYPES:
            raise ValueError(f"Column '{column_object.name}' "
                             f"is not comparable with '{operator}'.")

def translate_table(metadata, raw_table):
    for table in metadata.tables.keys():
        if table.lower() == raw_table.lower():
            return table

    raise ValueError(f"Table '{raw_table}' requested to be filtered "
                     f"is not in selected database.")

def translate_column(metadata, raw_column):
    parsed_column = parse_column(raw_column)

    table = translate_table(metadata, parsed_column[0])
    table_obj = metadata.tables[table]
    for column in table_obj.columns.keys():
        if column.lower() == parsed_column[1].lower():
            return column

    raise ValueError(f"Field '{parsed_column[1]}' requested to be filtered"
                     f" is not in '{table_obj.name}'")

