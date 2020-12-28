from datetime import datetime
from decimal import Decimal
import re

from sqlalchemy import Column
#----------------------------------------------------------------------------
#GLOBAL VARIABLES
VALUE_CHARACTERS      = "\w\d\-_%/(),'"
VC = VALUE_CHARACTERS
VALUE_EXPRESSION      = f"'[ {VC}]*'|[{VC}]+"
VE = VALUE_EXPRESSION
FILTER_FORMAT         = (" *\w+\.\w+ *([=<>!]+ *| +LIKE +| +IS NOT +){1} *"
                        f"({VE}) *")

IN_FORMAT             = (" *\w+\.\w+ *( +IN +| +NOT IN +)"
                        f"\(([\s\w\W\d]+)\)")



NUMERIC_OPERATORS     = [">", ">=", "<", "<="]
NONNUMERIC_OPERATORS  = ["=", "!=", "IS NOT", "LIKE", "IN", "NOT IN"]
OPERATORS             = NUMERIC_OPERATORS + NONNUMERIC_OPERATORS
COMPARABLE_TYPES      = [int, Decimal, float, datetime]
TYPES                 = [str, bytes] + COMPARABLE_TYPES
GROUP_OPTIONS         = ["limited_set", "num_set", "str_set"]

#----------------------------------------------------------------------------
#SPACE HANDLING

def parse_out_ends(unparsed_string):
    """Parse and remove beginning and end whitespace of a string.

    :param unparsed_string: String with variable terminal whitespaces.
    :type unparsed_string: str
    :returns: String with parsed and removed beginning and ending whitespace.
    :rtype: str
    """
    if not isinstance(unparsed_string, str):
        raise TypeError(f"Parameter type required is of type list.")

    beginning_trimmed = unparsed_string.lstrip()
    trimmed_string = beginning_trimmed.rstrip()

    return trimmed_string

def parse_in_spaces(unparsed_string_list):
    """Convert a list of strings to a single space separated string.

    :param unparsed_string_list: String list to be concatenated
    :type unparsed_string_list: list[str]
    :returns: String with parsed in whitespace.
    :rtype: str
    """
    if not isinstance(unparsed_string_list, list):
        raise TypeError(f"Parameter type required is of type list.")

    joined_string = unparsed_string_list[0]
    unparsed_string_list = unparsed_string_list[1:]

    for string in unparsed_string_list:
        joined_string = (f"{joined_string} {string}")

    return joined_string

def parse_out_spaces(unparsed_string):
    """Parse and remove beginning and internal white space of a string.

    :param unparsed_string: String with variable terminal whitespaces.
    :type unparsed_string: str
    :returns: String with parsed and removed beginning and ending whitespace.
    :rtype: str
    """
    if not isinstance(unparsed_string, str):
        raise TypeError(f"Paremeter type required is of type list.")

    trimmed_string = parse_out_ends(unparsed_string)
    split_string = unparsed_string.split(" ")
    joined_string = "".join(split_string)

    return joined_string

#----------------------------------------------------------------------------
#MYSQL TABLE/COLUMN HANDLING

def parse_column(unparsed_column):
    """Recognizes and parses a MySQL structured column.

    :param unparsed_column: Formatted MySQL column.
    :type unparsed_column: str
    :returns: List containing segments of a MySQL column.
    :rtype: list[str]
    """
    column_format = re.compile(" *\w+\.\w+ *", re.IGNORECASE)

    if re.match(column_format, unparsed_column) != None:
        trimmed_column = parse_out_ends(unparsed_column)
        parsed_column = re.split("\.", trimmed_column)

        parsed_column[0] = parse_out_ends(parsed_column[0])
        parsed_column[1] = parse_out_ends(parsed_column[1])
    else:
        raise ValueError(f"Unsupported table/column format: "
                         f"'{unparsed_column}'")

    return parsed_column

def translate_table(metadata, raw_table):
    """Converts a case-insensitive table name to a case-sensitive str.

    :param metadata: Reflected SQLAlchemy MetaData object.
    :type metadata: MetaData
    :param raw_table: Case-insensitive table name.
    :type_table: str
    :returns: Case-sensitive table name.
    :rtype: str
    """
    for table in metadata.tables.keys():
        if table.lower() == raw_table.lower():
            return table

    raise ValueError(f"Table '{raw_table}' requested to be filtered "
                     f"is not in selected database.")

def translate_column(metadata, raw_column):
    """Converts a case-insensitve {table}.{column} str to a case-sensitive str.

    :param metadata: Reflected SQLAlchemy MetaData object.
    :type metadata: MetaData
    :param raw_column: Case-insensitive {table}.{column}.
    :type raw_column: str
    :returns: Case-sensitive column name.
    :rtype: str
    """
    parsed_column = parse_column(raw_column)

    table = translate_table(metadata, parsed_column[0])
    table_obj = metadata.tables[table]
    for column in table_obj.columns.keys():
        if column.lower() == parsed_column[1].lower():
            return column

    raise ValueError(f"Field '{parsed_column[1]}' requested to be filtered"
                     f" is not in '{table_obj.name}'")

#----------------------------------------------------------------------------
#MYSQL CONDITIONAL STRING HANDLING

def parse_filter(unparsed_filter):
    """Recognizes and parses a MySQL structured WHERE clause.

    :param unparsed_filter: Formatted MySQL WHERE clause.
    :type unparsed_filters: str
    :returns: List containing segments of a MySQL WHERE clause.
    :rtype: list[str]
    """

    filter_format = re.compile(FILTER_FORMAT) 
    in_format = re.compile(IN_FORMAT) 

    value_format = re.compile(f"({VE})")
    quote_value_format = re.compile(f"('[ {VC}]*'|'[ {VC}]+'s')")
    in_delimiter = re.compile("' *,")
    whitespace = re.compile(" +")

    if not re.match(filter_format, unparsed_filter) is None:
        operators = re.compile("( *[=<>!]+ *| +LIKE +| +IS NOT +)")

        operator_split = re.split(operators, unparsed_filter)
        parsed_column = parse_column(operator_split[0])
       
        operator_split[1] = parse_out_ends(operator_split[1])
        operator_split[2] = parse_out_ends(operator_split[2])

        if re.match(quote_value_format, operator_split[2]) != None:
            operator_split[2] = operator_split[2].rstrip("'")
            operator_split[2] = operator_split[2].lstrip("'")

        parsed_filter = parsed_column + operator_split[1:] 

    elif not re.match(in_format, unparsed_filter) is None:
        operators = re.compile("( +IN +| +NOT IN +)")

        operator_split = re.split(operators, unparsed_filter)
        parsed_column = parse_column(operator_split[0])

        operator_split[1] = parse_out_ends(operator_split[1]) 
        operator_split[2] = parse_out_ends(operator_split[2])  
 
        parenthesized = operator_split[2]
        parenthesized = parenthesized.lstrip("(")
        parenthesized = parenthesized.rstrip(")")

        in_values = []
        if not re.search(in_delimiter, parenthesized) is None:
            parenthesized = re.split(in_delimiter, parenthesized)     
            for value in parenthesized:
                value = parse_out_ends(value)
                value = value.lstrip("'")
                value = value.rstrip("'")
                in_values.append(parse_out_ends(value))
        else:
            parenthesized = parenthesized.split(",")
            for value in parenthesized:
                value = parse_out_ends(value)
                if value == "" or not re.match(whitespace, value) is None:
                    continue
                in_values.append(value)

        operator_split[2] = in_values
        parsed_filter = parsed_column + operator_split[1:]

    else:
        raise ValueError(f"Unsupported filtering format: '{unparsed_filter}'")
     
    return parsed_filter

def create_filter_key(unparsed_filter):
    """Creates a standardized filter string from a valid unparsed_filter.

    :param unparsed_filter: Formatted MySQL WHERE clause.
    :type unparsed_filters: str
    :returns: Standardized MySQL conditional string.
    :rtype: str
    """
    parsed_filter = parse_filter(unparsed_filter)
    column = ".".join(parsed_filter[0:2])
    
    if isinstance(parsed_filter[3], list):
        parsed_filter[3] = ",".join(parsed_filter[3])
        parsed_filter[3] = "".join(["(", parsed_filter[3], ")"])
    filter_right = "".join(parsed_filter[2:])

    filter_key = column + filter_right
    return filter_key

def check_operator(operator, column_object):
    """Validates an operator's application on a MySQL column.

    :param operator: Accepted MySQL operator.
    :type operator: str
    :param column_object: A SQLAlchemy Column object.
    :type column_object: Column
    """
    if operator not in OPERATORS:
        raise ValueError(f"Operator {operator} is not supported.")
    
    if not isinstance(column_object, Column):
        raise TypeError(f"Type {type(column_object)} is not a Column.")

    column_type = column_object.type.python_type

    if column_type not in TYPES:
        raise ValueError(f"Column '{column_object.name}' "
                         f"has an unsupported type, {column_type}.")
    if operator in NUMERIC_OPERATORS and \
       column_type not in COMPARABLE_TYPES:
        raise ValueError(f"Column '{column_object.name}' "
                         f"is not comparable with '{operator}'.")

def parse_cmd_string(unparsed_cmd_string): 
    """Recognizes and parses MySQL WHERE clause structures.

    :param unparsed_cmd_string: Formatted MySQL WHERE clause string.
    :type unparsed_cmd_string: str
    :returns: 2-D array containing lists of statements joined by ORs.
    :rtype: list[list]
    """
    cmd_line_format = re.compile(".+ + AND|.+ + OR|.+")
   
    or_splits = []
    if unparsed_cmd_string == "":
        pass
    elif re.match(cmd_line_format, unparsed_cmd_string) != None:
        or_splits = []

        or_split = (re.split(" +OR", unparsed_cmd_string))
        for split in or_split: 
            and_split = re.split(" +AND", split)

            for index in range(len(and_split)):
                 and_split[index] = parse_out_ends(and_split[index])

                 if and_split[index] == "":
                     raise ValueError (f"Cmd filter string contains invalid "
                                        "OR or AND conjuctions: "
                                       f"'{unparsed_cmd_string}'")

            or_splits.append(and_split)
    else:
        raise ValueError(f"Unsupported cmd filter string: "
                         f"'{unparsed_cmd_string}'")

    return or_splits

def parse_cmd_list(unparsed_string_list):
    """Recognizes and parses MySQL WHERE clause structures from cmd lists.

    :param unparsed_string_list: Formatted MySQL WHERE clause arguments.
    :type unparsed_string_list: list[str]
    :returns: 2-D array containing lists of statements joined by ORs.
    :rtype: list[list]
    """
    joined_string = parse_in_spaces(unparsed_string_list)
    parsed_cmd_filters = parse_cmd_string(joined_string)

    return parsed_cmd_filters



