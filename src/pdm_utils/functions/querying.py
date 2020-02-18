from sqlalchemy import Column
from sqlalchemy import join
from sqlalchemy import MetaData
from sqlalchemy import select
from sqlalchemy.sql import distinct
from sqlalchemy.sql import func
from sqlalchemy.sql import functions
from sqlalchemy.sql.elements import BinaryExpression
from sqlalchemy.sql.elements import UnaryExpression
from datetime import datetime
from decimal import Decimal
import re

#Global file constants
COMPARATIVE_OPERATORS = [">", ">=", "<", "<="]
OPERATORS             = ["=", "!="] + COMPARATIVE_OPERATORS
COMPARABLE_TYPES      = [int, Decimal, float, datetime]
TYPES                 = [str, bytes] + COMPARABLE_TYPES
GROUP_OPTIONS = ["limited_set", "num_set", "str_set"]

def translate_table(db_graph, raw_table):
    """Parses a case-insensitive string to match a case-sensitive 
    table_node id string.

    :param raw_table:
        Input a case-insensitive string for a TableNode id.
    :type raw_table: str
    :param verbose:
        Set a boolean to control terminal output.
    :type verbose: Boolean
    :return table:
        Returns a case-sensitive string for a TableNode id.
    :type table: str
    """
    for table in db_graph.show_tables():
        if table.lower() == raw_table.lower():
            return table

    raise ValueError(f"Table '{raw_table}' requested to be filtered "
                     f"is not in '{self.sql_handle.database}'")

def translate_column(db_graph, raw_column, table_object):
    """Parses a case-insensitive string to match a case-sensitive
    column_node id string.

    :param raw_field:
        Input a raw string for a ColumnNode id.
    :type raw_field: str
    :param table:
        Input a case-sensitive string for a TableNode id.
    :type table: str
    :param verbose:
        Set a boolean to control the terminal output.
    :type verbose: Boolean
    :return field:
        Returns a case-sensitive string for a ColumnNode id.
    """
    for column in table_object.columns.keys():
        if column.lower() == raw_column.lower():
            return column
    
    raise ValueError(f"Field '{raw_column}' requested to be filtered"
                     f" is not in '{table_object.name}'")

def check_operator(operator, column_object):
    """Parses a operator string to match a MySQL query operators.

    :param operator:
        Input a raw operator string for an accepted MySQL query operator.
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
                         f"is not comparable with operator '{operator}'.")

def parse_filter(unparsed_filter):
    """Helper function to return a two-dimensional array of filter parameters.

    :param unparsed_filters:
        Input a list of filter expressions to parse and split.
    :type unparsed_filters: List[str]
    :return filters:
        Returns a two-dimensional array of filter parameters.
    :type filters: List[List[str]]
    """
    filter_format = re.compile("\w+\.\w+[=<>!]+\w+", re.IGNORECASE)

    if re.match(filter_format, unparsed_filter) != None:
        filter = (re.split("\W+", unparsed_filter) +\
                        re.findall("[!=<>]+", unparsed_filter))
    else:
        raise ValueError(f"Unsupported filtering format: '{unparsed_filter}'")
                
    return filter

def build_whereclause(db_graph, filter_expression): 
    filter_params = parse_filter(filter_expression)

    table = translate_table(db_graph, filter_params[0])
    table_object = db_graph.get_table(table)

    column = translate_column(db_graph, filter_params[1], table_object)  
    column_object = table_object.columns[column]
 
    check_operator(filter_params[3], column_object)

    whereclause = None

    if filter_params[3] == "=":
        whereclause = column_object == filter_params[2]
    elif filter_params[3] == "!=":
        whereclause = column_object != filter_params[2]
    elif filter_params[3] == ">" :
        whereclause = column_object >  filter_params[2]
    elif filter_params[3] == ">=":
        whereclause = column_object >= filter_params[2]
    elif filter_params[3] == "<" :
        whereclause = column_object <  filter_params[2]
    elif filter_params[3] == "<=":
        whereclause = column_object <= filter_params[2]

    return whereclause 

def build_onclause(db_graph, source_table, adjacent_table):
    edge = db_graph[source_table][adjacent_table]
    foreign_key = edge["key"]

    referent_column = foreign_key.column
    referenced_column = foreign_key.parent
    onclause = referent_column == referenced_column

    return onclause
 
def get_table_list(columns):
    table_set = set()
    for column in  columns:
        if isinstance(column, Column): 
            table_set.add(column.table)
        elif isinstance(column, functions.count):
            for column_clause in column.clauses.clauses:
                table_set.add(column_clause.table)
        elif isinstance(column, UnaryExpression):
            table_set.add(column.element.table)
        elif isinstance(column, BinaryExpression):
            expression = column.left
            if isinstance(expression, UnaryExpression):
                table_set.add(unary_expression.element.table)
            elif isinstance(expression, Column):
                table_set.add(expression.table)

        else:
            raise TypeError(f"Column {column} is not a SqlAlchemy Column.")
                            
    table_list = list(table_set)
    return table_list

def get_table_pathing(db_graph, table_list, center_table=None):
    if not center_table:
        center_table = table_list[0]

    table_list = table_list[1:]

    table_paths = []
    for table in table_list:
        path = db_graph.traverse(center_table.name, table.name)

        if not path:
            raise ValueError(f"Table {table_node} is not connected by any "
                             f"foreign keys to table {center_table}.")
        table_paths.append(path)

    table_pathing = [center_table, table_paths]
    return table_pathing

def join_pathed_tables(db_graph, table_pathing):
    center_table = table_pathing[0]
    joined_tables = center_table

    joined_table_set = set()
    joined_table_set.add(center_table.name)

    for path in table_pathing[1]:
        for index in range(len(path)):
            table = path[index]
            previous_table = path[index-1]
            if table not in joined_table_set:
                joined_table_set.add(table)
                table_object = db_graph.get_table(table)

                onclause = build_onclause(db_graph, previous_table, table)
                joined_tables = join(joined_tables, table_object, 
                                                            isouter=True, 
                                                            onclause=onclause)

    return joined_tables

def build_fromclause(db_graph, columns):
    table_list = get_table_list(columns)
    table_pathing = get_table_pathing(db_graph, table_list)
    joined_table = join_pathed_tables(db_graph, table_pathing)

    return joined_table

def build_select(db_graph, columns, where=None, order_by=None, 
                                                    from_=None, in_=None):
    if from_ == None:
        where_columns = []
        if where != None:
            for clause in where:
                where_columns.append(clause.left) 

        order_by_columns = []

        total_columns = columns + where_columns + order_by_columns
        fromclause = build_fromclause(db_graph, total_columns) 
    else:
        fromclause = from_

    select_query = select(columns).select_from(fromclause)

    if where != None:
        for clause in where:
            select_query = select_query.where(clause)

    if order_by != None:
        for clause in order_by:
            select_query = select_query.order_by(clause)

    return select_query

def build_count(db_graph, columns, where=None, order_by=None, 
                                                    from_=None, in_=None):    
    if from_ == None:
        where_columns = []

        if where != None:
            for clause in where:
                where_columns.append(clause.left) 

        order_by_columns = []

        total_columns = columns + where_columns + order_by_columns
        fromclause = build_fromclause(db_graph, total_columns)
    else:
        fromclause = from_

    column_params = []
    for column_param in columns:
        column_params.append(func.count(column_param))

    count_query = select(column_params).select_from(fromclause)

    if where != None:
        for clause in where:
            count_query = count_query.where(clause)

    if order_by != None:
        for clause in order_by:
            count_query = count_query.order_by(clause)
    
    return count_query

def build_distinct(db_graph, columns, where=None, order_by=None, 
                                                    from_=None, in_=None):
    query = build_select(db_graph, columns, where=where, 
                                                    order_by=order_by,
                                                    from_=from_,
                                                    in_=in_)

    distinct_query = query.distinct()
    return distinct_query

