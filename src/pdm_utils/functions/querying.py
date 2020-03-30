from networkx import Graph
from networkx import shortest_path
from pdm_utils.functions import parsing
from sqlalchemy import Column
from sqlalchemy import join
from sqlalchemy import MetaData
from sqlalchemy import select
from sqlalchemy.sql import distinct
from sqlalchemy.sql import func
from sqlalchemy.sql import functions
from sqlalchemy.sql.elements import BinaryExpression
from sqlalchemy.sql.elements import Grouping
from sqlalchemy.sql.elements import UnaryExpression
from datetime import datetime
from decimal import Decimal
import re

#GLOBAL CONSTANTS
COLUMN_TYPES = [Column, functions.count, BinaryExpression, UnaryExpression]

def get_table(metadata, table):
    table = parsing.translate_table(metadata, table)
    table_obj = metadata.tables[table]
    return table_obj

def get_column(metadata, column):
    parsed_column = parsing.parse_column(column)
    table_obj = get_table(metadata, parsed_column[0])
    column = parsing.translate_column(metadata, column)

    columns_dict = dict(table_obj.columns)
    column = columns_dict[parsed_column[1]]

    return column

def build_graph(metadata):
    if not isinstance(metadata, MetaData):
        raise TypeError("Graph requires MetaData object, "
                       f"object passed was of type {type(metadata)}.")

    graph = Graph(metadata=metadata)
    for table in metadata.tables.keys():
        table_object = metadata.tables[table]
        graph.add_node(table_object.name, table=table_object) 

    for target_table in metadata.tables.keys():
        target_table_obj  = metadata.tables[target_table]
    
        for foreign_key in target_table_obj.foreign_keys:   
            referent_column = foreign_key.column
            referent_table_obj = referent_column.table
            referent_table = referent_table_obj.name

            graph.add_edge(target_table, referent_table, 
                                         key=foreign_key)

    return graph

def build_where_clause(db_graph, filter_expression): 
    filter_params = parsing.parse_filter(filter_expression)

    table_object = db_graph.nodes[filter_params[0]]["table"]
    column_object = table_object.columns[filter_params[1]]
 
    parsing.check_operator(filter_params[2], column_object)

    right = filter_params[3]

    if column_object.type.python_type == bytes:
        right = right.encode("utf-8")

    where_clause = None

    if filter_params[2] == "=":
        where_clause = (column_object  ==  right)
    elif filter_params[2] == "LIKE":
        where_clause = (column_object.like(right))
    elif filter_params[2] == "!=" or filter_params[2] == "IS NOT":
        where_clause = (column_object  !=  right)
    elif filter_params[2] == ">" :
        where_clause = (column_object  >   right)
    elif filter_params[2] == ">=":
        where_clause = (column_object  >=  right)
    elif filter_params[2] == "<" :
        where_clause = (column_object  <   right)
    elif filter_params[2] == "<=":
        where_clause = (column_object  <=  right)

    return where_clause 

def build_onclause(db_graph, source_table, adjacent_table):
    source_table = parsing.translate_table(db_graph.graph["metadata"], 
                                           source_table)
    adjacent_table = parsing.translate_table(db_graph.graph["metadata"],
                                           adjacent_table)

    if source_table == adjacent_table:
        source_table = db_graph.nodes[source_table]["table"]
        return source_table

    edge = db_graph[source_table][adjacent_table]
    foreign_key = edge["key"]

    referent_column = foreign_key.column
    referenced_column = foreign_key.parent
    onclause = referent_column == referenced_column

    return onclause

def extract_column(column, check=None):
    if check != None:
        if not isinstance(column, check):
            raise TypeError(f"Type {type(column)} object passed as column "
                             "failed a specified type check.")

    if isinstance(column, Column): 
        pass
    #For handling SQLAlchemy count(Column) expressions
    elif isinstance(column, functions.count):
        column = column.clauses.clauses[0]
    #For handling SQLAlchemy Column.distinct expressions
    elif isinstance(column, UnaryExpression):
        column = column.element
    #For handling SQLAlchemy comparison expressions
    elif isinstance(column, BinaryExpression):
        expression = column.left
        #For handling SQLAlchemy Column.in_() expressions
        if isinstance(expression, UnaryExpression):        
            column = expression.element
        #For handling SQLAlchemy count(Column) comparisons
        elif isinstance(expression, functions.count):
            column = expression.clauses.clauses[0]
        #For handling SQLAlchemy Column.distinct() comparisons
        elif isinstance(expression, Grouping):
            column = expression.element.element
        elif isinstance(expression, Column):
            column = expression
        else:
            raise TypeError(f"BinaryExpression type {type(expression)} "
                            f"of expression {expression} is not supported.")

    else:
        raise TypeError(f"Input type {column} is not a derivative "
                         "of a SqlAlchemy Column.")
    
    return column

def extract_columns(columns, check=None):
    extracted_columns = []

    if isinstance(columns, list):
        for column in columns:
            extracted_columns.append(extract_column(column, check=check))

    else:
        extracted_columns.append(extract_column(columns, check=check))

    return extracted_columns

def extract_where_clauses(where_clauses):
    where_columns = []
    if not where_clauses is None:
        where_columns = extract_columns(where_clauses, check=BinaryExpression)

    return where_columns

def extract_order_by_clauses(order_by_clauses):
    order_by_columns = []
    if not order_by_clauses is None:
        order_by_columns = extract_columns(order_by_clauses, check=Column)

    return order_by_columns

def get_table_list(columns):
    table_set = set()
   
    if not isinstance(columns, list):
        columns = [columns]

    for column in columns:
        table_set.add(extract_column(column).table)
                            
    table_list = list(table_set)
    return table_list

def get_table_pathing(db_graph, table_list, center_table=None):
    table_list = table_list.copy()

    if not center_table is None:
        table_list.remove(center_table)

    else:
        center_table = table_list[0]
        table_list.remove(center_table)


    table_paths = []
    for table in table_list:
        path = shortest_path(db_graph, center_table.name, table.name)

        if not path:
            raise ValueError( "Operation cannot be performed. "
                             f"Table {table.name} is not connected by any "
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
                table_object = db_graph.nodes[table]["table"]

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

def append_where_clauses(executable, where_clauses):
    if where_clauses is None:
        return executable

    if isinstance(where_clauses, list):
        for clause in where_clauses:
            executable = executable.where(clause)
    else:
        executable = executable.where(where_clauses)

    return executable
    
def append_order_by_clauses(executable, order_by_clauses):
    if order_by_clauses is None:
        return executable

    if isinstance(order_by_clauses, list):
        for clause in order_by_clauses:
            executable = executable.order_by(clause)
    else:
        executable = executable.order_by(order_by_clauses)

    return executable

def build_select(db_graph, columns, where=None, order_by=None):
    where_columns = extract_where_clauses(where)
    order_by_columns = extract_order_by_clauses(order_by)

    if not isinstance(columns, list):
        columns = [columns]

    total_columns = columns + where_columns + order_by_columns
    fromclause = build_fromclause(db_graph, total_columns) 

    select_query = select(columns).select_from(fromclause)
    select_query = append_where_clauses(select_query, where)
    select_query = append_order_by_clauses(select_query, order_by)

    return select_query

def build_count(db_graph, columns, where=None, order_by=None):
    where_columns = extract_where_clauses(where)
    order_by_columns = extract_order_by_clauses(order_by)
   
    if not isinstance(columns, list):
        columns = [columns]

    total_columns = columns + where_columns + order_by_columns
    fromclause = build_fromclause(db_graph, total_columns)

    column_params = []
    for column_param in columns:
        column_params.append(func.count(column_param))

    count_query = select(column_params).select_from(fromclause)
    count_query = append_where_clauses(count_query, where)
    count_query = append_order_by_clauses(count_query, order_by)

    return count_query

def build_distinct(db_graph, columns, where=None, order_by=None):
    query = build_select(db_graph, columns, where=where, 
                                                    order_by=order_by)

    distinct_query = query.distinct()
    return distinct_query

