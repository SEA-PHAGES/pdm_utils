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
from sqlalchemy.sql.elements import UnaryExpression
from datetime import datetime
from decimal import Decimal
import re

def get_table(metadata, table):
    table = parsing.translate_table(metadata, table)
    table_obj = metadata.tables[table]
    return table_obj

def get_column(metadata, column):
    parsed_column = parsing.parse_column(column)
    table = get_table(metadata, parsed_column[0])
    column = parsing.translate_column(metadata, column)

    columns_dict = dict(table.columns)
    column = columns_dict[parsed_column[1]]

    return column

def build_graph(metadata):
    graph = Graph()
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

def build_whereclause(db_graph, filter_expression): 
    filter_params = parsing.parse_filter(filter_expression)

    table_object = db_graph.nodes[filter_params[0]]["table"]
    column_object = table_object.columns[filter_params[1]]
 
    parsing.check_operator(filter_params[2], column_object)

    whereclause = None

    if filter_params[2] == "=":
        whereclause = (column_object  ==  filter_params[3])
    elif filter_params[2] == "LIKE":
        whereclause = (column_object.like(filter_params[3]))
    elif filter_params[2] == "!=" or filter_params[2] == "IS NOT":
        whereclause = (column_object  !=  filter_params[3])
    elif filter_params[2] == ">" :
        whereclause = (column_object  >   filter_params[3])
    elif filter_params[2] == ">=":
        whereclause = (column_object  >=  filter_params[3])
    elif filter_params[2] == "<" :
        whereclause = (column_object  <   filter_params[3])
    elif filter_params[2] == "<=":
        whereclause = (column_object  <=  filter_params[3])

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
        path = shortest_path(db_graph, center_table.name, table.name)

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

def build_select(db_graph, columns, where=None, order_by=None):
    where_columns = []
    if where != None:
        for clause in where:
            where_columns.append(clause.left) 

    order_by_columns = []

    total_columns = columns + where_columns + order_by_columns
    fromclause = build_fromclause(db_graph, total_columns) 

    select_query = select(columns).select_from(fromclause)

    if where != None:
        for clause in where:
            select_query = select_query.where(clause)

    if order_by != None:
        for clause in order_by:
            select_query = select_query.order_by(clause)

    return select_query

def build_count(db_graph, columns, where=None, order_by=None):
    where_columns = []

    if where != None:
        for clause in where:
            where_columns.append(clause.left) 

    order_by_columns = []

    total_columns = columns + where_columns + order_by_columns
    fromclause = build_fromclause(db_graph, total_columns)

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

def build_distinct(db_graph, columns, where=None, order_by=None):
    query = build_select(db_graph, columns, where=where, 
                                                    order_by=order_by)

    distinct_query = query.distinct()
    return distinct_query

