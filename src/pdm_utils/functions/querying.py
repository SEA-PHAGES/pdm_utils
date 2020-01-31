from sqlalchemy import select, MetaData, join, Column
from sqlalchemy.sql import func, distinct

def get_table_set(columns):
    table_set = set()
    for column in  columns:
        if isinstance(column, Column): 
            table_set.add(column.table)
        elif isinstance(column, functions.count):
            for column_clause in column.clauses.clauses:
                table_set.add(column_clause.table)
        else:
            raise TypeError(f"Column {column} is not a SqlAlchemy Column.")
                            

    return table_set

def get_table_nodes(db_graph, table_set):
    table_nodes = []
    for table in table_set:
        table_node = db_graph.get_table(table.name)

        if not table_node:
            raise ValueError(f"Table {table} is not in {db_graph.id}")

        table_nodes.append(table_node)

    return table_nodes

def get_table_pathing(db_graph, table_nodes, center_table=None):
    if not center_table:
        center_table = table_nodes[0]

    table_nodes = table_nodes[1:]

    table_paths = []
    for table_node in table_nodes:
        path = db_graph.traverse(center_table, table_node)

        if not path:
            raise ValueError(f"Table {table_node} is not connected by any "
                             f"foreign keys to table {center_table}.")

        truncated_path = []
        for stop in path:
            truncated_path.append(stop[0])
        truncated_path.remove(truncated_path[0])

        table_paths.append(truncated_path)
   
    pathing = [center_table, table_paths]
    return pathing

def join_pathed_tables(db_graph, table_pathing):
    joined_table = table_pathing[0].table

    joined_table_set = set()
    for path in table_pathing[1]:
        for table in path:
            if table not in joined_table_set:
                joined_table_set.add(table)
                table_node = db_graph.get_table(table)

                if not table_node:
                    raise ValueError(f"Table {table} is not in {db_graph.id}")
    
                joined_table = join(joined_table, table_node.table, isouter=True)

    return joined_table

def build_fromclause(db_graph, columns):
    table_set = get_table_set(columns)
    table_nodes = get_table_nodes(db_graph, table_set)
    table_pathing = get_table_pathing(db_graph, table_nodes)
    joined_table = join_pathed_tables(db_graph, table_pathing)

    return joined_table

def build_select(db_graph, columns, where_clause=None, order_by_clause=None):
    fromclause = build_fromclause(db_graph, columns) 

    select_query = select(columns).\
                                select_from(fromclause).\
                                where(where_clause).\
                                order_by(order_by_clause)

    return select_query

def build_count(db_graph, columns, where_clause=None, order_by_clause=None):    
    fromclause = build_fromclause(db_graph, columns) 


    column_params = []
    for column_param in columns:
        column_params.append(func.count(column_param))

    count_query = select(column_params).\
                                select_from(fromclause).\
                                where(where_clause).\
                                order_by(order_by_clause)

    return count_query

def build_distinct(db_graph, columns, where_clause=None, order_by_clause=None):
    query = build_select(db_graph, columns, where_clause=where_clause, 
                                            order_by_clause=order_by_clause)

    distinct_query = query.distinct()
    return distinct_query

