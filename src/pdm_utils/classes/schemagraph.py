""" Module containing various data structure objects for the handling
and mapping of the characteristics and relationships of a MySQL database.
"""
from sqlalchemy import MetaData
from networkx import Graph
from networkx import shortest_path
import re

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

def setup_graph_nodes(db_graph, metadata):
    for table in metadata.tables.keys():
        table_object = metadata.tables[table]
        db_graph.add_node(table_object.name, table=table_object)
 
def setup_graph_edges(db_graph, metadata):
    for target_table in metadata.tables.keys():
        target_table_obj  = metadata.tables[target_table]
        
        for foreign_key in target_table_obj.foreign_keys:   
            referent_column = foreign_key.column
            referent_table_obj = referent_column.table
            referent_table = referent_table_obj.name

            db_graph.add_edge(target_table, referent_table, key=foreign_key)

def traverse(db_graph, source_table, target_table):
    paths = shortest_path(db_graph, source_table, target_table)

    return paths

class SchemaGraph(Graph):
    def __init__(self, metadata=None):
        super(SchemaGraph, self).__init__() 
    
        if metadata:
            self.setup(metadata)

        self.metadata = metadata
        self.id = ""

    def show_tables(self):
        tables = []
        for table in self.nodes:
            tables.append(table) 

        return tables

    def has_table(self, table):
        tables = self.show_tables()
        if table in tables:
            return True

        return False

    def get_table(self, table):
        if self.has_table(table):
            return self.nodes[table]["table"]

        return None

    def get_column(self, column):
        parsed_column = parse_column(column)
        table = self.get_table(parsed_column[0])
        
        column = None
        if table != None:
            columns_dict = dict(table.columns)
            column = columns_dict[parsed_column[1]]
        
        return column

    def traverse(self, source_table, target_table):
        path = traverse(self, source_table, target_table)
        return path

    def setup(self, metadata):
        if metadata:
            if not isinstance(metadata, MetaData):
                raise TypeError(
                "Metadata object passed is not of SQLAlchemy MetaData type.")

        setup_graph_nodes(self, metadata)
        setup_graph_edges(self, metadata)

