""" Module containing various data structure objects for the handling
and mapping of the characteristics and relationships of a MySQL database.
"""

from sqlalchemy import MetaData
import re

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
        filter = (re.split("\W+", filter) +\
                        re.findall("[!=<>]+", filter))
    else:
        raise ValueError(f"Unsupported filtering format: '{filter}'")
                
    return filter

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

def setup_column_node_group(column_node):
    parsed_type = column_node.parse_type()

    str_set = ["BLOB", "MEDIUMBLOB", "VARCHAR"]
    num_set = ["MEDIUMINT", "INTEGER", "FLOAT", "DECIMAL", "DOUBLE"]
    limited_set = ["TINYINT", "ENUM", "CHAR"]
    datetime_set = ["DATETIME"]

    if parsed_type in str_set:
        column_node.group = "string"

        if parsed_type == "VARCHAR":
            if column_node.type.length <= 5:
                column_node.group = "limited"

    elif parsed_type in num_set:
        column_node.group = "numeric"

    elif parsed_type in limited_set:
        column_node.group = "limited"

    elif parsed_type in datetime_set:
        column_node.group = "datetime"

def setup_graph_nodes(db_graph, metadata):
    for table in metadata.tables.keys():
        table_node = TableNode(table, table=metadata.tables[table])
        db_graph.add_table(table_node)

        for column in table_node.table.columns.keys():
            column_object = table_node.table.columns[column]
            type = column_object.type

            column_node = ColumnNode(column, type=type,
                                        nullable=column_object.nullable, 
                                        primary_key=column_object.primary_key,
                                        column=column_object)
            setup_column_node_group(column_node)

            table_node.add_column(column_node)
            if column_node.column.primary_key:
                table_node.primary_key = column_node

def setup_graph_edges(db_graph, metadata):
    for table_node in db_graph.children:
        table_object = table_node.table

        for column in table_object.columns.keys():
            foreign_key_set = table_object.columns[column].foreign_keys

            for foreign_key in foreign_key_set:
                reference = foreign_key.target_fullname
                parsed_reference = parse_column(reference)
                
                ref_table_node = db_graph.get_table(parsed_reference[0])
                ref_column_node = ref_table_node.get_column(parsed_reference[1])          
                column_node = table_node.get_column(column)

                table_node.remove_child(column_node)
                table_node.add_column(ref_column_node)
                if table_node.primary_key == column_node:
                    table_node.primary_key = ref_column_node

def traverse(curr_table, target_table, current_path=None):
    """Recursive function that finds a path between two TableNodes.

    :param curr_table:
        Input a starting TableNode under a DatabaseNode.
    :type curr_table: TableNode
    :param target_table:
        Input a target TableNode under a DatabaseNode.
    :type target_table: TableNode
    :param current_path:
        Input a 2-D array containing the traversed nodes in a path.
    :returns current_path:
        Returns a 2-D array containing the tables and keys in the path.
    :type current_path: List[List[str]] 
    """
    if current_path == None:
        current_path = [[curr_table.id, curr_table.primary_key.id]]

    previous_tables = []

    for previous in current_path:
        previous_tables.append(previous[0])

    if curr_table == target_table:
        return current_path

    table_links = curr_table.get_linked_keys()
    if len(curr_table.primary_key.parents) > 1:
        table_links.insert(0, curr_table.primary_key)

    for link in table_links:
        if target_table in link.parents:
            path = current_path.copy()
            path.append([target_table.id, link.id])

            return path
        for table in link.parents:
            if table.id not in previous_tables and table != curr_table:
                path = current_path.copy()
                path.append([table.id, link.id])

                path = traverse(table, target_table, path)
                if path != None:
                    if (path[-1])[0] == target_table.id:
                        return path

class Node:
    """Object which can store other Node objects and an id."""

    def __init__(self, id, parents=None, children=None):
        """Initializes a Node object.

        :param id:
            Input a ID for a Node object as a string.
        :type id: str
        :param parents:
            Input parents for a Node object as a list of Node objects.
        :type parents: List[Node]
        :param children:
            Input children for a Node object as a list of Node objects.
        :type children: List[Node]
        """
        if parents == None:
            self.parents = []
        else:
            self.parents = parents

        if children == None:
            self.children = []
        else:
            self.children = children

        if not isinstance(id, str):
            raise TypeError("id parameter must be a string object.")
        self.id = id

    def has_parent(self, parent):
        """Function that returns a boolean of if a parent Node exists.

        :param parent:
            Input the name of a parent Node as a string.
        :type parent: str
        :returns (parent in parents)
            Returns a boolean expression of if a ID matches any parent Node.
        :type (parent in parents): Boolean
        """
        if not isinstance(parent, str):
            raise TypeError("parent parameter must be a string object")

        parents = self.show_parents()
        return (parent in parents)

    def has_child(self, child):
        """Function that returns a boolean of if a child Node exists.

        :param child:
            Input the name of a child Node as a string.
        :type child: str
        :returns (child in children)
            Returns a boolean expression of if a ID matches any child Node.
        :type (child in children): Boolean
        """
        if not isinstance(child, str):
            raise TypeError("child paremeter must be a string object")
        children = self.show_children()
        return (child in children)

    def get_parent(self, parent_id):
        """Function that returns a connected parent Node.

        :param parent_id:
            Input the ID of a parent Node as a string.
        :type child_id: str
        :returns parent_node
            Returns a parent Node with a matching ID.
        :type parent_node: Node
        """
        if not isinstance(parent_id, str):
            raise TypeError("parent_id parameter must be a string object")

        parent_node = None
        for parent in self.parents:
            if parent.id == parent_id:
                parent_node = parent

        return parent_node

    def get_child(self, child_id):
        """Function that returns a connected child Node.

        :param child_id:
            Input the ID of a child Node as a string.
        :type child_id: str
        :returns child_node
            Returns a child Node with a matching ID.
        :type child_node: Node
        """
        if not isinstance(child_id, str):
            raise TypeError("parent_id parameter must be a string object")

        child_node = None
        for child in self.children:
            if child.id == child_id:
                child_node = child

        return child_node

    def show_parents(self):
        """Function that returns a list representing the all parent Nodes.

        :returns parents:
            Returns a list of IDs representing the parent Nodes.
        :type show_parents: List[str]
        """
        parents = []
        for parent in self.parents:
            parents.append(parent.id)

        return parents

    def show_children(self):
        """Function that returns a list representing all the children Nodes.
        
        :returns children:
            Returns a list of IDs representing the children Nodes.
        :type children: List[str]
        """
        children = []  
        for child in self.children:
            children.append(child.id)

        return children

    def add_parent(self, parent_node):
        """Function to link an existing Node object as a parent.

        :param parent_node:
            Input a Node object.
        :type parent_node: Node
        """
        if not isinstance(parent_node, Node):
            raise TypeError("parent_node parameter must be a Node object.")
        
        if parent_node in self.parents:
            raise ValueError(
            "parent_node cannot be a duplicate of an existing parent Node.")

        if self.has_parent(parent_node.id):
            raise ValueError(
            "parent_node cannot have the id of an existing parent Node.")

        parent_node.children.append(self)
        self.parents.append(parent_node)

    def add_child(self, child_node):
        """Function to link an existing Node object as a child.

        :param child_node:
            Input a Node object.
        :type child_node: Node
        """
        if not isinstance(child_node, Node):
            raise TypeError("child_node parameter must be a Node object.")

        if child_node in self.parents:
            raise ValueError(
            "child_node cannot be a duplicate of an existing child Node.")

        if self.has_child(child_node.id):
            raise ValueError(
            "child_node cannot have the id of and existing child Node.")

        child_node.parents.append(self)
        self.children.append(child_node)

    def create_parent(self, parent_id):
        """Function to create and link a new Node object as a parent.

        :param parent_id:
            Input an ID for the new Node object as a string.
        :type parent_id: str
        :returns parent_node:
            Returns the created and linked Node object.
        :type parent_node: Node
        """
        if not isinstance(parent_id, str):
            raise TypeError("parent_id parameter must be a string object")

        parent_node = Node(parent_id)
        self.add_parent(parent_node)
        return parent_node

    def create_child(self, child_id):
        """Function to create and link a new Node object as a child.
        
        :param child_id:
            Input an ID for the new Node object as a string.
        :type child_id: str
        :returns child_node:
            Returns the created and linked Node object.
        :type child_id: Node
        """
        if not isinstance(child_id, str):
            raise TypeError("child_id parameter must be a string object")

        child_node = Node(child_id)
        self.add_child(child_node)
        return child_node
  
    def remove_parent(self, parent_node):
        """Function that disconnects a connected parent Node.

        :param parent_node:
            Input the parent Node.
        :type parent_node: Node
        :returns parent_node
            Returns the removed parent_node.
        :type parent_node: Node
        """
        if not isinstance(parent_node, Node):
            raise TypeError

        if parent_node != None:
            self.parents.remove(parent_node)
            parent_node.children.remove(self)

        return parent_node

    def remove_child(self, child_node):
        """Function that disconnects a connected child Node.

        :param child_node:
            Input the child Node.
        :type child_node: Node
        :returns child_node
            Returns the removed child_node.
        :type child_node: Node
        """
        if not isinstance(child_node, Node):
            raise TypeError

        if child_node != None:
            self.children.remove(child_node)
            child_node.parents.remove(self)

        return child_node

class SchemaGraph(Node):
    """Object which can store TableNode objects and an id """
    def __init__(self, id, table_nodes=None):
        super(SchemaGraph, self).__init__(id, parents=None, 
                                               children=table_nodes)

    def add_table(self, table_node):
        """Function to link an existing TableNode object as a child.

        :param table_node:
            Input a Node object.
        :type table_node: Node
        """
        if not isinstance(table_node, TableNode):
            raise TypeError("table_node parameter must be a TableNode object.")
        self.add_child(table_node)

    def create_table(self, table):
        """Function to create and link a new TableNode object as a child.
        
        :param table:
            Input an ID for the new TableNode object as a string.
        :type table: str
        :returns child_id:
            Returns the created and linked Node object.
        :type child_id: Node
        """
        if not isinstance(table, str):
            raise TypeError("table parameter must be a string object.")

        table_node = TableNode(table)
        self.add_table(table_node)

        return table_node
 
    def show_tables(self):
        """Function that returns a list representing all the TableNodes.
        
        :returns children:
            Returns a list of IDs representing the TableNodes.
        :type children: List[str]
        """
        return self.show_children()

    def has_table(self, table):
        """Function that returns a boolean of if a TableNode exists.

        :param table:
            Input the name of a TableNode as a string.
        :type table: str
        :returns (child in children)
            Returns a boolean expression of if a ID matches any TableNode.
        :type (child in children): Boolean
        """
        return self.has_child(table)

    def get_table(self, table):
        """Function that returns a connected TableNode.

        :param table:
            Input the ID of a TableNode as a string.
        :type table: str
        :returns child_node
            Returns a child Node with a matching ID.
        :type child_node: Node
        """
        return self.get_child(table)

    def print_info(self):
        """Function to display information about a MySQL SchemaGraph.

        :param db_node:
            Input a DatabaseNode root of a SchemaGraph.
        :type db_node: DatabaseNode
        """
        for table in self.children:
            print("| " + table.id + " |")
            table.print_columns_info()
            print("\n")
    
    def traverse(self, curr_table, target_table):
        path = traverse(curr_table, target_table)
        return path

    def setup(self, metadata):
        setup_graph_nodes(self, metadata)
        setup_graph_edges(self, metadata)

class TableNode(Node):
    """Object which can store ColumnNode objects, an id, and table info."""

    def __init__(self, id, database_nodes=None, column_nodes=None, table=None):
        super(TableNode, self).__init__(id, parents=database_nodes, 
                                            children=column_nodes)

        self.primary_key = None

        self.table = table

    def add_column(self, column_node):
        """Function to link an existing Node object as a child.

        :param column_node:
            Input a Node object.
        :type column_node: ColumnNode
        """
        if not isinstance(column_node, ColumnNode):
            raise TypeError("column_node parameter must be a ColumnNode object.")
        self.add_child(column_node)

    def create_column(self, column):
        """Function to create and link a ColumnNode object.

        :param column:
            Input an ID for the ColumnNode object as a string.
        :type column: str
        :returns column_node:
            Returns the created and linked ColumnNode object.
        :type column_node: ColumnNode
        """
        if not isinstance(column, str):
            raise TypeError("column parameter must be a string object.")
        column_node = ColumnNode(column)
        self.add_column(column_node)

        return column_node

    def show_columns(self):
        """Function that returns a list representing all the ColumnNodes.
        
        :returns children:
            Returns a list of IDs representing the ColumnNodes.
        :type children: List[str]
        """
        return self.show_children()

    def show_columns_info(self):
        """Function that returns a list of info of attatched ColumnNodes.

        :returns columns_info:
            Returns a list of lists containing info about the MySQL column.
        :type columns_info: List[List[str]]
        """
        columns_info = []
        for columns in self.children:
            columns_info.append(columns.show_info())

        return columns_info
    
    def print_columns_info(self):
        """"Function that formats and prints a list of Columns' info."""

        print(" " + "_"*61 + " ")
        print("|" + " "*61 + "|")
        print("| %-16s | %-10s | %-12s | %-5s | %-4s |" % \
              ("Field", "Type", "Grouping", "Null", "Pri"))
     
        print("|" + "-"*61 + "|")
        for column in self.show_columns_info():
            primary = "    "
            if column[4]:
                primary = True

            print("| %-16s | %-10s | %-12s | %-5s | %-4s |" % \
            (column[0], column[1], column[2], column[3], primary))

        print("|" + "_"*61 + "|")

    def has_column(self, column):
        """Function that returns a boolean of if a child Node exists.

        :param column:
            Input the name of a ColumnNode as a string.
        :type column: str
        :returns (child in children)
            Returns a boolean expression of if a ID matches any ColumnNode.
        :type (child in children): Boolean
        """
        return self.has_child(column) 

    def get_column(self, column):
        """Function that returns a connected ColumnNode.

        :param column:
            Input the name of a ColumnNode as a string.
        :type column: str
        :returns child_node
            Returns a child Node with a matching ID.
        :type child_node: Node
        """

        return self.get_child(column)

    def show_foreign_keys(self):
        """Function that returns a list representing Foreign Key ColumnNodes.
        
        :returns foreign_keys:
            Returns a list of Foreign Key ColumnNode names as strings.
        :type foreign_keys: List[str]
        """
        foreign_keys = []
        for column in self.children:
            if column != self.primary_key:
                if not isinstance(column, ColumnNode):
                    raise TypeError(
                    "Object in TableNode.children is not a ColumnNode object.")
                if len(column.parents) > 1:
                    foreign_keys.append(column.id)

        return foreign_keys

    def get_linked_keys(self):
        """Function that returns a list of ColumnNodes connected
        to other tables.

        :returns foreign_keys:
            Returns a list of attatched Foreign Key ColumnNodes.
        :type foreign_keys: List[ColumnNode]
        """
        foreign_keys = []
        for column in self.children:
            if column != self.primary_key:
                if not isinstance(column, ColumnNode):
                    raise TypeError(
                    "Object in TableNode.children is not a ColumnNode object.")
                if len(column.parents) > 1:
                    foreign_keys.append(column)

        return foreign_keys

    def show_primary_key(self):
        """Function that returns the name of the Primary Key ColumnNode.

        :returns id:
            Returns the name of the primary key ColumnNode of the TableNode.
        :type id: str
        """
        if self.primary_key:
            return self.primary_key.id

        return ""

class ColumnNode(Node):
    def __init__(self, id, table_nodes=None, column=None,
                 type="", nullable=None, primary_key=False, group="Undefined"):
        super(ColumnNode, self).__init__(id, parents=table_nodes, 
                                             children=None)
        
        self.type = type
        self.nullable = nullable
        self.primary_key = primary_key

        self.group = group

        self.column = column
    
    def add_table(self, table_node):
        """Function to link an existing TableNode object as a parent.

        :param table_node:
            Input a TableNode object.
        :type table_node: TableNode
        """
        if not isinstance(table_node, TableNode):
            raise TypeError("table_node parameter must be a TableNode object.")

        self.add_parent(table_node)
    
    def create_table(self, table):
        """Function to create and link a new TableNode object as a parent.

        :param table:
            Input an ID for the new TableNode object as a string.
        :type table_id: str
        :returns table_node:
            Returns the created and linked TableNode object.
        :type table_node: TableNode
        """
        if not isinstance(table, str):
            raise TypeError("table parameter must be a string object.")
        table_node = TableNode(table)
        self.add_table(table_node) 

        return table_node

    def get_type(self):
        """Function that returns the raw column type.

        :returns type:
            Returns the column type as a string.
        :type type: str
        """
        return self.type
   
    def parse_type(self):
        """Function that returns a truncated version of the column type.

        :returns type:
            Returns the a truncated version of the column type as a string.
        :type type: str
        """
        column_type = str(self.type)

        type_format = re.compile("\w+")
        if re.match(type_format, column_type) != None:
            parsed_type = (re.split("\W+", column_type))


        return parsed_type[0]
       
    def show_info(self):
        """Function that returns a list of the info in the ColumnNode.

        :returns info:
            Returns a list of all the info in the ColumnNode as strings.
        :type info: List[str]
        """
        info = [self.id, self.parse_type(), self.group, 
                         self.nullable, self.primary_key]
        return info        

