from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler


def setup_db_tree_leaves(sql_handle, db_node):
    table_query = (
        "SELECT DISTINCT TABLE_NAME FROM INFORMATION_SCHEMA.COLUMNS "
       f"WHERE TABLE_SCHEMA='{sql_handle.database}'")
    table_dicts = sql_handle.execute_query(table_query)

    for data_dict in table_dicts:
        table_node = TableNode(data_dict["TABLE_NAME"])
        db_node.add_table(table_node)
        
        column_query = f"SHOW columns IN {data_dict['TABLE_NAME']}"
        column_dicts = sql_handle.execute_query(column_query)
        for column_dict in column_dicts:
            if column_dict["Null"] == "YES":
                Null = True
            else:
                Null = False

            column_node = ColumnNode(column_dict["Field"],
                                     type=column_dict["Type"],
                                     Null=Null,
                                     key=column_dict["Key"])
            table_node.add_column(column_node)
            if(column_dict["Key"] == "PRI"):
                table_node.primary_key = column_node

def setup_db_tree_webs(sql_handle, db_node):
    for table in db_node.children:
        #print(table.id + ": ")
        foreign_key_query = (
              f"SELECT TABLE_NAME, COLUMN_NAME, CONSTRAINT_NAME, "
               "REFERENCED_TABLE_NAME, REFERENCED_COLUMN_NAME "
               "FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE "
              f"WHERE REFERENCED_TABLE_SCHEMA='{sql_handle.database}' AND "
                    f"REFERENCED_TABLE_NAME='{table.id}'")
        foreign_key_results = sql_handle.execute_query(foreign_key_query)
        if foreign_key_results != ():
            for dict in foreign_key_results:
                primary_key_node = table.get_column(
                                             dict["COLUMN_NAME"])
                ref_table = db_node.get_table(
                                             dict["TABLE_NAME"])
                foreign_key_node = ref_table.get_column(
                                             dict["COLUMN_NAME"])

                primary_key_node.add_table(ref_table)
                
                if foreign_key_node == ref_table.primary_key:
                    ref_table.primary_key = primary_key_node

                for parent in foreign_key_node.parents:
                    if parent not in primary_key_node.parents:
                        primary_key_node.add_table(parent)
                        parent.remove_column(foreign_key_node)

                ref_table.remove_column(foreign_key_node)

def setup_grouping_options(sql_handle, db_node):
    limited_sets = ["varchar", "enum", "char", "tinyint"]
    str_sets     = ["blob", "mediumblob", "varchar"]
    num_sets     = ["int", "decimal", "mediumint", "float", "datetime",
                    "double"]

    for table in db_node.children:
        for column in table.children:
            if column == table.primary_key:
                continue
            if column.parse_type() in str_sets:
                column.group = "large_str_set"

            elif column.parse_type() in num_sets:
                column.group = "large_num_set"

            elif column.parse_type() in limited_sets:
                query = f"SELECT COUNT(DISTINCT {column.id}) FROM {table.id}"
                results_dict = sql_handle.execute_query(query)[0]

                if results_dict[f"COUNT(DISTINCT {column.id})"] < 150:
                    column.group = "limited_set"

class DatabaseTree: 
    def __init__(self, sql_handle):
        self.sql_handle = sql_handle 
        self.db_node = DatabaseNode(sql_handle.database)
        
        setup_db_tree_leaves(sql_handle, self.db_node)                
        setup_db_tree_webs(sql_handle, self.db_node) 

        setup_grouping_options(sql_handle, self.db_node)

    def get_table(self, table):
        return self.db_node.get_table(table)

    def has_table(self, table):
        return self.db_node.has_table(table)

    def show_tables(self):
        return self.db_node.show_tables()

    def get_root(self):
        return self.db_node

    def find_path(self, curr_table, target_table, current_path=None):
            if current_path == None:
                #current_path = [[curr_table.id, curr_table.primary_key.id]]
                current_path = []

            previous_tables = []

            for previous in current_path:
                previous_tables.append(previous[0])

            if curr_table == target_table:
                return current_path

            table_links = curr_table.get_foreign_keys()
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
                        path = self.find_path(table, target_table, path)
                        if path != None:
                            if (path[-1])[0] == target_table.id:
                                return path

    def build_values(self, table, column, values_column=None,
                     queries=[], values=[]):
        table_node = self.get_table(table)
        if table_node == None:
            raise ValueError

        if not table_node.has_column(column):
            raise ValueError
        query = f"SELECT {column} FROM {table}"

        if queries or values:
            query = query + " WHERE "

        if queries:
            query = query + " and ".join(queries)
            if values:
                query = query + " and "

        if values and values_column == None:
            query = query +f"{table_node.primary_key.id} IN ('" +\
                              "','".join(values) + "')"
        elif values and values_column != None: 
            if not table_node.has_column(values_column):
                raise ValueError
            query = query +f"{values_column} IN ('" +\
                            "','".join(values) + "')"

        values = []
        for result in self.sql_handle.execute_query(query):
            values.append(result[column])
  
        return values

class Node:
    def __init__(self, id, parents=None, children=None):
        if parents == None:
            self.parents = []
        else:
            self.parents = parents

        if children == None:
            self.children = []
        else:
            self.children = children

        self.id = id

    def add_parent(self, parent_node):
        parent_node.children.append(self)
        self.parents.append(parent_node)

    def add_child(self, child_node):
        child_node.parents.append(self)
        self.children.append(child_node)

    def create_parent(self, parent_id):
        parent_node = Node(parent_id)
        self.add_parent(parent_node)
        return parent_node

    def create_child(self, child_id):
        child_node = Node(child_id)
        self.add_child(child_node)
        return child_node

    def show_parents(self):
        parents = []
        for parent in self.parents:
            parents.append(parent.id)

        return parents

    def show_children(self):
        children = []  
        for child in self.children:
            children.append(child.id)

        return children

    def has_parent(self, parent):
        parents = self.show_parents()
        return (parent in parents)

    def has_child(self, child):
        children = self.show_children()
        return (child in children)

    def get_parent(self, parent_id):
        parent_node = None
        for parent in self.parents:
            if parent.id == parent_id:
                parent_node = parent

        return parent_node

    def get_child(self, child_id):
        child_node = None
        for child in self.children:
            if child.id == child_id:
                child_node = child

        return child_node

    def remove_parent(self, parent_node):
        if parent_node != None:
            self.parents.remove(parent_node)

        return parent_node

    def remove_child(self, child_node):
        if child_node != None:
            self.children.remove(child_node)

        return child_node

class DatabaseNode(Node):
    def __init__(self, id, parents=None, children=None):
        super(DatabaseNode, self).__init__(
                                    id, parents=parents, children=children)
   
    def create_table(self, table):
        return self.create_child(table)

    def add_table(self, table_node):
        return self.add_child(table_node)

    def remove_table(self, table):
        return self.remove_child(table)

    def show_tables(self):
        return self.show_children()

    def has_table(self, table):
        return self.has_child(table)

    def get_table(self, table):
        return self.get_child(table)

class TableNode(Node):
    def __init__(self, id, parents=None, children=None):
        super(TableNode, self).__init__(
                                    id, parents=parents, children=children)

        self.primary_key = None

    def create_column(self, column):
        return self.create_child(column)

    def add_column(self, column_node):
        return self.add_child(column_node)

    def remove_column(self, column):
        return self.remove_child(column)

    def show_columns(self):
        return self.show_children()

    def show_columns_info(self):
        columns_info = []
        for columns in self.children:
            columns_info.append(columns.show_info())

        return columns_info
    
    def print_columns_info(self):
        print(" " + "_"*46 + " ")
        print("|" + " "*46 + "|")
        print("| %-16s | %-10s | %-6s | %s |" % \
                         ("Field", "Type", "Null", "Key"))
     
        print("|" + "-"*46 + "|")
        for column in self.show_columns_info():
            key = column[3]
            if key == "":
                key = "   "
            print("| %-16s | %-10s | %-6s | %s |" % \
                         (column[0], column[1], column[2], key))

        print("|" + "_"*46 + "|")


    def has_column(self, column):
        return self.has_child(column) 

    def get_column(self, column):
        return self.get_child(column)

    def show_foreign_keys(self):
        foreign_keys = []
        for column in self.children:
            if column != self.primary_key:
                if len(column.parents) > 1:
                    foreign_keys.append(column.id)

        return foreign_keys

    def get_foreign_keys(self):
        foreign_keys = []
        for column in self.children:
            if column != self.primary_key:
                if len(column.parents) > 1:
                    foreign_keys.append(column)

        return foreign_keys

    def show_primary_key(self):
        return self.primary_key.id

class ColumnNode(Node):
    def __init__(self, id, parents=None, children=None, 
                 type=None, Null=None, key=None):
        super(ColumnNode, self).__init__(
                                    id, parents=parents, children=children)
        
        self.type = type
        self.null = Null
        self.key = key

        self.group = None

    def create_table(self, table):
        return self.create_parent(table)

    def add_table(self, table_node):
        return self.add_parent(table_node)

    def remove_table(self, table):
        return self.remove_parent(table)

    def get_type(self):
        return self.type
   
    def parse_type(self):
        type = self.type.split("(")
        type = type[0].split(" ")
        
        return type[0]

    def show_info(self):
        return [self.id, self.parse_type(), self.null, self.key]

def print_database_tables(db_node):
    for table in db_node.children:
        print(table.id + ":")
        table.print_columns_info()
        print("\n")
        

if __name__ == "__main__":
    sql_handle = MySQLConnectionHandler()
    sql_handle.database = input("Please enter database name: ")
    sql_handle.get_credentials()
    sql_handle.validate_credentials()
    db_tree = DatabaseTree(sql_handle)
    db_node = db_tree.get_root()
    print_database_tables(db_node)
