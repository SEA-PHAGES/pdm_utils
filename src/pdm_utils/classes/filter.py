"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import mysqldb
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from typing import List
import cmd, readline, argparse, os, sys, re, string

#Global file constants
GROUP_OPTIONS = ["cluster2", "subcluster2",
                 "status",  "hoststrain"]

GROUP_OPTIONS_DICT = {"phage" : ["cluster2", "subcluster2",
                                 "status",  "hoststrain"],
                      "gene"  : []}

def build_tables(sql_handle):
    query = (
        "SELECT DISTINCT TABLE_NAME FROM INFORMATION_SCHEMA.COLUMNS "
       f"WHERE TABLE_SCHEMA='{sql_handle.database}'")
    table_dicts = sql_handle.execute_query(query)

    tables = {}
    for data_dict in table_dicts:
        table = data_dict["TABLE_NAME"]
        column_dicts = sql_handle.execute_query(f"SHOW columns IN {table}")
        column_list = []
        for column_dict in column_dicts:
            column_list.append((column_dict["Field"]))
        tables.update({table.lower(): column_list})

    del tables["version"]
    return tables

def build_primary_keys(table_list, sql_handle):
    query = (
        "SELECT "
        "KCU.TABLE_NAME AS Table_Name, "
        "KCU.CONSTRAINT_NAME AS Constraint_Name, "
        "KCU.COLUMN_NAME AS Column_Name "
        "FROM INFORMATION_SCHEMA.TABLE_CONSTRAINTS TC "
        "JOIN INFORMATION_SCHEMA.KEY_COLUMN_USAGE KCU "
        "ON KCU.CONSTRAINT_SCHEMA = TC.CONSTRAINT_SCHEMA "
        "AND KCU.CONSTRAINT_NAME = TC.CONSTRAINT_NAME "
        "AND KCU.TABLE_NAME = TC.TABLE_NAME "
        "WHERE TC.CONSTRAINT_TYPE = 'PRIMARY KEY' "
        "and TC.TABLE_NAME in ('" +\
        "','".join(table_list) + "') "
        "ORDER BY "
        "KCU.TABLE_SCHEMA, KCU.TABLE_NAME, KCU.CONSTRAINT_NAME"
        )
    primary_key_dicts = sql_handle.execute_query(query)

    primary_keys = {}
    for dict in primary_key_dicts:
        primary_keys.update({dict["Table_Name"] : dict["Column_Name"]})

    return primary_keys

def build_foreign_keys(tables, sql_handle):
    query = (
        "SELECT "
        "INFORMATION_SCHEMA.KEY_COLUMN_USAGE.REFERENCED_TABLE_NAME, "
        "INFORMATION_SCHEMA.KEY_COLUMN_USAGE.REFERENCED_COLUMN_NAME, "
        "INFORMATION_SCHEMA.KEY_COLUMN_USAGE.TABLE_NAME, "
        "INFORMATION_SCHEMA.KEY_COLUMN_USAGE.COLUMN_NAME "
        "FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE "
        "WHERE INFORMATION_SCHEMA.KEY_COLUMN_USAGE.REFERENCED_TABLE_NAME IS NOT NULL "
        "AND INFORMATION_SCHEMA.KEY_COLUMN_USAGE.TABLE_SCHEMA= "
        f"'{sql_handle.database}' "
        "AND REFERENCED_TABLE_NAME IN ('" +\
        "','".join(tables) + "') "
        "AND TABLE_NAME IN ('" +\
        "','".join(tables) + "') "
        "ORDER BY INFORMATION_SCHEMA.KEY_COLUMN_USAGE.REFERENCED_TABLE_NAME")
    foreign_key_dicts = sql_handle.execute_query(query)

    foreign_keys = {}
    for table in tables:
        foreign_keys.update({table: []})

    for dict in foreign_key_dicts:
        foreign_keys[dict["REFERENCED_TABLE_NAME"]].append([
            dict["REFERENCED_COLUMN_NAME"],
            dict["TABLE_NAME"],
            dict["COLUMN_NAME"]])

    return foreign_keys

def build_values_set(table, attribute, sql_handle):
    query = f"SELECT {attribute} FROM {table}"
    value_dicts = sql_handle.execute_query(query)

    values = []
    for dict in value_dicts:
        values.append(dict[attribute])

    return values

class Filter:
    def __init__(self, sql_handle, values_list=[], table='phage'):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.sql_handle = sql_handle

        self.db_tree = DatabaseForeignKeyTree(sql_handle)

        self.tables = self.db_tree.tables_dict
        self.primary_keys = self.db_tree.primary_keys

        if table not in self.tables.keys():
            print(f"Table passed to filter is not in {sql_handle.database}")
            raise ValueError

        if values_list == []:
           self.values = list(build_values_set(table,
                                              self.primary_keys[table],
                                              sql_handle))
        else:
           self.values = values_list
        self.key = self.primary_keys[table]

        self.table = table
        self.history = []
        self.filters = {}

    def translate_field(self, raw_field, verbose=False):
        for table in self.tables.keys():
            for field in self.tables[table]:
                if field.lower() == raw_field.lower():
                    return field

        if verbose:
                    print(
                    f"Field '{raw_field}' requested to be filtered"
                    f" is not in '{sql_handle.database}'")

    def add_filter(self, table, raw_field, value, verbose=False):
        """
        Adds a filter to database queries
        """
        field = self.translate_field(raw_field)
        if field == self.primary_keys[table]:
            if verbose:
                print(f"Primary key {field} in {table} "
                      "cannot be a field requested for filtering")
        elif table in self.tables.keys():
            if field in self.tables[table]:
                if table.lower() in self.filters.keys():
                    if field in (self.filters[table]).keys():
                        (self.filters[table])[field].append(value)
                    else:
                        (self.filters[table])[field] = [value]
                else:
                    self.filters.update({table: {field: [value]}})
            else:
                if verbose:
                    print(
                    f"Field '{field}' requested to be filtered"
                    f" is not in table '{table}'")
        else:
            if verbose:
                print(
                f"Table '{table}' requested to be filtered "
                f"is not in '{self.sql_handle.database}'")

    def refresh(self):
        """
        Refreshes values
        """
        query = (f"SELECT {self.key} "
                 f"FROM {self.table} WHERE {self.key} IN "
                 "('" + "','".join(self.values) + "')")

        results_dicts = self.sql_handle.execute_query(query)
        results = []
        for result in results_dicts:
            results.append(result[self.key])

        self.values = results

    def update(self, verbose=False):
        """
        Updates results list for the Filter object
        """
        queries = []
        query_results = self.values
        if len(query_results) == 0:
            return

        for table in self.filters.keys():
            for field in self.filters[table].keys():
                queries.append(
                   f"{field}='"+ \
                    "','".join((self.filters[table])[field]) + "'")
            if table == self.table:
                if queries:
                    query = (f"SELECT {self.key} FROM {table} WHERE " +\
                              " and ".join(queries) + \
                             f" and {self.key} IN ('" + \
                              "','".join(query_results) + "')")

                    if verbose:
                        print(
                        f"Filtering {table} in "
                        f"{self.sql_handle.database} for "+\
                               " and ".join(queries) + "...")

                    query_results=[]
                    for result in self.sql_handle.execute_query(query):
                        query_results.append(result[self.key])

            else:
                    transposed_values = self.transpose(
                                                query_results, table, field,
                                                queries, verbose=verbose)
                    if transposed_values:
                        query_results = transposed_values

            queries = []
            if not query_results:
                break

        self.values = query_results

        if len(query_results) > 0:
            self.history.append(query_results)

    def sort(self, sort_field, verbose=False):
        sort_field = self.translate_field(sort_field)
        if sort_field in self.tables[self.table]:
            if verbose:
                print(f"Sorting by '{sort_field}'...")
            query = (f"SELECT {self.key} "
                     f"FROM {self.table} WHERE {self.key} IN "
                      "('" + "','".join(self.values) + "') "
                     f"ORDER BY {sort_field}")
            query_results = []
            for result in self.sql_handle.execute_query(query):
                query_results.append(result[self.key])

            self.values = query_results

        else:
            if verbose:
                print(f"Sort key '{sort}' not supported.")

    def transpose(self, query_results, table, field, queries, verbose=False):
        traversal_path = self.db_tree.find_path(table, self.table)
        if not traversal_path:
            if verbose:
                print(f"Table {table} is not connected to {self.table} "
                       "by foreign keys")
            return []

        elif (traversal_path[-1])[0] == self.table:
            initial_query = (
                    f"SELECT {self.primary_keys[table]} FROM {table} "
                    f"WHERE " +\
                     " and ".join(queries))
            current_key = self.primary_keys[table]
            current_table = table
            current_results = []
            results_dicts = self.sql_handle.execute_query(initial_query)
            for dict in results_dicts:
                current_results.append(dict[self.primary_keys[table]])

            for connection_key in traversal_path:
                if current_key != connection_key[1]:
                    query = (
                        f"SELECT {connection_key[1]} FROM "
                        f"{current_table} "
                        f"WHERE {current_key} IN ('" +\
                         "','".join(current_results) + "')")
                    results_dicts = self.sql_handle.execute_query(query)
                    current_results = []
                    for dict in results_dicts:
                        current_results.append(dict[connection_key[1]])
                current_key = connection_key[1]
                current_table = connection_key[0]

            return list(set(query_results) & set(current_results))

    def results(self, verbose=False):
        """
        Sets results according to current filters
        """
        if verbose:
            print("Results:")
            for row in range(0, len(self.values), 3):
                result_row = self.values[row:row+3]
                if len(result_row) == 3:
                    print("%-20s %-20s %s" % \
                         (result_row[0], result_row[1], result_row[2]))
                elif len(result_row) == 2:
                    print("%-20s %-20s" % \
                         (result_row[0], result_row[1]))

                elif len(result_row) == 1:
                    print("%s" % (result_row[0]))

        return self.values

    def undo(self, verbose=False):
        """
        Undos last filter option
        """
        if len(self.history) == 1:
            if verbose:
                print("No results history.")
        else:
            current = self.history.pop()
            self.values = self.history.pop()
            if verbose:
                print("Returned last valid results.")

    def reset(self, verbose=False):
        """
        Resets created queries and filters
        for the Filter object
        """
        self.values = list(mysqldb.create_phage_id_set(self.sql_handle))
        self.history = []
        if verbose:
            print("Results and results history cleared.")

    def hits(self, verbose=False):
        """
        Returns length of current results
        """
        if verbose:
            print(f"Database hits: {len(self.values)}")
        return len(self.values)

    def group(self, field, verbose=False):
        """
        Function that creates a two-dimensional array of
        values from the results list
        separated by a characteristic.
        """
        if field in GROUP_OPTIONS_DICT[self.table]:
            return self.group_str_qualitative(field)
        else:
            if verbose:
                print(f"Grouping option by {field} is not supported.")

    def group_str_qualitative(self, field_name):
        """
        Helper function that groups qualitative
        values into a two-dimensional array
        """
        group_dicts = self.sql_handle.execute_query(
            f"SELECT DISTINCT {field_name} FROM {self.table}")
        group_set = []
        for group_dict in group_dicts:
            group_set.append(group_dict[f"{field_name}"])

        groups = {}
        for key in group_set:
            groups.update({key: []})

        for group in group_set:
            query =(f"SELECT {self.key} FROM {self.table} "
                    f"WHERE {field_name}='{group}'"+\
                    f" and {self.key} IN " + "('" + \
                     "','".join(self.values) + "')")
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                     results_list.append(result_dict[self.key])
                groups[group] = results_list

        groups.update({"None": []})

        query =(f"SELECT {self.key} FROM {self.table} "
                f"WHERE {field_name} is NULL "
                f"and {self.key} IN" + "('" + \
                 "','".join(self.values) + "')")
        results_dicts = self.sql_handle.execute_query(query)
        if results_dicts != ():
            results_list = []
            for result_dict in results_dicts:
                results_list.append(result_dict[self.key])
            groups["None"] = results_list

        return groups

    def interactive(self, sql_handle = None):
        """
        Function to start interactive filtering.
        """
        pass

class DatabaseForeignKeyTree:
    def __init__(self, sql_handle):
        self.sql_handle = sql_handle

        self.tables_dict = build_tables(sql_handle)
        self.foreign_keys_dict = build_foreign_keys(self.tables_dict.keys(),
                                                    sql_handle)
        self.primary_keys = build_primary_keys(self.tables_dict.keys(),
                                               sql_handle)


        self.root = Node(sql_handle.database)
        self.tables = []
        table_node_dictionary = {}
        for table in self.tables_dict.keys():
            table_node = self.root.create_child(table)
            table_node_dictionary.update({table: table_node})
        self.table_refs = table_node_dictionary
        for table in self.foreign_keys_dict.keys():
            for foreign_key in self.foreign_keys_dict[table]:
                foreign_key_node = (table_node_dictionary[table]).\
                                        create_child(foreign_key[0])
                table_node_dictionary[foreign_key[1]].add_child(
                                        foreign_key_node)

    def find_path(self, current_table, target_table, current_path=None):
        if current_path == None:
            current_path = [[current_table, self.primary_keys[current_table]]]

        previous_tables = []

        for previous in current_path:
            previous_tables.append(previous[0])

        if current_table == target_table:
            return current_path

        current_node = self.table_refs[current_table]
        for child_node in current_node.children:
            for parent_node in child_node.parents:
                if parent_node.id not in previous_tables and\
                                            parent_node.id != current_table:
                    path = current_path.copy()
                    path.append([parent_node.id, child_node.id])
                    path = self.find_path(parent_node.id, target_table,
                                          path)
                    if path != None:
                        if (path[-1])[0] == target_table:
                            return path

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

class Cmd_Filter(cmd.Cmd):
    def __init__(self, db_filter=None, sql_handle=None):
        super(Cmd_Filter, self).__init__()
        self.filter = db_filter
        self.sql_handle = sql_handle
        self.intro =\
        """---------------Database Search and Query---------------
        Type help or ? to list commands.\n"""
        self.prompt = "(database) (filter)user@localhost: "
        self.data = None

    def preloop(self):

        if self.filter == None:
            if self.sql_handle == None:
                self.sql_handle = \
                        mysqlconnectionhandler.MySQLConnectionHandler()
                print(\
                    "---------------------Database Login---------------------")
                self.sql_handle.database = input("MySQL database: ")
                self.sql_handle.get_credentials()
                try:
                    self.sql_handle.validate_credentials()
                except:
                    print("Unable to create a mysql connection.")
                    exit(1)
            self.filter = Filter(self.sql_handle.database, self.sql_handle)

        self.prompt = "({}) (export){}@localhost: ".\
                format(self.sql_handle.database, self.sql_handle._username)

    def do_filter(self, *unparsed_args):
        """
        Formats a MySQL database query for a given attribute
        FILTER OPTIONS: Accession, Name, Cluster
        Subcluster, Annotation_Status, Retrieve_Record,
        Annotation_QC, Annotation_Author, Notes, GeneID
        """

        filter = unparsed_args[0]

        if filter == "":
            print(""" Please select a valid attribute:
            FILTER OPTIONS: Accession, Name, Cluster
            Subcluster, Annotation_Status, Retrieve_Record,
            Annotation_QC, Annotation_Author, Notes, GeneID
            """)

        elif not (filter.lower() in self.filter.phage_attributes\
                or filter.lower() in self.filter.gene_attributes):

            print("Attribute not in database.\n ")

        else:

            filter_value = input("Enter {} value: ".format(filter))
            self.filter.add_filter(filter, filter_value)
            self.retrieve_hits()

    def retrieve_hits(self, *args):
        "Function to retrieve the hits for filtering functions"
        hits = self.filter.hits()
        if hits <= 0:
            print("        No Database Hits.")
            print("        Reloading last filtering options...\n")
            self.filter.undo()
            self.filter.retrieve_results()
        else:
            print("        Database Hits: {}\n".format(hits))

    def do_undo(self, *args):
        """Reverts back to queries generated by previous filters
        USAGE: undo
        """
        self.filter.undo()
        print("        Reloaded last filtering options")

        self.do_show_filters()

    def do_show(self,*args):
        """Displays information on current filtering
        SHOW OPTIONS: Results, Filters, Hits"""

        options = ["results", "filters", "hits"]
        option = args[0].lower()

        if option in options:
            if option == "results":
                self.show_results()
            elif option == "filters":
                self.show_filters()
            elif option == "hits":
                self.retrieve_hits()

        else:
            print("""Please input a valid option
            SHOW OPTIONS: Results, Filter,Hits""")

    def show_results(self):
        """Shows results for current database filtering
        """

        print("        Results:\n")
        for row in range(0, len(self.filter.results), 3):
            result_row = self.filter.results[row:row+3]
            if len(result_row) == 3:
                print("%-20s %-20s %s" % \
                     (result_row[0], result_row[1], result_row[2]))
            elif len(result_row) == 2:
                print("%-20s %-20s" % \
                     (result_row[0], result_row[1]))

            elif len(result_row) == 1:
                print("%s" % (result_row[0]))

        print("\n")

    def show_filters(self, *args):
        """Displays current fitlers
        """
        if self.filter.phage_filters:
            print("        Phage table filters:")
            for phage_filter in self.filter.phage_filters.values():
                print("    " + phage_filter)
            print("\n")

        if self.filter.gene_filters:
            print("        Gene table filters:")
            for gene_filter in self.filter.gene_filters.values():
                print("    " + gene_filter)
            print("\n")

        if not self.filter.phage_filters and not self.filter.gene_filters:
            print("        No current filters applied.")

    def do_reset(self, *args):
        """Resets results history and current filters
        USAGE: reset
        """
        print("        Resetting Filters and Results History...\n")
        self.filter.reset()

    def do_clear(self, *args):
        """Clears display terminal
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')
        print(self.intro)

    def do_return(self, *args):
        """Finish quering and return results
        USAGE: return
        """

        print("        Exiting Filtering...\n")

        return True

    def do_exit(self, *args):
        """Exits program entirely without returning values
        USAGE: exit
        """

        sys.exit(1)

    def postloop(self):

        self.data = self.filter.results

if __name__ == "__main__":
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle._username = "root"
    sql_handle._password = "phage"
    sql_handle.database = "Actino_Draft"

    #db_filter = Filter(sql_handle, table="phage")
    #db_filter.add_filter("phage", "hoststrain", "Gordonia", verbose=True)
    #db_filter.add_filter("gene", "notes", "antirepressor", verbose=True)
    #db_filter.update()
    #db_filter.sort("PhageID")
    #db_filter.results(verbose=True)

    #db_filter = Filter(sql_handle, table="gene")
    #db_filter.add_filter("gene", "notes", "antirepressor",verbose=True)
    #db_filter.add_filter("phage", "hoststrain", "gordonia", verbose=True)
    #db_filter.update()
    #db_filter.sort("GeneID", verbose=True)
    #db_filter.results(verbose=True)
