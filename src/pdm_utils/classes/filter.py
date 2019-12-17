"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome, cds, mysqlconnectionhandler, databasetree
from typing import List
import cmd, readline, argparse, os, sys, re, string

#Global file constants
OPERATOR_OPTIONS      = ["=", "!=", ">", "<"]
COMPARATIVE_OPERATORS = [">", "<"]
COMPARABLE_TYPES      = ["int", "decimal",
                         "mediumint", "float",
                         "datetime", "double"]

class Filter:
    def __init__(self, sql_handle, values_list=None, table='phage'):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.sql_handle = sql_handle

        self.db_tree = databasetree.DatabaseTree(sql_handle)

        if table not in self.db_tree.show_tables():
            print(f"Table passed to filter is not in {sql_handle.database}")
            raise ValueError

        self.values = values_list

        self.table = table

        table_node = self.db_tree.get_table(table)
        self.key = table_node.primary_key.id

        self.history = []
        self.filters = {}

        self.values_valid = False
        self.updated = False

    def translate_table(self, raw_table, verbose=False):
        for table in self.db_tree.db_node.show_tables():
            if table.lower() == raw_table.lower():
                return table

        print(f"Table '{raw_table}' requested to be filtered "
              f"is not in '{self.sql_handle.database}'")
        raise ValueError

    def translate_field(self, raw_field, table, verbose=False):
        table_node = self.db_tree.get_table(table)
        if table_node == None:
            print(
              f"Table '{table}' requested to be filtered "
              f"is not in '{self.sql_handle.database}'")
            raise ValueError

        for field in table_node.show_columns():
            if field.lower() == raw_field.lower():
                return field
        
        print(f"Field '{raw_field}' requested to be filtered"
                  f" is not in '{table_node.id}'")
        raise ValueError

    def check_operator(self, operator, table, field, verbose=False):
        if operator not in OPERATOR_OPTIONS:
            return None

        table_node = self.db_tree.get_table(table) 
        if table_node == None:
            print(f"Table '{table}' requested to be filtered "
              f"is not in '{self.sql_handle.database}'")
            raise ValueError

        column_node = table_node.get_column(field)
        if column_node == None:
            print(f"Field '{field}' requested to be filtered "
                  f"is not in '{table_node.id}'")
            raise ValueError
        
        type = column_node.parse_type()
        if operator in COMPARATIVE_OPERATORS and \
                type not in COMPARABLE_TYPES:
            print(f"Field ' {field}' requested to be filtered "
                  f"is not comparable (Operator: {operator})")
            raise ValueError

    def add_filter(self, raw_table, raw_field, value, operator, verbose=False):
        """
        Adds a filter to database queries
        """
        table = self.translate_table(raw_table)
        table_node = self.db_tree.get_table(table)

        field = self.translate_field(raw_field, table)
        self.check_operator(operator, table, field)


        if field == table_node.show_primary_key():
            print(f"Primary key {field} in {table} "
                   "cannot be a field requested for filtering")
            raise ValueError

        elif table.lower() in self.filters.keys():
            if field in (self.filters[table]).keys():
                if operator in (self.filters[table])[field].keys():
                    (self.filters[table])[field][operator].append(value)
                else:
                    (self.filters[table])[field][operator] =  [value]
            else:
                (self.filters[table])[field] = {operator : [value]}
        else:
            self.filters[table] = {field: {operator : [value]}}

    def refresh(self):
        """
        Refreshes values
        """
        if self.values == None:
            self.values_valid = True

        else:
            query = (f"SELECT {self.key} "
                     f"FROM {self.table} WHERE {self.key} IN "
                     "('" + "','".join(self.values) + "')")
            
            results_dicts = self.sql_handle.execute_query(query)
            results = []
            for result in results_dicts:
                results.append(result[self.key])

            self.values = results
            self.values_valid = True

    def update(self, verbose=False):
        """
        Updates results list for the Filter object
        """
        self.refresh()
        
        if not self.filters:
            self.updated = True
            return
        if self.values != None and len(self.values) == 0:
            self.updated = True
            return
       
        if self.values == None:
            query_values = []
        else:
            query_values = self.values 
       
        for table in self.filters.keys():
            queries = self.build_queries(table)
            if verbose:
                print(f"Filtering {table} in "
                      f"{self.sql_handle.database} for "+\
                            " and ".join(queries) + "...")

            if table == self.table:
                query_values = self.db_tree.build_values(
                                                table, self.key,
                                                queries=queries,
                                                values=query_values)
            else:
                key = self.db_tree.get_table(table).primary_key.id
                untransposed = self.db_tree.build_values(
                                                table, key,
                                                queries=queries)
                query_values = list(set(query_values) &\
                                    set(self.transpose(table, untransposed)))

            if not query_values:
                break

        self.values = query_values
        if len(query_values) > 0:
            self.history.append(query_values)

    def build_queries(self, table):
        queries = []
        for field in self.filters[table].keys():
            for operator in self.filters[table][field].keys():
                queries.append(
                    f"{field}{operator}'"+ \
                        "','".join((self.filters[table])\
                                                [field]\
                                                [operator]) + "'")
        return queries

    def sort(self, sort_field, verbose=False):
        sort_field = self.translate_field(sort_field, self.table) 
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

    def transpose(self, table, values, verbose=False):
        table_node = self.db_tree.get_table(table)
        target_table_node = self.db_tree.get_table(self.table)
        traversal_path = self.db_tree.find_path(table_node, target_table_node)

        if not traversal_path:
            print(f"Table {table} is not connected to {self.table} "
                       "by foreign keys")
            raise ValueError 

        current_key = table_node.primary_key.id
        current_table = table
        current_values = values

        for connection_key in traversal_path:
            if current_key != connection_key[1]:
                current_values = self.db_tree.build_values(
                                    current_table, connection_key[1],
                                    values_column=current_key,
                                    values=current_values)
            current_key = connection_key[1]
            current_table = connection_key[0] 
       
        if current_key != self.key:
            current_values = self.db_tree.build_values(
                                    current_table, self.key,
                                    values_column=current_key,
                                    values=current_values)

        return current_values

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
        self.values = list(phamerator.create_phage_id_set(self.sql_handle))
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

    def group(self, table, field, verbose=False):
        """
        Function that creates a two-dimensional array of
        values from the results list
        separated by a characteristic.
        """
        table = self.translate_table(raw_table)
        table_node = self.db_tree.get_table(table)
        field = self.translate_field(raw_field, table)
        field_node = table_node.get_column(field) 

        group_options = {"limited_set"    : self.group_limited_set,
                         "large_num_set"  : self.group_large_num_set,
                         "large_str_set"  : self.group_large_str_set}

        if table_node.group in group_options.keys():
            return group_options[table_node.group]
        else:
            print(f"Grouping option by {field} is not supported.")
            raise ValueError
        
    def group_limited_set(self, table_node, field_node):
        """
        Helper function that groups limited set
        values into a two-dimensional array
        """
        group_dicts = self.sql_handle.execute_query(
            f"SELECT DISTINCT {field_node.id} FROM {table_node.id}")
        group_set = []
        for group_dict in group_dicts:
            group_set.append(group_dict[f"{field_node.id}"])

        groups = {}
        for key in group_set:
            groups.update({key: []})

        for group in group_set:
            query =(f"SELECT {table_node.primary_key.id} FROM {table_node.id} "
                    f"WHERE {field_node.id}='{group}'"+\
                    f" and {self.key} IN " + "('" + \
                     "','".join(self.values) + "')")
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                     results_list.append(result_dict[self.key])
                groups[group] = results_list

        if field_node.null:
            groups.update({"None": []})

            query =(f"SELECT {table_node.primary_key.id} FROM {table_node.id} "
                    f"WHERE {field_node.id} is NULL "
                    f"and {self.key} IN" + "('" + \
                     "','".join(self.values) + "')")
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                    results_list.append(result_dict[self.key])
                groups["None"] = results_list
    
        return groups
   
    def group_large_num_set(self, table, field):
        pass

    def group_large_str_set(self, table, field):
        pass

def debug_filter(db_filter):
    db_filter.add_filter("phage", "hoststrain", "Gordonia", "=", verbose=True)
    db_filter.add_filter("gene", "notes", "antirepressor", "=", verbose=True)
    db_filter.add_filter("phage", "sequencelength", "50000", ">", verbose=True)
    db_filter.update(verbose=True)
    db_filter.sort("PhageID")
    db_filter.results(verbose=True)
    db_filter.hits(verbose=True)

if __name__ == "__main__":
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = input("Please enter database name: ")
    sql_handle.get_credentials()
    sql_handle.validate_credentials()
    db_filter = Filter(sql_handle) 
    #debug_filter(db_filter)

    
