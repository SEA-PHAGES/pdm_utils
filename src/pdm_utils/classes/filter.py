"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome, cds, mysqlconnectionhandler, databasetree
from typing import List
import cmd, readline, argparse, os, sys, re, string, math

#Global file constants
OPERATOR_OPTIONS      = ["=", "!=", ">", "<"]
COMPARATIVE_OPERATORS = [">", "<"]
COMPARABLE_TYPES      = ["int", "decimal",
                         "mediumint", "float",
                         "datetime", "double"]
GROUP_OPTIONS = ["limited_set", "num_set", "str_set"]

class Filter:
    def __init__(self, sql_handle, table='phage'):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.sql_handle = sql_handle

        self.db_tree = databasetree.DatabaseTree(sql_handle)

        if table not in self.db_tree.show_tables():
            print(f"Table passed to filter is not in {sql_handle.database}")
            raise ValueError

        self.values = None

        self.table = table

        table_node = self.db_tree.get_table(table)
        self.key = table_node.primary_key.id

        self.history = []
        self.filters = {}

        self.values_valid = True
        self.updated = True

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
            raise ValueError

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

    def switch_table(self, raw_table, verbose=False):
        table = self.translate_table(raw_table)
        table_node = self.db_tree.get_table(table)
        
        self.values = self.transpose(self.table, self.values, 
                                     target_table=table)
        self.table = table
        self.key = table_node.primary_key.id
        self.update()
        
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

        self.updated = False

    def set_values(self, values):
        if values:
            if not isinstance(values[0], str):
                raise ValueError

        self.values = values
        self.valid_values = False
        self.updated = False

    def refresh(self):
        """
        Refreshes values
        """
        if not self.values or self.values_valid:
            self.values_valid = True
            return

        self.values = self.db_tree.build_values(self.table, self.key, 
                                                values=self.values)
        self.values_valid = True

    def update(self, verbose=False):
        """
        Updates results list for the Filter object
        """
        self.refresh()
        if self.updated:
            return

        if not self.filters:
            self.values = self.db_tree.build_values(self.table, self.key,
                                                    values=self.values)
            self.updated = True
            return

        query_values = self.values 
       
        for table in self.filters.keys():
            queries = self.build_queries(table)
            if verbose:
                print(f"Filtering {self.sql_handle.database}.{table} for "
                       + " and ".join(queries) + "...")

            if table == self.table:
                query_values = self.db_tree.build_values(
                                                table, self.key,
                                                queries=queries,
                                                values=query_values)

            else:
                key = self.db_tree.get_table(table).primary_key.id
                untransposed_values = self.db_tree.build_values(
                                                table, key,
                                                queries=queries)
                transposed_values = self.transpose(table, untransposed_values)
                if query_values:
                    query_values = list(set(query_values) &\
                                        set(transposed_values))
                else:
                    query_values = transposed_values
            if not query_values:
                self.values = query_values
                return

        self.values = query_values
        self.history.append(query_values)

    def build_queries(self, table):
        queries = []
        for field in self.filters[table].keys():
            for operator in self.filters[table][field].keys():
                operator_queries = []
                for value in self.filters[table][field][operator]:
                    operator_queries.append(f"{field}{operator}'{value}'")
                query = ("(" + " or ".join(operator_queries) + ")")
                queries.append(query)

        return queries

    def sort(self, sort_field, verbose=False):
        if not self.values:
            return

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

    def transpose(self, table, values, target_table=None):
        if not values:
            return []
        if target_table == None:
            target_table = self.table
        table_node = self.db_tree.get_table(table)
        target_table_node = self.db_tree.get_table(target_table)
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
       
        if current_key != target_table_node.primary_key.id:
            current_values = self.db_tree.build_values(
                                    current_table, 
                                    target_table_node.primary_key.id,
                                    values_column=current_key,
                                    values=current_values)

        return current_values

    def results(self, verbose=False):
        """
        Sets results according to current filters
        """
        if not self.values:
            return []

        if not self.values_valid:
            self.refresh()

        if verbose:
            print(" " + "_"*57 + " ")
            print("| %-20s" % \
                 (f"{self.key} Results:") \
                 + " "*36 + "|")
            print("|" + "-"*57 + "|")
            for row in range(0, len(self.values), 3):
                result_row = self.values[row:row+3]
                if len(result_row) == 3:
                    print("| %-20s %-20s %-13s |" % \
                         (result_row[0], result_row[1], result_row[2]))
                elif len(result_row) == 2:
                    print("| %-20s %-20s %-13s" % \
                         (result_row[0], result_row[1], "") \
                         + " " + "|")

                elif len(result_row) == 1:
                    print("| %-20s %-20s %-13s" % \
                         (result_row[0], "", "") \
                         + " " + "|")
            print("|" + "_"*57 + "|")

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
        self.values = None
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
        table = self.translate_table(table)
        table_node = self.db_tree.get_table(table)
        field = self.translate_field(field, table)
        field_node = table_node.get_column(field) 

        if field_node.group in GROUP_OPTIONS:
            if verbose:
                print(f"Grouping by {field} in {table}...")

            if field_node.group == "limited_set":
                distinct_field = field

            elif field_node.group == "num_set":
                distinct_field = self.large_num_set_distinct_query(
                                                             table, field)
            elif field_node.group == "str_set":
                distinct_field = f"substring({field}, 1, 1)"
            return self.build_groups(table_node, field_node,
                                     distinct_field=distinct_field,
                                     verbose=verbose)
        else:
            print(f"Grouping option by {field} is not supported.")
            raise ValueError
       
    def build_groups(self, table_node, field_node, distinct_field=None,
                     verbose=False):
        if not self.values:
            return {}

        if distinct_field == None:
            distinct_field = field_node.id
        query = f"SELECT DISTINCT {distinct_field} FROM {table_node.id}"

        if table_node.id != self.table:
            group_assist_filter = self.copy()
            group_assist_filter.switch_table(table_node.id)
        
        groups = {}
        for dict in self.sql_handle.execute_query(query):
            query = f"{distinct_field}='{str(dict[distinct_field])}'"

            if table_node.id == self.table:
                values = self.db_tree.build_values(table_node.id, 
                                                   table_node.primary_key.id,
                                                   queries=[query],
                                                   values=self.values)

            else:
                values = self.group_transpose(group_assist_filter,
                                              table_node, query)
           
            if values:
                groups[str((dict[distinct_field]))] = values
                if verbose:
                    print(f"...Group found for {field_node.id}="
                          f"'{str(dict[distinct_field])}' "
                          f"in {table_node.id}...")
                    print("......Database hits in group: " + str(len(values)))
            else:
                if verbose:
                    print(f"...No group found for {field_node.id}="
                          f"'{str(dict[distinct_field])}' "
                          f"in {table_node.id}...")
 
        if field_node.null:
            null_query = f"{field_node.id} IS NULL"
            values = self.db_tree.build_values(table_node.id, field_node.id,
                                               values_column=distinct_field,
                                               queries=[null_query],
                                               values=self.values)
            if values:
                groups["Null"] = values

        return groups

    def group_transpose(self, assist_filter, table_node, query):
        values = self.db_tree.build_values(table_node.id,
                                           table_node.primary_key.id,
                                           queries=[query],
                                           values=assist_filter.values)

        values = self.transpose(table_node.id, values)
        if values:
            values = self.db_tree.build_values(self.table, self.key, 
                                               values=values)

        return values
        
    def large_num_set_distinct_query(self, table, field):
        range_query = (
                f"SELECT round(log10(Max({field}) - Min({field}))) as pow "
                f"FROM {table}")
        range_pow = int(self.sql_handle.execute_query(range_query)[0]["pow"])
        range_pow = 10**(range_pow-2)
        return f"round({field}/{range_pow})*{range_pow}"

    def copy(self):
        duplicate_filter = Filter(self.sql_handle, table=self.table)
        duplicate_filter.db_tree = self.db_tree
        if self.values != None:
            duplicate_filter.values = self.values.copy()

        duplicate_filter.key = self.key

        duplicate_filter.values_valid = self.values_valid
        duplicate_filter.updated = self.updated

        duplicate_filter.filters = self.copy_filters()

        return duplicate_filter

    def copy_filters(self):
        duplicate_filters = {}
        for table in self.filters.keys():
            duplicate_filters[table] = {}
            for field in self.filters[table].keys():
                duplicate_filters[table].update({field: {}})
                for operator in self.filters[table][field].keys():
                    values = []
                    for value in self.filters[table][field][operator]:
                        values.append(value)
                    duplicate_filters[table][field].update({operator : values})

        return duplicate_filters


def debug_filter(db_filter):
    db_filter.add_filter("phage", "hoststrain", "Gordonia", "=", verbose=True)
    db_filter.add_filter("gene", "notes", "antirepressor", "=", verbose=True)
    db_filter.add_filter("phage", "sequencelength", "50000", ">", verbose=True)
    db_filter.update(verbose=True)
    db_filter.sort("PhageID")
    db_filter.results(verbose=True)
    db_filter.hits(verbose=True)
    print("")
    db_filter.group("phage", "Subcluster2", verbose=True)
    db_filter.group("phage", "PhageID", verbose=True)
    db_filter.group("gene", "Orientation", verbose=True)

if __name__ == "__main__":
    #sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    #sql_handle.database = input("Please enter database name: ")
    #sql_handle.get_credentials()
    #sql_handle.validate_credentials()
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = "Actino_Draft"
    sql_handle.username = "root"
    sql_handle.password = "phage"
    sql_handle.validate_credentials()
    db_filter = Filter(sql_handle) 
    debug_filter(db_filter)

    
