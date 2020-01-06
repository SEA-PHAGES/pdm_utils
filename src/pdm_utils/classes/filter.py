"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""

import cmd
import readline
import os
import sys
import re
import string
import math
import time
import csv
from pathlib import Path
from pdm_utils.classes import mysqlconnectionhandler, databasetree
from pdm_utils.pipelines import export_db

#Global file constants
OPERATOR_OPTIONS      = ["=", "!=", ">", "<"]
COMPARATIVE_OPERATORS = [">", "<"]
COMPARABLE_TYPES      = ["int", "decimal",
                         "mediumint", "float",
                         "datetime", "double"]
GROUP_OPTIONS = ["limited_set", "num_set", "str_set"]

class Filter:
    "MySQL database filtering object."
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

        self.history = HistoryNode("head", None)
        self.history_count = 0

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

    def add_history(self, function):
        if self.history_count == 100:
            current = self.history
            for i in range(0, 50):
                current = current.get_next()
            
            current.next = None
            self.history_count = 50

        if function == "switch":
            self.history.create_next(function, 
                                 [[self.table, self.key, self.copy_values()], \
                                  self.values_valid, self.updated])
        elif function == "add":
            self.history.create_next(function, 
                                  [self.copy_filters(), \
                                   self.values_valid, self.updated])

        elif function in ["set", "update", "sort"]:
            self.history.create_next(function, 
                                  [self.copy_values(), \
                                   self.values_valid, self.updated])

        elif function == "reset":
            self.history.create_next(function, 
                                 [[self.copy_values(), self.copy_filters()], \
                                   self.values_valid, self.updated])

        self.history_count += 1

    def undo(self, verbose=False):
        """
        Undos last filter option 
        """
        if self.history_count == 0:
            if verbose:
                print("No results history.")
        else:
            current = self.history.remove_next()
            if current.id in ["set", "update", "sort"]:
                self.values = current.history[0]

            elif current.id == "add":
                self.filters = current.history[0]

            elif current.id == "reset":
                self.values = current.history[0][0]
                self.filters = current.history[0][1]

            elif current.id == "switch":
                self.table = current.history[0][0]
                self.key = current.history[0][1]
                self.values = current.history[0][2]

            else:
                pass

            self.values_valid = current.history[1]
            self.updated = current.history[2]

            self.history_count -= 1

            if verbose:
                print(f"Undid {current.id} operation.")

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

    def switch_table(self, raw_table, verbose=False): 
        self.add_history("switch")

        table = self.translate_table(raw_table)
        table_node = self.db_tree.get_table(table)
        
        self.values = self.transpose(self.table, self.values, 
                                     target_table=table)
        if not self.values:
            self.values = None 

        self.table = table
        self.key = table_node.primary_key.id
        self.update()
        
    def add_filter(self, raw_table, raw_field, operator, value, verbose=False):
        """
        Adds a filter to database queries
        """
        self.add_history("add")

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

        self.add_history("set")

        self.values = values
        self.values_valid = False
        self.updated = False
 
    def update(self, verbose=False):
        """
        Updates results list for the Filter object
        """
        self.add_history("update")

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

    def sort(self, sort_field, verbose=False): 
        if not self.values:
            return

        self.add_history("sort")

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

    def reset(self, verbose=False):
        """
        Resets created queries and filters
        for the Filter object
        """ 
        self.add_history("reset")

        self.values = None
        self.filters = {}
        self.values_valid = True
        self.updated = True

        if verbose:
            print("Results and filters cleared.")

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
   
    def hits(self, verbose=False):
        """
        Returns length of current results
        """
        if self.values == None:
            if verbose:
                print("Database results currently empty.")
            return 0

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
            groups = self.build_groups(table_node, field_node,
                                     distinct_field=distinct_field,
                                     verbose=verbose)
            if verbose:
                for group in groups.keys():
                    values = groups[group]
                    print(" " + "_"*57 + " ")
                    print("| %-40s" % \
                        (f"{self.key} Group {group} Results:") \
                         + " "*16 + "|")
                    print("|" + "-"*57 + "|")
                    for row in range(0, len(values), 3):
                        result_row = values[row:row+3]
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
                    print(f"Database Hits: {len(values)}")

            return groups
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
                    #print("......Database hits in group: " + str(len(values)))
            #else:
                #if verbose:
                #    print(f"...No group found for {field_node.id}="
                #          f"'{str(dict[distinct_field])}' "
                #          f"in {table_node.id}...")
 
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
        if self.values:
            values = list(set(values) & set(self.values))

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
        duplicate_filter.values = self.copy_values()
        
        duplicate_filter.key = self.key

        duplicate_filter.values_valid = self.values_valid
        duplicate_filter.updated = self.updated

        duplicate_filter.filters = self.copy_filters()

        return duplicate_filter

    def copy_values(self):
        values = None
        if self.values != None:
            values = self.values.copy()

        return values

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

    def write_csv(self, output_path, csv_name=None, verbose=False):
        if csv_name == None:
            date = time.strftime("%Y%m%d")
            csv_name = f"{date}_filter"

        csv_path = output_path.joinpath(f"{csv_name}.csv")
        csv_version = 1

        while(csv_path.exists()):
            csv_version += 1
            csv_path = output_path.joinpath(f"{csv_name}{csv_version}.csv")
    
        if self.values == None:
            if verbose: 
                print("Database results currently empty.")
            return

        if verbose:
                print(f"...Writing {csv_path.name}.csv...")

        csv_path.touch()
        with open(csv_path, 'w', newline="") as csv_file:
            csvwriter=csv.writer(csv_file, delimiter=",",
                                 quotechar="|",
                                 quoting=csv.QUOTE_MINIMAL)
            for value in self.values:
                csvwriter.writerow([value])

def parse_filters(unparsed_filters):
    """
    Helper function to return a two-dimensional
    array of filter parameters
    """
    filter_format = re.compile("\w+\.\w+[=<>!]+\w+", re.IGNORECASE)
    filters = []
    for filter in unparsed_filters:
        if re.match(filter_format, filter) != None:
            filters.append(re.split("\W+", filter) +\
                           re.findall("[!=<>]+", filter))
        else:
            raise ValueError(f"Unsupported filtering format: '{filter}'")
                
    return filters

def parse_groups(unparsed_groups):
    group_format = re.compile("\w+\.\w+", re.IGNORECASE)
    groups = []
    for group in unparsed_groups:
        if re.match(group_format, group) != None:
            groups.append(re.split("\W+", group))
        else:
            raise ValueError(f"Unsupported grouping format: '{group}'")

    return groups

class HistoryNode:
    "Filter history storage object."
    def __init__(self, id, history):
        if not isinstance(id, str):
            raise TypeError("id parameter must be a string object.")
        self.id = id

        self.history = history

        self.next = None

    def get_id(self):
        return self.id

    def get_history(self):
        return self.history
    
    def has_next(self):
        has_next = (self.next!=None)
        return has_next
    
    def get_next(self):
        return self.next

    def add_next(self, history_node):
        if not isinstance(history_node, HistoryNode):
            raise TypeError("history_node parameter must be a string object.")
        if not self.has_next():
            self.next = history_node 
        else:
            history_node.add_next(self.next)
            self.next = history_node

    def create_next(self, id, history):
        if not isinstance(id, str):
            raise TypeError("id parameter must be a string object.")
        next_node = HistoryNode(id, history)
        self.add_next(next_node)

        return next_node

    def remove_next(self):
        next_node = None
        if self.has_next():
            next_node = self.get_next()
            self.next=None
            if next_node.has_next():
                new_next = next_node.get_next()
                self.add_next(new_next) 

        return next_node

class CmdFilter(cmd.Cmd):
    "Filtering CmdLoop object."
    def __init__(self, sql_handle):
        super(CmdFilter, self).__init__()

        self.sql_handle = sql_handle
        db_filter = Filter(sql_handle)
        self.db_filter = db_filter

        self.intro =\
       f"""---------------Filtering in {sql_handle.database}---------------
        Type help or ? to list commands.\n"""
        self.prompt = (f"({sql_handle.database}) "
                       f"({db_filter.table})"
                       f"{sql_handle.username}@localhost: ")

    def preloop(self):
        self.prompt = (f"({self.sql_handle.database}) "
                       f"({self.db_filter.table})"
                       f"{self.sql_handle.username}@localhost: ")

    def do_add(self, *args):
        """
        Adds a filter to the current filtering object.
        USAGE: add {table}.{field}{(=/!=/>/<)}{value}"""

        try:
            filters = parse_filters(args)
            for filter in filters:
                self.db_filter.add_filter(filter[0], filter[1],
                                           filter[3], filter[2],
                                           verbose=True)
        except ValueError:
            print("Filter not accepted.") 

    def do_update(self, *args):
        """
        Updates filtering object results list with current filters.
        USAGE: update"""
        self.db_filter.update(verbose=True)

    def do_group(self, *args):
        """
        Groups current filtering object results.
        USAGE: group {table}.{field}"""
        try:
            group = parse_groups([args[0]])
            self.db_filter.group(group[0][0], group[0][1], verbose=True)

        except ValueError:
            print(f"Unsupported grouping format: '{args[0]}'")

    def do_results(self, *args):
        """
        Shows current filter results.
        USAGE: results"""
        self.db_filter.results(verbose=True)

    def do_hits(self, *args):
        """
        Displays the number of filtering object results.
        USAGE: hits"""
        self.db_filter.hits(verbose=True)

    def do_switch(self, *args):
        """
        Switches the current table of the filtering object.
        USAGE: switch {table}"""
        self.db_filter.switch_table(args[0])

    def do_reset(self, *args): 
        """
        Resets filtering object results, filters, and history.
        USAGE: reset"""
        self.db_filter.reset()

    def do_dump(self, *args):
        """
        Dumps results in a csv file at a given path.
        USAGE: dump path/to/csv_folder"""
        if args[0] != "":
            output_path = Path(args[0])
        else:
            output_path = Path.cwd()

        if output_path.is_dir():
            try:
                self.db_filter.write_csv(output_path, verbose=True) 

            except:
                print(f"Csv-file writing in {output_path} failed.")

        else:
            print(f"{args[0]} is not a valid folder path.")

    def do_set(self, *args):
        """
        Retrieves and sets results in a csv file at a given path.
        USAGE: set path/to/existing/csv_file.csv"""
        import_path = Path(args[0])
        import_path = import_path.expanduser()
        import_path = import_path.resolve()

        if import_path.is_file():
            try:
                values = export_db.parse_value_list_input(import_path)
                self.db_filter.set_values(values)
                self.db_filter.refresh()

            except:
                print(f"Csv-file reading of {args[0]} failed.")
        else:
            print(f"{args[0]} is not a valid csv-file path.")

    def do_show(self, *args):
        """
        Prints current characteristics of a table(s). 
        Leave empty to print the entire database.
        USAGE: show {table}"""
        
        if args[0] == "":
            self.db_filter.db_tree.print_info()
            return
        
        try:
            table = self.db_filter.translate_table(args[0]) 
            table_node = self.db_filter.db_tree.get_table(table)
            table_node.print_columns_info()

        except ValueError:
            pass

    def do_undo(self, *args):
        """
        Undos last status change.
        USAGE undo
        """
        self.db_filter.undo(verbose=True)

    def do_clear(self, *args):
        """
        Clears display terminal.
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')
        print(self.intro)

    def do_exit(self, *args):
        """
        Exits program entirely.
        USAGE: exit
        """
        print("       Exiting...\n")

        sys.exit(1)

if __name__ == "__main__":
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = input("Please enter database name: ")
    sql_handle.get_credentials()
    sql_handle.validate_credentials()
    cmd_filter = CmdFilter(sql_handle)
    cmd_filter.cmdloop()

    
