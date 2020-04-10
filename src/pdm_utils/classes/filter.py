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
from collections import OrderedDict
from networkx import Graph
from pathlib import Path
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import parsing
from pdm_utils.functions import querying as q
from sqlalchemy import Column
from sqlalchemy.engine.base import Engine

class Filter:
    def __init__(self, alchemist=None, key=None):
        self._engine=None
        self._graph=None

        self._updated = True 
        self._values_valid = True
        self._connected = False 

        if isinstance(alchemist, AlchemyHandler):
            if not alchemist.connected:
                alchemist.connect(ask_database=True)

            if alchemist.graph == None:
                alchemist.build_graph()
                
            self._engine = alchemist.engine
            self._graph = alchemist.graph

            self._connected = True
       
        if isinstance(key, Column):
            self._key = key

        else:
            self._key = None

        self._values = []
        self._filters = {}

        self.verbose = False

    @property
    def updated(self):
        updated = self._updated
        return updated

    @property
    def values_valid(self):
        values_valid = self._values_valid
        return values_valid
    
    @property
    def connected(self):
        connected = self._connected
        return connected

    @property
    def filters(self):
        filters = self.copy_filters()
        return filters

    @property
    def engine(self):
        engine = self._engine
        return engine

    @engine.setter
    def engine(self, engine):
        if not isinstance(engine, Engine):
            raise TypeError
        
        self._engine = engine
        self._connected = False

    @property
    def graph(self):
        graph = self._graph
        return graph

    @graph.setter
    def graph(self, graph):
        if not isinstance(graph, Graph):
            raise TypeError

        self._graph = graph
        self._connected = False

    @property 
    def values(self):
        if self._values != None:
            values = self._values.copy()  
            return values
        else:
            return []
 
    @values.setter
    def values(self, values):
        if not isinstance(values, list):
            raise TypeError("Filter object values must be of type list.")
        self._values = values
        self._values_valid = False

    @property
    def key(self):
        key = self._key
        return key 

    @key.setter
    def key(self, key):
        if isinstance(key, Column):
            self._key = key
        elif isinstance(key, str):
            self._key = q.get_column(self.graph.graph["metadata"], key)
        else:
            raise TypeError("Filter key value must be of type SqlAlchemy Column.")
        
    def connect(self):
        """Connect Filter object to a database via an AlchemyHandler.
        """
        if self._connected:
            return

        alchemist = AlchemyHandler()
        alchemist.connect(ask_database=True)

        alchemist.build_graph()
        
        self._engine = alchemist.engine
        self._graph = alchemist.graph

        self._connected = True

    def check(self):
        """Check Filter object contains valid essential objects.
        Filter object requires a SQLAlchemy Engine and Column as well as
        a NetworkX Graph.
        """
        if not self._connected:
            self.connect()
        else: 
            if not isinstance(self._engine, Engine): 
                raise AttributeError("Filter object is missing valid engine.")
            if not isinstance(self._graph, Graph):
                raise AttributeError("Filter object is missing valid graph.")

        if not isinstance(self._key, Column):
            raise AttributeError("Filter object is missing valid column key.")

    def add(self, filter):
        """Add a filter to the Filter object class.

        :param filter: Formatted MySQL WHERE clause.
        :type filter: str
        """
        where_clause = q.build_where_clause(self.graph, filter)
        parsed_filter = parsing.parse_filter(filter)
        filter_left = parsed_filter[0] + "." + parsed_filter[1]\
                    + parsed_filter[2]
      
        #Stores the BinaryExpression with a key of the Column/Operator pairing.
        if filter_left not in self._filters.keys():
            self._filters.update({filter_left : [where_clause]}) 
           
        #Appends the BinaryExpression to the existing Column/Operator key.
        else:  
            self._filters[filter_left].append(where_clause)

        self._updated = False

    def remove(self, filter):
        """Remove a filter from the filter object class.
        :param filter: Formatted MySQL WHERE clause.
        :type filter: str
        """
        parsed_filter = parsing.parse_filter(filter)
        filter_left = parsed_filter[0] + "." + parsed_filter[1]\
                    + parsed_filter[2]

        #Indexes into the dictionary using the Column/Operator pairing.
        if filter_left in self._filters.keys():
            filters = self._filters[filter_left]

            #Matches the string value to the BinaryExpressions stored value
            for clause in filters:
                if clause.right.value == parsed_filter[3]:
                    filters.remove(clause)

            if len(self._filters[filter_left]) == 0:
                self._filters.pop(filter_left)

        self._updated = False

    def convert_column_input(self, raw_column):
        """Converts a column input, string or Column, to a Column.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        """
        self.check()

        if isinstance(raw_column, str):
            column = q.get_column(self.graph.graph["metadata"], raw_column)
        elif isinstance(raw_column, Column):
            column = raw_column
        else:
            raise TypeError("Column must be either a string or a Column object")



        return column

    def build_where_clauses(self):
        """Builds BinaryExpression objects from stored Filter object filters.
        """
        where_clauses = []
        for filter in self._filters.keys():
            where_clauses = where_clauses + self._filters[filter]

        return where_clauses

    def build_values(self, where=None, order_by=None, column=None, limit=8000):
        """Queries for values from stored WHERE clauses and Filter key.
        
        :param where: MySQL WHERE clause_related SQLAlchemy object(s).
        :type where: BinaryExpression
        :type where: list
        :param order_by: MySQL ORDER BY clause-related SQLAlchemy object(s).
        :type order_by: Column
        :type order_by: list
        :param column: SQLAlchemy Column object or object name.
        :type column: Column
        :type column: str
        :param limit: SQLAlchemy IN clause query length limiter.
        :type limit: int
        :returns: Distinct values fetched from given and innate constraints.
        :rtype: list
        """
        self.check()

        if column is None:
            select_column = self._key
        else:
            select_column = self.convert_column_input(column)
 
        if not where is None:
            if isinstance(where, list):
                base_clauses = where
            else:
                base_clauses = [where]
        else:
            base_clauses = []

        query = q.build_distinct(self.graph, select_column, where=base_clauses, 
                                                            order_by=order_by,
                                                            add_in=self._key)

        values = []
        if self._values != []:
            values = q.first_column_value_subqueries(
                        self._engine, query, self._key, self._values, 
                                                        limit=limit)

        else:
            proxy = self.engine.execute(query)
            results = proxy.fetchall()

            for result in results:
                values.append(result[0])

        if self._key.type.python_type == bytes:
            parsing.convert_from_encoded(values)

        return values
  
    def transpose(self, raw_column, return_dict=False, set_values=False):
        """Queries for distinct values from stored values and a MySQL Column.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        :param return_dict: Toggle whether to return data as a dictionary.
        :type return_dict: Boolean
        :param set_values: Toggle whether to replace Filter key and values.
        :type set_values: Boolean
        :returns: Distinct values fetched from given and innate constraints.
        :rtype: list
        :rtype: dict
        """
        self.check()

        if not self._values:
            return []
    
        column = self.convert_column_input(raw_column)
        name = column.name
   
        values = self.build_values(column=column)
        
        if set_values:
            self._key = column
            self._values = values
            values_valid = True

        if return_dict:
            values = {name : values}

        return values

    def mass_transpose(self, columns):
        """Queries for sets of distinct values, using self.transpose()

        :param columns: SQLAlchemy Column object(s)
        :type columns: Column
        :type columns: list
        :returns: Distinct values fetched from given and innate restraints.
        :rtype: dict
        """
        self.check()

        if not self._values:
            return {}

        if not isinstance(columns, list):
            columns = [columns]
        
        values={}
        for column in columns:
            column_values = self.transpose(column, return_dict=True)
            values.update(column_values)

        return values

    def retrieve(self, raw_columns):
        """Queries for distinct data for each value in the Filter object.

        :param columns: SQLAlchemy Column object(s)
        :type columns: Column
        :type columns: list
        :type columns: str
        :returns: Distinct values for each Filter value.
        :rtype: dict{dict}
        """
        self.check()

        if not self._values:
            return {}
        
        if not isinstance(raw_columns, list):
            raw_columns = [raw_columns]

        columns = []
        for raw_column in raw_columns:
            columns.append(self.convert_column_input(raw_column))

        values = {}
        for value in self._values:
            where_clause = (self._key == value)
            values.update({value : {}}) 

            #For each column for each value, add the respective data to a dict.
            for i in range(len(columns)):
                query = q.build_distinct(self._graph, columns[i], 
                                                        where=where_clause)
                value_data = q.first_column(self._engine, query)

                values[value].update({columns[i].name : value_data})

        return values

    def refresh(self):
        """Re-queries for the Filter's values.
        """
        self.check()

        if self._values_valid:
            return

        values = self.build_values()
        self._values = values
        self._values_valid = True

    def update(self):
        """Queries using the Filter's key and its stored BinaryExpressions.
        """
        self.check()

        if not self._values_valid:
            self.refresh()

        if self._updated:
            return

        where_clauses = self.build_where_clauses()
        values = self.build_values(where=where_clauses)
        self._values = values

        self._updated = True
        self._values_valid = True

    def sort(self, raw_column):
        """Re-queries for the Filter's values, applying a ORDER BY clause.
       
        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        """
        order_by_clause = self.convert_column_input(raw_column)

        values = self.build_values(order_by=order_by_clause)
        self._values = values
        self._values_valid = True

    def reset(self):
        """Resets all filters, values, and Filter state conditions.
        """
        self._filters = {}
        self._values = []
        self._values_valid = True
        self._updated = True

    def print_results(self): 
        """Prints the Filter object's values in a formatted way.
        """
        if not self._values:
            return

        if not self._values_valid:
            self.refresh()

        print(" " + "_"*57 + " ")
        print("| %-20s" % \
             (f"{self.key.name} Results:") \
             + " "*36 + "|")
        print("|" + "-"*57 + "|")
        for row in range(0, len(self._values), 3):
            result_row = self._values[row:row+3]
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
   
    def hits(self):
        """Gets the number of a Filter object's values.
        """
        if self.values == None:
            return 0

        return len(self._values)    

    def group(self, raw_column): 
        """Queries and separates Filter object's values based on a Column.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        """
        self.check()

        column = self.convert_column_input(raw_column)

        groups = self.transpose(column)
        
        group_results = {}
        for group in groups:
            where_clause = (column == group) 
            values = self.build_values(where=[where_clause])
            group_results.update({group : values})

        return group_results

    def copy(self):
        """Returns a copy of a Filter object.
        """
        copy = Filter()
        copy._connected = self._connected
        copy._updated = self.updated
        copy._values_valid = self.values_valid
        copy._filters = self.copy_filters()
        copy._engine = self.engine
        copy._graph = self.graph
        copy._key = self.key
        copy._values = self.values

        return copy
 
    def copy_filters(self):
        """Returns a copy of a Filter object's filter dictionary.
        """
        filters = {}
        #Copies of each list are required so that multiple filter copies
        #are not adding and removing from the same space in memory.
        for filter in self._filters.keys():
            clauses = self._filters[filter].copy()
            filters.update({filter : clauses})

        return filters

