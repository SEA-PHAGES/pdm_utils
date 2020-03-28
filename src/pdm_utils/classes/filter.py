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

        if isinstance(alchemist, AlchemyHandler):
            if not alchemist.connected:
                alchemist.connect(ask_database=True)

            if alchemist.graph == None:
                alchemist.build_graph()
                
            self._engine = alchemist.engine
            self._graph = alchemist.graph
       
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

    @property
    def graph(self):
        graph = self._graph
        return graph

    @graph.setter
    def graph(self, graph):
        if not isinstance(graph, Graph):
            raise TypeError

        self._graph = graph

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
        if not isinstance(key, Column):
            raise TypeError("Filter key value must be of type SqlAlchemy Column.")
        
        table = self._graph.nodes[key.table.name]["table"]
        column = table.columns[key.name]

        self._key = key

    def setup(loader): 
        load_filter(self, loader)
    
    def check(self):
        if not isinstance(self.engine, Engine): 
            raise AttributeError("Filter object is missing valid engine.")

        if not isinstance(self.graph, Graph):
            raise AttributeError("Filter object is missing valid graph.")

        if not isinstance(self._key, Column):
            raise AttributeError("Filter object is missing valid column key.")

    def add(self, filter):
        where_clause = q.build_where_clause(self.graph, filter)
        parsed_filter = parsing.parse_filter(filter)
        filter_left = parsed_filter[0] + "." + parsed_filter[1]\
                    + parsed_filter[2]
       
        if filter_left not in self._filters.keys():
            self._filters.update({filter_left : [where_clause]}) 
            
        else:  
            self._filters[filter_left].append(where_clause)

        self._updated = False

    def remove(self, filter):
        parsed_filter = parsing.parse_filter(filter)
        filter_left = parsed_filter[0] + "." + parsed_filter[1]\
                    + parsed_filter[2]

        if filter_left in self._filters.keys():
            filters = self._filters[filter_left]
            if len(filters) == 1:
                self._filters.pop(filter_left)

            else:
                for clause in filters:
                    if clause.right.value == parsed_filter[3]:
                        filters.remove(clause)
        
        self._updated = False

    def build_where_clauses(self):
        where_clauses = []
        for filter in self._filters.keys():
            where_clauses = where_clauses + self._filters[filter]

        return where_clauses

    def build_values(self, where=None, order_by=None):
        self.check()

        values = []

        if where == None:
            where_clauses = []
        else:
            where_clauses = where

        if self.values != []:
            where_clauses.append(self.key.in_(self.values))

        query = q.build_select(self.graph, [self._key], where=where_clauses)

        if isinstance(order_by, Column):
            query = query.order_by(order_by)

        proxy = self.engine.execute(query)
        results = proxy.fetchall()

        for result in results:
            values.append(result[0])

        return values
  
    def transpose(self, column, return_dict=False):
        self.check()

        if not self.values:
            return []
       
        name = ""

        if isinstance(column, str):
            column = q.get_column(self.graph.graph["metadata"], column)
            name = column.name

        elif isinstance(column, Column):
            name = column.name         
    
        else:
            raise TypeError
    
        where_clause = (self.key.in_(self.values))
        query = q.build_distinct(self.graph, [column], where=[where_clause])
        
        proxy = self.engine.execute(query)
        results = proxy.fetchall()

        values = []
        for result in results:
            values.append(result[0])

        if return_dict:
            values = {name : values}

        return values

    def retrieve(self, selectables):
        self.check()

        if not self.values:
            return {}

        where_clause = (self.key.in_(self.values))
    
        results = {}
        for selectable in selectables:
            values = self.transpose(selectable, return_dict=True) 
            results.update(values)            

        return results

    def refresh(self):
        self.check()

        if self._values_valid:
            return

        values = self.build_values()
        self._values = values
        self._values_valid = True

    def update(self):
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

    def sort(self, column):
        if isinstance(column, Column):
            order_by_clause = column

        elif isinstance(column, str):
            order_by_clause = self._graph.get_column(column)
        else:
            raise TypeError("Sort column must be either a string or a Column object")

        values = self.build_values(order_by=order_by_clause)
        self._values = values
        self._values_valid = True

    def reset(self):
        self._filters = {}
        self._values = []
        self._values_valid = True
        self._updated = True

    def print_results(self): 
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
        if self.values == None:
            return 0

        return len(self._values)    

    def group(self, column): 
        self.check()

        if isinstance(column, str):
            column = q.get_column(self.graph.graph["metadata"], column)

        elif isinstance(column, Column):
            pass
        else:
            raise TypeError

        groups = self.transpose(column)

        group_results = {}
        for group in groups:
            where_clause = (column == group) 
            values = self.build_values(where=[where_clause])
            group_results.update({group : values})

        return group_results

    def copy(self):
        copy = Filter()
        copy._updated = self.updated
        copy._values_valid = self.values_valid
        copy._filters = self.copy_filters()
        copy._engine = self.engine
        copy._graph = self.graph
        copy._key = self.key
        copy._values = self.values

        return copy
 
    def copy_filters(self):
        filters = {}
        for filter in self._filters.keys():
            clauses = self._filters[filter].copy()
            filters.update({filter : clauses})

        return filters

