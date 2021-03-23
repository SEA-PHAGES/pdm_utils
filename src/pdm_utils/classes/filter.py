"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""
from networkx import Graph
from sqlalchemy import Column
from sqlalchemy import and_
from sqlalchemy import or_
from sqlalchemy.engine.base import Engine
from sqlalchemy.orm.decl_api import DeclarativeMeta
from sqlalchemy.sql.elements import BooleanClauseList
from sqlalchemy.sql.elements import Label
from sqlalchemy.orm.session import Session

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import basic
from pdm_utils.functions import cartography
from pdm_utils.functions import parsing
from pdm_utils.functions import querying as q

# -----------------------------------------------------------------------------
# GLOBAL VARIABLES

COLUMN_TYPES = [Column, Label]


class Filter:
    def __init__(self, alchemist=None, key=None):
        self._engine = None
        self._graph = None
        self._session = None
        self._mapper = None

        self._connected = False

        if isinstance(alchemist, AlchemyHandler):
            if not alchemist.connected:
                alchemist.connect(ask_database=True)

            self._engine = alchemist.engine
            self._graph = alchemist.graph
            self._session = alchemist.session
            self._mapper = alchemist.mapper

            self._connected = True

        if isinstance(key, Column):
            self._key = key
        else:
            self._key = None

        self._values = []
        self._values_valid = True

        self._filters = []
        self._updated = True
        self._or_index = -1

        self.verbose = False

# -----------------------------------------------------------------------------
# FILTER CONNECTION HANDLING

    def connect(self, alchemist=None):
        """Connect Filter object to a database with an AlchemyHandler.

        :param alchemist: An AlchemyHandler object.
        :type alchemist: AlchemyHandler
        """
        if alchemist is not None:
            self.link(alchemist)
            return

        if self._connected:
            return

        alchemist = AlchemyHandler()
        alchemist.connect(ask_database=True)

        self._engine = alchemist.engine
        self._graph = alchemist.graph
        self._session = alchemist.session
        self._mapper = alchemist.mapper

        self._connected = True

    def link(self, alchemist):
        """Connect Filter object to a database with an existing AlchemyHandler.

        :param alchemist: An AlchemyHandler object.
        :type alchemist: AlchemyHandler
        """
        if not isinstance(alchemist, AlchemyHandler):
            raise TypeError("AlchemyHandler to be linked is not "
                            "of type AlchemyHandler.")

        if not alchemist.connected:
            alchemist.connect(ask_database=True)

        self._engine = alchemist.engine
        self._graph = alchemist.graph
        self._session = alchemist.session
        self._mapper = alchemist.mapper

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
                raise AttributeError("Filter object is missing valid Engine.")
            if not isinstance(self._graph, Graph):
                raise AttributeError("Filter object is missing valid Graph.")
            if not isinstance(self._session, Session):
                raise AttributeError("Filter object is missing valid Session.")

        if not isinstance(self._key, Column):
            raise AttributeError("Filter object is missing valid column key.")

# -----------------------------------------------------------------------------
# FILTER CONNECTION HANDLING PROPERTIES

    @property
    def connected(self):
        connected = self._connected
        return connected

    @property
    def engine(self):
        engine = self._engine
        return engine

    @property
    def graph(self):
        graph = self._graph
        return graph

    @property
    def session(self):
        session = self._session
        return session

    @property
    def mapper(self):
        mapper = self._mapper
        return mapper

# -----------------------------------------------------------------------------
# FILTER VALUE HANDLING

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

    def sort(self, raw_columns):
        """Re-queries for the Filter's values, applying a ORDER BY clause.

        :param raw_column: SQLAlchemy Column object(s) or object name(s).
        :type raw_columns: Column
        :type raw_columns: str
        :type raw_columns: list[Column]
        :type raw_columns: list[str]
        """
        self.check()

        columns = self.get_columns(raw_columns)

        query = q.build_select(self._graph, self._key, order_by=columns)

        values = q.first_column(self._engine, query, in_column=self._key,
                                values=self._values)
        self._values = values
        self._values_valid = True

# -----------------------------------------------------------------------------
# FILTER VALUE HANDLING PROPERTIES

    @property
    def values(self):
        if self._values is not None:
            values = self._values.copy()
            return values
        else:
            return []

    @values.setter
    def values(self, values):
        if not (isinstance(values, list) or values is None):
            raise TypeError("Filter object values must be of type list.")
        self._values = values
        self._values_valid = False
        self._updated = False

    @property
    def key(self):
        key = self._key
        return key

    @key.setter
    def key(self, key):
        if isinstance(key, Column):
            self._key = key
        elif isinstance(key, str):
            if self.graph is None:
                raise ValueError("String key input requires MySQL connection.")

            metadata = self.graph.graph["metadata"]

            try:
                self._key = q.get_column(self.graph.graph["metadata"], key)
            except:
                try:
                    table_obj = q.get_table(metadata, key)
                except:
                    raise ValueError("Inputted string key is neither a valid "
                                     "MySQL column or table.")

                self._key = list(table_obj.primary_key.columns)[0]

        else:
            raise TypeError("Filter key value is invalid."
                            "Filter key must be one of the following: \n"
                            "SQLAlchemy Column\n"
                            "MySQL column string\n"
                            "MySQL table string\n")

    @property
    def values_valid(self):
        values_valid = self._values_valid
        return values_valid

# FILTER CLAUSE HANDLING
# -----------------------------------------------------------------------------

    def reset(self):
        """Resets all filters, values, and Filter state conditions.
        """
        self.reset_filters()

        self._values = []
        self._values_valid = True

    def reset_filters(self):
        """Resets all filters and relevant Filter state condition.
        """
        self._filters = []
        self._updated = True
        self._or_index = -1

    def parenthesize(self):
        """Condense current filters into an isolated clause"""
        if self._or_index > 0:
            conditionals = self.build_where_clauses()

            self.reset_filters()
            self.new_or_()
            or_block = self._filters[self._or_index]

            or_block.update({"parenthetical": conditionals})

            self._updated = False

    def and_(self, filter):
        """Add an and conditional to the Filter object class.

        :param_filter: Formatted MySQL WHERE clause.
        :type_filter: str
        """
        if self._or_index < 0:
            self.new_or_()

        where_clause = q.build_where_clause(self.graph, filter)
        filter_key = parsing.create_filter_key(filter)

        or_block = self._filters[self._or_index]

        # Stores the BinaryExpression with a key of the Column/Operator pairing
        or_block.update({filter_key: where_clause})

        self._updated = False

    def new_or_(self):
        """Create a new conditional block to the Filter object class.
        """
        if self._or_index >= 0:
            if not self._filters[self._or_index]:
                return

        self._filters.append({})
        self._or_index += 1

    def remove(self, filter):
        """Remove an and filter from the current block of and conditionals.

        :param filter: Formatted MySQL WHERE clause.
        :type filter: str
        """
        filter_key = parsing.create_filter_key(filter)

        or_block = self._filters[self._or_index]

        # Indexes into the dictionary using the Column/Operator pairing.
        if filter_key in or_block.keys():
            or_block.pop(filter_key)
        else:
            raise IndexError(f"Filter {filter} is not within the current "
                             "block of and conditionals.")

        self._updated = False

    def add(self, filter_string):
        """Add a MySQL where filter(s) to the Filter object class.

        :param filter: Formatted MySQL WHERE clause.
        :type filter: str
        """
        filters = parsing.parse_cmd_string(filter_string)
        while(filters):
            and_filters = filters[0]
            for filter in and_filters:
                self.and_(filter)

            filters.pop(0)

            if filters:
                self.new_or_()

    def build_where_clauses(self):
        """Builds BinaryExpression objects from stored Filter object filters.

        :returns: A list of SQLAlchemy WHERE conditionals.
        :rtype: list
        """
        and_conditionals = []

        for or_block in self._filters:
            where_clauses = []
            for filter_key in or_block.keys():
                where_clauses.append(or_block[filter_key])

            if where_clauses:
                and_conditionals.append(and_(*where_clauses))

        if len(and_conditionals) > 1:
            conditionals = or_(*and_conditionals)
        elif len(and_conditionals) == 1:
            conditionals = and_conditionals
        else:
            conditionals = []

        return conditionals

# -----------------------------------------------------------------------------
# FILTER CLAUSE HANDLING PROPERTIES

    @property
    def filters(self):
        filters = self.copy_filters()
        return filters

    @property
    def updated(self):
        updated = self._updated
        return updated

    @property
    def or_index(self):
        or_index = self._or_index
        return or_index

    @or_index.setter
    def or_index(self, index):
        or_count = len(self._filters)
        if index < 0 or index > (or_count - 1):
            raise IndexError("Or_index is out of bounds for "
                             "the current or_ blocks.\n"
                             f"Current or_ blocks stored are {or_count}.")

        self._or_index = index

# -----------------------------------------------------------------------------
# FILTER QUERYING
    def build_values(self, where=None, column=None, raw_bytes=False,
                     limit=8000):
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
            column_obj = self._key
        else:
            column_obj = self.get_column(column)

        if where is not None:
            if isinstance(where, list) or isinstance(where, BooleanClauseList):
                base_clauses = where
            else:
                base_clauses = [where]
        else:
            base_clauses = []

        query = q.build_distinct(self.graph, column_obj, where=base_clauses,
                                 add_in=self._key)

        values = q.first_column(self.engine, query, in_column=self._key,
                                values=self._values, limit=limit)

        if not raw_bytes:
            if column_obj.type.python_type is bytes:
                values = basic.convert_to_decoded(values)

        return values

    def select(self, raw_columns, return_dict=True):
        """Queries for data conditioned on the values in the Filter object.

        :param columns: SQLAlchemy Column object(s)
        :type columns: Column
        :type columns: str
        :type columns: list[Column]
        :type columns: list[str]
        :param return_dict: Toggle whether to return data as a dictionary.
        :type return_dict: Boolean
        :returns: SELECT data conditioned on the values in the Filter object.
        :rtype: dict
        :rtype: list[RowProxy]
        """
        self.check()

        columns = self.get_columns(raw_columns)

        query = q.build_select(self._graph, columns, add_in=self._key)
        results = q.execute(self._engine, query, in_column=self._key,
                            values=self._values, return_dict=return_dict)

        return results

    def query(self, table_map):
        """Queries for ORM object instances conditioned on Filter values.

        :param table_map: SQLAlchemy ORM map object.
        :type table_map: str
        :type table_map: DeclarativeMeta
        :returns: List of mapped object instances.
        :rtype: list
        """
        self.check()

        if not self._values:
            return []

        if isinstance(table_map, str):
            table_map = cartography.get_map(self._mapper, table_map)
        elif isinstance(table_map, DeclarativeMeta):
            pass
        else:
            raise TypeError("Table map object must be either "
                            "type str or DeclarativeMeta.")

        in_clause = self._key.in_(self._values)

        instances = q.query(self._session, self._graph, table_map,
                            where=in_clause)
        return instances

    def transpose(self, raw_column, return_dict=False, set_values=False,
                  raw_bytes=False, filter=False):
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
            if return_dict:
                return {}
            else:
                return []

        column = self.get_column(raw_column)
        name = column.name

        where_clauses = None
        if filter:
            where_clauses = self.build_where_clauses()

        values = self.build_values(column=column, where=where_clauses,
                                   raw_bytes=raw_bytes)

        if set_values:
            self._key = column
            self._values = values

        if return_dict:
            values = {name: values}

        return values

    def mass_transpose(self, raw_columns, raw_bytes=False, filter=False):
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

        if not isinstance(raw_columns, list):
            raw_columns = [raw_columns]

        values = {}
        for column in raw_columns:
            column_values = self.transpose(column, return_dict=True,
                                           raw_bytes=raw_bytes,
                                           filter=filter)
            values.update(column_values)

        return values

    def group(self, raw_column, raw_bytes=False, filter=False):
        """Queries and separates Filter object's values based on a Column.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        """
        self.check()

        column = self.get_column(raw_column)

        groups = self.transpose(column)

        where_clauses = []
        if filter:
            where_clauses = self.build_where_clauses()

        group_results = {}
        for group in groups:
            if column.type.python_type == bytes:
                group_clauses = where_clauses + \
                                            [(column == group.encode("utf_8"))]
            else:
                group_clauses = where_clauses + [(column == group)]

            values = self.build_values(where=group_clauses,
                                       raw_bytes=raw_bytes)
            group_results.update({group: values})

        return group_results

    def retrieve(self, raw_columns, raw_bytes=False, filter=False):
        """Queries for distinct data for each value in the Filter object.

        :param columns: SQLAlchemy Column object(s)
        :type columns: Column
        :type columns: str
        :type columns: list[Column]
        :type columns: list[str]
        :returns: Distinct values for each Filter value.
        :rtype: dict{dict}
        """
        self.check()

        if not self._values:
            return {}

        columns = self.get_columns(raw_columns)

        where_clauses = []
        if filter:
            where_clauses = self.build_where_clauses()

        values = {}
        for value in self._values:
            compare_value = value
            if self._key.type.python_type == bytes and value is not None:
                compare_value = value.encode("utf-8")

            value_clauses = where_clauses + [(self._key == compare_value)]
            values.update({value: {}})

            # For each column for each value, add the respective data to a dict
            for i in range(len(columns)):
                query = q.build_distinct(self._graph, columns[i],
                                         where=value_clauses)
                value_data = q.first_column(self._engine, query)

                if not raw_bytes:
                    if columns[i].type.python_type == bytes:
                        value_data = basic.convert_to_decoded(value_data)

                values[value].update({columns[i].name: value_data})

        return values

    def get_column(self, raw_column):
        """Converts a column input, string or Column, to a Column.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: Column
        :type raw_column: str
        """
        self.check()

        if isinstance(raw_column, str):
            column = q.get_column(self.graph.graph["metadata"], raw_column)
        elif type(raw_column) in COLUMN_TYPES:
            column = raw_column
        else:
            raise TypeError(
                        "Column must be either a string or a Column object")

        return column

    def get_columns(self, raw_columns):
        """Converts a column input list, string or Column, to a list of Columns.

        :param raw_column: SQLAlchemy Column object or object name.
        :type raw_column: list[Column]
        :type raw_column: list[str]
        :returns: Returns SQLAlchemy Columns
        :rtype: list[Column]
        """
        if not isinstance(raw_columns, list):
            raw_columns = [raw_columns]

        columns = []
        for raw_column in raw_columns:
            columns.append(self.get_column(raw_column))

        return columns

# -----------------------------------------------------------------------------
# FILTER QUALITY-OF-LIFE

    def hits(self):
        """Gets the number of a Filter object's values.
        """
        if self.values is None:
            return 0

        return len(self._values)

    def copy(self):
        """Returns a copy of a Filter object.
        """
        copy = Filter()
        copy._connected = self._connected
        copy._updated = self.updated
        copy._values_valid = self.values_valid
        copy._filters = self.copy_filters()
        copy._engine = self._engine
        copy._graph = self._graph
        copy._session = self._session
        copy._key = self.key
        copy._values = self.values

        return copy

    def copy_filters(self):
        """Returns a copy of a Filter object's filter dictionary.
        """
        filters_copy = []
        # Copies of each list are required so that multiple filter copies
        # are not adding and removing from the same space in memory.
        for or_block in self._filters:
            new_block = {}

            for clause in or_block.keys():
                new_block.update({clause: or_block[clause]})

            filters_copy.append(new_block)

        return filters_copy

    def print_results(self):
        """Prints the Filter object's values in a formatted way.
        """
        if not self._values:
            return

        if not self._values_valid:
            self.refresh()

        print(" " + "_"*57 + " ")
        print("| %-20s" %
              (f"{self.key.name} Results:")
              + " "*36 + "|")
        print("|" + "-"*57 + "|")
        for row in range(0, len(self._values), 3):
            result_row = self._values[row:row+3]
            if len(result_row) == 3:
                print("| %-20s %-20s %-13s |" %
                      (result_row[0], result_row[1], result_row[2]))
            elif len(result_row) == 2:
                print("| %-20s %-20s %-13s" %
                      (result_row[0], result_row[1], "")
                      + " " + "|")

            elif len(result_row) == 1:
                print("| %-20s %-20s %-13s" %
                      (result_row[0], "", "")
                      + " " + "|")
        print("|" + "_"*57 + "|")
