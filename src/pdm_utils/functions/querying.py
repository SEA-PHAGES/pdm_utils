from collections import OrderedDict

from networkx import Graph
from networkx import shortest_path
from sqlalchemy import and_
from sqlalchemy import Column
from sqlalchemy import join
from sqlalchemy import MetaData
from sqlalchemy import select
from sqlalchemy import Table
from sqlalchemy.orm.decl_api import DeclarativeMeta
from sqlalchemy.sql import func
from sqlalchemy.sql import functions
from sqlalchemy.sql.elements import BinaryExpression
from sqlalchemy.sql.elements import BooleanClauseList
from sqlalchemy.sql.elements import Grouping
from sqlalchemy.sql.elements import Label
from sqlalchemy.sql.elements import UnaryExpression

from pdm_utils.functions import basic
from pdm_utils.functions import parsing

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

COLUMN_TYPES = [Column, Table, functions.count, BinaryExpression,
                UnaryExpression, Label, DeclarativeMeta, Grouping]


# SQLALCHEMY OBJECT RETRIEVAL
# Functions that functionalize retrieval of SqlAlchemy objects.
# -----------------------------------------------------------------------------
def get_table(metadata, table):
    """Get a SQLAlchemy Table object, with a case-insensitive input.

    :param metadata: Reflected SQLAlchemy MetaData object.
    :type metadata: MetaData
    :param table: Case-insensitive table name.
    :type table: str
    :returns:
        Corresponding SQLAlchemy Table object.
    :rtype: Table
    """
    # Ensures table name is accurate before indexing metadata.
    table = parsing.translate_table(metadata, table)
    table_obj = metadata.tables[table]
    return table_obj


def get_column(metadata, column):
    """Get a SQLAlchemy Column object, with a case-insensitive input.
    Input must be formatted {Table_name}.{Column_name}.

    :param metadata: Reflected SQLAlchemy MetaData object.
    :type metadata: MetaData
    :param table: Case-insensitive column name.
    :type table: str
    :returns: Corresponding SQLAlchemy Column object.
    :rtype: Column
    """
    parsed_column = parsing.parse_column(column)

    table_obj = get_table(metadata, parsed_column[0])  # Retrieves SQL Table.
    # Ensures column name is accurate before indexing Table.
    column = parsing.translate_column(metadata, column)

    columns_dict = dict(table_obj.columns)  # Converts to indexable object.
    column_obj = columns_dict[column]

    return column_obj


# MySQL DATABASE MAPPING AND TRAVERSAL
# Functions to traverse MySQL database structures and map valid JOIN clauses.
# -----------------------------------------------------------------------------
def build_graph(metadata):
    """Get a NetworkX Graph object populated from a SQLAlchemy MetaData object.

    :param metadata: Reflected SQLAlchemy MetaData object.
    :type metadata: MetaData
    :returns: Populated and structured NetworkX Graph object.
    :rtype: Column
    """
    if not isinstance(metadata, MetaData):
        raise TypeError("Graph requires MetaData object, "
                        f"object passed was of type {type(metadata)}.")

    graph = Graph(metadata=metadata)  # Stores metadata in Graph object.
    for table in metadata.tables.keys():  # Creates nodes from all Tables.
        table_object = metadata.tables[table]
        graph.add_node(table_object.name, table=table_object)

    # Stores Table objects in Nodes
    for target_table in metadata.tables.keys():
        target_table_obj = metadata.tables[target_table]

        # Uses the Table object foreign keys to create Graph edges.
        for foreign_key in target_table_obj.foreign_keys:
            referent_column = foreign_key.column  # Retreives source Column
            referent_table_obj = referent_column.table
            referent_table = referent_table_obj.name

            # Stores the foreign key as a table name key pair [Target][Source]
            graph.add_edge(target_table, referent_table, key=foreign_key)

    return graph


def build_fromclause(db_graph, columns):
    """Get a joined table from pathing instructions for joining MySQL Tables.
    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param columns: SQLAlchemy Column object(s).
    :type columns: Column
    :type columns: list
    :returns: SQLAlchemy Table object containing left outer-joined tables.
    :rtype: Table
    """
    table_list = get_table_list(columns)
    table_pathing = get_table_pathing(db_graph, table_list)
    joined_table = join_pathed_tables(db_graph, table_pathing)

    return joined_table


def build_onclause(db_graph, source_table, adjacent_table):
    """Creates a SQLAlchemy BinaryExpression object for a MySQL ON clause
       expression

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param source_table: Case-insensitive MySQL table name.
    :type source_table: str
    :param adjacent_table: Case-insensitive MySQL table name.
    :type source_table: str
    :returns: MySQL foreign key related SQLAlchemy BinaryExpression object.
    :rtype: BinaryExpression
    """
    # Retrieves accurate table names
    source_table = parsing.translate_table(db_graph.graph["metadata"],
                                           source_table)
    adjacent_table = parsing.translate_table(db_graph.graph["metadata"],
                                             adjacent_table)

    if source_table == adjacent_table:
        source_table = db_graph.nodes[source_table]["table"]
        return source_table

    edge = db_graph[source_table][adjacent_table]  # Indexes for graph edge.
    foreign_key = edge["key"]  # Retreives the stored foreign key object

    referent_column = foreign_key.column  # Gets the original primary key
    referenced_column = foreign_key.parent  # Gets the foreign key

    onclause = referent_column == referenced_column

    return onclause


def get_table_list(columns):
    """Get a nonrepeating list SQLAlchemy Table objects from Column objects.

    :param columns: SQLAlchemy Column object(s).
    :type columns: Column
    :type columns: list
    :returns: List of corresponding SQLAlchemy Table objects.
    :rtype: list
    """
    table_list = []

    if not isinstance(columns, list):
        columns = [columns]

    for column in columns:
        table_list.append(extract_column(column).table)

    # OrderedDict is used to remove repeats, trading space for runtime
    table_list = list(OrderedDict.fromkeys(table_list))

    return table_list


def get_table_pathing(db_graph, table_list, center_table=None):
    """Get pathing instructions for joining MySQL Table objects.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param table_list: List of SQLAlchemy Table objects.
    :type table_list: list[Table]
    :param center_table: SQLAlchemy Table object to begin traversals from.
    :type center_table: Table
    :returns: 2-D list containing the center table and pathing instructions.
    :rtype: list
    """
    table_list = table_list.copy()

    # Removes the center_table from the list, if specified
    if center_table is not None:
        table_list.remove(center_table)
    # Creates a random center_table
    else:
        center_table = table_list[0]
        table_list.remove(center_table)

    table_paths = []
    for table in table_list:
        # NetworkX traversal algorithm finds path
        path = shortest_path(db_graph, center_table.name, table.name)

        if not path:
            raise ValueError("Operation cannot be performed. "
                             f"Table {table.name} is not connected by any "
                             f"foreign keys to table {center_table}.")

        table_paths.append(path)

    table_pathing = [center_table, table_paths]
    return table_pathing


def join_pathed_tables(db_graph, table_pathing):
    """Get a joined table from pathing instructions for joining MySQL Tables.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param table_pathing: 2-D list containing a Table and pathing lists.
    :type table_pathing: list
    :returns: SQLAlchemy Table object containing left outer-joined tables.
    :rtype: Table
    """
    center_table = table_pathing[0]  # The center_table from table pathing
    joined_tables = center_table  # Creates a Table 'base' to add to

    joined_table_set = set()
    joined_table_set.add(center_table.name)  # Creates running table history

    for path in table_pathing[1]:
        for index in range(len(path)):
            table = path[index]
            previous_table = path[index-1]

            # If the table has not already been added:
            if table not in joined_table_set:
                joined_table_set.add(table)
                table_object = db_graph.nodes[table]["table"]
                # Build a MySQL ON clause to left outer-join the table
                onclause = build_onclause(db_graph, previous_table, table)
                joined_tables = join(joined_tables, table_object,
                                     isouter=True, onclause=onclause)

    return joined_tables


# SQLALCHEMY COLUMN HANDLING
# Functions to deal with column formats and constructing column expressions.
# -----------------------------------------------------------------------------
def build_where_clause(db_graph, filter_expression):
    """Creates a SQLAlchemy BinaryExpression object from a MySQL WHERE
       clause expression.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param filter_expression: MySQL where clause expression.
    :type filter_expression: str
    :returns: MySQL expression-related SQLAlchemy BinaryExpression object.
    :rtype: BinaryExpression
    """
    filter_params = parsing.parse_filter(filter_expression)

    # Retrieve Table and Column objects from split filter expression strings
    table_object = db_graph.nodes[filter_params[0]]["table"]

    column_name = parsing.translate_column(
                                db_graph.graph["metadata"],
                                f"{table_object.name}.{filter_params[1]}")
    column_object = table_object.columns[column_name]

    # Checks the operator and Column type compatability
    parsing.check_operator(filter_params[2], column_object)

    right = filter_params[3]  # Stores the expressions 'right' value

    if column_object.type.python_type == bytes:
        if isinstance(right, list):
            for i in range(len(right)):
                right[i] = right[i].encode("utf-8")
        else:
            right = right.encode("utf-8")

    if right == "None":
        right = None

    where_clause = None

    if filter_params[2] == "=":
        where_clause = (column_object == right)
    elif filter_params[2] == "LIKE":
        where_clause = (column_object.like(right))
    elif filter_params[2] == "!=" or filter_params[2] == "IS NOT":
        where_clause = (column_object != right)
    elif filter_params[2] == ">":
        where_clause = (column_object > right)
    elif filter_params[2] == ">=":
        where_clause = (column_object >= right)
    elif filter_params[2] == "<":
        where_clause = (column_object < right)
    elif filter_params[2] == "<=":
        where_clause = (column_object <= right)
    elif filter_params[2] == "IN":
        where_clause = (column_object.in_(right))
    elif filter_params[2] == "NOT IN":
        where_clause = (column_object.notin_(right))

    return where_clause


def extract_column(column, check=None):
    """Get a column from a supported SQLAlchemy Column-related object.

    :param column: SQLAlchemy Column-related object.
    :type column: BinaryExpression
    :type column: Column
    :type column: count
    :type column: UnaryExpression
    :param check: SQLAlchemy Column-related object type.
    :type check: <type Column>
    :type check: <type count>
    :type check: <type UnaryExpression>
    :type check: <type BinaryExpression>
    :returns: Corresponding SQLAlchemy Column object.
    :rtype: Column
    """
    if check is not None:
        if not isinstance(column, check):
            raise TypeError(f"Type {type(column)} object passed as column "
                            "failed a specified type check.")

    if isinstance(column, Column):
        pass
    # For handling SQLAlchemy Table input
    elif isinstance(column, Table):
        column = list(column.primary_key.columns)[0]
    # For handling SQLAlchemy count(Column) expressions
    elif isinstance(column, functions.count):
        column = column.clauses.clauses[0]
    # For handling SQLAlchemy Column.distinct expressions
    elif isinstance(column, UnaryExpression):
        column = column.element
    elif isinstance(column, Label):
        column = column.element
    elif isinstance(column, Grouping):
        column = column.element.element
    # For handling SQLAlchemy comparison expressions
    elif isinstance(column, BinaryExpression):
        expression = column.left
        # For handling SQLAlchemy Column.in_() expressions
        if isinstance(expression, UnaryExpression):
            column = expression.element
        # For handling SQLAlchemy count(Column) comparisons
        elif isinstance(expression, functions.count):
            column = expression.clauses.clauses[0]
        # For handling SQLAlchemy Column.distinct() comparisons
        elif isinstance(expression, Grouping):
            column = expression.element.element
        elif isinstance(expression, Column):
            column = expression
        else:
            raise TypeError(f"BinaryExpression type {type(expression)} "
                            f"of expression {expression} is not supported.")
    elif isinstance(column, DeclarativeMeta):
        try:
            table = column.__table__
        except:
            raise ValueError("SQLAlchemy Map-related object is not supported.")
        finally:
            column = list(table.primary_key.columns)[0]
    else:
        raise TypeError(f"Input type {type(column)} is not a derivative "
                        "of a SqlAlchemy Column.")

    return column


def extract_columns(columns, check=None):
    """Get a column from a supported SQLAlchemy Column-related object(s).

    :param column: SQLAlchemy Column-related object.
    :type column: BinaryExpression
    :type column: Column
    :type column: count
    :type column: list
    :type column: UnaryExpression
    :param check: SQLAlchemy Column-related object type.
    :type check: <type Column>
    :type check: <type count>
    :type check: <type UnaryExpression>
    :type check: <type BinaryExpression>
    :returns: List of SQLAlchemy Column objects.
    :rtype: list[Column]
    """
    extracted_columns = []

    if isinstance(columns, list) or isinstance(columns, BooleanClauseList):
        for column in columns:
            if isinstance(column, BooleanClauseList):
                extracted_columns = extracted_columns + extract_columns(column)
            elif isinstance(column, Grouping):
                extracted_columns = extracted_columns + \
                                    extract_columns(column.element)
            else:
                extracted_columns.append(extract_column(column, check=check))
    elif columns is None:
        pass
    else:
        extracted_columns.append(extract_column(columns, check=check))

    return extracted_columns


# SQLALCHEMY SELECT HANDLING
# Functions that build and modify SqlAlchemy executable objects.
# -----------------------------------------------------------------------------
def build_select(db_graph, columns, where=None, order_by=None, add_in=None,
                 having=None, group_by=None):
    """Get MySQL SELECT expression SQLAlchemy executable.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param columns: SQLAlchemy Column object(s).
    :type columns: Column
    :type columns: list
    :param where: MySQL WHERE clause-related SQLAlchemy object(s).
    :type where: BinaryExpression
    :type where: list
    :param order_by: MySQL ORDER BY clause-related SQLAlchemy object(s).
    :type order_by: Column
    :type order_by: list
    :param add_in: MySQL Column-related inputs to be considered for joining.
    :type add_in: Column
    :type add_in: list
    :param having: MySQL HAVING clause-related SQLAlchemy object(s).
    :type having: BinaryExpression
    :type having: list
    :param group_by: MySQL GROUP BY clause-related SQLAlchemy object(s).
    :type group_by: Column
    :type group_by: list
    :returns: MySQL SELECT expression-related SQLAlchemy executable.
    :rtype: Select
    """
    where_columns = extract_columns(where)
    having_columns = extract_columns(having)
    add_in_columns = extract_columns(add_in, check=Column)
    order_by_columns = extract_columns(order_by, check=Column)
    group_by_columns = extract_columns(group_by, check=Column)

    if not isinstance(columns, list):
        columns = [columns]

    # Gathers all the tables from the SELECT, WHERE, and ORDER BY clauses
    total_columns = (columns + where_columns + order_by_columns
                     + add_in_columns + having_columns + group_by_columns)
    fromclause = build_fromclause(db_graph, total_columns)

    select_query = select(columns).select_from(fromclause)
    select_query = append_group_by_clauses(select_query, group_by)
    select_query = append_order_by_clauses(select_query, order_by)
    select_query = append_where_clauses(select_query, where)
    select_query = append_having_clauses(select_query, having)

    return select_query


def build_count(db_graph, columns, where=None, add_in=None):
    """Get MySQL COUNT() expression SQLAlchemy executable.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param columns: SQLAlchemy Column object(s).
    :type columns: Column
    :type columns: list
    :param where: MySQL WHERE clause-related SQLAlchemy object(s).
    :type where: BinaryExpression
    :type where: list
    :param add_in: MySQL Column-related inputs to be considered for joining.
    :type add_in: Column
    :type add_in: list
    :returns: MySQL COUNT() expression-related SQLAlchemy executable.
    :rtype: Select
    """
    where_columns = extract_columns(where)
    add_in_columns = extract_columns(add_in, check=Column)

    if not isinstance(columns, list):
        columns = [columns]

    # Gathers all the tables from the SELECT, WHERE, and ORDER BY clauses
    total_columns = columns + where_columns + add_in_columns
    fromclause = build_fromclause(db_graph, total_columns)
    # Converts SQLAlchemy Columns into SQLAlchemy COUNT-related objects
    column_params = []
    for column_param in columns:
        column_params.append(func.count(column_param))

    count_query = select(column_params).select_from(fromclause)
    count_query = append_where_clauses(count_query, where)

    return count_query


def build_distinct(db_graph, columns, where=None, order_by=None, add_in=None):
    """Get MySQL DISTINCT expression SQLAlchemy executable.

    :param db_graph: SQLAlchemy structured NetworkX Graph object.
    :type db_graph: Graph
    :param columns: SQLAlchemy Column object(s).
    :type columns: Column
    :type columns: list
    :param where: MySQL WHERE clause-related SQLAlchemy object(s).
    :type where: BinaryExpression
    :type where: list
    :param order_by: MySQL ORDER BY clause-related SQLAlchemy object(s).
    :type order_by: Column
    :type order_by: list
    :param add_in: MySQL Column-related inputs to be considered for joining.
    :type add_in: Column
    :type add_in: list
    :returns: MySQL DISTINCT expression-related SQLAlchemy executable.
    :rtype: Select
    """
    query = build_select(db_graph, columns, where=where,
                         order_by=order_by, add_in=add_in)

    distinct_query = query.distinct()  # Converts a SELECT to a DISTINCT
    return distinct_query


def append_where_clauses(executable, where_clauses):
    """Add WHERE SQLAlchemy BinaryExpression objects to a Select object.

    :param executable: SQLAlchemy executable query object.
    :type executable: Select
    :param where_clauses: MySQL WHERE clause-related SQLAlchemy object(s).
    :type where_clauses: BinaryExpression
    :type where_clauses: list
    :returns: MySQL expression-related SQLAlchemy exectuable.
    :rtype: Select
    """
    if where_clauses is None:
        return executable

    if isinstance(where_clauses, list):
        executable = executable.where(and_(*where_clauses))
    else:
        executable = executable.where(where_clauses)

    return executable


def append_having_clauses(executable, having_clauses):
    """Add HAVING SQLAlchemy Column objects to a Select object.

    :param executable: SQLAlchemy executable query object.
    :type executable: Select
    :param having_clauses: MySQL HAVING clause-related SQLAlchemy object(s).
    :type order_by_clauses: Column
    :type order_by_clauses: List
    :returns MySQL expression-related SQLAlchemy executable.
    :rtype: Select
    """
    if having_clauses is None:
        return executable

    if isinstance(having_clauses, list):
        executable = executable.having(and_(*having_clauses))
    else:
        executable = executable.having(having_clauses)

    return executable


def append_order_by_clauses(executable, order_by_clauses):
    """Add ORDER BY SQLAlchemy Column objects to a Select object.

    :param executable: SQLAlchemy executable query object.
    :type executable: Select
    :param order_by_clauses: MySQL ORDER BY clause-related SQLAlchemy object(s)
    :type order_by_clauses: Column
    :type order_by_clauses: list
    :returns: MySQL expression-related SQLAlchemy exectuable.
    :rtype: Select
    """
    if order_by_clauses is None:
        return executable

    if isinstance(order_by_clauses, list):
        for clause in order_by_clauses:
            executable = executable.order_by(clause)
    else:
        executable = executable.order_by(order_by_clauses)

    return executable


def append_group_by_clauses(executable, group_by_clauses):
    """Add GROUP BY SQLAlchemy Column objects to a Select object.

    :param executable: SQLAlchemy executable query object.
    :type executable: Select
    :param order_by_clauses: MySQL GROUP BY clause-related SQLAlchemy object(s)
    :type order_by_clauses: Column
    :type order_by_clauses: list
    :returns: MySQL expression-related SQLAlchemy exectuable.
    :rtype: Select
    """
    if group_by_clauses is None:
        return executable

    if isinstance(group_by_clauses, list):
        for clause in group_by_clauses:
            executable = executable.group_by(clause)
    else:
        executable = executable.group_by(group_by_clauses)

    return executable


# SQLALCHEMY EXECUTE QUERY FUNCTIONS
# Functions that execute SqlAlchemy select statements and handle outputs.
# -----------------------------------------------------------------------------
def execute(engine, executable, in_column=None, values=[], limit=8000,
            return_dict=True):
    """Use SQLAlchemy Engine to execute a MySQL query.

    :param engine: SQLAlchemy Engine object used for executing queries.
    :type engine: Engine
    :param executable: Input a executable MySQL query.
    :type executable: Select
    :type executable: str
    :param return_dict: Toggle whether execute returns dict or tuple.
    :type return_dict: Boolean
    :returns: Results from execution of given MySQL query.
    :rtype: list[dict]
    :rtype: list[tuple]
    """
    if values:
        if in_column is None:
            raise ValueError("Column input is required to condition "
                             "SQLAlchemy select for a set of values.")

        results = execute_value_subqueries(engine, executable,
                                           in_column, values,
                                           return_dict=return_dict,
                                           limit=limit)

    else:
        proxy = engine.execute(executable)

        results = proxy.fetchall()

        if return_dict:
            results_dicts = []
            for result in results:
                results_dicts.append(dict(result))

            results = results_dicts

    return results


def first_column(engine, executable, in_column=None, values=[], limit=8000):
    """Use SQLAlchemy Engine to execute and return the first column of fields.

    :param engine: SQLAlchemy Engine object used for executing queries.
    :type engine: Engine
    :param executable: Input an executable MySQL query.
    :type executable: Select
    :type executable: str
    :returns: A column for a set of MySQL values.
    :rtype: list[str]
    """
    if values:
        if in_column is None:
            raise ValueError("Column input is required to condition "
                             "SQLAlchemy select for a set of values.")

        values = first_column_value_subqueries(engine, executable,
                                               in_column, values,
                                               limit=limit)
    else:
        proxy = engine.execute(executable)
        results = proxy.fetchall()

        values = []
        for result in results:
            values.append(result[0])

    return values


def execute_value_subqueries(engine, executable, in_column, source_values,
                             return_dict=True, limit=8000):
    """Query with a conditional on a set of values using subqueries.

    :param engine: SQLAlchemy Engine object used for executing queries.
    :type engine: Engine
    :param executable: Input a executable MySQL query.
    :type executable: Select
    :type executable: str
    :param in_column: SQLAlchemy Column object.
    :type in_column: Column
    :param source_values: Values from specified MySQL column.
    :type source_values: list[str]
    :param return_dict: Toggle whether to return data as a dictionary.
    :type return_dict: Boolean
    :param limit: SQLAlchemy IN clause query length limiter.
    :type limit: int
    :returns: List of grouped data for each value constraint.
    :rtype: list
    """
    if not isinstance(in_column, Column):
        raise ValueError("Inputted column to conditional values against "
                         "is not a SqlAlchemy Column."
                         f"Object is instead type {type(in_column)}.")

    if not executable.is_derived_from(in_column.table):
        raise ValueError("Inputted column to conditional values against "
                         "must be a column from the table(s) joined in the "
                         "SQLAlchemy select.")

    values = []
    if in_column.type.python_type == bytes:
        source_values = basic.convert_to_encoded(source_values)

    chunked_values = basic.partition_list(source_values, limit)

    for value_chunk in chunked_values:
        subquery = executable.where(in_column.in_(value_chunk))

        proxy = engine.execute(subquery)
        results = proxy.fetchall()

        for result in results:
            if return_dict:
                result = dict(result)

            values.append(result)

    return values


def first_column_value_subqueries(engine, executable, in_column, source_values,
                                  limit=8000):
    """Query with a conditional on a set of values using subqueries.

    :param engine: SQLAlchemy Engine object used for executing queries.
    :type engine: Engine
    :param executable: Input a executable MySQL query.
    :type executable: Select
    :type executable: str
    :param in_column: SQLAlchemy Column object.
    :type in_column: Column
    :param source_values: Values from specified MySQL column.
    :type source_values: list[str]
    :param return_dict: Toggle whether to return data as a dictionary.
    :type return_dict: Boolean
    :param limit: SQLAlchemy IN clause query length limiter.
    :type limit: int
    :returns: Distinct values fetched from value constraints.
    :rtype: list
    """
    if not isinstance(in_column, Column):
        raise ValueError("Inputted column to conditional values against "
                         "is not a SqlAlchemy Column."
                         f"Object is instead type {type(in_column)}.")

    if not executable.is_derived_from(in_column.table):
        raise ValueError("Inputted column to conditional values against "
                         "must be a column from the table(s) joined in the "
                         "SQLAlchemy select.")

    values = []
    if in_column.type.python_type == bytes:
        source_values = basic.convert_to_encoded(source_values)

    chunked_values = basic.partition_list(source_values, limit)

    for value_chunk in chunked_values:
        subquery = executable.where(in_column.in_(value_chunk))

        proxy = engine.execute(subquery)
        results = proxy.fetchall()

        for result in results:
            values.append(result[0])

    values = list(OrderedDict.fromkeys(values))
    return values


def query(session, db_graph, table_map, where=None):
    """Use SQLAlchemy session to retrieve ORM objects from a mapped object.

    :param session: Bound and connected SQLAlchemy Session object.
    :type session: Session
    :param table_map: SQLAlchemy ORM map object.
    :param where: MySQL WHERE clause-related SQLAlchemy object(s).
    :type where: BinaryExpression
    :type where: list
    :param order_by: MySQL ORDER BY clause-related SQLAlchemy object(s).
    :type order_by: Column
    :type order_by: list
    :returns: List of mapped object instances.
    :rtype: list
    """
    where_columns = extract_columns(where)

    total_columns = [table_map] + where_columns
    from_clause = build_fromclause(db_graph, total_columns)

    query_obj = session.query(table_map).select_from(from_clause)

    if where is not None:
        if isinstance(where, list):
            query_obj = query_obj.filter(and_(*where))
        else:
            query_obj = query_obj.filter(where)

    instances = query_obj.all()
    return instances
