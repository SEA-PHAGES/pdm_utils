.. _filtering_querying:

Filtering and dynamic querying
==============================

By combining the :sqlalchemy:`SQLAlchemy <>` Core Metadata object with graphing tools from the :networkx:`NetworkX package <>`, ``pdm_utils`` can support dynamic querying of the database and filtering of data.

Dynamic MySQL joining and query construction
********************************************

The Filter can dynamically accept a column and table inputs since it combines SQLAlchemy Metadata with a :networkx:`NetworkX <>` graph. The graph describes the relationships and constraints between the tables of a database and allows for identifying the relational path between the tables. The graph thus acts as a guide for joining database tables during MySQL query construction. To create this graph, you can use a ``pdm_utils`` AlchemyHandler object connected to a MySQL database::

    >>> from pdm_utils.classes import alchemyhandler
    >>> alchemist = alchemyhandler.AlchemyHandler()
    >>> alchemist.username = "user123"
    >>> alchemist.password = "p@ssword"
    >>> alchemist.database = "Actinobacteriophage"
    >>> alchemist.connect()
    >>> graph = alchemist.graph

The graph can be used in conjunction with SQLAlchemy Column objects to dynamically create SQLAlchemy Select objects that represent a MySQL statement.  SQLAlchemy Column and Table objects can be retrieved from a SQLAlchemy Metadata object::

    >>> metadata = alchemist.metadata
    >>> phage_table = metadata.tables["phage"]
    >>> phageid_column = phage_table.columns["PhageID"]

Alternatively, SQLAlchemy Columns can be retrieved from case-insensitive MySQL formatted inputs using the ``pdm_utils`` 'querying' module::

    >>> from pdm_utils.functions import querying
    >>> metadata = alchemist.metadata
    >>> phageid = querying.get_column(metadata, "PHAGE.phageID")

Columns can be used with the ``pdm_utils`` 'querying' module and the graph object to create SQLAlchemy Select statements that dynamically join the required tables::

    >>> from pdm_utils.functions import querying
    >>> cluster_column = querying.get_column(alchemist.metadata, "phage.cluster")
    >>> domain_name_column = querying.get_column(alchemist.metadata, "domain.name")
    >>> gene_notes_column = querying.get_column(alchemist.metadata, "gene.notes")
    >>> columns = [cluster_column, domain_name_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns)
    >>> results_proxy = alchemist.engine.execute(select_statement)

Columns can be used to create conditionals to narrow the results of a select statement, populating the MySQL WHERE clause of the statement.  The ``pdm_utils`` 'querying' module incorporates the tables of the columns used in the conditionals when joining tables to allow for easier addition of conditionals with indirect relationships::

    >>> domain_name_column = querying.get_column(alchemist.metadata, "domain.name")
    >>> domain_name_conditional = (domain_name_column == "HTH_26")
    >>> columns = [cluster_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns, where=domain_name_conditional)
    >>> results_proxy = alchemist.engine.execute(select_statement)

Alternatively, SQLAlchemy ORM mapped objects can be retrieved from MySQL formatted conditional inputs using the ``pdm_utils`` 'querying' module build_where_clause() function and a graph object::

    >>> domain_name_conditional = querying.build_where_clause(alchemist.graph, "domain.name = 'HTH_26'")
    >>> columns = [cluster_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns, where=domain_name_conditional)
    >>> results_proxy = alchemist.engine.execute(select_statement)

Select statements can be executed by the execute() function of the ``pdm_utils`` 'querying' module. This allows for input of values to condition the query on and automatically retrieves the results from the results proxy. populating the MySQL IN clause of the statement.  The execute() comes with safer procedures for IN clause statement and automatically subqueries if the number of values exceeds the average query length limit for MySQL statements::

    >>> domain_name_column = querying.get_column(alchemist.metadata, "domain.name")
    >>> gene_notes_column = querying.get_column(alchemist.metadata, "gene.notes")
    >>> phageid_column = querying.get_column(alchemist.metadata, "phage.PhageID")
    >>> columns = [gene_notes_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns)
    >>> results = querying.execute(alchemist.engine, select_statement, in_column=phageid_column, values=["Trixie", "D29", "Myrna", "Alice"])


The ``pdm_utils`` Filter object wraps much of the querying functionality to add a layer of abstraction (see below for more details)::

    >>> db_filter.key = "phage.PhageID"
    >>> db_filter.values = ["Trixie", "D29", "Myrna", "Alice"]
    >>> results = db_filter.select(["domain.name", "gene.notes"])


Access subsets of data using a ``pdm_utils`` Filter
***************************************************

The ``pdm_utils`` Filter object enables iterative querying, sorting, and retrieving of subsets of data::

     >>> from pdm_utils.classes import filter
     >>> db_filter = filter.Filter()

A filter object can connect to a MySQL database by calling a method Filter.connect()::

     >>> db_filter = filter.Filter()
     >>> db_filter.connect()

Alternatively, a valid AlchemyHandler can be provided when the filter is instantiated::

     >>> db_filter = filter.Filter(alchemist=alchemy_handler_obj)

A filter object requires a key to retrieve values from, which can be set using the attribute Filter.key and passing in a MySQL formatted column::

    >>> db_filter.key = "phage.PhageID"

Alternatively, a SQLAlchemy Metadata Column object can be passed::

    >>> db_filter.key = PhageID_Column_obj
    >>> db_filter.key
    Column('PhageID', VARCHAR(length=25), table=<phage>, primary_key=True, nullable=False)


A filter object key can also be set by passing in the name of a table as a string, and the filter key will be the primary key of that table::

    >>> db_filter.key = "phage"
    >>> db_filter.key
    Column('PhageID', VARCHAR(length=25), table=<phage>, primary_key=True, nullable=False)


The filter retrieves values from the database depending on previously-retrieved values and given conditionals.  New conditionals can be applied using a method Filter.add() and MySQL syntax.
Ex. Creating the Subcluster filter identifies 90 phages in Subcluster A2::

     >>> db_filter.add("phage.Subcluster = 'A2'")
     >>> db_filter.update()
     >>> db_filter.hits()
     90

The filter results are stored in the values attribute, and can be sorted and accessed::

     >>> db_filter.sort("phage.PhageID")
     >>> len(db_filter.values)
     90
     >>> db_filter.values[:10]
     ['20ES', 'AbbyPaige', 'Acolyte', 'Adzzy', 'AN3', 'AN9', 'ANI8', 'AnnaL29', 'Anselm', 'ArcherNM']


This list of PhageIDs can now be passed to other functions, such as mysqldb.parse_genome_data(). The filtered results can be subsequently filtered if needed. Suppose that only Subcluster A2 phages that contain at least
one gene that is annotated as the 'repressor' are needed. This filter can be added, resulting in a list of only 4 phages::

     >>> db_filter.add("gene.Notes = 'repressor'")
     >>> db_filter.update()
     >>> db_filter.hits()
     4
     >>> db_filter.values
     ['Pukovnik', 'RedRock', 'Odin', 'Adzzy']

The same list of PhageIDs can be retrieved by adding conjunctions to the added conditionals::

    >>> db_filter.add("gene.Notes = 'repressor' AND phage.Subcluster = 'A2'")
    >>> db_filter.update()
    >>> db_filter.hits()
    4
    >>> db_filter.values
    ['Pukovnik', 'RedRock', 'Odin', 'Adzzy']

OR conjunctions can be used with these statements, and OR conjunctions work in the same way that MySQL OR's do: AND conjunctions take precedence and are processed first before OR conjunctions::

    >>> db_filter.add("gene.Notes = 'repressor' AND phage.Subcluster = 'A2' OR gene.Notes = 'repressor' AND phage.Subcluster = 'C1'")
    >>> db_filter.update()
    >>> db_filter.hits()
    5
    >>> db_filter.values
    ['Pukovnik', 'RedRock', 'Odin', 'Adzzy', 'LRRHood']

To get the distinct set of values from another column related to the current set of values, the method Filter.transpose() accepts a column input and returns the distinct values in a list::

    >>> db_filter.values
    ['Trixie', 'D29', 'Myrna', 'Alice']
    >>> db_filter.transpose("phage.Cluster")
    ['A', 'C']

Filter.transpose() has an optional parameter that allow you to switch a Filter's key and values::

    >>> db_filter.values
    ['Trixie', 'D29', 'Myrna', 'Alice']
    >>> db_filter.key
    Column('PhageID', VARCHAR(length=25), table=<phage>, primary_key=True, nullable=False)
    >>> db_filter.transpose("phage.Cluster", set_values=True)
    ['A', 'C']
    >>> db_filter.values
    ['A', 'C']
    >>> db_filter.key
    Column('Cluster', VARCHAR(length=5), table=<phage>)

To get the distinct set of values from multiple columns related to the current set of values, the method Filter.mass_transpose() accepts column input(s) and returns the values in a dictionary where the keys are the names of the entered columns, and the values are the respective distinct values in a list::

    >>> db_filter.values
    ['Trixie', 'D29', 'Myrna', 'Alice']
    >>> db_filter.mass_transpose(["phage.Cluster",  "phage.Subcluster"])
    {'Cluster' : ['A', 'C'], 'Subcluster' : ['A2', 'C1', 'C2']}

To get the distinct set of values from multiple columns for each individual value in the set of values, the method Filter.retrieve() accepts a column input and returns a dictionary where the keys are the set of values and the values are dictionaries where the keys are the names of the entered columns, and the values are the respective distinct values in a list::

    >>> db_filter.values
    ['Trixie', 'Myrna', 'Alice']
    >>> db_filter.retrieve(["phage.Cluster", "phage.Subcluster"])
    {'Trixie' : {'Cluster' : ['A'], 'Subcluster' : ['A2']}, 'Myrna' : {'Cluster' : ['C'], 'Subcluster' : ['C2']}, 'Alice' : {'Cluster' : ['C'], 'Subcluster' : ['C1']}}

To group the current set of values based on the distinct set of values related to the current set of values, the method Filter.group() accepts a column input and returns the values in a dictionary where the keys are the set of related distinct columns and the values are the respective subset of the current set of values in a list::

    >>> db_filter.values
    ['Trixie', 'D29', 'Myrna', 'Alice']
    >>> db_filter.group("phage.Cluster")
    {'A' : ['Trixie', 'D29'], 'C' : ['Myrna', 'Alice']}
