.. _mysql_connection:

Set up connection to MySQL
==========================

The primary tool to connect to MySQL and access databases is the ``pdm_utils`` AlchemyHandler. It stores login credentials, manages connections, and stores different types of information about the database and MySQL environment.

Create an AlchemyHandler object::

    >>> from pdm_utils.classes import alchemyhandler
    >>> alchemist = alchemyhandler.AlchemyHandler()

The AlchemyHandler can take credentials by setting its various attributes::


    >>> alchemist.username = "user123"
    >>> alchemist.password = "p@ssword"

Alternatively, the connect() method prompts for the credentials::

    >>> alchemist.connect()
    MySQL username:
    MySQL password:

Similarly, a specific database can be directly set::

    >>> alchemist.database = "Actinobacteriophage"

or indirectly using the connect() method::

    >>> alchemist.connect(ask_database=True)
    MySQL database:

Engine
------

If login credentials are valid, a SQLAlchemy engine is created and stored in the engine attribute. The engine object provides the core interface between Python and MySQL. It stores information about the database of interest, manages connections, and contains methods to directly interact with the database. With the engine, the database can be queried using pure SQL syntax::

    >>> engine = alchemist.engine
    >>> result = engine.execute("SELECT PhageID, HostGenus FROM phage WHERE Subcluster = 'A2'")
    >>> phages = result.fetchall()
    >>> len(phages)
    90
    >>> dict(phages[0])
    {'PhageID': '20ES', 'HostGenus': 'Mycobacterium'}


MySQL transactions can also be executed using the engine. It returns 0 if successful, or 1 if unsuccessful::

    >>> engine.execute("UPDATE phage SET HostGenus = 'Arthrobacter' WHERE PhageID = '20ES'")
    >>> result = engine.execute("SELECT PhageID,HostGenus FROM phage WHERE Subcluster = 'A2'")
    >>> phages = result.fetchall()
    >>> dict(phages[0])
    {'PhageID': '20ES', 'HostGenus': 'Arthrobacter'}


The ``pdm_utils`` 'mysqldb_basic' module provides several functions that rely on the engine, such as retrieving the list of tables in the database, or the list of columns in one of the tables::

    >>> from pdm_utils.functions import mysqldb_basic
    >>> tables = mysqldb_basic.get_tables(engine, alchemist.database)
    >>> tables
    {'gene_domain', 'phage', 'gene', 'pham', 'tmrna', 'version', 'trna', 'domain'}
    >>> columns = mysqldb_basic.get_columns(engine, alchemist.database, 'phage')
    >>> columns
    {'Length', 'Notes', 'Subcluster', 'Sequence', 'HostGenus', 'DateLastModified', 'RetrieveRecord', 'Cluster', 'Accession', 'AnnotationAuthor', 'GC', 'Status', 'Name', 'PhageID'}




Metadata
--------

A SQLAlchemy MetaData object can also be created from a AlchemyHandler object with a
valid credentials and database.  A SQLAlchemy MetaData object allows for access to 
SQLAlchemy Base Table and Object classes, where the Tables can be directly accessed 
from the metadata object and the Columns from a Table::

    >>> metadata = alchemist.metadata
    >>> phage = metadata.tables["phage"]
    >>> type(phage)
    <class 'sqlalchemy.sql.schema.Table'>

    >>> PhageID = phage.columns.PhageID
    >>> type(PhageID)
    <class 'sqlalchemy.sql.schema.Column'>

Table objects retrieved this way can be used to retrieve information about the 
structure of a table in the database the AlchemyHandler is connected to::

    >>> metadata = alchemist.metadata
    >>> list(metadata.tables)
    ['domain', 'gene', 'phage', 'pham', 'gene_domain', 'tmrna', 'trna', 'version']

    >>> phage = metadata.tables["phage"]
    >>> phage.name
    'phage'

    >>> primary_keys = phage.primary_key
    >>> type(primary_keys)
    <class 'sqlalchemy.sql.schema.PrimaryKeyConstraint'>

    >>> dict(primary_keys.columns).keys()
    dict_keys(['PhageID'])

Column objects retrieved from Table objects can be used to retrieve information about 
the characteristics of a Column within the table represented in the connected database::

    >>> PhageID = phage.columns.PhageID 
    >>> PhageID.name
    'PhageID'
    >>> str(PhageID)
    'phage.PhageID'

    >>> PhageID.type
    VARCHAR(length=25)
    >>> PhageID.nullable
    False 
    >>> PhageID.primary_key
    True

These Column and Table objects can be used to manually select, insert, or update data 
in a more object-oriented way when paired with an Engine object::

    >>> phage = alchemist.metadata.tables["phage"]
    >>> HostGenus = phage.columns.HostGenus
    >>> PhageID = phage.columns.PhageID
    >>> Subcluster = phage.columns.Subcluster

    >>> query = phage.select(Subcluster == 'A2')
    >>> result = alchemist.engine.execute(query)
    >>> phages = result.fetchall()
    >>> len(phages)
    90

    >>> from sqlalchemy import select

    >>> query = select([PhageID, Subcluster]).where(Subcluster == 'A2')
    >>> result = alchemist.engine.execute(query)
    >>> phages = result.fetchall()
    >>> len(phages)
    90

To query for information by indirect relationship conditionals, Tables and Columns can
be used to join tables to select from::

    >>> phage = alchemist.metadata.tables["phage"]
    >>> gene = alchemist.metadata.tables["gene"]

    >>> Cluster = phage.columns.Cluster
    >>> PhamID = gene.columns.PhamID

    >>> from sqlalchemy import join
    >>> joined_table = join(phage, gene, isouter=True)

    >>> from sqlalchemy import select
    >>> query = select([Cluster.distinct()]).select_from(joined_tables).where(PhamID == 2002)
    >>> result = alchemist.engine.execute(query)
    >>> clusters = result.fetchall()
    >>> dict(clusters[0])
    {'Cluster' : 'A'}

Database graph
--------------

An AlchemyHandler also has the ability to generate and store a graphical representation
 of the SQLAlchemy MetaData object as a ``NetworkX`` Graph object.  The graph object 
 has access to the same Table and objects as the MetaData as well as similar basic 
 information::

    >>> db_graph = alchemist.graph
    >>> list(db_graph.nodes)
    ['domain', 'gene', 'phage', 'pham', 'gene_domain', 'tmrna', 'trna', 'version']

    >>> phage_node = db_graph.nodes["phage"]
    >>> phage = phage_node["table"]
    >>> phage.name
    'phage'

The graph object also stores information about the relationships between two 
tables, specifically the foreign key constraints between tables (and if joining
 two tables is possible)::

    >>> from networkx import shortest_graph
    >>> db_graph = alchemist.graph
    >>> shortest_path(db_graph, 'phage', 'domain')
    ['phage', 'gene', 'gene_domain', 'domain']

    >>> foreign_key_edge = db_graph['phage']['gene']
    >>> foreign_key = foreign_key_edge["key"]

    >>> type(foreign_key)
    <class 'sqlalchemy.sql.schema.ForeignKey'>
    >>> foreign_key
    ForeignKey('phage.PhageID')


Mapper
------

The AlchemyHandler provides support for using the SQLAlchemy ORM module and 
SQLAlchemy ORM objects based on the schema of the connected database.  Access 
to the SQLAlchemy ORM objects is possible though the Automap Base object 
generated by the AlchemyHandler::
    
    >>> mapper = alchemist.mapper

    >>> Phage = mapper.classes["phage"]
    >>> type(Phage)
    <class 'sqlalchemy.ext.declarative.api.DeclarativeMeta'>

SQLAlchemy ORM objects have attributes that directly correspond to columns 
in the table that they represent, and these columns can be used in a similar 
way to the Base SQLAlchemy Column objects::

    >>> mapper = alchemist.mapper
    >>> Phage = mapper.classes["phage"]
    >>> conditional = Phage.Subcluster == 'A2'

    >>> Phage = alchemist.metadata.tables["phage"]
    >>> query = phage.select(conditional)
    >>> result = alchemist.engine.execute(query)
    >>> phages = result.fetchall()
    >>> len(phages)
    90

Session
-------

SQLAlchemy ORM objects are extremely powerful when used in combination with the
 SQLAlchemy session object.  The session object can be used with basic queries 
 to create objects that represent entries in the database that store information
 as attributes named after the columns in the table which the ORM object 
 represents.  In addition, the ORM object instances can be created, updated, 
 or deleted in a python environment, and the session object will manage and track the 
 changes::

    >>> session = alchemist.session
    >>> Phage = alchemist.mapper.classes["phage"]

    >>> phages = session.query(Phage).filter(Phage.Subcluster == 'A2')
    >>> len(phages)
    90
    
    >>> phage[0].PhageID
    '20ES'
    >>> phage[0].HostGenus
    'Mycobacterium'
