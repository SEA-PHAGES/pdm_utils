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




Database graph
--------------




Mapper
------


Session
-------
