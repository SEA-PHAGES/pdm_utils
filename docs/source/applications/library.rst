Library tutorial
================

``pdm_utils`` provides a library of functions, classes, and methods that can be used to develop customized data analysis pipelines. Below is a brief introduction to how the library can be used.

In the shell terminal, activate the Conda environment containing the installed
``pdm_utils`` package (if needed) to ensure all dependencies are present. Then open a Python IDE::

    > conda activate pdm_utils
    (pdm_utils)>
    (pdm_utils)> python3
    >>>



Connect to the MySQL database
*****************************

In order to retrieve and explore data stored within the MySQL database, create a SQLAlchemy 'engine'. This object provides the core interface between Python and MySQL. It stores information about the database of interest, login credentials and connection status so that they do not need to be provided multiple times, and it contains methods to directly interact with the database. The ``pdm_utils`` 'mysqldb' module provides several functions that rely on the engine. To start, create an engine to the primary database, and provide the username and password when prompted::

    >>> from pdm_utils.functions import mysqldb
    >>> engine = mysqldb.connect_to_db(database='Actinobacteriophage')
    MySQL username:
    MySQL password:

MySQL queries can be executed using the engine. In the following example, a list of 90 phages in Subcluster A2 are retrieved. For each phage, a dictionary of data is returned::

    >>> result = engine.execute("SELECT PhageID,HostGenus FROM phage WHERE Subcluster = 'A2'")
    >>> phages = result.fetchall()
    >>> len(phages)
    90
    >>> dict(phages[0])
    {'PhageID': '20ES', 'HostGenus': 'Mycobacterium'}

MySQL transactions can also be executed using the engine. It returns 0 if successful, or 1 if unsuccessful::

    >>> txn_result = mysqldb.execute_transaction(engine, ["UPDATE phage SET HostGenus = 'Arthrobacter' WHERE PhageID = '20ES'"])
    >>> txn_result
    0
    >>> result = engine.execute("SELECT PhageID,HostGenus FROM phage WHERE Subcluster = 'A2'")
    >>> phages = result.fetchall()
    >>> dict(phages[0])
    {'PhageID': '20ES', 'HostGenus': 'Arthrobacter'}



Access ``pdm_utils`` Genome data
********************************

Data can also be retrieved in an object-oriented structure. First, create a list of phages for which data should be retrieved. These are expected to be stored in the PhageID column of the *phage* table::

    >>> phage_id_list = ['L5', 'Trixie', 'D29']

Construct the MySQL query to retrieve the specific types of data from the *phage* table and the *gene* table::

    >>> phage_query = 'SELECT PhageID, Name, Sequence, Cluster, Subcluster, Status, HostGenus FROM phage'
    >>> gene_query = 'SELECT GeneID, Start, Stop, Orientation, Translation, Notes FROM gene'

The parse_genome_data function retrieves the data and constructs ``pdm_utils`` Genome and Cds objects from the data. In the example below, there are three Genome objects created, each corresponding to a different phage in phage_id_list::

    >>> phage_data = mysqldb.parse_genome_data(engine, phage_id_list=phage_id_list, phage_query=phage_query, gene_query=gene_query)
    >>> len(phage_data)
    3


Data for each phage can be directly accessed::

    >>> phage_data[0].id
    'D29'
    >>> d29 = phage_data[0]
    >>> d29.host_genus
    'Mycobacterium'
    >>> d29.cluster
    'A'
    >>> d29.subcluster
    'A2'
    >>> d29.annotation_status
    'final'

The genome sequence is stored in the seq attribute as a Biopython Seq object,
so Biopython Seq attributes and methods (such as 'lower' or 'reverse_complement') can also be directly accessed::

    >>> len(d29.seq)
    49136
    >>> d29.seq[:10]
    Seq('GGTCGGTTAT')
    >>> d29.seq[:10].lower()
    Seq('ggtcggttat')
    >>> d29.seq[:10].reverse_complement()
    Seq('ATAACCGACC')



Access ``pdm_utils`` Cds data
*****************************

Data from the *gene* table is retrieved and parsed into Cds objects.
For each phage, all Cds objects for are stored in the Genome object's 'cds_features' attribute as a list. Data for each CDS feature can be directly accessed::

    >>> len(d29.cds_features)
    77
    >>> cds54 = d29.cds_features[54]
    >>> cds54.description
    'DNA primase'
    >>> cds54.start
    38737
    >>> cds54.stop
    39127
    >>> cds54.orientation
    'R'
    >>> cds54.coordinate_format
    '0_half_open'


Similar to the nucleotide sequence in the Genome object, the CDS translation is stored in the translation attribute as a Biopython Seq object::

    >>> cds54.translation
    Seq('MTATGIAEVIQRYYPDWDPPPDHYEWNKCLCPFHGDETPSAAVSYDLQGFNCLA...PWS', IUPACProtein())


The nucleotide sequence for each Cds feature is not explicitly stored in the MySQL database. The sequence can be extracted from the parent genome, but this relies on the Cds object containing a Biopython SeqFeature object stored in the seqfeature attribute, but this is also empty at first::

    >>> cds54.seq
    Seq('', IUPACAmbiguousDNA())
    >>> cds54.seqfeature



To extract the sequence, first construct the Biopython SeqFeature object::

    >>> cds54.set_seqfeature()
    >>> cds54.seqfeature
    SeqFeature(FeatureLocation(ExactPosition(38737), ExactPosition(39127), strand=-1), type='CDS')

With the SeqFeature constructed, the 390 bp nucleotide sequence can be retrieved from the parent genome::

    >>> cds54.set_nucleotide_sequence(parent_genome_seq=d29.seq)
    >>> cds54.seq
    Seq('TTGACAGCCACCGGCATCGCGGAGGTCATCCAGCGGTACTACCCGGACTGGGAT...TGA')
    >>> len(cds54.seq)
    390


Connect to the MySQL database using a ``pdm_utils`` AlchemyHandler
******************************************************************
Many of the tools in pdm_utils require a ``pdm_utils`` AlchemyHandler object to manage connections and to abstract the relationships between the various SQLAlchemy package objects that allow for easy database access::
    
    >>> from pdm_utils.classes import alchemyhandler
    >>> alchemist = alchemyhandler.AlchemyHandler()

The AlchemyHandler can take credentials by setting its various attributes or by using the AlchemyHandler.connect() method to prompt for the credentials::

    >>> alchemist.username = "user123"
    >>> alchemist.password = "p@ssword"

    or

    >>> alchemist.connect()
    MySQL username:
    MySQL password:

The AlchemyHandler can connect to a specific database in a similar fashion, either through setting an attribute or using the AlchemyHandler.connect() method::

    >>> alchemist.database = "Actinobacteriophage"

    or

    >>> alchemist.connect(ask_database=True)
    MySQL database:

The AlchemyHandler then allows for access to a variety of SQLAlchemy objects related to the established connection::

    >>> engine = alchemist.engine
    >>> result = engine.execute("SELECT PhageID,HostGenus FROM phage WHERE Subcluster = 'A2'")
    >>> phages = result.fetchall()
    >>> len(phages)
    90
    >>> dict(phages[0])
    {'PhageID': '20ES', 'HostGenus': 'Mycobacterium'}

Access subsets of data using a ``pdm_utils`` Filter
***************************************************

Sometimes data pertaining to a large set of phages (for instance, all Subcluster A2 phages) is needed. Manually constructing the list of PhageIDs is time intensive and error prone, but can be automatically generated using a ``pdm_utils`` Filter object::

     >>> from pdm_utils.classes import filter
     >>> db_filter = filter.Filter()

A filter object can connect to a MySQL database by calling a method Filter.connect() or by passing the filter initialization a valid AlchemyHandler.

     >>> db_filter = filter.Filter()
     >>> db_filter.connect()

     or

     >>> db_filter = filter.Filter(alchemist=alchemy_handler_obj)

A filter object requires a key to retrieve values from, which can be set using the attribute Filter.key and passing in a MySQL formatted column, or a SQLAlchemy Column object.

    >>> db_filter.key = "phage.PhageID"

    or

    >>> db_filter.key = PhageID_Column_obj
    >>> db_filter.key
    Column('PhageID', VARCHAR(length=25), table=<phage>, primary_key=True, nullable=False)


A filter object key can also be set by passing in the name of a table as a string, and the filter key will be the primary key of that table.

    >>> db_filter.key = "phage"
    >>> db_filter.key
    Column('PhageID', VARCHAR(length=25), table=<phage>, primary_key=True, nullable=False)


The filter retrieves values from the database depending on previous retrieved values and given conditionals.  New conditionals can be applied using a method Filter.add() and MySQL syntax.
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


This list of PhageIDs can now be passed to other functions, such as mysqldb.parse_genome_data(). The filtered results can be filtered further if needed. Suppose that only Subcluster A2 phages that contain at least
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

OR conjunctions can be used with these statements, and OR conjunctions work in the same way that MySQL OR's do: AND conjunctions take precendance and are processed first before OR conjunctions::

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

To get the distinct set of values from multiple columns for each individual value in the set of values, the method Filter.retrieve() accepts a column input and returns a dictionary where the keys are the set of values and the valeus are dictionaries where the keys are the names of the entered columns, and the values are the respective distinct values in a list::

    >>> db_filter.values
    ['Trixie', 'Myrna', 'Alice']
    >>> db_filter.retrieve(["phage.Cluster", "phage.Subcluster"])
    {'Trixie' : {'Cluster' : ['A'], 'Subcluster' : ['A2']}, 'Myrna' : {'Cluster' : ['C'], 'Subcluster' : ['C2']}, 'Alice' : {'Cluster' : ['C'], 'Subcluster' : ['C1']}}

To group the current set of values based on the distinct set of values related to the current set of values, the method Filter.group() accepts a column input and returns the values in a dictionary where the keys are the set of related distinct columns and the values are the respective subset of the current set of values in a list::

    >>> db_filter.values
    ['Trixie', 'D29', 'Myrna', 'Alice']
    >>> db_filter.group("phage.Cluster")
    {'A' : ['Trixie', 'D29'], 'C' : ['Myrna', 'Alice']}


Dynamic MySQL joining and query construction
********************************************

The base of the ability of the Filter to dynamically accept a variety of column inputs from a variety of tables comes from a manipulation of SQLAlchemy base. This is possible though a NetworkX graph that describes the relationships and constraints between the tables of a database and allow for pathings between the tables to act as guides for joining database tables during MySQL query construction.  To create this graph, you can use a ``pdm_utils`` AlchemyHandler object connected to a MySQL database::

    >>> from pdm_utils.classes import alchemyhandler
    >>> alchemist = alchemyhandler.AlchemyHandler()

    >>> alchemist.username = "user123"
    >>> alchemist.password = "p@ssword"
    >>> alchemist.database = "Actinobacteriophage"

    >>> alchemist.connect()

    >>> graph = alchemist.graph

The graph can be used in conjunction with SQLAlchemy Column objects to dynamically create SQLAlchemy select objects that represent a MySQL statement.  SQLAlchemy Column and Table objects can be retrieved from a SQLAlchemy metadata object::

    >>> metadata = alchemist.metadata
    >>> phage_table = metadata.tables["phage"]
    >>> phageid_column = phage_table.columns["PhageID"]

Alternatively, SQLAlchemy Columns can be retrieved from case-insensitive MySQL formatted inputs using the ``pdm_utils`` 'querying' module::

    >>> from pdm_utils.functions import querying
    >>> metadata = alchemist.metadata
    >>> phageid = querying.get_column(metadata, "PHAGE.phageID")

Columns can be used with the ``pdm_utils`` 'querying' module and the populated graph object to create SQLAlchemy select statements that dynamically join the required tables::
    
    >>> from pdm_utils.functions import querying

    >>> cluster_column = querying.get_column(alchemist.metadata, "phage.cluster")
    >>> domain_name_column = querying.get_column(alchemist.metadata, "domain.name")
    >>> gene_notes_column = querying.get_column(alchemist.metadata, "gene.notes")

    >>> columns = [cluster_column, domain_name_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns)
    >>> results_proxy = alchemist.engine.execute(select_statement)

Columns can be used to create conditionals to narrow the results of a select statement, populating the MySQL WHERE clause of the statement.  The ``pdm_utils`` 'querying' module incoorporates the tables of the columns used in the conditionals when joining tables to allow for easier addition of conditionals with indirect relationships::

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

Select statements can be executed by the execute() function of the ``pdm_utils`` 'querying' module that allows for input of values to condition the query on and automatically retrieves the results from the results proxy. populating the MySQL IN clause of the statement.  The execute() comes with safer procedures for IN clause statement and automatically subqueries if the number of values exceeds the average query length limit for MySQL statements::

    >>> domain_name_column = querying.get_column(alchemist.metadata, "domain.name")
    >>> gene_notes_column = querying.get_column(alchemist.metadata, "gene.notes")
    >>> phageid_column = querying.get_column(alchemist.metadata, "phage.PhageID")

    >>> columns = [gene_notes_column, gene_notes_column]
    >>> select_statement = querying.build_select(alchemist.graph, columns)

    >>> results = querying.execute(alchemist.engine, select_statement, in_column=phageid_column, values=["Trixie", "D29", "Myrna", "Alice"])

Alternatively, since the ``pdm_utils`` Filter object wraps much of this functionality, the Filter object can replicate many of the same tasks, adding a layer of abstraction to the querying module::

    >>> db_filter.key = "phage.PhageID"
    >>> db_filter.values = ["Trixie", "D29", "Myrna", "Alice"]

    >>> results = db_filter.select(["domain.name", "gene.notes"])

Working with the SQLAlchemy ORM
*******************************

The SQLAlchemy ORM is a powerful tool for modifying and exploring a database.  Access to ORM mapped objects and mapped object instances is available through use of an ``pdm_utils`` AlchemyHandler object::
    
    >>> phage = alchemist.mapper.classes["phage"]

Alternatively, SQLAlchemy ORM mapped objects can be retrieved from case-insensitive MySQL formatted inputs using the ``pdm_utils`` 'cartography' module::

    >>> from pdm_utils.functions import cartography

    >>> phage = cartography.get_map(alchemist.mapper, "PHAGE")

These SQLAlchemy ORM mapped objects can be used for adding new entries, updating entries, and removing entries through the use of a session object.  A session object can be used to retrieve instances that are tracked by the session.  A session can be created using the ``pdm_utils`` AlchemyHandler object::

    >>> session = alchemist.session

Python objects, ORM map instances, that reflect entries in the database can be retrieved using the session object, and conditionals to filter the results of the query can be formed using attributes of the retrieved ORM map.::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")

    >>> Trixie = alchemist.session.query(phage).filter(phage.PhageID == 'Trixie').scalar()\

The retrieved instance is a reflection of a single entry tied to the primary_key from the mapped table specified, and has attributes that reflect the related columns of that table::
    >>> Trixie.PhageID
    'Trixie'
    >>> Trixie.Cluster
    'A'
    >>> Trixie.Length
    53526

The dynamic querying from ``pdm_utils`` 'querying' module can be applied to SQLAlchemy ORM queries using the query() function, and SQLAlchemy base objects and conditionals can be incoorporated from the querying module into ORM queries to generate ORM objects::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")

    >>> subcluster_conditional = phage.Subcluster == 'A2'
    >>> notes_conditional = querying.build_where_clause(alchemist.graph, "phage.Notes = 'antirepressor'")

    >>> conditionals = [subcluster_conditional, notes_conditional]

    >>> mapped_obj_instances = querying.query(alchemist.session, alchemist.graph, phage, where=conditionals)
    >>> phage_instance = mapped_obj_instances[0]
    >>> phage_instance.PhageID
    'IronMan'

Additionally, the ``pdm_utils`` Filter object can be used to retrieve these mapped instances.  The filter object can apply filters and retrieve a list of values that can be used to retrieve a similar set of mapped obj instances::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")

    >>> db_filter.add("phage.Subcluster = 'A2' AND gene.Notes = 'antirepressor'")
    >>> db_filter.update()

    >>> mapped_obj_instances = db_filter.query(phage)
    >>> phage_instance = mapped_obj_instances[0]
    >>> phage_instance.PhageID
    'IronMan'

The SQLAlchemy session tracks instances generated through these queries, and can be used to manually manage database information, and can be used to manually update entries in the database::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")

    >>> IronMan = alchemist.session.query(phage).filter(phage.PhageID == 'IronMan').scalar()
    >>> IronMan.DateLastModified
    datetime.datetime(2020, 3, 13, 0, 0)

    >>> from datetime import datetime
    >>> today = datetime(2020, 5, 22, 0, 0)

    >>> IronMan.DateLastModified = today
    >>> alchemist.session.commit()
    
    >>> IronMan = alchemist.session.query(phage).filter(phage.PhageID == 'IronMan').scalar()
    >>> IronMan.DateLastModified
    datetime.datetime(2020, 5, 22, 0, 0)

Once references to instances have been acquired using the session, entries in the database can also be deleted::

    >>> IronMan = alchemist.session.query(phage).filter(phage.PhageID == 'IronMan').scalar()
    >>> alchemist.session.delete(IronMan)

    >>> alchemist.session.query(phage).filter_by(PhageID='IronMan').count()
    0

The SQLAlchemy map can also be used to instantiate new objects that are then added as entries to the database::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")
    >>> Phabulous = phage(PhageID='Phabulous', Cluster='A', Subcluster='A2', Length=52342)

    >>> alchemist.session.commit()

    >>> alchemist.session.query(phage).filter_by(PhageID='Phabulous').count()
    1

SQLAlchemy mapped instances generated from a session also have access to the data that the relevant entry is has a relationship with::

    >>> phage = cartography.get_map(alchemist.mapper, "phage")
    >>> IronMan = alchemist.session.query(phage).filter_by(PhageID='IronMan').scalar()

    >>> IronMan_genes = IronMan.gene_collection
    >>> IronMan_gene1 = IronMan_genes[0]

    >>> IronMan_gene1.PhageID
    'IronMan' 
    >>> IronMan_gene1.Name
    1
    >>> IronMan_gene1.GeneID
    IronMan_CDS_1

These instances retrieved from the relationship attributes of another mapped instance can likewise be updated or deleted with use of the SQLAlchemy session::
    
    >>> IronMan = alchemist.session.query(phage).filter_by(PhageID='IronMan').scalar()
    >>> IronMan_gene1 = IronMan.gene_collection[0]

    >>> IronMan_gene1.PhamID
    42415
    >>> IronMan_gene1.PhamID = 54326

    >>> alchemist.session.commit()

When all interaction with MySQL is complete, the DBAPI connections can be closed::

    >>> engine.dispose()

For more information on how different Genome and Cds object attributes map to the MySQL database, refer to the :ref:`object attribute maps <attributemap>`.
