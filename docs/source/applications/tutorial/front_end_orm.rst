.. _front_end_orm:

Front-end ORM
=============


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
