.. _back_end_orm:

Back-end ORM
============

A collection of customized classes represent the 'back-end' biology-centric ORM that is used to parse, evaluate, and exchange biological data between different data sources, including GenBank, PhagesDB, and MySQL. The objects leverage the highly-refined and well-supported :biopython:`BioPython package <>` to coordinate data exchange as well as perform biology-related tasks (sequence extraction, translation, etc.). In contrast to the SQLAlchemy 'front-end' ORM, this ORM is only compatible with the most current ``pdm_utils`` database schema.

For information on how different classes in the bio-centric ORM map to different types of data sources, including the MySQL database, refer to the :ref:`bio-centric object relational mappings <attributemap>`.


Genome data
***********

Data for specific phage genomes can be retrieved in an object-oriented structure using the 'mysqldb' module. First, create a list of phages for which data should be retrieved. These are expected to be stored in the PhageID column of the *phage* table::

    >>> from pdm_utils.functions import mysqldb
    >>> phage_id_list = ['L5', 'Trixie', 'D29']

Construct the MySQL query to retrieve the specific types of data from the *phage* table and the *gene* table::

    >>> phage_query = 'SELECT PhageID, Name, Sequence, Cluster, Subcluster, Status, HostGenus FROM phage'
    >>> gene_query = 'SELECT GeneID, Start, Stop, Orientation, Translation, Notes FROM gene'

The mysqldb.parse_genome_data() function retrieves the data and constructs ``pdm_utils`` Genome, Cds, Trna, and Tmrna objects from the data. In the example below, there are three Genome objects created, each corresponding to a different phage in phage_id_list::

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



Cds data
********

Data from the *gene* table is retrieved and parsed into Cds objects.
For each phage, all Cds objects are stored in the Genome object's 'cds_features' attribute as a list. Data for each CDS feature can be directly accessed::

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



Trna and Tmrna data
*******************

Similar to CDS data, data from the *trna* and *tmrna* tables are retrieved and parsed into Trna and Tmrna objects, and stored in the Genome.trna_features and Genome.tmrna_features attributes, respectively. Each class contains a variety of methods to validate and manipulate this type of data.



Source data
***********

Similar to CDS data, data from source features in GenBank-formatted flat files are parsed into Source objects. There is no equivalent source table in the database, but the class contains a variety of methods to validate and manipulate this type of data from flat files.
