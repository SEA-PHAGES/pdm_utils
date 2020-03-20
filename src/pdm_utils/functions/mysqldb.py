"""Functions to interact with MySQL."""

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import getpass
import subprocess
import sys
import sqlalchemy
from pdm_utils.classes import genome
from pdm_utils.classes import genomepair
from pdm_utils.classes import cds
from pdm_utils.constants import constants
from pdm_utils.functions import basic



# TODO unittest.
def query_dict_list(engine, query):
    """Get the results of a MySQL query as a list of dictionaries.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param query: MySQL query statement.
    :type query: str
    :returns:
        List of dictionaries, where each dictionary represents a row of data.
    :rtype: list
    """
    result_list = engine.execute(query).fetchall()
    result_dict_list = []
    for row in result_list:
        row_as_dict = dict(row)
        result_dict_list.append(row_as_dict)
    return result_dict_list

def parse_phage_table_data(data_dict, trans_table=11, gnm_type=""):
    """Parse a MySQL database dictionary to create a Genome object.

    :param data_dict:
        Dictionary of data retrieved from the phage table.
    :type data_dict: dict
    :param trans_table:
        The translation table that can be used to translate CDS features.
    :type trans_table: int
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :returns: A pdm_utils genome object.
    :rtype: genome
    """

    gnm = genome.Genome()
    try:
        gnm.id = data_dict["PhageID"]
    except:
        pass

    try:
        gnm.accession = data_dict["Accession"]
    except:
        pass

    try:
        gnm.name = data_dict["Name"]
    except:
        pass

    try:
        gnm.host_genus = data_dict["HostGenus"]
    except:
        pass

    try:
        # Sequence data is stored as MEDIUMBLOB, so decode to string.
        gnm.set_sequence(data_dict["Sequence"].decode("utf-8"))
    except:
        pass

    try:
        gnm.length = int(data_dict["Length"])
    except:
        pass

    try:
        # DateLastModified gets returned as a datetime.datetime object.
        # TODO some phages have no date, so it will be returned NULL.
        gnm.date = data_dict["DateLastModified"]
    except:
        pass

    try:
        gnm.description = data_dict["Notes"].decode("utf-8")
    except:
        pass

    try:
        gnm.gc = float(data_dict["GC"])
    except:
        pass

    try:
        # Singletons are stored in the MySQL database as NULL, which gets
        # returned as None.
        gnm.set_cluster(data_dict["Cluster"])
    except:
        pass


    # Non-subclustered phages are stored in the MySQL database as NULL, which gets
    # returned as None.
    try:
        gnm.set_subcluster(data_dict["Subcluster"])
    except:
        pass

    try:
        gnm.annotation_status = data_dict["Status"]
    except:
        pass

    try:
        gnm.set_retrieve_record(data_dict["RetrieveRecord"])
    except:
        pass

    try:
        gnm.set_annotation_author(data_dict["AnnotationAuthor"])
    except:
        pass

    gnm.translation_table = trans_table
    gnm.type = gnm_type
    return gnm


def parse_gene_table_data(data_dict, trans_table=11):
    """Parse a MySQL database dictionary to create a Cds object.

    :param data_dict:
        Dictionary of data retrieved from the gene table.
    :type data_dict: dict
    :param trans_table:
        The translation table that can be used to translate CDS features.
    :type trans_table: int
    :returns: A pdm_utils Cds object.
    :rtype: Cds
    """

    cds_ftr = cds.Cds()
    cds_ftr.type = "CDS"
    try:
        cds_ftr.id = data_dict["GeneID"]
    except:
        pass

    try:
        cds_ftr.genome_id = data_dict["PhageID"]
    except:
        pass

    try:
        cds_ftr.start = int(data_dict["Start"])
    except:
        pass

    try:
        cds_ftr.stop = int(data_dict["Stop"])
    except:
        pass

    try:
        cds_ftr.parts = int(data_dict["Parts"])
    except:
        pass

    cds_ftr.coordinate_format = "0_half_open"

    try:
        cds_ftr.length = int(data_dict["Length"])
    except:
        pass

    try:
        cds_ftr.name = data_dict["Name"]
    except:
        pass

    try:
        cds_ftr.set_translation(data_dict["Translation"])
    except:
        pass

    try:
        cds_ftr.orientation = data_dict["Orientation"]
    except:
        pass

    try:
        cds_ftr.description = data_dict["Notes"].decode("utf-8")
    except:
        pass

    try:
        cds_ftr.set_locus_tag(data_dict["LocusTag"])
    except:
        pass

    try:
        cds_ftr.pham_id = int(data_dict["PhamID"])
    except:
        pass

    try:
        cds_ftr.domain_status = int(data_dict["DomainStatus"])
    except:
        pass

    try:
        cds_ftr.translation_table = trans_table
    except:
        pass
    return cds_ftr


def retrieve_data(engine, column=None, query=None, phage_id_list=None):
    """Retrieve genome data from a MySQL database for a single genome.

    The query is modified to include one or more PhageIDs

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param query:
        A MySQL query that selects valid, specific columns
        from the a valid table without conditioning on a PhageID
        (e.g. 'SELECT PhageID, Cluster FROM phage').
    :type query: str
    :param column:
        A valid column in the table upon which the query can be conditioned.
    :type column: str
    :param phage_id_list:
        A list of valid PhageIDs upon which the query can be conditioned.
        In conjunction with the 'column' parameter, the 'query' is
        modified (e.g. "WHERE PhageID IN ('L5', 'Trixie')").
    :type phage_id_list: list
    :returns:
        A list of items, where each item is a dictionary of
        SQL data for each PhageID.
    :rtype: list
    """
    if (phage_id_list is not None and len(phage_id_list) > 0):
        query = query \
                + f" WHERE {column} IN ('" \
                + "','".join(phage_id_list) \
                + "')"
    query = query + ";"
    result_dict_list = query_dict_list(engine, query)
    return result_dict_list


def parse_cds_data(engine, column=None, phage_id_list=None, query=None):
    """Returns Cds objects containing data parsed from a
    MySQL database.

    :param engine:
        This parameter is passed directly to the 'retrieve_data' function.
    :type engine: Engine
    :param query:
        This parameter is passed directly to the 'retrieve_data' function.
    :type query: str
    :param column:
        This parameter is passed directly to the 'retrieve_data' function.
    :type column: str
    :param phage_id_list:
        This parameter is passed directly to the 'retrieve_data' function.
    :type phage_id_list: list
    :returns: A list of pdm_utils Cds objects.
    :rtype: list
    """
    cds_list = []
    result_list = retrieve_data(
                    engine, column=column, query=query,
                    phage_id_list=phage_id_list)
    for data_dict in result_list:
        cds_ftr = parse_gene_table_data(data_dict)
        cds_list.append(cds_ftr)
    return cds_list


def parse_genome_data(engine, phage_id_list=None, phage_query=None,
                      gene_query=None, trna_query=None, gnm_type=""):
    """Returns a list of Genome objects containing data parsed from a MySQL
    database.

    :param engine:
        This parameter is passed directly to the 'retrieve_data' function.
    :type engine: Engine
    :param phage_query:
        This parameter is passed directly to the 'retrieve_data' function
        to retrieve data from the phage table.
    :type phage_query: str
    :param gene_query:
        This parameter is passed directly to the 'parse_cds_data' function
        to retrieve data from the gene table.
        If not None, pdm_utils Cds objects for all of the phage's
        CDS features in the gene table will be constructed
        and added to the Genome object.
    :type gene_query: str
    :param trna_query:
        This parameter is passed directly to the '' function
        to retrieve data from the tRNA table. Note: not yet implemented.
        If not None, pdm_utils Trna objects for all of the phage's
        CDS features in the gene table will be constructed
        and added to the Genome object.
    :type trna_query: str
    :param phage_id_list:
        This parameter is passed directly to the 'retrieve_data' function.
        If there is at at least one valid PhageID, a pdm_utils genome
        object will be constructed only for that phage. If None, or an
        empty list,  genome objects for all phages in the
        database will be constructed.
    :type phage_id_list: list
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :returns: A list of pdm_utils Genome objects.
    :rtype: list
    """
    genome_list = []
    result_list1 = retrieve_data(engine, column="PhageID",
                                 phage_id_list=phage_id_list,
                                 query=phage_query)
    for data_dict in result_list1:
        gnm = parse_phage_table_data(data_dict, gnm_type=gnm_type)
        if gene_query is not None:
            cds_list = parse_cds_data(engine, column="PhageID",
                                      phage_id_list=[gnm.id],
                                      query=gene_query)

            for x in range(len(cds_list)):
                cds_list[x].genome_length = gnm.length
            gnm.cds_features = cds_list
        if trna_query is not None:
            # TODO develop this step once tRNA table and objects are built.
            pass
        genome_list.append(gnm)
    return genome_list




def get_distinct_data(engine, table, column, null=None):
    """Get set of distinct values currently in a MySQL database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param table: A valid table in the database.
    :type table: str
    :param column: A valid column in the table.
    :type column: str
    :param null: Replacement value for NULL data.
    :type null: misc
    :returns: A set of distinct values from the database.
    :rtype: set
    """
    query = f"SELECT DISTINCT({column}) FROM {table}"
    result_set = query_set(engine, query)

    if None in result_set:
        result_set.remove(None)
        result_set.add(null)

    return result_set


def create_seq_set(engine):
    """Create set of genome sequences currently in a MySQL database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: A set of unique values from phage.Sequence.
    :rtype: set
    """
    query = "SELECT Sequence FROM phage"

    # Returns a list of items, where each item is a tuple of
    # SQL data for each row in the table.
    result_list = engine.execute(query).fetchall()

    # Convert to a set of sequences.
    # Sequence data is stored as MEDIUMBLOB, so data is returned as bytes
    # "b'AATT", "b'TTCC", etc.
    result_set = set()
    for tup in result_list:
        gnm_seq = tup[0].decode("utf-8")
        gnm_seq = Seq(gnm_seq, IUPAC.ambiguous_dna).upper()
        result_set.add(gnm_seq)
    return result_set


def convert_for_sql(value, check_set=set(), single=True):
    """Convert a value for inserting into MySQL.

    :param value: Value that should be checked for conversion.
    :type value: misc
    :param check_set: Set of values to check against.
    :type check_set: set
    :param single: Indicates whether single quotes should be used.
    :type single: bool
    :returns:
        Returns either "NULL" or the value encapsulated in quotes
        ("'value'" or '"value"')
    :rtype: str
    """
    if value in check_set:
        value = "NULL"
    else:
        if single == True:
            value = f"'{value}'"
        else:
            value = f'"{value}"'
    return value


def create_update(table, field2, value2, field1, value1):
    """Create MySQL UPDATE statement.

    "'UPDATE <table> SET <field2> = '<value2' WHERE <field1> = '<data1>'."

    When the new value to be added is 'singleton' (e.g. for Cluster
    fields), or an empty value (e.g. None, "none", etc.),
    the new value is set to NULL.

    :param table: The database table to insert information.
    :type table: str
    :param field1: The column upon which the statement is conditioned.
    :type field1: str
    :param value1:
        The value of 'field1' upon which the statement is conditioned.
    :type value1: str
    :param field2: The column that will be updated.
    :type field2: str
    :param value2:
        The value that will be inserted into 'field2'.
    :type value2: str
    :returns: A MySQL UPDATE statement.
    :rtype: set
    """
    check_set = constants.EMPTY_SET | {"Singleton"}
    part1 = f"UPDATE {table} SET {field2} = "
    part3 = f" WHERE {field1} = '{value1}';"
    part2 = convert_for_sql(value2, check_set=check_set, single=True)
    statement = part1 + part2 + part3
    return statement


def create_delete(table, field, data):
    """Create MySQL DELETE statement.

    "'DELETE FROM <table> WHERE <field> = '<data>'."

    :param table: The database table to insert information.
    :type table: str
    :param field: The column upon which the statement is conditioned.
    :type field: str
    :param data: The value of 'field' upon which the statement is conditioned.
    :type data: str
    :returns: A MySQL DELETE statement.
    :rtype: str
    """
    statement = f"DELETE FROM {table} WHERE {field} = '{data}';"
    return statement


def create_gene_table_insert(cds_ftr):
    """Create a MySQL gene table INSERT statement.

    :param cds_ftr: A pdm_utils Cds object.
    :type cds_ftr: Cds
    :returns:
        A MySQL statement to INSERT a new row in the 'gene' table
        with data for several fields.
    :rtype: str
    """
    locus_tag = convert_for_sql(cds_ftr.locus_tag, check_set={""}, single=False)

    # cds_ftr.translation is a BioPython Seq object.
    # It is coerced to string by default.
    # Statement should be triple quoted.
    # Descriptions may contain a "'", so at least description attribute needs
    # to be encapsulated with double quotes.
    # locus_tag can be NULL or can be a string, so it should not be
    # encapsulated with '"'.
    statement = ("""INSERT INTO gene """
                 """(GeneID, PhageID, Start, Stop, Length, Name, """
                 """Translation, Orientation, Notes, LocusTag, Parts) """
                 """VALUES """
                 """("{}", "{}", {}, {}, {}, "{}", "{}", "{}", "{}", {}, {});""")
    statement = statement.format(cds_ftr.id, cds_ftr.genome_id,
                                 cds_ftr.start, cds_ftr.stop,
                                 cds_ftr.translation_length,
                                 cds_ftr.name, cds_ftr.translation,
                                 cds_ftr.orientation, cds_ftr.description,
                                 locus_tag, cds_ftr.parts)
    return statement


def create_phage_table_insert(gnm):
    """Create a MySQL phage table INSERT statement.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :returns:
        A MySQL statement to INSERT a new row in the 'phage' table
        with data for several fields.
    :rtype: str
    """
    cluster = convert_for_sql(gnm.cluster, check_set={"Singleton"}, single=True)
    subcluster = convert_for_sql(gnm.subcluster, check_set={"none"}, single=True)

    # gnm.seq is a BioPython Seq object.
    # It is coerced to string by default.
    statement = ("INSERT INTO phage "
                 "(PhageID, Accession, Name, HostGenus, Sequence, "
                 "Length, GC, Status, DateLastModified, "
                 "RetrieveRecord, AnnotationAuthor, "
                 "Cluster, Subcluster) "
                 "VALUES "
                 f"('{gnm.id}', '{gnm.accession}', '{gnm.name}', "
                 f"'{gnm.host_genus}', '{gnm.seq}', "
                 f"{gnm.length}, {gnm.gc}, '{gnm.annotation_status}', "
                 f"'{gnm.date}', {gnm.retrieve_record}, "
                 f"{gnm.annotation_author}, "
                 f"{cluster}, {subcluster});"
                 )
    return statement


def create_genome_statements(gnm, tkt_type=""):
    """Create list of MySQL statements based on the ticket type.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :param tkt_type: 'add' or 'replace'.
    :type tkt_type: str
    :returns:
        List of MySQL statements to INSERT all data
        from a genome into the database
        (DELETE FROM genome, INSERT INTO phage, INSERT INTO gene, ...).
    :rtype: list
    """

    sql_statements = []
    if tkt_type == "replace":
        statement1 = create_delete("phage", "PhageID", gnm.id)
        sql_statements.append(statement1)
    statement2 = create_phage_table_insert(gnm)
    sql_statements.append(statement2)
    for cds_ftr in gnm.cds_features:
        statement3 = create_gene_table_insert(cds_ftr)
        sql_statements.append(statement3)

    # TODO add steps to insert tRNA and tmRNA data.

    return sql_statements


def get_phage_table_count(engine):
    """Get the current number of genomes in the database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: Number of rows from the phage table.
    :rtype: int
    """
    query = "SELECT COUNT(*) FROM phage"
    result_list = engine.execute(query).fetchall()
    count = result_list[0][0]
    return count


def change_version(engine, amount=1):
    """Change the database version number.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param amount: Amount to increment/decrement version number.
    :type amount: int
    """
    result = get_version_table_data(engine)
    current = result["Version"]
    new = current + amount
    print(f"Updating version from {current} to {new}.")
    statement = (f"UPDATE version SET Version = {new}")
    engine.execute(statement)


# TODO originally coded in export pipeline, so ensure that function is removed.
# TODO unittest.
def get_version_table_data(engine):
    """Retrieves data from the version table.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: Dictionary containing keys "Version" and "SchemaVersion".
    :rtype: dict
    """
    query = "SELECT * FROM version"
    result_dict_list = query_dict_list(engine, query)
    return result_dict_list[0]


# TODO unittest.
def get_mysql_dbs(engine):
    """Retrieve database names from MySQL.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: Set of database names.
    :rtype: set
    """
    query = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA"
    databases = query_set(engine, query)
    return databases


# TODO unittest.
def get_db_tables(engine, database):
    """Retrieve tables names from the database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: Set of table names.
    :rtype: set
    """
    query = ("SELECT table_name FROM information_schema.tables "
             f"WHERE table_schema = '{database}'")
    db_tables = query_set(engine, query)
    return db_tables


# TODO unittest.
def get_table_columns(engine, database, table_name):
    """Retrieve columns names from a table.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param database: Name of the database to query.
    :type database: str
    :param table_name: Name of the table to query.
    :type table_name: str
    :returns: Set of column names.
    :rtype: set
    """
    query = ("SELECT column_name FROM information_schema.columns WHERE "
              f"table_schema = '{database}' AND "
              f"table_name = '{table_name}'")
    columns = query_set(engine, query)
    return columns



def query_set(engine, query):
    """Retrieve set of data from MySQL query.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param query: MySQL query statement.
    :type query: str
    :returns: Set of queried data.
    :rtype: set
    """
    result_list = engine.execute(query).fetchall()
    set_of_data = basic.get_values_from_tuple_list(result_list)
    return set_of_data


# TODO unittest.
def get_schema_version(engine):
    """Identify the schema version of the database_versions_list.

    Schema version data has not been persisted in every schema version,
    so if schema version data is not found, it is deduced from other
    parts of the schema:
    1. If the version table does not exist, schema_version = 0.
    2. If there is no schema_version or SchemaVersion field,
        it is either schema_version = 1 or 2.
    3. If AnnotationAuthor, Program, AnnotationQC, and RetrieveRecord
        columns are in phage table, schema_version = 2.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: The version of the pdm_utils database schema.
    :rtype: int
    """
    db_tables = get_db_tables(engine, engine.url.database)
    if "version" in db_tables:
        version_table = True
    else:
        version_table = False

    if version_table == True:
        version_columns = get_version_table_data(engine)
        if "schema_version" in version_columns.keys():
            schema_version = version_columns["schema_version"]
        elif "SchemaVersion" in version_columns.keys():
            schema_version = version_columns["SchemaVersion"]
        else:
            phage_columns = get_table_columns(
                                engine, engine.url.database, "phage")
            expected = {"AnnotationAuthor", "Program",
                        "AnnotationQC", "RetrieveRecord"}
            diff = expected - phage_columns
            if len(diff) == 0:
                schema_version = 2
            else:
                schema_version = 1
    else:
        schema_version = 0
    return schema_version

# TODO unittest.
def check_schema_compatibility(engine, pipeline, code_version=None):
    """Confirm database schema is compatible with code.

    If schema version is not compatible, sys.exit is called.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param pipeline: Description of the pipeline checking compatibility.
    :type pipeline: str
    :param code_version:
        Schema version on which the pipeline operates. If no schema version is
        provided, the package-wide schema version value is used.
    :type code_version: int
    """
    schema_version = get_schema_version(engine)
    if code_version is None:
        code_version = constants.CODE_SCHEMA_VERSION
    if code_version != schema_version:
        print(f"The database schema version is {schema_version}, but "
              f"{pipeline} is compatible with "
              f"schema version {code_version}.")
        sys.exit(1)


# TODO unittest.
def drop_create_db(engine, database):
    """Creates a new, empty database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param database: Name of the database to drop and create.
    :type database: str
    :returns: Indicates if drop/create was successful (0) or failed (1).
    :rtype: int
    """
    # First, test if the database already exists within mysql.
    # If there is, delete it so that a new database is installed.
    databases = get_mysql_dbs(engine)
    if database in databases:
        result = drop_db(engine, database)
    else:
        result = 0
    if result == 0:
        result = create_db(engine, database)
    return result


# TODO unittest.
def drop_db(engine, database):
    """Delete a database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param database: Name of the database to drop.
    :type database: str
    :returns: Indicates if drop was successful (0) or failed (1).
    :rtype: int
    """
    statement = f"DROP DATABASE {database}"
    try:
        engine.execute(statement)
        return 0
    except:
        return 1

# TODO unittest.
def create_db(engine, database):
    """Create a new, empty database.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param database: Name of the database to create.
    :type database: str
    :returns: Indicates if create was successful (0) or failed (1).
    :rtype: int
    """
    statement = f"CREATE DATABASE {database}"
    try:
        engine.execute(statement)
        return 0
    except:
        return 1


# TODO unittest.
def copy_db(engine, new_database):
    """Copies a database.

    :param engine:
        SQLAlchemy Engine object able to connect to a MySQL database, which
        contains the name of the database that will be copied into
        the new database.
    :type engine: Engine
    :param new_database: Name of the new copied database.
    :type new_database: str
    :returns: Indicates if copy was successful (0) or failed (1).
    :rtype: int
    """
    #mysqldump -u root -pPWD database1 | mysql -u root -pPWD database2
    command_string1 = ("mysqldump "
                      f"-u {engine.url.username} "
                      f"-p{engine.url.password} "
                      f"{engine.url.database}")
    command_string2 = ("mysql -u "
                      f"{engine.url.username} "
                      f"-p{engine.url.password} "
                      f"{new_database}")
    command_list1 = command_string1.split(" ")
    command_list2 = command_string2.split(" ")
    try:
        print("Copying database...")

        # Per subprocess documentation:
        # 1. For pipes, use Popen instead of check_call.
        # 2. Call p1.stdout.close() to allow p1 to receive a SIGPIPE if p2 exits.
        #    which gets called when used as a context manager.
        # communicate() waits for the process to complete.
        with subprocess.Popen(command_list1, stdout=subprocess.PIPE) as p1:
            with subprocess.Popen(command_list2, stdin=p1.stdout) as p2:
                p2.communicate()
        print("Copy complete.")
        result = 0
    except:
        print(f"Unable to copy {engine.url.database} to {new_database} in MySQL.")
        result = 1
    return result


def connect_to_db(database):
    """Connect to a MySQL database.

    :param database: Name of the database to connect to.
    :type database: str
    :returns:
        SQLAlchemy Engine object able to connect to a MySQL database.
        If unable to connect to the database, sys.exit is called.
    :rtype: Engine
    """
    engine, msg = get_engine(database=database, echo=False)
    if engine is None:
        print(msg)
        sys.exit(1)
    else:
        return engine


# TODO unittest.
def install_db(engine, schema_filepath):
    """Install a MySQL file into the indicated database.

    :param engine:
        SQLAlchemy Engine object able to connect to a MySQL databas.
    :type engine: Engine
    :param schema_filepath: Path to the MySQL database file.
    :type schema_filepath: Path
    """
    command_string = (f"mysql -u {engine.url.username} "
                      f"-p{engine.url.password} {engine.url.database}")
    command_list = command_string.split(" ")
    with schema_filepath.open("r") as fh:
        try:
            print("Installing database...")
            subprocess.check_call(command_list, stdin=fh)
            print("Installation complete.")
        except:
            print(f"Unable to install {schema_filepath.name} in MySQL.")


# TODO function is to replace MySQLConnectionHandler usage.
def get_engine(username=None, password=None, database=None,
               echo=True, attempts=5):
    """Create SQLAlchemy Engine object.

    If username, password, or database is None, user is
    prompted to provide values. To connect to MySQL without specifying
    a database, set database="".

    :param username: Username to login to MySQL.
    :type username: str
    :param password: Password to login to MySQL.
    :type password: str
    :param database: Name of the database to connect to.
    :type database: str
    :param echo:
        This parameter is passed directly to sqlalchemy.create_engine().
    :type echo: bool
    :param attempts: Number of login attempts allowed, before returning.
    :type attempts: int
    :returns:
        tuple (SQLAlchemy Engine, message)
        WHERE
        SQLAlchemy Engine(Engine) is a SQLAlchemy Engine object
        able to connect to a MySQL database.
        message(str) is the result of the attempt to create the engine.
    :rtype: tuple
    """
    attempt = 0
    valid = False
    msg = "Setting up MySQL connection. "
    while (attempt < attempts and valid == False):
        # The value provided by getpass needs to be distinct from
        # the function parameter. With this setup, if the user inputs
        # an incorrect value using getpass, they can try again.
        # But if the user passes an incorrect value as a parameter,
        # getpass is not called.
        if username is None:
            input_username = getpass.getpass(prompt="MySQL username: ")
        else:
            input_username = username
        if password is None:
            input_password = getpass.getpass(prompt="MySQL password: ")
        else:
            input_password = password
        if database is None:
            input_database = input("Database: ")
        else:
            input_database = database
        engine_string = construct_engine_string(
                            username=input_username,
                            password=input_password,
                            database=input_database)
        engine = sqlalchemy.create_engine(engine_string, echo=echo)
        try:
            conn = engine.connect()
            conn.close()
            valid = True
            msg = msg + "Valid MySQL login credentials. "
        except:
            valid = False
            print("Invalid login credentials.")
        attempt += 1

    if valid != True:
        engine = None
        msg = msg + "Invalid MySQL login credentials. "
        if attempt == attempts:
            msg = msg + (f"For security purposes, only {attempts} MySQL "
                         "login attempt(s) are permitted at once. "
                         "Please verify your login credentials and "
                         "database name, and try again.")
    return (engine, msg)


def construct_engine_string(db_type="mysql", driver="pymysql",
                            username="", password="", database=""):
    """Construct a SQLAlchemy engine URL.

    :param db_type: Type of SQL database.
    :type db_type: str
    :param driver: Name of the Python DBAPI used to connect to the SQL database.
    :type driver: str
    :param username: Username to login to SQL database.
    :type username: str
    :param password: Password to login to SQL database.
    :type password: str
    :param database: Name of the database to connect to.
    :type database: str
    :returns: URL string to create SQLAlchemy engine.
    :rtype: str
    """
    engine_string = f"{db_type}+{driver}://{username}:{password}@localhost/{database}"
    return engine_string


# Handles exceptions similar to find_domains pipeline.
# Copied and modeled from MySQLConnectionHandler.execute_transaction,
# but simplified since Engine does the work of connecting to the database.
def execute_transaction(engine, statement_list=[]):
    """Execute list of MySQL statements within a single defined transaction.

    :param engine:
        SQLAlchemy Engine object able to connect to a MySQL databas.
    :type engine: Engine
    :param statement_list:
        a list of any number of MySQL statements with
        no expectation that anything will return
    :returns:
        tuple (result, message)
        WHERE
        result (int) is 0 or 1 status code. 0 means no problems, 1 means problems
        message(str) is a description of the result.
    :rtype: tuple
    """
    msg = "Unable to execute MySQL statements. "
    connection = engine.connect()
    trans = connection.begin()
    try:
        for statement in statement_list:
            connection.execute(statement)
        trans.commit()

    except sqlalchemy.exc.DBAPIError as err:
        err_stmt = err.statement
        sqla_err_type = str(type(err))
        pymysql_err_type = str(type(err.orig))
        pymysql_err_code = err.orig.args[0]
        pymysql_err_msg = err.orig.args[1]
        msg = msg + ("A MySQL error was encountered. "
                     f"SQLAlchemy Error type: {sqla_err_type}. "
                     f"PyMySQL Error type: {pymysql_err_type}. "
                     f"PyMYSQL Error code: {pymysql_err_code}. "
                     f"PyMySQL Error message: {pymysql_err_msg}. "
                     f"Statement: {err_stmt}")
        result = 1
    except TypeError as err:
        msg = msg + "A Python TypeError was encountered. "
        result = 1
    except:
        # TODO not sure how to test this block. Would need to construct a
        # statement that causes a Python built-in exception other than TypeError.
        msg = msg + "An unidentified error was encountered. "
        result = 1
    else:
        msg = "MySQL statements were successfully executed. "
        result = 0

    if result == 1:
        print(msg)
        print("Rolling back transaction...")
        trans.rollback()

    return result, msg

    # Code block below does same thing, but doesn't return a value based if there is an error.
    # with engine.begin() as connection:
    #     for statement in statement_list:
    #         r1 = connection.execute(statement)
