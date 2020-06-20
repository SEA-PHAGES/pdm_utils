"""Functions to interact with MySQL."""

import getpass
import subprocess
import sys
import re

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import sqlalchemy

from pdm_utils.classes import cds, trna, tmrna
from pdm_utils.classes import genome
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import mysqldb_basic

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

    ftr = cds.Cds()
    ftr.type = "CDS"
    try:
        ftr.id = data_dict["GeneID"]
    except:
        pass

    try:
        ftr.genome_id = data_dict["PhageID"]
    except:
        pass

    try:
        ftr.start = int(data_dict["Start"])
    except:
        pass

    try:
        ftr.stop = int(data_dict["Stop"])
    except:
        pass

    try:
        ftr.parts = int(data_dict["Parts"])
    except:
        pass

    ftr.coordinate_format = "0_half_open"

    try:
        ftr.length = int(data_dict["Length"])
    except:
        pass

    try:
        ftr.name = data_dict["Name"]
    except:
        pass

    try:
        ftr.set_translation(data_dict["Translation"].decode("utf-8"))
    except:
        pass

    try:
        ftr.orientation = data_dict["Orientation"]
    except:
        pass

    try:
        ftr.description = data_dict["Notes"].decode("utf-8")
    except:
        pass

    try:
        ftr.set_locus_tag(data_dict["LocusTag"])
    except:
        pass

    try:
        ftr.pham_id = int(data_dict["PhamID"])
    except:
        pass

    try:
        ftr.domain_status = int(data_dict["DomainStatus"])
    except:
        pass

    try:
        ftr.translation_table = trans_table
    except:
        pass

    return ftr


# TODO Christian review
def parse_trna_table_data(data_dict):
    """Parse a MySQL database dictionary to create a Trna object.

    :param data_dict:
        Dictionary of data retrieved from the gene table.
    :type data_dict: dict
    :returns: A pdm_utils Trna object.
    :rtype: Trna
    """

    ftr = trna.Trna()
    ftr.type = "tRNA"

    try:
        ftr.id = data_dict["GeneID"]
    except:
        pass

    try:
        ftr.genome_id = data_dict["PhageID"]
    except:
        pass

    try:
        ftr.start = int(data_dict["Start"])
    except:
        pass

    try:
        ftr.stop = int(data_dict["Stop"])
    except:
        pass

    try:
        ftr.length = int(data_dict["Length"])
    except:
        pass

    try:
        ftr.name = data_dict["Name"]
    except:
        pass

    try:
        ftr.orientation = data_dict["Orientation"]
    except:
        pass

    try:
        ftr.set_locus_tag(data_dict["LocusTag"])
    except:
        pass

    try:
        ftr.note = data_dict["Note"].decode("utf-8")
    except:
        pass

    try:
        ftr.amino_acid = data_dict["AminoAcid"]
    except:
        pass

    try:
        ftr.anticodon = data_dict["Anticodon"]
    except:
        pass

    try:
        ftr.structure = data_dict["Structure"].decode("utf-8")
    except:
        pass

    try:
        ftr.use = data_dict["Source"]
    except:
        pass

    try:
        ftr.parts = 1
    except:
        pass

    ftr.coordinate_format = "0_half_open"

    return ftr

# TODO Christian review
def parse_tmrna_table_data(data_dict):
    """Parse a MySQL database dictionary to create a Tmrna object.

    :param data_dict:
        Dictionary of data retrieved from the gene table.
    :type data_dict: dict
    :returns: A pdm_utils Tmrna object.
    :rtype: Tmrna
    """

    ftr = tmrna.Tmrna()
    ftr.type = "tmRNA"

    try:
        ftr.id = data_dict["GeneID"]
    except:
        pass

    try:
        ftr.genome_id = data_dict["PhageID"]
    except:
        pass

    try:
        ftr.start = int(data_dict["Start"])
    except:
        pass

    try:
        ftr.stop = int(data_dict["Stop"])
    except:
        pass

    try:
        ftr.length = int(data_dict["Length"])
    except:
        pass

    try:
        ftr.name = data_dict["Name"]
    except:
        pass

    try:
        ftr.orientation = data_dict["Orientation"]
    except:
        pass

    try:
        ftr.set_locus_tag(data_dict["LocusTag"])
    except:
        pass

    try:
        ftr.note = data_dict["Note"].decode("utf-8")
    except:
        pass

    try:
        ftr.peptide_tag = data_dict["PeptideTag"]
    except:
        pass

    try:
        ftr.parts = 1
    except:
        pass

    ftr.coordinate_format = "0_half_open"

    return ftr


def parse_feature_data(engine, ftr_type, column=None, phage_id_list=None,
                       query=None):
    """Returns Cds objects containing data parsed from a
    MySQL database.

    :param engine:
        This parameter is passed directly to the 'retrieve_data' function.
    :type engine: Engine
    :param query:
        This parameter is passed directly to the 'retrieve_data' function.
    :type query: str
    :param ftr_type:
        Indicates the type of features retrieved.
    :type ftr_type: str
    :param column:
        This parameter is passed directly to the 'retrieve_data' function.
    :type column: str
    :param phage_id_list:
        This parameter is passed directly to the 'retrieve_data' function.
    :type phage_id_list: list
    :returns: A list of pdm_utils Cds objects.
    :rtype: list
    """
    result_list = mysqldb_basic.retrieve_data(engine, column=column,
                                              query=query,
                                              id_list=phage_id_list)
    ftrs = []
    for data_dict in result_list:
        if ftr_type == "cds":
            ftr = parse_gene_table_data(data_dict)
        elif ftr_type == "trna":
            ftr = parse_trna_table_data(data_dict)
        elif ftr_type == "tmrna":
            ftr = parse_tmrna_table_data(data_dict)
        else:
            # If the ftr_type is invalid, just take the data dictionary.
            # Alternatively it could raise an error.
            ftr = data_dict
        ftrs.append(ftr)
    return ftrs


def parse_genome_data(engine, phage_id_list=None, phage_query=None,
                      gene_query=None, trna_query=None, tmrna_query=None,
                      gnm_type=""):
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
        This parameter is passed directly to the 'parse_feature_data' function
        to retrieve data from the gene table.
        If not None, pdm_utils Cds objects for all of the phage's
        CDS features in the gene table will be constructed
        and added to the Genome object.
    :type gene_query: str
    :param trna_query:
        This parameter is passed directly to the 'parse_feature_data' function
        to retrieve data from the trna table.
        If not None, pdm_utils Trna objects for all of the phage's
        tRNA features in the trna table will be constructed
        and added to the Genome object.
    :type trna_query: str
    :param tmrna_query:
        This parameter is passed directly to the 'parse_feature_data' function
        to retrieve data from the tmrna table.
        If not None, pdm_utils Tmrna objects for all of the phage's
        tmRNA features in the tmrna table will be constructed
        and added to the Genome object.
    :type tmrna_query: str
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
    COLUMN = "PhageID"
    genome_list = []
    result_list1 = mysqldb_basic.retrieve_data(engine, column=COLUMN,
                                               id_list=phage_id_list,
                                               query=phage_query)
    for data_dict in result_list1:
        gnm = parse_phage_table_data(data_dict, gnm_type=gnm_type)

        if gene_query is not None:
            cds_list = parse_feature_data(engine, "cds", column=COLUMN,
                                          phage_id_list=[gnm.id],
                                          query=gene_query)

            for x in range(len(cds_list)):
                cds_list[x].genome_length = gnm.length
            gnm.cds_features = cds_list

        if trna_query is not None:
            trna_list = parse_feature_data(engine, "trna", column=COLUMN,
                                          phage_id_list=[gnm.id],
                                          query=trna_query)

            for x in range(len(trna_list)):
                trna_list[x].genome_length = gnm.length
            gnm.trna_features = trna_list

        if tmrna_query is not None:
            tmrna_list = parse_feature_data(engine, "tmrna", column=COLUMN,
                                          phage_id_list=[gnm.id],
                                          query=tmrna_query)

            for x in range(len(tmrna_list)):
                tmrna_list[x].genome_length = gnm.length
            gnm.tmrna_features = tmrna_list

        genome_list.append(gnm)
    return genome_list


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
    part2 = mysqldb_basic.convert_for_sql(value2, check_set=check_set, single=True)
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
    locus_tag = mysqldb_basic.convert_for_sql(cds_ftr.locus_tag,
                                              check_set={""}, single=False)

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
                                 cds_ftr.length,
                                 cds_ftr.name, cds_ftr.translation,
                                 cds_ftr.orientation, cds_ftr.description,
                                 locus_tag, cds_ftr.parts)
    return statement


def create_trna_table_insert(trna_ftr):
    """
    Create a MySQL trna table INSERT statement.
    :param trna_ftr: a pdm_utils Trna object
    :type trna_ftr: Trna
    :returns: a MySQL statement to INSERT a new row in the 'trna' table
    with all of trna_ftr's relevant data
    :rtype: str
    """
    # Any values that could be none should be run through mysqldb_basic's
    # convert_for_sql method so that the statement is does not cause errors
    # when evaluated
    geneid, phageid = trna_ftr.id, trna_ftr.genome_id
    name = trna_ftr.name
    locus_tag = mysqldb_basic.convert_for_sql(
        trna_ftr.locus_tag, check_set={""}, single=False)
    start, stop, length = trna_ftr.start, trna_ftr.stop, trna_ftr.length
    orientation = trna_ftr.orientation
    note = mysqldb_basic.convert_for_sql(
        trna_ftr.note, check_set={""}, single=False)
    amino_acid, anticodon = trna_ftr.amino_acid, trna_ftr.anticodon
    structure = mysqldb_basic.convert_for_sql(
        trna_ftr.structure, check_set={""}, single=False)
    source = mysqldb_basic.convert_for_sql(
        trna_ftr.use, check_set={None}, single=False)

    statement = ("""INSERT INTO trna """
                 """(GeneID, PhageID, Start, Stop, Length, """
                 """Name, Orientation, Note, LocusTag, AminoAcid, Anticodon, """
                 """Structure, Source) VALUES ("{}", "{}", {}, {}, {}, """
                 """"{}", "{}", {}, {}, "{}", "{}", {}, {});""")
    statement = statement.format(geneid, phageid, start, stop, length,
                                 name, orientation, note, locus_tag,
                                 amino_acid, anticodon, structure, source)
    return statement


def create_tmrna_table_insert(tmrna_ftr):
    """

    :param tmrna_ftr:
    :return:
    """
    geneid, phageid = tmrna_ftr.id, tmrna_ftr.genome_id
    name = tmrna_ftr.name
    locus_tag = mysqldb_basic.convert_for_sql(
        tmrna_ftr.locus_tag, check_set={""}, single=False)
    start, stop, length = tmrna_ftr.start, tmrna_ftr.stop, tmrna_ftr.length
    orientation = tmrna_ftr.orientation
    note = mysqldb_basic.convert_for_sql(
        tmrna_ftr.note, check_set={""}, single=False)
    peptide_tag = mysqldb_basic.convert_for_sql(
        tmrna_ftr.peptide_tag, check_set={""}, single=False)

    statement = ("""INSERT INTO tmrna """
                 """(GeneID, PhageID, Start, Stop, Length, """
                 """Name, Orientation, Note, LocusTag, PeptideTag) VALUES """
                 """("{}", "{}", {}, {}, {}, "{}", "{}", {}, {}, {});""")
    statement = statement.format(geneid, phageid, start, stop, length,
                                 name, orientation, note, locus_tag,
                                 peptide_tag)
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
    cluster = mysqldb_basic.convert_for_sql(gnm.cluster,
                                            check_set={"Singleton"},
                                            single=True)
    subcluster = mysqldb_basic.convert_for_sql(gnm.subcluster,
                                               check_set={"none"},
                                               single=True)

    # gnm.seq is a BioPython Seq object.
    # It is coerced to string by default.
    statement = ("""INSERT INTO phage """
                 """(PhageID, Accession, Name, HostGenus, Sequence, """
                 """Length, GC, Status, DateLastModified, RetrieveRecord, """
                 """AnnotationAuthor, Cluster, Subcluster) """
                 """VALUES """
                 """('{}', '{}', '{}', '{}', '{}', {}, {}, '{}', """
                 """'{}', {}, {}, {}, {});""")
    statement = statement.format(
                    gnm.id, gnm.accession, gnm.name, gnm.host_genus, gnm.seq,
                    gnm.length, gnm.gc, gnm.annotation_status, gnm.date,
                    gnm.retrieve_record, gnm.annotation_author, cluster,
                    subcluster)

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

    stmts = []
    if tkt_type == "replace":
        stmt1 = create_delete("phage", "PhageID", gnm.id)
        stmts.append(stmt1)
    stmt2 = create_phage_table_insert(gnm)
    stmts.append(stmt2)
    for cds_ftr in gnm.cds_features:
        stmt3 = create_gene_table_insert(cds_ftr)
        stmts.append(stmt3)
    for trna_ftr in gnm.trna_features:
        stmt4 = create_trna_table_insert(trna_ftr)
        stmts.append(stmt4)
    for tmrna_ftr in gnm.tmrna_features:
        stmt5 = create_tmrna_table_insert(tmrna_ftr)
        stmts.append(stmt5)

    return stmts



def change_version(engine, amount=1):
    """Change the database version number.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :param amount: Amount to increment/decrement version number.
    :type amount: int
    """
    result = mysqldb_basic.get_first_row_data(engine, "version")
    current = result["Version"]
    new = current + amount
    print(f"Updating version from {current} to {new}.")
    statement = (f"UPDATE version SET Version = {new}")
    engine.execute(statement)



# TODO unittest.
def get_schema_version(engine):
    """Identify the schema version of the database_versions_list.

    Schema version data has not been persisted in every schema version,
    so if schema version data is not found, it is deduced from other
    parts of the schema.

    :param engine: SQLAlchemy Engine object able to connect to a MySQL database.
    :type engine: Engine
    :returns: The version of the pdm_utils database schema.
    :rtype: int
    """
    # 1. If the version table does not exist, schema_version = 0.
    # 2. If there is no schema_version or SchemaVersion field,
    #    it is either schema_version = 1 or 2.
    # 3. If AnnotationAuthor, Program, AnnotationQC, and RetrieveRecord
    #    columns are in phage table, schema_version = 2.


    db_tables = mysqldb_basic.get_tables(engine, engine.url.database)
    if "version" in db_tables:
        version_table = True
    else:
        version_table = False

    if version_table == True:
        version_columns = mysqldb_basic.get_first_row_data(engine, "version")
        if "schema_version" in version_columns.keys():
            schema_version = version_columns["schema_version"]
        elif "SchemaVersion" in version_columns.keys():
            schema_version = version_columns["SchemaVersion"]
        else:
            phage_columns = mysqldb_basic.get_columns(
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

    # Code block below does same thing, but doesn't return a value
    # based on if there is an error.
    # with engine.begin() as connection:
    #     for statement in statement_list:
    #         r1 = connection.execute(statement)

