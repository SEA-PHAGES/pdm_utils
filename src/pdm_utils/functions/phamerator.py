"""Functions to interact with PhameratorDB."""


from pdm_utils.classes import genome
from pdm_utils.classes import genomepair
from pdm_utils.classes import cds
from pdm_utils.functions import basic


def parse_phage_table_data(data_dict, trans_table=11, gnm_type=""):
    """Parse a Phamerator database dictionary to create a Genome object.

    :param data_dict:
        Dictionary of data retrieved from the phage table of PhameratorDB.
    :type data_dict: dict
    :param trans_table:
        The translation table that can be used to translate CDS features.
    :type trans_table: int
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
        gnm.host_genus = data_dict["HostStrain"]
    except:
        pass

    try:
        # Sequence data is stored as MEDIUMBLOB, so decode to string.
        gnm.set_sequence(data_dict["Sequence"].decode("utf-8"))
    except:
        pass

    try:
        gnm.length = int(data_dict["SequenceLength"])
    except:
        pass

    try:
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
        gnm.set_cluster_subcluster(data_dict["Cluster"])
    except:
        pass

    try:
        # Singletons are stored in PhameratorDB as NULL, which gets
        # returned as None.
        gnm.set_cluster(data_dict["Cluster2"])
    except:
        pass

    try:
        gnm.set_subcluster(data_dict["Subcluster2"])
    except:
        pass

    try:
        gnm.annotation_status = data_dict["status"]
    except:
        pass

    try:
        gnm.set_retrieve_record(data_dict["RetrieveRecord"])
    except:
        pass

    try:
        gnm.set_annotation_qc(data_dict["AnnotationQC"])
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
    """Parse a Phamerator database dictionary to create a Cds object.

    :param data_dict:
        Dictionary of data retrieved from the gene table of PhameratorDB.
    :type data_dict: dict
    :param trans_table:
        The translation table that can be used to translate CDS features.
    :type trans_table: int
    :returns: A pdm_utils cds object.
    :rtype: cds
    """

    cds_ftr = cds.Cds()
    try:
        cds_ftr.id = data_dict["GeneID"]
    except:
        pass

    try:
        cds_ftr.genome_id = data_dict["PhageID"]
    except:
        pass

    try:
        cds_ftr.left = int(data_dict["Start"])
    except:
        pass

    try:
        cds_ftr.right = int(data_dict["Stop"])
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
        cds_ftr.type = data_dict["TypeID"]
    except:
        pass

    try:
        cds_ftr.set_translation(data_dict["translation"])
    except:
        pass

    try:
        cds_ftr.strand = data_dict["Orientation"]
    except:
        pass

    try:
        cds_ftr.description = data_dict["Notes"].decode("utf-8")
    except:
        pass

    try:
        cds_ftr.locus_tag = data_dict["LocusTag"]
    except:
        pass

    try:
        cds_ftr.translation_table = trans_table
    except:
        pass
    return cds_ftr


def retrieve_data(sql_handle, column=None, query=None, phage_id_list=None):
    """Retrieve genome data from Phamerator for a single genome.

    The query is modified to include one or more PhageIDs

    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
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
                + " WHERE %s IN ('" % column \
                + "','".join(phage_id_list) \
                + "')"
    query = query + ";"
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    return result_list


def parse_cds_data(sql_handle, column=None, phage_id_list=None, query=None):
    """Returns Cds objects containing data parsed from a
    Phamerator database.

    :param sql_handle:
        This parameter is passed directly to the 'retrieve_data' function.
    :type sql_handle: MySQLConnectionHandler
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
                    sql_handle, column=column, query=query,
                    phage_id_list=phage_id_list)
    for data_dict in result_list:
        cds_ftr = parse_gene_table_data(data_dict)
        cds_list.append(cds_ftr)
    return cds_list


def parse_genome_data(sql_handle, phage_id_list=None, phage_query=None,
                      gene_query=None, trna_query=None, gnm_type=""):
    """Returns a list of Genome objects containing data parsed from MySQL
    Phamerator database.

    :param sql_handle:
        This parameter is passed directly to the 'retrieve_data' function.
    :type sql_handle: MySQLConnectionHandler
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
    :returns: A list of pdm_utils Genome objects.
    :rtype: list
    """
    genome_list = []
    result_list1 = retrieve_data(sql_handle, column="PhageID",
                                             phage_id_list=phage_id_list,
                                             query=phage_query)
    for data_dict in result_list1:
        gnm = parse_phage_table_data(data_dict, gnm_type=gnm_type)
        if gene_query is not None:
            cds_list = parse_cds_data(sql_handle, column="PhageID",
                                      phage_id_list=[gnm.id],
                                      query=gene_query)
            x = 0
            while x < len(cds_list):
                cds_list[x].genome_length = gnm.length
                x += 1
            gnm.cds_features = cds_list
        if trna_query is not None:
            # TODO develop this step once tRNA table and objects are built.
            pass
        genome_list.append(gnm)
    return genome_list


# TODO this can be improved if the MCH.execute_query() method
# is able to switch to a standard cursor instead of only using
# dictcursor.
def create_phage_id_set(sql_handle):
    """Create set of phage_ids currently in PhameratorDB.

    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :returns: A set of PhageIDs.
    :rtype: set
    """
    query = "SELECT PhageID FROM phage"
    # Returns a list of items, where each item is a dictionary of
    # SQL data for each row in the table.
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    # Convert to a set of PhageIDs.
    result_set = set([])
    for dict in result_list:
        result_set.add(dict["PhageID"])
    return result_set


def create_seq_set(sql_handle):
    """Create set of genome sequences currently in PhameratorDB.

    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :returns: A set of genome sequences.
    :rtype: set
    """
    query = "SELECT Sequence FROM phage"
    # Returns a list of items, where each item is a dictionary of
    # SQL data for each row in the table.
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    # Convert to a set of sequences.
    # Sequence data is stored as MEDIUMBLOB, so data is returned as bytes
    # "b'AATT", "b'TTCC", etc.
    result_set = set([])
    for dict in result_list:
        result_set.add(dict["Sequence"].decode("utf-8"))
    return result_set


def create_accession_set(sql_handle):
    """Create set of accessions currently in PhameratorDB.

    :param sql_handle:
        A pdm_utils MySQLConnectionHandler object containing
        information on which database to connect to.
    :type sql_handle: MySQLConnectionHandler
    :returns: A set of accessions.
    :rtype: set
    """
    query = "SELECT Accession FROM phage"
    # Returns a list of items, where each item is a dictionary of
    # SQL data for each row in the table.
    result_list = sql_handle.execute_query(query)
    sql_handle.close_connection()
    # Convert to a set of accessions.
    result_set = set([])
    for dict in result_list:
        result_set.add(dict["Accession"])
    return result_set


def create_update(table, field2, value2, field1, value1):
    """Create MySQL UPDATE statement.

    "'UPDATE <table> SET <field2> = '<value2' WHERE <field1> = '<data1>'."

    When the new value to be added is 'singleton' (e.g. for Cluster and
    Cluster2 fields), or an empty value (e.g. None, "none", etc.),
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
    part1 = "UPDATE %s SET %s = " % (table, field2)
    part3 = " WHERE %s = '%s';" % (field1, value1)
    part2a = "NULL"
    part2b = "'%s'" % value2
    if (basic.check_empty(value2) == True or \
        value2.lower() == "singleton"):
        part2 = part2a
    else:
        part2 = part2b
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
    :returns: A MySQL UPDATE statement.
    :rtype: str
    """
    statement = "DELETE FROM %s WHERE %s = '%s';" % (table, field, data)
    return statement


# TODO should this function first check datatypes and formats?
# Strand should be either "F" or "R"
# Start, stop, length = int
# TypeID = 'CDS'
def create_gene_table_insert(cds_ftr):
    """Create a MySQL gene table INSERT statement.

    :param cds_ftr: A pdm_utils Cds object.
    :type cds_ftr: Cds
    :returns:
        A MySQL statement to INSERT a new row in the 'gene' table
        with data for several fields.
    :rtype: str
    """
    statement = ("INSERT INTO gene "
                 "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, "
                 "translation, Orientation, Notes, LocusTag) "
                 "VALUES "
                 "('%s', '%s', %s, %s, %s, '%s', '%s', '%s', '%s', '%s', '%s');"
                 % (cds_ftr.id,
                    cds_ftr.genome_id,
                    cds_ftr.left,
                    cds_ftr.right,
                    cds_ftr.translation_length,
                    cds_ftr.name,
                    cds_ftr.type,
                    cds_ftr.translation,
                    cds_ftr.strand,
                    cds_ftr.processed_description,
                    cds_ftr.locus_tag)
                 )
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
    statement = ("INSERT INTO phage (PhageID, Accession, Name, "
                 "HostStrain, Sequence, SequenceLength, GC, status, "
                 "DateLastModified, RetrieveRecord, AnnotationQC, "
                 "AnnotationAuthor) VALUES ("
                 "'%s', '%s', '%s', '%s', '%s',"
                 " %s, %s, '%s', '%s', %s, %s, %s);"
                 % (gnm.id,
                 gnm.accession,
                 gnm.name,
                 gnm.host_genus,
                 gnm.seq,
                 gnm.length,
                 gnm.gc,
                 gnm.annotation_status,
                 gnm.date,
                 gnm.retrieve_record,
                 gnm.annotation_qc,
                 gnm.annotation_author)
                 )
    return statement


# TODO integration test? Not sure if this is really needed.
def create_genome_insert(gnm):
    """Create a collection of MySQL INSERT and UPDATE statements.

    :param gnm: A pdm_utils Genome object.
    :type gnm: Genome
    :returns:
        A list of MySQL INSERT and UPDATE statement for the
        phage, gene, and trna tables.
    :rtype: list
    """
    table = "phage"
    field1 = "PhageID"
    value1 = gnm.id
    statements = []
    statements.append(create_phage_table_insert(gnm))
    statements.append(create_update(
        table, field1, value1, "Cluster", gnm.cluster_subcluster))
    statements.append(create_update(
        table, field1, value1, "Cluster2", gnm.cluster))
    statements.append(create_update(
        table, field1, value1, "Subcluster2", gnm.subcluster))
    return statements


def copy_data(bndl, from_type, to_type, flag="retain"):
    """Copy data from a 'phamerator' genome object.

    If a genome object stored in the Bundle object has
    attributes that are set to be 'retained' from Phamerator,
    copy any necessary data from the genome with 'type' attribute
    set to 'phamerator' to the new genome.

    :param bndl: A pdm_utils Bundle object.
    :type bndl: Bundle
    :param from_type:
        Indicates the value of the source genome's 'type',
        indicating the genome from which data will be copied.
    :type from_type: str
    :param to_type:
        Indicates the value of the target genome's 'type',
        indicating the genome to which data will be copied.
    :type to_type: str
    :param flag:
        Indicates the value that attributes of the target genome object
        must have in order be updated from the 'phamerator' genome object.
    :type flag: str
    """
    if to_type in bndl.genome_dict.keys():
        to_gnm = bndl.genome_dict[to_type]
        to_gnm.set_value_flag(flag)
        if to_gnm._value_flag:
            if from_type in bndl.genome_dict.keys():
                from_gnm = bndl.genome_dict[from_type]

                # Copy all data that is set to be copied and
                # add to Bundle object.
                genome_pair = genomepair.GenomePair()
                genome_pair.genome1 = to_gnm
                genome_pair.genome2 = from_gnm
                genome_pair.copy_data("type", from_gnm.type, to_gnm.type, flag)
                bndl.set_genome_pair(genome_pair, to_gnm.type, from_gnm.type)
        to_gnm.set_value_flag(flag)
























# TODO unit test below.

















# TODO need to work on this.
# TODO implement.
# TODO unit test.
def implement_update_statements():

    #If it looks like there is a problem with some of the genomes on the list,
    #cancel the transaction, otherwise proceed
    if updated == update_total:
        if run_type == "production":
            con = mdb.connect(mysqlhost, username, password, database)
            con.autocommit(False)
            cur = con.cursor()

            try:
                cur.execute("START TRANSACTION")
                for statement in update_statements:
                    cur.execute(statement)
                    write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
                cur.execute("COMMIT")
                write_out(output_file,"\nAll update statements committed.")
                cur.close()
                con.autocommit(True)

            except:
                success_action_file_handle.close()
                mdb_exit("\nError: problem updating genome information.\nNo changes have been made to the database.")

            con.close()
        else:
            write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)
    else:
        write_out(output_file,"\nError: problem processing data list to update genomes. Check input table format.\nNo changes have been made to the database.")
        write_out(output_file,"\nExiting import script.")
        output_file.close()
        success_action_file_handle.close()
        sys.exit(1)

    #Document the update actions
    for element in update_data_list:
        if element[7] == "":
            element[7] = "none"

        update_output_list = [element[0],\
                                element[1],\
                                element[2],\
                                element[3],\
                                element[8],\
                                element[4],\
                                author_dictionary[element[9]],\
                                element[5],\
                                element[7],\
                                element[10],\
                                element[6]]
        success_action_file_writer.writerow(update_output_list)

    write_out(output_file,"\nAll field update actions have been implemented.")
    raw_input("\nPress ENTER to proceed to next import stage.")




# TODO need to work on this.
# TODO implement.
# TODO unit test.
def implement_remove_statements():

    #If it looks like there is a problem with some of the genomes on the list,
    #cancel the transaction, otherwise proceed

    if removed == remove_total:

        if run_type == "production":
            con = mdb.connect(mysqlhost, username, password, database)
            con.autocommit(False)
            cur = con.cursor()

            try:
                cur.execute("START TRANSACTION")
                for statement in removal_statements:
                    cur.execute(statement)
                    write_out(output_file,"\n" + statement + " executed successfully, but not yet committed.")
                cur.execute("COMMIT")
                write_out(output_file,"\nAll remove statements committed.")
                cur.close()
                con.autocommit(True)

            except:
                success_action_file_handle.close()
                mdb_exit("\nError: problem removing genomes with no replacements.\nNo remove actions have been implemented.")
            con.close()
        else:
            write_out(output_file,"\nRUN TYPE IS %s, SO NO CHANGES TO THE DATABASE HAVE BEEN IMPLEMENTED.\n" % run_type)
    else:
        write_out(output_file,"\nError: problem processing data list to remove genomes. Check input table format.\nNo remove actions have been implemented.")
        output_file.close()
        success_action_file_handle.close()
        sys.exit(1)

    #Document the remove actions
    for element in remove_data_list:
        remove_output_list = [element[0],\
                                element[1],\
                                element[2],\
                                element[3],\
                                element[8],\
                                element[4],\
                                element[9],\
                                element[5],\
                                element[7],\
                                element[10],\
                                element[6]]
        success_action_file_writer.writerow(remove_output_list)
    write_out(output_file,"\nAll genome remove actions have been implemented.")
    raw_input("\nPress ENTER to proceed to next import stage.")








# # TODO probably no longer needed.
# # TODO this needs to be revamped. It was originally constructed to handle
# # 'update' tickets that contain several pieces of data.
# def create_genome_update_statements(gnm):
#     """Create a collection of genome-level UPDATE statements using data
#     in a Genome object.
#
#     :param gnm: A pdm_utils Genome object.
#     :type gnm: Genome
#     :returns:
#         A list of MySQL INSERT and UPDATE statement for the
#         phage, gene, and trna tables.
#     :rtype: list
#
#     """
#
#     table = "phage"
#     field1 = "PhageID"
#     value1 = gnm.id
#     statements = []
#     statements.append(create_update(
#         table, field1, value1, "HostStrain", gnm.host_genus))
#     statements.append(create_update(
#         table, field1, value1, "status", gnm.annotation_status))
#     statements.append(create_update(
#         table, field1, value1, "Accession", gnm.accession))
#     statements.append(create_update(
#         table, field1, value1, "AnnotationAuthor", gnm.author))
#     statements.append(create_update(
#         table, field1, value1, "Cluster", gnm.cluster_subcluster))
#     statements.append(create_update(
#         table, field1, value1, "Cluster2", gnm.cluster))
#     statements.append(create_update(
#         table, field1, value1, "Subcluster2", gnm.subcluster))
#     return statements






###
