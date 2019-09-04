"""Functions to interact with PhameratorDB."""


from classes import Genome
from classes import GenomePair
from classes import cds
from functions import basic
import pymysql


def parse_phage_table_data(data_dict, trans_table=11):
    """Parse a Phamerator database dictionary to create a Genome object."""

    genome = Genome.Genome()
    try:
        genome.id = data_dict["PhageID"]
    except:
        pass

    try:
        genome.accession = data_dict["Accession"]
    except:
        pass

    try:
        genome.name = data_dict["Name"]
    except:
        pass

    try:
        genome.host_genus = data_dict["HostStrain"]
    except:
        pass

    try:
        # Sequence data is stored as MEDIUMBLOB, so decode to string.
        genome.set_sequence(data_dict["Sequence"].decode("utf-8"))
    except:
        pass

    try:
        genome._length = data_dict["SequenceLength"]
    except:
        pass

    try:
        genome.date = data_dict["DateLastModified"]
    except:
        pass

    try:
        genome.description = data_dict["Notes"].decode("utf-8")
    except:
        pass

    try:
        genome._gc = data_dict["GC"]
    except:
        pass

    try:
        genome.set_cluster_subcluster(data_dict["Cluster"])
    except:
        pass

    try:
        # Singletons are stored in PhameratorDB as NULL, which gets
        # returned as None.
        genome.set_cluster(data_dict["Cluster2"])
    except:
        pass

    try:
        genome.set_subcluster(data_dict["Subcluster2"])
    except:
        pass

    try:
        genome.annotation_status = data_dict["status"]
    except:
        pass

    try:
        genome.retrieve_record = data_dict["RetrieveRecord"]
    except:
        pass

    try:
        genome.annotation_qc = data_dict["AnnotationQC"]
    except:
        pass

    try:
        genome.annotation_author = data_dict["AnnotationAuthor"]
    except:
        pass

    genome.translation_table = trans_table
    genome.type = "phamerator"
    return genome


def parse_gene_table_data(data_dict, trans_table=11):
    """Parse a Phamerator database dictionary to create a Cds object."""

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
        cds_ftr.left = data_dict["Start"]
    except:
        pass

    try:
        cds_ftr.right = data_dict["Stop"]
    except:
        pass

    try:
        cds_ftr._length = data_dict["Length"]
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


# TODO revamp so that sql handler makes the query instead of
# simply passing data to pymysql.
def retrieve_data(sql_handle, column=None, query=None, phage_id_list=None):
    """Retrieve genome data from Phamerator for a single genome.

    The function expects a query that selects valid, specific columns
    from the Phage table but does not condition on a PhageID.
    (e.g. 'SELECT PhageID,Cluster FROM phage')
    """

    if (phage_id_list is not None and len(phage_id_list) > 0):
        query = query \
                + " WHERE %s IN ('" % column \
                + "','".join(phage_id_list) \
                + "')"
    query = query + ";"
    # Create the connection.
    connection = pymysql.connect(host = "localhost",
                                    user = sql_handle.username,
                                    password = sql_handle.password,
                                    database = sql_handle.database,
                                    cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(query)

    # Data is returned as a list of items, where each item is a
    # dictionary of SQL data for each PhageID.
    result_list = cur.fetchall()
    connection.close()
    return result_list


def parse_cds_data(sql_handle, column=None, phage_id_list=None, query=None):
    """Returns Cds objects containing data parsed from a
    Phamerator database."""

    cds_list = []
    result_list = retrieve_data(
                    sql_handle, column=column, query=query,
                    phage_id_list=phage_id_list)
    for data_dict in result_list:
        cds_ftr = parse_gene_table_data(result_list[0])
        cds_list.append(cds_ftr)
    return cds_list


def parse_genome_data(sql_handle, phage_id_list=None, phage_query=None,
                      gene_query=None, trna_query=None):
    """Returns a list of Genome objects containing data parsed from MySQL
    Phamerator database.

    If the 'phage_id' parameter contains a list of at least
    one valid PhageID, a Genome object will be constructed
    only for that phage. If the 'phage_id' parameter is None,
    Genome objects for all phages in the database will be constructed.
    If the 'gene_query' parameter is not None, Cds objects for all of
    the phage's CDS features in the gene table will be constructed
    and added to the Genome object using the provided query.
    If the 'trna_query' parameter is True, Trna objects for all of
    the phage's tRNA features in the tRNA table will be constructed
    add added to the Genome object using the provided query.
    """
    genome_list = []
    result_list1 = retrieve_data(sql_handle, column="PhageID",
                                             phage_id_list=phage_id_list,
                                             query=phage_query)
    for data_dict in result_list1:
        genome = parse_phage_table_data(data_dict)
        if gene_query is not None:
            cds_list = parse_cds_data(sql_handle, column="PhageID",
                                      phage_id_list=[genome.id],
                                      query=gene_query)
            genome.cds_features = cds_list
        if trna_query is not None:
            # TODO develop this step once tRNA table and objects are built.
            pass
        genome_list.append(genome)
    return genome_list


# TODO revamp so that sql handler makes the query instead of
# simply passing data to pymysql.
def create_phage_id_set(sql_handle):
    """Create set of phage_ids currently in PhameratorDB."""

    # Create the connection.
    connection = pymysql.connect(host = "localhost",
                                    user = sql_handle.username,
                                    password = sql_handle.password,
                                    database = sql_handle.database)
    cur = connection.cursor()
    cur.execute("SELECT PhageID FROM phage")

    # Data is returned as a tuple of tuples: (("L5",), ("Trixie",), etc.)
    result_tuple = cur.fetchall()
    connection.close()

    # Convert to a set of PhageIDs.
    result_set = set([])
    for tup in result_tuple:
        result_set.add(tup[0])
    return result_set


# TODO revamp so that sql handler makes the query instead of
# simply passing data to pymysql.
def create_seq_set(sql_handle):
    """Create set of genome sequences currently in PhameratorDB."""


    # Create the connection.
    connection = pymysql.connect(host = "localhost",
                                    user = sql_handle.username,
                                    password = sql_handle.password,
                                    database = sql_handle.database)
    cur = connection.cursor()
    cur.execute("SELECT Sequence FROM phage")

    # Data is returned as a tuple of tuples
    result_tuple = cur.fetchall()
    connection.close()

    # Convert to a set of sequences.
    # Sequence data is stored as MEDIUMBLOB, so data is returned as bytes
    # (("b'AATT",), ("b'TTCC",), etc.)
    result_set = set([])
    for tup in result_tuple:
        result_set.add(tup[0].decode("utf-8"))
    return result_set


def create_update_statement(table, field1, value1, field2, value2):
    """Create MySQL UPDATE statement. When:
    The new value to be added is 'singleton' (e.g. for Cluster and
    Cluster2 fields), or
    The new value to be added is an empty value (e.g. None, "none", etc.),
    the new value is set to NULL."""


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



def create_genome_update_statements(genome):
    """Create a collection of genome-level UPDATE statements using data
    in a Genome object."""

    table = "phage"
    field1 = "PhageID"
    value1 = genome.id
    statements = []
    statements.append(create_update_statement( \
        table, field1, value1, "HostStrain", genome.host_genus))
    statements.append(create_update_statement( \
        table, field1, value1, "status", genome.annotation_status))
    statements.append(create_update_statement( \
        table, field1, value1, "Accession", genome.accession))
    statements.append(create_update_statement( \
        table, field1, value1, "AnnotationAuthor", genome.author))
    statements.append(create_update_statement( \
        table, field1, value1, "Cluster", genome.cluster_subcluster))
    statements.append(create_update_statement( \
        table, field1, value1, "Cluster2", genome.cluster))
    statements.append(create_update_statement( \
        table, field1, value1, "Subcluster2", genome.subcluster))
    return statements


def create_delete_statement(table, field1, data1):
    """Create MySQL DELETE statement."""
    statement = "DELETE FROM %s WHERE %s = '%s';" % (table, field1, data1)
    return statement


# TODO this may no longer be needed.
def create_genome_delete_statement(genome):
    """Create a genome-level DELETE statements using data
    in a Genome object."""

    table = "phage"
    field1 = "PhageID"
    value1 = genome.id
    statements = []
    statements.append(create_delete_statement(table, field1, value1))
    return statements


def create_cds_insert_statement(cds_feature):
    """Create a CDS-level INSERT statement using data in a CDS object."""

    statement = "INSERT INTO gene " + \
        "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
        "translation, Orientation, Notes, LocusTag) " + \
        "VALUES " + \
        "('%s', '%s', %s, %s, %s, '%s', '%s', '%s', '%s', '%s', '%s');" % \
        (cds_feature.id, \
        cds_feature.genome_id, \
        cds_feature.left, \
        cds_feature.right, \
        cds_feature._translation_length, \
        cds_feature.name, \
        cds_feature.type, \
        cds_feature.translation, \
        cds_feature.strand, \
        cds_feature.processed_description, \
        cds_feature.locus_tag)
    return statement


# TODO this function could also receive a genome object.
def create_cds_insert_statements(list_of_features):
    """Create a collection of CDS-level INSERT statements using data
    in a list of CDS objects."""

    statements = []
    for cds_feature in list_of_features:
        statements.append(create_cds_insert_statement(cds_feature))
    return statements


def create_genome_insert_statement(genome):
    """Create a genome-level INSERT statements using data
    in a Genome object."""

    statement = \
        "INSERT INTO phage (PhageID, Accession, Name, " + \
        "HostStrain, Sequence, SequenceLength, GC, status, " + \
        "DateLastModified, RetrieveRecord, AnnotationQC, " + \
        "AnnotationAuthor) " + \
        "VALUES (" + \
        "'%s', '%s', '%s', '%s', '%s', %s, %s, '%s', '%s', '%s', '%s', '%s');" \
        % (genome.id, \
        genome.accession, \
        genome.name, \
        genome.host_genus, \
        genome.seq, \
        genome._length, \
        genome._gc, \
        genome.annotation_status, \
        genome.date, \
        genome.retrieve_record, \
        genome.annotation_qc, \
        genome.annotation_author)
    return statement


def create_genome_insert_statements(genome):
    """Create a collection of genome-level INSERT statements using data
    in a Genome object."""

    table = "phage"
    field1 = "PhageID"
    value1 = genome.id
    statements = []
    statements.append(create_genome_insert_statement(genome))
    statements.append(create_update_statement( \
        table, field1, value1, "Cluster", genome.cluster_subcluster))
    statements.append(create_update_statement( \
        table, field1, value1, "Cluster2", genome.cluster))
    statements.append(create_update_statement( \
        table, field1, value1, "Subcluster2", genome.subcluster))
    return statements


def copy_data_from(bndl, type, flag="retain"):
    """Copy data from a 'phamerator' genome object.

    If a genome object stored in the Bundle object has
    attributes that are set to be 'retained' from Phamerator,
    copy any necessary data from the Phamerator genome to the new genome.
    The 'type' parameter indicates the type of genome that may need
    to be populated from Phamerator."""

    if type in bndl.genome_dict.keys():

        genome1 = bndl.genome_dict[type]
        genome1.set_value_flag(flag)

        if genome1._value_flag:

            if "phamerator" in bndl.genome_dict.keys():

                genome2 = bndl.genome_dict["phamerator"]

                # Copy all data that is set to be copied and
                # add to Bundle object.
                genome_pair = GenomePair.GenomePair()
                genome_pair.genome1 = genome1
                genome_pair.genome2 = genome2
                genome_pair.copy_data("type", genome2.type, genome1.type, flag)
                bndl.set_genome_pair(genome_pair, genome1.type, genome2.type)

        # Now record an error if there are still fields
        # that need to be retained.
        genome1.set_value_flag(flag)
        genome1.check_value_flag()




























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















#TODO below: functions that may no longer be needed.







# # TODO this may no longer be needed now that
# # parse_phage_table_data() is available.
# def parse_phamerator_data(genome, data_tuple):
#     """Parses tuple of data derived from a Phamerator database
#     and populates a genome object.
#     Expected data structure:
#     0 = PhageID
#     1 = Name
#     2 = HostStrain
#     3 = Sequence
#     4 = status
#     5 = Cluster2
#     6 = DateLastModified
#     7 = Accession
#     8 = Subcluster2
#     9 = AnnotationAuthor
#     10 = AnnotationQC
#     11 = RetrieveRecord
#     """
#
#     genome.set_id(data_tuple[0])
#     genome.name = data_tuple[1]
#     genome.set_host_genus(data_tuple[2])
#     genome.set_sequence(data_tuple[3])
#     genome.annotation_status = data_tuple[4]
#     genome.set_cluster(data_tuple[5])
#     genome.set_subcluster(data_tuple[8])
#     genome.set_date(data_tuple[6])
#     genome.set_accession(data_tuple[7])
#     genome.annotation_author = str(data_tuple[9])
#     genome.annotation_qc = str(data_tuple[10])
#     genome.retrieve_record = str(data_tuple[11])
#     genome.type = "phamerator"





# TODO this may no longer be needed now that
# parse_phage_table_data() is available.
# def create_phamerator_dict(phamerator_data_tuples):
#     """
#     Returns a dictionary of Phamerator data retrieved from MySQL query.
#     Key = PhageID.
#     Value = Genome object containing parsed MySQL data.
#     """
#
#     genome_dict = {}
#     for genome_tuple in phamerator_data_tuples:
#         genome = Genome.Genome()
#         parse_phamerator_data(genome,genome_tuple)
#         genome_dict[genome.id] = genome
#
#     return genome_dict

#
# # TODO this may no longer be needed.
# def create_data_sets(genome_dict):
#     """
#     Create sets of all unique values for several fields in the Phamerator data.
#     """
#     phage_id_set = set()
#     host_genus_set = set()
#     status_set = set()
#     cluster_set = set()
#     accession_set = set()
#     subcluster_set = set()
#
#     for genome_id in genome_dict.keys():
#
#         genome = genome_dict[genome_id]
#         phage_id_set.add(genome.id)
#         host_genus_set.add(genome.host_genus)
#         status_set.add(genome.annotation_status)
#         cluster_set.add(genome.cluster)
#
#         # TODO this was not implemented in original import script,
#         # so maybe the subcluster 'empty check' is not needed.
#         # Only add to the accession set if there was an accession,
#         # and not if it was empty.
#         if basic.check_empty(genome.subcluster) == False:
#             subcluster_set.add(genome.subcluster)
#
#         # Only add to the accession set if there was an accession,
#         # and not if it was empty.
#         if basic.check_empty(genome.accession) == False:
#             accession_set.add(genome.accession)
#
#     dictionary_of_sets = {}
#     dictionary_of_sets["phage_id"] = phage_id_set
#     dictionary_of_sets["host_genus"] = host_genus_set
#     dictionary_of_sets["annotation_status"] = status_set
#     dictionary_of_sets["cluster"] = cluster_set
#     dictionary_of_sets["subcluster"] = subcluster_set
#     dictionary_of_sets["accession"] = accession_set
#
#     return dictionary_of_sets

###
