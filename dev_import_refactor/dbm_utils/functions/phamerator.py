"""Functions to interact with PhameratorDB."""


from classes import Genome
from classes import GenomePair
from functions import basic


def parse_phamerator_data(genome, data_tuple):
    """Parses tuple of data derived from a Phamerator database
    and populates a genome object.
    Expected data structure:
    0 = PhageID
    1 = Name
    2 = HostStrain
    3 = Sequence
    4 = status
    5 = Cluster2
    6 = DateLastModified
    7 = Accession
    8 = Subcluster2
    9 = AnnotationAuthor
    10 = AnnotationQC
    11 = RetrieveRecord
    """

    genome.set_phage_id(data_tuple[0])
    genome.phage_name = data_tuple[1]
    genome.set_host(data_tuple[2])
    genome.set_sequence(data_tuple[3])
    genome.status = data_tuple[4]
    genome.set_cluster(data_tuple[5])
    genome.set_subcluster(data_tuple[8])
    genome.set_date_last_modified(data_tuple[6])
    genome.set_accession(data_tuple[7])
    genome.annotation_author = str(data_tuple[9])
    genome.annotation_qc = str(data_tuple[10])
    genome.retrieve_record = str(data_tuple[11])
    genome.type = "phamerator"




def create_phamerator_dict(phamerator_data_tuples):
    """
    Returns a dictionary of Phamerator data retrieved from MySQL query.
    Key = PhageID.
    Value = Genome object containing parsed MySQL data.
    """

    genome_dict = {}
    for genome_tuple in phamerator_data_tuples:
        genome = Genome.Genome()
        parse_phamerator_data(genome,genome_tuple)
        genome_dict[genome.phage_id] = genome

    return genome_dict


def create_data_sets(genome_dict):
    """
    Create sets of all unique values for several fields in the Phamerator data.
    """
    phage_id_set = set()
    host_genus_set = set()
    status_set = set()
    cluster_set = set()
    accession_set = set()
    subcluster_set = set()

    for genome_id in genome_dict.keys():

        genome = genome_dict[genome_id]
        phage_id_set.add(genome.phage_id)
        host_genus_set.add(genome.host_genus)
        status_set.add(genome.status)
        cluster_set.add(genome.cluster)

        # TODO this was not implemented in original import script,
        # so maybe the subcluster 'empty check' is not needed.
        # Only add to the accession set if there was an accession,
        # and not if it was empty.
        if basic.check_empty(genome.subcluster) == False:
            subcluster_set.add(genome.subcluster)

        # Only add to the accession set if there was an accession,
        # and not if it was empty.
        if basic.check_empty(genome.accession) == False:
            accession_set.add(genome.accession)

    dictionary_of_sets = {}
    dictionary_of_sets["phage_id"] = phage_id_set
    dictionary_of_sets["host_genus"] = host_genus_set
    dictionary_of_sets["status"] = status_set
    dictionary_of_sets["cluster"] = cluster_set
    dictionary_of_sets["subcluster"] = subcluster_set
    dictionary_of_sets["accession"] = accession_set

    return dictionary_of_sets


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
    value1 = genome.phage_id
    statements = []

    statements.append(create_update_statement( \
        table, field1, value1, "HostStrain", genome.host_genus))

    statements.append(create_update_statement( \
        table, field1, value1, "status", genome.status))

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





# TODO I may not need this function.
def create_genome_delete_statement(genome):
    """Create a genome-level DELETE statements using data
    in a Genome object."""

    table = "phage"
    field1 = "PhageID"
    value1 = genome.phage_id
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
        (cds_feature.gene_id, \
        cds_feature.phage_id, \
        cds_feature.left_boundary, \
        cds_feature.right_boundary, \
        cds_feature._translation_length, \
        cds_feature.gene_name, \
        cds_feature.type_id, \
        cds_feature.translation, \
        cds_feature.strand, \
        cds_feature.processed_primary_description, \
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
        % (genome.phage_id, \
        genome.accession, \
        genome.phage_name, \
        genome.host_genus, \
        genome.sequence, \
        genome._length, \
        genome._gc, \
        genome.status, \
        genome.date_last_modified, \
        genome.retrieve_record, \
        genome.annotation_qc, \
        genome.annotation_author)

    return statement


def create_genome_insert_statements(genome):
    """Create a collection of genome-level INSERT statements using data
    in a Genome object."""

    table = "phage"
    field1 = "PhageID"
    value1 = genome.phage_id

    statements = []
    statements.append(create_genome_insert_statement(genome))


    statements.append(create_update_statement( \
        table, field1, value1, "Cluster", genome.cluster_subcluster))

    statements.append(create_update_statement( \
        table, field1, value1, "Cluster2", genome.cluster))

    statements.append(create_update_statement( \
        table, field1, value1, "Subcluster2", genome.subcluster))

    return statements



def copy_data_from_phamerator(matched_data_obj, type):
    """If a genome object stored in the DataGroup object has
    attributes that are set to be 'retained' from Phamerator,
    copy any necessary data from the Phamerator genome to the new genome.
    The 'type' parameter indicates the type of genome that may need
    to be populated from Phamerator."""

    if type in matched_data_obj.genome_dict.keys():

        genome1 = matched_data_obj.genome_dict[type]
        genome1.set_retain()

        if genome1._retain:

            if "phamerator" in matched_data_obj.genome_dict.keys():

                genome2 = matched_data_obj.genome_dict["phamerator"]

                # Copy all data that is set to be retained and
                # add to DataGroup object.
                genome_pair = GenomePair.GenomePair()
                genome_pair.genome1 = genome1
                genome_pair.genome2 = genome2
                genome_pair.copy_data("type", genome2.type, genome1.type, "retain")
                matched_data_obj.set_genome_pair(genome_pair, genome1.type, genome2.type)

        # Now record an error if there are still fields
        # that need to be retained.
        genome1.set_retain()
        genome1.check_fields_retained()













# TODO unit test below.









#TODO delete below once confirmed it is functionalized

def retrieve_sql_data(sql_obj):


    #Retrieve database version
    #Retrieve current data in database
    #0 = PhageID
    #1 = Name
    #2 = HostStrain
    #3 = Sequence
    #4 = status
    #5 = Cluster2
    #6 = DateLastModified
    #7 = Accession
    #8 = Subcluster2
    #9 = AnnotationAuthor
    #10 = AnnotationQC
    #11 = RetrieveRecord
    try:
        con = mdb.connect(mysqlhost, username, password, database)
        con.autocommit(False)
        cur = con.cursor()
    except:
        print("Unsuccessful attempt to connect to the database. Please verify the database, username, and password.")
        output_file.close()
        sys.exit(1)

    try:

        cur.execute("START TRANSACTION")
        cur.execute("SELECT version FROM version")
        db_version = str(cur.fetchone()[0])
        cur.execute("SELECT PhageID, \
                            Name, \
                            HostStrain, \
                            Sequence,status, \
                            Cluster2, \
                            DateLastModified, \
                            Accession,\
                            Subcluster2, \
                            AnnotationAuthor,\
                            AnnotationQC, \
                            RetrieveRecord FROM phage")
        current_genome_data_tuples = cur.fetchall()
        cur.execute("COMMIT")
        cur.close()
        con.autocommit(True)
    except:
        mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")


    con.close()




    return list_of_tuples
    # TODO delete above once confirmed it is functionalized








# TODO need to work on this.
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


























###
