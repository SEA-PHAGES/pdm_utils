"""Functions to interact with PhameratorDB."""


from classes import Genome
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
    host_set = set()
    status_set = set()
    cluster_set = set()
    accession_set = set()
    subcluster_set = set()

    for genome_id in genome_dict.keys():

        genome = genome_dict[genome_id]
        phage_id_set.add(genome.phage_id)
        host_set.add(genome.host)
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
    dictionary_of_sets["host"] = host_set
    dictionary_of_sets["status"] = status_set
    dictionary_of_sets["cluster"] = cluster_set
    dictionary_of_sets["subcluster"] = subcluster_set
    dictionary_of_sets["accession"] = accession_set

    return dictionary_of_sets









#If phage Cluster is Singleton, make sure MySQL statement is created correctly
def create_cluster_statement(phage_id,cluster):
    """Create MySQL statement to update Cluster data.
    If Cluster is 'singleton', Cluster is set to NULL."""

    part1 = """UPDATE phage SET Cluster = """
    part3 = """ WHERE PhageID = '""" + phage_id + """';"""

    if cluster.lower() == "singleton":
        part2 = """NULL"""
    else:
        part2 = """'""" + cluster + """'"""

    cluster_statement = part1 + part2 + part3

    return cluster_statement
















# TODO unit test below.

#TODO delete below once confirmed it is functionalized

# def prepare_sql_data(sql_obj):
#
#
#     #Retrieve database version
#     #Retrieve current data in database
#     #0 = PhageID
#     #1 = Name
#     #2 = HostStrain
#     #3 = Sequence
#     #4 = status
#     #5 = Cluster2
#     #6 = DateLastModified
#     #7 = Accession
#     #8 = Subcluster2
#     #9 = AnnotationAuthor
#     #10 = AnnotationQC
#     #11 = RetrieveRecord
#     try:
#         con = mdb.connect(mysqlhost, username, password, database)
#         con.autocommit(False)
#         cur = con.cursor()
#     except:
#         print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
#         output_file.close()
#         sys.exit(1)
#
#     try:
#
#         cur.execute("START TRANSACTION")
#         cur.execute("SELECT version FROM version")
#         db_version = str(cur.fetchone()[0])
#         cur.execute("SELECT PhageID, \
#                             Name, \
#                             HostStrain, \
#                             Sequence,status, \
#                             Cluster2, \
#                             DateLastModified, \
#                             Accession,\
#                             Subcluster2, \
#                             AnnotationAuthor,\
#                             AnnotationQC, \
#                             RetrieveRecord FROM phage")
#         current_genome_data_tuples = cur.fetchall()
#         cur.execute("COMMIT")
#         cur.close()
#         con.autocommit(True)
#     except:
#         mdb_exit("\nUnable to access the database to retrieve genome information.\nNo changes have been made to the database.")
#
#
#     con.close()

    #TODO delete above once confirmed it is functionalized



































#If phage Cluster is Singleton, make sure MySQL statement is created correctly
def create_cluster2_statement(phage_name,cluster):
    cluster2_statement = ""
    if cluster == "singleton":
        cluster2_statement = "UPDATE phage SET Cluster2 = NULL" + " WHERE PhageID = '" + phage_name + "';"
    else:
        cluster2_statement = "UPDATE phage SET Cluster2 = '" + cluster + "' WHERE PhageID = '" + phage_name + "';"
    return cluster2_statement


#If phage Subcluster is empty ("none"), make sure MySQL statement is created correctly
def create_subcluster2_statement(phage_name,subcluster):
    subcluster2_statement = ""
    if subcluster == "none":
        subcluster2_statement = "UPDATE phage SET Subcluster2 = NULL" + " WHERE PhageID = '" + phage_name + "';"
    else:
        subcluster2_statement = "UPDATE phage SET Subcluster2 = '" + subcluster + "' WHERE PhageID = '" + phage_name + "';"
    return subcluster2_statement















def create_update_statement1(table_name, field1, data1, field2, data2):

        update_statement = "UPDATE %s SET %s = '%s' WHERE %s = '%s';" \
                            % (table_name, field1, data1, field2, data2)
        return update_statement

def create_list_of_update_statements(matched_update_object):


        #TODO will need to ensure this is implemented at some point.
        # if genome_data[7] == "none":
        # 	genome_data[7] = ""
        update_statements = []
        update_statements.append(\
                create_update_statement1("phage","HostStrain",host_data,"PhageID",phage_id))
        update_statements.append(\
                create_update_statement1("phage","status",status_data,"PhageID",phage_id))
        update_statements.append(\
                create_update_statement1("phage","Accession",accession_data,"PhageID",phage_id))
        update_statements.append(\
                create_update_statement1("phage","AnnotationAuthor",author_data,"PhageID",phage_id))
        update_statements.append(\
                create_cluster2_statement(genome_data[1],genome_data[3]))
        update_statements.append(\
                create_subcluster2_statement(genome_data[1],genome_data[8]))
        update_statements.append(\
                create_cluster_statement(genome_data[1],\
                assign_cluster_field(genome_data[8],genome_data[3])))

        return update_statements









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













def create_remove_statement1(table_name, field1, data1):

    remove_statement = "DELETE FROM %s WHERE %s = '%s';" \
                        % (table_name, field1, data1)
    return remove_statement





def create_list_of_remove_statements(matched_update_object):


    #Remove actions implemented
    #Prepare to remove any genomes that are not accompanied by a new genome

    removal_statements = []
    for genome_data in remove_data_list:
        removal_statements.append( \
                create_remove_statement1("phage","PhageID",phage_id))
    return removal_statements




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











def create_insert_statement1(cds_feature):

    insert_statement = """INSERT INTO gene (GeneID, PhageID, Start, Stop, Length, Name, TypeID, translation, Orientation, Notes, LocusTag) VALUES ("%s","%s",%s,%s,%s,"%s","%s","%s","%s","%s","%s");""" \
                            % (cds_feature.gene_id,\
                            cds_feature.phage_id,\
                            cds_feature.left_boundary,\
                            cds_feature.right_boundary,\
                            cds_feature.gene_length,\
                            cds_feature.gene_name,\
                            cds_feature.type_id,\
                            cds_feature.translation,\
                            cds_feature.strand,\
                            cds_feature.notes,\
                            cds_feature.locus_tag)
    return insert_statement



def create_list_of_gene_insert_statements(genome_object):

    insert_statements = []

    #Add all updated gene feature data to the add_replace_statements list
    for cds_feature in genome_object.cds_features:

        insert_statements.append(\
            create_insert_statement1(cds_feature))

    return insert_statements














###
