"""Retrieve data from Phamerator and parse into Genome objects

"""










#TODO delete below once confirmed it is functionalized



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
    print "Unsuccessful attempt to connect to the database. Please verify the database, username, and password."
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

#TODO delete above once confirmed it is functionalized





#TODO re-factor this?
# write_out(output_file,"\nDatabase: " + database)
# write_out(output_file,"\nDatabase version: " + db_version)
# write_out(output_file,"\nTotal phages in database before changes: " + str(len(current_genome_data_tuples)))




#TODO delete?
# modified_genome_data_lists = []
# print "Preparing genome data sets from the database..."





#Parse Phamerator genome data
def parse_phamerator_genomes():

    #Now that sets and dictionaries have been made, create a phamerator_data_dict
    #Key = phageID
    #Value = Genome object of Phamerator data
    phamerator_genome_dict = {}

    for genome_tuple in current_genome_data_tuples:

        phamerator_genome = Genome()
        phamerator_genome = parse_phamerator_data(phamerator_genome,genome_tuple)
        phamerator_genome_dict[phamerator_genome.phage_id] = phamerator_genome

    return phamerator_genome_dict



def create_data_sets(phamerator_genome_dict):

    #Create data sets to compare new data with
    #Originally, a phageName_set, phageSequence_set, phageGene_set, and phageGene_dict were implemented, but I ended up not using them here.
    #The phageGene_set get implemented later in the script.
    #The SQL query still returns these values so if I need to re-implement those, I am able to.
    phage_id_set = set()
    phage_host_set = set()
    phage_status_set = set()
    phage_cluster_set = set()
    phage_accession_set = set()
    phage_subcluster_set = set()


    for genome in phamerator_genome_dict.keys():


        phage_id_set.add(genome.phage_id)
        phage_host_set.add(genome.host)
        phage_status_set.add(genome.status)
        phage_cluster_set.add(genome.cluster)
        phage_subcluster_set.add(genome.subcluster)


        #TODO I will need to modify this, since accessions are now "none" instead of ""
        #Only add to the accession set if there was an accession, and not if it was empty "".
        phage_accession_set.add(genome.accession)


    dictionary_of_sets = {}
    dictionary_of_sets["id"] = phage_id_set
    dictionary_of_sets["host"] = phage_host_set
    dictionary_of_sets["status"] = phage_status_set
    dictionary_of_sets["cluster"] = phage_cluster_set
    dictionary_of_sets["subcluster"] = phage_subcluster_set


    return dictionary_of_sets






def main(sql_obj):

    retrieved_phamerator_data = prepare_sql_data(sql_obj)

    phamerator_genome_dict = parse_phamerator_genomes(retrieved_phamerator_data)

    phamerator_genome_sets = create_data_sets(phamerator_genome_dict)

    return phamerator_genome_dict,phamerator_genome_sets

###
