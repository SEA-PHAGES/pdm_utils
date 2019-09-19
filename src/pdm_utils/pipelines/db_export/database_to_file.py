"""Functions for converting a local SQL database query for a selction of phages to formatted files"""
"""Pipeline for converting a database, filtered for some phages, and writing appropriate output files"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from pdm_utils.functions import flat_files, phamerator, basic
import os, sys

# TODO Owen unittest.
def database_to_file(database_name, file_export_format, export_folder_path, phage_name_filter_list = []):
    """Use SQL database to export files of the desired format of selected phage data
    :param database_name:
        Input SQL database name.
    :type database_name: str
    :param file_export_format:
        Input SeqIO file export format.
    :type file_export_format: str
    :param export_folder_path:
        Input the path for the placement of the directory
        of exported files.
    :type export_folder_path:
    :param phage_name_filter_list:
        Input a list of phage names within the selected
        SQL database.
    "type phage_name_filter_list: str[]
    """

    sql_handle = establish_database_connection(database_name)
    seqfeature_file_output\
            (retrieve_seqfeature_from_database\
            (sql_handle, phage_name_filter_list),\
            file_export_format, export_folder_path)

# TODO Owen unittest.
def establish_database_connection(database_name):
    """Creates a mysqlconnectionhandler object and populates its credentials

    :param tag database_name:
        Input SQL database name.
    "type database_name: str
    """

    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle.database = database_name
    sql_handle.get_credentials()
    try:
        sql_handle.validate_credentials
    except:
        print("SQL connection to database {} with username\
                and password failed".format(database_name))
    return sql_handle

# TODO Owen unittest.
def retrieve_seqfeature_from_database (sql_database_handle, phage_name_filter_list = []):
    """Reads a local SQL database and converts it to a SeqRecord list

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :param phage_name_filter_list:
        Input a list of phage names within the selected
        SQL database.
    :type phage_name_filter_list: str[]
    """

    genome_query = "SELECT * FROM phage"
    cds_query = "SELECT * FROM gene"
    genome_list = phamerator.parse_genome_data\
            (sql_database_handle,\
                    phage_id_list = phage_name_filter_list\
                    ,phage_query = genome_query,\
                    gene_query = cds_query)
    database_versions = retrieve_database_version\
            (sql_database_handle)
    seq_record_list = []
    for genome in genome_list:
        set_cds_seqfeatures(genome)
        seqrecord = flat_files.genome_to_seqrecord(genome)
        append_database_version(seqrecord, database_versions)
        seq_record_list.append(seqrecord)


    return seq_record_list

# TODO Owen unittest.
def set_cds_seqfeatures(phage_genome):
    """Helper function that queries for and returns cds data from a SQL database for a specific phage

    :param phage_genome:
        Input a genome object to query cds data for.
    :type phage_genome: genome
    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    """

    try:
        def _sorting_key(cds): return cds.left
        phage_genome.cds_features.sort(key=_sorting_key)
    except:
        print("Genome cds features unable to be sorted")
        pass
    for cds in phage_genome.cds_features:
        cds.set_seqfeature()


# TODO Owen unittest.
def retrieve_database_version(sql_database_handle):
    """Helper function that queries a SQL database for the database version and schema version

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    """

    database_versions_list = phamerator.retrieve_data(sql_database_handle, query='SELECT * FROM version')
    return database_versions_list[0]

# TODO Owen unittest.
def append_database_version(genome_seqrecord, version_data):
    """Helper function that appends the working database version in a comment within a SeqFeature annotation

    :param genome_seqfeature:
        Input a SeqRecord object generated from the working
        SQL database.
    :type genome_seqfeature: SeqRecord
    :param version_data:
        Input a version data dictionary parsed from a SQL database.
    :type version_data: dictionary
    """

    assert len(version_data) >= 2, "Version data dictionary\
        containing SQL database version\
        data does not contain enough values"
    try:
        genome_seqrecord.annotations["comment"] =\
                genome_seqrecord.annotations["comment"] +\
                ("Database Version: {};\
                Schema Version: {}".format\
                (version_data["version"],\
                version_data["schema_version"]),)
    except:
        raise


def seqfeature_file_output(seq_record_list, file_format, input_path):
    """Outputs files with a particuar format from a SeqRecord list

    :param seq_record_list:
        Input a list of SeqRecords.
    :type seq_record_list: SeqRecord[]
    :param file_format:
        Input SeqIO file output format
    :type file_format: str
    :param input_path:
        Input the path for the placement of the directory
        of exported files
    :type input_path: str
    """


    output_dir="database_export_output"
    try:
        os.mkdir(os.path.join(input_path, output_dir))
    except:
        if os.path.exists(os.path.join(input_path, output_dir)):
            pass
        else:
            print("Mkdir function failed to \
                    create database_export_output\
                    directory in {}".format(input_path))
    for record in seq_record_list:
        output_dir="database_export_output/{}.{}".format\
                (record.name, file_format)
        output_path=os.path.join(input_path, output_dir)
        output_handle=open(output_path, "w+")
        SeqIO.write(record, output_handle, file_format)
        output_handle.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Please try again with the following arguments:")
        print("database_to_file.py database_name file_type PhageID1 PhageID2 etc..")
    else:
        if len(sys.argv) == 3:
            database_to_file(database_name = sys.argv[1], file_export_format = sys.argv[2], export_folder_path =os.getcwd())
        else:
            phage_id_list = []
            for args in sys.argv[2:]:
                phage_id_list.append(args)
            database_to_file(database_name = sys.argv[1],\
                    file_export_format = sys.argv[2],\
                    export_folder_path = os.getcwd(),\
                    phage_name_filter_list = phage_id_list)
