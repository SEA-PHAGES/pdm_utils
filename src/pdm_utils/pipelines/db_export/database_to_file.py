"""Functions for converting a local SQL database query for a selction of phages to formatted files"""
"""Pipeline for converting a database, filtered for some phages, and writing appropriate output files"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from pdm_utils.functions import flat_files, phamerator, basic
from functools import singledispatch
from typing import List, Dict
import os, sys, typing

def database_to_file(database_name: str, file_export_format: str, export_folder_path: str, phage_list_input):
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
    :type phage_name_filter_list: list[str]:
    """
   

        phage_name_filter_list =\
            parse_phage_list_input(phage_list_input)

    sql_handle = establish_database_connection(database_name)
    seqfeature_file_output\
            (retrieve_seqrecord_from_database\
            (sql_handle, phage_name_filter_list),\
            file_export_format, export_folder_path, export_dir_name = database_name)

@singledispatch
def parse_phage_list_input(phage_list_input): 
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a list of phage names.
    :type phage_list_input: list[str]:
    :return phage_list_input:
        returns the same phage_list_input:
    """


    return phage_list_input

@parse_phage_list_input.register(str)
def _(phage_list_input):
    """Helper function to populate the filter list for a SQL query
    :param phage_list_input:
        Input a csv file path.
    :type phage_list_input: str:
    :return phage_list:
        Returns a list of phage names
    :type phage_list: list[str]:
    """


    if not os.path.isfile(phage_list_input):
        raise ValueError("File {} is not found".\
                format(phage_list_input))

    phage_list = []
    with open(phage_list_input) as csv:
        csv_reader = csv.reader(csv, delimiter = ",")
        for name in csv_reader[1:]:
            phage_list.append(name[0])

    return phage_list

def establish_database_connection(database_name: str):
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

def retrieve_seqrecord_from_database\
        (sql_database_handle: mysqlconnectionhandler.MySQLConnectionHandler\
        ,phage_name_filter_list: List[str] = []):
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

def set_cds_seqfeatures(phage_genome: genome.Genome):
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
        if phage_genome == None:
            raise TypeError
        print("Genome cds features unable to be sorted")
        pass
    for cds in phage_genome.cds_features:
        cds.set_seqfeature()


def retrieve_database_version\
        (sql_database_handle: mysqlconnectionhandler.MySQLConnectionHandler):
    """Helper function that queries a SQL database for the database version and schema version

    :param sql_database_handle:
        Input a mysqlconnectionhandler object.
    :type sql_database_handle: mysqlconnectionhandler
    :returns:
        database_versions_list(dictionary) is a dictionary
        of size 2 that contains values tied to keys
        "version" and "schema_version"
    """

    database_versions_list = phamerator.retrieve_data(sql_database_handle, query='SELECT * FROM version')
    return database_versions_list[0]

def append_database_version(genome_seqrecord: SeqRecord,\
        version_data: Dict):
    """Helper function that appends the working database version in a comment within a SeqFeature annotation

    :param genome_seqfeature:
        Input a SeqRecord object generated from the working
        SQL database.
    :type genome_seqfeature: SeqRecord
    :param version_data:
        Input a version data dictionary parsed from a SQL database.
    :type version_data: dictionary
    """

    if len(version_data) < 2:
        print("Version data dictionary\
        containing SQL database version\
        data does not contain enough values")
        raise ValueError
    try:
        genome_seqrecord.annotations["comment"] =\
                genome_seqrecord.annotations["comment"] +\
                ("Database Version: {};\
                Schema Version: {}".format\
                (version_data["version"],\
                version_data["schema_version"]),)
    except:
        if genome_seqrecord == None:
            raise TypeError
        raise


def seqfeature_file_output(seq_record_list: List[SeqRecord], file_format: str, input_path: str, export_dir_name: str = "Database"):
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

    if not os.path.exists(input_path):
        print("Path parameter passed to seqfeature_file_output\
            is not a valid path")
        raise ValueError

    try:
        os.mkdir(os.path.join(input_path, export_dir_name))
    except:
        if os.path.exists(os.path.join(input_path, export_dir_name)):
            pass
        else:
            print("Mkdir function failed to \
                    create database_export_output\
                    directory in {}".format(input_path))
            raise ValueError 
    for record in seq_record_list:
        output_dir="{}/{}.{}".format\
                (export_dir_name,\
                record.name, file_format)
        output_path=os.path.join(input_path, output_dir)
        output_handle=open(output_path, "w+")
        SeqIO.write(record, output_handle, file_format)
        output_handle.close()

#To be removed and replaced with run.py functionality
def main(args):
    if len(args) < 2:
        print("Please try again with the following arguments:")
        print("database_to_file.py database_name file_type PhageID1 Phage ID2 etc..")
    else:
        if len(args) == 2:
            database_to_file(database_name = args[0],\
                    file_export_format = args[1],\
                    export_folder_path = os.getcwd())
        else:
            phage_id_list = []
            for arg in args[2:]:
                phage_id_list.append(arg)
            database_to_file(database_name = args[0],\
                        file_export_format = args[1],\
                        export_folder_path = os.getcwd(),\
                        phage_list_input = phage_id_list)

if __name__ == "__main__":
    main_args = ()
    for arg in sys.argv[1:]:
        main_args = main_args + (arg,)
    main(main_args)
    
