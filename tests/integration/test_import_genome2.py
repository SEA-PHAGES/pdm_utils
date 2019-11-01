"""Integration tests for the entire import pipeline."""

import csv
# import time
from datetime import datetime
import unittest
import pymysql
import shutil
from unittest.mock import patch
from pathlib import Path
from pdm_utils import run
import subprocess
from Bio import SeqIO
from pdm_utils.classes import mysqlconnectionhandler as mch



# import argparse
# from pdm_utils.pipelines.db_import import import_genome
# from pdm_utils.constants import constants
# from pdm_utils.functions import basic
# from pdm_utils.classes import bundle, genome, ticket, eval
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
# from Bio.SeqFeature import ExactPosition, Reference
# from unittest.mock import patch
# import getpass

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"
unittest_file = Path(__file__)
unittest_dir = unittest_file.parent
schema_file = "test_schema5.sql"
schema_filepath = Path(unittest_dir, "test_files/", schema_file)
base_flat_file = "test_flat_file_10.gb"
base_flat_file_path = Path(unittest_dir, "test_files/", base_flat_file)


def create_new_db(schema_file, db, user, pwd):
    """Creates a new, empty database."""

    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()

    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    sql = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
          f"WHERE SCHEMA_NAME = '{db}'")
    cur.execute(sql)
    result = cur.fetchall()
    if len(result) != 0:
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()

    # Next, create the database within mysql.
    cur.execute(f"CREATE DATABASE {db}")
    connection.commit()
    connection.close()

    # Now import the empty schema from file.
    # Seems like pymysql has trouble with this step, so use subprocess.
    handle = open(schema_filepath, "r")
    command_string = f"mysql -u {user} -p{pwd} {db}"
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list, stdin = handle)
    handle.close()




def remove_db(db, user, pwd):
    """Remove the MySQL database created for the test."""
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(f"DROP DATABASE {db}")
    connection.commit()
    connection.close()


def insert_data_into_phage_table(db, user, pwd, data_dict):
    """Insert data into the phage table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO phage "
        "(PhageID, Accession, Name, "
        "HostStrain, Sequence, SequenceLength, GC, status, "
        "DateLastModified, RetrieveRecord, AnnotationAuthor, "
        "Cluster, Cluster2, Subcluster2) "
        "VALUES ("
        f"'{data_dict['phage_id']}', '{data_dict['accession']}', "
        f"'{data_dict['name']}', '{data_dict['host_strain']}', "
        f"'{data_dict['sequence']}', {data_dict['sequence_length']}, "
        f"{data_dict['gc']}, '{data_dict['status']}', "
        f"'{data_dict['date_last_modified']}', "
        f"{data_dict['retrieve_record']}, "
        f"{data_dict['annotation_author']}, '{data_dict['cluster']}', "
        f"'{data_dict['cluster2']}', '{data_dict['subcluster2']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


def insert_data_into_gene_table(db, user, pwd, data_dict):
    """Insert data into the gene table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO gene "
        "(GeneID, PhageID, Start, Stop, Length, Name, "
        "translation, Orientation, Notes, LocusTag) "
        "VALUES ("
        f"'{data_dict['gene_id']}', '{data_dict['phage_id']}', "
        f"{data_dict['start']}, {data_dict['stop']}, "
        f"{data_dict['length']}, '{data_dict['name']}', "
        f"'{data_dict['translation']}', '{data_dict['orientation']}', "
        f"'{data_dict['notes']}', '{data_dict['locus_tag']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


phage_table_query = (
    "SELECT "
    "PhageID, Accession, Name, "
    "HostStrain, Sequence, SequenceLength, GC, status, "
    "DateLastModified, RetrieveRecord, AnnotationAuthor, "
    "Cluster, Cluster2, Subcluster2 "
    "FROM phage;"
    )

gene_table_query = (
    "SELECT "
    "GeneID, PhageID, Start, Stop, Length, Name, "
    "translation, Orientation, Notes, LocusTag "
    "FROM gene;"
)



def get_sql_data(db, user, pwd, query):
    """Get data from the phage table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(query)
    results = cur.fetchall()
    cur.close()
    connection.close()
    return results


l5_phage_table_data = {
    "phage_id": "L5",
    "accession": "ABC123",
    "name": "L5",
    "host_strain": "Mycobacterium",
    "sequence": "AAAAAAAAAA",
    "sequence_length": 10,
    "gc": 0,
    "status": "final",
    "date_last_modified": datetime.strptime('1/1/0001', '%m/%d/%Y'),
    "retrieve_record": 1,
    "annotation_author": 1,
    "cluster": "A",
    "cluster2": "A",
    "subcluster2": "A2"
    }



trixie_phage_table_data = {
    "phage_id": "Trixie",
    "accession": "BCD456",
    "name": "Trixie",
    "host_strain": "Gordonia",
    "sequence": "GGGGGGGGGGGGGGGGGGGG",
    "sequence_length": 20,
    "gc": 1,
    "status": "final",
    "date_last_modified": datetime.strptime('1/1/2000', '%m/%d/%Y'),
    "retrieve_record": 1,
    "annotation_author": 1,
    "cluster": "A",
    "cluster2": "A",
    "subcluster2": "A3"
    }

redrock_phage_table_data = {
    "phage_id": "RedRock",
    "accession": "BCD456",
    "name": "RedRock_Draft",
    "host_strain": "Arthrobacter",
    "sequence": "CCCCCCCCCCCCCCC",
    "sequence_length": 15,
    "gc": 1,
    "status": "draft",
    "date_last_modified": datetime.strptime('1/1/2010', '%m/%d/%Y'),
    "retrieve_record": 1,
    "annotation_author": 1,
    "cluster": "B",
    "cluster2": "B",
    "subcluster2": "B1"
    }


d29_phage_table_data = {
    "phage_id": "D29",
    "accession": "XYZ123",
    "name": "D29",
    "host_strain": "Microbacterium",
    "sequence": "ATGCATGCATGCATGC",
    "sequence_length": 16,
    "gc": 0.5,
    "status": "unknown",
    "date_last_modified": datetime.strptime('1/1/0001', '%m/%d/%Y'),
    "retrieve_record": 0,
    "annotation_author": 0,
    "cluster": "C",
    "cluster2": "C",
    "subcluster2": "C1"
    }

l5_gene_table_data_1 = {
    "gene_id": "L5_0001",
    "phage_id": "L5",
    "start": 100,
    "stop": 1100,
    "length": 1000,
    "name": "1",
    "translation": "ACTGC",
    "orientation": "F",
    "notes": "int",
    "locus_tag": "SEA_L5_0001"
    }



trixie_gene_table_data_1 = {
    "gene_id": "TRIXIE_0001",
    "phage_id": "Trixie",
    "start": 100,
    "stop": 1100,
    "length": 1000,
    "name": "1",
    "translation": "ACTGC",
    "orientation": "F",
    "notes": "int",
    "locus_tag": "SEA_TRIXIE_0001"
    }


def get_seq(filepath):
    """Get genome sequence from a flat file so that it can be added
    to a MySQL database for testing."""
    seqrecord = SeqIO.read(filepath, "genbank")
    return seqrecord.seq


def create_min_tkt_dict(old_tkt):

    min_keys = set(["id", "type", "phage_id"])
    new_tkt = {}
    for key in min_keys:
        new_tkt[key] = old_tkt[key]
    return new_tkt

def set_data(old_tkt, keys, value):
    """Set selected keys to certain values."""
    new_tkt = old_tkt.copy()
    for key in keys:
        new_tkt[key] = value
    return new_tkt



l5_ticket_data_complete = {
    "id": 1,
    "type": "add",
    "phage_id": "L5",
    "host_genus": "Mycobacterium",
    "cluster": "A",
    "subcluster": "A2",
    "accession": "ABC123",
    "description_field": "product",
    "annotation_status": "draft",
    "annotation_author": 1,
    "retrieve_record": 1,
    "run_mode": "phagesdb",
    }

l5_ticket_data_min = create_min_tkt_dict(l5_ticket_data_complete)

valid_retrieve_keys = set(["host_genus", "cluster", "subcluster", "accession"])
l5_ticket_data_retrieve = set_data(l5_ticket_data_complete, valid_retrieve_keys, "retrieve")

valid_retain_keys = set(["host_genus", "cluster", "subcluster",
                         "accession", "annotation_author", "retrieve_record"])






def create_import_table(list_of_data_dicts, file_path):
    """Create an import table."""

    headers = set()
    for dict in list_of_data_dicts:
        headers = headers | dict.keys()

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    with open(file_path, "w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        file_writer.writerow(headers_dict)
        for data_dict in list_of_data_dicts:
            file_writer.writerow(data_dict)






class TestImportGenomeMain1(unittest.TestCase):


    def setUp(self):

        create_new_db(schema_file, db, user, pwd)


        self.base_dir = Path(unittest_dir, "test_wd/test_import")
        self.base_dir.mkdir()

        self.import_table = Path(self.base_dir, "import_table.csv")
        self.genome_folder = Path(self.base_dir, "genome_folder")
        self.genome_folder.mkdir()

        self.output_folder = Path(self.base_dir, "output_folder")
        self.output_folder.mkdir()

        self.log_file = Path(self.output_folder, "import_log.txt")

        self.test_flat_file1 = Path(unittest_dir,
                                    "test_files/test_flat_file_1.gb")

        self.test_import_table1 = Path(unittest_dir,
                                    "test_files/test_import_table_1.csv")


        self.sql_handle = mch.MySQLConnectionHandler()
        self.sql_handle.database = db
        self.sql_handle.username = user
        self.sql_handle.password = pwd


    def tearDown(self):
        shutil.rmtree(self.base_dir)
        remove_db(db, user, pwd)



    @patch("pdm_utils.pipelines.db_import.import_genome.setup_sql_handle")
    def test_import_pipeline_1(self, setup_sql_mock):
        """Test pipeline with:
        valid add ticket,
        valid flat file."""
        setup_sql_mock.return_value = self.sql_handle
        insert_data_into_phage_table(db, user, pwd, d29_phage_table_data)
        insert_data_into_phage_table(db, user, pwd, redrock_phage_table_data)
        insert_data_into_phage_table(db, user, pwd, trixie_phage_table_data)
        # insert_data_into_gene_table(db, user, pwd, trixie_gene_table_data_1)

        # results = get_sql_data(db, user, pwd, phage_table_query)
        # print(len(results))
        # results2 = get_sql_data(db, user, pwd, gene_table_query)
        # print(len(results2))
        # input("paused")
        # seq = get_seq(self.test_flat_file1)
        # l5_phage_table_data["sequence"] = seq
        # l5_phage_table_data["sequence_length"] = len(seq)
        # print(trixie_phage_table_data["sequence"][:25])
        # print(trixie_phage_table_data["sequence_length"])
        # input("pause")
        # insert_data_into_phage_table(db, user, pwd, trixie_phage_table_data)
        # input("pause2")
        # keys = set(["host_genus", "cluster", "subcluster"])
        # l5_ticket_data_retrieve = set_data(l5_ticket_data_complete, keys, "retrieve")


        list_of_data_dicts = [l5_ticket_data_complete]
        create_import_table(list_of_data_dicts, self.import_table)
        shutil.copy(self.test_flat_file1, self.genome_folder)
        unparsed_args = ["run.py", "import_dev", db,
                         str(self.genome_folder),
                         str(self.import_table),
                         "-g", "_organism_name",
                         "-p",
                         "-r", "phagesdb",
                         "-d", "product",
                         "-o", str(self.output_folder),
                         "-l", str(self.log_file)
                         ]
        # input("paused for run")
        run.main(unparsed_args)

        phage_table_results = get_sql_data(db, user, pwd, phage_table_query)
        gene_table_results = get_sql_data(db, user, pwd, gene_table_query)

        # TODO in progress below.
        # with self.subTest():
        #     self.assertEqual(len(phage_table_results), 4)
        # with self.subTest():
        #     self.assertTrue(len(gene_table_results) > 0)











if __name__ == '__main__':
    unittest.main()
