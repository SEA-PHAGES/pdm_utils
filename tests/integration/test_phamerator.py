"""Integration tests for misc. functions that interact with a MySQL database."""


import unittest
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome
from pdm_utils.classes import cds
from pdm_utils.constants import constants
from pdm_utils.classes import mysqlconnectionhandler
import subprocess, os
import pymysql
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from pathlib import Path
from unittest.mock import patch, Mock
import shutil

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"

unittest_file = Path(__file__)
unittest_dir = unittest_file.parent
schema_file = "test_schema6.sql"
schema_filepath = Path(unittest_dir, "test_files/", schema_file)



class TestPhameratorFunctions1(unittest.TestCase):

    def setUp(self):
        """In order to test MySQL database-related functions, create a
        new, empty 'test' database and structure it using the
        expected schema.
        Each unittest will populate the empty database as needed."""



        self.sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        self.sql_handle.username = user
        self.sql_handle.password = pwd
        self.sql_handle.database = db


        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()

        # First, test if a test database already exists within mysql.
        # If there is, delete it so that a fresh test database is installed.
        sql = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA " + \
              f"WHERE SCHEMA_NAME = '{db}'"
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


    def tearDown(self):
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()
        connection.close()


    def test_create_phage_id_set_1(self):
        """Retrieve a set of all data from PhageID column."""

        input_phage_ids = ["L5", "Trixie", "D29"]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id in input_phage_ids:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id}', '', '', '', '', 1, 1, 'final', " + \
                f"'{constants.EMPTY_DATE}', 1, 1);"

            cur.execute(sql)
        connection.commit()
        connection.close()

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result = phamerator.create_phage_id_set(sql_handle)
        self.assertEqual(len(result), 3)




    def test_create_accession_set_1(self):
        """Retrieve a set of all data from PhageID column."""
        input_phage_ids_and_accs = [["L5", "ABC123"],
                                    ["Trixie", "XYZ456"],
                                    ["D29", "MNO789"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_accs in input_phage_ids_and_accs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_accs[0]}', '{id_and_accs[1]}', '', '', " + \
                f"'', 1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()
        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result = phamerator.create_accession_set(sql_handle)
        self.assertEqual(len(result), 3)




    def test_create_seq_set_1(self):
        """Retrieve a set of all data from Sequence column."""


        input_phage_ids_and_seqs = [["L5", "atcg"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = (
                "INSERT INTO phage (PhageID, Accession, Name, "
                "HostStrain, Sequence, SequenceLength, GC, Status, "
                "DateLastModified, RetrieveRecord, AnnotationAuthor) "
                "VALUES ("
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', "
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);")
            cur.execute(sql)
        connection.commit()
        connection.close()

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result = phamerator.create_seq_set(sql_handle)
        with self.subTest():
            self.assertEqual(len(result), 3)
        with self.subTest():
            self.assertTrue(Seq("ATCG", IUPAC.ambiguous_dna) in result)





    def test_retrieve_data_1(self):
        """Verify that a dictionary of data is retrieved for a valid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, RetrieveRecord" \
                + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["L5"], query=query)
        with self.subTest():
            self.assertEqual(len(result_list[0].keys()), 11)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "L5")

    def test_retrieve_data_2(self):
        """Verify that an empty list is retrieved
        for an invalid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, RetrieveRecord" \
                + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["EagleEye"], query=query)
        self.assertEqual(len(result_list), 0)


    def test_retrieve_data_3(self):
        """Verify that dictionaries of data are retrieved for a list of two
        valid PhageIDs."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, RetrieveRecord" \
                + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["L5","Trixie"], query=query)
        self.assertEqual(len(result_list), 2)


    def test_retrieve_data_4(self):
        """Verify that dictionaries of data are retrieved for a list of three
        valid PhageIDs and one invalid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, RetrieveRecord" \
                + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["L5","Trixie","EagleEye","D29"],
                        query=query)
        self.assertEqual(len(result_list), 3)


    def test_retrieve_data_5(self):
        """Verify that dictionaries of data are retrieved for multiple
        valid PhageIDs when no list is provided."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, RetrieveRecord" \
                + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, query=query)
        self.assertEqual(len(result_list), 3)


    def test_retrieve_data_6(self):
        """Verify that a list of CDS data is retrieved for a valid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, " + \
                "1000, '', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["L5"], query=query)
        with self.subTest():
            self.assertEqual(len(result_list), 3)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "L5")


    def test_retrieve_data_7(self):
        """Verify that an empty list of CDS data is retrieved
        for an invalid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"]]
        input_cds_data = [["L5_001", "L5"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, " + \
                "1000, '', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["Trixie"], query=query)
        self.assertEqual(len(result_list), 0)


    def test_retrieve_data_8(self):
        """Verify that a list of all CDS data is retrieved when no
        PhageID is provided."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()


        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', " + \
                f"1, 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, " + \
                "1000, '', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        result_list = phamerator.retrieve_data(
                        sql_handle, query=query)
        self.assertEqual(len(result_list), 5)




    def test_parse_genome_data_1(self):
        """Verify that a Genome object is constructed correctly for a
        valid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATTC"],
                                    ["D29", "GGCCATT"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                " VALUES (" + \
                f"'{id_and_seq[0]}', 'ABC123', '', 'Mycobacterium', " + \
                f" '{id_and_seq[1]}', " + \
                f" 1, 10.10, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, " + \
                "1000, '', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, RetrieveRecord" \
                      + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        genome_list = phamerator.parse_genome_data(sql_handle,
                        phage_id_list=["L5"], phage_query=phage_query,
                        gnm_type="phamerator")
        with self.subTest():
            self.assertEqual(len(genome_list), 1)
        with self.subTest():
            self.assertEqual(genome_list[0].id, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "ATCG")
        with self.subTest():
            self.assertEqual(genome_list[0].type, "phamerator")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 0)


    def test_parse_genome_data_2(self):
        """Verify that an empty Genome object list is constructed for an
        invalid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATTC"],
                                    ["D29", "GGCCATT"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{id_and_seq[0]}', 'ABC123', '', 'Mycobacterium', " + \
                f"'{id_and_seq[1]}', " + \
                f" 1, 10.10, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql)
        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, RetrieveRecord" \
                      + " FROM phage"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        genome_list = phamerator.parse_genome_data(
                          sql_handle, phage_id_list=["EagleEye"],
                          phage_query=phage_query)
        self.assertEqual(len(genome_list), 0)


    def test_parse_genome_data_3(self):
        """Verify that a Genome object with CDS features
        is constructed correctly for a valid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATTC"],
                                    ["D29", "GGCCATT"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{id_and_seq[0]}', 'ABC123', '', 'Mycobacterium', " + \
                f"'{id_and_seq[1]}', " + \
                f"1, 10.10, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, 1000, " + \
                "'', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, RetrieveRecord" \
                      + " FROM phage"
        gene_query = "SELECT GeneID, PhageID, Start, Stop, Length, Name," \
                     + " Translation, Orientation, Notes, LocusTag" \
                     + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        genome_list = phamerator.parse_genome_data(
                        sql_handle, phage_id_list=["L5"], phage_query=phage_query,
                        gene_query=gene_query)
        with self.subTest():
            self.assertEqual(len(genome_list), 1)
        with self.subTest():
            self.assertEqual(genome_list[0].id, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "ATCG")
        with self.subTest():
            self.assertEqual(genome_list[0].type, "")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 3)
        with self.subTest():
            self.assertEqual(genome_list[0].cds_features[0].genome_length, 4)


    def test_parse_genome_data_4(self):
        """Verify that multiple Genome objects with CDS features
        are constructed correctly for multiple valid PhageIDs."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATTC"],
                                    ["D29", "GGCCATT"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{id_and_seq[0]}', 'ABC123', '', 'Mycobacterium', " + \
                f"'{id_and_seq[1]}', " + \
                f" 1, 10.10, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, 1000, " + \
                "'', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, Status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, RetrieveRecord" \
                      + " FROM phage"
        gene_query = "SELECT GeneID, PhageID, Start, Stop, Length, Name," \
                     + " Translation, Orientation, Notes, LocusTag" \
                     + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        genome_list = phamerator.parse_genome_data(
                        sql_handle, phage_query=phage_query,
                        gene_query=gene_query)

        genome_dict = {}
        for gnm in genome_list:
            genome_dict[gnm.id] = gnm

        with self.subTest():
            self.assertEqual(len(genome_list), 3)
        with self.subTest():
            self.assertEqual(genome_dict["L5"].seq, "ATCG")
        with self.subTest():
            self.assertEqual(len(genome_dict["L5"].cds_features), 3)
        with self.subTest():
            self.assertEqual(
                genome_dict["L5"].cds_features[0].genome_length, 4)
        with self.subTest():
            self.assertEqual(
                genome_dict["L5"].cds_features[1].genome_length, 4)
        with self.subTest():
            self.assertEqual(len(genome_dict["Trixie"].cds_features), 2)
        with self.subTest():
            self.assertEqual(
                genome_dict["Trixie"].cds_features[0].genome_length, 5)
        with self.subTest():
            self.assertEqual(len(genome_dict["D29"].cds_features), 0)




    def test_parse_cds_data_1(self):
        """Verify that a Cds object is constructed correctly for a
        valid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"], ["Trixie", "GGGG"]]
        input_cds_data = [["L5_001", "L5"], ["Trixie_001", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', 1, " + \
                f" 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, 1000, " + \
                "'', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        cds_list = phamerator.parse_cds_data(sql_handle, column="PhageID",
                                             phage_id_list=["L5"], query=query)
        with self.subTest():
            self.assertEqual(len(cds_list), 1)
        with self.subTest():
            self.assertEqual(cds_list[0].id, "L5_001")
        with self.subTest():
            self.assertEqual(cds_list[0].genome_id, "L5")


    def test_parse_cds_data_2(self):
        """Verify that an empty Cds object list is constructed for an
        invalid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"]]
        input_cds_data = [["L5_001", "L5"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', 1, " + \
                f" 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, 1000, " + \
                "'', '', 'F', '', '');"
            cur.execute(sql2)
        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        cds_list = phamerator.parse_cds_data(sql_handle,
                                             column="PhageID",
                                             phage_id_list=["Trixie"],
                                             query=query)
        self.assertEqual(len(cds_list), 0)


    def test_parse_cds_data_3(self):
        """Verify that Cds objects are constructed correctly for multiple
        valid PhageID when the phage_id parameter is not specified."""

        input_phage_ids_and_seqs = [["L5", "ATCG"], ["Trixie", "GGGG"]]
        input_cds_data = [["L5_001", "L5"], ["Trixie_001", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, AnnotationAuthor) " + \
                "VALUES (" + \
                f"'{id_and_seq[0]}', '', '', '', '{id_and_seq[1]}', 1, " + \
                f" 1, 'final', '{constants.EMPTY_DATE}', 1, 1);"
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, " + \
                "Translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                f"('{cds_data[0]}', '{cds_data[1]}', 1, 100, 1000, " + \
                "'', '', 'F', '', '');"
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " Translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db
        cds_list = phamerator.parse_cds_data(sql_handle, query=query)

        with self.subTest():
            self.assertEqual(len(cds_list), 2)




    def test_create_phage_table_insert_1(self):
        """Verify phage table INSERT statement is created correctly."""
        # Note: even though this function returns a string and doesn't
        # actually utilize a MySQL database, this test ensures
        # that the returned statement will function properly in MySQL.
        gnm = genome.Genome()
        gnm.id = "L5"
        gnm.name = "L5_Draft"
        gnm.host_genus = "Mycobacterium"
        gnm.annotation_status = "final"
        gnm.accession = "ABC123"
        gnm.seq = "ATCG"
        gnm.length = 4
        gnm.gc = 0.5001
        gnm.date = constants.EMPTY_DATE
        gnm.retrieve_record = 1
        gnm.annotation_author = 1
        gnm.cluster_subcluster = ""
        gnm.cluster = "singleton"
        gnm.subcluster = "A2"
        statement = phamerator.create_phage_table_insert(gnm)
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()
        query =  ("SELECT PhageID, Accession, Name, "
                 "HostStrain, Sequence, SequenceLength, GC, Status, "
                 "DateLastModified, RetrieveRecord, "
                 "AnnotationAuthor, Cluster, Cluster2, Subcluster2 "
                 "FROM phage WHERE PhageID = 'L5'")
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = ("INSERT INTO phage "
               "(PhageID, Accession, Name, HostStrain, Sequence, "
               "SequenceLength, GC, Status, DateLastModified, RetrieveRecord, "
               "AnnotationAuthor, Cluster, Cluster2, Subcluster2) "
               "VALUES "
               "('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', 'ATCG', "
               f"4, 0.5001, 'final', '{constants.EMPTY_DATE}', 1, "
               "1, NULL, NULL, 'A2');")
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["PhageID"], "L5")
        with self.subTest():
            self.assertEqual(results["Accession"], "ABC123")
        with self.subTest():
            self.assertEqual(results["Name"], "L5_Draft")
        with self.subTest():
            self.assertEqual(results["HostStrain"], "Mycobacterium")
        with self.subTest():
            self.assertEqual(results["Sequence"].decode("utf-8"), "ATCG")
        with self.subTest():
            self.assertEqual(results["SequenceLength"], 4)
        with self.subTest():
            self.assertEqual(results["GC"], 0.5001)
        with self.subTest():
            self.assertEqual(results["Status"], "final")
        with self.subTest():
            self.assertEqual(results["DateLastModified"], constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(results["RetrieveRecord"], 1)
        with self.subTest():
            self.assertEqual(results["AnnotationAuthor"], 1)
        with self.subTest():
            self.assertIsNone(results["Cluster"])
        with self.subTest():
            self.assertIsNone(results["Cluster2"])
        with self.subTest():
            self.assertEqual(results["Subcluster2"], "A2")




    def test_get_phage_table_count_1(self):
        """Verify the correct number of phages is returned when
        the database is empty."""
        count = phamerator.get_phage_table_count(self.sql_handle)
        self.assertEqual(count, 0)


    def test_get_phage_table_count_2(self):
        """Verify the correct number of phages is returned when
        the database contains one genome."""

        # Add the L5 genome to the phage table.
        insert1 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostStrain, Sequence, "
                   "SequenceLength, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster2, Subcluster2) "
                   "VALUES ('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(insert1)
        connection.commit()
        connection.close()
        count = phamerator.get_phage_table_count(self.sql_handle)
        self.assertEqual(count, 1)




    def test_change_version_1(self):
        """Verify the version is incremented by 1."""
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        statement = ("INSERT INTO version (Version, SchemaVersion) "
                    "VALUES (10, 0);")
        cur.execute(statement)
        connection.commit()
        connection.close()
        phamerator.change_version(self.sql_handle)
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        query = "SELECT Version from version;"
        cur.execute(query)
        result = cur.fetchall()
        connection.close()
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 11)

    def test_change_version_2(self):
        """Verify the version is incremented by 5."""
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        statement = ("INSERT INTO version (Version, SchemaVersion) "
                    "VALUES (10, 0);")
        cur.execute(statement)
        connection.commit()
        connection.close()
        phamerator.change_version(self.sql_handle, amount=5)
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        query = "SELECT Version from version;"
        cur.execute(query)
        result = cur.fetchall()
        connection.close()
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 15)

    def test_change_version_3(self):
        """Verify the version is decremented by 5."""
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        statement = ("INSERT INTO version (Version, SchemaVersion) "
                    "VALUES (10, 0);")
        cur.execute(statement)
        connection.commit()
        connection.close()
        phamerator.change_version(self.sql_handle, amount=-5)
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        query = "SELECT Version from version;"
        cur.execute(query)
        result = cur.fetchall()
        connection.close()
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 5)




class TestPhameratorFunctions2(unittest.TestCase):

    def setUp(self):
        """In order to test MySQL database-related functions, create a
        new, empty 'test' database and structure it using the
        expected schema.
        Each unittest will populate the empty database as needed."""


        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()

        # First, test if a test database already exists within mysql.
        # If there is, delete it so that a fresh test database is installed.
        sql = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA " + \
              f"WHERE SCHEMA_NAME = '{db}'"
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

        # Add the L5 genome to the phage table.
        insert1 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostStrain, Sequence, "
                   "SequenceLength, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster2, Subcluster2) "
                   "VALUES ('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(insert1)
        connection.commit()
        connection.close()

        self.std_phage_query = \
            ("SELECT PhageID, Accession, Name, "
             "HostStrain, Sequence, SequenceLength, GC, Status, "
             "DateLastModified, RetrieveRecord, "
             "AnnotationAuthor, Cluster2, Subcluster2 "
             "FROM phage WHERE PhageID = 'L5'")

        self.std_cds_query = \
            ("SELECT GeneID, PhageID, Start, Stop, Length, Name, "
             "Translation, Orientation, Notes, LocusTag FROM gene "
             "WHERE PhageID = 'L5'")


    def tearDown(self):
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()
        connection.close()




    def test_create_gene_table_insert_1(self):
        """Verify gene table INSERT statement is created correctly."""
        # Note: even though this function returns a string and doesn't
        # actually utilize a MySQL database, this test ensures
        # that the returned statement will function properly in MySQL.
        cds1 = cds.Cds()
        cds1.id = "SEA_L5_123"
        cds1.genome_id = "L5"
        cds1.left = 5
        cds1.right = 10
        cds1.translation_length = 20
        cds1.name = "Int"
        cds1.type = "CDS"
        cds1.translation = "ACKLG"
        cds1.strand = "F"
        cds1.processed_description = "integrase"
        cds1.locus_tag = "TAG1"
        insert2 = phamerator.create_gene_table_insert(cds1)
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(insert2)
        connection.commit()
        connection.close()
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_cds_query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = ("INSERT INTO gene "
               "(GeneID, PhageID, Start, Stop, Length, Name, "
               "Translation, Orientation, Notes, LocusTag) "
               "VALUES "
               "('SEA_L5_123', 'L5', 5, 10, 20, 'Int', "
               "'ACKLG', 'F', 'integrase', 'TAG1');")
        with self.subTest():
            self.assertEqual(insert2, exp)
        with self.subTest():
            self.assertEqual(results["GeneID"], "SEA_L5_123")
        with self.subTest():
            self.assertEqual(results["PhageID"], "L5")
        with self.subTest():
            self.assertEqual(results["Start"], 5)
        with self.subTest():
            self.assertEqual(results["Stop"], 10)
        with self.subTest():
            self.assertEqual(results["Length"], 20)
        with self.subTest():
            self.assertEqual(results["Name"], "Int")
        with self.subTest():
            self.assertEqual(results["Translation"], "ACKLG")
        with self.subTest():
            self.assertEqual(results["Orientation"], "F")
        with self.subTest():
            self.assertEqual(results["Notes"].decode("utf-8"), "integrase")
        with self.subTest():
            self.assertEqual(results["LocusTag"], "TAG1")




    def test_create_update_1(self):
        """Verify correct Cluster2 statement is created for a non-singleton."""
        statement = phamerator.create_update(
            "phage", "Cluster2", "B", "PhageID", "L5")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_phage_query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = "UPDATE phage SET Cluster2 = 'B' WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Cluster2"], "B")

    def test_create_update_2(self):
        """Verify correct Cluster2 statement is created for a singleton."""
        statement = phamerator.create_update(
            "phage", "Cluster2", "SINGLETON", "PhageID", "L5")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_phage_query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = "UPDATE phage SET Cluster2 = NULL WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertIsNone(results["Cluster2"])

    def test_create_update_3(self):
        """Verify correct Subcluster2 statement is created for a
        non-empty value."""
        statement = phamerator.create_update(
            "phage", "Subcluster2", "A2", "PhageID", "L5")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_phage_query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = "UPDATE phage SET Subcluster2 = 'A2' WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Subcluster2"], "A2")

    def test_create_cluster_statement_4(self):
        """Verify Gene table statement is created correctly."""
        # First add gene SEA_L5_123 to the database.
        input = ("INSERT INTO gene "
               "(GeneID, PhageID, Start, Stop, Length, Name, "
               "Translation, Orientation, Notes, LocusTag) "
               "VALUES "
               "('SEA_L5_123', 'L5', 5, 10, 20, 'Int', "
               "'ACKLG', 'F', 'integrase', 'TAG1');")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(input)
        connection.commit()
        connection.close()

        # Second run the update statement.
        statement = phamerator.create_update(
            "gene", "Notes", "Repressor", "GeneID", "SEA_L5_123")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()

        # Third, retrieve the current state of the gene table.
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_cds_query)
        results = cur.fetchall()[0]
        cur.close()
        connection.close()
        exp = "UPDATE gene SET Notes = 'Repressor' WHERE GeneID = 'SEA_L5_123';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Notes"].decode("utf-8"), "Repressor")




    def test_create_delete_1(self):
        """Verify correct DELETE statement is created
        for a PhageID in the phage table."""
        # First, retrieve the current state of the phage table.
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_phage_query)
        results1 = cur.fetchall()
        cur.close()
        connection.close()
        results1_phageids = set()
        for dict in results1:
            results1_phageids.add(dict["PhageID"])

        # Second, execute the DELETE statement.
        statement = phamerator.create_delete("phage", "PhageID", "L5")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()

        # Third, retrieve the current state of the phage table.
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_phage_query)
        results2 = cur.fetchall()
        cur.close()
        connection.close()
        results2_phageids = set()
        for dict in results2:
            results2_phageids.add(dict["PhageID"])

        exp = "DELETE FROM phage WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(len(results1_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results2_phageids), 0)






    def test_create_delete_2(self):
        """Verify correct DELETE statement is created
        for a single GeneID in the gene table."""

        # First add two genes to the database.
        input1 = ("INSERT INTO gene "
                  "(GeneID, PhageID, Start, Stop, Length, Name, "
                  "Translation, Orientation, Notes, LocusTag) "
                  "VALUES "
                  "('SEA_L5_1', 'L5', 5, 10, 20, 'LysA', "
                  "'ABCDEF', 'F', 'lysin', 'TAG1');")
        input2 = ("INSERT INTO gene "
                  "(GeneID, PhageID, Start, Stop, Length, Name, "
                  "Translation, Orientation, Notes, LocusTag) "
                  "VALUES "
                  "('SEA_L5_123', 'L5', 5, 10, 20, 'Int', "
                  "'ACKLG', 'F', 'integrase', 'TAG123');")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(input1)
        cur.execute(input2)
        connection.commit()
        connection.close()

        # Second, retrieve the current state of the gene table.
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_cds_query)
        results1 = cur.fetchall()
        cur.close()
        connection.close()
        results1_phageids = set()
        results1_geneids = set()
        for dict in results1:
            results1_phageids.add(dict["PhageID"])
            results1_geneids.add(dict["GeneID"])

        # Third, execute the DELETE statement.
        statement = phamerator.create_delete("gene", "GeneID", "SEA_L5_1")
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(statement)
        connection.commit()
        connection.close()

        # Fourth, retrieve the current state of the gene table.
        connection = pymysql.connect(host="localhost",user=user,
                                     password=pwd, database=db,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(self.std_cds_query)
        results2 = cur.fetchall()
        cur.close()
        connection.close()
        results2_phageids = set()
        results2_geneids = set()
        for dict in results2:
            results2_phageids.add(dict["PhageID"])
            results2_geneids.add(dict["GeneID"])

        exp = "DELETE FROM gene WHERE GeneID = 'SEA_L5_1';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(len(results1_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results1_geneids), 2)
        with self.subTest():
            self.assertEqual(len(results2_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results2_geneids), 1)




class TestPhameratorFunctions3(unittest.TestCase):
    def setUp(self):
        self.database = "Actino_Draft"

    # setup_sql_handle() calls MySQLConnectionHandler.open_connection(),
    # which uses getpass.getpass() to get the username and password.
    @patch("getpass.getpass")
    def test_setup_sql_handle_1(self, getpass_mock):
        """Verify that handle returned with valid info
        when database is provided."""
        getpass_mock.side_effect = [user, pwd]
        sql_handle, msg = phamerator.setup_sql_handle("Actino_Draft")
        with self.subTest():
            self.assertTrue(getpass_mock.called)
        with self.subTest():
            self.assertIsNotNone(sql_handle)
        with self.subTest():
            self.assertTrue(sql_handle.credential_status)
        with self.subTest():
            self.assertTrue(sql_handle._database_status)

    # setup_sql_handle() calls MySQLConnectionHandler.open_connection(),
    # which uses getpass.getpass() to get the username and password.
    @patch("getpass.getpass")
    def test_setup_sql_handle_2(self, getpass_mock):
        """Verify that handle returned with valid info
        when database is not provided."""
        getpass_mock.side_effect = [user, pwd]
        sql_handle, msg = phamerator.setup_sql_handle()
        with self.subTest():
            self.assertTrue(getpass_mock.called)
        with self.subTest():
            self.assertIsNotNone(sql_handle)
        with self.subTest():
            self.assertTrue(sql_handle.credential_status)
        with self.subTest():
            self.assertTrue(sql_handle._database_status)

    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_3(self, mch_mock):
        """Verify that no handle is returned when database is provided and
        credential_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = self.database
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = True
        mch_mock.return_value = sql_handle_mock
        sql_handle, msg = phamerator.setup_sql_handle(self.database)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertIsNone(sql_handle)

    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_4(self, mch_mock):
        """Verify that no handle is returned when database is provided and
        _database_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = self.database
        sql_handle_mock.credential_status = True
        sql_handle_mock._database_status = False
        mch_mock.return_value = sql_handle_mock
        sql_handle, msg = phamerator.setup_sql_handle(self.database)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertIsNone(sql_handle)

    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_5(self, mch_mock):
        """Verify that no handle is returned when database is provided and
        credential_status = False and _database_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = self.database
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = False
        mch_mock.return_value = sql_handle_mock
        sql_handle, msg = phamerator.setup_sql_handle(self.database)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertIsNone(sql_handle)

    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_6(self, mch_mock):
        """Verify that no handle is returned when database is not provided and
        credential_status = False and _database_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = True
        mch_mock.return_value = sql_handle_mock
        sql_handle, msg = phamerator.setup_sql_handle()
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertIsNone(sql_handle)



class TestPhameratorFunctions4(unittest.TestCase):
    def setUp(self):

        self.sql_handle_1 = mysqlconnectionhandler.MySQLConnectionHandler()
        self.sql_handle_1.database = "Actino_Draft"
        self.sql_handle_1.username = user
        self.sql_handle_1.password = pwd

        self.sql_handle_2 = mysqlconnectionhandler.MySQLConnectionHandler()
        self.sql_handle_2.database = "Actino_Draft"
        self.sql_handle_2.username = user
        self.sql_handle_2.password = pwd
        self.sql_handle_2.credential_status = False
        self.sql_handle_2._database_status = False


    @patch("pdm_utils.functions.phamerator.setup_sql_handle")
    def test_connect_to_db_1(self, setup_sql_mock):
        """Verify that sql_handle returned when valid info provided."""
        setup_sql_mock.return_value = (self.sql_handle_1, "")
        sql_handle = phamerator.connect_to_db("Actino_Draft")
        with self.subTest():
            self.assertTrue(setup_sql_mock.called)
        with self.subTest():
            self.assertIsNotNone(sql_handle)

    @patch("sys.exit")
    @patch("pdm_utils.functions.phamerator.setup_sql_handle")
    def test_connect_to_db_2(self, setup_sql_mock, sys_exit_mock):
        """Verify that sys exit is called when sql_handle is None."""
        setup_sql_mock.return_value = (None, "")
        sql_handle = phamerator.connect_to_db("Actino_Draft")
        with self.subTest():
            self.assertTrue(setup_sql_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)



if __name__ == '__main__':
    unittest.main()
