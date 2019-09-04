"""Integration tests for misc. functions that interact with PhameratorDB."""


import unittest
from functions import phamerator
from classes import Genome
from classes import cds
from constants import constants
from classes import MySQLConnectionHandler
import subprocess, os
import pymysql
# import getpass
# import paramiko
# paramiko.util.log_to_file("/tmp/paramiko.log")

# user = getpass.getpass(prompt="mysql username: ")
# pwd = getpass.getpass(prompt="mysql password: ")


# The following integration tests user the 'tester' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = 'tester'
pwd = 'tester'


class TestPhameratorFunctions(unittest.TestCase):



    def setUp(self):
        """In order to test MySQL database-related functions, create a
        new, empty 'test' database and structure it using the
        same schema as in PhameratorDB.
        Each unittest will populate the empty database as needed."""

        # self.db = "test_schema4"
        # self.schema_file = self.db + ".sql"

        self.db = "test_db"
        self.schema_file = "test_schema4.sql"


        self.sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        self.sql_handle.username = user
        self.sql_handle.password = pwd
        self.sql_handle.database = self.db


        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()

        # First, test if a test database already exists within mysql.
        # If there is, delete it so that a fresh test database is installed.
        sql = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA " + \
              "WHERE SCHEMA_NAME = '%s'" % self.db
        cur.execute(sql)
        result = cur.fetchall()

        if len(result) != 0:
            cur.execute("DROP DATABASE %s" % self.db)
            connection.commit()

        # Next, create the database within mysql.
        cur.execute("CREATE DATABASE %s" % self.db)
        connection.commit()
        connection.close()

        # Now import the empty schema from file.
        # Seems like pymysql has trouble with this step, so use subprocess.
        schema_filepath = \
            os.path.join(os.path.dirname(__file__),
                        "test_files/",
                        self.schema_file)

        handle = open(schema_filepath, "r")
        command_string = "mysql -u %s -p%s %s" % (user, pwd, self.db)
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list, stdin = handle)
        handle.close()





    def test_create_phage_id_set_1(self):
        """Retrieve a set of all data from PhageID column."""

        input_phage_ids = ["L5", "Trixie", "D29"]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id in input_phage_ids:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '', 1, 1, '', '%s', 1, 1, 1);" % \
                (id, constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result = phamerator.create_phage_id_set(sql_handle)
        self.assertEqual(len(result), 3)










    def test_create_seq_set_1(self):
        """Retrieve a set of all data from Sequence column."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result = phamerator.create_seq_set(sql_handle)
        with self.subTest():
            self.assertEqual(len(result), 3)
        with self.subTest():
            self.assertTrue("ATCG" in result)





    def test_retrieve_data_1(self):
        """Verify that a dictionary of data is retrieved for a valid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_data(
                        sql_handle, column="PhageID",
                        phage_id_list=["L5"], query=query)
        with self.subTest():
            self.assertEqual(len(result_list[0].keys()), 12)
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql)
        connection.commit()
        connection.close()

        query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                + " Cluster2, DateLastModified, Accession, Subcluster2," \
                + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()


        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_data(
                        sql_handle, query=query)
        self.assertEqual(len(result_list), 5)




    def test_parse_genome_data_1(self):
        """Verify that a Genome object is constructed correctly for a
        valid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', 'ABC123', '', 'Mycobacterium', '%s', " \
                 % (id_and_seq[0], id_and_seq[1]) + \
                " 1, 10.10, 'final', '%s', 1, 1, 1);" \
                % constants.EMPTY_DATE
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                      + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        genome_list = phamerator.parse_genome_data(
                        sql_handle, phage_id_list=["L5"], phage_query=phage_query)
        with self.subTest():
            self.assertEqual(len(genome_list), 1)
        with self.subTest():
            self.assertEqual(genome_list[0].id, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "ATCG")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 0)


    def test_parse_genome_data_2(self):
        """Verify that an empty Genome object list is constructed for an
        invalid PhageID."""


        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]
        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', 'ABC123', '', 'Mycobacterium', '%s', " \
                 % (id_and_seq[0], id_and_seq[1]) + \
                " 1, 10.10, 'final', '%s', 1, 1, 1);" \
                % constants.EMPTY_DATE
            cur.execute(sql)
        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                      + " FROM phage"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        genome_list = phamerator.parse_genome_data(
                          sql_handle, phage_id_list=["EagleEye"],
                          phage_query=phage_query)
        self.assertEqual(len(genome_list), 0)


    def test_parse_genome_data_3(self):
        """Verify that a Genome object with CDS features
        is constructed correctly for a valid PhageID."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', 'ABC123', '', 'Mycobacterium', '%s', " \
                 % (id_and_seq[0], id_and_seq[1]) + \
                " 1, 10.10, 'final', '%s', 1, 1, 1);" \
                % constants.EMPTY_DATE
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                      + " FROM phage"
        gene_query = "SELECT GeneID, PhageID, Start, Stop, Length, Name," \
                     + " TypeID, translation, Orientation, Notes, LocusTag" \
                     + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 3)


    def test_parse_genome_data_4(self):
        """Verify that multiple Genome objects with CDS features
        are constructed correctly for multiple valid PhageIDs."""

        input_phage_ids_and_seqs = [["L5", "ATCG"],
                                    ["Trixie", "AATT"],
                                    ["D29", "GGCC"]]

        input_cds_data = [["L5_001", "L5"],
                          ["L5_002", "L5"],
                          ["L5_003", "L5"],
                          ["TRIXIE_001", "Trixie"],
                          ["TRIXIE_002", "Trixie"]]

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', 'ABC123', '', 'Mycobacterium', '%s', " \
                 % (id_and_seq[0], id_and_seq[1]) + \
                " 1, 10.10, 'final', '%s', 1, 1, 1);" \
                % constants.EMPTY_DATE
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        phage_query = "SELECT PhageID, Name, HostStrain, Sequence, status," \
                      + " Cluster2, DateLastModified, Accession, Subcluster2," \
                      + " AnnotationAuthor, AnnotationQC, RetrieveRecord" \
                      + " FROM phage"
        gene_query = "SELECT GeneID, PhageID, Start, Stop, Length, Name," \
                     + " TypeID, translation, Orientation, Notes, LocusTag" \
                     + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        genome_list = phamerator.parse_genome_data(
                        sql_handle, phage_query=phage_query,
                        gene_query=gene_query)

        genome_dict = {}
        for genome in genome_list:
            genome_dict[genome.id] = genome

        with self.subTest():
            self.assertEqual(len(genome_list), 3)
        with self.subTest():
            self.assertEqual(genome_dict["L5"].seq, "ATCG")
        with self.subTest():
            self.assertEqual(len(genome_dict["L5"].cds_features), 3)
        with self.subTest():
            self.assertEqual(len(genome_dict["Trixie"].cds_features), 2)
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)
        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
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
                                        database = self.db,
                                        cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for id_and_seq in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) " + \
                "VALUES (" + \
                "'%s', '', '', '', '%s', 1, 1, '', '%s', 1, 1, 1);" % \
                (id_and_seq[0], id_and_seq[1], constants.EMPTY_DATE)
            cur.execute(sql1)

        for cds_data in input_cds_data:
            sql2 = "INSERT INTO gene " + \
                "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
                "translation, Orientation, Notes, LocusTag) " + \
                "VALUES " + \
                "('%s', '%s', 1, 100, 1000, '', '', '', 'F', '', '');" % \
                (cds_data[0], cds_data[1])
            cur.execute(sql2)

        connection.commit()
        connection.close()

        query = "SELECT" \
                + " GeneID, PhageID, Start, Stop, Length, Name," \
                + " TypeID, translation, Orientation, Notes, LocusTag" \
                + " FROM gene"

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        cds_list = phamerator.parse_cds_data(sql_handle, query=query)

        with self.subTest():
            self.assertEqual(len(cds_list), 2)




    def tearDown(self):

        connection = pymysql.connect(host = "localhost",
                                        user = user,
                                        password = pwd,
                                        cursorclass = pymysql.cursors.DictCursor)

        cur = connection.cursor()
        cur.execute("DROP DATABASE %s" % self.db)
        connection.commit()
        connection.close()







if __name__ == '__main__':
    unittest.main()
