"""Integration tests for misc. functions that interact with PhameratorDB."""


import unittest
from functions import phamerator
from classes import Genome
from classes import Cds
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























    def test_retrieve_genome_data_1(self):
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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_genome_data(sql_handle, "L5")
        with self.subTest():
            self.assertEqual(len(result_list[0].keys()), 12)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "L5")


    def test_retrieve_genome_data_2(self):
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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_genome_data(sql_handle, "EagleEye")
        self.assertEqual(len(result_list), 0)





    def test_create_phamerator_genome_1(self):
        """Verify that a Genome object is constructed correctly for a
        valid PhageID."""


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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        genome_list = phamerator.create_phamerator_genome(sql_handle, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].id, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "ATCG")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)


    def test_create_phamerator_genome_2(self):
        """Verify that an empty Genome object is constructed for an
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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        genome_list = phamerator.create_phamerator_genome(sql_handle, "EagleEye")
        self.assertEqual(len(genome_list), 0)
        # with self.subTest():
        #     self.assertEqual(genome.seq, "")
        # with self.subTest():
        #     self.assertEqual(genome.type, "phamerator")







    def test_retrieve_cds_data_1(self):
        """Verify that a list of data is retrieved for a valid PhageID."""

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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_cds_data(sql_handle, "L5")
        with self.subTest():
            self.assertEqual(len(result_list), 3)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "L5")





    def test_retrieve_cds_data_2(self):
        """Verify that an empty list of data is retrieved
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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_cds_data(sql_handle, "Trixie")
        self.assertEqual(len(result_list), 0)







    def test_retrieve_cds_data_3(self):
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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        result_list = phamerator.retrieve_cds_data(sql_handle)
        self.assertEqual(len(result_list), 5)




    def test_create_cds_1(self):
        """Verify that a Cds object is constructed correctly for a
        valid PhageID."""

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

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        cds_list = phamerator.create_cds(sql_handle, "L5")

        with self.subTest():
            self.assertEqual(len(cds_list), 1)
        with self.subTest():
            self.assertEqual(cds_list[0].id, "L5_001")
        with self.subTest():
            self.assertEqual(cds_list[0].genome_id, "L5")







    def test_create_cds_2(self):
        """Verify that an empty Cds object is constructed for an
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
        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = self.db
        cds_list = phamerator.create_cds(sql_handle, "Trixie")
        self.assertEqual(len(cds_list), 0)





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
