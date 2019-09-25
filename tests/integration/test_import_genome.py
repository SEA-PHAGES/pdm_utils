"""Integration tests for the main import pipeline."""



import unittest
import pymysql
import os
import shutil
import subprocess
from unittest.mock import patch, Mock
import argparse
from pdm_utils.pipelines.db_import import import_genome
from pdm_utils.constants import constants
from pdm_utils.classes import bundle, genome, ticket
from pdm_utils.classes import mysqlconnectionhandler as mch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference



# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"


class TestImportGenomeMain1(unittest.TestCase):


    def setUp(self):
        self.schema_file = "test_schema4.sql"
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()

        # First, test if a test database already exists within mysql.
        # If there is, delete it so that a fresh test database is installed.
        sql = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
              "WHERE SCHEMA_NAME = '%s'" % db)
        cur.execute(sql)
        result = cur.fetchall()
        if len(result) != 0:
            cur.execute("DROP DATABASE %s" % db)
            connection.commit()

        # Next, create the database within mysql.
        cur.execute("CREATE DATABASE %s" % db)
        connection.commit()
        connection.close()

        # Now import the empty schema from file.
        # Seems like pymysql has trouble with this step, so use subprocess.
        schema_filepath = \
            os.path.join(os.path.dirname(__file__),
                        "test_files/",
                        self.schema_file)

        handle = open(schema_filepath, "r")
        command_string = "mysql -u %s -p%s %s" % (user, pwd, db)
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list, stdin = handle)
        handle.close()





        self.base_dir = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_import")
        os.mkdir(self.base_dir)

        self.genome_folder = \
            os.path.join(self.base_dir, "genome_folder")
        os.mkdir(self.genome_folder)

        self.test_import_table = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_import_table_1.csv")

        self.test_flat_file1 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_1.gb")

        self.test_flat_file2 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_2.gb")

        self.sql_handle = mch.MySQLConnectionHandler()
        self.sql_handle.database = db
        self.sql_handle.username = user
        self.sql_handle.password = pwd

        self.eval_flags = {
            "check_locus_tag":True,
            "check_description_field":True,
            "check_replace":True,
            "check_trna":True,
            "import_locus_tag":True,
            "check_id_typo":True,
            "check_host_typo":True,
            "check_author":True,
            "check_description":True,
            "check_gene":True
            }
        self.tkt1 = ticket.GenomeTicket()
        self.tkt1.id = 1
        self.tkt1.type = "add"
        self.tkt1.phage_id = "L5"
        self.tkt1.run_mode = "phagesdb"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.host_genus = "Mycobacterium"
        self.tkt1.cluster = "A"
        self.tkt1.subcluster = "A2"
        self.tkt1.annotation_status = "draft"
        self.tkt1.annotation_author = 1
        self.tkt1.annotation_qc = 1
        self.tkt1.retrieve_record = 1
        self.tkt1.accession = "ABC123"

        self.tkt2 = ticket.GenomeTicket()









    def test_prepare_tickets_1(self):
        """Verify dictionary is returned from a correct import table file."""
        eval_flags=constants.RUN_MODE_PHAGESDB
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table,
                        eval_flags=constants.RUN_MODE_PHAGESDB,
                        description_field="product")
        self.assertEqual(len(tkt_dict.keys()), 2)










    def test_prepare_bundle_1(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no phagesdb data."""
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1, genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["add"]
        bndl_tkt = bndl.ticket
        ff_tkt_pair = bndl.genome_pair_dict["flat_file_add"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(bndl.id, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "L5")






    def test_prepare_bundle_2(self):
        """Verify bundle is returned from a flat file with:
        no record."""
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file2,
                    ticket_dict=tkt_dict, id=1, genome_id_field="organism_name")
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)




    def test_prepare_bundle_3(self):
        """Verify bundle is returned from a flat file with:
        one record, no ticket."""
        tkt_dict = {"L5x":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1, genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)




    def test_prepare_bundle_4(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, with phagesdb data."""
        self.tkt1.host_genus = "retrieve"
        self.tkt1.cluster = "B"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1, genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        pdb_gnm = bndl.genome_dict["phagesdb"]
        ff_tkt_pair = bndl.genome_pair_dict["flat_file_add"]
        ff_pdb_pair = bndl.genome_pair_dict["flat_file_phagesdb"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")






    def test_prepare_bundle_5(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with phamerator data,
        and no phagesdb data."""

        l5_data = ["L5", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', '%s', '%s', '%s', '%s', " \
                 % (data[0], data[1], data[2], data[3], data[4] ) + \
                " %s, %s, '%s', '%s', %s, %s, %s);" \
                % (data[5], data[6], data[7], data[8], data[9],
                data[10], data[11])
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.cluster = "B"
        self.tkt1.annotation_author = "retain"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, sql_handle=self.sql_handle, id=1,
                    genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        pmr_gnm = bndl.genome_dict["phamerator"]
        ff_tkt_pair = bndl.genome_pair_dict["flat_file_add"]
        ff_pmr_pair = bndl.genome_pair_dict["flat_file_phamerator"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.annotation_author, 1)






    def test_prepare_bundle_6(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, no phamerator data,
        and no phagesdb data."""

        l5_data = ["L5x", "ABC123", "L5_Draft", "Mycobacterium", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Gordonia", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', '%s', '%s', '%s', '%s', " \
                 % (data[0], data[1], data[2], data[3], data[4] ) + \
                " %s, %s, '%s', '%s', %s, %s, %s);" \
                % (data[5], data[6], data[7], data[8], data[9],
                data[10], data[11])
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.annotation_author = "retain"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, sql_handle=self.sql_handle, id=1,
                    genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.annotation_author, "retain")






    def test_prepare_bundle_7(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with phamerator data,
        no MySQL connection handle, and no phagesdb data."""

        l5_data = ["L5", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = \
                "INSERT INTO phage (PhageID, Accession, Name, " + \
                "HostStrain, Sequence, SequenceLength, GC, status, " + \
                "DateLastModified, RetrieveRecord, AnnotationQC, " + \
                "AnnotationAuthor) VALUES (" + \
                "'%s', '%s', '%s', '%s', '%s', " \
                 % (data[0], data[1], data[2], data[3], data[4] ) + \
                " %s, %s, '%s', '%s', %s, %s, %s);" \
                % (data[5], data[6], data[7], data[8], data[9],
                data[10], data[11])
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.annotation_author = "retain"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1, genome_id_field="organism_name")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.annotation_author, "retain")




    def tearDown(self):

        # Remove all contents in the directory created for the test.
        shutil.rmtree(self.base_dir)

        # Remove the MySQL database created for the test.
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute("DROP DATABASE %s" % db)
        connection.commit()
        connection.close()






class TestImportGenomeMain2(unittest.TestCase):

    def setUp(self):
        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")

        self.test_directory1 = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/input_folder")
        os.mkdir(self.test_directory1)

        self.sql_handle_1 = mch.MySQLConnectionHandler()
        # self.sql_handle_1.database = "Actino_Draft"
        self.sql_handle_1.username = user
        self.sql_handle_1.password = pwd

        self.sql_handle_2 = mch.MySQLConnectionHandler()
        self.sql_handle_2.database = "Actino_Draft"
        self.sql_handle_2.username = user
        self.sql_handle_2.password = pwd
        self.sql_handle_2.credential_status = False
        self.sql_handle_2._database_status = False

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("pipeline")
        self.args1 = self.parser.parse_args(["import"])


    @patch("pdm_utils.pipelines.db_import.import_genome.input_output")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_run_import_1(self, mch_mock, input_output_mock):
        """Verify that correct args calls input_output."""
        args_list = ["run.py",
                     "import",
                     "Actino_Draft",
                     self.test_directory1,
                     self.test_filepath1]
        mch_mock.return_value = self.sql_handle_1
        import_genome.run_import(args_list)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(input_output_mock.called)

    @patch("pdm_utils.pipelines.db_import.import_genome.input_output")
    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_run_import_2(self, mch_mock, input_output_mock,
                          sys_exit_mock):
        """Verify that input_output is not called due to
        invalid sql credential_status."""
        args_list = ["run.py",
                     "import",
                     "Actino_Draft",
                     self.test_directory1,
                     self.test_filepath1]

        sql_handle_mock = Mock()
        sql_handle_mock.database = "Actino_Draft"
        sql_handle_mock.credential_status = True
        sql_handle_mock._database_status = False
        mch_mock.return_value = sql_handle_mock

        import_genome.run_import(args_list)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(input_output_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.db_import.import_genome.input_output")
    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_run_import_3(self, mch_mock, input_output_mock,
                          sys_exit_mock):
        """Verify that input_output is not called due to
        invalid sql database_status."""
        args_list = ["run.py",
                     "import",
                     "Actino_Draft",
                     self.test_directory1,
                     self.test_filepath1]

        sql_handle_mock = Mock()
        sql_handle_mock.database = "Actino_Draft"
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = True
        mch_mock.return_value = sql_handle_mock

        import_genome.run_import(args_list)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(input_output_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)

    @patch("pdm_utils.pipelines.db_import.import_genome.input_output")
    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_run_import_4(self, mch_mock, input_output_mock,
                          sys_exit_mock):
        """Verify that input_output is not called due to invalid folder."""
        args_list = ["run.py",
                     "import",
                     "Actino_Draft",
                     self.test_directory1 + "asdf",
                     self.test_filepath1]
        mch_mock.return_value = self.sql_handle_1
        import_genome.run_import(args_list)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(input_output_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.db_import.import_genome.input_output")
    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_run_import_5(self, mch_mock, input_output_mock,
                          sys_exit_mock):
        """Verify that input_output is not called due to invalid file."""
        args_list = ["run.py",
                     "import",
                     "Actino_Draft",
                     self.test_directory1,
                     self.test_filepath1 + "asdf"]
        mch_mock.return_value = self.sql_handle_1
        import_genome.run_import(args_list)
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(input_output_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)



    def tearDown(self):
        shutil.rmtree(self.test_directory1)




if __name__ == '__main__':
    unittest.main()
