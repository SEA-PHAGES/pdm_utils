"""Integration tests for the main import pipeline."""


import csv
import time
import unittest
import pymysql
import os
import shutil
import subprocess
from unittest.mock import patch, Mock
import argparse
from pdm_utils.pipelines.db_import import import_genome
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.classes import bundle, genome, genomepair, ticket, eval, cds, source
from pdm_utils.classes import mysqlconnectionhandler as mch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference
from unittest.mock import patch
from pathlib import Path
import getpass


def get_errors_from_dict(dict):
    """Iterate through a dictionary produce from bundle.get_evaluations()."""
    errors = 0
    for key in dict:
        sub_dict = dict[key]
        for sub_key in sub_dict:
            evl_list = sub_dict[sub_key]
            for evl in evl_list:
                if evl.status == "error":
                    errors += 1
    return errors

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"


class TestImportGenomeMain1(unittest.TestCase):


    def setUp(self):
        self.schema_file = "test_schema6.sql"
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
        self.unittest_file = Path(__file__)
        self.unittest_dir = self.unittest_file.parent
        self.test_files_path = Path(self.unittest_dir, "test_files")
        schema_filepath = Path(self.test_files_path, self.schema_file)

        handle = open(schema_filepath, "r")
        command_string = f"mysql -u {user} -p{pwd} {db}"
        command_list = command_string.split(" ")
        proc = subprocess.check_call(command_list, stdin = handle)
        handle.close()



        self.base_dir = Path(self.unittest_dir,"test_wd/test_import")
        self.base_dir.mkdir()

        self.genome_folder = Path(self.base_dir,"genome_folder")
        self.genome_folder.mkdir()

        self.test_flat_file1 = Path(self.test_files_path, "test_flat_file_1.gb")
        self.test_flat_file2 = Path(self.test_files_path, "test_flat_file_2.gb")

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

        self.data_dict = {}
        self.data_dict["host_genus"] = "Arthrobacter"
        self.data_dict["cluster"] = "B"
        self.data_dict["subcluster"] = "B2"
        self.data_dict["annotation_status"] = "draft"
        self.data_dict["annotation_author"] = 1
        self.data_dict["retrieve_record"] = 1
        self.data_dict["accession"] = "ABC123"

        self.tkt1 = ticket.GenomeTicket()
        self.tkt1.id = 1
        self.tkt1.type = "add"
        self.tkt1.phage_id = "L5"
        self.tkt1.run_mode = "phagesdb"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.data_dict = self.data_dict

        self.tkt2 = ticket.GenomeTicket()

        self.id_dict = constants.PHAGE_ID_DICT


    def tearDown(self):

        # Remove all contents in the directory created for the test.
        shutil.rmtree(self.base_dir)

        # Remove the MySQL database created for the test.
        connection = pymysql.connect(host="localhost",
                                     user=user,
                                     password=pwd,
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = connection.cursor()
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()
        connection.close()

    def test_prepare_bundle_1(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no phagesdb data."""
        # Omit cluster to verify that only attributes in the data_add set
        # are copied.
        self.tkt1.data_add = set(["host_genus", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])

        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket")
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(bndl.id, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Arthrobacter")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "")
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "L5")
        with self.subTest():
            self.assertTrue(ff_gnm._cds_processed_descriptions_tally > 0)
        with self.subTest():
            self.assertEqual(ff_gnm._cds_processed_descriptions_tally,
                             ff_gnm._cds_processed_products_tally)


    def test_prepare_bundle_2(self):
        """Verify bundle is returned from a flat file with:
        no record."""
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file2,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name")
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)


    def test_prepare_bundle_3(self):
        """Verify bundle is returned from a flat file with:
        one record, no ticket."""
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        tkt_dict = {"L5x":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertIsNone(bndl.ticket)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)


    def test_prepare_bundle_4(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no data_add, no phagesdb data."""
        self.tkt1.data_add = set()
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)


    def test_prepare_bundle_5(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, with phagesdb data."""
        # Use cluster and host_genus to confirm that only attributes
        # within the data_retrieve set are copied.
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_dict["host_genus"] = "retrieve"
        self.tkt1.data_retrieve = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    retrieve_ref="phagesdb")
        ff_gnm = bndl.genome_dict["flat_file"]
        pdb_gnm = bndl.genome_dict["phagesdb"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")
        with self.subTest():
            self.assertEqual(pdb_gnm.cluster, "A")


    def test_prepare_bundle_6(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with phamerator data,
        and no phagesdb data."""
        # Use host_genus and accession to confirm that only attributes
        # in the data_retain set are copied.
        l5_data = ["L5", "EFG789", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
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
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{data[0]}', '{data[1]}', '{data[2]}', " + \
                f"'{data[3]}', '{data[4]}', " + \
                f" {data[5]}, {data[6]}, '{data[7]}', " + \
                f"'{data[8]}', {data[9]}, {data[10]});"
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    retain_ref="phamerator",
                    sql_handle=self.sql_handle)
        ff_gnm = bndl.genome_dict["flat_file"]
        pmr_gnm = bndl.genome_dict["phamerator"]
        ff_pmr_pair = bndl.genome_pair_dict["flat_file_phamerator"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(ff_gnm.accession, "ABC123")
        with self.subTest():
            self.assertEqual(pmr_gnm.accession, "EFG789")


    def test_prepare_bundle_7(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, no phamerator data,
        and no phagesdb data."""

        l5_data = ["L5x", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Gordonia", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
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
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{data[0]}', '{data[1]}', '{data[2]}', " + \
                f"'{data[3]}', '{data[4]}', " + \
                f" {data[5]}, {data[6]}, '{data[7]}', " + \
                f"'{data[8]}', {data[9]}, {data[10]});"
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    sql_handle=self.sql_handle)
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "")


    def test_prepare_bundle_8(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with phamerator data,
        no MySQL connection handle, and no phagesdb data."""

        l5_data = ["L5", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
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
                "HostStrain, Sequence, SequenceLength, GC, Status, " + \
                "DateLastModified, RetrieveRecord, " + \
                "AnnotationAuthor) VALUES (" + \
                f"'{data[0]}', '{data[1]}', '{data[2]}', " + \
                f"'{data[3]}', '{data[4]}', " + \
                f" {data[5]}, {data[6]}, '{data[7]}', " + \
                f"'{data[8]}', {data[9]}, {data[10]});"
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)


    def test_prepare_bundle_9(self):
        """Verify bundle is returned with the genome id converted using
        the id dictionary and with the cluster_subcluster set."""
        self.id_dict["L5"] = "new_id"
        self.tkt1.phage_id = "new_id"
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])

        tkt_dict = {"new_id":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket",
                    id_conversion_dict = self.id_dict)
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(ff_gnm.id, "new_id")
        with self.subTest():
            self.assertEqual(ff_gnm.name, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")
        with self.subTest():
            self.assertEqual(ff_gnm.subcluster, "B2")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster_subcluster, "B2")
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "new_id")


    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_prepare_bundle_10(self, ask_mock, choose_mock):
        """Verify bundle is returned with the CDS descriptions derived
        from 'function' instead of 'product', after it is
        interactively selected using interactive = True."""
        self.tkt1.data_add = set(["host_genus", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        ask_mock.side_effect = [False]
        choose_mock.side_effect = ["function"]
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket",
                    interactive=True)
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertTrue(ff_gnm._cds_processed_descriptions_tally != \
                            ff_gnm._cds_processed_products_tally)
        with self.subTest():
            self.assertEqual(ff_gnm._cds_processed_descriptions_tally,
                             ff_gnm._cds_processed_functions_tally)





class TestImportGenomeMain2(unittest.TestCase):

    def setUp(self):
        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")
        self.test_filepath1 = Path(self.test_filepath1)

        self.test_directory1 = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_dir")
        self.test_directory1 = Path(self.test_directory1)
        self.test_directory1.mkdir()

        # Minimum args list
        self.args_list = ["run.py",
                          "import",
                          "Actino_Draft",
                          str(self.test_directory1),
                          str(self.test_filepath1)]


    def tearDown(self):
        shutil.rmtree(self.test_directory1)




    def test_parse_args_1(self):
        """Verify args when minimum args_list is provided."""
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.database, "Actino_Draft")
        with self.subTest():
            self.assertEqual(args.input_folder, self.test_directory1)
        with self.subTest():
            self.assertEqual(args.import_table, self.test_filepath1)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "_organism_name")
        with self.subTest():
            self.assertFalse(args.prod_run)
        with self.subTest():
            self.assertEqual(args.run_mode, "phagesdb")
        with self.subTest():
            self.assertEqual(args.description_field, "product")
        with self.subTest():
            self.assertEqual(args.output_folder, Path("/tmp/"))
        with self.subTest():
            self.assertEqual(args.log_file, "import.log")
        with self.subTest():
            self.assertFalse(args.interactive)


    @patch("sys.exit")
    def test_parse_args_2(self, sys_exit_mock):
        """Verify sys exit when too few arguments are provided."""
        self.args_list.pop()
        args = import_genome.parse_args(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    def test_parse_args_3(self):
        """Verify parsed args when all args are explicitly provided
        using short name."""
        output_folder = "/path/to/output/"
        log_file = "logfile.txt"
        self.args_list.extend(["-g", "FILENAME",
                               "-p",
                               "-r", "PECAAN",
                               "-d", "FUNCTION",
                               "-o", output_folder,
                               "-l", log_file,
                               "-i"
                               ])
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "filename")
        with self.subTest():
            self.assertTrue(args.prod_run)
        with self.subTest():
            self.assertEqual(args.run_mode, "pecaan")
        with self.subTest():
            self.assertEqual(args.description_field, "function")
        with self.subTest():
            self.assertEqual(args.output_folder, Path(output_folder))
        with self.subTest():
            self.assertEqual(args.log_file, log_file)
        with self.subTest():
            self.assertTrue(args.interactive)


    def test_parse_args_4(self):
        """Verify parsed args when all args are explicitly provided
        using long name."""
        output_folder = "/path/to/output/"
        log_file = "logfile.txt"
        self.args_list.extend(["--genome_id_field", "FILENAME",
                               "--prod_run",
                               "--run_mode", "PECAAN",
                               "--description_field", "FUNCTION",
                               "--output_folder", output_folder,
                               "--log_file", log_file,
                               "--interactive"
                               ])
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "filename")
        with self.subTest():
            self.assertTrue(args.prod_run)
        with self.subTest():
            self.assertEqual(args.run_mode, "pecaan")
        with self.subTest():
            self.assertEqual(args.description_field, "function")
        with self.subTest():
            self.assertEqual(args.output_folder, Path(output_folder))
        with self.subTest():
            self.assertEqual(args.log_file, log_file)
        with self.subTest():
            self.assertTrue(args.interactive)




class TestImportGenomeMain3(unittest.TestCase):

    def setUp(self):
        self.test_directory1 = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_dir")
        os.mkdir(self.test_directory1)
        self.file = Path(self.test_directory1, "new_file.txt")
        self.dir = Path(self.test_directory1, "new_dir")
        self.parser = argparse.ArgumentParser()


    def tearDown(self):
        shutil.rmtree(self.test_directory1)


    def test_set_path_1(self):
        """Verify output when file exists and is expected to exist."""
        self.file.touch()
        output = import_genome.set_path(self.file, kind="file", expect=True)
        with self.subTest():
            self.assertIsInstance(output, Path)
        with self.subTest():
            self.assertEqual(str(self.file), str(output))

    @patch("sys.exit")
    def test_set_path_2(self, sys_exit_mock):
        """Verify script exits when file does not exist,
        and is not expected to exist."""
        output = import_genome.set_path(self.file, kind="file", expect=True)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_3(self, verify_path2_mock):
        """Verify home directory expansion."""
        home = Path("~")
        home = home.expanduser()
        test_file = Path("~/path/to/file.txt")
        verify_path2_mock.return_value = (True, None)
        output = import_genome.set_path(test_file, kind="file", expect=True)
        exp = Path(home, "path/to/file.txt")
        self.assertEqual(output, exp)


    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_4(self, verify_path2_mock):
        """Verify '..' directory resolution."""
        test_file = Path("/dir1/dir2/../file.txt")
        verify_path2_mock.return_value = (True, None)
        output = import_genome.set_path(test_file, kind="file", expect=True)
        exp = Path("/dir1/file.txt")
        self.assertEqual(output, exp)


    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_5(self, verify_path2_mock):
        """Verify home directory expansion and '..' directory resolution."""
        home = Path("~")
        home = home.expanduser()
        test_file = Path("~/dir1/dir2/../file.txt")
        verify_path2_mock.return_value = (True, None)
        output = import_genome.set_path(test_file, kind="file", expect=True)
        exp = Path(home, "dir1/file.txt")
        self.assertEqual(output, exp)




class TestImportGenomeMain4(unittest.TestCase):

    def setUp(self):

        self.import_table = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_import_table_1.csv")
        self.import_table = Path(self.import_table)

        self.base_dir = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_folder")
        self.base_dir = Path(self.base_dir)
        self.base_dir.mkdir()

        self.input_folder = Path(self.base_dir, "input_folder")
        self.output_folder = Path(self.base_dir, "output_folder")
        self.log_file = Path(self.base_dir, "test_log.txt")


        self.sql_handle_1 = mch.MySQLConnectionHandler()
        self.sql_handle_1.database = "Actino_Draft"
        self.sql_handle_1.username = user
        self.sql_handle_1.password = pwd

        self.sql_handle_2 = mch.MySQLConnectionHandler()
        self.sql_handle_2.database = "Actino_Draft"
        self.sql_handle_2.username = user
        self.sql_handle_2.password = pwd
        self.sql_handle_2.credential_status = False
        self.sql_handle_2._database_status = False

        self.args_list = ["run.py",
                          "import",
                          "Actino_Draft",
                          str(self.input_folder),
                          str(self.import_table),
                          "-g", "FILENAME",
                          "-p",
                          "-r", "PECAAN",
                          "-d", "FUNCTION",
                          "-o", str(self.output_folder),
                          "-l", str(self.log_file)
                          ]

    def tearDown(self):
        shutil.rmtree(self.base_dir)

    # import_genome.setup_sql_handle() calls
    # MySQLConnectionHandler.open_connection(),
    # which uses getpass.getpass() to get the username and password.
    @patch("getpass.getpass")
    def test_setup_sql_handle_1(self, getpass_mock):
        """Verify that sql_handle returned with valid info."""
        getpass_mock.side_effect = [user, pwd]
        sql_handle = import_genome.setup_sql_handle("Actino_Draft")
        with self.subTest():
            self.assertTrue(getpass_mock.called)
        with self.subTest():
            self.assertIsNotNone(sql_handle)


    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_2(self, mch_mock, sys_exit_mock):
        """Verify that sys exit is called when credential_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = "Actino_Draft"
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = True
        mch_mock.return_value = sql_handle_mock
        sql_handle = import_genome.setup_sql_handle("Actino_Draft")
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_3(self, mch_mock, sys_exit_mock):
        """Verify that sys exit is called when _database_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = "Actino_Draft"
        sql_handle_mock.credential_status = True
        sql_handle_mock._database_status = False
        mch_mock.return_value = sql_handle_mock
        sql_handle = import_genome.setup_sql_handle("Actino_Draft")
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("sys.exit")
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler")
    def test_setup_sql_handle_4(self, mch_mock, sys_exit_mock):
        """Verify that sys exit is called when credential_status = False and
        _database_status = False."""
        sql_handle_mock = Mock()
        sql_handle_mock.database = "Actino_Draft"
        sql_handle_mock.credential_status = False
        sql_handle_mock._database_status = False
        mch_mock.return_value = sql_handle_mock
        sql_handle = import_genome.setup_sql_handle("Actino_Draft")
        with self.subTest():
            self.assertTrue(mch_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)




    @patch("pdm_utils.pipelines.db_import.import_genome.data_io")
    @patch("pdm_utils.pipelines.db_import.import_genome.setup_sql_handle")
    def test_main_1(self, setup_sql_mock, data_io_mock):
        """Verify that correct args calls data_io."""
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        setup_sql_mock.return_value = self.sql_handle_1
        import_genome.main(self.args_list)
        self.assertTrue(data_io_mock.called)


    @patch("pdm_utils.pipelines.db_import.import_genome.data_io")
    @patch("pdm_utils.pipelines.db_import.import_genome.setup_sql_handle")
    @patch("sys.exit")
    def test_main_2(self, sys_exit_mock, setup_sql_mock, data_io_mock):
        """Verify that invalid input folder calls sys exit."""
        self.output_folder.mkdir()
        setup_sql_mock.return_value = self.sql_handle_1
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.db_import.import_genome.data_io")
    @patch("pdm_utils.pipelines.db_import.import_genome.setup_sql_handle")
    @patch("sys.exit")
    def test_main_3(self, sys_exit_mock, setup_sql_mock, data_io_mock):
        """Verify that invalid import file calls sys exit."""
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        self.args_list[4] = ""
        setup_sql_mock.return_value = self.sql_handle_1
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.db_import.import_genome.data_io")
    @patch("sys.exit")
    @patch("getpass.getpass")
    def test_main_4(self, getpass_mock, sys_exit_mock, data_io_mock):
        """Verify that invalid database calls sys exit."""
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        self.args_list[2] = "Actino_Draft_x"
        getpass_mock.side_effect = [user, pwd]
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)




class TestImportGenomeMain5(unittest.TestCase):

    def setUp(self):

        self.run_mode_eval_dict = {"run_mode": "custom_run_mode",
                                   "eval_flag_dict": {"a":1}}

        self.table_structure_dict = constants.IMPORT_TABLE_STRUCTURE

        self.unittest_file = Path(__file__)
        self.unittest_dir = self.unittest_file.parent
        self.test_import_table_1 = Path(self.unittest_dir,
                                     "test_files/test_import_table_1.csv")

        # Valid data dictionary.
        self.data_dict1 = {}
        self.data_dict1["type"] = "replace"
        self.data_dict1["phage_id"] = "Trixie"
        self.data_dict1["description_field"] = "product"
        self.data_dict1["run_mode"] = "phagesdb"
        self.data_dict1["host_genus"] = "retrieve"
        self.data_dict1["cluster"] = "retain"

        # Valid data dictionary.
        self.data_dict2 = {}
        self.data_dict2["type"] = "replace"
        self.data_dict2["phage_id"] = "L5"
        self.data_dict2["description_field"] = "product"
        self.data_dict2["run_mode"] = "phagesdb"
        self.data_dict2["host_genus"] = "retrieve"
        self.data_dict2["cluster"] = "retain"


        self.dict1_dict2 = [self.data_dict1, self.data_dict2]




    def test_prepare_tickets_1(self):
        """Verify dictionary is returned from a correct import table file."""

        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        run_mode_eval_dict=self.run_mode_eval_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertEqual(len(tkt_dict.keys()), 2)


    # Patch so that a variety of different types of files don't need
    # to be created just to test this function.
    @patch("pdm_utils.functions.tickets.retrieve_ticket_data")
    def test_prepare_tickets_2(self, mock_retrieve_tickets):
        """Verify dictionary is returned from two correct
        import data dictionaries."""

        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        run_mode_eval_dict=self.run_mode_eval_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertEqual(len(tkt_dict.keys()), 2)


    @patch("pdm_utils.functions.tickets.retrieve_ticket_data")
    def test_prepare_tickets_3(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        and one incorrect import data dictionaries."""

        self.data_dict2.pop("phage_id")
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        run_mode_eval_dict=self.run_mode_eval_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)




    @patch("pdm_utils.functions.tickets.retrieve_ticket_data")
    def test_prepare_tickets_4(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        data dictionary and one correct data dictionary with
        duplicated phage_id."""

        self.data_dict2["phage_id"] = "Trixie"
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        run_mode_eval_dict=self.run_mode_eval_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)


    @patch("pdm_utils.functions.tickets.retrieve_ticket_data")
    def test_prepare_tickets_5(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        data dictionary and one incorrect data dictionary with
        invalid ticket type."""

        self.data_dict2["type"] = "invalid"
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        run_mode_eval_dict=self.run_mode_eval_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)





class TestImportGenomeMain6(unittest.TestCase):


    def setUp(self):

        self.flat_file_l5 = os.path.join(os.path.dirname(__file__),
                            "test_files/test_flat_file_1.gb")
        self.flat_file_l5 = Path(self.flat_file_l5)

        self.flat_file_trixie = os.path.join(os.path.dirname(__file__),
                            "test_files/test_flat_file_6.gb")
        self.flat_file_trixie = Path(self.flat_file_trixie)

        self.sql_handle = mch.MySQLConnectionHandler()
        self.sql_handle.database = "Actino_Draft"
        self.sql_handle.username = user
        self.sql_handle.password = pwd

        self.eval_flags = {}
        self.eval_flags["check_seq"] = True
        self.eval_flags["check_id_typo"] = True
        self.eval_flags["check_host_typo"] = True
        self.eval_flags["check_author"] = True
        self.eval_flags["check_trna"] = True
        self.eval_flags["check_gene"] = True
        self.eval_flags["check_locus_tag"] = True
        self.eval_flags["check_description"] = True
        self.eval_flags["check_description_field"] = True
        self.eval_flags["import_locus_tag"] = True

        self.data_dict1 = {}
        self.data_dict1["type"] = "replace"
        self.data_dict1["phage_id"] = "L5"
        self.data_dict1["description_field"] = "product"
        self.data_dict1["run_mode"] = "phagesdb"
        self.data_dict1["host_genus"] = "Mycobacterium"
        self.data_dict1["cluster"] = "A"

        self.tkt1 = ticket.GenomeTicket()
        self.tkt1.id = 1
        self.tkt1.phage_id = "L5"
        self.tkt1.run_mode = "phagesdb"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.data_dict = self.data_dict1

        self.data_dict2 = {}
        self.data_dict2["type"] = "replace"
        self.data_dict2["id"] = 1
        self.data_dict2["phage_id"] = "Trixie"
        self.data_dict2["description_field"] = "product"
        self.data_dict2["run_mode"] = "phagesdb"
        self.data_dict2["host_genus"] = "Gordonia"
        self.data_dict2["cluster"] = "B"

        self.tkt2 = ticket.GenomeTicket()
        self.tkt2.id = 2
        self.tkt2.phage_id = "Trixie"
        self.tkt2.run_mode = "phagesdb"
        self.tkt2.description_field = "product"
        self.tkt2.eval_flags = self.eval_flags
        self.tkt2.data_dict = self.data_dict2




    def test_process_files_and_tickets_1(self):
        """Verify correct output using:
        no files,
        no unmatched tickets."""
        ticket_dict = {}
        files = []
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(evaluation_dict.keys()), 0)


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_into_db")
    def test_process_files_and_tickets_2(self, import_into_db_mock):
        """Verify correct output using:
        two files with no tickets,
        unsuccessful import,
        no unmatched tickets."""
        ticket_dict = {}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [False, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")

        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 2)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_into_db")
    def test_process_files_and_tickets_3(self, import_into_db_mock):
        """Verify correct output using:
        two files with matched tickets,
        successful import,
        no unmatched tickets."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [True, True]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 2)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_into_db")
    def test_process_files_and_tickets_4(self, import_into_db_mock):
        """Verify correct output using:
        two files with matched tickets,
        unsuccessful import,
        no unmatched tickets."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [False, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 2)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))


    def test_process_files_and_tickets_5(self):
        """Verify correct output using:
        no files,
        two unmatched tickets."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = []
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_into_db")
    def test_process_files_and_tickets_6(self, import_into_db_mock):
        """Verify correct output using:
        one file matched to ticket with successful import,
        one file unmatched to ticket with unsuccessful import,
        one unmatched ticket."""
        self.tkt2.phage_id = "Trixie_x"
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [True, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 1)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 1)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 1)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 1)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2, 3]))


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_7(self, ask_mock, execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = False."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]

        # self.tkt2.evaluations = [eval.Eval(status="error"),
        #                          eval.Eval(status="error")]
        # ticket_dict = {self.tkt1.phage_id: self.tkt1,
        #                self.tkt2.phage_id: self.tkt2}
        # files = [self.flat_file_l5, self.flat_file_trixie]

        # ask_mock.side_effect = [True, True]

        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=False)
        evaluation_dict = results_tuple[4]
        errors = get_errors_from_dict(evaluation_dict)
        with self.subTest():
            self.assertFalse(ask_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertTrue(errors > 0)


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_8(self, ask_mock, execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = True,
        with no 'error' status corrections."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ask_mock.return_value = False
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        evaluation_dict = results_tuple[4]
        errors = get_errors_from_dict(evaluation_dict)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertTrue(errors > 0)


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_9(self, ask_mock, execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = True,
        with all 'error' status corrections to 'warning'."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ask_mock.return_value = True
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        evaluation_dict = results_tuple[4]
        errors = get_errors_from_dict(evaluation_dict)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertTrue(execute_mock.called)
        with self.subTest():
            self.assertEqual(errors, 0)


    # Patching so avoid an attempt to add data to the database.
    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    @patch("pdm_utils.pipelines.db_import.import_genome.run_checks")
    @patch("pdm_utils.pipelines.db_import.import_genome.prepare_bundle")
    def test_process_files_and_tickets_10(self, prep_mock, run_checks_mock,
                                          ask_mock, execute_mock):
        """Verify correct output using:
        two files with matched tickets,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = True,
        with all 'error' status corrections to 'warning' in first file
        but no corrections in second file."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        self.tkt2.evaluations = [eval.Eval(status="error")]
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        gnm1 = genome.Genome()
        bndl1 = bundle.Bundle()
        bndl1.id = 1
        bndl1.ticket = self.tkt1
        bndl1.genome_dict["flat_file"] = gnm1

        gnm2 = genome.Genome()
        bndl2 = bundle.Bundle()
        bndl2.id = 2
        bndl2.ticket = self.tkt2
        bndl2.genome_dict["flat_file"] = gnm2

        prep_mock.side_effect = [bndl1, bndl2]
        ask_mock.side_effect = [True, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, sql_handle=self.sql_handle,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertTrue(execute_mock.called)
        with self.subTest():
            self.assertEqual(bndl1._errors, 0)
        with self.subTest():
            self.assertEqual(bndl2._errors, 1)




class TestImportGenomeMain7(unittest.TestCase):

    def setUp(self):

        self.base_dir = os.path.join(os.path.dirname(__file__),
                            "test_wd/test_import")
        self.base_dir = Path(self.base_dir)
        self.base_dir.mkdir()

        self.genome_folder = Path(self.base_dir, "genomes")

        self.flat_file1 = Path(self.genome_folder, "flat_file1.txt")
        self.flat_file2 = Path(self.genome_folder, "flat_file2.txt")

        self.valid_import_table_file = os.path.join(os.path.dirname(__file__),
                            "test_files/test_import_table_1.csv")
        self.valid_import_table_file = Path(self.valid_import_table_file)

        self.invalid_import_table_file = os.path.join(os.path.dirname(__file__),
                            "test_files/test_import_table_2.csv")
        self.invalid_import_table_file = Path(self.invalid_import_table_file)

        self.output_folder = Path(self.base_dir, "output_folder")

        self.sql_handle = mch.MySQLConnectionHandler()
        self.sql_handle.database = db
        self.sql_handle.username = user
        self.sql_handle.password = pwd


        self.tkt1 = ticket.GenomeTicket()
        self.tkt1.phage_id = "L5"
        self.tkt1.host_genus = "Mycobacterium"

        self.tkt2 = ticket.GenomeTicket()
        self.tkt2.phage_id = "Trixie"
        self.tkt2.host_genus = "Mycobacterium"

        self.tkt_dict1 = {"phage_id": "L5", "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"phage_id": "Trixie", "host_genus": "Mycobacterium"}

        self.date = time.strftime("%Y%m%d")
        self.results_folder1 = "{}_results".format(self.date)
        self.results_folder2 = "{}_results_1".format(self.date)
        self.results_folder3 = "{}_results_2".format(self.date)

        self.exp_success = Path(self.output_folder, self.results_folder1, "success")
        self.exp_success_tkt_table = Path(self.exp_success, "import_tickets.csv")
        self.exp_success_genomes = Path(self.exp_success, "genomes")

        self.exp_fail = Path(self.output_folder, self.results_folder1, "fail")
        self.exp_fail_tkt_table = Path(self.exp_fail, "import_tickets.csv")
        self.exp_fail_genomes = Path(self.exp_fail, "genomes")

    def tearDown(self):
        shutil.rmtree(self.base_dir)


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_1(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io runs correctly when there are no errors."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        success_ticket_list = [self.tkt_dict1]
        failed_ticket_list = [self.tkt_dict1, self.tkt_dict2]
        success_filename_list = [self.flat_file1]
        failed_filename_list = [self.flat_file2]
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        exp_success_tkts = []
        with open(self.exp_success_tkt_table,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_success_tkts.append(dict)

        exp_fail_tkts = []
        with open(self.exp_fail_tkt_table,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_fail_tkts.append(dict)

        input_genomes_count = 0
        for item in self.genome_folder.iterdir():
            input_genomes_count += 1

        success_genomes_count = 0
        for item in self.exp_success_genomes.iterdir():
            success_genomes_count += 1
        fail_genomes_count = 0
        for item in self.exp_fail_genomes.iterdir():
            fail_genomes_count += 1

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)
        with self.subTest():
            self.assertEqual(len(exp_success_tkts), 1)
        with self.subTest():
            self.assertEqual(len(exp_success_tkts[0].keys()), 11)
        with self.subTest():
            self.assertEqual(len(exp_fail_tkts), 2)
        with self.subTest():
            self.assertEqual(len(exp_fail_tkts[0].keys()), 11)
        with self.subTest():
            self.assertEqual(input_genomes_count, 0)
        with self.subTest():
            self.assertEqual(success_genomes_count, 1)
        with self.subTest():
            self.assertEqual(fail_genomes_count, 1)


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_2(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is successful with
        success tickets but no success files, and
        fail tickets but no fail files."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        success_ticket_list = [self.tkt_dict1]
        failed_ticket_list = [self.tkt_dict1, self.tkt_dict2]
        success_filename_list = []
        failed_filename_list = []
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)
        with self.subTest():
            self.assertTrue(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_genomes.exists())


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_3(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is successful with
        success files but no success tickets, and
        fail files but no fail tickets."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()
        success_ticket_list = []
        failed_ticket_list = []
        success_filename_list = [self.flat_file1]
        failed_filename_list = [self.flat_file2]
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)
        with self.subTest():
            self.assertFalse(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertTrue(self.exp_fail_genomes.exists())


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_4(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is successful with two previously-existing
        results folders and no items in any output dataset."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        self.results_folder1 = Path(self.output_folder, self.results_folder1)
        self.results_folder1.mkdir()
        self.results_folder2 = Path(self.output_folder, self.results_folder2)
        self.results_folder2.mkdir()
        self.results_folder3 = Path(self.output_folder, self.results_folder3)
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        input_genomes_count = 0
        for item in self.genome_folder.iterdir():
            input_genomes_count += 1

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)
        with self.subTest():
            self.assertTrue(self.results_folder3.exists())
        with self.subTest():
            self.assertFalse(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_genomes.exists())
        with self.subTest():
            self.assertEqual(input_genomes_count, 2)


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_5(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is not successful when results folder is invalid."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        self.results_folder1 = Path(self.output_folder, self.results_folder1)
        self.results_folder1.mkdir()
        self.results_folder2 = Path(self.output_folder, self.results_folder2)
        self.results_folder2.mkdir()
        self.results_folder3 = Path(self.output_folder, self.results_folder3)
        self.results_folder3.mkdir()
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)
        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_6(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is not successful when there are no files to process."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)


    # TODO remove log_evaluations patch if it is no longer called in data_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.log_evaluations")
    @patch("sys.exit")
    @patch("pdm_utils.pipelines.db_import.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.phamerator.get_phage_table_count")
    def test_data_io_7(self, get_count_mock, pft_mock, sys_exit_mock, log_eval_mock):
        """Verify data_io is not successful when there are no tickets to process."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(sql_handle=self.sql_handle,
            genome_folder=self.genome_folder,
            import_table_file=self.invalid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        # with self.subTest():
        #     self.assertTrue(log_eval_mock.called)




class TestImportGenomeMain8(unittest.TestCase):


    def setUp(self):

        self.eval_error1 = eval.Eval(status="error")
        self.eval_error2 = eval.Eval(status="error")
        self.eval_error3 = eval.Eval(status="error")
        self.eval_error4 = eval.Eval(status="error")
        self.eval_error5 = eval.Eval(status="error")
        self.eval_error6 = eval.Eval(status="error")
        self.eval_error7 = eval.Eval(status="error")
        self.eval_error8 = eval.Eval(status="error")
        self.eval_correct1 = eval.Eval(status="correct")
        self.eval_correct2 = eval.Eval(status="correct")
        self.eval_correct3 = eval.Eval(status="correct")

        self.tkt = ticket.GenomeTicket()

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()
        self.cds3 = cds.Cds()
        self.cds4 = cds.Cds()

        self.src1 = source.Source()
        self.src2 = source.Source()
        self.src3 = source.Source()

        self.gnm1 = genome.Genome()
        self.gnm2 = genome.Genome()
        self.gnm3 = genome.Genome()

        self.genome_pair1 = genomepair.GenomePair()
        self.genome_pair2 = genomepair.GenomePair()
        self.genome_pair3 = genomepair.GenomePair()

        self.bndl = bundle.Bundle()




    def test_review_evaluation_list_1(self):
        """Verify no need for user input if there are no 'error' evals."""
        eval_list = [self.eval_correct1, self.eval_correct2, self.eval_correct3]
        exit = import_genome.review_evaluation_list(eval_list)
        with self.subTest():
            self.assertEqual(self.eval_correct1.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_correct2.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_correct3.status, "correct")
        with self.subTest():
            self.assertFalse(exit)

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_2(self, ask_mock):
        """Verify only one call for user input, and
        verify no changes to status."""
        ask_mock.side_effect = [False]
        eval_list = [self.eval_correct1, self.eval_error1, self.eval_correct2]
        exit = import_genome.review_evaluation_list(eval_list)
        with self.subTest():
            self.assertEqual(self.eval_correct1.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_error1.status, "error")
        with self.subTest():
            self.assertEqual(self.eval_correct3.status, "correct")
        with self.subTest():
            self.assertFalse(exit)




    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_4(self, ask_mock):
        """Verify two calls for user input, and
        verify change to status."""
        ask_mock.side_effect = [True, True]
        eval_list = [self.eval_error1, self.eval_correct1, self.eval_error2]
        exit = import_genome.review_evaluation_list(eval_list)
        with self.subTest():
            self.assertEqual(self.eval_error1.status, "warning")
        with self.subTest():
            self.assertEqual(self.eval_correct1.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_error2.status, "warning")
        with self.subTest():
            self.assertFalse(exit)


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_5(self, ask_mock):
        """Verify only two calls for user input, with
        change to first status, no change for second, and
        exit before third."""
        ask_mock.side_effect = [True, None, True]
        eval_list = [self.eval_error1, self.eval_error2, self.eval_error3]
        exit = import_genome.review_evaluation_list(eval_list)
        with self.subTest():
            self.assertEqual(self.eval_error1.status, "warning")
        with self.subTest():
            self.assertEqual(self.eval_error2.status, "error")
        with self.subTest():
            self.assertEqual(self.eval_error3.status, "error")
        with self.subTest():
            self.assertTrue(exit)




    def test_review_evaluations_1(self):
        """Verify results when bundle has:
        no bundle evaluations,
        no ticket,
        no genomes in genome_dict,
        no genome_pairs in genome_pair_dict."""
        import_genome.review_evaluations(self.bndl)


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_2(self, ask_mock):
        """Verify results when bundle has:
        one error evaluation."""
        ask_mock.side_effect = [True]
        self.bndl.evaluations = [self.eval_error1]
        import_genome.review_evaluations(self.bndl)
        self.assertEqual(self.bndl.evaluations[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_3(self, ask_mock):
        """Verify results when bundle has:
        a ticket with one error evaluation."""
        ask_mock.side_effect = [True]
        self.tkt.evaluations = [self.eval_error1]
        self.bndl.ticket = self.tkt
        import_genome.review_evaluations(self.bndl)
        self.assertEqual(self.bndl.ticket.evaluations[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_4(self, ask_mock):
        """Verify results when bundle has:
        three genomes, two with one error evaluation."""
        ask_mock.side_effect = [True, True]
        self.gnm1.evaluations = [self.eval_error1]
        self.gnm2.evaluations = [self.eval_correct1]
        self.gnm3.evaluations = [self.eval_error2]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        self.bndl.genome_dict["gnm2"] = self.gnm2
        self.bndl.genome_dict["gnm3"] = self.gnm3
        import_genome.review_evaluations(self.bndl)
        evl_list1 = self.bndl.genome_dict["gnm1"].evaluations
        evl_list2 = self.bndl.genome_dict["gnm2"].evaluations
        evl_list3 = self.bndl.genome_dict["gnm3"].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_5(self, ask_mock):
        """Verify results when bundle has:
        one genome with three CDS features, two with one error evaluation."""
        ask_mock.side_effect = [True, True]
        self.cds1.evaluations = [self.eval_error1]
        self.cds2.evaluations = [self.eval_correct1]
        self.cds3.evaluations = [self.eval_error2]
        self.gnm1.cds_features = [self.cds1, self.cds2, self.cds3]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        import_genome.review_evaluations(self.bndl)
        evl_list1 = self.bndl.genome_dict["gnm1"].cds_features[0].evaluations
        evl_list2 = self.bndl.genome_dict["gnm1"].cds_features[1].evaluations
        evl_list3 = self.bndl.genome_dict["gnm1"].cds_features[2].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_6(self, ask_mock):
        """Verify results when bundle has:
        one genome with three source features, two with one error evaluation."""
        ask_mock.side_effect = [True, True]
        self.src1.evaluations = [self.eval_error1]
        self.src2.evaluations = [self.eval_correct1]
        self.src3.evaluations = [self.eval_error2]
        self.gnm1.source_features = [self.src1, self.src2, self.src3]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        import_genome.review_evaluations(self.bndl)
        evl_list1 = self.bndl.genome_dict["gnm1"].source_features[0].evaluations
        evl_list2 = self.bndl.genome_dict["gnm1"].source_features[1].evaluations
        evl_list3 = self.bndl.genome_dict["gnm1"].source_features[2].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_7(self, ask_mock):
        """Verify results when bundle has:
        three genome_pairs, two with one error evaluation."""
        ask_mock.side_effect = [True, True]
        self.genome_pair1.evaluations = [self.eval_error1]
        self.genome_pair2.evaluations = [self.eval_correct1]
        self.genome_pair3.evaluations = [self.eval_error2]
        self.bndl.genome_pair_dict["gnm_pr1"] = self.genome_pair1
        self.bndl.genome_pair_dict["gnm_pr2"] = self.genome_pair2
        self.bndl.genome_pair_dict["gnm_pr3"] = self.genome_pair3
        import_genome.review_evaluations(self.bndl)
        evl_list1 = self.bndl.genome_pair_dict["gnm_pr1"].evaluations
        evl_list2 = self.bndl.genome_pair_dict["gnm_pr2"].evaluations
        evl_list3 = self.bndl.genome_pair_dict["gnm_pr3"].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluations_8(self, ask_mock):
        """Verify results when bundle has:
        one ticket with four evaluations, three that are errors,
        three CDS features, each with one error evaluation, and
        one source feature with one error,
        but with user quitting ticket review after first ticket error review
        and quitting after first CDS feature error review."""
        self.tkt.evaluations = [self.eval_error1, self.eval_correct1,
                                self.eval_error2, self.eval_error3]
        self.cds1.evaluations = [self.eval_error4]
        self.cds2.evaluations = [self.eval_error5]
        self.cds3.evaluations = [self.eval_error6]
        self.src1.evaluations = [self.eval_error7]
        self.gnm1.cds_features = [self.cds1, self.cds2, self.cds3]
        self.gnm1.source_features = [self.src1]
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["gnm1"] = self.gnm1
        ask_mock.side_effect = [True, None, True, None, True]
        import_genome.review_evaluations(self.bndl)
        tkt_list = self.bndl.ticket.evaluations
        cds1_list = self.bndl.genome_dict["gnm1"].cds_features[0].evaluations
        cds2_list = self.bndl.genome_dict["gnm1"].cds_features[1].evaluations
        cds3_list = self.bndl.genome_dict["gnm1"].cds_features[2].evaluations
        src1_list = self.bndl.genome_dict["gnm1"].source_features[0].evaluations
        with self.subTest():
            self.assertEqual(tkt_list[0].status, "warning")
        with self.subTest():
            self.assertEqual(tkt_list[1].status, "correct")
        with self.subTest():
            self.assertEqual(tkt_list[2].status, "error")
        with self.subTest():
            self.assertEqual(tkt_list[3].status, "error")
        with self.subTest():
            self.assertEqual(cds1_list[0].status, "warning")
        with self.subTest():
            self.assertEqual(cds2_list[0].status, "error")
        with self.subTest():
            self.assertEqual(cds3_list[0].status, "error")
        with self.subTest():
            self.assertEqual(src1_list[0].status, "warning")







class TestImportGenomeMain9(unittest.TestCase):

    def setUp(self):

        self.cds1 = cds.Cds()
        self.cds1.id = "L5_001"
        self.cds1.left = 10
        self.cds1.right = 20
        self.cds1.strand = "F"
        self.cds1.processed_product = "int"
        self.cds1.processed_function = ""
        self.cds1.processed_note = ""

        self.cds2 = cds.Cds()
        self.cds2.id = "L5_002"
        self.cds2.left = 12345
        self.cds2.right = 12445
        self.cds2.strand = "F"
        self.cds2.processed_product = "capsid"
        self.cds2.processed_function = ""
        self.cds2.processed_note = "random"

        self.cds3 = cds.Cds()
        self.cds3.id = "L5_001"
        self.cds3.left = 150000
        self.cds3.right = 151000
        self.cds3.strand = "R"
        self.cds3.processed_product = "rep"
        self.cds3.processed_function = "lysB"
        self.cds3.processed_note = ""

        self.field = "product"

        self.gnm = genome.Genome()
        self.tkt = ticket.GenomeTicket()


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_1(self, ask_mock):
        """Verify no change if description_field is correct."""
        ask_mock.return_value = True
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_2(self, ask_mock):
        """Verify no change if description_field is not correct,
        and new field is not selected."""
        ask_mock.side_effect = [False, None]
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_3(self, ask_mock):
        """Verify no change if description_field is not correct,
        and no new field is selected."""
        ask_mock.side_effect = [False, False, False]
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_4(self, ask_mock):
        """Verify no change if description_field is not correct,
        and new field is selected."""
        ask_mock.side_effect = [False, False, True]
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertTrue(new_field in {"function", "note"})




    def test_set_cds_descriptions_1(self):
        """Verify descriptions are set from product field,
        with no interactivity."""
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3]
        self.tkt.description_field = "product"
        import_genome.set_cds_descriptions(self.gnm, self.tkt, interactive=False)
        with self.subTest():
            self.assertEqual(self.gnm._cds_processed_descriptions_tally, 3)
        with self.subTest():
            self.assertEqual(self.cds1.processed_description, "int")
        with self.subTest():
            self.assertEqual(self.cds2.processed_description, "capsid")
        with self.subTest():
            self.assertEqual(self.cds3.processed_description, "rep")

    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_set_cds_descriptions_2(self, ask_mock, choose_mock):
        """Verify descriptions are set from product field,
        with interactivity = True and user changes description_field."""
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3]
        self.tkt.description_field = "product"
        ask_mock.side_effect = [False]
        choose_mock.side_effect = ["function"]
        import_genome.set_cds_descriptions(self.gnm, self.tkt, interactive=True)
        with self.subTest():
            self.assertEqual(self.gnm._cds_processed_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(self.cds1.processed_description, "")
        with self.subTest():
            self.assertEqual(self.cds2.processed_description, "")
        with self.subTest():
            self.assertEqual(self.cds3.processed_description, "lysB")
        with self.subTest():
            self.assertEqual(self.tkt.description_field, "function")

if __name__ == '__main__':
    unittest.main()
