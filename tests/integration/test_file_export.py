"""Tests the functionality of the unique functions in the file_export pipeline"""

from pdm_utils.classes import genome
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.pipelines.db_export import file_export
from pdm_utils.classes import mysqlconnectionhandler
from pathlib import Path
import os, sys, shutil, unittest

class TestFileExport(unittest.TestCase):

    def setUp(self):
        self.test_sql_handle = mysqlconnectionhandler.MySQLConnectionHandler\
               (username = 'pdm_anon',\
                password = 'pdm_anon',\
                database = 'Actino_Draft')
        self.test_sql_handle._credential_status = True
        self.test_sql_handle._database_status = True
        self.test_phage_name_filter_list = [] 
        self.test_seqrecord_list = []
        self.test_format = "gb"
        self.database_name = "Test_db"

    def test_retrieve_seqrecord_from_database_1(self):
        self.test_phage_name_filter_list = ["Trixie", "Alice", "D29"]
        test_seqrecord_list = file_export.\
                retrieve_seqrecord_from_database\
                (self.test_sql_handle, self.test_phage_name_filter_list)
        self.assertEqual(test_seqrecord_list[0].name, "Alice")
        self.assertEqual(test_seqrecord_list[1].name, "D29")
        self.assertEqual(test_seqrecord_list[2].name, "Trixie")
        self.assertEqual(len(test_seqrecord_list), 3)

    def test_retrieve_seqrecord_from_database_2(self):
        self.test_phage_name_filter_list = ["iNvAlID"]
        self.test_seqrecord_list = file_export.retrieve_seqrecord_from_database(self.test_sql_handle, self.test_phage_name_filter_list)
        self.assertFalse(self.test_seqrecord_list)

    def test_retrieve_database_version_1(self):
        test_version_data = \
                file_export.retrieve_database_version\
                (self.test_sql_handle)

        self.assertEqual(len(test_version_data), 2)

    def test_seqfeature_file_output_1(self):
        file_export.seqfeature_file_output\
                (self.test_seqrecord_list, self.test_format,\
                Path(os.getcwd()))
        self.assertTrue(os.path.exists(os.path.join\
                (os.getcwd(), "file_export")))
        self.assertFalse(os.listdir(os.path.join\
                (os.getcwd(), "file_export")))
        shutil.rmtree(os.path.join(os.getcwd(), "file_export"))

    def test_seqfeature_file_output_2(self):
        file_export.seqfeature_file_output\
                (self.test_seqrecord_list,\
                self.test_format, Path(os.getcwd()),\
                export_dir_name = self.database_name)
        self.assertTrue(os.path.exists(os.path.join\
                (os.getcwd(), "Test_db")))
        self.assertFalse(os.listdir(os.path.join\
                (os.getcwd(), "Test_db")))

if __name__ == "__main__":
    unittest.main()
