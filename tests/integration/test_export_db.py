"""Tests the functionality of the unique functions in the export_db pipeline"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import *
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from pdm_utils.classes import genome, filter
from pdm_utils.pipelines import export_db
import sqlalchemy
from pathlib import Path
from unittest.mock import patch, Mock, call, MagicMock
import os, sys, unittest, shutil, csv, pymysql

class TestFileExport(unittest.TestCase):

    def setUp(self):

        #Creates a test database
        self.connection = pymysql.connect(host="localhost",
                                     user="pdm_anon",
                                     password="pdm_anon",
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = (self.connection).cursor()
        cur.execute("SELECT SCHEMA_NAME FROM "
                    "INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = 'test_db'")
        result = cur.fetchall()
        if len(result) != 0:
            cur.execute("DROP DATABASE test_db")
            (self.connection).commit()

        cur.execute("CREATE DATABASE test_db;")
        (self.connection).commit()
        (self.connection).close

        engine_string = f"mysql+pymysql://pdm_anon:pdm_anon@localhost/test_db"
        self.engine = sqlalchemy.create_engine(engine_string, echo=False)

        #Creates Genome object
        gnm = genome.Genome()
        gnm.id = "TestID"
        gnm.accession = "TestAccession"
        gnm.name = "Test"
        gnm.host_genus = "TestHost"
        gnm.length = "TestLength"
        gnm.date = "TestDate"
        gnm.description = "TestDescription"
        gnm.gc = "TestGC"
        gnm.cluster = "TestCluster"
        gnm.subcluster = "TestSubcluster"
        gnm.annotation_status = "TestStatus"
        gnm.retrieve_record = "TestRecord"
        gnm.annotation_author = "TestAuthor"
        self.genome = gnm
        #Creates SeqRecord object
        seqrecord = SeqRecord(Seq("ATGC"))
        seqrecord.seq.alphabet = IUPACAmbiguousDNA()
        seqrecord.id = "Test_Accession"
        seqrecord.name = "Test"
        self.test_record = seqrecord
        #Creates working test directory
        self.test_cwd = (Path.cwd()).joinpath("DELETE_ME")
        (self.test_cwd).mkdir()

    def test_convert_path(self):
        """
        Unittest of export_db.convert_path()
            -Tests for functionality with valid
             and invalid path inputs
        """
        #Sub test to test functionality with a valid path input
        with self.subTest(input_path="valid"):
            #Creates test directory
            test_path = (self.test_cwd).joinpath("test")
            test_path.mkdir()
            #Creates string path to test directory
            input_path = str(self.test_cwd)
            input_path = os.path.join(input_path, "test")
            returned_path = export_db.convert_path(input_path)
            #Asserts that convert_path() functions as expected
            self.assertEqual(returned_path, test_path)
            #Breaks down test_directory
            test_path.rmdir()
        #Sub test to test functionality with an invalid path input
        with self.subTest(input_path="invalid"):
            #Creates path to invalid directory
            input_path = str(self.test_cwd)
            input_path = os.path.join(input_path, "iNvAlId_dIrEcToRy")
            #Asserts that convert_path() raises an exception
            with self.assertRaises(ValueError):
                export_db.convert_path(input_path)

    def test_convert_file_path(self):
        """
        Unittest to test export_db.convert_file_path() helper function
            -Tests for paths to a current or nonexisting file
        """
        #Sub test to test functionality with an existing file path input
        with self.subTest(file_path=True):
            #Creates a path to a file
            file_path = (self.test_cwd).joinpath("test_file")
            file_path.touch()
            #Asserts that convert_file_path() correctly converted input path
            input_path = os.path.join(str(self.test_cwd), "test_file")
            returned_path = export_db.convert_file_path(input_path)
            self.assertEqual(returned_path, file_path)
            #Removes created file
            file_path.unlink()
            self.assertFalse(file_path.exists())
        #Sub test to test functionality with an invalid file path input
        with self.subTest(file_path=False):
            #Creates path to directory
            file_path = (self.test_cwd).joinpath("test_file")
            file_path.mkdir()
            #Asserts that convert_file_path() raises an exception
            input_path = os.path.join(str(self.test_cwd), "test_file")
            with self.assertRaises(ValueError):
                export_db.convert_file_path(input_path)
            #Removes created directory
            file_path.rmdir()
            self.assertFalse(file_path.exists())

        #Tests undefined singledispatch type
    #Tests list singledispatch type
    #Tests Path singledispatch type
    def test_parse_value_list_input(self):
        """
        Unittest that tests export_db.parse_phage_list_input()
            -Tests for Path object single dispatch handling
            and csv reader functionality
        """
        #Creates test csv file
        csv_path = (self.test_cwd).joinpath("test_csv")
        csv_path.touch()
        csv_path.write_text("Test, NotSeen, NotSeen")
        #Asserts that parse_phage_list_input correctly read csv
        test_value_list = export_db.parse_value_list_input(csv_path)
        self.assertEqual(test_value_list[0], "Test")
        #Removes test csv file
        csv_path.unlink()
        self.assertFalse(csv_path.exists())


    def test_write_csv(self):
        """
        Unittest that tests export_db.write_csv()
            -Patches print
            -Tests file writing functionality
                -Asserts print statement calls with MagicMock object calls
                -Asserts seqrecord_to_file() raises exceptions on bad
                 base input directories
                -Asserts files are created according to naming conventions
        """
        pass

    def test_write_database(self):
        """
        Unittest that tests export_db.write_database()
            -Tests file writing,
                -Asserts files are created according to naming conventions
        """
        #Create folder and file paths
        folder_path = (self.test_cwd).joinpath("export")
        # TODO how the preferred output will be.
        # db_path = folder_path.joinpath("test_db_v1.sql")
        # version_path = db_path.with_name("test_db_v1.version")

        # TODO how the temporary output is formatted.
        db_path = folder_path.joinpath("test_db.sql")
        version_path = db_path.with_name("test_db.version")

        #Assert write_database() correctly created files
        export_db.write_database(self.engine, 1, self.test_cwd)
        self.assertTrue(folder_path.exists())
        self.assertTrue(db_path.exists())
        self.assertTrue(version_path.exists())
        test_version = version_path.read_text()
        self.assertEqual(test_version, "1")
        #Remove test folder and files
        version_path.unlink()
        db_path.unlink()
        folder_path.rmdir()
        self.assertFalse(folder_path.exists())
        self.assertFalse(db_path.exists())
        self.assertFalse(version_path.exists())

    def tearDown(self):
        #Tears down current working test directory
        shutil.rmtree(str(self.test_cwd))
        #Tears down test database
        cur = (self.connection).cursor()
        cur.execute("DROP DATABASE test_db;")
        (self.connection).commit()
        (self.connection).close()

if __name__ == "__main__":
    unittest.main()
