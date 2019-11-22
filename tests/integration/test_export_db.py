"""Tests the functionality of the unique functions in the export_db pipeline"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import *
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from pdm_utils.classes import genome
from pdm_utils.pipelines import export_db
from pdm_utils.classes import mysqlconnectionhandler

from pathlib import Path
from unittest.mock import patch, Mock, call
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
        #Creates valid MySqlConnectionHandler
        mch = mysqlconnectionhandler.MySQLConnectionHandler()
        mch._username = "pdm_anon"
        mch._password = "pdm_anon"
        mch.database = "test_db"
        self.sql_handle = mch
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

    def test_convert_dir_path(self):
        """
        Unittest to test export_db.convert_file_path() helper function
            -Tests for paths to a current or nonexisting file
        """
        #Sub test to test functionality with an existing directory path input
        with self.subTest(dir_path=True):
            #Creates path to directory
            dir_path = (self.test_cwd).joinpath("test_dir")
            dir_path.mkdir()
            #Asserts that convert_file_path() correctly converted input path
            input_path = os.path.join(str(self.test_cwd), "test_dir")
            returned_path = export_db.convert_dir_path(input_path)
            self.assertEqual(returned_path, dir_path)
            #Removes created directory
            dir_path.rmdir()
            self.assertFalse(dir_path.exists())
        #Sub test to test functionality with an invalid directory path input
        with self.subTest(dir_path=False):
            #Creates a path to a file
            dir_path = (self.test_cwd).joinpath("test_dir")
            dir_path.touch()
            #Asserts that convert_file_path() raises exception
            input_path = os.path.join(str(self.test_cwd), "test_dir")
            with self.assertRaises(ValueError):
                export_db.convert_dir_path(input_path)
            #Removes created file
            dir_path.unlink()
            self.assertFalse(dir_path.exists())

    #Tests undefined singledispatch type
    #Tests list singledispatch type
    #Tests Path singledispatch type
    def test_parse_phage_list_input(self):
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
        test_phage_list = export_db.parse_phage_list_input(csv_path)
        self.assertEqual(test_phage_list[0], "Test")
        #Removes test csv file
        csv_path.unlink()
        self.assertFalse(csv_path.exists())

    #Patches get_credentials()
    #Tests mysqlconnectionhandler return type
    #Tests valid connection
    @patch("pdm_utils.classes.mysqlconnectionhandler."
                                                "MySQLConnectionHandler."
                                                "ask_username_and_password")
    def test_establish_database_connection(self, GetPasswordMock):
        """
        Unittest that tests export_db.establish_connection()
            -Tests for valid and invalid mysqlconnectionhandler objects
        """
        #Sub test to test functionality with a valid mysqlconnectionhandler
        with self.subTest(valid_mysqlconnectionhandler=True, input_type="Str"):
            #Patches mysqlconnectionhandler in context to return a valid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.export_db."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Asserts establish_database_connection functionality
                MCHMock.return_value = self.sql_handle
                test_sql_handle = export_db.establish_database_connection(
                                                            "test_db")
                self.assertEqual(test_sql_handle, self.sql_handle)
                GetPasswordMock.assert_called_once()
                GetPasswordMock.reset_mock()
        #Sub test to test functionality with an invalid mysqlconnectionhandler
        with self.subTest(valid_mysqlconnectionhandler=False, input_type="Str"):
            #Patches mysqlconnectionhandler in context to return an invalid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.export_db."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Creates faulty MySQLConnectionHandler object
                mch = mysqlconnectionhandler.MySQLConnectionHandler()
                mch._username = "invalid"
                mch._password = "invalid"
                MCHMock.return_value = mch
                #Asserts establish_database_connection() functionality
                export_db.establish_database_connection("test_db")
                #self.assertEqual(GetPasswordMock.call_count, 4)
                GetPasswordMock.assert_not_called()
                GetPasswordMock.reset_mock()
        #Sub test to test establish_database_connection raises exception
        with self.subTest(valid_mysqlconnectionhandler=True, input_type=None):
            #Patches mysqlconnectionhandler in context to return an invalid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.export_db."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Creates valid MySQLConnectionHandler object
                MCHMock.return_value = self.sql_handle
                #Asserts establish_database_connection() raises exception
                with self.assertRaises(TypeError):
                    export_db.establish_database_connection(None)
                GetPasswordMock.assert_not_called()
                GetPasswordMock.reset_mock()


    @patch("pdm_utils.pipelines.export_db.print")
    def test_write_seqrecord(self, PrintMock):
        """
        Unittest that tests export_db.write_seqrecord()
            -Patches print
            -Tests input error handling, file writing,
             and printing functionalities
                -Asserts print statement calls with MagicMock object calls
                -Asserts write_seqrecord() raises exceptions on bad
                 base input directories
                -Asserts files are created according to naming conventions
        """
        #Sub test to test the writing functionalities of write_seqrecord()
        with self.subTest(seqrecord_list=["test_record"], file_format="gb",
                          export_path=self.test_cwd,
                          export_directory_name="export_db", verbose=False):
            #Assert write_seqrecord() correctly creates directory and file
            export_db.write_seqrecord([self.test_record], "gb",
                                          self.test_cwd)
            test_path = (self.test_cwd).joinpath("export")
            self.assertTrue(test_path.is_dir())
            file_path = test_path.joinpath("Test.gb")
            self.assertTrue(file_path.is_file())
            PrintMock.assert_not_called()
            #Remove test file and directory and reset MagicMock object
            PrintMock.reset_mock()
            file_path.unlink()
            test_path.rmdir()
        #Sub test to test the functionality of optional directory parameters
        with self.subTest(seqrecord_list=["test_record"], file_format="fasta",
                          export_path=self.test_cwd,
                          export_directory_name="test_folder", verbose=False):
            #Assert write_seqrecord() correctly creates directly and file
            export_db.write_seqrecord([self.test_record], "fasta",
                                          self.test_cwd,
                                          export_dir_name="test_folder")
            test_path = (self.test_cwd).joinpath("test_folder")
            self.assertTrue(test_path.is_dir())
            file_path = test_path.joinpath("Test.fasta")
            self.assertTrue(file_path.is_file())
            PrintMock.assert_not_called()
            #Remove test file and directory and reset MagicMock object
            PrintMock.reset_mock()
            file_path.unlink()
            test_path.rmdir()
        #Sub test to test error handling of write_seqrecord()
        with self.subTest(seqrecord_list=[], file_format="gb",
                          export_path=Path("~/iNvAlId_dIrEcToRy"),
                          export_directory_name="export_db", verbose=False):
            #Asserts write_seqrecord() raises exception
            with self.assertRaises(ValueError):
                export_db.write_seqrecord([], "gb",
                                              Path("~/iNvAlId_dIrEcToRy"))
            PrintMock.assert_called_once()
            #Reset MagicMock object
            PrintMock.reset_mock()
        #Sub test to test verbose option print statements
        with self.subTest(seqrecord_list=["test_record"], file_format="gb",
                          export_path=self.test_cwd,
                          export_dir_name="export", verbose=True):
            #Assert write_seqrecord() correctly creates directly and file
            export_db.write_seqrecord([self.test_record], "gb",
                                          self.test_cwd, verbose=True)
            test_path = (self.test_cwd).joinpath("export")
            self.assertTrue(test_path.is_dir())
            file_path = test_path.joinpath("Test.gb")
            self.assertTrue(file_path.is_file())
            PrintMock.assert_has_calls([call("Resolving export path..."),
                                        call("Resolving current export "
                                             "directory status..."),
                                        call("Writing selected data to "
                                             "files..."),
                                        call("Writing Test")])
            #Remove test file and directory
            file_path.unlink()
            test_path.rmdir()

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
        #Sub test to test the writing functionalities of write_csv()
        with self.subTest(export_directory_name="export",
                          csv_file_name="database"):
            #Asserts write_csv correctly creates directory and file
            export_db.write_csv([self.genome], self.test_cwd)
            dir_path = (self.test_cwd).joinpath("export")
            self.assertTrue(dir_path.is_dir())
            csv_path = dir_path.joinpath("database.csv")
            self.assertTrue(csv_path.is_file())
            #Asserts write_csv correctly recognizes existing log
            export_db.write_csv([self.genome], self.test_cwd)
            second_csv_path = csv_path.with_name("database2.csv")
            self.assertTrue(second_csv_path.is_file())
            #Removes test files and directory
            second_csv_path.unlink()
            csv_path.unlink()
            dir_path.rmdir()
        #Sub test to test the directory and file naming functionality
        #of write_csv()
        with self.subTest(export_directory_name="test_directory",
                          csv_file_name="log"):
            #Asserts write_csv correctly creates directory and file
            export_db.write_csv([self.genome], self.test_cwd,
                                      export_dir_name="test_directory",
                                      csv_name="log")
            dir_path = (self.test_cwd).joinpath("test_directory")
            self.assertTrue(dir_path.is_dir())
            csv_path = dir_path.joinpath("log.csv")
            self.assertTrue(csv_path.is_file())
            #Asserts write_csv correctly recognizes existing log
            export_db.write_csv([self.genome], self.test_cwd,
                                      export_dir_name="test_directory",
                                      csv_name="log")
            second_csv_path = csv_path.with_name("log2.csv")
            self.assertTrue(second_csv_path.is_file())
            #Removes test files and directory
            second_csv_path.unlink()
            csv_path.unlink()
            dir_path.rmdir()
        #Test to test the format of write_csv()
        export_db.write_csv([self.genome], self.test_cwd,
                              export_dir_name="Test",
                              csv_name="Test")
        dir_path = (self.test_cwd).joinpath("Test")
        self.assertTrue(dir_path.is_dir())
        csv_path = dir_path.joinpath("Test.csv")
        self.assertTrue(csv_path.is_file())
        if(csv_path.is_file()):
            csv_contents = []
            #Reads csv file contents
            with open(csv_path, newline="") as test_csv:
                csv_reader = csv.reader(test_csv, delimiter=",", quotechar="|")
                for row in csv_reader:
                    csv_contents.append(row)
            self.assertEqual(csv_contents[0],["PhageID",
                                              "Accession",
                                              "Name",
                                              "HostStrain",
                                              "SequenceLength",
                                              "DateLastModified",
                                              "Notes",
                                              "GC",
                                              "Cluster",
                                              "Subcluster",
                                              "Status",
                                              "RetrieveRecord",
                                              "AnnotationAuthor"] )
            self.assertEqual(csv_contents[1],["TestID",
                                              "TestAccession",
                                              "Test",
                                              "TestHost",
                                              "TestLength",
                                              "TestDate",
                                              "TestDescription",
                                              "TestGC",
                                              "TestCluster",
                                              "TestSubcluster",
                                              "TestStatus",
                                              "TestRecord",
                                              "TestAuthor"])

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
        export_db.write_database(self.sql_handle, 1, self.test_cwd)
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
