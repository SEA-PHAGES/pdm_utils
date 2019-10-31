"""Tests the functionality of the unique functions in the file_export pipeline"""

from pdm_utils.classes import genome
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import * 
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.pipelines.db_export import file_export
from pdm_utils.classes import mysqlconnectionhandler
from pathlib import Path
from unittest.mock import patch, Mock, call
import os, sys, unittest

class TestFileExport(unittest.TestCase):

    def setUp(self):
        #Creates valid MySqlConnectionHandler
        mch = mysqlconnectionhandler.MySQLConnectionHandler()
        mch._username = "pdm_anon"
        mch._password = "pdm_anon"
        gnm = genome.Genome()
        gnm.name = "Test"
        self.genome = gnm
        self.sql_handle = mch
        seqrecord = SeqRecord(Seq("ATGC"))
        seqrecord.seq.alphabet = IUPACAmbiguousDNA()
        seqrecord.id = "Test_Accession"
        seqrecord.name = "Test"
        self.test_record = seqrecord

    def test_convert_path(self):
        """
        Unittest of file_export.convert_path()
            -Tests for functionality with valid
             and invalid path inputs
        """
        #Sub test to test functionality with a valid path input
        with self.subTest(input_path="valid"):
            #Creates test directory
            test_path = Path.cwd()
            test_path = test_path.joinpath("test")
            test_path.mkdir()
            #Creates string path to test directory
            input_path = os.getcwd()
            input_path = os.path.join(input_path, "test")
            returned_path = file_export.convert_path(input_path)
            #Asserts that convert_path() functions as expected
            self.assertEqual(returned_path, test_path)
            #Breaks down test_directory
            test_path.rmdir()
        #Sub test to test functionality with an invalid path input 
        with self.subTest(input_path="invalid"):
            #Creates path to invalid directory
            input_path = os.getcwd()
            input_path = os.path.join(input_path, "iNvAlId_dIrEcToRy")
            #Asserts that convert_path() raises an exception
            with self.assertRaises(ValueError):
                file_export.convert_path(input_path)
 
    def test_convert_file_path(self):
        """
        Unittest to test file_export.convert_file_path() helper function
            -Tests for paths to a current or nonexisting file
        """
        #Sub test to test functionality with an existing file path input
        with self.subTest(file_path=True):
            #Creates a path to a file
            file_path = Path.cwd()
            file_path = file_path.joinpath("test_file")
            file_path.touch()
            #Asserts that convert_file_path() correctly converted input path
            input_path = os.path.join(os.getcwd(), "test_file")
            returned_path = file_export.convert_file_path(input_path)
            self.assertEqual(returned_path, file_path)
            #Removes created file
            file_path.unlink()
            self.assertFalse(file_path.exists())
        #Sub test to test functionality with an invalid file path input
        with self.subTest(file_path=False):
            #Creates path to directory
            file_path=Path.cwd()
            file_path = file_path.joinpath("test_file")
            file_path.mkdir()
            #Asserts that convert_file_path() raises an exception
            input_path = os.path.join(os.getcwd(), "test_file")
            with self.assertRaises(ValueError):
                file_export.convert_file_path(input_path)
            #Removes created directory
            file_path.rmdir()
            self.assertFalse(file_path.exists())
            
    def test_convert_dir_path(self):
        """
        Unittest to test file_export.convert_file_path() helper function
            -Tests for paths to a current or nonexisting file
        """
        #Sub test to test functionality with an existing directory path input
        with self.subTest(dir_path=True):
            #Creates path to directory
            dir_path=Path.cwd()
            dir_path = dir_path.joinpath("test_dir")
            dir_path.mkdir()
            #Asserts that convert_file_path() correctly converted input path
            input_path = os.path.join(os.getcwd(), "test_dir")
            returned_path = file_export.convert_dir_path(input_path)
            self.assertEqual(returned_path, dir_path)
            #Removes created directory
            dir_path.rmdir()
            self.assertFalse(dir_path.exists())
        #Sub test to test functionality with an invalid directory path input 
        with self.subTest(dir_path=False):
            #Creates a path to a file
            dir_path = Path.cwd()
            dir_path = dir_path.joinpath("test_dir")
            dir_path.touch()
            #Asserts that convert_file_path() raises exception
            input_path = os.path.join(os.getcwd(), "test_dir")
            with self.assertRaises(ValueError):
                file_export.convert_dir_path(input_path)
            #Removes created file
            dir_path.unlink()
            self.assertFalse(dir_path.exists())

    #Tests undefined singledispatch type
    #Tests list singledispatch type
    #Tests Path singledispatch type
    def test_parse_phage_list_input(self):
        """
        Unittest that tests file_export.parse_phage_list_input()
            -Tests for Path object single dispatch handling
            and csv reader functionality
        """
        #Creates test csv file
        csv_path = Path.cwd()
        csv_path = csv_path.joinpath("test_csv")
        csv_path.touch()
        csv_path.write_text("Test")
        #Asserts that parse_phage_list_input correctly read csv
        test_phage_list = file_export.parse_phage_list_input(csv_path)
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
        Unittest that tests file_export.establish_connection()
            -Tests for valid and invalid mysqlconnectionhandler objects
        """
        #Sub test to test functionality with a valid mysqlconnectionhandler
        with self.subTest(valid_mysqlconnectionhandler=True, input_type="Str"):
            #Patches mysqlconnectionhandler in context to return a valid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.db_export.file_export."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Asserts establish_database_connection functionality
                MCHMock.return_value = self.sql_handle
                test_sql_handle = file_export.establish_database_connection(
                                                            "Actino_Draft")
                self.assertEqual(test_sql_handle, self.sql_handle)
                GetPasswordMock.assert_called_once()
                GetPasswordMock.reset_mock()
        #Sub test to test functionality with an invalid mysqlconnectionhandler
        with self.subTest(valid_mysqlconnectionhandler=False, input_type="Str"):
            #Patches mysqlconnectionhandler in context to return an invalid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.db_export.file_export."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Creates faulty MySQLConnectionHandler object
                mch = mysqlconnectionhandler.MySQLConnectionHandler()
                mch._username = "invalid"
                mch._password = "invalid"
                MCHMock.return_value = mch
                #Asserts establish_database_connection() functionality
                file_export.establish_database_connection("Actino_Draft")
                #self.assertEqual(GetPasswordMock.call_count, 4)
                GetPasswordMock.assert_not_called()
                GetPasswordMock.reset_mock()
        #Sub test to test establish_database_connection raises exception
        with self.subTest(valid_mysqlconnectionhandler=True, input_type=None):
            #Patches mysqlconnectionhandler in context to return an invalid
            #MySQLConnectionHandler object
            with patch("pdm_utils.pipelines.db_export.file_export."
                       "mysqlconnectionhandler.MySQLConnectionHandler") \
                                                        as MCHMock:
                #Creates valid MySQLConnectionHandler object
                MCHMock.return_value = self.sql_handle
                #Asserts establish_database_connection() raises exception
                with self.assertRaises(TypeError):
                    file_export.establish_database_connection(None)
                GetPasswordMock.assert_not_called()
                GetPasswordMock.reset_mock()
 

    @patch("pdm_utils.pipelines.db_export.file_export.print")
    def test_seqrecord_to_file(self, PrintMock):
        """
        Unittest that tests file_export.seqrecord_to_file()
            -Patches print
            -Tests input error handling, file writing,
             and printing functionalities
                -Asserts print statement calls with MagicMock object calls
                -Asserts seqrecord_to_file() raises exceptions on bad
                 base input directories
                -Asserts files are created according to naming conventions
        """
        #Sub test to test the writing functionalities of seqrecord_to_file()
        with self.subTest(seqrecord_list=["test_record"], file_format="gb",
                          export_path=Path.cwd(), 
                          export_directory_name="file_export", verbose=False):
            #Assert seqrecord_to_file() correctly creates directory and file
            file_export.seqrecord_to_file([self.test_record], "gb",
                                          Path.cwd()) 
            test_path = (Path.cwd()).joinpath("file_export")
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
                          export_path=Path.cwd(),
                          export_directory_name="test_folder", verbose=False): 
            #Assert seqrecord_to_file() correctly creates directly and file
            file_export.seqrecord_to_file([self.test_record], "fasta",
                                          Path.cwd(),
                                          export_dir_name="test_folder") 
            test_path = (Path.cwd()).joinpath("test_folder")
            self.assertTrue(test_path.is_dir())
            file_path = test_path.joinpath("Test.fasta")
            self.assertTrue(file_path.is_file())
            PrintMock.assert_not_called()
            #Remove test file and directory and reset MagicMock object 
            PrintMock.reset_mock()
            file_path.unlink()
            test_path.rmdir()
        #Sub test to test error handling of seqrecord_to_file()
        with self.subTest(seqrecord_list=[], file_format="gb", 
                          export_path=Path("~/iNvAlId_dIrEcToRy"),
                          export_directory_name="file_export", verbose=False):
            #Asserts seqrecord_to_file() raises exception
            with self.assertRaises(ValueError):
                file_export.seqrecord_to_file([], "gb", 
                                              Path("~/iNvAlId_dIrEcToRy"))
            PrintMock.assert_called_once()
            #Reset MagicMock object
            PrintMock.reset_mock()
        #Sub test to test verbose option print statements
        with self.subTest(seqrecord_list=["test_record"], file_format="gb",
                          export_path=Path.cwd(),
                          export_dir_name="file_export", verbose=True):
            #Assert seqrecord_to_file() correctly creates directly and file
            file_export.seqrecord_to_file([self.test_record], "gb",
                                          Path.cwd(), verbose=True)
            test_path = (Path.cwd()).joinpath("file_export")
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

    def test_write_csv_log(self):
        """
        Unittest that tests file_export.write_csv_log()
            -Patches print
            -Tests file writing functionality
                -Asserts print statement calls with MagicMock object calls
                -Asserts seqrecord_to_file() raises exceptions on bad
                 base input directories
                -Asserts files are created according to naming conventions
        """
        #Sub test to test the writing functionalities of write_csv_log()
        with self.subTest(export_directory_name="file_export"):
            #Asserts write_csv_log correctly creates directory and file
            file_export.write_csv_log([self.genome], Path.cwd())
            dir_path = (Path.cwd()).joinpath("file_export")
            self.assertTrue(dir_path.is_dir())
            log_path = dir_path.joinpath("log.csv")
            self.assertTrue(log_path.is_file())
            #Asserts write_csv_log correctly recognizes existing log
            file_export.write_csv_log([self.genome], Path.cwd())
            second_log_path = log_path.with_name("log2.csv")
            self.assertTrue(second_log_path.is_file())
            #Removes test files and directory
            second_log_path.unlink()
            log_path.unlink()
            dir_path.rmdir()
        #Sub test to test the writing functionalities of write_csv_log()
        with self.subTest(export_directory_name="test_directory"):
            #Asserts write_csv_log correctly creates directory and file
            file_export.write_csv_log([self.genome], Path.cwd(),
                                      export_dir_name="test_directory")
            dir_path = (Path.cwd()).joinpath("test_directory")
            self.assertTrue(dir_path.is_dir())
            log_path = dir_path.joinpath("log.csv")
            self.assertTrue(log_path.is_file())
            #Asserts write_csv_log correctly recognizes existing log
            file_export.write_csv_log([self.genome], Path.cwd(),
                                      export_dir_name="test_directory")
            second_log_path = log_path.with_name("log2.csv")
            self.assertTrue(second_log_path.is_file())
            #Removes test files and directory
            second_log_path.unlink()
            log_path.unlink()
            dir_path.rmdir()


if __name__ == "__main__":
    unittest.main()
