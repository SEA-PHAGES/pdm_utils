"""Tests the functionality of the unique functions in the file_export pipeline"""

import unittest, os
from argparse import ArgumentError
from unittest.mock import patch, Mock, PropertyMock, call
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.pipelines.db_export import file_export
from pdm_utils.functions import flat_files
from pathlib import Path

class TestFileExport(unittest.TestCase):

    def setUp(self):
        """
        Creates objects for unit testing of the file_export pipeline
        """
        #Creates test Genome objects
        test_phage1 = genome.Genome()
        test_phage1.name = "TestPhage_1"
        test_phage2 = genome.Genome()
        test_phage2.name = "TestPhage_2"
        test_phage3 = genome.Genome()
        test_phage3.name = "TestPhage_3"
        self.genomes = [ test_phage1, test_phage2, test_phage3 ]
        #Creates a test MySQLConnectionHandler object
        mch = mysqlconnectionhandler.MySQLConnectionHandler()
        mch._username = "pdm_anon"
        mch._password = "pdm_anon"
        mch.database = "Test"
        self.sql_handle = mch
        #Creates test Cds objects
        test_cds1 = cds.Cds()
        test_cds1.left = 1
        test_cds1.right = 2
        test_cds1.coordinate_format = "0_half_open"
        test_cds1.strand = 1
        test_cds2 = cds.Cds()
        test_cds2.left = 2
        test_cds2.right = 3
        test_cds2.coordinate_format = "0_half_open"
        test_cds2.strand = 1
        test_cds3 = cds.Cds()
        test_cds3.left = 3
        test_cds3.right = 4
        test_cds3.coordinate_format = "0_half_open"
        test_cds3.strand = 1
        self.cds_list = [test_cds1, test_cds2, test_cds3]
        #Creates a test SeqRecord object
        seq = Seq("ATGC")
        seqrecord = SeqRecord(seq)
        seqrecord.annotations.update({"comment" : ()})
        self.test_seqrecord = seqrecord
        self.test_version_dictionary = \
                {"version" : "Test", "schema_version": "Test"}

    @patch("pdm_utils.pipelines.db_export.file_export.print")
    @patch("pdm_utils.pipelines.db_export.file_export.Cmd_Export")
    @patch("pdm_utils.pipelines.db_export.file_export.execute_file_export")
    @patch(
    "pdm_utils.pipelines.db_export.file_export.establish_database_connection")
    @patch("pdm_utils.pipelines.db_export.file_export.parse_phage_list_input")
    def test_run_file_export(self, ParsePhageInputMock,   EstablishDBMock, 
                             ExecuteFileExportMock, InteractiveMockClass,
                             PrintMock):
        """
        Unittest for file_export.run_file_export function
            -Patches print()
            -Patches file_export.parse_phage_list_input()
            -Patches file_export.establish_database_connection()
            -Patches file_export.execute_file_export()
            -Patches file_export.Cmd_Export
            -Patches file_export.parse_file_export_args() in context
                -Within sub tests, a MagicMock object has predefined 
                 attributes that test the functionality of run_file_export()
                 to successfully branch on input options, verbose and
                 interactivity options, and export options.
        """
        #Sub test that isolates the import_table input branch
        with patch(
        "pdm_utils.pipelines.db_export.file_export.parse_file_export_args") \
            as ParseArgsMock:
            with self.subTest(import_table=["Test"], interactive=False,
                    verbose=False):
                #Establishes predefined MagicMock attributes
                type(ParseArgsMock.return_value).import_table = \
                                            PropertyMock(return_value=["Test"])
                type(ParseArgsMock.return_value).interactive = \
                                            PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).verbose = \
                                            PropertyMock(return_value=False)
                # 
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                PrintMock.assert_not_called()
                ParseArgsMock.assert_called_once_with("Test")
                ParsePhageInputMock.assert_called_once_with("Test")
                ExecuteFileExportMock.assert_called_once()
                InteractiveMockClass.assert_not_called()
        #Resets MagicMock object calls
        PrintMock.reset_mock()
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        InteractiveMockClass.reset_mock()
        #
        #Sub test that isolates the single_genomes input branch
        with patch(
        "pdm_utils.pipelines.db_export.file_export.parse_file_export_args") \
            as ParseArgsMock:
            with self.subTest(import_table=False, single_genomes=True,
                              interactive=False, verbose=False):
                type(ParseArgsMock.return_value).import_table = \
                                            PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).single_genomes = \
                                            PropertyMock(return_value=True)
                type(ParseArgsMock.return_value).interactive = \
                                            PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).verbose = \
                                            PropertyMock(return_value=False)
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                PrintMock.assert_not_called()
                ParseArgsMock.assert_called_once_with("Test")
                ParsePhageInputMock.assert_called_once()
                ExecuteFileExportMock.assert_called_once()
                InteractiveMockClass.assert_not_called()
        #Resets MagicMock object calls
        PrintMock.reset_mock()
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        InteractiveMockClass.reset_mock()
        #
        #Sub test that isolates the interactivity option
        with patch(
        "pdm_utils.pipelines.db_export.file_export.parse_file_export_args") \
            as ParseArgsMock:
            with self.subTest(import_table=False, single_genomes=False,
                              interactive=True, verbose=False):
                #Establishes predefined MagicMock attributes
                type(ParseArgsMock.return_value).import_table = \
                                            PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).single_genomes = \
                                            PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).interactive = \
                                            PropertyMock(return_value=True)
                type(ParseArgsMock.return_value).verbose = \
                                            PropertyMock(return_value=False)
                #Asserts option selection through MagicMock object calls
                file_export.run_file_export("Test")
                PrintMock.assert_not_called()
                ParseArgsMock.assert_called_once_with("Test")
                ParsePhageInputMock.assert_not_called()
                ExecuteFileExportMock.assert_not_called()
                InteractiveMockClass.assert_called_once()
        #Resets MagicMock object calls
        PrintMock.reset_mock()
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        InteractiveMockClass.reset_mock()                                            
        #
        #Sub test that isolates the verbosity option
        with patch(
        "pdm_utils.pipelines.db_export.file_export.parse_file_export_args") \
            as ParseArgsMock:
            with self.subTest(import_table=False, single_genomes=False,
                                  interactive=False, verbose=True):
                #Establishes predefined MagicMock attributes
                type(ParseArgsMock.return_value).import_table = \
                        PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).single_genomes = \
                        PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).interactive = \
                        PropertyMock(return_value=False)
                type(ParseArgsMock.return_value).database = \
                        PropertyMock(return_value="Test")
                type(ParseArgsMock.return_value).verbose = \
                     PropertyMock(return_valie=True)
                #Asserts option selection through MagicMock object calls
                file_export.run_file_export("Test")
                PrintMock.assert_called_with("Please input credentials")
                ParseArgsMock.assert_called_once_with("Test")
                ParsePhageInputMock.assert_not_called()
                ExecuteFileExportMock.assert_called_once()
                InteractiveMockClass.assert_not_called()

    @patch("pdm_utils.pipelines.db_export.file_export.write_database")
    @patch("pdm_utils.pipelines.db_export.file_export.write_csv")
    @patch("pdm_utils.pipelines.db_export.file_export.write_seqrecord")
    @patch(
    "pdm_utils.pipelines.db_export.file_export.phamerator.parse_genome_data")
    @patch(
    "pdm_utils.pipelines.db_export.file_export.flat_files.genome_to_seqrecord")
    @patch("pdm_utils.pipelines.db_export.file_export.set_cds_seqfeatures")
    @patch("pdm_utils.pipelines.db_export.file_export.append_database_version")
    @patch(
    "pdm_utils.pipelines.db_export.file_export.retrieve_database_version")
    @patch("pdm_utils.pipelines.db_export.file_export.print")
    def test_execute_file_export(self, PrintMock, RetrieveDBVersionMock, 
                                 AppendDBVersionMock, SetCdsMock, 
                                 GenomeToSeqRecordMock, ParseGenomeMock, 
                                 WriteFFileMock, WriteCsvMock, WriteDBMock):
        """ 
        Unittest for file_export.run_file_export
            -Patches print()
            -Patches file_export.retrieve_database_version()
            -Patches file_export.append_database_version()
            -Patches file_export.set_cds_seqfeatures()
            -Patches file_export.write_seqrecord()
            -Patches file_export.write_csv()
            -Patches file_export.write_database()
            -Patches flat_files.genome_to_seqrecord()
            -Patches phamerator.parse_genome_data()
            -Executes file_export.execute_file_export() with test parameters
                -Within sub tests, execute_file_export() is called
                 with varying test parameters that test handling and 
                 helper function calling depending on file format,
                 filter list, directory path, and directory name export inputs
                 while also testing for verbose and csv log optional
                 functionality.

        """    
        #Establishes MagicMock return values
        ParseGenomeMock.return_value = self.genomes
        GenomeToSeqRecordMock.return_value = "Test"
        #
        #Sub test that tests for formatted file export without verbose options
        with self.subTest(verbose=False, phage_filter_list="test_filter_list",
                                             export_path="test_path", 
                                             folder_name="test_folder_name",
                                             ff_export="gb",
                                             csv_export=False,
                                             db_export=False):
            file_export.execute_file_export(self.sql_handle, 
                                            ["test_filter_list"], 
                                            "test_path",
                                            "test_folder_name", 
                                            verbose=False, 
                                            ff_export="gb", csv_export=False,
                                            db_export=False)
            #Asserts helper function calls MagicMock object calls
            PrintMock.assert_not_called()
            RetrieveDBVersionMock.assert_called_once()
            AppendDBVersionMock.assert_called()
            SetCdsMock.assert_has_calls([call(self.genomes[0]),
                                         call(self.genomes[1]),
                                         call(self.genomes[2])])
            ParseGenomeMock.assert_called_once()
            WriteFFileMock.assert_called()
            WriteCsvMock.assert_not_called()
            WriteDBMock.assert_not_called()
        #Resets MagicMock object calls
        PrintMock.reset_mock()
        RetrieveDBVersionMock.reset_mock()
        AppendDBVersionMock.reset_mock()
        SetCdsMock.reset_mock()
        ParseGenomeMock.reset_mock()
        WriteFFileMock.reset_mock()
        WriteCsvMock.reset_mock()
        WriteDBMock.reset_mock()
        #Sub test that tests for csv file export with verbose options
        with self.subTest(verbose=False, phage_filter_list="",
                                             export_path="test_path", 
                                             folder_name="test_folder_name",
                                             ff_export=None,
                                             csv_export=True,
                                             db_export=False):
            file_export.execute_file_export(self.sql_handle, [], 
                                            "test_path",
                                            "test_folder_name", 
                                            verbose=True, 
                                            ff_export=None,
                                            csv_export=True,
                                            db_export=False)
            #Asserts helper function calls through MagicMock object calls
            PrintMock.assert_has_calls(
                [call("Retrieving database version..."),
                 call("Retrieving genomic data from Test..."),
                 call("Writing csv file...")], any_order=True)
            RetrieveDBVersionMock.assert_called_once()
            AppendDBVersionMock.assert_not_called()
            SetCdsMock.assert_has_calls([]) 
            ParseGenomeMock.assert_called_once()
            WriteFFileMock.assert_not_called()
            WriteCsvMock.assert_called_once()
            WriteDBMock.assert_not_called()
        #Resets MagicMock object calls 
        PrintMock.reset_mock()
        RetrieveDBVersionMock.reset_mock()
        AppendDBVersionMock.reset_mock()
        SetCdsMock.reset_mock()
        ParseGenomeMock.reset_mock()
        WriteFFileMock.reset_mock()
        WriteCsvMock.reset_mock()
        WriteDBMock.reset_mock()
        #Sub test that tests for csv file export and database file export
        #with verbose options
        with self.subTest(verbose=True, phage_filter_list="test_filter_list",
                                             export_path="test_path", 
                                             folder_name="test_folder_name",
                                             ffile_export=True,
                                             csv_export=False,
                                             db_export=True):
            file_export.execute_file_export(self.sql_handle, 
                                            ["test_filter_list"], 
                                            "test_path",
                                            "test_folder_name", 
                                            verbose=True, 
                                            ff_export="gb",
                                            csv_export=True,
                                            db_export=True)
            #Asserts helper function calls through MagicMock object calls
            PrintMock.assert_has_calls(
                [call("Retrieving database version..."), 
                 call("Writing SQL database file..."),
                 call("Retrieving genomic data from Test..."), 
                 call("Converting genomic data to SeqRecord format..."),
                 call("Converting TestPhage_1"),
                 call("Converting TestPhage_2"),
                 call("Converting TestPhage_3"),
                 call("Appending database version..."),
                 call("Writing csv file...")], any_order=True)
            RetrieveDBVersionMock.assert_called_once()
            AppendDBVersionMock.assert_called()
            SetCdsMock.assert_has_calls([call(self.genomes[0]),
                                         call(self.genomes[1]),
                                         call(self.genomes[2])])
            ParseGenomeMock.assert_called_once()
            WriteFFileMock.assert_called_once()
            WriteCsvMock.assert_called_once()
            WriteDBMock.assert_called_once()
        #Resets MagicMock object calls
        PrintMock.reset_mock()
        RetrieveDBVersionMock.reset_mock()
        AppendDBVersionMock.reset_mock()
        SetCdsMock.reset_mock()
        ParseGenomeMock.reset_mock()
        WriteFFileMock.reset_mock()
        WriteCsvMock.reset_mock()
        WriteDBMock.reset_mock()
        #
        #Sub test that tests for database file export without verbose options
        with self.subTest(verbose=False, phage_filter_list="test_filter_list",
                                             export_path="test_path", 
                                             folder_name="test_folder_name",
                                             ffile_export="gb",
                                             csv_export=False,
                                             db_export=False):
            file_export.execute_file_export(self.sql_handle, 
                                            ["test_filter_list"], 
                                            "test_path",
                                            "test_folder_name", 
                                            verbose=False, 
                                            ff_export=None,
                                            csv_export=False,
                                            db_export=True)
            #Asserts helper function calls through MagicMock object calls
            PrintMock.assert_not_called()
            RetrieveDBVersionMock.assert_called_once()
            AppendDBVersionMock.assert_not_called()
            SetCdsMock.assert_not_called()
            ParseGenomeMock.assert_not_called()
            WriteFFileMock.assert_not_called()
            WriteCsvMock.assert_not_called()
            WriteDBMock.assert_called_once()

    @patch("pdm_utils.pipelines.db_export.file_export.exit")
    @patch("pdm_utils.pipelines.db_export.file_export.print")
    def test_parse_file_export_args(self, PrintMock, ExitMock):
        """ 
        Unittest for file_export.run_file_export()
            -Tests parse_file_export_args() for argparse attribute creation
             depending on various command line inputs for:
                -database 
                -import_table
                -single_genomes
                -all
                -verbose
                -interactive
                -export_directory
                -folder_name
                -ffile_export
                -csv_export
                -database_export
        """         
        #Sub test that tests database export argument parsing
        with self.subTest(import_table=None, single_genomes=None, all=False, 
                          verbose=False, interactive=False, 
                          export_directory=None, folder_name=None, 
                          ffile_export="gb", 
                          csv_export=False, database_export=False):
            #Creates and parses argument list
            unparsed_args = ["blank", "blank", "Test", "-dbx"]
            args = file_export.parse_file_export_args(unparsed_args)
            #Asserts attribute population
            self.assertEqual(args.database, "Test")
            self.assertEqual(args.import_table, None)
            self.assertEqual(args.single_genomes, None)
            self.assertFalse(args.all)
            self.assertFalse(args.verbose)
            self.assertFalse(args.interactive) 
            self.assertEqual(args.export_path, Path.cwd())
            self.assertEqual(args.folder_name, None)
            self.assertEqual(args.ffile_export, None)
            self.assertFalse(args.csv_export)
            self.assertTrue(args.database_export)
        #Sub test that tests csv export argument parsing
        with patch(
        "pdm_utils.pipelines.db_export.file_export.convert_file_path")\
                                                        as ConvertFilePathMock:
            with self.subTest(file_format=None, 
                              import_table="Test", single_genomes=None, 
                              all=False, verbose=False, interactive=False, 
                              export_directory=None, folder_name=None, 
                              ffile_export=None, csv_export=False,
                              database_export=False):
                #Creates and parses argument list
                unparsed_args = ["blank", "blank",
                                 "Test", "--import_table", "Test",
                                 "--csv_export"]
                args = file_export.parse_file_export_args(unparsed_args)
                #Asserts attribute population
                self.assertEqual(args.database, "Test")
                self.assertEqual(args.single_genomes, None)
                self.assertFalse(args.all)
                self.assertFalse(args.verbose)
                self.assertFalse(args.interactive)
                self.assertEqual(args.export_path, Path.cwd())
                self.assertEqual(args.folder_name, None)
                self.assertEqual(args.ffile_export, None)
                self.assertTrue(args.csv_export)
                self.assertFalse(args.database_export)
                PrintMock.assert_not_called()
                ConvertFilePathMock.assert_called_with("Test")
                PrintMock.reset_mock()
                ConvertFilePathMock.reset_mock()
        #Sub test that tests formatted file export, csv export,
        # and database export argument parsing
        with patch(
        "pdm_utils.pipelines.db_export.file_export.convert_dir_path")\
                                                        as ConvertDirPathMock:
            ConvertDirPathMock.return_value = Path.cwd()
            with self.subTest(import_table=None, single_genomes="Test", 
                              all=False, verbose=True, interactive=True,
                              export_directory=os.getcwd(), 
                              folder_name="Test", ffile_export="fasta",
                              csv_export=True, database_export=True):
                #Creates and parses argument list
                unparsed_args = ["blank", "blank",
                                 "Test", "--single_genomes", "Test",
                                 "--verbose", "--interactive",
                                 "--export_path", os.getcwd(),
                                 "--folder_name", "Test", "--ffile_export",
                                 "fasta", "--csv_export", "--database_export"]
                args = file_export.parse_file_export_args(unparsed_args)
                #Asserts attribute population
                self.assertEqual(args.database, "Test")
                self.assertEqual(args.import_table, None)
                self.assertTrue(args.single_genomes, ["Test"])
                self.assertFalse(args.all)
                self.assertTrue(args.verbose)
                self.assertTrue(args.interactive)
                self.assertEqual(args.folder_name, "Test")
                self.assertEqual(args.export_path, Path.cwd())
                self.assertEqual(args.ffile_export, "fasta")
                self.assertTrue(args.csv_export)
                self.assertTrue(args.database_export)
                PrintMock.assert_not_called()
                ConvertDirPathMock.assert_called_with(os.getcwd())
                ConvertDirPathMock.reset_mock()
                PrintMock.reset_mock()
                #Creates and parses argument list with shorthand arguments
                unparsed_args = ["blank", "blank",
                                 "Test", "-sgin", "Test", "-v", "-i",
                                 "-pth", os.getcwd(),
                                 "-name", "Test", 
                                 "-ffx", "fasta", "-csvx", "-dbx"]
                args = file_export.parse_file_export_args(unparsed_args)
                #Asserts attribute population
                self.assertEqual(args.database, "Test")
                self.assertEqual(args.import_table, None)
                self.assertEqual(args.single_genomes, ["Test"])
                self.assertFalse(args.all)
                self.assertTrue(args.verbose)
                self.assertTrue(args.interactive) 
                self.assertEqual(args.export_path, Path.cwd())
                self.assertEqual(args.folder_name, "Test")
                self.assertEqual(args.ffile_export, "fasta")
                self.assertTrue(args.csv_export)
                self.assertTrue(args.database_export)
                PrintMock.assert_not_called()
                ConvertDirPathMock.assert_called_with(os.getcwd())
                ConvertDirPathMock.reset_mock()
                PrintMock.reset_mock()
        #Sub test that tests for no export option handling
        with self.subTest(import_table=None, single_genomes="Test", all=False,
                          verbose=True, interactive=True,
                          export_directory=os.getcwd(), 
                          folder_name="Test", ffile_export=None,
                          csv_export=False, database_export=False):
            #Creates and parses argument list
            unparsed_args = ["blank", "blank", "Test", "--all"]
            args = file_export.parse_file_export_args(unparsed_args)
            #Asserts attribute population
            self.assertEqual(args.database, "Test")
            self.assertEqual(args.import_table, None)
            self.assertEqual(args.single_genomes, None)
            self.assertTrue(args.all)
            self.assertFalse(args.verbose)
            self.assertFalse(args.interactive)
            self.assertEqual(args.export_path, Path.cwd())
            self.assertEqual(args.folder_name, None)
            self.assertFalse(args.csv_export) 
            self.assertEqual(args.ffile_export, None)
            self.assertFalse(args.database_export)
            PrintMock.assert_called()
            PrintMock.reset_mock()
            #Creates and parses argument list with shorthand arguments
            unparsed_args = ["blank", "blank", "Test", "-a"]
            args = file_export.parse_file_export_args(unparsed_args)
            #Asserts attribute population
            self.assertEqual(args.database, "Test")
            self.assertEqual(args.import_table, None)
            self.assertEqual(args.single_genomes, None)
            self.assertTrue(args.all)
            self.assertFalse(args.verbose)
            self.assertFalse(args.interactive)
            self.assertEqual(args.export_path, Path.cwd())
            self.assertEqual(args.folder_name, None) 
            self.assertEqual(args.ffile_export, None)
            self.assertFalse(args.csv_export)
            self.assertFalse(args.database_export)
            PrintMock.assert_called()

    def parse_phage_list_input(self): 
        """
        Unittest for file_export.parse_phage_list_input()
            -Tests for single dispatch handling of list parameter type
             and undefined parameter type
        """
        #Sub test that tests for single dispatch handling of list type
        with self.subTest(input_type="list"):
            phage_list = parse_phage_list_input(["Test", "Test", "Test"])
            self.assertEqual(phage_list, ["Test", "Test", "Test"])
        #Sub test that tests for single dispatch handling of
        with self.subTest(input_type=None):
            with self.assertRaises(TypeError):
                phage_list = parse_phage_list_input(None)          

    def test_set_cds_seqfeatures(self):
        """
        Unittest for file_export.parse_set_cds_seqfeatures()
            -Tests for ability to order cds_features of a given 
             Genome object.
            -Tests for None parameter error handling.
        """
        #Sub test that tests for correct ordering of an ordered list
        with self.subTest(genome_cds_features="Ordered List"):
            test_genome = self.genomes[0]
            test_genome.cds_features = self.cds_list
            file_export.set_cds_seqfeatures(test_genome)
            self.assertEqual(self.cds_list, test_genome.cds_features)
        #Sub test that tests for correct handling of an unordered list
        with self.subTest(genome_cds_features="Unordered List"):
            test_genome = self.genomes[1]
            test_genome.cds_features = [self.cds_list[1], self.cds_list[0],
                                                          self.cds_list[2]]
            file_export.set_cds_seqfeatures(test_genome)
            self.assertEqual(self.cds_list, test_genome.cds_features)
        #Sub test that tests for error handling of None object
        with self.subTest(genome_cds_features="Faulty Genome Object"):
            test_genome = None
            with self.assertRaises(TypeError):
                file_export.set_cds_seqfeatures(test_genome)

    def test_append_database_version(self):
        """
        Unittest for file_export.test_append_database_version()
            -Tests for handling of valid and invalid inputs
             of SeqRecord objects and database version dictionaries
             for append_database_version()
        """
        #Sub test that tests for appending a formatted database version comment
        with self.subTest(test_seqrecord = "Valid Seqrecord",
                          test_version_dictionary = "Valid Dictionary"):
            file_export.append_database_version(self.test_seqrecord, 
                                                self.test_version_dictionary)
            self.assertEqual((self.test_seqrecord.annotations["comment"])[0],
                             "Database Version: Test; Schema Version: Test")
        #Sub test that tests for handling of an invalid version dictionary
        with self.subTest(test_seqrecord = "Valid Seqrecord",
                          test_version_dictionary = "Invalid Dictionary"):
            with self.assertRaises(ValueError):
                file_export.append_database_version(self.test_seqrecord, {})
        #Sub test that tests for handling of an invalid SeqRecord object
        with self.subTest(test_seqrecord = "Invalid Seqrecord",
                          test_version_dictionary = "Valid Dictionary"):
            with self.assertRaises(TypeError):
                file_export.append_database_version(
                                            None, self.test_version_dictionary)

if __name__ == "__main__":
    unittest.main()
