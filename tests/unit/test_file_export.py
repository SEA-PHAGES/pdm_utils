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
        #Creates test phage name list
        phage_name_list = []
        phage_name_list.append("TestPhage_1")
        phage_name_list.append("TestPhage_2")
        phage_name_list.append("TestPhage_3")
        self.names = phage_name_list
        #Creates test Genome objects
        test_phage1 = genome.Genome()
        test_phage1.name = self.names[0]
        test_phage2 = genome.Genome()
        test_phage2.name = self.names[1]
        test_phage3 = genome.Genome()
        test_phage3.name = self.names[2]
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
    @patch("pdm_utils.pipelines.db_export.file_export.parse_export_select")
    @patch("pdm_utils.pipelines.db_export.file_export.parse_phage_list_input")
    @patch("pdm_utils.pipelines.db_export.file_export.execute_file_export")
    @patch(
    "pdm_utils.pipelines.db_export.file_export.establish_database_connection")
    def test_run_file_export(self, EstablishDBMock, 
                             ExecuteFileExportMock,
                             ParsePhageInputMock,
                             ParseArgsMock):
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
        #Establishes predefined MagicMock attributes
        ParsePhageInputMock.return_value = self.names
        EstablishDBMock.return_value = self.sql_handle
        type(ParseArgsMock.return_value).csv_export = \
                                PropertyMock(return_value=False)
        type(ParseArgsMock.return_value).db_export = \
                                PropertyMock(return_value=False)
        type(ParseArgsMock.return_value).ffile_export = \
                                PropertyMock(return_value=False)
        type(ParseArgsMock.return_value).multi_export = \
                                PropertyMock(return_value=False)
        type(ParseArgsMock.return_value).interactive = \
                                PropertyMock(return_value=False)
        type(ParseArgsMock.return_value).database = \
                                PropertyMock(return_value="TestDatabase")

        with patch("pdm_utils.pipelines.db_export.file_export.parse_csvx") \
                as ParseCsvxMock:
            type(ParseCsvxMock.return_value).list_input = \
                                    PropertyMock(return_value=self.names)
            type(ParseCsvxMock.return_value).export_path = \
                                    PropertyMock(return_value=Path.cwd())
            type(ParseCsvxMock.return_value).folder_name = \
                                    PropertyMock(return_value="TestName")
            type(ParseCsvxMock.return_value).verbose = \
                                    PropertyMock(return_value=False)

            #Sub test that tests the csv export branch
            with self.subTest(export_option="csv_export"):
                #Changes predefined MagicMock attributes
                type(ParseArgsMock.return_value).csv_export = \
                                    PropertyMock(return_value=True)
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                ParseArgsMock.assert_called_once_with("Test")
                ExecuteFileExportMock.assert_called_once_with(
                                        self.sql_handle, Path.cwd(),
                                        "TestName", 
                                        phage_filter_list=self.names,
                                        csv_export=True, verbose=False)
                EstablishDBMock.assert_called_once_with("TestDatabase") 
                ParsePhageInputMock.assert_called_once_with(self.names)
        #Resets MagicMock object calls
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        EstablishDBMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        #Restores MagicMock attributes
        type(ParseArgsMock.return_value).csv_export = \
                                PropertyMock(return_value=False)

        with patch("pdm_utils.pipelines.db_export.file_export.parse_dbx") \
                                        as ParseDbxMock:
            #Sub test that isolates the database export branch
            with self.subTest(export_option="db_export"):
                type(ParseDbxMock.return_value).export_path = \
                                        PropertyMock(return_value=Path.cwd())
                type(ParseDbxMock.return_value).folder_name = \
                                        PropertyMock(return_value="TestName")
                type(ParseDbxMock.return_value).verbose = \
                                        PropertyMock(return_value=False)

                #Changes predefined MagicMock attributes
                type(ParseArgsMock.return_value).db_export = \
                                    PropertyMock(return_value=True)
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                ParseArgsMock.assert_called_with("Test")
                ExecuteFileExportMock.assert_called_once_with(
                                        self.sql_handle, Path.cwd(),
                                        "TestName",
                                        db_export=True, verbose=False)
                EstablishDBMock.assert_called_once_with("TestDatabase") 
                ParsePhageInputMock.assert_not_called()
        #Resets MagicMock object calls
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        EstablishDBMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        #Restores MagicMock attributes
        type(ParseArgsMock.return_value).db_export = \
                                PropertyMock(return_value=False)

        with patch("pdm_utils.pipelines.db_export.file_export.parse_ffx") \
                                    as ParseFfxMock:             
            #Sub test that tests the formatted file export branch
            with self.subTest(export_option="ffile_export"):
                type(ParseFfxMock.return_value).file_format = \
                                        PropertyMock(return_value="Format")
                type(ParseFfxMock.return_value).list_input = \
                                        PropertyMock(return_value=self.names)
                type(ParseFfxMock.return_value).export_path = \
                                        PropertyMock(return_value=Path.cwd())
                type(ParseFfxMock.return_value).folder_name = \
                                        PropertyMock(return_value="TestName")
                type(ParseFfxMock.return_value).verbose = \
                                        PropertyMock(return_value=False)


                #Changes predefined MagicMock attributes
                type(ParseArgsMock.return_value).ffile_export = \
                                    PropertyMock(return_value=True)
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                ParseArgsMock.assert_called_once_with("Test")
                ExecuteFileExportMock.assert_called_once_with(
                                        self.sql_handle, Path.cwd(),
                                        "TestName",
                                        phage_filter_list=self.names,
                                        ffile_export="Format", verbose=False)
                EstablishDBMock.assert_called_once_with("TestDatabase") 
                ParsePhageInputMock.assert_called_once_with(self.names)
        #Resets MagicMock object calls
        ParseArgsMock.reset_mock()
        ParsePhageInputMock.reset_mock()
        EstablishDBMock.reset_mock()
        ExecuteFileExportMock.reset_mock()
        #Restores MagicMock attributes
        type(ParseArgsMock.return_value).ffile_export = \
                                PropertyMock(return_value=False)

        with patch("pdm_utils.pipelines.db_export.file_export.parse_multix") \
                                    as ParseMultixMock:
            #Sub test that tests the multiple pipeline export branch
            with self.subTest(export_option="multi_export"):
                type(ParseMultixMock.return_value).list_input = \
                                        PropertyMock(return_value=self.names)
                type(ParseMultixMock.return_value).export_path = \
                                        PropertyMock(return_value=Path.cwd())
                type(ParseMultixMock.return_value).folder_name = \
                                        PropertyMock(return_value="TestName")
                type(ParseMultixMock.return_value).verbose = \
                                        PropertyMock(return_value=False)
                type(ParseMultixMock.return_value).csv_export = \
                                        PropertyMock(return_value=True)
                type(ParseMultixMock.return_value).db_export = \
                                        PropertyMock(return_value=True)
                type(ParseMultixMock.return_value).ffile_export = \
                                        PropertyMock(return_value="Format")

                #Changes predefined MagicMock attributes
                type(ParseArgsMock.return_value).multi_export = \
                                    PropertyMock(return_value=True)
                #Asserts branching through MagicMock object calls
                file_export.run_file_export("Test")
                ParseArgsMock.assert_called_once_with("Test")
                ExecuteFileExportMock.assert_called_once_with(
                                        self.sql_handle, Path.cwd(),
                                        "TestName", 
                                        phage_filter_list=self.names,
                                        csv_export=True, 
                                        db_export=True, ffile_export="Format",
                                        verbose=False) 
                EstablishDBMock.assert_called_once_with("TestDatabase") 
                ParsePhageInputMock.assert_called_once_with(self.names) 
        
    def test_parse_export_select(self):
        """
        Unittest that tests functionality of a created argparse object
        for export selection opttions.
        """
        unparsed_args = ["blank", "blank", "database", "-csvx"]
        args = file_export.parse_export_select(unparsed_args)
        self.assertEqual(args.database, "database")
        self.assertTrue(args.csv_export)
        self.assertFalse(args.db_export)
        self.assertFalse(args.ffile_export)
        self.assertFalse(args.multi_export)
        self.assertFalse(args.interactive)
       
        unparsed_args = ["blank", "blank", "database", "-dbx"]
        args = file_export.parse_export_select(unparsed_args)
        self.assertEqual(args.database, "database") 
        self.assertFalse(args.csv_export)
        self.assertTrue(args.db_export)
        self.assertFalse(args.ffile_export)
        self.assertFalse(args.multi_export)
        self.assertFalse(args.interactive)

        unparsed_args = ["blank", "blank", "database", "-ffx"]
        args = file_export.parse_export_select(unparsed_args)
        self.assertEqual(args.database, "database")
        self.assertFalse(args.csv_export)
        self.assertFalse(args.db_export)
        self.assertTrue(args.ffile_export)
        self.assertFalse(args.multi_export)
        self.assertFalse(args.interactive)

        unparsed_args = ["blank", "blank", "database", "-multix"]
        args = file_export.parse_export_select(unparsed_args)
        self.assertEqual(args.database, "database")
        self.assertFalse(args.csv_export)
        self.assertFalse(args.db_export)
        self.assertFalse(args.ffile_export)
        self.assertTrue(args.multi_export)
        self.assertFalse(args.interactive)

    @patch("pdm_utils.pipelines.db_export.file_export.convert_file_path")
    @patch("pdm_utils.pipelines.db_export.file_export.convert_dir_path")
    def test_parse_multix(self, ConvertDirMock, ConvertFileMock):
        """
        Unittest that tests functionality of a created argparse object
        for multiple export.
        """
    
        unparsed_args = ["blank", "blank", "blank", "blank", "-tin", os.getcwd(), "-v", "-ffx", "gb", "-csvx", "-dbx"] 
        args = file_export.parse_multix(unparsed_args)
        self.assertTrue(args.verbose)
        self.assertEqual(args.ffile_export, "gb")
        self.assertTrue(args.csv_export)
        self.assertTrue(args.db_export)
        ConvertFileMock.assert_called_once()
        ConvertDirMock.assert_not_called()
        ConvertFileMock.reset_mock()
        ConvertDirMock.reset_mock()

        unparsed_args = ["blank", "blank", "blank", "blank", "-sgin", "Test", "-pth", os.getcwd(), "-name", "Test"]
        args = file_export.parse_multix(unparsed_args)
        self.assertEqual(args.list_input, ["Test"])
        self.assertFalse(args.verbose)
        self.assertFalse(args.csv_export)
        self.assertFalse(args.db_export)
        self.assertEqual(args.folder_name, "Test")
        self.assertFalse(args.ffile_export)
        ConvertDirMock.assert_called_once()
        ConvertFileMock.assert_not_called()

    @patch("pdm_utils.pipelines.db_export.file_export.convert_file_path")
    @patch("pdm_utils.pipelines.db_export.file_export.convert_dir_path")
    def test_parse_csvx(self, ConvertDirMock, ConvertFileMock):
        """
        Unittest that tests functionality of a created argparse object
        for csv export.
        """
        unparsed_args = ["blank", "blank", "blank", "blank", "-tin", os.getcwd(), "-v"] 
        args = file_export.parse_csvx(unparsed_args)
        self.assertTrue(args.verbose)
        ConvertFileMock.assert_called_once()
        ConvertDirMock.assert_not_called()
        ConvertFileMock.reset_mock()
        ConvertDirMock.reset_mock()

        unparsed_args = ["blank", "blank", "blank", "blank", "-sgin", "Test", "-pth", os.getcwd(), "-name", "Test"]
        args = file_export.parse_csvx(unparsed_args)
        self.assertFalse(args.verbose)
        self.assertEqual(args.list_input, ["Test"])
        ConvertDirMock.assert_called_once()
        ConvertFileMock.assert_not_called()

    @patch("pdm_utils.pipelines.db_export.file_export.convert_file_path")
    @patch("pdm_utils.pipelines.db_export.file_export.convert_dir_path")
    def test_parse_dbx(self, ConvertDirMock, ConvertFileMock):
        """
        Unittest that tests functionality of a created argparse object
        for sql file export.
        """
        unparsed_args = ["blank", "blank", "blank", "blank", "-v"] 
        args = file_export.parse_dbx(unparsed_args)
        self.assertTrue(args.verbose)
        ConvertFileMock.assert_not_called()
        ConvertDirMock.assert_not_called()
        ConvertFileMock.reset_mock()
        ConvertDirMock.reset_mock()

        unparsed_args = ["blank", "blank", "blank", "blank", "-pth", os.getcwd(), "-name", "Test"]
        args = file_export.parse_dbx(unparsed_args)
        self.assertFalse(args.verbose)
        self.assertEqual(args.folder_name, "Test")
        ConvertDirMock.assert_called_once()
        ConvertFileMock.assert_not_called()
    
    @patch("pdm_utils.pipelines.db_export.file_export.convert_file_path")
    @patch("pdm_utils.pipelines.db_export.file_export.convert_dir_path")
    def test_parse_ffx(self, ConvertDirMock, ConvertFileMock):
        """
        Unittest that tests functionality of a created argparse object
        for formatted file export.
        """
        unparsed_args = ["blank", "blank", "blank", "blank", "gb", "-tin", os.getcwd(), "-v"] 
        args = file_export.parse_ffx(unparsed_args)
        self.assertTrue(args.verbose)
        self.assertEqual(args.file_format, "gb")
        ConvertFileMock.assert_called_once()
        ConvertDirMock.assert_not_called()
        ConvertFileMock.reset_mock()
        ConvertDirMock.reset_mock()

        unparsed_args = ["blank", "blank", "blank", "blank", "fasta", "-sgin", "Test", "-pth", os.getcwd(), "-name", "Test"]
        args = file_export.parse_ffx(unparsed_args)
        self.assertEqual(args.file_format, "fasta")
        self.assertEqual(args.list_input, ["Test"])
        self.assertFalse(args.verbose)
        self.assertEqual(args.folder_name ,"Test")
        ConvertDirMock.assert_called_once()
        ConvertFileMock.assert_not_called()

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
        with self.subTest(verbose=False, export_path="test_path", 
                          folder_name="test_folder_name",
                          phage_filter_list="test_filter_list",
                          ffile_export="gb", csv_export=False,
                          db_export=False):
            file_export.execute_file_export(self.sql_handle, 
                                            "test_path",
                                            "test_folder_name", 
                                            phage_filter_list=\
                                                    ["test_filter_list"],
                                            verbose=False, 
                                            ffile_export="gb", csv_export=False,
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
                                             ffile_export=None,
                                             csv_export=True,
                                             db_export=False):
            file_export.execute_file_export(self.sql_handle, 
                                            "test_path",
                                            "test_folder_name", 
                                            phage_filter_list=[],
                                            verbose=True, 
                                            ffile_export=None,
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
                                            "test_path",
                                            "test_folder_name", 
                                            phage_filter_list=\
                                                    ["test_filter_list"], 
                                            verbose=True, 
                                            ffile_export="gb",
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
            ParseGenomeMock.assert_called()
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
                                            "test_path",
                                            "test_folder_name", 
                                            phage_filter_list=\
                                                    ["test_filter_list"], 
                                            verbose=False, 
                                            ffile_export=None,
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
