"""Tests the functionality of the unique functions in the export_db pipeline"""

import unittest, os
from argparse import ArgumentError
from unittest.mock import patch, Mock, PropertyMock, call
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.pipelines import export_db
from pdm_utils.functions import flat_files
from pathlib import Path

class TestFileExport(unittest.TestCase):

    def setUp(self):
        """
        Creates objects for unit testing of the export_db pipeline
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
                {"Version" : "Test", "SchemaVersion": "Test"}

    @patch("pdm_utils.pipelines.export_db.execute_file_export")
    @patch("pdm_utils.pipelines.export_db.parse_filters")
    @patch("pdm_utils.pipelines.export_db.parse_phage_list_input")
    @patch(
    "pdm_utils.pipelines.export_db.establish_database_connection")
    @patch("pdm_utils.pipelines.export_db.print")
    @patch("pdm_utils.pipelines.export_db.parse_file_export")
    def test_run_file_export(self, ArgParseMock, PrintMock, EstablishDBMock,
                             ParsePhageListMock, ParseFiltersMock,
                             ExecuteExportMock):
        """
        Unittest for export_db.run_file_export()
            -Tests branching based on returned
            argparse attributes
        """
        pass

    def test_parse_file_export(self):
        """
        Unittest for export_db.parse_file_export()
        """
        pass

    def test_execute_file_export(self):
        """
        Unittest for export_db.execute_file_export()
        """
        pass

    def test_execute_ffx_export(self):
        """
        Unittest for export_db.execute_ffx_export()
        """
        pass

    def test_execute_csvx_export(self):
        """
        Unittest for export_db.execute_csvx_export()
        """
        pass

    def test_csvx_grouping(self):
        """
        Unittest for export_db.csvx_grouping()
        """
        pass

    def test_ffx_grouping(self):
        """
        Unittest for export_db.ffx_grouping()
        """
        pass

    def parse_filters(self):
        """
        Unittest for export_db.parse_filters()
        """
        pass

    def test_parse_phage_list_input(self):
        """
        Unittest for export_db.parse_phage_list_input()
            -Tests for single dispatch handling of list parameter type
             and undefined parameter type
        """
        #Sub test that tests for single dispatch handling of list type
        with self.subTest(input_type="list"):
            phage_list = export_db.parse_phage_list_input(
                                         ["Test", "Test", "Test"])
            self.assertEqual(phage_list, ["Test", "Test", "Test"])
        #Sub test that tests for single dispatch handling of
        with self.subTest(input_type=None):
            with self.assertRaises(TypeError):
                phage_list = export_db.parse_phage_list_input(None)

    def test_set_cds_seqfeatures(self):
        """
        Unittest for export_db.parse_set_cds_seqfeatures()
            -Tests for ability to order cds_features of a given
             Genome object.
            -Tests for None parameter error handling.
        """
        #Sub test that tests for correct ordering of an ordered list
        with self.subTest(genome_cds_features="Ordered List"):
            test_genome = self.genomes[0]
            test_genome.cds_features = self.cds_list
            export_db.set_cds_seqfeatures(test_genome)
            self.assertEqual(self.cds_list, test_genome.cds_features)
        #Sub test that tests for correct handling of an unordered list
        with self.subTest(genome_cds_features="Unordered List"):
            test_genome = self.genomes[1]
            test_genome.cds_features = [self.cds_list[1], self.cds_list[0],
                                                          self.cds_list[2]]
            export_db.set_cds_seqfeatures(test_genome)
            self.assertEqual(self.cds_list, test_genome.cds_features)
        #Sub test that tests for error handling of None object
        with self.subTest(genome_cds_features="Faulty Genome Object"):
            test_genome = None
            with self.assertRaises(TypeError):
                export_db.set_cds_seqfeatures(test_genome)

    def test_append_database_version(self):
        """
        Unittest for export_db.test_append_database_version()
            -Tests for handling of valid and invalid inputs
             of SeqRecord objects and database version dictionaries
             for append_database_version()
        """
        #Sub test that tests for appending a formatted database version comment
        with self.subTest(test_seqrecord = "Valid Seqrecord",
                          test_version_dictionary = "Valid Dictionary"):
            export_db.append_database_version(self.test_seqrecord,
                                                self.test_version_dictionary)
            self.assertEqual((self.test_seqrecord.annotations["comment"])[0],
                             "Database Version: Test; Schema Version: Test")
        #Sub test that tests for handling of an invalid version dictionary
        with self.subTest(test_seqrecord = "Valid Seqrecord",
                          test_version_dictionary = "Invalid Dictionary"):
            with self.assertRaises(ValueError):
                export_db.append_database_version(self.test_seqrecord, {})
        #Sub test that tests for handling of an invalid SeqRecord object
        with self.subTest(test_seqrecord = "Invalid Seqrecord",
                          test_version_dictionary = "Valid Dictionary"):
            with self.assertRaises(TypeError):
                export_db.append_database_version(
                                            None, self.test_version_dictionary)

if __name__ == "__main__":
    unittest.main()
