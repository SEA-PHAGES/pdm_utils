"""Tests the functionality of the unique functions in the export_db pipeline"""

from argparse import ArgumentError
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pdm_utils.classes import genome, cds
from pdm_utils.functions import flat_files
from pdm_utils.pipelines import export_db
from unittest.mock import call
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock
import unittest
import os

class TestExportMain(unittest.TestCase):
    def setUp(self):
        self.args_mock = Mock()  

        self.pipeline_mock = Mock()
        self.alchemist_mock = Mock() 

        self.output_path_mock = Mock()
        self.output_name_mock = Mock() 
        self.values_mock = Mock()
        self.verbose_mock = Mock() 

        self.table_mock = Mock()
        self.filters_mock = Mock()
        self.groups_mock = Mock()

        type(self.args_mock).pipeline = \
                            PropertyMock(return_value=self.pipeline_mock) 
        type(self.args_mock).alchemist = \
                            PropertyMock(return_value=self.alchemist_mock)

        type(self.args_mock).output_path = \
                            PropertyMock(return_value=self.output_path_mock)
        type(self.args_mock).output_name = \
                            PropertyMock(return_value=self.output_name_mock) 
        type(self.args_mock).values = \
                            PropertyMock(return_value=self.values_mock)
        type(self.args_mock).verbose = \
                            PropertyMock(return_value=self.verbose_mock)

        type(self.args_mock).table = \
                            PropertyMock(return_value=self.table_mock)
        type(self.args_mock).fitlers = \
                            PropertyMock(return_value=self.filters_mock)
        type(self.args_mock).groups = \
                            PropertyMock(return_value=self.groups_mock)

    def test_main_1(self):
        """Verify function structure of main(). 
        """ 
        pass

class TestExecuteExport(unittest.TestCase):
    def setUp(self):
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

    def test_parse_value_list_input(self):
        """
        Unittest for export_db.parse_phage_list_input()
            -Tests for single dispatch handling of list parameter type
             and undefined parameter type
        """
        pass

    def test_set_cds_seqfeatures(self):
        """
        Unittest for export_db.parse_set_cds_seqfeatures()
          -Tests for ability to order cds_features of a given
             Genome object.
            -Tests for None parameter error handling.
        """
        pass

    def test_append_database_version(self):
        """
        Unittest for export_db.test_append_database_version()
            -Tests for handling of valid and invalid inputs
             of SeqRecord objects and database version dictionaries
             for append_database_version()
        """
        pass
        
if __name__ == "__main__":
    unittest.main()
