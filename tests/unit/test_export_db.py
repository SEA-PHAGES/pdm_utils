"""Tests the functionality of the unique functions in the export_db pipeline"""
import os
import unittest
from argparse import ArgumentError
from pathlib import Path
from unittest.mock import call
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from pdm_utils.classes import genome, cds
from pdm_utils.functions import flat_files
from pdm_utils.pipelines import export_db

class TestExportMain(unittest.TestCase):
    def setUp(self):
        self.mock_args = Mock()  

        self.mock_pipeline = "gb"
        self.mock_alchemist = Mock() 
        self.mock_engine = Mock() 
        self.mock_database = Mock()

        type(self.mock_alchemist).engine = \
                            PropertyMock(return_value=self.mock_engine)

        self.mock_folder_path = Mock()
        self.mock_folder_name = Mock() 
        self.mock_input = Mock()
        self.mock_values = Mock()
        self.mock_verbose = Mock() 
        self.mock_config = Mock()

        self.mock_table = Mock()
        self.mock_filters = Mock()
        self.mock_groups = Mock()
        self.mock_sort = Mock()

        self.mock_include_columns = Mock()
        self.mock_exclude_columns = Mock()
        self.mock_sequence_columns = Mock()
        self.mock_raw_bytes = Mock()

        self.mock_concatenate = Mock()

        type(self.mock_args).pipeline = \
                            PropertyMock(return_value=self.mock_pipeline) 
        type(self.mock_args).database = \
                            PropertyMock(return_value=self.mock_database)

        type(self.mock_args).folder_path = \
                            PropertyMock(return_value=self.mock_folder_path)
        type(self.mock_args).folder_name = \
                            PropertyMock(return_value=self.mock_folder_name) 
        type(self.mock_args).input = \
                            PropertyMock(return_value=self.mock_input)
        type(self.mock_args).verbose = \
                            PropertyMock(return_value=self.mock_verbose)

        type(self.mock_args).table = \
                            PropertyMock(return_value=self.mock_table)
        type(self.mock_args).filters = \
                            PropertyMock(return_value=self.mock_filters)
        type(self.mock_args).groups = \
                            PropertyMock(return_value=self.mock_groups)
        type(self.mock_args).sort = \
                            PropertyMock(return_value=self.mock_sort)

        type(self.mock_args).include_columns = \
                        PropertyMock(return_value=self.mock_include_columns)
        type(self.mock_args).exclude_columns = \
                        PropertyMock(return_value=self.mock_exclude_columns)
        type(self.mock_args).sequence_columns = \
                        PropertyMock(return_value=self.mock_sequence_columns)
        type(self.mock_args).raw_bytes = \
                            PropertyMock(return_value=self.mock_raw_bytes)

        type(self.mock_args).concatenate = \
                            PropertyMock(return_value=self.mock_concatenate)
       
    @patch("pdm_utils.pipelines.revise.configfile.build_complete_config")
    @patch("pdm_utils.pipelines.export_db.execute_export")
    @patch("pdm_utils.pipelines.export_db.pipelines_basic.parse_value_input")
    @patch("pdm_utils.pipelines.export_db.mysqldb.check_schema_compatibility")
    @patch("pdm_utils.pipelines.export_db.pipelines_basic.build_alchemist")
    @patch("pdm_utils.pipelines.export_db.parse_export")
    def test_main_1(self, parse_export_mock, build_alchemist_mock,
                                             check_schema_compatibility_mock,
                                             parse_value_input_mock,
                                             execute_export_mock,
                                             build_complete_config_mock):
        """Verify function structure of main(). 
        """ 
        parse_export_mock.return_value = self.mock_args
        parse_value_input_mock.return_value = self.mock_values
        build_alchemist_mock.return_value = self.mock_alchemist
        build_complete_config_mock.return_value = self.mock_config
        
        export_db.main(["cmd", "args"])

        parse_export_mock.assert_called_with(["cmd", "args"])
        build_alchemist_mock.assert_called_with(
                                   self.mock_database, config=self.mock_config)
        check_schema_compatibility_mock.assert_called_with(self.mock_engine,
                                                           "export")
        parse_value_input_mock.assert_called_with(self.mock_input)

        execute_export_mock.assert_called_with(
                                    self.mock_alchemist,
                                    self.mock_folder_path,
                                    self.mock_folder_name,
                                    self.mock_pipeline,
                                    table=self.mock_table,
                                    values=self.mock_values,
                                    filters=self.mock_filters,
                                    groups=self.mock_groups,
                                    sort=self.mock_sort, 
                                    include_columns=self.mock_include_columns,
                                    exclude_columns=self.mock_exclude_columns,
                                    raw_bytes=self.mock_raw_bytes,
                                    sequence_columns=self.mock_sequence_columns,
                                    concatenate=self.mock_concatenate,
                                    verbose=self.mock_verbose)



                
        
if __name__ == "__main__":
    unittest.main()
