"""Tests the functionality of unique functions in the review pipeline
"""
import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from pdm_utils.pipelines import review
class TestMain(unittest.TestCase):
    def setUp(self):
        self.test_args_list = ["test", "args", "list"]
        self.args = Mock()

        self.mock_alchemist = Mock()
        self.mock_values = Mock()

        self.mock_database = Mock()
        self.mock_folder_path = Mock()
        self.mock_folder_name = Mock()
        self.mock_no_review = Mock()
        self.mock_input = Mock()
        self.mock_filters = Mock()
        self.mock_groups = Mock()
        self.mock_sort = Mock()
        self.mock_gene_report = Mock()
        self.mock_summary_report = Mock()
        self.mock_pham_summary_report = Mock()
        self.mock_verbose = Mock()
        self.mock_config = Mock()
        self.mock_force = Mock()


        type(self.args).database = PropertyMock(
                                        return_value=self.mock_database)
        type(self.args).folder_path = PropertyMock(
                                        return_value=self.mock_folder_path)
        type(self.args).folder_name = PropertyMock(
                                        return_value=self.mock_folder_name)
        type(self.args).no_review = PropertyMock(
                                        return_value=self.mock_no_review)
        type(self.args).input = PropertyMock(
                                        return_value=self.mock_input)
        type(self.args).filters = PropertyMock(
                                        return_value=self.mock_filters)
        type(self.args).groups = PropertyMock(
                                        return_value=self.mock_groups)
        type(self.args).sort = PropertyMock(
                                        return_value=self.mock_sort)
        type(self.args).gene_reports = PropertyMock(
                                        return_value=self.mock_gene_report)
        type(self.args).summary_report = PropertyMock(
                                        return_value=self.mock_summary_report)
        type(self.args).verbose = PropertyMock(
                                        return_value=self.mock_verbose)
        type(self.args).force = PropertyMock(
                                        return_value=self.mock_force)
        type(self.args).pham_summary_reports = PropertyMock(
                                        return_value=\
                                                self.mock_pham_summary_report)
        type(self.args).all_reports = PropertyMock(
                                        return_value=False)
   
    @patch("pdm_utils.pipelines.revise.configfile.build_complete_config")
    @patch("pdm_utils.pipelines.review.execute_review")
    @patch("pdm_utils.pipelines.review.parse_review")
    @patch("pdm_utils.pipelines.review.pipelines_basic.build_alchemist")
    @patch("pdm_utils.pipelines.review.pipelines_basic.parse_value_input")
    def test_main_1(self, parse_value_input_mock, build_alchemist_mock, 
                                        parse_review_mock, execute_review_mock,
                                        build_complete_config_mock):
        """Verify the function structure of main().
        """
        parse_review_mock.return_value = self.args
        build_alchemist_mock.return_value = self.mock_alchemist
        parse_value_input_mock.return_value = self.mock_values
        build_complete_config_mock.return_value = self.mock_config

        review.main(self.test_args_list)     

        parse_review_mock.assert_called_with(self.test_args_list)
        build_alchemist_mock.assert_called_with(self.mock_database, 
                                                config=self.mock_config)
        parse_value_input_mock.assert_called_with(self.mock_input)
        execute_review_mock.assert_called_with(
                            self.mock_alchemist, 
                            folder_path=self.mock_folder_path, 
                            folder_name=self.mock_folder_name, 
                            no_review=self.mock_no_review,
                            force=self.mock_force,
                            values=self.mock_values, filters=self.mock_filters,
                            groups=self.mock_groups, sort=self.mock_sort, 
                            gr_reports=self.mock_gene_report, 
                            psr_reports=self.mock_pham_summary_report,
                            s_report=self.mock_summary_report,
                            verbose=self.mock_verbose)

if __name__ == "__main__":
    unittest.main()
