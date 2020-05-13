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
        self.mock_review = Mock()
        self.mock_input = Mock()
        self.mock_filters = Mock()
        self.mock_groups = Mock()
        self.mock_sort = Mock()
        self.mock_pham_gene_report = Mock()
        self.mock_verbose = Mock()


        type(self.args).database = PropertyMock(
                                        return_value=self.mock_database)
        type(self.args).folder_path = PropertyMock(
                                        return_value=self.mock_folder_path)
        type(self.args).folder_name = PropertyMock(
                                        return_value=self.mock_folder_name)
        type(self.args).review = PropertyMock(
                                        return_value=self.mock_review)
        type(self.args).input = PropertyMock(
                                        return_value=self.mock_input)
        type(self.args).filters = PropertyMock(
                                        return_value=self.mock_filters)
        type(self.args).groups = PropertyMock(
                                        return_value=self.mock_groups)
        type(self.args).sort = PropertyMock(
                                        return_value=self.mock_sort)
        type(self.args).pham_gene_report = PropertyMock(
                                        return_value=self.mock_pham_gene_report)
        type(self.args).verbose = PropertyMock(
                                        return_value=self.mock_verbose)

    @patch("pdm_utils.pipelines.review.execute_review")
    @patch("pdm_utils.pipelines.review.parse_review")
    @patch("pdm_utils.pipelines.review.AlchemyHandler")
    @patch("pdm_utils.pipelines.review.export_db.parse_value_input")
    def test_main_1(self, parse_value_input_mock, alchemyhandler_mock, 
                                        parse_review_mock, execute_review_mock):
        """Verify the function structure of main().
        """
        parse_review_mock.return_value = self.args
        alchemyhandler_mock.return_value = self.mock_alchemist
        parse_value_input_mock.return_value = self.mock_values

        review.main(self.test_args_list)     

        parse_review_mock.assert_called_with(self.test_args_list)
        alchemyhandler_mock.assert_called_with(database=self.mock_database)
        parse_value_input_mock.assert_called_with(self.mock_input)
        execute_review_mock.assert_called_with(
                            self.mock_alchemist, self.mock_folder_path, 
                            self.mock_folder_name, review=self.mock_review,
                            values=self.mock_values, filters=self.mock_filters,
                            groups=self.mock_groups, sort=self.mock_sort, 
                            pg_report=self.mock_pham_gene_report, 
                            verbose=self.mock_verbose)

if __name__ == "__main__":
    unittest.main()
