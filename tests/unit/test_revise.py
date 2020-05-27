import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from pdm_utils.pipelines import revise

class TestReviseMain(unittest.TestCase):
    def setUp(self):
        self.mock_args = Mock()
        self.test_args_list = ["test", "args", "list"]

        self.mock_alchemist = Mock()
        self.mock_revisions_data_dicts = Mock()

        self.mock_database = Mock()
        self.mock_revisions_file = Mock()

        self.mock_folder_path = Mock()
        self.mock_folder_name = Mock()
        self.mock_filters = Mock()
        self.mock_groups = Mock()
        self.mock_verbose = Mock()

        type(self.mock_args).database = PropertyMock(
                                    return_value=self.mock_database)
        type(self.mock_args).revisions_file = PropertyMock(
                                    return_value=self.mock_revisions_file)
        type(self.mock_args).folder_path = PropertyMock(
                                    return_value=self.mock_folder_path)
        type(self.mock_args).folder_name = PropertyMock(
                                    return_value=self.mock_folder_name)
        type(self.mock_args).filters = PropertyMock(
                                    return_value=self.mock_filters)
        type(self.mock_args).groups = PropertyMock(
                                    return_value=self.mock_groups)
        type(self.mock_args).verbose = PropertyMock(
                                    return_value=self.mock_verbose)

    @patch("pdm_utils.pipelines.revise.execute_revise")
    @patch("pdm_utils.pipelines.revise.basic.retrieve_data_dict")
    @patch("pdm_utils.pipelines.revise.AlchemyHandler")
    @patch("pdm_utils.pipelines.revise.parse_revise")
    def test_main_1(self, parse_revise_mock, alchemyhandler_mock,
                            retrieve_data_dict_mock, execute_revise_mock):
        """Verfiy function structure of main().
        """
        parse_revise_mock.return_value = self.mock_args
        alchemyhandler_mock.return_value = self.mock_alchemist
        retrieve_data_dict_mock.return_value = self.mock_revisions_data_dicts

        revise.main(self.test_args_list)

        parse_revise_mock.assert_called_with(self.test_args_list)
        alchemyhandler_mock.assert_called_with(database=self.mock_database)
        retrieve_data_dict_mock.assert_called_with(self.mock_revisions_file)

        execute_revise_mock.assert_called_with(
                            self.mock_alchemist, self.mock_revisions_data_dicts,
                            self.mock_folder_path, self.mock_folder_name,
                            filters=self.mock_filters, groups=self.mock_groups,
                            verbose=self.mock_verbose)

if __name__ == "__main__":
    unittest.main()
