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

        self.mock_input_type = Mock()
        self.mock_output_type = Mock()

        self.mock_folder_path = Mock()
        self.mock_folder_name = Mock()
        self.mock_filters = Mock()
        self.mock_groups = Mock()
        self.mock_verbose = Mock()
        self.mock_config = Mock()

        type(self.mock_args).database = PropertyMock(
                                    return_value=self.mock_database)
        type(self.mock_args).revisions_file = PropertyMock(
                                    return_value=self.mock_revisions_file)
        type(self.mock_args).folder_path = PropertyMock(
                                    return_value=self.mock_folder_path)
        type(self.mock_args).folder_name = PropertyMock(
                                    return_value=self.mock_folder_name)
        type(self.mock_args).input_type = PropertyMock(
                                    return_value=self.mock_input_type)
        type(self.mock_args).output_type = PropertyMock(
                                    return_value=self.mock_output_type)
        type(self.mock_args).filters = PropertyMock(
                                    return_value=self.mock_filters)
        type(self.mock_args).groups = PropertyMock(
                                    return_value=self.mock_groups)
        type(self.mock_args).verbose = PropertyMock(
                                    return_value=self.mock_verbose)

    @patch("pdm_utils.pipelines.revise.configfile.build_complete_config")
    @patch("pdm_utils.pipelines.revise.execute_local_revise")
    @patch("pdm_utils.pipelines.revise.pipelines_basic.build_alchemist")
    @patch("pdm_utils.pipelines.revise.parse_revise")
    def test_main_1(self, parse_revise_mock, build_alchemist_mock,
                                             execute_local_revise_mock,
                                             build_complete_config_mock):
        """Verfiy function structure of main().
        """
        type(self.mock_args).pipeline = PropertyMock(
                                    return_value="local")

        parse_revise_mock.return_value = self.mock_args
        build_alchemist_mock.return_value = self.mock_alchemist
        build_complete_config_mock.return_value = self.mock_config

        revise.main(self.test_args_list)

        parse_revise_mock.assert_called_with(self.test_args_list)
        build_alchemist_mock.assert_called_with(self.mock_database, 
                                                config=self.mock_config)

        execute_local_revise_mock.assert_called_with(
                            self.mock_alchemist, self.mock_revisions_file,
                            self.mock_folder_path, self.mock_folder_name,
                            config=self.mock_config,
                            input_type=self.mock_input_type,
                            output_type=self.mock_output_type,
                            filters=self.mock_filters, groups=self.mock_groups,
                            verbose=self.mock_verbose)

if __name__ == "__main__":
    unittest.main()
