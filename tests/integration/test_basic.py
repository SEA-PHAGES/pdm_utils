"""Integration tests for general functions."""


from pdm_utils.functions import basic
from datetime import datetime
import unittest
from unittest.mock import patch
import os
import shutil
from pathlib import Path




class TestBasicFunctions1(unittest.TestCase):

    def setUp(self):
        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")

        self.test_directory1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_directory1")




    def test_close_files_1(self):
        """Verify one handle is closed."""
        handle1 = open(self.test_filepath1, "r")
        handle_list = [handle1]
        basic.close_files(handle_list)
        self.assertTrue(handle1.closed)

    def test_close_files_2(self):
        """Verify multiple handles are closed."""
        handle1 = open(self.test_filepath1, "r")
        handle2 = open(self.test_filepath1, "r")

        handle_list = [handle1, handle2]
        basic.close_files(handle_list)
        with self.subTest():
            self.assertTrue(handle1.closed)
        with self.subTest():
            self.assertTrue(handle2.closed)




    # Using the patch decorator to bypass the user-required input() value.
    @patch("pdm_utils.functions.basic.get_input", return_value = "yes")
    def test_get_input(self, input):
        """Verify the function retrieves input."""
        value = basic.get_input()
        self.assertEqual(value, "yes")




    @patch("pdm_utils.functions.basic.get_input", return_value = "YES")
    def test_ask_yes_no_1(self, input):
        """Verify user-supplied 'YES' results in True."""
        response = basic.ask_yes_no()
        self.assertTrue(response)

    @patch("pdm_utils.functions.basic.get_input", return_value = "NO")
    def test_ask_yes_no_2(self, input):
        """Verify user-supplied 'NO' results in False."""
        response = basic.ask_yes_no()
        self.assertFalse(response)

    @patch("pdm_utils.functions.basic.get_input", return_value = "invalid")
    def test_ask_yes_no_3(self, input):
        """Verify invalid user-supplied response results in None response."""
        response = basic.ask_yes_no()
        self.assertIsNone(response)












    def test_expand_path_1(self):
        """Verify that home path is expanded."""
        partial_path = "/fake/path"
        input_path = "~" + partial_path
        expanded_path = basic.expand_path(input_path)
        home_dir = os.path.expanduser("~")
        expected_path = home_dir + partial_path
        self.assertEqual(expanded_path, expected_path)

    def test_expand_path_2(self):
        """Verify no change when an absoluate path is provided."""
        input_path = "/fake/path"
        expanded_path = basic.expand_path(input_path)
        home_dir = os.path.expanduser("~")
        expected_path = input_path
        self.assertEqual(expanded_path, expected_path)

    def test_expand_path_3(self):
        """Verify absolute path is expanded when a relative path is provided."""
        partial_path = "/fake/path"
        input_path = "." + partial_path
        expanded_path = basic.expand_path(input_path)
        local_path = os.path.abspath(".")
        expected_path = local_path + partial_path
        self.assertEqual(expanded_path, expected_path)




    def test_verify_path_1(self):
        """Verify a true file produces no error when 'file' is indicated."""
        result = basic.verify_path(self.test_filepath1, "file")
        self.assertTrue(result)

    def test_verify_path_2(self):
        """Verify a fake file produces an error when 'file' is indicated."""
        result = basic.verify_path(self.test_filepath1 + "abcxyz", "file")
        self.assertFalse(result)

    def test_verify_path_3(self):
        """Verify a true directory produces no error when
        'directory' is indicated."""
        result = basic.verify_path(self.test_directory1, "dir")
        self.assertTrue(result)

    def test_verify_path_4(self):
        """Verify a fake directory produces an error when
        'directory' is indicated."""
        result = basic.verify_path(self.test_directory1 + "abcxyz", "dir")
        self.assertFalse(result)

    def test_verify_path_5(self):
        """Verify a true directory produces no error when
        'None' is indicated."""
        result = basic.verify_path(self.test_directory1)
        self.assertTrue(result)

    def test_verify_path_6(self):
        """Verify a fake directory produces an error when
        'None' is indicated."""
        result = basic.verify_path(self.test_directory1 + "abcxyz")
        self.assertFalse(result)

    def test_verify_path_7(self):
        """Verify a true directory produces an error when
        'invalid' is indicated."""
        result = basic.verify_path(self.test_directory1, "invalid")
        self.assertFalse(result)




    def test_parse_flag_file_1(self):
        """Verify that all three rows are parsed correctly."""
        flag_file = os.path.join(os.path.dirname(__file__),
                                 "test_files/test_flag_file_1.csv")
        flag_dict = basic.parse_flag_file(flag_file)
        with self.subTest():
            self.assertEqual(len(flag_dict.keys()), 3)
        with self.subTest():
            self.assertTrue(flag_dict["check_a"])
        with self.subTest():
            self.assertFalse(flag_dict["check_b"])
        with self.subTest():
            self.assertFalse(flag_dict["check_c"])

    def test_parse_flag_file_2(self):
        """Verify that only one row is parsed correctly."""
        flag_file = os.path.join(os.path.dirname(__file__),
                                 "test_files/test_flag_file_2.csv")
        flag_dict = basic.parse_flag_file(flag_file)
        self.assertEqual(len(flag_dict.keys()), 1)



    # Using the patch decorator to bypass the user-required input() value.
    @patch("getpass.getpass")
    def test_get_user_pwd_1(self, getpass_mock):
        """."""
        mock_results = ["user", "pwd"]
        getpass_mock.side_effect = mock_results
        results = basic.get_user_pwd()
        with self.subTest():
            self.assertEqual(results[0], "user")
        with self.subTest():
            self.assertEqual(results[1], "pwd")



class TestBasicFunctions2(unittest.TestCase):

    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/new_dir")
        os.mkdir(self.base_dir)




    def test_make_new_dir_1(self):
        """Verify new directory is created."""
        base = Path(self.base_dir)
        test_dir = Path("test_dir")
        output_path = basic.make_new_dir(base, test_dir)
        exp_dir = "test_dir"
        exp_path = Path(base, exp_dir)
        with self.subTest():
            self.assertTrue(exp_path.is_dir())
        with self.subTest():
            self.assertEqual(exp_dir, output_path.stem)

    def test_make_new_dir_2(self):
        """Verify no new directory is created when a directory already
        exists and attempt = 1."""
        base = Path(self.base_dir)
        new_dir = Path("test_dir")
        Path(base, new_dir).mkdir()
        output_path = basic.make_new_dir(base, new_dir)
        self.assertIsNone(output_path)

    def test_make_new_dir_3(self):
        """Verify new directory is created when a directory already
        exists and attempt = 2."""
        base = Path(self.base_dir)
        new_dir = Path("test_dir")
        Path(base, new_dir).mkdir()
        output_path = basic.make_new_dir(base, new_dir, attempt=2)
        exp_dir = Path("test_dir_1")
        exp_path = Path(base, exp_dir)
        with self.subTest():
            self.assertTrue(exp_path.is_dir())
        with self.subTest():
            self.assertEqual(exp_dir.stem, output_path.stem)

    def test_make_new_dir_4(self):
        """Verify new directory is created when two directories already
        exists and attempt = 3."""
        base = Path(self.base_dir)
        new_dir = Path("test_dir")
        Path(base, new_dir).mkdir()
        Path(base, new_dir.stem + "_1").mkdir()
        output_path = basic.make_new_dir(base, new_dir, attempt=3)
        exp_dir = Path("test_dir_2")
        exp_path = Path(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(os.path.isdir(exp_path))
        with self.subTest():
            self.assertEqual(exp_dir.stem, output_path.stem)



    def test_identify_files_1(self):
        """Verify the correct number of elements are returned when
        no ignore set is provided."""
        base = Path(self.base_dir)
        Path(base, "new_dir").mkdir()
        Path(base, "file1.txt").touch()
        Path(base, ".DS_Store").touch()
        list_of_files = basic.identify_files(base)
        exp_num_files = 2
        self.assertEqual(len(list_of_files), exp_num_files)

    def test_identify_files_2(self):
        """Verify the correct number of elements are returned when
        an ignore set is provided."""
        base = Path(self.base_dir)
        Path(base, "new_dir").mkdir()
        Path(base, "file1.txt").touch()
        Path(base, ".DS_Store").touch()
        suffix_set = set([".DS_Store"])
        list_of_files = basic.identify_files(base, suffix_set)
        exp_num_files = 1
        self.assertEqual(len(list_of_files), exp_num_files)

    def tearDown(self):
        shutil.rmtree(self.base_dir)




if __name__ == '__main__':
    unittest.main()
