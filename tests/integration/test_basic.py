"""Integration tests for general functions."""


from pdm_utils.functions import basic
from datetime import datetime
import unittest
from unittest.mock import patch
import os
import shutil




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




    # Using the patch decorator to bypass the user-requried input() value.
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




    def test_identify_files_1(self):
        """Verify only correct number of elements are returned.
        There is one folder in the test directory, so the total number
        of files should be 1 less than the total contents of the
        directory."""
        list_of_directory_contents = os.listdir(self.test_directory1)
        expected = len(list_of_directory_contents) - 1
        list_of_files = basic.identify_files(self.test_directory1)
        self.assertEqual(len(list_of_files), expected)




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




class TestBasicFunctions2(unittest.TestCase):

    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_make_new_dir")
        os.mkdir(self.base_dir)




    def test_make_new_dir_1(self):
        """Verify new directory is created."""
        new_dir = "test_dir"
        output_dir = basic.make_new_dir(self.base_dir, new_dir)
        exp_dir = "test_dir"
        exp_path = os.path.join(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(os.path.isdir(exp_path))
        with self.subTest():
            self.assertEqual(exp_dir, output_dir)

    def test_make_new_dir_2(self):
        """Verify no new directory is created when a directory already
        exists and attempt = 1."""
        new_dir = "test_dir"
        os.mkdir(os.path.join(self.base_dir, new_dir))
        output_dir = basic.make_new_dir(self.base_dir, new_dir)
        self.assertEqual(output_dir, "")



    def test_make_new_dir_3(self):
        """Verify new directory is created when a directory already
        exists and attempt = 2."""
        new_dir = "test_dir"
        os.mkdir(os.path.join(self.base_dir, new_dir))
        output_dir = basic.make_new_dir(self.base_dir, new_dir, attempt=2)
        exp_dir = "test_dir_1"
        exp_path = os.path.join(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(os.path.isdir(exp_path))
        with self.subTest():
            self.assertEqual(output_dir, exp_dir)

    def test_make_new_dir_3(self):
        """Verify new directory is created when two directories already
        exists and attempt = 3."""
        new_dir = "test_dir"
        os.mkdir(os.path.join(self.base_dir, new_dir))
        os.mkdir(os.path.join(self.base_dir, new_dir + "_1"))
        output_dir = basic.make_new_dir(self.base_dir, new_dir, attempt=3)
        exp_dir = "test_dir_2"
        exp_path = os.path.join(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(os.path.isdir(exp_path))
        with self.subTest():
            self.assertEqual(output_dir, exp_dir)




    def tearDown(self):
        shutil.rmtree(self.base_dir)




if __name__ == '__main__':
    unittest.main()
