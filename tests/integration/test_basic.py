"""Integration tests for general functions."""

import argparse
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from pdm_utils.functions import basic

unittest_file = Path(__file__)
unittest_dir = unittest_file.parent

TMPDIR_PREFIX = "pdm_utils_basic_tests_"
# Can set TMPDIR_BASE to string such as "/tmp/" to track tmp directory location.
TMPDIR_BASE = "/tmp"


class TestBasicFunctions1(unittest.TestCase):

    def setUp(self):

        self.test_filepath1 = Path(unittest_dir,
                                "test_files/test_flat_file_1.gb")

        self.test_directory1 = Path(unittest_dir,
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




    @patch("builtins.input", return_value = "YES")
    def test_ask_yes_no_1(self, input_mock):
        """Verify output when user answers provides YES response."""
        response = basic.ask_yes_no()
        self.assertTrue(response)

    @patch("builtins.input", return_value = "NO")
    def test_ask_yes_no_2(self, input_mock):
        """Verify output when user answers provides NO response."""
        response = basic.ask_yes_no()
        self.assertFalse(response)

    @patch("builtins.input", return_value = "EXIT")
    def test_ask_yes_no_3(self, input_mock):
        """Verify output when user answers provides exit response."""
        response = basic.ask_yes_no(response_attempt=3)
        self.assertIsNone(response)

    @patch("builtins.input", return_value = "invalid")
    def test_ask_yes_no_4(self, input_mock):
        """Verify output when user answers provides invalid response."""
        response = basic.ask_yes_no()
        self.assertIsNone(response)

    @patch("builtins.input")
    def test_ask_yes_no_5(self, input_mock):
        """Verify output when attempts = 3 and
        user provides valid response on last attempt."""
        input_mock.side_effect = ["invalid", "invalid", "Y"]
        response = basic.ask_yes_no(response_attempt=3)
        self.assertTrue(response)




    def test_expand_path_1(self):
        """Verify that home path is expanded."""
        partial_path = "/fake/path"
        input_path = "~" + partial_path
        expanded_path = basic.expand_path(input_path)
        home_dir = Path("~").expanduser()
        expected_path = str(home_dir) + partial_path
        self.assertEqual(expanded_path, expected_path)

    def test_expand_path_2(self):
        """Verify no change when an absoluate path is provided."""
        input_path = "/fake/path"
        expanded_path = basic.expand_path(input_path)
        expected_path = input_path
        self.assertEqual(expanded_path, expected_path)

    def test_expand_path_3(self):
        """Verify absolute path is expanded when a relative path is provided."""
        partial_path = "/fake/path"
        input_path = "." + partial_path
        expanded_path = basic.expand_path(input_path)
        local_path = Path(".").resolve()
        expected_path = str(local_path) + partial_path
        self.assertEqual(expanded_path, expected_path)




    def test_verify_path_1(self):
        """Verify a true file produces no error when 'file' is indicated."""
        result = basic.verify_path(self.test_filepath1, "file")
        self.assertTrue(result)

    def test_verify_path_2(self):
        """Verify a fake file produces an error when 'file' is indicated."""
        result = basic.verify_path(str(self.test_filepath1) + "abcxyz", "file")
        self.assertFalse(result)

    def test_verify_path_3(self):
        """Verify a true directory produces no error when
        'directory' is indicated."""
        result = basic.verify_path(str(self.test_directory1), "dir")
        self.assertTrue(result)

    def test_verify_path_4(self):
        """Verify a fake directory produces an error when
        'directory' is indicated."""
        result = basic.verify_path(str(self.test_directory1) + "abcxyz", "dir")
        self.assertFalse(result)

    def test_verify_path_5(self):
        """Verify a true directory produces no error when
        'None' is indicated."""
        result = basic.verify_path(str(self.test_directory1))
        self.assertTrue(result)

    def test_verify_path_6(self):
        """Verify a fake directory produces an error when
        'None' is indicated."""
        result = basic.verify_path(str(self.test_directory1) + "abcxyz")
        self.assertFalse(result)

    def test_verify_path_7(self):
        """Verify a true directory produces an error when
        'invalid' is indicated."""
        result = basic.verify_path(str(self.test_directory1), "invalid")
        self.assertFalse(result)




    def test_parse_flag_file_1(self):
        """Verify that all three rows are parsed correctly."""
        flag_file = Path(unittest_dir, "test_files/test_flag_file_1.csv")
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
        flag_file = Path(unittest_dir, "test_files/test_flag_file_2.csv")
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
        self.tmpdir = tempfile.TemporaryDirectory(prefix=TMPDIR_PREFIX,
                                                  dir=TMPDIR_BASE)
        self.base_dir = Path(self.tmpdir.name)

    def tearDown(self):
        print(self.tmpdir.name)
        input("check")
        self.tmpdir.cleanup()




    def test_make_new_dir_1(self):
        """Verify new directory is created."""
        test_dir = Path("test_dir")
        output_path = basic.make_new_dir(self.base_dir, test_dir)
        exp_dir = "test_dir"
        exp_path = Path(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(exp_path.is_dir())
        with self.subTest():
            self.assertEqual(exp_dir, output_path.stem)

    def test_make_new_dir_2(self):
        """Verify no new directory is created when a directory already
        exists and attempt = 1."""
        new_dir = Path("test_dir")
        Path(self.base_dir, new_dir).mkdir()
        output_path = basic.make_new_dir(self.base_dir, new_dir)
        self.assertIsNone(output_path)

    def test_make_new_dir_3(self):
        """Verify new directory is created when a directory already
        exists and attempt = 2."""
        new_dir = Path("test_dir")
        Path(self.base_dir, new_dir).mkdir()
        output_path = basic.make_new_dir(self.base_dir, new_dir, attempt=2)
        exp_dir = Path("test_dir_1")
        exp_path = Path(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(exp_path.is_dir())
        with self.subTest():
            self.assertEqual(exp_dir.stem, output_path.stem)

    def test_make_new_dir_4(self):
        """Verify new directory is created when two directories already
        exists and attempt = 3."""
        new_dir = Path("test_dir")
        Path(self.base_dir, new_dir).mkdir()
        Path(self.base_dir, new_dir.stem + "_1").mkdir()
        output_path = basic.make_new_dir(self.base_dir, new_dir, attempt=3)
        exp_dir = Path("test_dir_2")
        exp_path = Path(self.base_dir, exp_dir)
        with self.subTest():
            self.assertTrue(exp_path.is_dir())
        with self.subTest():
            self.assertEqual(exp_dir.stem, output_path.stem)



    def test_identify_files_1(self):
        """Verify the correct number of elements are returned when
        no ignore set is provided."""
        Path(self.base_dir, "new_dir").mkdir()
        Path(self.base_dir, "file1.txt").touch()
        Path(self.base_dir, ".DS_Store").touch()
        list_of_files = basic.identify_files(self.base_dir)
        exp_num_files = 2
        self.assertEqual(len(list_of_files), exp_num_files)

    def test_identify_files_2(self):
        """Verify the correct number of elements are returned when
        an ignore set is provided."""
        Path(self.base_dir, "new_dir").mkdir()
        Path(self.base_dir, "file1.txt").touch()
        Path(self.base_dir, ".DS_Store").touch()
        suffix_set = set([".DS_Store"])
        list_of_files = basic.identify_files(self.base_dir, suffix_set)
        exp_num_files = 1
        self.assertEqual(len(list_of_files), exp_num_files)




class TestBasicFunctions3(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory(prefix=TMPDIR_PREFIX,
                                                  dir=TMPDIR_BASE)
        self.test_directory1 = Path(self.tmpdir.name)
        self.file = Path(self.test_directory1, "new_file.txt")
        self.dir = Path(self.test_directory1, "new_dir")

    def tearDown(self):
        self.tmpdir.cleanup()




    def test_verify_path2_1(self):
        """Verify output when file exists, is expected to exist,
        and kind = 'file'."""
        self.file.touch()
        result, msg = basic.verify_path2(self.file, kind="file", expect=True)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_2(self):
        """Verify output when file exists, is not expected to exist,
        and kind = 'file'."""
        self.file.touch()
        result, msg = basic.verify_path2(self.file, kind="file", expect=False)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_3(self):
        """Verify output when file does not exist, is not expected to exist,
        and kind = 'file'."""
        result, msg = basic.verify_path2(self.file, kind="file", expect=False)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_4(self):
        """Verify output when file does not exist, is expected to exist,
        and kind = 'file'."""
        result, msg = basic.verify_path2(self.file, kind="file", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_5(self):
        """Verify output when file exists, is expected to exist,
        and kind = 'dir'."""
        self.file.touch()
        result, msg = basic.verify_path2(self.dir, kind="dir", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_6(self):
        """Verify output when file exists, is expected to exist,
        and kind = None."""
        self.file.touch()
        result, msg = basic.verify_path2(self.file, kind=None, expect=True)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_7(self):
        """Verify output when file exists, is not expected to exist,
        and kind = None."""
        self.file.touch()
        result, msg = basic.verify_path2(self.file, kind=None, expect=False)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_8(self):
        """Verify output when file does not exist, is not expected to exist,
        and kind = None."""
        result, msg = basic.verify_path2(self.file, kind=None, expect=False)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_9(self):
        """Verify output when file does not exist, is expected to exist,
        and kind = None."""
        result, msg = basic.verify_path2(self.file, kind=None, expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_10(self):
        """Verify output when file exists, is expected to exist,
        and kind = 'invalid'."""
        result, msg = basic.verify_path2(self.file, kind="invalid", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_11(self):
        """Verify output when directory exists, is expected to exist,
        and kind = 'dir'."""
        self.dir.mkdir()
        result, msg = basic.verify_path2(self.dir, kind="dir", expect=True)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_12(self):
        """Verify output when directory exists, is not expected to exist,
        and kind = 'dir'."""
        self.dir.mkdir()
        result, msg = basic.verify_path2(self.dir, kind="dir", expect=False)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_13(self):
        """Verify output when directory does not exist,
        is not expected to exist, and kind = 'dir'."""
        result, msg = basic.verify_path2(self.dir, kind="dir", expect=False)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_14(self):
        """Verify output when directory does not exist,
        is expected to exist, and kind = 'dir'."""
        result, msg = basic.verify_path2(self.dir, kind="dir", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_15(self):
        """Verify output when directory exists, is expected to exist,
        and kind = 'file'."""
        self.dir.mkdir()
        result, msg = basic.verify_path2(self.dir, kind="file", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)


    def test_verify_path2_16(self):
        """Verify output when directory exists, is expected to exist, and
        kind = None."""
        self.dir.mkdir()
        result, msg = basic.verify_path2(self.dir, kind=None, expect=True)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertIsNone(msg)


    def test_verify_path2_17(self):
        """Verify output when directory exists, is expected to exist,
        and kind = 'invalid'."""
        self.dir.mkdir()
        result, msg = basic.verify_path2(self.dir, kind="invalid", expect=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertIsNotNone(msg)



class TestBasicFunctions4(unittest.TestCase):

    def setUp(self):
        self.options = ["note", "function", "product"]




    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_1(self, ask_mock):
        """Verify no choice given with no field_set options."""
        field = basic.choose_from_list(set())
        with self.subTest():
            self.assertFalse(ask_mock.called)
        with self.subTest():
            self.assertIsNone(field)

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_2(self, ask_mock):
        """Verify first field is selected."""
        ask_mock.side_effect = [True]
        field = basic.choose_from_list(self.options)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertEqual(field, "note")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_3(self, ask_mock):
        """Verify second field is selected after one "no" response."""
        ask_mock.side_effect = [False, True]
        field = basic.choose_from_list(self.options)
        self.assertEqual(field, "function")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_4(self, ask_mock):
        """Verify exit after one "exit" response."""
        ask_mock.side_effect = [None]
        field = basic.choose_from_list(self.options)
        self.assertIsNone(field)

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_5(self, ask_mock):
        """Verify exit after one "no" response and one "exit" response."""
        ask_mock.side_effect = [False, None]
        field = basic.choose_from_list(self.options)
        self.assertIsNone(field)

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_choose_from_list_6(self, ask_mock):
        """Verify exit after three "no" responses."""
        ask_mock.side_effect = [False, False, False]
        field = basic.choose_from_list(self.options)
        self.assertIsNone(field)




class TestBasicFunctions5(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory(prefix=TMPDIR_PREFIX,
                                                  dir=TMPDIR_BASE)
        self.test_directory1 = Path(self.tmpdir.name)
        self.file = Path(self.test_directory1, "new_file.txt")
        self.dir = Path(self.test_directory1, "new_dir")
        self.parser = argparse.ArgumentParser()

    def tearDown(self):
        self.tmpdir.cleanup()


    def test_set_path_1(self):
        """Verify output when file exists and is expected to exist."""
        self.file.touch()
        # Since using tempfile, there is an added quirk.
        # the tempfile path may be a symlink, so passing it through set path
        # will resolve the symlink, changing the path, and breaking the test.
        self.file = self.file.resolve()
        output = basic.set_path(self.file, kind="file", expect=True)
        with self.subTest():
            self.assertIsInstance(output, Path)
        with self.subTest():
            self.assertEqual(str(self.file), str(output))

    @patch("sys.exit")
    def test_set_path_2(self, sys_exit_mock):
        """Verify script exits when file does not exist,
        and is not expected to exist."""
        output = basic.set_path(self.file, kind="file", expect=True)
        self.assertTrue(sys_exit_mock.called)

    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_3(self, verify_path2_mock):
        """Verify home directory expansion."""
        home = Path("~")
        home = home.expanduser()
        test_file = Path("~/path/to/file.txt")
        verify_path2_mock.return_value = (True, None)
        output = basic.set_path(test_file, kind="file", expect=True)
        exp = Path(home, "path/to/file.txt")
        self.assertEqual(output, exp)

    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_4(self, verify_path2_mock):
        """Verify '..' directory resolution."""
        test_file = Path("/dir1/dir2/../file.txt")
        verify_path2_mock.return_value = (True, None)
        output = basic.set_path(test_file, kind="file", expect=True)
        exp = Path("/dir1/file.txt")
        self.assertEqual(output, exp)

    @patch("pdm_utils.functions.basic.verify_path2")
    def test_set_path_5(self, verify_path2_mock):
        """Verify home directory expansion and '..' directory resolution."""
        home = Path("~")
        home = home.expanduser()
        test_file = Path("~/dir1/dir2/../file.txt")
        verify_path2_mock.return_value = (True, None)
        output = basic.set_path(test_file, kind="file", expect=True)
        exp = Path(home, "dir1/file.txt")
        self.assertEqual(output, exp)




class TestBasicFunctions6(unittest.TestCase):


    def setUp(self):
        self.test_import_table_1 = Path(unittest_dir,
                                     "test_files/test_import_table_1.csv")
        self.tkt_dict1 = {"phage_id": "L5", "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"phage_id": "Trixie", "host_genus": "Mycobacterium"}
        self.tmpdir = tempfile.TemporaryDirectory(prefix=TMPDIR_PREFIX,
                                                  dir=TMPDIR_BASE)
        self.base_dir = Path(self.tmpdir.name)
        self.export_file = Path(self.base_dir, "table.csv")

    def tearDown(self):
        self.tmpdir.cleanup()




    def test_retrieve_data_dict_1(self):
        """Verify a correctly structured file can be opened."""
        list_of_data_dicts = \
            basic.retrieve_data_dict(self.test_import_table_1)
        self.assertEqual(len(list_of_data_dicts), 2)




    def test_export_data_dict_1(self):
        """Verify data is exported correctly."""

        list_of_data = [self.tkt_dict1, self.tkt_dict2]
        headers = ["type", "phage_id", "host_genus", "cluster"]
        basic.export_data_dict(list_of_data, self.export_file,
                                    headers, include_headers=True)

        exp_success_tkts = []
        with open(self.export_file,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_success_tkts.append(dict)

        with self.subTest():
            self.assertEqual(len(exp_success_tkts), 2)
        with self.subTest():
            self.assertEqual(set(exp_success_tkts[0].keys()), set(headers))





if __name__ == '__main__':
    unittest.main()
