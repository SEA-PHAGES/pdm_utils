"""Integration tests for configfile module."""

import configparser
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

from pdm_utils.functions import configfile

# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_config")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder that will get created/removed for each test.
test_folder = Path(test_root_dir, "output")
config_filename = "config"
config_filepath = Path(test_folder, config_filename + ".txt")

USER = "pdm_anon"
PWD = "pdm_anon"
TOOL = "test"

def build_parser(mysql=False, password="", ncbi=False):
    """Creates a mock config to save to file."""
    parser = configparser.ConfigParser(allow_no_value=True)
    if mysql:
        parser["mysql"] = {}
        parser["mysql"]["user"] = USER
        parser["mysql"]["password"] = password
    if ncbi:
        parser["ncbi"] = {}
        parser["ncbi"]["tool"] = TOOL
        parser["ncbi"]["email"] = None
    return parser

def create_config_file(parser, filepath):
    with filepath.open("w") as fh:
        parser.write(fh)

def open_config_file(filepath):
    parser = configparser.ConfigParser()
    parser.read(filepath)
    return parser


class TestConfigFile(unittest.TestCase):

    def setUp(self):
        test_folder.mkdir()

    def tearDown(self):
        shutil.rmtree(test_folder)




    def test_parse_config_1(self):
        """Confirm that parser is NOT complete when NO input parser is provided."""
        parser1 = build_parser(mysql=True)
        create_config_file(parser1, config_filepath)
        parser2 = configfile.parse_config(config_filepath, parser=None)
        with self.subTest():
            self.assertEqual(parser2["mysql"]["user"], USER)
        with self.subTest():
            self.assertEqual(parser2["mysql"]["password"], "")
        with self.subTest():
            self.assertTrue("invalid" not in parser2["mysql"].keys())
        with self.subTest():
            self.assertTrue("ncbi" not in parser2.keys())

    def test_parse_config_2(self):
        """Confirm that parser is complete when input parser is provided,
        and confirm that password value in first parser is replaced with
        new value from file."""
        parser1 = build_parser(mysql=True, password=PWD)
        create_config_file(parser1, config_filepath)
        parser2 = build_parser(mysql=True, password="old", ncbi=True)
        parser3 = configfile.parse_config(config_filepath, parser=parser2)
        with self.subTest():
            self.assertEqual(parser3["mysql"]["user"], USER)
        with self.subTest():
            self.assertEqual(parser3["mysql"]["password"], PWD)
        with self.subTest():
            self.assertEqual(parser3["ncbi"]["tool"], TOOL)
        with self.subTest():
            self.assertIsNone(parser3["ncbi"]["email"])

    def test_build_complete_config_1(self):
        """Confirm that ConfigParser is constructed with default None
        using valid file."""
        parser1 = build_parser(mysql=True)
        create_config_file(parser1, config_filepath)
        parser2 = configfile.build_complete_config(config_filepath)
        with self.subTest():
            self.assertEqual(parser2["mysql"]["user"], USER)
        with self.subTest():
            self.assertIsNone(parser2["ncbi"]["tool"])
        with self.subTest():
            self.assertIsNone(parser2["ncbi"]["email"])

    def test_build_complete_config_2(self):
        """Confirm that ConfigParser is constructed with default None
        using no file."""
        parser2 = configfile.build_complete_config(None)
        with self.subTest():
            self.assertIsNone(parser2["mysql"]["user"])
        with self.subTest():
            self.assertIsNone(parser2["ncbi"]["tool"])
        with self.subTest():
            self.assertIsNone(parser2["ncbi"]["email"])

    def test_write_config_1(self):
        """Confirm that ConfigParser is written to file."""
        parser1 = build_parser(mysql=True)
        configfile.write_config(parser1, config_filepath)
        parser2 = open_config_file(config_filepath)
        with self.subTest():
            self.assertEqual(parser2["mysql"]["user"], USER)
        with self.subTest():
            self.assertTrue("ncbi" not in parser2.keys())

    def test_create_empty_config_file_1(self):
        """Confirm that ConfigParser is written to file with "empty" values."""
        value = "empty"
        configfile.create_empty_config_file(test_folder, config_filename, value)
        parser = open_config_file(config_filepath)
        with self.subTest():
            self.assertEqual(parser["mysql"]["user"], value)
        with self.subTest():
            self.assertEqual(parser["ncbi"]["tool"], value)



if __name__ == '__main__':
    unittest.main()
