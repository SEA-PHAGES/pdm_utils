"""Unit tests for the configfile module."""

import unittest

from pdm_utils.functions import configfile

MYSQL = "mysql"

class TestConfigFile(unittest.TestCase):

    def test_default_sections_keys_1(self):
        """Confirm that dictionary of sections and keys is constructed."""
        dict = configfile.default_sections_keys()
        self.assertTrue(MYSQL in dict.keys())

    def test_setup_section_1(self):
        """Confirm that default string settings are added."""
        keys = {"a", "b", "c"}
        value = "default"
        dict = configfile.setup_section(keys, value)
        for key in keys:
            with self.subTest():
                self.assertEqual(dict[key], value)

    def test_setup_section_2(self):
        """Confirm that default None settings are added."""
        keys = {"a", "b", "c"}
        value = None
        dict = configfile.setup_section(keys, value)
        for key in keys:
            with self.subTest():
                self.assertIsNone(dict[key])


    def test_default_parser_1(self):
        """Confirm that ConfigParser is constructed with default string."""
        section = "mysql"
        key = "user"
        null_value = "default"
        parser = configfile.default_parser(null_value)
        for key in parser.keys():
            section = parser[key]
            for key in section.keys():
                with self.subTest():
                    self.assertEqual(section[key], null_value)

    def test_default_parser_2(self):
        """Confirm that ConfigParser is constructed with default None."""
        section = "mysql"
        key = "user"
        null_value = None
        parser = configfile.default_parser(null_value)
        for key in parser.keys():
            section = parser[key]
            for key in section.keys():
                with self.subTest():
                    self.assertIsNone(section[key])



if __name__ == '__main__':
    unittest.main()
