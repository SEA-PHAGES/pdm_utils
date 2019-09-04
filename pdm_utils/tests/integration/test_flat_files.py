"""Integration tests for misc. functions that interact with
GenBank-formatted flat files."""


import unittest
from functions import flat_files
from classes import genome
import os

class TestFlatFileFunctions(unittest.TestCase):


    def setUp(self):

        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")

        self.test_filepath2 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_2.gb")

        self.test_filepath3 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_8.ggb")



    # TODO this may no longer be needed.
    # def test_create_parsed_flat_file_1(self):
    #     """Verify file with one GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     gnm = flat_files.create_parsed_flat_file(self.test_filepath1)
    #     self.assertEqual(gnm.id, "L5")
    #
    #
    # def test_create_parsed_flat_file_2(self):
    #     """Verify file with no GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     gnm = flat_files.create_parsed_flat_file(self.test_filepath2)
    #     self.assertEqual(gnm.id, "")
    #
    #
    # def test_create_parsed_flat_file_3(self):
    #     """Verify file with GenBank-formatted record but incorrect
    #     file extension is correctly parsed."""
    #
    #     gnm = flat_files.create_parsed_flat_file(self.test_filepath3)
    #     self.assertEqual(gnm.id, "")
    #
    #
    #
    #
    # def test_create_parsed_flat_file_list_1(self):
    #     """Verify file with one GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     files = [self.test_filepath1]
    #     genomes, successes, fails = \
    #         flat_files.create_parsed_flat_file_list(files)
    #
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 0)
    #
    #
    # def test_create_parsed_flat_file_list_2(self):
    #     """Verify file with no GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     files = [self.test_filepath2]
    #     genomes, successes, fails = \
    #         flat_files.create_parsed_flat_file_list(files)
    #
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 0)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)
    #
    #
    # def test_create_parsed_flat_file_list_3(self):
    #     """Verify one file with GenBank-formatted flat file
    #     record and one with no record is correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath2]
    #     genomes, successes, fails = \
    #         flat_files.create_parsed_flat_file_list(files)
    #
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)
    #
    #
    # def test_create_parsed_flat_file_list_4(self):
    #     """Verify two files with GenBank-formatted flat file
    #     records are correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath1]
    #     genomes, successes, fails = \
    #         flat_files.create_parsed_flat_file_list(files)
    #
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 0)
    #
    #
    # def test_create_parsed_flat_file_list_5(self):
    #     """Verify two files with GenBank-formatted flat file
    #     record and one with no record is correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath2, self.test_filepath1]
    #
    #     genomes, successes, fails = \
    #         flat_files.create_parsed_flat_file_list(files)
    #
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 3)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)





    def test_parse_files_1(self):
        """Verify file with one GenBank-formatted flat file
        record is correctly parsed."""

        files = [self.test_filepath1]
        genomes, successes, fails = flat_files.parse_files(files)
        with self.subTest():
            self.assertEqual(len(genomes), 1)
        with self.subTest():
            self.assertEqual(len(successes), 1)
        with self.subTest():
            self.assertEqual(len(fails), 0)


    def test_parse_files_2(self):
        """Verify file with no GenBank-formatted flat file
        record is correctly parsed."""

        files = [self.test_filepath2]
        genomes, successes, fails = flat_files.parse_files(files)
        with self.subTest():
            self.assertEqual(len(genomes), 0)
        with self.subTest():
            self.assertEqual(len(successes), 0)
        with self.subTest():
            self.assertEqual(len(fails), 1)


    def test_parse_files_3(self):
        """Verify one file with GenBank-formatted flat file
        record and one with no record is correctly parsed."""

        files = [self.test_filepath1, self.test_filepath2]
        genomes, successes, fails = flat_files.parse_files(files)
        with self.subTest():
            self.assertEqual(len(genomes), 1)
        with self.subTest():
            self.assertEqual(len(successes), 1)
        with self.subTest():
            self.assertEqual(len(fails), 1)


    def test_parse_files_4(self):
        """Verify two files with GenBank-formatted flat file
        records are correctly parsed."""

        files = [self.test_filepath1, self.test_filepath1]
        genomes, successes, fails = flat_files.parse_files(files)
        with self.subTest():
            self.assertEqual(len(genomes), 2)
        with self.subTest():
            self.assertEqual(len(successes), 2)
        with self.subTest():
            self.assertEqual(len(fails), 0)


    def test_parse_files_5(self):
        """Verify two files with GenBank-formatted flat file
        record and one with no record is correctly parsed."""

        files = [self.test_filepath1, self.test_filepath2, self.test_filepath1]
        genomes, successes, fails = flat_files.parse_files(files)
        with self.subTest():
            self.assertEqual(len(genomes), 2)
        with self.subTest():
            self.assertEqual(len(successes), 2)
        with self.subTest():
            self.assertEqual(len(fails), 1)





if __name__ == '__main__':
    unittest.main()
