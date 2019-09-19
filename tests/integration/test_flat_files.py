"""Integration tests for misc. functions that interact with
GenBank-formatted flat files."""


import unittest
from pdm_utils.functions import flat_files
from pdm_utils.classes import genome
import os
from Bio.SeqRecord import SeqRecord

class TestFlatFileFunctions(unittest.TestCase):


    def setUp(self):

        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_1.gb")

        self.test_filepath2 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_2.gb")

        self.test_filepath3 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_5.gb")




    def test_retrieve_genome_data_1(self):
        """Verify file with one GenBank-formatted flat file
        record is correctly parsed."""
        record = flat_files.retrieve_genome_data(self.test_filepath1)
        self.assertIsInstance(record, SeqRecord)

    def test_retrieve_genome_data_2(self):
        """Verify file with no GenBank-formatted flat file
        record is correctly parsed."""
        record = flat_files.retrieve_genome_data(self.test_filepath2)
        self.assertIsNone(record)

    def test_retrieve_genome_data_3(self):
        """Verify file with two GenBank-formatted flat file
        records is correctly parsed."""
        record = flat_files.retrieve_genome_data(self.test_filepath3)
        self.assertIsNone(record)




    # def test_parse_files_1(self):
    #     """Verify file with one GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     files = [self.test_filepath1]
    #     genomes, successes, fails = flat_files.parse_files(files)
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 0)
    #
    #
    # def test_parse_files_2(self):
    #     """Verify file with no GenBank-formatted flat file
    #     record is correctly parsed."""
    #
    #     files = [self.test_filepath2]
    #     genomes, successes, fails = flat_files.parse_files(files)
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 0)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 0)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)
    #
    #
    # def test_parse_files_3(self):
    #     """Verify one file with GenBank-formatted flat file
    #     record and one with no record is correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath2]
    #     genomes, successes, fails = flat_files.parse_files(files)
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 1)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)
    #
    #
    # def test_parse_files_4(self):
    #     """Verify two files with GenBank-formatted flat file
    #     records are correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath1]
    #     genomes, successes, fails = flat_files.parse_files(files)
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 0)
    #
    #
    # def test_parse_files_5(self):
    #     """Verify two files with GenBank-formatted flat file
    #     record and one with no record is correctly parsed."""
    #
    #     files = [self.test_filepath1, self.test_filepath2, self.test_filepath1]
    #     genomes, successes, fails = flat_files.parse_files(files)
    #     with self.subTest():
    #         self.assertEqual(len(genomes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(successes), 2)
    #     with self.subTest():
    #         self.assertEqual(len(fails), 1)





if __name__ == '__main__':
    unittest.main()
