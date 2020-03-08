"""Integration tests for misc. functions that interact with
GenBank-formatted flat files."""


from pathlib import Path
import unittest
from pdm_utils.functions import flat_files
from pdm_utils.classes import genome
import os
from Bio.SeqRecord import SeqRecord


unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
test_file_dir = Path(test_dir, "test_files")


class TestFlatFileFunctions(unittest.TestCase):


    def setUp(self):
        self.test_filepath1 = Path(test_file_dir, "test_flat_file_1.gb")
        self.test_filepath2 = Path(test_file_dir, "test_flat_file_2.gb")
        self.test_filepath3 = Path(test_file_dir, "test_flat_file_5.gb")


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

if __name__ == '__main__':
    unittest.main()
