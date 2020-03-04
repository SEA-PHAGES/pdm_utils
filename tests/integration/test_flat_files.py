"""Integration tests for misc. functions that interact with
GenBank-formatted flat files."""


from pathlib import Path
import unittest
from pdm_utils.functions import flat_files
from pdm_utils.classes import genome
import os
from Bio.SeqRecord import SeqRecord

class TestFlatFileFunctions(unittest.TestCase):


    def setUp(self):

        self.unittest_file = Path(__file__)
        self.unittest_dir = self.unittest_file.parent
        self.base_dir = Path(self.unittest_dir,"test_files")
        self.test_filepath1 = Path(self.base_dir, "test_flat_file_1.gb")
        self.test_filepath2 = Path(self.base_dir, "test_flat_file_2.gb")
        self.test_filepath3 = Path(self.base_dir, "test_flat_file_5.gb")



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
