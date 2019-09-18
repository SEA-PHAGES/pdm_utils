"""Tests the functionality of the unique functions in the database_to_file pipeline"""

import unittest
from classes import genome
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pipelines.db_export import database_to_file
import os, sys, shutil

class TestDatabaseToFile(unittest.TestCase):

    def setUp(self):
        self.test_seq_record_list = []
        self.test_format = "gb"

    def test_seqfeature_file_output_1(self):
        database_to_file.seqfeature_file_output(self.test_seq_record_list, self.test_format, os.getcwd())
        self.assertTrue(os.path.exists(os.path.join(os.getcwd(), "database_export_output")))
        self.assertFalse(os.listdir(os.path.join(os.getcwd(), "database_export_output")))
        shutil.rmtree(os.path.join(os.getcwd(), "database_export_output"))

