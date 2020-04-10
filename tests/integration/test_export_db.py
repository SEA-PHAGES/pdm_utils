"""Tests the functionality of the unique functions in the export_db pipeline"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import *
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from pdm_utils.classes import genome, filter
from pdm_utils.pipelines import export_db

from pathlib import Path
from unittest.mock import patch, Mock, call, MagicMock
import os, sys, unittest, shutil, csv, pymysql

class TestFileExport(unittest.TestCase):

    def setUp(self):
        pass
 
    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
