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

        #Creates a test database
        self.connection = pymysql.connect(host="localhost",
                                     user="pdm_anon",
                                     password="pdm_anon",
                                     cursorclass=pymysql.cursors.DictCursor)
        cur = (self.connection).cursor()
        cur.execute("SELECT SCHEMA_NAME FROM "
                    "INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = 'test_db'")
        result = cur.fetchall()
        if len(result) != 0:
            cur.execute("DROP DATABASE test_db")
            (self.connection).commit()

        cur.execute("CREATE DATABASE test_db;")
        (self.connection).commit()
        (self.connection).close
        #Creates valid alchemist
        #Creates Genome object
        gnm = genome.Genome()
        gnm.id = "TestID"
        gnm.accession = "TestAccession"
        gnm.name = "Test"
        gnm.host_genus = "TestHost"
        gnm.length = "TestLength"
        gnm.date = "TestDate"
        gnm.description = "TestDescription"
        gnm.gc = "TestGC"
        gnm.cluster = "TestCluster"
        gnm.subcluster = "TestSubcluster"
        gnm.annotation_status = "TestStatus"
        gnm.retrieve_record = "TestRecord"
        gnm.annotation_author = "TestAuthor"
        self.genome = gnm
        #Creates SeqRecord object
        seqrecord = SeqRecord(Seq("ATGC"))
        seqrecord.seq.alphabet = IUPACAmbiguousDNA()
        seqrecord.id = "Test_Accession"
        seqrecord.name = "Test"
        self.test_record = seqrecord
        #Creates working test directory
        self.test_cwd = (Path.cwd()).joinpath("DELETE_ME")
        (self.test_cwd).mkdir()
 
    def tearDown(self):
        #Tears down current working test directory
        shutil.rmtree(str(self.test_cwd))
        #Tears down test database
        cur = (self.connection).cursor()
        cur.execute("DROP DATABASE test_db;")
        (self.connection).commit()
        (self.connection).close()

if __name__ == "__main__":
    unittest.main()
