import csv
import os
import shutil 
import sys 
import unittest
from pathlib import Path
from unittest.mock import call
from unittest.mock import MagicMock
from unittest.mock import Mock 
from unittest.mock import patch 

from Bio.Alphabet.IUPAC import *
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.classes import genome
from pdm_utils.pipelines import export_db


# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

# pdm_anon, pdm_anon, and pdm_test_db
USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
TEST_DIR = "/tmp/pdm_utils_tests_export"

class TestFileExport(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

        self.test_dir = Path(TEST_DIR)
        if self.test_dir.is_dir():
            shutil.rmtree(TEST_DIR)

        self.test_dir.mkdir()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()
        shutil.rmtree(TEST_DIR)

    def setUp(self):
        self.alchemist = AlchemyHandler()
        self.alchemist.username=USER
        self.alchemist.password=PWD
        self.alchemist.database=DB
        self.alchemist.connect(ask_database=True, login_attempts=0)
        self.alchemist.build_graph()

        self.db_filter = Filter(alchemist=self.alchemist)
        
        self.export_test_dir = self.test_dir.joinpath("export_test_dir")

    def tearDown(self):
        if self.export_test_dir.is_dir():
            shutil.rmtree(str(self.export_test_dir))

    def test_execute_export_1(self):
        """Verify execute_export() creates new directory as expected.
        """
        for pipeline in export_db.PIPELINES:
            with self.subTest(pipeline=pipeline):
                export_db.execute_export(self.alchemist, self.test_dir, 
                                         self.export_test_dir.name, pipeline)
                self.assertTrue(self.export_test_dir.is_dir())
                shutil.rmtree(str(self.export_test_dir))

    def test_execute_export_2(self):
        """Verify execute_export() 'sql' pipeline functions as expected.
        """
        export_db.execute_export(self.alchemist, self.test_dir,
                                  self.export_test_dir.name, "sql")

        self.assertTrue(self.export_test_dir.is_dir())
        
        sql_file_path = self.export_test_dir.joinpath(
                                            f"{self.alchemist.database}.sql")
        self.assertTrue(sql_file_path.is_file())

    def test_execute_export_3(self):
        """Verify execute_export() 'csv' pipeline functions as expected.
        """
        for table in export_db.TABLES:
            with self.subTest(table=table):
                export_db.execute_export(self.alchemist, self.test_dir,
                                          self.export_test_dir.name, "csv",
                                          table=table)
                self.assertTrue(self.export_test_dir.is_dir())

                csv_file_path = self.export_test_dir.joinpath(
                                            f"{table}.csv")

                self.assertTrue(csv_file_path.is_file())

                shutil.rmtree(str(self.export_test_dir))

    def test_execute_export_4(self):
        """Verify execute_export() SeqRecord pipelines function as expected.
        """
        for file_type in export_db.BIOPYTHON_PIPELINES:
            with self.subTest(file_type=file_type):
                export_db.execute_export(self.alchemist, self.test_dir,
                                         self.export_test_dir.name, file_type)
                self.assertTrue(self.export_test_dir.is_dir())

                flat_file_path = self.export_test_dir.joinpath(
                                            f"Trixie.{file_type}")
                self.assertTrue(flat_file_path.is_file())

                shutil.rmtree(str(self.export_test_dir))
   
    def test_execute_export_5(self):
        """Verify execute_export() filter parameter functions as expected.
        """
        filters = "phage.PhageID!=Trixie AND phage.Cluster=A"
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "fasta",
                                 filters=filters)

        D29_file_path = self.export_test_dir.joinpath("D29.fasta")
        Trixie_file_path = self.export_test_dir.joinpath("Trixie.fasta")

        self.assertTrue(D29_file_path.is_file())
        self.assertFalse(Trixie_file_path.is_file())
    
    def test_execute_export_6(self):
        """Verify execute_export() group parameter functions as expected.
        """
        groups = ["phage.Cluster", "phage.Subcluster"]
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "fasta",
                                 groups=groups)

        A_path = self.export_test_dir.joinpath("A")
        C_path = self.export_test_dir.joinpath("C")

        A2_path = A_path.joinpath("A2")
        C1_path = C_path.joinpath("C1")
        C2_path = C_path.joinpath("C2")

        Trixie_path = A2_path.joinpath("Trixie.fasta")
        D29_path = A2_path.joinpath("D29.fasta")
        Alice_path = C1_path.joinpath("Alice.fasta")
        Myrna_path = C2_path.joinpath("Myrna.fasta")

                            
        self.assertTrue(A_path.is_dir())
        self.assertTrue(C_path.is_dir())

        self.assertTrue(A2_path.is_dir())
        self.assertTrue(C1_path.is_dir())
        self.assertTrue(C2_path.is_dir())

        self.assertTrue(Trixie_path.is_file())
        self.assertTrue(D29_path.is_file())
        self.assertTrue(Alice_path.is_file())
        self.assertTrue(Myrna_path.is_file())
        
    def test_execute_export_7(self):
        """Verify execute_export() sort parameter is functional.
        """
        sort_columns = ["phage.Subcluster"]
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "csv",
                                 sort=sort_columns)
    
    def test_execute_export_8(self):
        """Verify execute_export() concatenate parameter functions as expected.
        """
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "fasta",
                                 concatenate=True)

        fasta_path = self.test_dir.joinpath(
                                        f"{self.export_test_dir.name}.fasta")

        self.assertTrue(fasta_path.is_file())
        self.assertFalse(self.export_test_dir.is_file())

    def test_execute_export_9(self):
        """Verify execute_export() include_columns functions as expected.
        """
        include_columns = ["phage.Cluster"]
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "csv", table="gene",
                                 include_columns=include_columns)

        csv_path = self.export_test_dir.joinpath(
                                        f"gene.csv")

        with open(csv_path) as csv_handle:
            reader = csv.reader(csv_handle)
            headers = next(reader)

        self.assertTrue("Cluster" in headers)
        self.assertEqual("GeneID", headers[0])
        self.assertFalse("Translation" in headers)
    
    def test_execute_export_10(self):
        """Verify execute_export() exclude_columns functions as expected.
        """
        exclude_columns = ["phage.Subcluster"]
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "csv",
                                 exclude_columns=exclude_columns)
        
        csv_path = self.export_test_dir.joinpath(
                                        f"phage.csv")

        with open(csv_path) as csv_handle:
            reader = csv.reader(csv_handle)
            headers = next(reader)

        self.assertTrue("Cluster" in headers)
        self.assertEqual("PhageID", headers[0])
        self.assertFalse("Subcluster" in headers)
        self.assertFalse("Sequence" in headers)

    def test_execute_export_11(self):
        """Verify execute_export() sequence_columns functions as expected.
        """
        export_db.execute_export(self.alchemist, self.test_dir,
                                 self.export_test_dir.name, "csv",
                                 sequence_columns=True)

        csv_path = self.export_test_dir.joinpath(
                                        f"phage.csv")

        with open(csv_path) as csv_handle:
            reader = csv.reader(csv_handle)
            headers = next(reader)

        self.assertTrue("Cluster" in headers)
        self.assertEqual("PhageID", headers[0])
        self.assertTrue("Sequence" in headers)

    def test_execute_export_12(self):
        """Verify execute_export() SeqRecord pipeline functions as expected.
        """
        for file_type in export_db.BIOPYTHON_PIPELINES:
            with self.subTest(file_type=file_type):
                export_db.execute_export(self.alchemist, self.test_dir,
                                         self.export_test_dir.name, file_type,
                                         table="gene")

                self.assertTrue(self.export_test_dir.is_dir())

                flat_file_path = self.export_test_dir.joinpath(
                                            f"Trixie_CDS_70.{file_type}")
                self.assertTrue(flat_file_path.is_file())

                shutil.rmtree(str(self.export_test_dir))
 
if __name__ == "__main__":
    unittest.main()
