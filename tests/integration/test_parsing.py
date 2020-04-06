from pathlib import Path
from pdm_utils.functions import parsing
from sqlalchemy import Column
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import Table
from unittest.mock import Mock, patch, PropertyMock
import unittest
import re
import sys

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

class TestParsing(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    def setUp(self): 
        self.engine = create_engine(test_db_utils.create_engine_string())

        self.metadata = MetaData(bind=self.engine)
        self.metadata.reflect()

        self.phage = self.metadata.tables["phage"]
        self.gene = self.metadata.tables["gene"]

        self.PhageID = self.phage.c.PhageID
        self.PhamID = self.gene.c.PhamID
        self.Notes = self.gene.c.Notes

    def test_check_operator_1(self):
        """Verify check_operator() recognizes all operators.
        """
        for operator in parsing.OPERATORS:
            with self.subTest(operator=operator):
                parsing.check_operator(operator, self.PhamID)

    def test_check_operator_2(self):
        """Verify check_operator() raises ValueError from invalid pairing.
        check_operator() should recognize when an operator is inappropriate
        for a given Column.
        """
        for operator in parsing.NUMERIC_OPERATORS:
            with self.subTest(operator=operator):
                with self.assertRaises(ValueError):
                    parsing.check_operator(operator, self.PhageID)

        for operator in parsing.NUMERIC_OPERATORS:
            with self.subTest(operator=operator):
                with self.assertRaises(ValueError):
                    parsing.check_operator(operator, self.Notes)

    def test_check_operator_3(self):
        """Verify check_operator() recognizes valid pairings.
        """
        for operator in parsing.NONNUMERIC_OPERATORS:
            with self.subTest(operator=operator):
                parsing.check_operator(operator, self.PhageID)

        for operator in parsing.NONNUMERIC_OPERATORS:
            with self.subTest(operator=operator):
                parsing.check_operator(operator, self.Notes)

    def test_translate_table_1(self):
        """Verify translate_table() conserves case-sensitive table names.
        """
        table = parsing.translate_table(self.metadata, "phage")
  
        self.assertEqual(table, "phage")

    def test_translate_table_2(self):
        """Verify translate_table() returns case-sensitive table names.
        """
        table = parsing.translate_table(self.metadata, "GENE")

        self.assertEqual(table, "gene")

    def test_translate_table_3(self):
        """Verify translate_table() returns case-sensitive table names.
        """
        table = parsing.translate_table(self.metadata, "tRNA")
        self.assertEqual(table, "trna")

    def test_translate_table_4(self):
        """Verify translate_table() raises ValueError from invalid table name.
        """
        with self.assertRaises(ValueError):
            parsing.translate_table(self.metadata, "not_table")

    def test_translate_column_1(self):
        """Verify translate_column() conserves case-sensitive table names.
        """
        column = parsing.translate_column(self.metadata, "phage.PhageID")

        self.assertEqual(column, "PhageID")

    def test_translate_column_2(self):
        """Verify translate_column() returns case-sensitive column names.
        """
        column = parsing.translate_column(self.metadata, "GENE.GENEID")

        self.assertEqual(column, "GeneID")

    def test_translate_column_3(self):
        """Verify translate_column() returns case-sensitive column names.
        """
        column = parsing.translate_column(self.metadata, "tRNA.tRNAid")
        
        self.assertEqual(column, "TrnaID")

    def test_translate_column_4(self):
        """Verify translate_column() raises ValueError from invalid column name.
        """
        with self.assertRaises(ValueError):
            parsing.translate_column(self.metadata, "phage.not_column")

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

if __name__ == "__main__":
    unittest.main()
