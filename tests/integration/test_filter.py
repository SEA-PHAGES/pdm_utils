"""Integration unittests for the filter module"""

from networkx import Graph
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import querying
from sqlalchemy.engine.base import Engine
from sqlalchemy.sql.elements import BinaryExpression
from unittest.mock import Mock, patch
import unittest

class TestFilter(unittest.TestCase):
     def setUp(self):
         mock_alchemist = AlchemyHandler()
         mock_alchemist.username="pdm_anon"
         mock_alchemist.password="pdm_anon"
         mock_alchemist.database="Actinobacteriophage"
         mock_alchemist.connect()
         mock_alchemist.build_graph()
         self.alchemist = mock_alchemist
         self.db_filter = Filter(loader=self.alchemist)

     def test_constructor_1(self):
         db_filter = Filter(loader=self.alchemist)

     def test_constructor_2(self):
         db_filter = Filter(loader=self.alchemist.engine)

     def test_constructor_2(self):
         self.assertTrue(isinstance(self.db_filter.graph, Graph))
         self.assertTrue(isinstance(self.db_filter.engine, Engine))

     def test_values_1(self):
         self.db_filter.values = ["Test1", "Test2"]
         self.assertFalse(self.db_filter.values_valid)

     def test_values_2(self):
         with self.assertRaises(TypeError):
             self.db_filter.values = "Hello"

     def test_add_1(self):
         self.db_filter.add("phage.PhageID=Trixie")
         self.assertTrue("phage.PhageID=" in self.db_filter.filters.keys())

     def test_add_2(self):
         self.db_filter.add("phage.PhageID=Trixie")
         where_clauses = self.db_filter.filters["phage.PhageID="]
         self.assertTrue(isinstance(where_clauses, list))
         self.assertTrue(isinstance(where_clauses[0], BinaryExpression))

     def test_add_3(self):
         self.db_filter.add("phage.PhageID=Trixie")
         self.db_filter.add("phage.PhageID=D29")
         where_clauses = self.db_filter.filters["phage.PhageID="]
         self.assertEqual(len(where_clauses), 2)

     def test_remove_1(self):
         self.db_filter.add("phage.PhageID=Trixie")
         self.db_filter.remove("phage.PhageID=Trixie")
         self.assertEqual(self.db_filter.filters, {})

     def test_remove_2(self):
         self.db_filter.add("phage.PhageID=Trixie")
         self.db_filter.add("phage.PhageID=D29")
         self.db_filter.remove("phage.PhageID=Trixie")
         where_clauses = self.db_filter.filters["phage.PhageID="]
         self.assertTrue(len(where_clauses) == 1)
         self.assertEqual(where_clauses[0].right.value, "D29")

     def test_build_where_clauses_1(self):
         self.db_filter.add("phage.PhageID=Trixie")
         self.db_filter.add("phage.PhageID=D29")

         queries = self.db_filter.build_where_clauses()
         self.assertEqual(len(queries), 2)

     def test_build_values(self):
         phageid = self.alchemist.get_column("phage.PhageID")
         self.db_filter.key = phageid

         self.db_filter.values = ["Trixie", "D29"]
         values = self.db_filter.build_values()

         self.assertTrue("Trixie" in values)
         self.assertTrue("D29" in values)

     def test_refresh(self):
         phageid = self.alchemist.get_column("phage.PhageID")
         self.db_filter.key = phageid
         self.db_filter.values = ["Trixie", "D29", "Sheetz"]
         self.db_filter.refresh()

         self.assertTrue("Trixie" in self.db_filter.values)
         self.assertTrue("D29" in self.db_filter.values)
         self.assertFalse("Sheetz" in self.db_filter.values)

     def test_sort(self):
         phageid = self.alchemist.get_column("phage.PhageID")
         self.db_filter.key = phageid
         self.db_filter.values = ["Trixie", "D29"]
         self.db_filter.sort(phageid)

         self.assertTrue("Trixie" in self.db_filter.values)
         self.assertTrue("D29" in self.db_filter.values)
         self.assertEqual(self.db_filter.values[0], "D29")

     def test_random(self):
         pass

if __name__ == "__main__":
    unittest.main()
