"""Integration unittests for the filter module"""

from networkx import Graph
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import querying
from pdm_utils.functions.test_db import setup_test_db
from pdm_utils.functions.test_db import teardown_test_db
from sqlalchemy.engine.base import Engine
from sqlalchemy.sql.elements import BinaryExpression
from unittest.mock import Mock, patch
import unittest

class TestFilter(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        setup_test_db()

    def setUp(self):
        alchemist = AlchemyHandler()
        alchemist.username="pdm_anon"
        alchemist.password="pdm_anon"
        alchemist.database="test_db"
        alchemist.connect()
        alchemist.build_graph()
        self.alchemist = alchemist

        self.db_filter = Filter(alchemist=self.alchemist)

        phageid = self.alchemist.get_column("phage.PhageID") 
        self.phageid = phageid

    def test_constructor_1(self):
        db_filter = Filter(alchemist=self.alchemist)

    def test_constructor_2(self):
        self.assertTrue(isinstance(self.db_filter.graph, Graph))
        self.assertTrue(isinstance(self.db_filter.engine, Engine))

    def test_add_1(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.assertTrue("phage.PhageID=" in self.db_filter.filters.keys())

    def test_add_2(self):
        self.db_filter.add("phage.PhageID=Myrna")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertTrue(isinstance(where_clauses, list))
        self.assertTrue(isinstance(where_clauses[0], BinaryExpression))

    def test_add_3(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertEqual(len(where_clauses), 2)

    def test_remove_1(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.remove("phage.PhageID=Myrna")
        self.assertEqual(self.db_filter.filters, {})

    def test_remove_2(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")
        self.db_filter.remove("phage.PhageID=Myrna")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertTrue(len(where_clauses) == 1)
        self.assertEqual(where_clauses[0].right.value, "D29")

    def test_build_where_clauses_1(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")

        queries = self.db_filter.build_where_clauses()
 
        self.assertEqual(len(queries), 2)

    def test_build_where_clauses_2(self):
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")

        queries = self.db_filter.build_where_clauses()
 
        for query in queries:
            self.assertTrue(isinstance(query, BinaryExpression))

    def test_build_values_1(self):
        self.db_filter.key = self.phageid

        self.db_filter.values = ["Myrna", "D29"]
        values = self.db_filter.build_values()

        self.assertTrue("Myrna" in values)
        self.assertTrue("D29" in values)
  
    def test_transpose_1(self):
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.phageid

        self.db_filter.refresh()
        
        host_genera = self.db_filter.transpose("phage.HostGenus")
        
        self.assertEqual(host_genera, ["Mycobacterium"])
    
    def test_retrieve_1(self):
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.phageid
        
        self.db_filter.refresh()

        myrna_data = self.db_filter.retrieve(["phage.HostGenus",
                                              "phage.Cluster",
                                              "gene.Notes"])

        self.assertTrue(len(myrna_data) == 3 )
        self.assertTrue(isinstance(myrna_data, dict))

        self.assertEqual(myrna_data["HostGenus"], ["Mycobacterium"])
        self.assertEqual(myrna_data["Cluster"], ["C"])

    def test_refresh_2(self):
        self.db_filter.key = self.phageid
        self.db_filter.values = ["Myrna", "D29", "Sheetz"]
        self.db_filter.refresh()

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertTrue("D29" in self.db_filter.values)
        self.assertFalse("Sheetz" in self.db_filter.values)

    def test_update_1(self):
        self.db_filter.key = self.phageid
        self.db_filter.values = ["Myrna", "D29"]
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.update() 

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertFalse("D29" in self.db_filter.values)

    def test_sort_1(self):
        self.db_filter.key = self.phageid
        self.db_filter.values = ["Myrna", "D29"]
        self.db_filter.sort(self.phageid)

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertTrue("D29" in self.db_filter.values)
        self.assertEqual(self.db_filter.values[0], "D29")

    def test_group_1(self):
        self.db_filter.key = self.phageid
        self.db_filter.values = ["Myrna", "D29"]
        group_results = self.db_filter.group(self.phageid)

        self.assertTrue("Myrna" in group_results.keys())
        self.assertTrue("Myrna" in group_results["Myrna"])

        self.assertTrue("D29" in group_results.keys())
        self.assertTrue("D29" in group_results.keys())

    def test_group_2(self):
        self.db_filter.key = self.phageid
        self.db_filter.values = ["Myrna", "D29"]
        group_results = self.db_filter.group("phage.HostGenus")

        self.assertTrue("Mycobacterium" in group_results.keys())

        self.assertTrue("Myrna" in group_results["Mycobacterium"])
        self.assertTrue("D29" in group_results["Mycobacterium"]) 

    @classmethod
    def tearDownClass(self):
        teardown_test_db()

if __name__ == "__main__":
    unittest.main()
