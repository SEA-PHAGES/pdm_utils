"""Integration unittests for the filter module"""

from networkx import Graph
from pathlib import Path
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import querying
from sqlalchemy.engine.base import Engine
from sqlalchemy.sql.elements import BinaryExpression
from unittest.mock import Mock, patch
import sys
import unittest

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

class TestFilter(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    def setUp(self):
        alchemist = AlchemyHandler()
        alchemist.username="pdm_anon"
        alchemist.password="pdm_anon"
        alchemist.database="test_db"
        alchemist.connect()
        alchemist.build_graph()
        self.alchemist = alchemist

        self.db_filter = Filter(alchemist=self.alchemist)

        self.phage = self.alchemist.metadata.tables["phage"]
        self.gene = self.alchemist.metadata.tables["gene"] 
        self.trna = self.alchemist.metadata.tables["trna"]

        self.PhageID = self.phage.c.PhageID
        self.Cluster = self.phage.c.Cluster
        self.Subcluster = self.phage.c.Subcluster
        self.GeneID = self.gene.c.GeneID
        self.PhamID = self.gene.c.PhamID
        self.TrnaID = self.trna.c.TrnaID
        self.Product = self.trna.c.Product

    def test_add_1(self):
        """Verify that add() creates a dictionary key as expected.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.assertTrue("phage.PhageID=" in self.db_filter.filters.keys())

    def test_add_2(self):
        """Verify that add() stores BinaryExpression data as expected.
        """

        self.db_filter.add("phage.PhageID=Myrna")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertTrue(isinstance(where_clauses, list))
        self.assertTrue(isinstance(where_clauses[0], BinaryExpression))

    def test_add_3(self):
        """Verify that add() recognizes previous add() data.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertEqual(len(where_clauses), 2)

    def test_remove_1(self):
        """Verify that remove() removes dictionary entry after depleted.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.remove("phage.PhageID=Myrna")
        self.assertEqual(self.db_filter.filters, {})

    def test_remove_2(self):
        """Verify that remove() conserves dictionary entry if not depleted.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")
        self.db_filter.remove("phage.PhageID=Myrna")
        where_clauses = self.db_filter.filters["phage.PhageID="]
        self.assertTrue(len(where_clauses) == 1)
        self.assertEqual(where_clauses[0].right.value, "D29")

    def test_convert_column_input_1(self):
        """Verify that convert_column_input() converts string column input.
        """
        self.db_filter.key = self.Cluster

        column = self.db_filter.convert_column_input("phage.PhageID") 

        self.assertEqual(column, self.PhageID)

    def test_convert_column_input_2(self):
        """Verify that convert_column_input() conserves Column input.
        """
        self.db_filter.key = self.Cluster

        column = self.db_filter.convert_column_input(self.PhageID)

        self.assertEqual(column, self.PhageID)

    def test_convert_column_input_3(self):
        """Verify that convert_column_input() raises TypeError.
        convert_column_input() should raise TypeError when column input is
        neither a string or a Column.
        """
        self.db_filter.key = self.Cluster

        with self.assertRaises(TypeError):
            self.db_filter.convert_column_input(None)

    def test_build_where_clauses_1(self):
        """Verify that build_where_clauses() forms list of expected length.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")

        queries = self.db_filter.build_where_clauses()
 
        self.assertEqual(len(queries), 2)

    def test_build_where_clauses_2(self):
        """Verify that build_where_clauses() forms list of BinaryExpressions.
        """
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.add("phage.PhageID=D29")

        queries = self.db_filter.build_where_clauses()
 
        for query in queries:
            self.assertTrue(isinstance(query, BinaryExpression))

    def test_build_values_1(self):
        """Verify that build_values() does not exclude values as expected.
        """
        self.db_filter.key = self.PhageID

        values = self.db_filter.build_values()

        self.assertTrue("Myrna" in values)
        self.assertTrue("D29" in values)
        self.assertTrue("Alice" in values)
        self.assertTrue("Trixie" in values)

    def test_build_values_2(self):
        """Verify that build_values() utilizes WHERE clauses as expected.
        """
        self.db_filter.key = self.PhageID

        where_clause = (self.Cluster == "A")
        values = self.db_filter.build_values(where=where_clause)

        self.assertTrue("D29" in values)
        self.assertTrue("Trixie" in values)
        self.assertFalse("Myrna" in values)
        self.assertFalse("Alice" in values)

    def test_build_values_3(self):
        """Verify that build_values() creates DISTINCT values as expected.
        """
        self.db_filter.key = self.Cluster

        where_clause = (self.Subcluster == "A2")
        values = self.db_filter.build_values(where=where_clause)

        self.assertEqual(len(values), 1)
        self.assertEqual(values, ["A"])
  
    def test_transpose_1(self):
        """Verify that transpose() utilizes Filter values as expected.
        """
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.PhageID

        self.db_filter.refresh()
        
        clusters = self.db_filter.transpose("phage.Cluster")
        
        self.assertEqual(clusters, ["C"])

    def test_transpose_2(self):
        """Verify that transpose() can optionally create dict return value.
        """
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.PhageID

        self.db_filter.refresh()
        
        clusters_dict = self.db_filter.transpose(self.Cluster, return_dict=True)
        
        self.assertEqual(clusters_dict["Cluster"], ["C"])

    def test_transpose_3(self):
        """Verify that transpose() can alter Filter properties as expected.
        """
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.PhageID

        self.db_filter.refresh()

        self.db_filter.transpose("phage.Cluster", set_values=True)

        self.assertEqual(self.db_filter.key, self.Cluster)
        self.assertEqual(self.db_filter.values, ["C"])
        
    def test_mass_transpose_1(self):
        """Verify that mass_tranpose() returns DISTINCT values as expected.
        """
        self.db_filter.values = ["Myrna"]
        self.db_filter.key = self.PhageID
        
        self.db_filter.refresh()

        myrna_data = self.db_filter.mass_transpose(["phage.HostGenus",
                                              "phage.Cluster",
                                              "gene.Notes"])

        self.assertTrue(len(myrna_data) == 3 )
        self.assertTrue(isinstance(myrna_data, dict))

        self.assertEqual(myrna_data["HostGenus"], ["Mycobacterium"])
        self.assertEqual(myrna_data["Cluster"], ["C"])

    def test_mass_transpose_2(self):
        """Verify that mass_tranpose() utilizes all values as expected.
        """
        self.db_filter.values = ["Myrna", "Trixie"]
        self.db_filter.key = self.PhageID

        self.db_filter.refresh()

        data = self.db_filter.mass_transpose(["phage.HostGenus",
                                        "phage.Cluster",
                                        "gene.Notes"])

        self.assertTrue(len(data) == 3)
        self.assertTrue(isinstance(data, dict))

        self.assertEqual(data["HostGenus"], ["Mycobacterium"])
        self.assertEqual(data["Cluster"], ["C", "A"])

    def test_retrieve_1(self):
        """Verify that retrieve() separates data as expected.
        """
        self.db_filter.values = ["Myrna", "Trixie"]
        self.db_filter.key = self.PhageID

        self.db_filter.refresh()

        data = self.db_filter.retrieve(["phage.HostGenus",
                                        "phage.Cluster"])
        
        myrna_data = data["Myrna"]
        self.assertEqual(myrna_data["HostGenus"], ["Mycobacterium"])
        self.assertEqual(myrna_data["Cluster"], ["C"])

        trixie_data = data["Trixie"]
        self.assertEqual(trixie_data["HostGenus"], ["Mycobacterium"])
        self.assertEqual(trixie_data["Cluster"], ["A"])

    def test_retrieve_2(self):
        """Verify that retrieve() separates data as expected.
        """
        self.db_filter.values = ["A", "C"]
        self.db_filter.key = self.Cluster

        self.db_filter.refresh()

        data = self.db_filter.retrieve(["phage.Cluster", "phage.PhageID"])

        a_data = data["A"]
        self.assertEqual(a_data["Cluster"], ["A"])
        self.assertTrue("Trixie" in a_data["PhageID"])
        self.assertFalse("Myrna" in a_data["PhageID"])

        c_data = data["C"]
        self.assertEqual(c_data["Cluster"], ["C"])
        self.assertFalse("Trixie" in c_data["PhageID"])
        self.assertTrue("Myrna" in c_data["PhageID"])

    def test_refresh_1(self):
        """Verify that refresh() eliminates invalid data.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29", "Sheetz"]
        self.db_filter.refresh()

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertTrue("D29" in self.db_filter.values)
        self.assertFalse("Sheetz" in self.db_filter.values)

    def test_update_1(self):
        """Verify that update() filters out values.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29"]
        self.db_filter.add("phage.PhageID=Myrna")
        self.db_filter.update() 

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertFalse("D29" in self.db_filter.values)

    def test_sort_1(self):
        """Verify that sort() orders values as expected.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29"]
        self.db_filter.sort(self.PhageID)

        self.assertTrue("Myrna" in self.db_filter.values)
        self.assertTrue("D29" in self.db_filter.values)
        self.assertEqual(self.db_filter.values[0], "D29")

    def test_group_1(self):
        """Verify that group() creates separate groups as expected.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29"]
        group_results = self.db_filter.group(self.PhageID)

        self.assertTrue("Myrna" in group_results.keys())
        self.assertTrue("Myrna" in group_results["Myrna"])

        self.assertTrue("D29" in group_results.keys())
        self.assertTrue("D29" in group_results["D29"])

    def test_group_2(self):
        """Verify that group() recognizes similarities in values as expected.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29"]
        group_results = self.db_filter.group("phage.HostGenus")

        self.assertTrue("Mycobacterium" in group_results.keys())

        self.assertTrue("Myrna" in group_results["Mycobacterium"])
        self.assertTrue("D29" in group_results["Mycobacterium"]) 

    def test_group_3(self):
        """Verify that group() recognizes differences in values as expected.
        """
        self.db_filter.key = self.PhageID
        self.db_filter.values = ["Myrna", "D29", "Trixie"]
        group_results = self.db_filter.group("phage.Cluster")

        self.assertTrue("A" in group_results.keys())
        self.assertTrue("C" in group_results.keys())

        self.assertTrue("Myrna" in group_results["C"])
        self.assertTrue("D29" in group_results["A"])
        self.assertTrue("Trixie" in group_results["A"])

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

if __name__ == "__main__":
    unittest.main()
