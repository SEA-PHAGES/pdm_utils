from datetime import datetime
from decimal import Decimal
from networkx import Graph
from networkx import shortest_path
from pathlib import Path
from pdm_utils.functions import querying
from pdm_utils.functions.mysqldb import query_dict_list
from sqlalchemy import Column
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import select
from sqlalchemy.exc import InternalError
from sqlalchemy.sql import func
from sqlalchemy.sql import functions
from sqlalchemy.sql.elements import BinaryExpression
from sqlalchemy.sql.elements import Grouping
from sqlalchemy.sql.elements import UnaryExpression
from unittest.mock import Mock
from unittest.mock import patch
import re
import sys
import unittest

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

class TestQuerying(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    def setUp(self):
        self.engine = create_engine(test_db_utils.create_engine_string())

        self.metadata = MetaData(bind=self.engine)
        self.metadata.reflect()

        self.graph = querying.build_graph(self.metadata)
   
        self.phage = self.metadata.tables["phage"]
        self.gene = self.metadata.tables["gene"]
        self.trna = self.metadata.tables["trna"]

        self.PhageID = self.phage.c.PhageID
        self.Cluster = self.phage.c.Cluster
        self.Subcluster = self.phage.c.Subcluster
        self.GeneID = self.gene.c.GeneID
        self.PhamID = self.gene.c.PhamID
        self.TrnaID = self.trna.c.TrnaID
        self.Product = self.trna.c.Product

    def test_get_table_1(self):
        """Verify get_table() retrieves correct Table.
        """
        self.assertEqual(querying.get_table(self.metadata, "phage"), self.phage)

    def test_get_table_2(self):
        """Verify get_table() operates case insensitive.
        """
        self.assertEqual(querying.get_table(self.metadata, "pHAgE"), self.phage)

    def test_get_table_3(self):
        """Verify get_table() raises ValueError from invalid table name.
        """
        with self.assertRaises(ValueError):
            querying.get_table(self.metadata, "not_a_table")

    def test_get_column_1(self):
        """Verify get_column() retrieves correct Column.
        """
        self.assertEqual(querying.get_column(self.metadata, "gene.GeneID"),
                                                             self.GeneID) 

    def test_get_column_2(self):
        """Verify get_column() retrieves correct Column.
        """
        self.assertEqual(querying.get_column(self.metadata, "GENE.GENEID"), 
                                                             self.GeneID)

    def test_get_column_3(self):
        """Verify get_column() raises ValueError from invalid column name.
        """
        with self.assertRaises(ValueError):
            querying.get_column(self.metadata, "gene.not_a_column")

    def test_build_graph_1(self):
        """Verify build_graph() returns a Graph object.
        """
        self.assertTrue(isinstance(self.graph, Graph))

    def test_build_graph_2(self):
        """Verify build_graph() correctly structures Graph nodes.
        """
        self.assertTrue("phage" in self.graph.nodes)
        self.assertTrue("gene" in self.graph.nodes)
        self.assertTrue("trna" in self.graph.nodes)

    def test_build_graph_3(self):
        """Verify build_graph correctly stores Table objects.
        """
        phage_node = self.graph.nodes["phage"]
        gene_node = self.graph.nodes["gene"]
        trna_node = self.graph.nodes["trna"]

        self.assertEqual(phage_node["table"], self.phage)
        self.assertEqual(gene_node["table"], self.gene)
        self.assertEqual(trna_node["table"], self.trna)

    def test_build_graph_4(self):
        """Verify build_graph() correctly stores its MetaData object.
        """
        self.assertEqual(self.graph.graph["metadata"], self.metadata)

    def test_build_where_clause_1(self):
        """Verify build_where_clause() returns a BinaryExpression object.
        """
        where_clause = querying.build_where_clause(self.graph, 
                                                        "phage.PhageID=Trixie") 
        
        self.assertTrue(isinstance(where_clause, BinaryExpression))

    def test_build_where_clause_2(self):
        """Verify build_where_clause() builds from specified Column.
        """
        where_clause = querying.build_where_clause(self.graph, 
                                                        "phage.PhageID=Trixie")

        self.assertEqual(where_clause.left, self.PhageID)

    def test_build_where_clause_3(self):
        """Verify build_where_clause() builds from specified value.
        """
        where_clause = querying.build_where_clause(self.graph,
                                                        "phage.PhageID=Trixie")

        self.assertEqual(where_clause.right.value, "Trixie")

    def test_build_onclause_1(self):
        """Verify build_onclause() returns a BinaryExpression object.
        """
        onclause = querying.build_onclause(self.graph, "gene", "phage")

        self.assertTrue(isinstance(onclause, BinaryExpression))

    def test_build_onclause_2(self):
        """Verify build_onclause() builds onclause with conserved direction.
        build_onclause() should build from primary key to foreign key.
        """
        onclause = querying.build_onclause(self.graph, "gene", "phage")

        self.assertEqual(onclause.left.table, self.phage)
        self.assertEqual(onclause.right.table, self.gene)

    def test_build_onclause_3(self):
        """Verify build_onclause() builds onclause with conserved direction.
        build_onclause() should build form primary key to foreign key.
        """
        onclause = querying.build_onclause(self.graph, "phage", "gene")

        self.assertEqual(onclause.left.table, self.phage)
        self.assertEqual(onclause.right.table, self.gene)

    def test_build_onclause_4(self):
        """Verify build_onclause() raises KeyError when Tables are not linked.
        """
        with self.assertRaises(KeyError):
            querying.build_onclause(self.graph, "gene", "trna")

    def test_get_table_list_1(self):
        """Verify get_table_list() returns correct ordered Tables.
        """
        columns = [self.GeneID, self.PhageID, self.TrnaID]

        table_list = querying.get_table_list(columns)

        self.assertEqual(table_list[0], self.gene)
        self.assertEqual(table_list[1], self.phage)
        self.assertEqual(table_list[2], self.trna)

    def test_get_table_list_2(self):
        """Verify get_table_list() accepts Column input.
        """
        table_list = querying.get_table_list(self.Product)

        self.assertEqual(table_list[0], self.trna)

    def test_get_table_pathing_1(self):
        """Verify get_table_pathing() builds from common center table.
        """
        table_list = [self.phage, self.gene, self.trna]

        table_pathing = querying.get_table_pathing(self.graph, table_list)

        self.assertEqual(table_pathing[0], self.phage)
        self.assertEqual(table_pathing[1][0], ["phage", "gene"])
        self.assertEqual(table_pathing[1][1], ["phage", "trna"])

    def test_get_table_pathing_2(self):
        """Verify get_table_pathing() builds from distal center table.
        """
        table_list = [self.gene, self.phage, self.trna]

        table_pathing = querying.get_table_pathing(self.graph, table_list)

        self.assertEqual(table_pathing[0], self.gene)
        self.assertEqual(table_pathing[1][0], ["gene", "phage"])
        self.assertEqual(table_pathing[1][1], ["gene", "phage", "trna"])

    def test_get_table_pathing_3(self):
        """Verify get_table_pathing() correctly utilizes center_table input.
        """
        table_list = [self.gene, self.phage, self.trna]

        table_pathing = querying.get_table_pathing(self.graph, table_list,
                                                   center_table=self.phage)

        self.assertEqual(table_pathing[0], self.phage)
        self.assertEqual(table_pathing[1][0], ["phage", "gene"])
        self.assertEqual(table_pathing[1][1], ["phage", "trna"])

    def test_join_pathed_tables_1(self):
        """Verify join_pathed_tables() joins Tables in specified order.
        """
        table_pathing = [self.phage, [["phage", "gene"],
                                      ["phage", "trna"]]]

        joined_tables = querying.join_pathed_tables(self.graph, table_pathing)

        last_table = joined_tables.right
        joined_tables = joined_tables.left
        self.assertEqual(last_table, self.trna)

        second_table = joined_tables.right
        first_table = joined_tables.left
        self.assertEqual(second_table, self.gene)

        self.assertEqual(first_table, self.phage)

    def test_join_pathed_tables_2(self):
        """Verify join_pathed_tables() recognizes already pathed Tables.
        """
        table_pathing = [self.gene, [["gene", "phage"],
                                     ["gene", "phage", "trna"]]]

        joined_tables = querying.join_pathed_tables(self.graph, table_pathing)

        last_table = joined_tables.right
        joined_tables = joined_tables.left
        self.assertEqual(last_table, self.trna)

        second_table = joined_tables.right
        first_table = joined_tables.left
        self.assertEqual(second_table, self.phage)

        self.assertEqual(first_table, self.gene)

    def test_build_select_1(self):
        """Verify build_select() creates valid SQLAlchemy executable.
        """
        select_query = querying.build_select(self.graph, self.PhageID)

        phage_ids = []
        dict_list = query_dict_list(self.engine, select_query)
        for dict in dict_list:
            phage_ids.append(dict["PhageID"]) 

        self.assertTrue("Myrna" in phage_ids)
        self.assertTrue("D29" in phage_ids)
        self.assertTrue("Trixie" in phage_ids)

    def test_build_select_2(self):
        """Verify build_select() appends WHERE clauses to executable.
        """
        where_clause = (self.Cluster == "A")
        select_query = querying.build_select(self.graph, self.PhageID,
                                                        where=where_clause)
        
        phage_ids = []
        dict_list = query_dict_list(self.engine, select_query)
        for dict in dict_list:
            phage_ids.append(dict["PhageID"])

        self.assertTrue("Trixie" in phage_ids)
        self.assertTrue("D29" in phage_ids)
        self.assertFalse("Myrna" in phage_ids)

    def test_build_select_3(self):
        """Verify build_select() appends ORDER BY clauses to executable.
        """
        select_query = querying.build_select(self.graph, self.PhageID,
                                                        order_by=self.PhageID)

        dict_list = query_dict_list(self.engine, select_query)

        phage_ids = []
        dict_list = query_dict_list(self.engine, select_query)
        for dict in dict_list:
            phage_ids.append(dict["PhageID"]) 

        self.assertEqual("Alice", phage_ids[0])
        self.assertTrue("Myrna" in phage_ids)
        self.assertTrue("D29" in phage_ids)
        self.assertTrue("Trixie" in phage_ids)

    def test_build_select_4(self):
        """Verify build_select() handles many-to-one relations as expected.
        build_select() queries should duplicate 'one' when filtering 'many'
        """
        where_clause = (self.Subcluster == "A2")
        select_query = querying.build_select(self.graph, self.Cluster,
                                                        where=where_clause)

        dict_list = query_dict_list(self.engine, select_query)

        self.assertTrue(len(dict_list) > 1)

    def test_build_count_1(self):
        """Verify build_count() creates valid SQLAlchemy executable.
        """
        count_query = querying.build_count(self.graph, self.PhageID)

        dict_list = query_dict_list(self.engine, count_query)
        count_dict = dict_list[0] 

        self.assertTrue(isinstance(count_dict["count_1"], int))

    def test_build_count_2(self):
        """Verify build_count() appends WHERE clauses to executable.
        """
        where_clause = (self.PhageID == "Trixie")
        count_query = querying.build_count(self.graph, self.PhageID,
                                                        where=where_clause)

        dict_list = query_dict_list(self.engine, count_query)
        count_dict = dict_list[0]

        self.assertEqual(count_dict["count_1"], 1)

    def test_build_count_3(self):
        """Verify build_count() recognizes multiple inputs as expected.
        """
        where_clause = (self.Cluster == "A")
        count_query = querying.build_count(self.graph, 
                                                [self.PhageID, 
                                                 self.Cluster.distinct()],
                                                        where=where_clause)

        dict_list = query_dict_list(self.engine, count_query)
        count_dict = dict_list[0]

        self.assertTrue(count_dict["count_1"] > 1)
        self.assertEqual(count_dict["count_2"], 1)

    def test_build_distinct_1(self):
        """Verify build_distinct() handles many-to-one relations as expected.
        build_distinct() should not duplicate 'one' when handling 'many.'
        """
        where_clause = (self.Subcluster == "A2")
        distinct_query = querying.build_distinct(self.graph, self.Cluster,
                                                        where=where_clause)

        dict_list = query_dict_list(self.engine, distinct_query)
        self.assertEqual(len(dict_list), 1)

        distinct_dict = dict_list[0] 
        self.assertEqual(distinct_dict["Cluster"], "A")

    def test_build_distinct_2(self):
        """Verify build_distinct() cannot variable aggregated columns.
        MySQL does not accept DISTINCT queries with aggregated
        and non-aggregated columns.
        """
        distinct_query = querying.build_distinct(self.graph, 
                                                    [self.PhageID,
                                                     func.count(self.Cluster)])

 
        with self.assertRaises(InternalError):
            dict_list = query_dict_list(self.engine, distinct_query)


    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

if __name__ == "__main__":
    unittest.main()
