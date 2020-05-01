import unittest
from unittest.mock import Mock
from unittest.mock import MagicMock
from unittest.mock import patch
from unittest.mock import PropertyMock

from networkx import Graph
from sqlalchemy import create_engine
from sqlalchemy import Column
from sqlalchemy import MetaData
from sqlalchemy import Table
from sqlalchemy.sql import distinct
from sqlalchemy.sql import func
from sqlalchemy.sql import functions
from sqlalchemy.sql.elements import BinaryExpression
from sqlalchemy.sql.elements import BindParameter
from sqlalchemy.sql.elements import UnaryExpression
from sqlalchemy.sql.schema import ForeignKey

from pdm_utils.functions import querying

class TestUseMetadata(unittest.TestCase):
    def setUp(self):
        self.metadata = Mock(spec=MetaData)

        self.phage = Mock(spec=Table) 
        self.gene  = Mock(spec=Table) 
        self.trna  = Mock(spec=Table)

        type(self.phage).name = PropertyMock(return_value="phage")
        type(self.gene).name  = PropertyMock(return_value="gene")
        type(self.trna).name  = PropertyMock(return_value="trna")

        tables = {"phage" : self.phage,
                  "gene"  : self.gene,
                  "trna"  : self.trna}

        self.PhageID = Mock(spec=Column)
        self.GeneID  = Mock(spec=Column)

        type(self.phage).columns = {"PhageID" : self.PhageID}
        type(self.gene).columns  = {"PhageID" : self.PhageID,
                                    "GeneID"  : self.GeneID}
        type(self.trna).columns  = {"PhageID" : self.PhageID}

        type(self.PhageID).table = PropertyMock(return_value=self.phage)
        type(self.GeneID).table  = PropertyMock(return_value=self.gene)

        self.PhageGeneKey = Mock(spec=ForeignKey)
        self.PhageTrnaKey = Mock(spec=ForeignKey)

        type(self.PhageGeneKey).column = PropertyMock(
                                        return_value=self.PhageID)
        type(self.PhageGeneKey).parent = PropertyMock(
                                        return_value=self.PhageID)
        type(self.PhageTrnaKey).column = PropertyMock(
                                        return_value=self.PhageID)
        type(self.PhageTrnaKey).parent = PropertyMock(
                                        return_value=self.PhageID)

        type(self.phage).foreign_keys = PropertyMock(
                                        return_value=set())
        type(self.gene).foreign_keys = PropertyMock(
                                        return_value=set([self.PhageGeneKey]))
        type(self.trna).foreign_keys = PropertyMock(
                                        return_value=set([self.PhageTrnaKey]))

        type(self.metadata).tables = PropertyMock(return_value=tables)

    @patch("pdm_utils.functions.querying.parsing.translate_table")
    def test_get_table_1(self, translate_table_mock):
        """Verify translate_table() is called with the correct parameters.
        """
        translate_table_mock.return_value = "phage"

        querying.get_table(self.metadata, "phage") 

        translate_table_mock.assert_called_with(self.metadata, "phage")

    def test_get_table_2(self):
        """Verify get_table() returns the correct Table object.
        """
        table_obj = querying.get_table(self.metadata, "phage")

        self.assertEqual(table_obj, self.phage)
   
    @patch("pdm_utils.functions.querying.parsing.parse_column")
    def test_get_column_1(self, ParseColumn):
        """Verify parse_column() is called with the correct parameters.
        """
        ParseColumn.return_value = ["phage", "PhageID"]

        querying.get_column(self.metadata, "phage.PhageID")

        ParseColumn.assert_called_with("phage.PhageID")

    @patch("pdm_utils.functions.querying.get_table")
    def test_get_column_2(self, GetTable):
        """Verify get_table() is called with the correct parameters.
        """
        GetTable.return_value = self.phage

        querying.get_column(self.metadata, "phage.PhageID")

        GetTable.assert_called_with(self.metadata, "phage")

    @patch("pdm_utils.functions.querying.parsing.translate_column")
    def test_get_column_3(self, TranslateColumn):
        """Verify translate_column() is called with the correct parameters.
        """
        TranslateColumn.return_value = "PhageID"

        querying.get_column(self.metadata, "phage.PhageID")

        TranslateColumn.assert_called_with(self.metadata, "phage.PhageID")

    def test_get_column_4(self):
        """Verify get_column() returns the correct Column object.
        """
        column_obj = querying.get_column(self.metadata, "phage.PhageID")

        self.assertEqual(column_obj, self.PhageID)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_get_table_list_1(self, extract_column_mock):  
        """Verify extract_column() is called with the correct parameters. 
        extract_column() should take Column objects as well as lists.
        """
        extract_column_mock.return_value = self.PhageID

        querying.get_table_list([self.PhageID])

        extract_column_mock.assert_called_with(self.PhageID)
    
    @patch("pdm_utils.functions.querying.extract_column")
    def test_get_table_list_2(self, extract_column_mock):
        """Verify extract_column() is called with the correct parameters.
        extract_column() should take Column objects as well as lists.
        """
        extract_column_mock.return_value = self.PhageID

        querying.get_table_list(self.PhageID)
    
        extract_column_mock.assert_called_with(self.PhageID)
    
    def test_get_table_list_3(self):
        """Verify get_table_list() retrieves specified Table objects.
        get_table_list() should exclude unspecified Tables.
        """
        table_list = querying.get_table_list(self.PhageID)

        self.assertTrue(self.phage in table_list)
        self.assertFalse(self.gene in table_list)

    def test_get_table_list_4(self):
        """Verify get_table_list() retrieves specified Table objects.
        """
        columns = [self.PhageID, self.GeneID]
        table_list = querying.get_table_list(columns)

        self.assertTrue(self.phage in table_list)
        self.assertTrue(self.gene in table_list)

    def test_build_graph_1(self):
        """Verify build_graph() creates a NetworkX Graph object.
        """
        graph = querying.build_graph(self.metadata)

        self.assertTrue(isinstance(graph, Graph))
 
    def test_build_graph_2(self):
        """Verify build_graph() stores Table objects with correct structure.
        """
        graph = querying.build_graph(self.metadata)

        self.assertTrue("phage" in graph.nodes)
        self.assertTrue("gene"  in graph.nodes)
        self.assertTrue("trna"  in graph.nodes)

        self.assertEqual(graph.nodes["phage"]["table"], self.phage)
        self.assertEqual(graph.nodes["gene"]["table"],  self.gene)
        self.assertEqual(graph.nodes["trna"]["table"],  self.trna)

    def test_build_graph_3(self):
        """Verify build_graph() stores foreign keys with correct structure.
        """
        graph = querying.build_graph(self.metadata)
        
        phage = graph.nodes["phage"]["table"]
        gene = graph.nodes["gene"]["table"]
        trna = graph.nodes["trna"]["table"]

        phage_gene_edge = graph["gene"]["phage"]
        phage_trna_edge = graph["trna"]["phage"]

        self.assertEqual(phage_gene_edge["key"], self.PhageGeneKey)
        self.assertEqual(phage_trna_edge["key"], self.PhageTrnaKey)
   
    @patch("pdm_utils.functions.querying.parsing.translate_table")
    def test_build_onclause_1(self, translate_table_mock):
        """Verify translate_table() is called with the correct parameters.
        """
        translate_table_mock.return_value = "phage"
        graph = querying.build_graph(self.metadata)
        
        table = querying.build_onclause(graph, "phage", "phage")

        translate_table_mock.assert_called_with(self.metadata, "phage")

    def test_build_onclause_2(self):
        """Verify build_onclause() detects table redundencies
        """
        graph = querying.build_graph(self.metadata) 
       
        table = querying.build_onclause(graph, "phage", "phage")

        self.assertEqual(table, self.phage)

class TestExtract(unittest.TestCase):
    def setUp(self):
        self.column = Mock(spec=Column)
        self.columns = [self.column]
        self.check_type = Mock()
    
    @patch("pdm_utils.functions.querying.isinstance")
    def test_extract_column_1(self, is_instance_mock):
        """Verify extract_column() raises a TypeError from innate type checks.
        extract_column() should throw an exception when the object is not
        any expected Column-related SQLAlchemy object.
        """
        with self.assertRaises(TypeError):
            is_instance_mock.return_value = False
            querying.extract_column(self.column)

        is_instance_mock.assert_any_call(self.column, Column)
        is_instance_mock.assert_any_call(self.column, functions.count)
        is_instance_mock.assert_any_call(self.column, UnaryExpression)
        is_instance_mock.assert_any_call(self.column, BinaryExpression)

    def test_extract_column_2(self):
        """Verify extract_column() returns received Column object.
        """
        self.assertEqual(querying.extract_column(self.column), self.column)
        
    def test_extract_column_3(self):
        """Verify extract_column() converts received func.count object.
        """
        column_list = PropertyMock(return_value=self.columns)
        clause_list = Mock()
        count = Mock(spec=functions.count)
        
        type(count).clauses = clause_list
        type(clause_list).clauses = column_list

        self.assertEqual(querying.extract_column(count), self.column)

    def test_extract_column_4(self):
        """Verify extract_column() converts received UnaryExpression object.
        """
        unary_expression = Mock(spec=UnaryExpression)

        type(unary_expression).element = self.column

        self.assertEqual(querying.extract_column(unary_expression), 
                                                 self.column)

    def test_extract_column_5(self):
        """Verify extract_column() converts received BinaryExpression object.
        """
        binary_expression = Mock(spec=BinaryExpression)

        type(binary_expression).left = self.column

        self.assertEqual(querying.extract_column(binary_expression), 
                                                 self.column)

    def test_extract_column_6(self): 
        """Verify extract_column() converts received BinaryExpression object.
        extract_column() should detect embedded UnaryExpression objects.
        """
        unary_expression = Mock(spec=UnaryExpression)
        binary_expression = Mock(spec=BinaryExpression)

        type(unary_expression).element = self.column
        type(binary_expression).left = unary_expression

        self.assertEqual(querying.extract_column(binary_expression), 
                                                 self.column)
    
    @patch("pdm_utils.functions.querying.isinstance")
    def test_extract_column_7(self, is_instance_mock):
        """Verify isinstance() is called with correct parameters.
        """
        is_instance_mock.return_value = True

        querying.extract_column(self.column, check=self.check_type)

        is_instance_mock.assert_any_call(self.column, self.check_type)

    def test_extract_column_8(self):
        """Verify extract_column() raises TypeError from imposed type check.
        """
        with self.assertRaises(TypeError):
            querying.extract_column(self.column, check=str)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_extract_columns_1(self, extract_column_mock):
        """Verify extract_column() is called with correct parameters.
        extract_columns() should accept both Column and list parameters.
        """
        querying.extract_columns(self.columns)

        extract_column_mock.assert_called_with(self.column, check=None)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_extract_columns_2(self, extract_column_mock): 
        """Verify extract_column() is called with correct parameters. 
        extract_columns() should accept both Column and list parameters.
        """
        querying.extract_columns(self.column)

        extract_column_mock.assert_called_with(self.column, check=None)

class TestUseGraph(unittest.TestCase): 
    def setUp(self):
        self.graph = Mock()
        self.center = Mock()
        self.table_object = Mock()

        self.table_1 = Mock()
        self.table_2 = Mock()
        self.table_3 = Mock()
        self.table_4 = Mock()

        type(self.center) .name = PropertyMock(return_value="center")
        type(self.table_1).name = PropertyMock(return_value="table_1")
        type(self.table_2).name = PropertyMock(return_value="table_2")
        type(self.table_3).name = PropertyMock(return_value="table_3")
        type(self.table_4).name = PropertyMock(return_value="table_4")

        self.table_list_1 = [self.center, 
                             self.table_1, 
                             self.table_2, 
                             self.table_3, 
                             self.table_4]

        self.table_names_1 = ["table_1",
                              "table_2",
                              "table_3",
                              "table_4"]

        self.table_list_2 = [self.table_1,
                             self.table_3,
                             self.table_4,
                             self.center,
                             self.table_2]

        self.table_names_2 = ["table_1",
                              "table_3",
                              "table_4",
                              "table_2"]
      
        nodes_dict = {"table_1" : {"table" : self.table_1},
                      "table_2" : {"table" : self.table_2},
                      "table_3" : {"table" : self.table_3},
                      "table_4" : {"table" : self.table_4}}
        mock_nodes = PropertyMock(return_value=nodes_dict)
     
        type(self.graph).nodes = mock_nodes
        type(self.center).name = PropertyMock(return_value="center")

        self.table_path_1 = ["center", "table_1"]
        self.table_path_2 = ["center", "table_2"]
        self.table_path_3 = ["center", "table_2", "table_3"]
        self.table_path_4 = ["center", "table_1", "table_4"]
        
        self.table_paths_1 = [self.table_path_1,
                              self.table_path_2,
                              self.table_path_3,
                              self.table_path_4]
        
        self.table_paths_2 = [self.table_path_1,
                              self.table_path_3,
                              self.table_path_4,
                              self.table_path_2]

        self.table_pathing_1 = [self.center, [self.table_path_1, 
                                              self.table_path_2, 
                                              self.table_path_3,
                                              self.table_path_4]]

        self.table_pathing_2 = [self.center, []]  

    @patch("pdm_utils.functions.querying.shortest_path")
    def test_get_table_pathing_1(self, shortest_path_mock):
        """Verify shortest_path() is called with correct parameters.
        """
        shortest_path_mock.return_value = self.table_path_1

        for index in range(len(self.table_names_1)):
            with self.subTest(table_index=index):
                querying.get_table_pathing(self.graph, self.table_list_1)

                shortest_path_mock.assert_any_call(self.graph, 
                                                "center",
                                                self.table_names_1[index])

                shortest_path_mock.clear()

    @patch("pdm_utils.functions.querying.shortest_path")
    def test_get_table_pathing_2(self, shortest_path_mock):
        """Verify center_table parameter influences shortest_path() parameters.
        """
        shortest_path_mock.return_value = self.table_path_1

        for index in range(len(self.table_names_2)):
            with self.subTest(table_index=index):
                querying.get_table_pathing(self.graph, self.table_list_2,
                                                center_table=self.center)

                shortest_path_mock.assert_any_call(self.graph,
                                                "center",
                                                self.table_names_2[index])

                shortest_path_mock.clear()

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_1(self, join_mock, build_on_clause_mock):
        """Verify build_onclause() influences join() parameters.
        """
        join_mock.return_value = self.center
        build_on_clause_mock.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        build_on_clause_mock.assert_any_call(self.graph, "center", "table_1")
        join_mock.assert_any_call(self.center, self.table_1,
                                          isouter=True,
                                          onclause="onclause")
    
    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_2(self, join_mock, build_on_clause_mock):
        """Verify join_pathed_tables() recognizes different Table connections.
        join_pathed_tables() should recognized multiple connections to center.
        """
        join_mock.return_value = self.center
        build_on_clause_mock.return_value = "onclause"
        
        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        build_on_clause_mock.assert_any_call(self.graph, "center", "table_2")
        join_mock.assert_any_call(self.center, self.table_2,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_3(self, join_mock, build_on_clause_mock):       
        """Verify join_pathed_tables() recognizes already pathed Tables.
        """
        join_mock.return_value = self.center
        build_on_clause_mock.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        build_on_clause_mock.assert_any_call(self.graph, "table_2", "table_3")
        join_mock.assert_any_call(self.center, self.table_3,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_4(self, join_mock, build_on_clause_mock):  
        """Verify join_pathed_tables() recognizes different Table connections.
        join_pathed_tables() should recognize connections to non-center Tables.
        """
        join_mock.return_value = self.center
        build_on_clause_mock.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        build_on_clause_mock.assert_any_call(self.graph, "table_1", "table_4")
        join_mock.assert_any_call(self.center, self.table_4,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_5(self, join_mock, build_on_clause_mock):
        """Verify join_pathed_tables() recognizes empty pathing.
        """
        join_mock.return_value = self.center
        build_on_clause_mock.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_2)

        build_on_clause_mock.assert_not_called()
        join_mock.assert_not_called()

class TestBuildClauses(unittest.TestCase):
    def setUp(self):
        self.graph = Mock()
        self.metadata = Mock()
        self.graph_properties = {"metadata" : self.metadata}
        type(self.graph).graph = PropertyMock(
                                 return_value=self.graph_properties)

        self.phage = Mock(spec=Table)
        self.gene  = Mock(spec=Table)
        self.trna  = Mock(spec=Table)
     
        self.tables = [self.phage, self.gene, self.trna]

        type(self.phage).name = PropertyMock(return_value="phage")
        type(self.gene).name  = PropertyMock(return_value="gene")
        type(self.trna).name  = PropertyMock(return_value="trna")
  
        self.table_names = ["phage", "gene", "trna"]

        nodes_dict = {"phage" : {"table" : self.phage},
                      "gene"  : {"table"  : self.gene},
                      "trna"  : {"table"  : self.trna}}
        tables_dict = {"phage" : self.phage,
                       "gene"  : self.gene,
                       "trna"  : self.trna}
    
        mock_nodes = PropertyMock(return_value=nodes_dict)
        mock_tables = PropertyMock(return_value=tables_dict)
        type(self.graph).nodes = mock_nodes
        type(self.metadata).tables = mock_tables

        self.pathing = [self.phage, [["phage", "gene"], ["phage", "trna"]]]

        self.Cluster = Mock(spec=Column)
        self.PhamID = Mock(spec=Column)
        self.Notes = Mock(spec=Column)

        self.columns = [self.Cluster, self.PhamID, self.Notes]

        type(self.Cluster).name = PropertyMock(return_value="Cluster")
        type(self.PhamID).name = PropertyMock(return_value="PhamID")
        type(self.Notes).name = PropertyMock(return_value="Notes")

        self.Cluster_type = Mock()
        self.PhamID_type = Mock()
        self.Notes_type = Mock()

        type(self.Cluster_type).python_type = PropertyMock(return_value=str)
        type(self.PhamID_type).python_type = PropertyMock(return_value=int)
        type(self.Notes_type).python_type = PropertyMock(return_value=bytes)

        type(self.Cluster).type = PropertyMock(return_value=self.Cluster_type)
        type(self.PhamID).type = PropertyMock(return_value=self.PhamID_type)
        type(self.Notes).type = PropertyMock(return_value=self.Notes_type)

        self.phage_columns = {"Cluster" : self.Cluster}
        self.gene_columns = {"PhamID" : self.PhamID}
        self.trna_columns = {"Notes" : self.Notes}

        type(self.phage).columns = PropertyMock(return_value=self.phage_columns)
        type(self.gene).columns  = PropertyMock(return_value=self.gene_columns)
        type(self.trna).columns  = PropertyMock(return_value=self.trna_columns)

        self.not_column = Mock()
        self.count_column = Mock()

        self.whereclause_1 = Mock(spec=BinaryExpression)
        self.whereclause_2 = Mock(spec=BinaryExpression)
        self.whereclause_3 = Mock(spec=BinaryExpression)

        self.not_whereclause = Mock()

        self.whereclauses = [self.whereclause_1,
                             self.whereclause_2,
                             self.whereclause_3]
 
    @patch("pdm_utils.functions.querying.parsing.parse_filter")
    def test_build_where_clause_1(self, parse_filter_mock): 
        """Verify parse_filter() is called with correct parameters.
        """
        parse_filter_mock.return_value = ["phage", "Cluster", "!=", "A"]

        querying.build_where_clause(self.graph, "phage.Cluster != A")

        parse_filter_mock.assert_called_with("phage.Cluster != A")

    @patch("pdm_utils.functions.querying.parsing.check_operator")
    def test_build_where_clause_2(self, check_operator_mock):
        """Verify check_operator() is called with correct parameters.
        """
        querying.build_where_clause(self.graph, "gene.PhamID = 2")

        check_operator_mock.assert_called_with("=", self.PhamID)
   
    @patch("pdm_utils.functions.querying.join_pathed_tables")
    @patch("pdm_utils.functions.querying.get_table_pathing")
    @patch("pdm_utils.functions.querying.get_table_list")
    def test_build_fromclause_1(self, get_table_list_mock, 
                                      get_table_pathing_mock, 
                                      join_pathed_tables_mock):
        """Verify function structure of build_fromclause().
        """
        get_table_list_mock.return_value = self.table_names
        get_table_pathing_mock.return_value = self.pathing
        join_pathed_tables_mock.return_value  = Mock()
        
        querying.build_fromclause(self.graph, self.tables)

        get_table_list_mock.assert_called_with(self.tables)
        get_table_pathing_mock.assert_called_with(self.graph, self.table_names)
        join_pathed_tables_mock.assert_called_with(self.graph, self.pathing)

    @patch("pdm_utils.functions.querying.append_order_by_clauses")
    @patch("pdm_utils.functions.querying.append_where_clauses")
    @patch("pdm_utils.functions.querying.select")
    @patch("pdm_utils.functions.querying.build_fromclause")
    @patch("pdm_utils.functions.querying.extract_columns")
    def test_build_select_1(self, extract_columns_mock,
                                  build_from_clause_mock, select_mock,
                                  append_where_clauses_mock, append_order_by_clauses_mock):
        """Verify function structure of build_select().
        """
        executable_mock = Mock()
        select_from_mock = Mock()
        type(executable_mock).select_from = select_from_mock
        select_from_mock.return_value = executable_mock

        extract_columns_mock.return_value = self.columns
        build_from_clause_mock.return_value = self.phage
        select_mock.return_value = executable_mock
        append_where_clauses_mock.return_value = executable_mock
        append_order_by_clauses_mock.return_value = executable_mock

        querying.build_select(self.graph, self.columns, 
                                            order_by=self.columns,
                                            add_in=self.columns)

        extract_columns_mock.assert_any_call(None)
        extract_columns_mock.assert_any_call(self.columns, check=Column)
        total_columns = self.columns * 4
        build_from_clause_mock.assert_called_once_with(self.graph, total_columns)
        
        select_mock.assert_called_once_with(self.columns)
        select_from_mock.assert_called_once_with(self.phage)
        append_where_clauses_mock.assert_called_once_with(executable_mock, None)
        append_order_by_clauses_mock.assert_called_once_with(executable_mock,
                                                   self.columns)

    @patch("pdm_utils.functions.querying.append_where_clauses")
    @patch("pdm_utils.functions.querying.select")
    @patch("pdm_utils.functions.querying.func.count")
    @patch("pdm_utils.functions.querying.build_fromclause")
    @patch("pdm_utils.functions.querying.extract_columns")
    def test_build_count_1(self, extract_columns_mock, 
                                 build_from_clause_mock, count_mock, select_mock,
                                 append_where_clauses_mock):
        """Verify function structure of build_count().
        """
        executable_mock = Mock()
        select_from_mock = Mock()
        type(executable_mock).select_from = select_from_mock
        select_from_mock.return_value = executable_mock

        extract_columns_mock.return_value = self.columns
        build_from_clause_mock.return_value = self.phage
        count_mock.return_value = self.count_column
        select_mock.return_value = executable_mock
        append_where_clauses_mock.return_value = executable_mock

        querying.build_count(self.graph, self.columns, 
                                            where=self.whereclauses)

        extract_columns_mock.assert_any_call(self.whereclauses)
        extract_columns_mock.assert_any_call(None, check=Column)
        total_columns = self.columns + self.columns + self.columns
        build_from_clause_mock.assert_called_once_with(self.graph, total_columns)
        
        for index in range(len(self.columns)):
            with self.subTest(columns_list_index=index):
                count_mock.assert_any_call(self.columns[index])

        count_param = [self.count_column, self.count_column, self.count_column]
        select_mock.assert_called_once_with(count_param)
        select_from_mock.assert_called_once_with(self.phage)
        append_where_clauses_mock.assert_called_once_with(executable_mock, 
                                                   self.whereclauses)

    @patch("pdm_utils.functions.querying.build_select")
    def test_build_distinct_1(self, build_select_mock):
        """Verify function structure of build_distinct()
        """
        executable_mock = Mock()
        DistinctMock = Mock()
        type(executable_mock).distinct = DistinctMock

        build_select_mock.return_value = executable_mock
        
        querying.build_distinct(self.graph, self.columns, 
                                            where=self.whereclauses, 
                                            order_by=self.columns)

        build_select_mock.assert_called_with(self.graph, self.columns,
                                            where=self.whereclauses, 
                                            order_by=self.columns,
                                            add_in=None)

class TestExecute(unittest.TestCase):
    def setUp(self):
        self.mock_engine = Mock()
        self.mock_proxy = Mock()

        self.mock_row_proxy = MagicMock()
        self.mock_results = [self.mock_row_proxy]

        self.data_tuple = ("Trixie", "Mycobacterium", "A")
        self.data_dict = {"PhageID"   : "Trixie",
                          "HostGenus" : "Mycobacterium",
                          "Cluster"   : "A"}
        self.mock_row_proxy.keys.return_value.__iter__.return_value = \
                                                    self.data_dict.keys()
        self.mock_row_proxy.__getitem__.side_effect = \
                                                    self.data_tuple

        self.mock_engine.execute.return_value = self.mock_proxy
        self.mock_proxy.fetchall.return_value = self.mock_results 

        self.mock_executable = Mock()
        self.mock_in_column = Mock(spec=Column)
        self.mock_executable.get_children.return_value = [self.mock_in_column]
        self.values = ["Trixie", "D29", "Myrna"]

    def test_execute_1(self):
        """Verify function structure of execute().
        """
        querying.execute(self.mock_engine, self.mock_executable)
        
        self.mock_engine.execute.assert_called()
        self.mock_proxy.fetchall.assert_called()

    def test_execute_2(self):
        """Verify execute() converts results to data dictionaries.
        """
        results = querying.execute(self.mock_engine, self.mock_executable)

        self.assertEqual(results, [self.data_dict])

    @patch("pdm_utils.functions.querying.dict")
    def test_execute_3(self, dict_mock):
        """Verify that execute() calls built-in function dict().
        """
        dict_mock.return_value = "dict_return_value"

        results = querying.execute(self.mock_engine, self.mock_executable)

        self.assertEqual(results, ["dict_return_value"])

    @patch("pdm_utils.functions.querying.dict")
    def test_execute_4(self, dict_mock):
        """Verify that parameter return_dict controls conversion with dict().
        """
        dict_mock.return_value = "dict_return_value"

        results = querying.execute(self.mock_engine, self.mock_executable,
                                                     return_dict=False)

        self.assertNotEqual(results, ["dict_return_value"])
        self.assertEqual(results, self.mock_results)

    @patch("pdm_utils.functions.querying.execute_value_subqueries")
    def test_execute_5(self, subqueries_mock):
        """Verify that execute() calls execute_value_subqueries().
        """
        querying.execute(self.mock_engine, self.mock_executable,
                         values=self.values, in_column=self.mock_in_column,
                         limit=8001, return_dict=False)

        subqueries_mock.assert_called_with(self.mock_engine,
                                           self.mock_executable,
                                           self.mock_in_column,
                                           self.values,
                                           limit=8001, 
                                           return_dict=False)

    def test_execute_6(self):
        """Verify that execute() raises ValueError with lacking instruction.
        """
        with self.assertRaises(ValueError):
            querying.execute(self.mock_engine, self.mock_executable,
                             values=self.values)

    def test_first_column_1(self):
        """Verify function structure of first_column().
        """
        querying.first_column(self.mock_engine, self.mock_executable)

        self.mock_engine.execute.assert_called()
        self.mock_proxy.fetchall.assert_called()

    def test_first_column_2(self):
        """Verify that first_column() retreives the first value from each row.
        """
        results = querying.first_column(self.mock_engine, self.mock_executable)

        self.assertEqual(results, [self.data_tuple[0]])

    @patch("pdm_utils.functions.querying.first_column_value_subqueries")
    def test_first_column_3(self, subqueries_mock):
        """Verify that first_column() calls first_column_value_subqueries().
        """
        querying.first_column(self.mock_engine, self.mock_executable,
                              in_column=self.mock_in_column, 
                              values=self.values, limit=8001)

        subqueries_mock.assert_called_with(self.mock_engine,
                                           self.mock_executable,
                                           self.mock_in_column,
                                           self.values,
                                           limit=8001)

    def test_first_column_4(self):
        """Verify first_column() raises ValueError with lacking instructions.
        """
        with self.assertRaises(ValueError):
            querying.execute(self.mock_engine, self.mock_executable,
                             values=self.values)

    def test_execute_value_subqueries_1(self):
        """Verify execute_value_subqueries() raises ValueError from bad column.
        """
        with self.assertRaises(ValueError):
            querying.execute_value_subqueries(
                                          self.mock_engine,
                                          self.mock_executable,
                                          None,
                                          self.values)

    def test_execute_value_subqueries_2(self):
        """Verify execute_value_subqueries() raises ValueError from bad column.
        """
        self.mock_executable.is_derived_from.return_value = False

        with self.assertRaises(ValueError):
            querying.execute_value_subqueries(
                                          self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

    @patch("pdm_utils.functions.querying.parsing.convert_to_encoded")
    def test_execute_value_subqueries_3(self, convert_to_encoded_mock):
        """Verify execute_value_subqueries() calls convert_to_encoded().
        """
        type_mock = Mock()
        type(type_mock).python_type = PropertyMock(return_value=bytes)
        type(self.mock_in_column).type = type_mock

        convert_to_encoded_mock.return_value = self.values

        querying.execute_value_subqueries(self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

        convert_to_encoded_mock.assert_called_with(self.values)

    def test_execute_value_subqueries_4(self):
        """Verify that execute_value_subqueries() calls Column.in_().
        """
        querying.execute_value_subqueries(self.mock_engine, 
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

        self.mock_in_column.in_.assert_called_with(self.values)

    def test_execute_value_subqueries_5(self):
        """Verify that execute_value_subqueries chunks values correctly.
        """
        querying.execute_value_subqueries(self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values,
                                          return_dict=False,
                                          limit=2)
        
        self.mock_in_column.in_.assert_any_call(self.values[:2])
        self.mock_in_column.in_.assert_any_call([self.values[2]])

    def test_first_column_value_subqueries_1(self):
        """Verify that first_column_value_subqueries() raises ValueError.
        First_column_value_subqueries should raise when inputted column
        is not a Column object.
        """
        with self.assertRaises(ValueError):
            querying.first_column_value_subqueries(
                                          self.mock_engine,
                                          self.mock_executable,
                                          None,
                                          self.values)

    def test_first_column_value_subqueries_2(self):
        """Verify that first_column_value_subqueries() raises ValueError.
        first_column_value_subqueries() should raise when column does not
        belong to any of the tables joined for the query.
        """
        self.mock_executable.is_derived_from.return_value = False
        with self.assertRaises(ValueError):
            querying.first_column_value_subqueries(
                                          self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

    @patch("pdm_utils.functions.querying.parsing.convert_to_encoded")
    def test_first_column_value_subqueries_3(self, convert_to_encoded_mock):
        """Verify first_column_value_subqueries() calls convert_to_encoded().
        """
        type_mock = Mock()
        type(type_mock).python_type = PropertyMock(return_value=bytes)
        type(self.mock_in_column).type = type_mock

        convert_to_encoded_mock.return_value = self.values

        querying.first_column_value_subqueries(self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

        convert_to_encoded_mock.assert_called_with(self.values)

    def test_first_column_value_subqueries_4(self):
        """Verify first_column_value_subqueries() calls Column.in_().
        """
        querying.first_column_value_subqueries(self.mock_engine, 
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values)

        self.mock_in_column.in_.assert_called_with(self.values)

    def test_first_column_value_subqueries_5(self):
        """Verify that first_column_value_subqueries() chunks values correctly.
        """
        querying.first_column_value_subqueries(self.mock_engine,
                                          self.mock_executable,
                                          self.mock_in_column,
                                          self.values,
                                          limit=2)
        
        self.mock_in_column.in_.assert_any_call(self.values[:2])
        self.mock_in_column.in_.assert_any_call([self.values[2]])

if __name__ == "__main__":
    unittest.main()
