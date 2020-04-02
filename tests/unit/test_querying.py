from networkx import Graph
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import querying
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
from unittest.mock import Mock, patch, PropertyMock
import unittest

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
    def test_get_table_1(self, TranslateTable):
        TranslateTable.return_value = "phage"

        querying.get_table(self.metadata, "phage") 

        TranslateTable.assert_called_with(self.metadata, "phage")

    def test_get_table_2(self):
        table_obj = querying.get_table(self.metadata, "phage")

        self.assertEqual(table_obj, self.phage)
    
    @patch("pdm_utils.functions.querying.parsing.parse_column")
    def test_get_column_1(self, TranslateColumn):
        TranslateColumn.return_value = ["phage", "PhageID"]

        querying.get_column(self.metadata, "phage.PhageID")

        TranslateColumn.assert_called_with("phage.PhageID")

    @patch("pdm_utils.functions.querying.get_table")
    def test_get_column_2(self, GetTable):
        GetTable.return_value = self.phage

        querying.get_column(self.metadata, "phage.PhageID")

        GetTable.assert_called_with(self.metadata, "phage")

    @patch("pdm_utils.functions.querying.parsing.translate_column")
    def test_get_column_3(self, TranslateColumn):
        TranslateColumn.return_value = "PhageID"

        querying.get_column(self.metadata, "phage.PhageID")

        TranslateColumn.assert_called_with(self.metadata, "phage.PhageID")

    def test_get_column_4(self):
        column_obj = querying.get_column(self.metadata, "phage.PhageID")

        self.assertEqual(column_obj, self.PhageID)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_get_table_list_1(self, ExtractColumn):  
        ExtractColumn.return_value = self.PhageID

        querying.get_table_list([self.PhageID])

        ExtractColumn.assert_called_with(self.PhageID)
    
    @patch("pdm_utils.functions.querying.extract_column")
    def test_get_table_list_2(self, ExtractColumn):
        ExtractColumn.return_value = self.PhageID

        querying.get_table_list(self.PhageID)
    
        ExtractColumn.assert_called_with(self.PhageID)
    
    def test_get_table_list_3(self):
        table_list = querying.get_table_list(self.PhageID)

        self.assertTrue(self.phage in table_list)
        self.assertFalse(self.gene in table_list)

    def test_get_table_list_4(self):
        columns = [self.PhageID, self.GeneID]
        table_list = querying.get_table_list(columns)

        self.assertTrue(self.phage in table_list)
        self.assertTrue(self.gene in table_list)

    def test_build_graph_1(self):
        graph = querying.build_graph(self.metadata)

        self.assertTrue(isinstance(graph, Graph))
 
    def test_build_graph_2(self):
        graph = querying.build_graph(self.metadata)

        self.assertTrue("phage" in graph.nodes)
        self.assertTrue("gene"  in graph.nodes)
        self.assertTrue("trna"  in graph.nodes)

        self.assertEqual(graph.nodes["phage"]["table"], self.phage)
        self.assertEqual(graph.nodes["gene"]["table"],  self.gene)
        self.assertEqual(graph.nodes["trna"]["table"],  self.trna)

    def test_build_graph_3(self):
        graph = querying.build_graph(self.metadata)
        
        phage = graph.nodes["phage"]["table"]
        gene = graph.nodes["gene"]["table"]
        trna = graph.nodes["trna"]["table"]

        phage_gene_edge = graph["gene"]["phage"]
        phage_trna_edge = graph["trna"]["phage"]

        self.assertEqual(phage_gene_edge["key"], self.PhageGeneKey)
        self.assertEqual(phage_trna_edge["key"], self.PhageTrnaKey)
   
    @patch("pdm_utils.functions.querying.parsing.translate_table")
    def test_build_onclause_1(self, TranslateTable):
        TranslateTable.return_value = "phage"
        graph = querying.build_graph(self.metadata)
        
        table = querying.build_onclause(graph, "phage", "phage")

        TranslateTable.assert_called_with(self.metadata, "phage")
        self.assertEqual(table, self.phage)

    def test_build_onclause_2(self):
        graph = querying.build_graph(self.metadata) 
       
        querying.build_onclause(graph, "phage", "gene")

class TestExtract(unittest.TestCase):
    def setUp(self):
        self.column = Mock(spec=Column)
        self.columns = [self.column]
        self.check_type = Mock()
    
    @patch("pdm_utils.functions.querying.isinstance")
    def test_extract_column_1(self, IsInstance):
        with self.assertRaises(TypeError):
            IsInstance.return_value = False
            querying.extract_column(self.column)

        IsInstance.assert_any_call(self.column, Column)
        IsInstance.assert_any_call(self.column, functions.count)
        IsInstance.assert_any_call(self.column, UnaryExpression)
        IsInstance.assert_any_call(self.column, BinaryExpression)

    def test_extract_column_2(self):
        self.assertEqual(querying.extract_column(self.column), self.column)
        
    def test_extract_column_3(self):
        column_list = PropertyMock(return_value=self.columns)
        clause_list = Mock()
        count = Mock(spec=functions.count)
        
        type(count).clauses = clause_list
        type(clause_list).clauses = column_list

        self.assertEqual(querying.extract_column(count), self.column)

    def test_extract_column_4(self):
        unary_expression = Mock(spec=UnaryExpression)

        type(unary_expression).element = self.column

        self.assertEqual(querying.extract_column(unary_expression), 
                                                 self.column)

    def test_extract_column_5(self):
        binary_expression = Mock(spec=BinaryExpression)

        type(binary_expression).left = self.column

        self.assertEqual(querying.extract_column(binary_expression), 
                                                 self.column)

    def test_extract_column_6(self):
        unary_expression = Mock(spec=UnaryExpression)
        binary_expression = Mock(spec=BinaryExpression)

        type(unary_expression).element = self.column
        type(binary_expression).left = unary_expression

        self.assertEqual(querying.extract_column(binary_expression), 
                                                 self.column)
    
    @patch("pdm_utils.functions.querying.isinstance")
    def test_extract_column_7(self, IsInstance):
        IsInstance.return_value = True

        querying.extract_column(self.column, check=self.check_type)

        IsInstance.assert_any_call(self.column, self.check_type)

    def test_extract_column_8(self):
        with self.assertRaises(TypeError):
            querying.extract_column(self.column, check=str)

    @patch("pdm_utils.functions.querying.isinstance")
    def test_extract_columns_1(self, IsInstance):
        querying.extract_columns(self.columns)

        IsInstance.assert_any_call(self.columns, list)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_extract_columns_2(self, ExtractColumn):
        querying.extract_columns(self.columns)

        ExtractColumn.assert_called_with(self.column, check=None)

    @patch("pdm_utils.functions.querying.extract_column")
    def test_extract_columns_3(self, ExtractColumn):
        querying.extract_columns(self.column)

        ExtractColumn.assert_called_with(self.column, check=None)

    def test_extract_where_clauses_1(self):
        where_columns = querying.extract_where_clauses(None)
        
        self.assertEqual(where_columns, []) 

    @patch("pdm_utils.functions.querying.extract_columns")
    def test_extract_where_clauses_2(self, ExtractColumns):
        querying.extract_where_clauses(self.columns)

        ExtractColumns.assert_called_with(self.columns, check=BinaryExpression)
    
    def test_extract_order_by_clauses_1(self):
        order_by_columns = querying.extract_order_by_clauses(None)

        self.assertEqual(order_by_columns, [])

    @patch("pdm_utils.functions.querying.extract_columns")
    def test_extract_order_by_clauses_2(self, ExtractColumns):
        querying.extract_order_by_clauses(self.columns)

        ExtractColumns.assert_called_with(self.columns, check=Column)
        
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
    def test_get_table_pathing_1(self, ShortestPath):
        ShortestPath.return_value = self.table_path_1

        for index in range(len(self.table_names_1)):
            with self.subTest(table_index=index):
                querying.get_table_pathing(self.graph, self.table_list_1)

                ShortestPath.assert_any_call(self.graph, 
                                                "center",
                                                self.table_names_1[index])

                ShortestPath.clear()

    @patch("pdm_utils.functions.querying.shortest_path")
    def test_get_table_pathing_2(self, ShortestPath):
        ShortestPath.return_value = self.table_path_1

        for index in range(len(self.table_names_2)):
            with self.subTest(table_index=index):
                querying.get_table_pathing(self.graph, self.table_list_2,
                                                center_table=self.center)

                ShortestPath.assert_any_call(self.graph,
                                                "center",
                                                self.table_names_2[index])

                ShortestPath.clear()

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_1(self, Join, BuildOnClause):
        Join.return_value = self.center
        BuildOnClause.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        BuildOnClause.assert_any_call(self.graph, "center", "table_1")
        Join.assert_any_call(self.center, self.table_1,
                                          isouter=True,
                                          onclause="onclause")
    
    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_2(self, Join, BuildOnClause):         
        Join.return_value = self.center
        BuildOnClause.return_value = "onclause"
        
        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        BuildOnClause.assert_any_call(self.graph, "center", "table_2")
        Join.assert_any_call(self.center, self.table_2,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_3(self, Join, BuildOnClause):       
        Join.return_value = self.center
        BuildOnClause.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        BuildOnClause.assert_any_call(self.graph, "table_2", "table_3")
        Join.assert_any_call(self.center, self.table_3,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_4(self, Join, BuildOnClause):       
        Join.return_value = self.center
        BuildOnClause.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_1)  

        BuildOnClause.assert_any_call(self.graph, "table_1", "table_4")
        Join.assert_any_call(self.center, self.table_4,
                                          isouter=True,
                                          onclause="onclause")

    @patch("pdm_utils.functions.querying.build_onclause")
    @patch("pdm_utils.functions.querying.join")
    def test_join_pathed_tables_5(self, Join, BuildOnClause):
        Join.return_value = self.center
        BuildOnClause.return_value = "onclause"

        querying.join_pathed_tables(self.graph, self.table_pathing_2)

        BuildOnClause.assert_not_called()
        Join.assert_not_called()

class TestBuildClauses(unittest.TestCase):
    def setUp(self):
        self.graph = Mock()

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
    
        mock_nodes = PropertyMock(return_value=nodes_dict)
        type(self.graph).nodes = mock_nodes

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
    def test_build_where_clause_1(self, ParseFilter): 
        ParseFilter.return_value = ["phage", "Cluster", "!=", "A"]

        querying.build_where_clause(self.graph, "phage.Cluster != A")

        ParseFilter.assert_called_with("phage.Cluster != A")

    @patch("pdm_utils.functions.querying.parsing.check_operator")
    def test_build_where_clause_2(self, CheckOperator):
        querying.build_where_clause(self.graph, "gene.PhamID = 2")

        CheckOperator.assert_called_with("=", self.PhamID)
   
    @patch("pdm_utils.functions.querying.join_pathed_tables")
    @patch("pdm_utils.functions.querying.get_table_pathing")
    @patch("pdm_utils.functions.querying.get_table_list")
    def test_build_fromclause_1(self, GetTableList, GetTablePathing, 
                                                            JoinPathedTables):
        GetTableList.return_value = self.table_names
        GetTablePathing.return_value = self.pathing
        JoinPathedTables.return_value  = Mock()
        
        querying.build_fromclause(self.graph, self.tables)

        GetTableList.assert_called_with(self.tables)
        GetTablePathing.assert_called_with(self.graph, self.table_names)
        JoinPathedTables.assert_called_with(self.graph, self.pathing)

    @patch("pdm_utils.functions.querying.append_order_by_clauses")
    @patch("pdm_utils.functions.querying.append_where_clauses")
    @patch("pdm_utils.functions.querying.select")
    @patch("pdm_utils.functions.querying.build_fromclause")
    @patch("pdm_utils.functions.querying.extract_order_by_clauses")
    @patch("pdm_utils.functions.querying.extract_where_clauses")
    def test_build_select_1(self, ExtractWhereClauses, ExtractOrderByClauses,
                                  BuildFromClause, Select,
                                  AppendWhereClauses, AppendOrderByClauses):
        ExecutableMock = Mock()
        SelectFromMock = Mock()
        type(ExecutableMock).select_from = SelectFromMock
        SelectFromMock.return_value = ExecutableMock

        ExtractWhereClauses.return_value = self.columns
        ExtractOrderByClauses.return_value = self.columns
        BuildFromClause.return_value = self.phage
        Select.return_value = ExecutableMock
        AppendWhereClauses.return_value = ExecutableMock
        AppendOrderByClauses.return_value = ExecutableMock

        querying.build_select(self.graph, self.columns, 
                                            where=self.whereclauses,
                                            order_by=self.columns)

        ExtractWhereClauses.assert_called_once_with(self.whereclauses)
        ExtractOrderByClauses.assert_called_once_with(self.columns)
        total_columns = self.columns + self.columns + self.columns
        BuildFromClause.assert_called_once_with(self.graph, total_columns)
        
        Select.assert_called_once_with(self.columns)
        SelectFromMock.assert_called_once_with(self.phage)
        AppendWhereClauses.assert_called_once_with(ExecutableMock,
                                                   self.whereclauses)
        AppendOrderByClauses.assert_called_once_with(ExecutableMock,
                                                   self.columns)

    @patch("pdm_utils.functions.querying.append_where_clauses")
    @patch("pdm_utils.functions.querying.select")
    @patch("pdm_utils.functions.querying.func.count")
    @patch("pdm_utils.functions.querying.build_fromclause")
    @patch("pdm_utils.functions.querying.extract_where_clauses")
    def test_build_count_1(self, ExtractWhereClauses, 
                                 BuildFromClause, Count, Select,
                                 AppendWhereClauses):
        ExecutableMock = Mock()
        SelectFromMock = Mock()
        type(ExecutableMock).select_from = SelectFromMock
        SelectFromMock.return_value = ExecutableMock

        ExtractWhereClauses.return_value = self.columns
        BuildFromClause.return_value = self.phage
        Count.return_value = self.count_column
        Select.return_value = ExecutableMock
        AppendWhereClauses.return_value = ExecutableMock

        querying.build_count(self.graph, self.columns, 
                                            where=self.whereclauses)

        ExtractWhereClauses.assert_called_once_with(self.whereclauses)
        total_columns = self.columns + self.columns
        BuildFromClause.assert_called_once_with(self.graph, total_columns)
        
        for index in range(len(self.columns)):
            with self.subTest(columns_list_index=index):
                Count.assert_any_call(self.columns[index])

        count_param = [self.count_column, self.count_column, self.count_column]
        Select.assert_called_once_with(count_param)
        SelectFromMock.assert_called_once_with(self.phage)
        AppendWhereClauses.assert_called_once_with(ExecutableMock, 
                                                   self.whereclauses)

    @patch("pdm_utils.functions.querying.build_select")
    def test_build_distinct_1(self, BuildSelect):
        ExecutableMock = Mock()
        DistinctMock = Mock()
        type(ExecutableMock).distinct = DistinctMock

        BuildSelect.return_value = ExecutableMock
        
        querying.build_distinct(self.graph, self.columns, 
                                            where=self.whereclauses, 
                                            order_by=self.columns)

        BuildSelect.assert_called_with(self.graph, self.columns,
                                            where=self.whereclauses, 
                                            order_by=self.columns)

if __name__ == "__main__":
    unittest.main()
