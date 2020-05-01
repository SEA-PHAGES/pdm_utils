from networkx import Graph
from sqlalchemy import Column
from sqlalchemy.engine.base import Engine
from sqlalchemy.sql.elements import BinaryExpression
from pdm_utils.classes.filter import Filter
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from unittest.mock import call
from unittest.mock import patch
from unittest.mock import Mock 
from unittest.mock import PropertyMock
import unittest

class TestFilter(unittest.TestCase):
    @patch("pdm_utils.classes.filter.isinstance")
    def setUp(self, is_instance_mock):
        is_instance_mock.return_value = True

        alchemist = Mock(spec=AlchemyHandler)
        engine = Mock(spec=Engine)
        graph = Mock(spec=Graph)
        key = Mock(spec=Column)
        proxy = Mock()
        
        type(alchemist).engine = PropertyMock(return_value=engine)
        type(alchemist).graph  = PropertyMock(return_value=graph)
        type(alchemist).connected = PropertyMock(return_value=True)

        engine.execute.return_value = proxy
        proxy.fetchall.return_value = []
      
        self.mock_isinstance = is_instance_mock
        self.mock_alchemist = alchemist
        self.mock_engine = engine
        self.mock_graph = graph
        self.mock_key = key

        self.db_filter = Filter(alchemist=alchemist, key=key) 
        self.assertEqual(self.db_filter.engine, self.mock_engine)
        self.assertEqual(self.db_filter.graph, self.mock_graph)

    def test_updated_1(self):
        """Verify that the upload property portrays Filter._upload.
        """
        self.db_filter._updated = True

        self.assertTrue(self.db_filter.updated)

    def test_updated_2(self):
        """Verify that upload property is immutable.
        """
        with self.assertRaises(AttributeError):
            self.db_filter.updated = True

    def test_values_valid_1(self):
        """Verify that the values_valid property portrays Filter._values_valid.
        """
        self.db_filter._values_valid = True

        self.assertTrue(self.db_filter.values_valid)

    def test_values_valid_2(self):
        """Verify that the values_valid property is immutable.
        """
        with self.assertRaises(AttributeError):
            self.db_filter.values_valid = True

    def test_connected_1(self):
        """Verify that the connected property portrays Filter._connected.
        """
        self.db_filter._connected = True

        self.assertTrue(self.db_filter.connected)

    def test_connected_2(self):
        """Verify that the connected property is immutable.
        """
        with self.assertRaises(AttributeError):
            self.db_filter.connected = True

    def test_engine_1(self):
        """Verify that the engine property portrays Filter._engine.
        """
        self.db_filter._engine = self.mock_engine

        self.assertEqual(self.db_filter.engine, self.mock_engine)

    def test_engine_2(self):
        """Verify that the engine property can set Filter._engine.
        """
        self.db_filter.engine = self.mock_engine

        self.assertEqual(self.db_filter._engine, self.mock_engine)

    def test_engine_3(self):
        """Verify that the engine property raises TypeError on invalid input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.engine = Mock()

    def test_graph_1(self):
        """Verify that the graph property portrays Filter._graph.
        """
        self.db_filter._graph = self.mock_graph

        self.assertEqual(self.db_filter.graph, self.mock_graph)

    def test_graph_2(self):
        """Verify that the graph property can set Filter._graph.
        """
        self.db_filter.graph = self.mock_graph

        self.assertEqual(self.db_filter._graph, self.mock_graph)

    def test_graph_3(self):
        """Verify that the graph property raises TypeError on invalid input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.graph = Mock()

    def test_values_1(self):
        """Verify that the values property portrays Filter._values.
        """
        self.db_filter._values = ["Test1", "Test2"]

        self.assertEqual(self.db_filter.values, ["Test1", "Test2"])
 
    def test_values_2(self):
        """Verify that the values property can set Filter._values.
        """
        self.db_filter.values = ["Test1", "Test2"]

        self.assertEqual(self.db_filter.values, ["Test1", "Test2"])

    def test_values_3(self):
        """Verify that the values property modifies the values_valid property.
        """
        self.db_filter.values = ["Test1", "Test2"]

        self.assertFalse(self.db_filter.values_valid)

    def test_values_2(self):
        """Verify that the values property raises TypeError on invalid input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.values = "Hello"

    def test_key_1(self):
        """Verify that the key property modifies the Filter._key
        """
        self.db_filter.key = self.mock_key  

        self.assertEqual(self.db_filter.key, self.mock_key)

    def test_key_2(self):
        """Verify that the key property raises TypeError on invalid input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.key = Mock()

    @patch("pdm_utils.classes.filter.AlchemyHandler")
    def test_connect_1(self, alchemyhandler_mock):
        """Verfiy that connect() returns when Filter is already connected.
        """
        alchemyhandler_mock.return_value = self.mock_alchemist

        self.db_filter._connected = True

        self.db_filter.connect()
        alchemyhandler_mock.assert_not_called()

    @patch("pdm_utils.classes.filter.AlchemyHandler")
    def test_connect_2(self, alchemyhandler_mock):
        """Verify that the Filter creates an AlchemyHandler.
        """
        alchemyhandler_mock.return_value = self.mock_alchemist
        
        self.db_filter._connected = False

        self.db_filter.connect()
        alchemyhandler_mock.assert_called()

    @patch("pdm_utils.classes.filter.AlchemyHandler")
    def test_connect_3(self, alchemyhandler_mock):
        """Verify that the Filter uses an AlchemyHandler to connect. 
        """
        alchemyhandler_mock.return_value = self.mock_alchemist

        self.db_filter._connected = False

        self.db_filter.connect()

        self.mock_alchemist.connect.assert_called()
        self.mock_alchemist.build_graph.assert_called()

    @patch("pdm_utils.classes.filter.isinstance")
    def test_link_1(self, isinstance_mock):
        """Verify link() calls isinstance().
        """
        self.db_filter.link(self.mock_alchemist)

        isinstance_mock.assert_called()

    def test_link_2(self):
        """Verify link() raises TypeError upon bad alchemist input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.link("Bad input")

    def test_link_3(self):
        """Verify function structure of link().
        """
        type(self.mock_alchemist).connected = PropertyMock(return_value=False)
        type(self.mock_alchemist).graph = PropertyMock(return_value=None)

        self.db_filter.link(self.mock_alchemist)


        self.mock_alchemist.connect.assert_called()
        self.mock_alchemist.build_graph.assert_called()

    @patch("pdm_utils.classes.filter.Filter.connect")
    def test_check_1(self, connect_mock):
        """Verify that the Filter will connect if not connected.
        """
        self.db_filter._connected = False
        
        self.db_filter.check()

        connect_mock.assert_called()

    @patch("pdm_utils.classes.filter.isinstance")
    def test_check_2(self, is_instance_mock):
        """Verify that the Filter calls isinstance() with correct paremeters.
        """
        is_instance_mock.return_value = True

        self.db_filter.check()

        is_instance_mock.assert_any_call(self.mock_engine, Engine)
        is_instance_mock.assert_any_call(self.mock_graph, Graph)
        is_instance_mock.assert_any_call(self.mock_key, Column)

    def test_check_3(self):
        """Verify that the Filter raises AttributeError with an invalid engine.
        """
        self.db_filter._engine = "Not a valid engine"

        with self.assertRaises(AttributeError):
            self.db_filter.check()

    def test_check_4(self):
        """Verify that the Filter raises AttributeError with an invalid graph.
        """
        self.db_filter._graph = "Not a valid graph"

        with self.assertRaises(AttributeError):
            self.db_filter.check()

    def test_check_5(self):
        """Verify that the Filter raises AttibuteError with an invalid key.
        """
        self.db_filter._key = "Not a valid key"

        with self.assertRaises(AttributeError):
            self.db_filter.check()

    @patch("pdm_utils.classes.filter.q.build_distinct")
    @patch("pdm_utils.classes.filter.isinstance")
    def test_build_values_1(self, is_instance_mock, build_distinct_mock):
        """Verify that build_distinct() is called with correct parameters.
        """
        is_instance_mock.return_value = True
        self.db_filter.build_values(where="Not a list")

        build_distinct_mock.assert_called_with(self.mock_graph, self.mock_key,
                                               where="Not a list", 
                                               add_in=self.mock_key)
 
    @patch("pdm_utils.classes.filter.Filter.check")
    def test_transpose_1(self, check_mock): 
        """Verify that transpose() returns without values.
        """
        self.db_filter.transpose("gene.Notes")

        check_mock.assert_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    def test_mass_transpose_1(self, check_mock): 
        """Verify that mass_tranpose() returns without values.
        """
        self.db_filter.mass_transpose("Column")
    
        check_mock.assert_called()
   
    @patch("pdm_utils.classes.filter.q.build_distinct")
    @patch("pdm_utils.classes.filter.Filter.check") 
    def test_retrieve_1(self, check_mock, build_distinct_mock):
        """Verify that retrieve() returns without values.
        """
        self.db_filter.retrieve("Column")
         
        check_mock.assert_called()
        build_distinct_mock.assert_not_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_refresh_1(self, build_values_mock, check_mock):
        """Verify that refresh() returns without values.
        """
        self.db_filter.refresh() 

        check_mock.assert_called()
        build_values_mock.assert_not_called()
        
    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_refresh_2(self, build_values_mock, check_mock):
        """Verify that refresh() calls build_values() and conserves values.
        """
        build_values_mock.return_value = ["Phage"]

        self.db_filter._values_valid = False
        self.db_filter.refresh()

        check_mock.assert_called()
        build_values_mock.assert_called()

        self.assertTrue(self.db_filter.values_valid)
        self.assertEqual(self.db_filter.values, ["Phage"])

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_1(self, build_values_mock, Refresh, check_mock):
        """Verify update() returns without values.
        """
        self.db_filter.update()

        check_mock.assert_called()
        Refresh.assert_not_called()
        build_values_mock.assert_not_called()
        
    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_2(self, build_values_mock, Refresh, check_mock):
        """Verify update() refreshes values before updating.
        """
        self.db_filter._values_valid = False

        self.db_filter.update()

        check_mock.assert_called()
        Refresh.assert_called()
        build_values_mock.assert_not_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_3(self, build_values_mock, Refresh, check_mock):
        """Verify function structure of update().
        """
        self.db_filter._values_valid = False
        self.db_filter._updated = False

        self.db_filter.update()

        check_mock.assert_called()
        Refresh.assert_called()
        build_values_mock.assert_called()

        self.assertTrue(self.db_filter._values_valid)
        self.assertTrue(self.db_filter._values_valid)
  
    @patch("pdm_utils.classes.filter.q.first_column")
    @patch("pdm_utils.classes.filter.q.build_select")
    @patch("pdm_utils.classes.filter.Filter.convert_column_input")
    def test_sort_1(self, convert_column_input_mock, build_select_mock,
                                                     first_column_mock):
        """Verify function structure of sort().
        """
        self.db_filter._values_valid = False

        first_column_mock.return_value = ["Phage"]

        self.db_filter.sort("column")

    
        convert_column_input_mock.assert_called()
        build_select_mock.assert_called()
        first_column_mock.assert_called()

        self.assertTrue(self.db_filter._values_valid)
        self.assertEqual(self.db_filter.values, ["Phage"])

    def test_sort_2(self):
        """Verify that sort() raises TypeError at bad ORDER BY input.
        """
        with self.assertRaises(TypeError):
            self.db_filter.sort(None)

    def test_reset_1(self):
        """Verify that reset() clears filters.
        """
        self.db_filter._filters = [{"Expression" : "Some Whereclause"}]

        self.db_filter.reset()

        self.assertEqual(self.db_filter.filters, [])

    def test_reset_2(self):
        """Verify that reset() clears values.
        """
        self.db_filter._values  = ["Phage1", "Phage2", "Phage3"]

        self.db_filter.reset()

        self.assertEqual(self.db_filter.values, [])

    def test_reset_3(self):
        """Verify that reset() sets values_valid Filter property.
        """
        self.db_filter._values_valid = False

        self.db_filter.reset()

        self.assertTrue(self.db_filter.values_valid)

    def test_reset_4(self):
        """Verify that reset() sets updated Filter property.
        """
        self.db_filter._updated = False

        self.db_filter.reset()

        self.assertTrue(self.db_filter.updated)

    def test_reset_5(self):
        """Verify that reset sets on_index property.
        """
        self.db_filter._on_index = 5

        self.db_filter.reset()

        self.assertEqual(self.db_filter.or_index, -1)

    def test_hits_1(self):
        """Verify that hits() accurately portrays no values.
        """
        self.assertEqual(self.db_filter.hits(), 0)

    def test_hits_2(self):
        """Verify that hits() accurately portrays a number of values.
        """
        self.db_filter._values = ["Phage1", "Phage2", "Phage3"]

        self.assertEqual(self.db_filter.hits(), 3)

    def test_copy_1(self):
        """Verify that copy() reflects a Filter's updated property.
        """
        self.db_filter._updated = False
 
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.updated, self.db_filter.updated) 

    def test_copy_2(self): 
        """Verify that copy() reflects a Filter's values_valid property.
        """
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.values_valid, self.db_filter.values_valid)

    def test_copy_3(self):
        """Verify that copy() reflects a Filter's filters property.
        """
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.filters, self.db_filter.filters)

    def test_copy_4(self):
        """Verify that copy() reflects a Filter's engine property.
        """
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.engine, self.db_filter.engine)

    def test_copy_5(self):
        """Verify that copy() reflects a Filter's graph property.
        """
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.graph, self.db_filter.graph)

    def test_copy_6(self): 
        """Verify that copy() reflects a Fitler's key property.
        """
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.key, self.db_filter.key)

    def test_copy_7(self):
        """Verify that copy() reflects a Filter's values property.
        """
        self.db_filter._values = ["Values"]

        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.values, self.db_filter.values)

    def test_copy_8(self):
        """Verify that copy() reflects a Filter's connected property.
        """
        self.db_filter._connected = True
        copy_filter = self.db_filter.copy()
        
        self.assertEqual(copy_filter._connected, self.db_filter._connected)

    def test_copy_filters_1(self):
        """Verify that copy_filters() replicates filters.
        """
        self.db_filter._filters = [{"Filter1" : "Filter1", "Filter2" : "Filter2", "Filter3" : "Filter3"}]

        filters_copy = self.db_filter.copy_filters()


        self.assertEqual(filters_copy, self.db_filter._filters)

    def test_copy_filters_2(self):
        """Verify that copy_filters() creates new address for copied filters.
        """
        self.db_filter._filters = [{"Filter1" : "Filter1", "Filter2" : "Filter2", "Filter3" : "Filter3"}]

        filters_copy = self.db_filter.copy_filters()

        self.db_filter._filters[0].update({"Filter2" : "filter2"})

        self.assertNotEqual(filters_copy, self.db_filter._filters)

if __name__ == "__main__":
     unittest.main()
