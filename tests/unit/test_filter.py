from sqlalchemy.engine.base import Engine
from sqlalchemy.sql.elements import BinaryExpression
from pdm_utils.classes.filter import Filter
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from sqlalchemy import Column
from unittest.mock import patch, Mock, PropertyMock, call
import unittest

class TestFilter(unittest.TestCase):
    @patch("pdm_utils.classes.filter.isinstance")
    def setUp(self, IsInstance):
        IsInstance.return_value = True

        alchemist = Mock()
        engine = Mock()
        graph = Mock()
        key = Mock()

        alchemist.engine.return_value =  engine
        alchemist.graph.return_value = graph
        alchemist.connected.return_value = True
      
        self.mock_isinstance = IsInstance
        self.mock_alchemist = alchemist
        self.mock_engine = engine
        self.mock_graph = graph
        self.mock_key = key

        self.db_filter = Filter(alchemist=alchemist, key=key) 


    def test_constructor_1(self):
        db_filter = Filter()
       
    def test_constructor_2(self):
        self.mock_isinstance.assert_any_call(self.mock_alchemist, 
                                             AlchemyHandler)
        self.mock_isinstance.assert_any_call(self.mock_key, Column)


    def test_values_1(self):
        self.db_filter.values = ["Test1", "Test2"]
        self.assertFalse(self.db_filter.values_valid)

    def test_values_2(self):
        with self.assertRaises(TypeError):
            self.db_filter.values = "Hello"

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.q.build_distinct")
    def test_transpose_2(self, BuildDistinct, Check):
        empty_values = self.db_filter.transpose("gene.Notes")

        self.assertEqual(empty_values, [])

        BuildDistinct.assert_not_called()
        Check.assert_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.transpose")
    def test_retrieve_1(self, Transpose, Check): 
        empty_data = self.db_filter.retrieve("Error")
    
        self.assertEqual(empty_data, {})

        Transpose.assert_not_called()
        Check.assert_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_refresh_1(self, BuildValues, Check):
        self.db_filter.refresh() 

        Check.assert_called()
        BuildValues.assert_not_called()
        
    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_refresh_2(self, BuildValues, Check):
        BuildValues.return_value = ["Phage"]

        self.db_filter._values_valid = False
        self.db_filter.refresh()

        Check.assert_called()
        BuildValues.assert_called()

        self.assertTrue(self.db_filter.values_valid)
        self.assertEqual(self.db_filter.values, ["Phage"])

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_1(self, BuildValues, Refresh, Check):
        self.db_filter.update()

        Check.assert_called()
        Refresh.assert_not_called()
        BuildValues.assert_not_called()
        
    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_2(self, BuildValues, Refresh, Check):
        self.db_filter._values_valid = False

        self.db_filter.update()

        Check.assert_called()
        Refresh.assert_called()
        BuildValues.assert_not_called()

    @patch("pdm_utils.classes.filter.Filter.check")
    @patch("pdm_utils.classes.filter.Filter.refresh")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_update_3(self, BuildValues, Refresh, Check):
        self.db_filter._values_valid = False
        self.db_filter._updated = False

        self.db_filter.update()

        Check.assert_called()
        Refresh.assert_called()
        BuildValues.assert_called()

        self.assertTrue(self.db_filter._values_valid)
        self.assertTrue(self.db_filter._values_valid)
    
    @patch("pdm_utils.classes.filter.isinstance")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_sort_1(self, BuildValues, IsInstance):
        IsInstance.return_value = True
        BuildValues.return_value = ["Phage"]

        self.db_filter._values_valid = False

        self.db_filter.sort("column")

        IsInstance.assert_called()
        BuildValues.assert_called()

        self.assertTrue(self.db_filter._values_valid)
        self.assertEqual(self.db_filter.values, ["Phage"])

    @patch("pdm_utils.classes.filter.isinstance")
    @patch("pdm_utils.classes.filter.Filter.build_values")
    def test_sort_1(self, BuildValues, IsInstance):
        IsInstance.return_value = False

        with self.assertRaises(TypeError):
            self.db_filter.sort()

    def test_reset_1(self):
        self.db_filter._filters = {"Expression" : "Some Whereclause"}

        self.db_filter.reset()

        self.assertEqual(self.db_filter.filters, {})

    def test_reset_2(self):
        self.db_filter._values  = ["Phage1", "Phage2", "Phage3"]

        self.db_filter.reset()

        self.assertEqual(self.db_filter.values, [])

    def test_reset_3(self):
        self.db_filter._values_valid = False

        self.db_filter.reset()

        self.assertTrue(self.db_filter.values_valid)

    def test_reset_4(self):
        self.db_filter._updated = False

        self.db_filter.reset()

        self.assertTrue(self.db_filter.updated)

    def test_hits_1(self):
        self.assertEqual(self.db_filter.hits(), 0)

    def test_hits_2(self):
        self.db_filter._values = ["Phage1", "Phage2", "Phage3"]

        self.assertEqual(self.db_filter.hits(), 3)

    def test_copy_1(self):
        self.db_filter._updated = False
 
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.updated, self.db_filter.updated) 

    def test_copy_2(self): 
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.values_valid, self.db_filter.values_valid)

    def test_copy_3(self):
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.filters, self.db_filter.filters)

    def test_copy_4(self):
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.engine, self.db_filter.engine)

    def test_copy_5(self):
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.graph, self.db_filter.graph)

    def test_copy_6(self): 
        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.key, self.db_filter.key)

    def test_copy_7(self):
        self.db_filter._values = ["Values"]

        copy_filter = self.db_filter.copy()

        self.assertEqual(copy_filter.values, self.db_filter.values)

    def test_copy_filters_1(self):
        self.db_filter._filters = {"Filter1" : ["Filter1"], "Filter2" : ["Filter2"], "Filter3" : ["Filter3"]}

        filters_copy = self.db_filter.copy_filters()


        self.assertEqual(filters_copy, self.db_filter._filters)

    def test_copy_filters_2(self):
        self.db_filter._filters = {"Filter1" : ["Filter1"], "Filter2" : ["Filter2"], "Filter3" : ["Filter3"]}

        filters_copy = self.db_filter.copy_filters()

        self.db_filter._filters.update({"Filter2" : []})

        self.assertNotEqual(filters_copy, self.db_filter._filters)

if __name__ == "__main__":
     unittest.main()
