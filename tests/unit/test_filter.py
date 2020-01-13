"""Unit unittests for the filter module."""

from pdm_utils.classes.filter import HistoryNode, Filter, CmdFilter
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from unittest.mock import patch, Mock, PropertyMock
import unittest

class TestFilter(unittest.TestCase):
    @patch("pdm_utils.classes.databasetree.DatabaseTree")
    def setUp(self, DatabaseTreeMock):
        DatabaseTreeMock.return_value.\
                show_tables.return_value = ["phage"]

        self.db_tree_mock = DatabaseTreeMock
        self.sql_handle = MySQLConnectionHandler()
        self.db_filter = Filter(self.sql_handle)
   
    def test_translate_table_1(self):
        "Verify translate_table() returns id match as expected."
        print(self.db_filter.db_tree.show_tables())
        table = self.db_filter.translate_table("phage")
        self.assertEqual(table, "phage")
        table = self.db_filter.translate_table("pHaGe")
        self.assertEqual(table, "phage")

    def test_translate_table_2(self):
        "Verify translate_table() raises ValueError as expected."
        with self.assertRaises(ValueError):
            self.db_filter.translate_table("test")

    def test_translate_field_1(self):
        "Verify translate_field() returns id match as expected."
        self.db_tree_mock.return_value.\
                get_table.return_value.\
                show_columns.return_value = ["phageID"]

        field = self.db_filter.translate_field("phageID", "phage")
        self.assertEqual(field, "phageID")
        field = self.db_filter.translate_field("pHaGeiD", "phage")
        self.assertEqual(field, "phageID")

    def test_translate_field_2(self):
        "Verify translate_field() raises ValueError as expected."
        self.db_tree_mock.return_value.\
                get_table.return_value.\
                show_columns.return_value = ["phageID"]

        with self.assertRaises(ValueError):
            self.db_filter.translate_field("gene", "phage")  

    def test_check_operator_1(self):
        "Verify as expected."
        self.db_tree_mock.return_value.\
                get_table.return_value.\
                show_columns.return_value = [""]

        self.db_tree_mock.return_value.\
                get_table.return_value.\
                show_columns.return_value = ["phageID"]  
    
    def test_build_queries_1(self):
        "Verify as expected."
        pass

    def test_transpose_1(self):
        "Verify as expected."
        pass

    def test_add_history_1(self):
        "Verify as expected."
        pass

    def test_undo_1(self):
        "Verify as expected."
        pass

    def test_refresh_1(self):
        "Verify as expected."
        pass

    def test_switch_table_1(self):
        "Verify as expected."
        pass
    
    def test_add_filter_1(self):
        "Verify as expected."
        pass

    def test_set_values_1(self):
        "Verify as expected."
        pass

    def test_update_1(self):
        "Verify as expected."
        pass

    def test_sort(self):
        "Verify as expected."
        pass

    def test_reset_1(self):
        "Verify as expected."
        pass

    def test_results_1(self):
        "Verify as expected."
        pass

    def test_hits_1(self):
        "Verify as expected."
        pass

    def test_group_1(self):
        "Verify as expected."
        pass

    def test_build_groups_1(self):
        "Verify as expected."
        pass

    def test_group_transpose_1(self):
        "Verify as expected."
        pass

    def test_large_num_set_distinct_query_1(self):
        "Verify as expected."
        pass

    def test_copy_1(self):
        "Verify as expected."
        pass

    def test_copy_values_1(self):
        "Verify as expected."
        pass

    def test_write_csv(self):
        "Verify as expected."
        pass

class TestHistoryNode(unittest.TestCase):
    def setUp(self):
        "Verify as expected."
        pass

    def test_init_1(self):
        "Verify as expected."
        pass

    def test_get_id_1(self):
        "Verify as expected."
        pass

    def test_get_history_1(self):
        "Verify as expected."
        pass

    def test_has_next_1(self):
        "Verify as expected."
        pass

    def test_get_next_1(self):
        "Verify as expected."
        pass

    def test_add_next_1(self):
        "Verify as expected."
        pass

    def test_create_next_1(self):
        "Verify as expected."
        pass

    def test_remove_next_1(self):
        "Verify as expected."
        pass

class TestCmdFilter(unittest.TestCase):
    def setUp(self):
        "Verify as expected."
        pass

    def test_preloop_1(self):
        "Verify as expected."
        pass

    def test_do_add_1(self):
        "Verify as expected."
        pass

    def test_do_update_1(self):
        "Verify as expected."
        pass

    def test_do_group_1(self):
        "Verify as expected."
        pass

    def test_do_results_1(self):
        "Verify as expected."
        pass

    def test_do_hits_1(self):
        "Verify as expected."
        pass

    def test_do_switch_1(self):
        "Verify as expected."
        pass

    def test_do_reset_1(self):
        "Verify as expected."
        pass

    def test_do_dump_1(self):
        "Verify as expected."
        pass

    def test_do_set_1(self):
        "Verify as expected."
        pass

    def test_do_show_1(self):
        "Verify as expected."
        pass

    def test_do_undo_1(self):
        "Verify as expected."
        pass

    def test_do_clear_1(self):
        "Verify as expected."
        pass

    def test_do_exit_1(self):
        "Verify as expected."
        pass

if __name__ == "__main__":
    unittest.main()
        

