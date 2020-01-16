"""Unit unittests for the schemagraph module"""

from pdm_utils.classes.schemagraph import  Node, \
                                           DatabaseNode, TableNode, ColumnNode
from pdm_utils.classes import schemagraph, mysqlconnectionhandler
from unittest.mock import Mock, patch
import unittest

class TestNode(unittest.TestCase):
    "Unittests for the Node class"
    def setUp(self):
        self.node = Node("test")

    def test_init_1(self):
        "Verify Node attributes are created as expected."

        self.assertEqual(self.node.id, "test")
        self.assertEqual(self.node.parents, [])
        self.assertEqual(self.node.children, [])
        
    def test_init_2(self):
        "Verify attribute boolean expressions function as expected."

        node = Node("test")

        self.assertTrue(not node.parents)
        self.assertTrue(not node.children)
        self.assertFalse(node == self.node) 

    def test_init_3(self):
        "Verify Node id type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            node = Node([])
 
    def test_show_parents_1(self):
        "Verify show_parents() creates return list as expected."

        parent_node_1 = Node("1")
        parent_node_2 = Node("2")
        parent_node_3 = Node("3")

        self.node.parents.append(parent_node_1)
        self.node.parents.append(parent_node_2)
        self.node.parents.append(parent_node_3)

        parents = self.node.show_parents()
        self.assertTrue("1" in parents)
        self.assertTrue("2" in parents)
        self.assertTrue("3" in parents)
    
    def test_show_children_1(self):
        "Verify show_children() creates return list as expected."

        child_node_1 = Node("1")
        child_node_2 = Node("2")
        child_node_3 = Node("3")

        self.node.children.append(child_node_1)
        self.node.children.append(child_node_2)
        self.node.children.append(child_node_3)

        children = self.node.show_children()
        self.assertTrue("1" in children)
        self.assertTrue("2" in children)
        self.assertTrue("3" in children)

    def test_has_parent_1(self):
        "Verify has_parent() identifies linked Nodes as expected."

        parent_node = Node("parent")
        random_node = Node("random")
        self.node.parents.append(parent_node) 

        self.assertTrue(self.node.has_parent("parent"))
        self.assertFalse(self.node.has_parent("random"))
   
    def test_has_parent_2(self):
        "Verify has_parent() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.has_parent([])

    def test_has_child_1(self):
        "Verify has_child() identifies linked Nodes as expected."

        child_node = Node("child")
        random_node = Node("random")
        self.node.children.append(child_node)

        self.assertTrue(self.node.has_child("child"))
        self.assertFalse(self.node.has_child("random"))

    def test_has_child_2(self):
        "Verify has_child() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.has_child([])

    def test_get_parent_1(self):
        "Verify get_parent() returns linked Nodes as expected."

        parent_node = Node("parent")
        random_node = Node("random")
        self.node.parents.append(parent_node)

        self.assertEqual(parent_node, 
                         self.node.get_parent("parent"))
        self.assertEqual(None,
                         self.node.get_parent("random"))

    def test_get_parent_2(self):
        "Verify get_parent() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            parent_node = Node([])

    def test_get_child_1(self):
        "Verify get_child() returns linked Nodes as expected."

        child_node = Node("child")
        random_node = Node("random")
        self.node.children.append(child_node)

        self.assertEqual(child_node,
                         self.node.get_child("child"))
        self.assertEqual(None,
                         self.node.get_parent("random"))

    def test_get_child_2(self):
        "Verify get_child() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.get_child([])

    def test_add_parent_1(self):
        "Verify add_parent() appends to node attributes as expected."

        parent_node = Node("parent")
        self.node.add_parent(parent_node)
        
        self.assertTrue(parent_node in self.node.parents)
        self.assertTrue(self.node in parent_node.children)

    def test_add_parent_2(self):
        "Verify add_parent() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.add_parent([])

    def test_add_parent_3(self):
        "Verify add_parent() type check raises ValueError as expected."
        parent_node = Node("parent")
        duplicate_parent = Node("parent")

        with self.assertRaises(ValueError):
            self.node.add_parent(parent_node)
            self.node.add_parent(parent_node)

        with self.assertRaises(ValueError):
            self.node.add_parent(parent_node)
            self.node.add_parent(duplicate_parent)

    @patch("pdm_utils.classes.schemagraph.Node.has_parent")
    def test_add_parent_4(self, HasParent):
        "Verify add_parent() calls has_parent() as expected."
        self.node.has_parent("parent")

        HasParent.assert_called_with("parent")

    def test_add_child_1(self):
        "Verify add_child() appends to node attributes as expected."

        child_node = Node("child")
        self.node.add_child(child_node)
        
        self.assertTrue(child_node in self.node.children)
        self.assertTrue(self.node in child_node.parents)
    
    def test_add_child_2(self):
        "Verify add_child() type check raises Type Error as expected."
        
        with self.assertRaises(TypeError):
            self.node.add_child([])

    def test_add_child_3(self):
        "Verify add_child() type check raises ValueError as expected."
        child_node = Node("child")
        duplicate_child = Node("child")

        with self.assertRaises(ValueError):
            self.node.add_child(child_node)
            self.node.add_child(child_node)

        with self.assertRaises(ValueError):
            self.node.add_child(child_node)
            self.node.add_child(duplicate_child)
    
    @patch("pdm_utils.classes.schemagraph.Node.has_child")
    def test_add_child_4(self, HasChild):
        "Verify add_child() calls has_child() as expected."
        self.node.has_child("child")

        HasChild.assert_called_with("child")

    def test_create_parent_1(self):
        "Verify create_parent() creates a Node object as expected."
        parent_node = self.node.create_parent("parent")

        self.assertTrue(isinstance(parent_node, Node))

    def test_create_parent_2(self):
        "Verify create_parent() appends to Node attributes as expected."

        parent_node = self.node.create_parent("parent")

        self.assertTrue(parent_node in self.node.parents)
        self.assertTrue(self.node in parent_node.children)

    def test_create_parent_3(self):
        "Verify create_parent() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.create_parent([])

    @patch("pdm_utils.classes.schemagraph.Node.add_parent")
    def test_create_parent_4(self, AddParent):
        "Verify create_parent() calls add_parent() as expected."

        parent_node = self.node.create_parent("parent")

        AddParent.assert_called_with(parent_node)

    def test_create_child_1(self):
        "Verify create_child() creates a Node object as expected."
        child_node = self.node.create_child("child")

        self.assertTrue(isinstance(child_node, Node))

    def test_create_child_2(self):
        "Verify create_child() appends to Node attributes as expected."

        child_node = self.node.create_child("child")

        self.assertTrue(child_node in self.node.children)
        self.assertTrue(self.node in child_node.parents)

    def test_create_child_3(self):
        "Verify create_child() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.create_child([])
    
    @patch("pdm_utils.classes.schemagraph.Node.add_child")
    def test_create_child_4(self, AddChild):
        "Verify create_child() calls add_child() as expected."

        child_node = self.node.create_child("child")

        AddChild.assert_called_with(child_node)

    def test_remove_parent_1(self):
        "Verify remove_parent() removes from Node attributes as expected."

        parent_node = Node("parent") 
        self.node.parents.append(parent_node)
        parent_node.children.append(self.node)
        
        self.assertEqual(parent_node, self.node.remove_parent(parent_node))

        self.assertFalse(parent_node in self.node.parents)
        self.assertFalse(self.node in parent_node.children)

    def test_remove_parent_2(self):
        "Verify remove_parent() type check raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.remove_parent([])

    def test_remove_child_1(self):
        "Verify remove_child() removes from Node attributes as expected."

        child_node = Node("child")
        self.node.children.append(child_node)
        child_node.parents.append(self.node)

        self.assertEqual(child_node, self.node.remove_child(child_node))

        self.assertFalse(child_node in self.node.children)
        self.assertFalse(self.node in child_node.parents)

    def test_remove_child_2(self):
        "Verify remove_child() type ccheck raises TypeError as expected."

        with self.assertRaises(TypeError):
            self.node.remove_parent([])

class TestDatabaseNode(unittest.TestCase):
    "Unittests for the DatabaseNode class."
    def setUp(self):
        self.db_node = DatabaseNode("test")

    @patch("pdm_utils.classes.schemagraph.DatabaseNode.add_child")
    def add_table_1(self, AddChild):
        "Verify add_table() calls add_child() as expected."
        table_node = TableNode("table")
        self.db_node.add_table(table_node)

        AddChild.assert_called_with(table_node)

    def add_table_3(self):
        "Verify add_table() type check raises TypeError as expected."
        with self.assertRaises(TypeError):
            self.db_node.add_table([])

    def create_table_1(self):
        "Verify create_table() creates a TableNode object as expected."
        table_node = self.db_node.create_table("table")
        self.assertTrue(isinstance(table_node, TableNode))

    def create_table_2(self):
        "Verify create_table() type check raises a TypeError as expected."    
        with self.assertRaises(TypeError):
            self.db_node.create_table([])

    @patch("pdm_utils.classes.schemagraph.TableNode.add_table")
    def create_table_3(self, AddTable):
        "Verify create_table() calls add_table() as expected."
        table_node = self.table_node.create_table("table")

        AddTable.assert_called_with(table_node)

    @patch("pdm_utils.classes.schemagraph.DatabaseNode.show_children")
    def show_tables_1(self, ShowChildren):
        "Verify show_tables() calls show_children() as expected."
        self.db_node.show_tables()

        ShowChildren.assert_called()

    @patch("pdm_utils.classes.schemagraph.DatabaseNode.has_child")       
    def has_table_1(self, HasChild):
        "Verify has_table() calls has_child() as expected."
        self.db_node.has_table("table")

        HasChild.assert_called_with("table")

    @patch("pdm_utils.classes.schemagraph.DatabaseNode.get_child")
    def get_table_1(self, GetChild):
        "Verify get_table() calls get_child() as expected."
        self.db_node.get_table("table")

        GetChild.assert_called_with("table")

class TestTableNode(unittest.TestCase):
    "Unittests for the TableNode class." 
    def setUp(self):
        self.table_node = TableNode("test")

    def test_init_1(self):
        "Verify TableNode attributes are created as expected."
        self.assertEqual(None, self.table_node.primary_key)

    @patch("pdm_utils.classes.schemagraph.TableNode.add_child")
    def test_add_column_1(self, AddChild):
        "Verify add_column() calls add_child() as expected."
        column_node = ColumnNode("column")
        self.table_node.add_column(column_node) 

        AddChild.assert_called_with(column_node)

    def test_add_column_2(self):
        "Verify add_column() type check raises TypeError as expected."
        with self.assertRaises(TypeError):
            self.table_node.add_column([])

    def test_create_column_1(self):
        "Verify create_column() creates a ColumnNode as expected."
        column_node = self.table_node.create_column("column")
        
        self.assertTrue(isinstance(column_node, ColumnNode))

    def test_create_column_2(self):
        "Verify create_column() type check raises TypeError as expected."
        with self.assertRaises(TypeError):
            self.table_node.create_column([])

    @patch("pdm_utils.classes.schemagraph.TableNode.add_column")
    def test_create_column_3(self, AddColumn):
        "Verify create_column() calls add_column() as expected."
        column_node = self.table_node.create_column("column") 

        AddColumn.assert_called_with(column_node)

    @patch("pdm_utils.classes.schemagraph.TableNode.show_children")     
    def test_show_columns_1(self, AddChildren):
        "Verify show_columns() calls show_children() as expected."
        self.table_node.show_columns()

        AddChildren.assert_called()

    @patch("pdm_utils.classes.schemagraph.ColumnNode.show_info")
    def test_show_columns_info_1(self, ShowColumnsInfo):
        "Verify show_columns_info() calls ColumnNode.show_info() as expected."
        column_node = ColumnNode("column")
        self.table_node.children.append(column_node)
        
        self.table_node.show_columns_info()

        ShowColumnsInfo.assert_called_once()

    def test_show_columns_info_2(self):
        "Verify show_columns_info() returns a list of lists as expected."
        column_node = ColumnNode("column")
        self.table_node.children.append(column_node)

        info = self.table_node.show_columns_info()
        self.assertTrue(isinstance(info, list))
        for info_list in info:
            self.assertTrue(isinstance(info_list, list))

    @patch("pdm_utils.classes.schemagraph.TableNode.has_child")
    def test_has_column_1(self, HasChild):
        "Verify has_column() calls has_child() as expected."
        self.table_node.has_column("column")

        HasChild.assert_called_with("column")

    @patch("pdm_utils.classes.schemagraph.TableNode.get_child")
    def test_get_column_1(self, GetChild):
        "Verify get_column() calls get_child() as expected."
        self.table_node.get_child("column")

        GetChild.assert_called_with("column")

    def test_show_foreign_keys_1(self):
        "Verify show_foreign_keys() returns a list as expected."
        other_table = TableNode("other")

        column_node_1 = ColumnNode("1", parents=[other_table, self.table_node]) 
        column_node_2 = ColumnNode("2", parents=[self.table_node])
        column_node_3 = ColumnNode("3", parents=[other_table, self.table_node])

        self.table_node.children.append(column_node_1)
        self.table_node.children.append(column_node_2)
        self.table_node.children.append(column_node_3)
        
        foreign_keys = self.table_node.show_foreign_keys()
        self.assertTrue(isinstance(foreign_keys, list))
        self.assertTrue("1" in foreign_keys)
        self.assertFalse("2" in foreign_keys)
        self.assertTrue("3" in foreign_keys)
    
    def test_show_foreign_keys_2(self):
        "Verify show_foreign_keys() type check raises TypeError as expected."
        self.table_node.children.append([])
        
        with self.assertRaises(TypeError):
            self.table_node.show_foreign_keys()

    def test_get_foreign_keys_1(self):
        "Verify get_foreign_keys() returns a list as expected."
        other_table = TableNode("other")

        column_node_1 = ColumnNode("1", parents=[other_table, self.table_node]) 
        column_node_2 = ColumnNode("2", parents=[self.table_node])
        column_node_3 = ColumnNode("3", parents=[other_table, self.table_node])

        self.table_node.children.append(column_node_1)
        self.table_node.children.append(column_node_2)
        self.table_node.children.append(column_node_3)
        
        foreign_keys = self.table_node.get_foreign_keys()
        self.assertTrue(isinstance(foreign_keys, list))
        self.assertTrue(column_node_1 in foreign_keys)
        self.assertFalse(column_node_2 in foreign_keys)
        self.assertTrue(column_node_3 in foreign_keys)

    def test_get_foreign_keys_2(self):
        "Verify get_foreign_keys() type check raises TypeError as expected."
        self.table_node.children.append([])
        
        with self.assertRaises(TypeError):
            self.table_node.get_foreign_keys()

    def test_show_primary_key_1(self):
        "Verify show_primary_key returns primary_key.id as expected."
        column_node = ColumnNode("column")

        self.table_node.primary_key = column_node

        self.assertEqual(self.table_node.show_primary_key(), "column")

    def test_show_primary_key_2(self):
        "Verify show_primary_key handles None primary_key as expected."
        self.assertEqual(self.table_node.show_primary_key(), "")
        
class TestColumnNode(unittest.TestCase):
    "Unittests for the ColumnNode class."
    def setUp(self):
        self.column_node = ColumnNode("test")

    def test_init_1(self):
        "Verify ColumnNode attribtues are created as expected."
        self.assertEqual(self.column_node.type, "")
        self.assertEqual(self.column_node.null, None)
        self.assertEqual(self.column_node.key, None)

    @patch("pdm_utils.classes.schemagraph.ColumnNode.add_parent")
    def test_add_table_1(self, AddParent):
        "Verify add_table() calls add_parent() as expected."
        table_node = TableNode("table")
        self.column_node.add_table(table_node)

        AddParent.assert_called_with(table_node)

    def test_add_table_2(self):
        "Verify add_table() type check raises TypeError as expected."
        with self.assertRaises(TypeError):
            self.column_node.add_table([])

    def test_create_table_1(self):
        "Verify create_table() returns TableNode object as expected."
        table_node = self.column_node.create_table("table")

        self.assertTrue(isinstance(table_node, TableNode))

    def test_create_table_2(self):
        "Verify create_table() type check raises TypeError as expected."
        with self.assertRaises(TypeError):
            self.column_node.create_table([])

    @patch("pdm_utils.classes.schemagraph.ColumnNode.add_table")
    def test_create_table_3(self, AddTable):
        "Verify create_table() calls add_table() as expected."
        table_node = self.column_node.create_table("table")

        AddTable.assert_called_with(table_node)

    def test_get_type_1(self):
        "Verify get_type() returns a string as expected."
        type = self.column_node.get_type()
        self.assertEqual(type, "")

        self.column_node.type = "varchar(5)"
        type = self.column_node.get_type()
        self.assertEqual(type, "varchar(5)")

    def test_parse_type_1(self):
        "Verify parse_type() returns a string as expected."
        type = self.column_node.parse_type() 
        self.assertEqual(type, "")

        self.column_node.type = "varchar (5)"
        type = self.column_node.parse_type()
        self.assertEqual(type, "varchar")

    def test_show_info_1(self):
        "Verify show_info() returns a list as expected."
        info = self.column_node.show_info()
        self.assertEqual(info, ["test", "", "Undefined", None, None])

        self.column_node.type = "varchar(5)"
        self.column_node.group = "limited_set"
        self.column_node.null = "True"
        self.column_node.key = "PRI"

        info = self.column_node.show_info()
        self.assertEqual(info, ["test", "varchar", "limited_set", "True", "PRI"])

if __name__ == "__main__":
    unittest.main()
