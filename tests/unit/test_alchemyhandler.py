from sqlalchemy import MetaData
from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from unittest.mock import patch, Mock, PropertyMock
import unittest

class TestAlchemyHandler(unittest.TestCase):
    def setUp(self):
        self.alchemist = AlchemyHandler()

    def test_constructor_1(self):
        self.assertEqual(self.alchemist._database, None)
        self.assertEqual(self.alchemist._username, None)
        self.assertEqual(self.alchemist._password, None)

    def test_constructor_2(self):
        self.assertEqual(self.alchemist._engine, None)
        self.assertEqual(self.alchemist.metadata, None)
        self.assertEqual(self.alchemist.graph, None)
        self.assertEqual(self.alchemist.session, None)

    def test_constructor_3(self):
        self.assertFalse(self.alchemist.connected)
        self.assertFalse(self.alchemist.has_database)
        self.assertFalse(self.alchemist.has_credentials)

    def test_database_1(self):
        self.alchemist.database = "Test"
        self.assertTrue(self.alchemist.has_database)
        self.assertFalse(self.alchemist.connected)

    def test_username_1(self):
        self.alchemist.username = "Test"
        self.assertFalse(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_username_2(self):
        self.alchemist.password = "Test"
        self.alchemist.username = "Test"
        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_password_1(self):
        self.alchemist.password ="Test"
        self.assertFalse(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_password_2(self):
        self.alchemist.username = "Test"
        self.alchemist.password = "Test"
        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)
   
    def test_engine_1(self):
        self.alchemist.connected = True
        self.alchemist.engine = None

        self.assertFalse(self.alchemist.connected)
       
    def test_engine_2(self):
        with self.assertRaises(TypeError):
            self.alchemist.engine = "Test"

    @patch("pdm_utils.classes.alchemyhandler.input")
    def test_ask_database_1(self, Input):
        self.alchemist.ask_database()
        Input.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.input")
    def test_ask_database_2(self, Input):
        self.alchemist.has_database = False
        self.alchemist.connected = True

        self.alchemist.ask_database()
 
        self.assertTrue(self.alchemist.has_database)
        self.assertFalse(self.alchemist.connected)

    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_ask_credentials_1(self, GetPass):
        self.alchemist.ask_credentials()

        GetPass.assert_called()
 
    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_ask_credentials_2(self, GetPass):
        self.alchemist.has_credentials = False
        self.alchemist.connected = True

        self.alchemist.ask_credentials()

        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_1(self, CreateEngine, AskCredentials):
        self.alchemist.engine = None
        self.alchemist.connected = True
        self.alchemist.build_engine()

        CreateEngine.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_2(self, CreateEngine, AskCredentials):
        self.alchemist.username = "user"
        self.alchemist.password = "pass"
        self.alchemist.has_credentials = False

        self.alchemist.build_engine()

        AskCredentials.assert_called()
        login_string = "mysql+pymysql://user:pass@localhost/"
        CreateEngine.assert_called_with(login_string)

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_3(self, CreateEngine, AskCredentials):
        self.alchemist.has_credentials = True
        self.alchemist.connected = False
        self.alchemist.metadata = "Test"
        self.alchemist.graph = "Test"

        self.alchemist.build_engine()

        self.alchemist.connected = True
        self.assertEqual(self.alchemist.metadata, None)
        self.assertEqual(self.alchemist.graph, None)

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_1(self, BuildEngine, AskDatabase, AskCredentials):
        self.alchemist.has_credentials = True
        self.alchemist.connect()
        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_2(self, BuildEngine, AskDatabase, AskCredentials):
        self.alchemist.connect(ask_database=True)
        BuildEngine.assert_called()
        AskDatabase.assert_called()
        AskCredentials.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.exit")
    @patch("pdm_utils.classes.alchemyhandler.print")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_3(self, BuildEngine, AskDatabase, AskCredentials,
                                                           Print, Exit):
        self.alchemist.connected = False
        BuildEngine.side_effect = OperationalError("", "", "")

        self.alchemist.connect()
        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_called()
        Print.assert_called()
        Exit.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_execute_1(self, BuildEngine):
        MockEngine = Mock()
        MockProxy  = Mock()

        MockEngine.execute.return_value = MockProxy 
        MockProxy.fetchall.return_value = []

        self.alchemist._engine = MockEngine

        self.alchemist.execute("Executable")

        MockEngine.execute.assert_called_with("Executable")
        MockProxy.fetchall.assert_called()
        BuildEngine.assert_not_called() 
   
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_execute_2(self, BuildEngine):
        pass
     
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_scalar_1(self, BuildEngine):
        MockEngine = Mock()
        MockProxy  = Mock()
        
        MockEngine.execute.return_value = MockProxy
        MockProxy.scalar.return_value = "Scalar"

        self.alchemist._engine = MockEngine
       
        self.alchemist.scalar("Executable")

        MockEngine.execute.assert_called_with("Executable")
        MockProxy.scalar.assert_called()
        BuildEngine.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_scalar_2(self, BuildEngine):
        pass
 
    @patch("pdm_utils.classes.alchemyhandler.MetaData")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_metadata_1(self, AskDatabase, BuildEngine, MetaData):
        self.alchemist.has_database = False
        self.alchemist.connected = False

        self.alchemist.build_metadata()

        AskDatabase.assert_called()
        BuildEngine.assert_called()
        MetaData.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.MetaData")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_metadata_2(self, AskDatabase, BuildEngine, MetaData):
        self.alchemist.has_database = True
        self.alchemist.connected = True
        
        self.alchemist.build_metadata()

        AskDatabase.assert_not_called()
        BuildEngine.assert_not_called()
        MetaData.assert_called() 

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_table_1(self, BuildMetadata, TranslateTable):
        self.alchemist.metadata = "Metadata"

        self.alchemist.translate_table("Test")

        TranslateTable.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_table_2(self, BuildMetadata, TranslateTable):
        self.alchemist.metadata = None

        self.alchemist.translate_table("Test")

        TranslateTable.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_column_1(self, BuildMetadata, TranslateColumn):
        self.alchemist.metadata = "Metadata"

        self.alchemist.translate_column("Test")

        TranslateColumn.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_column_2(self, BuildMetadata, TranslateColumn):
        self.alchemist.metadata = None

        self.alchemist.translate_column("Test")

        TranslateColumn.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_table_1(self, BuildMetadata, GetTable):
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_table("Test")

        GetTable.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_table_2(self, BuildMetadata, GetTable):
        self.alchemist.metadata = None

        self.alchemist.get_table("Test")

        GetTable.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_column_1(self, BuildMetadata, GetColumn):
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_column("Test")

        GetColumn.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()
        
    @patch("pdm_utils.classes.alchemyhandler.querying.get_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_column_2(self, BuildMetadata, GetColumn):
        self.alchemist.metadata = None

        self.alchemist.get_column("Test")

        GetColumn.assert_called_with(None, "Test")
        BuildMetadata.assert_called()
 
    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_1(self, BuildMetadata, BuildGraph):
        BuildGraph.return_value = "Graph"

        self.alchemist.metadata = "Metadata"

        self.alchemist.build_graph()

        BuildMetadata.assert_not_called()
        BuildGraph.assert_called_with("Metadata")
        
        self.assertEqual(self.alchemist.graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_2(self, BuildMetadata, BuildGraph):
        BuildGraph.return_value = "Graph"

        self.alchemist.metadata = None

        self.alchemist.build_graph()

        BuildMetadata.assert_called()
        BuildGraph.assert_called_with(None)
        
        self.assertEqual(self.alchemist.graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.cartography.get_map")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_map_1(self, BuildMetadata, GetMap): 
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_map("Test")

        BuildMetadata.assert_not_called()
        GetMap.assert_called_with("Metadata", "Test")

    @patch("pdm_utils.classes.alchemyhandler.cartography.get_map")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_map_2(self, BuildMetadata, GetMap):
        self.alchemist.metadata = None

        self.alchemist.get_map("Test")

        BuildMetadata.assert_called()
        GetMap.assert_called_with(None, "Test") 

if __name__ == "__main__":
    unittest.main()
