import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from sqlalchemy import MetaData
from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError

from pdm_utils.classes.alchemyhandler import AlchemyHandler

class TestAlchemyHandler(unittest.TestCase):
    def setUp(self):
        self.alchemist = AlchemyHandler()

    def test_constructor_1(self):
        """Verify AlchemyHandler credentials are initialized as None.
        """
        self.assertEqual(self.alchemist._database, None)
        self.assertEqual(self.alchemist._username, None)
        self.assertEqual(self.alchemist._password, None)

    def test_constructor_2(self):
        """Verify AlchemyHandler data objects are initialized as None.
        """
        self.assertEqual(self.alchemist._engine, None)
        self.assertEqual(self.alchemist.metadata, None)
        self.assertEqual(self.alchemist.graph, None)
        self.assertEqual(self.alchemist.session, None)

    def test_constructor_3(self):
        """Verify AlchemyHandler data booleans are initialized as False.
        """
        self.assertFalse(self.alchemist.connected)
        self.assertFalse(self.alchemist.has_database)
        self.assertFalse(self.alchemist.has_credentials)

    def test_database_1(self):
        """Verify database property sets has_database.
        """
        self.alchemist.database = "Test"
        self.assertTrue(self.alchemist.has_database)
        self.assertFalse(self.alchemist.connected_database)

    def test_username_1(self):
        """Verify username property conserves has_credentials and connected.
        """
        self.alchemist.username = "Test"
        self.assertFalse(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_username_2(self):
        """Verify username property sets has_credentials with valid password.
        """
        self.alchemist.password = "Test"
        self.alchemist.username = "Test"
        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_password_1(self):
        """Verify password property conserves has_credentials and connected.
        """
        self.alchemist.password ="Test"
        self.assertFalse(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_password_2(self):
        """Verify password property sets has_credentials with valid username.
        """
        self.alchemist.username = "Test"
        self.alchemist.password = "Test"
        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)
   
    def test_engine_1(self):
        """Verify engine property sets connected.
        """
        self.alchemist.connected = True
        self.alchemist.engine = None

        self.assertFalse(self.alchemist.connected)
       
    def test_engine_2(self):
        """Verify engine property raises TypeError on bad engine input.
        """
        with self.assertRaises(TypeError):
            self.alchemist.engine = "Test"

    @patch("pdm_utils.classes.alchemyhandler.input")
    def test_ask_database_1(self, Input):
        """Verify ask_database() calls input().
        """
        self.alchemist.ask_database()
        Input.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.input") 
    def test_ask_database_2(self, Input):
        """Verify ask_database() sets has_database.
        """
        self.alchemist.has_database = False
        self.alchemist.connected = True

        self.alchemist.ask_database()
 
        self.assertTrue(self.alchemist.has_database)
        self.assertFalse(self.alchemist.connected)

    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_ask_credentials_1(self, GetPass):
        """Verify ask_credentials() calls getpass().
        """
        self.alchemist.ask_credentials()

        GetPass.assert_called()
 
    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_ask_credentials_2(self, GetPass):
        """Verify ask_credentials() sets has_credentials.
        """
        self.alchemist.has_credentials = False
        self.alchemist.connected = True

        self.alchemist.ask_credentials()

        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)
    
    def test_validate_database_1(self):
        """Verify function structure of validate_database().
        """
        MockEngine = Mock()
        MockProxy = Mock()

        MockEngine.execute.return_value = MockProxy 
        MockProxy.fetchall.return_value = [("pdm_test_db",), 
                                           ("Actinobacteriophage",)]

        self.alchemist.connected = True
        self.alchemist.database = "pdm_test_db"
        self.alchemist._engine = MockEngine

        self.alchemist.validate_database()

        MockEngine.execute.assert_called_with("SHOW DATABASES")
        MockProxy.fetchall.assert_called()

    def test_validate_database_2(self):
        """Verify validate_database() raises IndexError without database.
        """
        self.alchemist.connected = True

        with self.assertRaises(AttributeError):
            self.alchemist.validate_database()

    def test_validate_database_3(self):
        """Verify validate_database() raises ValueError from bad database input.
        """
        MockEngine = Mock()
        MockProxy = Mock()

        MockEngine.execute.return_value = MockProxy
        MockProxy.fetchall.return_value = []

        self.alchemist.connected = True
        self.alchemist.database = "test db"
        self.alchemist._engine = MockEngine

        with self.assertRaises(ValueError):
            self.alchemist.validate_database()

        MockEngine.execute.assert_called_with("SHOW DATABASES")
        MockProxy.fetchall.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_1(self, create_engine_mock, ask_credentials_mock):
        """Verify build_engine() returns if connected already.
        """
        self.alchemist.engine = None
        self.alchemist.connected = True
        self.alchemist.build_engine()

        create_engine_mock.assert_not_called()
        ask_credentials_mock.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_2(self, create_engine_mock, ask_credentials_mock):
        """Verify build_engine() calls create_engine() with engine string.
        """
        self.alchemist.username = "user"
        self.alchemist.password = "pass"
        self.alchemist.has_credentials = False

        self.alchemist.build_engine()

        ask_credentials_mock.assert_called()
        login_string = "mysql+pymysql://user:pass@localhost/"
        create_engine_mock.assert_called_with(login_string, echo=False)
   
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.validate_database")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_3(self, create_engine_mock, validate_database_mock): 
        """Verify build_engine() calls create_engine() with db engine string.
        """
        self.alchemist.username = "user"
        self.alchemist.password = "pass"
        self.alchemist.database = "database"

        self.alchemist.build_engine()

        login_string = "mysql+pymysql://user:pass@localhost/"
        db_login_string = "mysql+pymysql://user:pass@localhost/database"

        create_engine_mock.assert_any_call(login_string, echo=False)
        create_engine_mock.assert_any_call(db_login_string, echo=False)
        
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_4(self, create_engine_mock, ask_credentials_mock):
        """Verify build_engine() sets has_credentials.
        """
        self.alchemist.has_credentials = True
        self.alchemist.connected = False
        self.alchemist.metadata = "Test"
        self.alchemist.graph = "Test"

        self.alchemist.build_engine()

        self.alchemist.connected = True
        self.assertEqual(self.alchemist.metadata, None)
        self.assertEqual(self.alchemist.graph, None)

    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_5(self, create_engine_mock):
        """Verify AlchemyHandler echo property controls create_engine() parameters.
        """
        self.alchemist.username = "user"
        self.alchemist.password = "pass"
        self.alchemist.build_engine()

        login_string = "mysql+pymysql://user:pass@localhost/"
    
        create_engine_mock.assert_any_call(login_string, echo=False)

        self.alchemist.echo = True
        self.alchemist.connected = False
        self.alchemist.build_engine()

        create_engine_mock.assert_any_call(login_string, echo=True)


    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_1(self, BuildEngine, AskDatabase, AskCredentials):
        """Verify connect() returns if build_engine() does not complain.
        """
        self.alchemist.has_credentials = True
        self.alchemist.connected = True
        self.alchemist.connect()
        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_2(self, BuildEngine, AskDatabase, AskCredentials):
        """Verify connect() AlchemyHandler properties control function calls.
        """ 
        self.alchemist.connected = True
        self.alchemist.connected_database = True
        self.alchemist.has_credentials = True
        self.alchemist.connect(ask_database=True)
        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
                                                        "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_3(self, BuildEngine, AskDatabase, AskCredentials):
        """Verify connect() depends on build_engine() to raise ValueError.
        """
        self.alchemist.connected = False
        BuildEngine.side_effect = OperationalError("", "", "")
        
        with self.assertRaises(ValueError):
            self.alchemist.connect()

        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_called()

    def build_engine_side_effect(self, mock_engine):
        """Helper function for side effect usage.
        """
        self.alchemist._engine = mock_engine

    @patch("pdm_utils.classes.alchemyhandler.querying.execute")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_execute_1(self, build_engine_mock, execute_mock):
        """Verify function structure of execute().
        """
        engine_mock = Mock()
        build_engine_mock.side_effect = self.build_engine_side_effect(
                                                                engine_mock)

        self.alchemist.execute("Executable")

        execute_mock.assert_called_with(engine_mock, "Executable", 
                                                            return_dict=True)
   
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_scalar_1(self, build_engine_mock):
        """Verify function structure of scalar().
        """
        engine_mock = Mock() 
        proxy_mock  = Mock()
        
        engine_mock.execute.return_value = proxy_mock
        proxy_mock.scalar.return_value = "Scalar"

        build_engine_mock.side_effect = self.build_engine_side_effect(
                                                               engine_mock)
       
        self.alchemist.scalar("Executable")

        engine_mock.execute.assert_called_with("Executable")
        proxy_mock.scalar.assert_called()
        build_engine_mock.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.first_column")
    def test_first_column_1(self, first_column_mock):
        """Verify function structure of first_row().
        """
        engine_mock = Mock()
        self.alchemist._engine = engine_mock

        self.alchemist.first_column("Executable")

        first_column_mock.assert_called_with(engine_mock, "Executable")

    @patch("pdm_utils.classes.alchemyhandler.MetaData")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_metadata_1(self, AskDatabase, BuildEngine, MetaData):
        """Verify build_metadata() relies on AlchemyHandler properties.
        """
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
        """Verify build_metadata() calls ask_database() and build_engine().
        """
        self.alchemist.has_database = True
        self.alchemist.connected = True
        
        self.alchemist.build_metadata()

        AskDatabase.assert_not_called()
        BuildEngine.assert_not_called()
        MetaData.assert_called() 

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_table_1(self, BuildMetadata, TranslateTable):
        """Verify translate_table() calls parsing.translate_table().
        """
        self.alchemist.metadata = "Metadata"

        self.alchemist.translate_table("Test")

        TranslateTable.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_table_2(self, BuildMetadata, TranslateTable):
        """Verify translate_table() calls build_metadata().
        """
        self.alchemist.metadata = None

        self.alchemist.translate_table("Test")

        TranslateTable.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_column_1(self, BuildMetadata, TranslateColumn):
        """Verify translate_column() calls parsing.translate_column().
        """
        self.alchemist.metadata = "Metadata"

        self.alchemist.translate_column("Test")

        TranslateColumn.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.parsing.translate_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_translate_column_2(self, BuildMetadata, TranslateColumn):
        """Verify translate_column() calls build_metadata().
        """
        self.alchemist.metadata = None

        self.alchemist.translate_column("Test")

        TranslateColumn.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_table_1(self, BuildMetadata, GetTable):
        """Verify get_table() calls querying.get_column().
        """
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_table("Test")

        GetTable.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_table")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_table_2(self, BuildMetadata, GetTable):
        """Verify get_table() calls build_metadata().
        """
        self.alchemist.metadata = None

        self.alchemist.get_table("Test")

        GetTable.assert_called_with(None, "Test")
        BuildMetadata.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.get_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_column_1(self, BuildMetadata, GetColumn):
        """Verify get_column() calls querying.get_column().
        """
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_column("Test")

        GetColumn.assert_called_with("Metadata", "Test")
        BuildMetadata.assert_not_called()
        
    @patch("pdm_utils.classes.alchemyhandler.querying.get_column") 
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_column_2(self, BuildMetadata, GetColumn):
        """Verify get_column() calls build_metadata()
        """
        self.alchemist.metadata = None

        self.alchemist.get_column("Test")

        GetColumn.assert_called_with(None, "Test")
        BuildMetadata.assert_called()
 
    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_1(self, BuildMetadata, BuildGraph):
        """Verify get_column() calls querying.build_graph().
        """
        BuildGraph.return_value = "Graph"

        self.alchemist.metadata = "Metadata"

        self.alchemist.build_graph()

        BuildMetadata.assert_not_called()
        BuildGraph.assert_called_with("Metadata")
        
        self.assertEqual(self.alchemist.graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_2(self, BuildMetadata, BuildGraph):
        """Verify get_column() calls build_metadata().
        """
        BuildGraph.return_value = "Graph"

        self.alchemist.metadata = None

        self.alchemist.build_graph()

        BuildMetadata.assert_called()
        BuildGraph.assert_called_with(None)
        
        self.assertEqual(self.alchemist.graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.cartography.get_map")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_map_1(self, BuildMetadata, GetMap): 
        """Verify get_map() calls cartography.get_map()
        """
        self.alchemist.metadata = "Metadata"

        self.alchemist.get_map("Test")

        BuildMetadata.assert_not_called()
        GetMap.assert_called_with("Metadata", "Test")

    @patch("pdm_utils.classes.alchemyhandler.cartography.get_map")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_get_map_2(self, BuildMetadata, GetMap):
        """Verify get_map() calls build_metadata().
        """
        self.alchemist.metadata = None

        self.alchemist.get_map("Test")

        BuildMetadata.assert_called()
        GetMap.assert_called_with(None, "Test") 

if __name__ == "__main__":
    unittest.main()
