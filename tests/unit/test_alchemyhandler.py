import unittest
from unittest.mock import Mock
from unittest.mock import patch

from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError

from pdm_utils.classes.alchemyhandler import (
                    AlchemyHandler, SQLCredentialsError, MySQLDatabaseError)


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
        self.assertEqual(self.alchemist._metadata, None)
        self.assertEqual(self.alchemist._graph, None)
        self.assertEqual(self.alchemist._session, None)

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

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.clear")
    def test_username_3(self, clear_mock):
        """Verify changing usrename property calls clear().
        """
        self.alchemist.username = "Test"

        clear_mock.assert_called()

    def test_password_1(self):
        """Verify password property conserves has_credentials and connected.
        """
        self.alchemist.password = "Test"
        self.assertFalse(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    def test_password_2(self):
        """Verify password property sets has_credentials with valid username.
        """
        self.alchemist.username = "Test"
        self.alchemist.password = "Test"
        self.assertTrue(self.alchemist.has_credentials)
        self.assertFalse(self.alchemist.connected)

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.clear")
    def test_password_3(self, clear_mock):
        """Verify changing password property calls clear().
        """
        self.alchemist.password = "Test"

        clear_mock.assert_called()

    def test_construct_engine_string_1(self):
        """Verify construct_engine_string generates an expected URI.
        """
        URI = self.alchemist.construct_engine_string(username="pdm_user",
                                                     password="pdm_pass")
        self.assertEqual(URI, "mysql+pymysql://pdm_user:pdm_pass@localhost/")

    def test_construct_engine_string_2(self):
        """Verify construct_engine_string accepts use of different drivers.
        """
        URI = self.alchemist.construct_engine_string(driver="mysqlconnector",
                                                     username="pdm_user",
                                                     password="pdm_pass")

        self.assertEqual(URI,
                         "mysql+mysqlconnector://pdm_user:pdm_pass@localhost/")

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

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_engine_3(self, build_engine_mock):
        """Verify engine property calls build_engine() selectively.
        """
        mock_engine = Mock()
        build_engine_mock.return_value = mock_engine

        self.alchemist._engine = "Test"
        self.assertEqual(self.alchemist.engine, "Test")

        build_engine_mock.assert_not_called()

        self.alchemist._engine = None
        self.alchemist.engine

        build_engine_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler"
           ".extract_engine_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.get_mysql_dbs")
    def test_engine_4(self, get_mysql_dbs_mock,
                      extract_engine_credentials_mock):
        """Verify call structure of engine property setter.
        """
        mock_engine = Mock(spec=Engine)

        self.alchemist.engine = mock_engine

        get_mysql_dbs_mock.assert_called()
        extract_engine_credentials_mock.assert_called_with(mock_engine)

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_metadata_1(self, build_metadata_mock):
        """Verify metadata property calls build_metadata() selectively.
        """
        self.alchemist._metadata = "Test"
        self.alchemist.metadata

        build_metadata_mock.assert_not_called()

        self.alchemist._metadata = None
        self.alchemist.metadata

        build_metadata_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_graph")
    def test_graph_1(self, build_graph_mock):
        """Verify graph property calls build_graph() selectively.
        """
        self.alchemist._graph = "Test"
        self.alchemist.graph

        build_graph_mock.assert_not_called()

        self.alchemist._graph = None
        self.alchemist.graph

        build_graph_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_session")
    def test_session_1(self, build_session_mock):
        """Verify session property calls build_session() selectively.
        """
        self.alchemist._session = "Test"
        self.alchemist.session

        build_session_mock.assert_not_called()

        self.alchemist._session = None
        self.alchemist.session

        build_session_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_mapper")
    def test_mapper_1(self, build_mapper_mock):
        """Verify mapper property calls build_mapper() selectively.
        """
        self.alchemist._mapper = "Test"
        self.alchemist.mapper

        build_mapper_mock.assert_not_called()

        self.alchemist._mapper = None
        self.alchemist.mapper

        build_mapper_mock.assert_called()

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
        mock_engine = Mock()
        mock_proxy = Mock()

        mock_engine.execute.return_value = mock_proxy
        mock_proxy.fetchall.return_value = [("pdm_test_db",),
                                            ("Actino_Draft",)]

        self.alchemist.connected = True
        self.alchemist.database = "pdm_test_db"
        self.alchemist._engine = mock_engine

        self.alchemist.validate_database()

        mock_engine.execute.assert_called_once()
        mock_proxy.fetchall.assert_called()

    def test_validate_database_2(self):
        """Verify validate_database() raises IndexError without database.
        """
        self.alchemist.connected = True

        with self.assertRaises(AttributeError):
            self.alchemist.validate_database()

    def test_validate_database_3(self):
        """Verify validate_database() raises ValueError from bad database input.
        """
        mock_engine = Mock()
        mock_proxy = Mock()

        mock_engine.execute.return_value = mock_proxy
        mock_proxy.fetchall.return_value = []

        self.alchemist.connected = True
        self.alchemist.database = "test db"
        self.alchemist._engine = mock_engine

        with self.assertRaises(MySQLDatabaseError):
            self.alchemist.validate_database()

        mock_engine.execute.assert_called_once()
        mock_proxy.fetchall.assert_called()

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
        """Verify build_engine() raises attribute error without credentials.
        """
        self.alchemist.username = "user"
        self.alchemist.password = "pass"
        self.alchemist.has_credentials = False

        with self.assertRaises(AttributeError):
            self.alchemist.build_engine()

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
        self.alchemist._metadata = "Test"
        self.alchemist._graph = "Test"

        self.alchemist.build_engine()

        self.alchemist.connected = True
        self.assertEqual(self.alchemist._metadata, None)
        self.assertEqual(self.alchemist._graph, None)

    @patch("pdm_utils.classes.alchemyhandler.sqlalchemy.create_engine")
    def test_build_engine_5(self, create_engine_mock):
        """Verify AlchemyHandler echo property controls create_engine()
        parameters.
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
    def test_connect_1(self, build_engine_mock, ask_database_mock,
                       AskCredentials):
        """Verify connect() returns if build_engine() does not complain.
        """
        self.alchemist.has_credentials = True
        self.alchemist.connected = True
        self.alchemist.connect()
        build_engine_mock.assert_called()
        ask_database_mock.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
           "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_2(self, build_engine_mock, ask_database_mock,
                       AskCredentials):
        """Verify connect() AlchemyHandler properties control function calls.
        """
        self.alchemist.connected = True
        self.alchemist.connected_database = True
        self.alchemist.has_credentials = True
        self.alchemist.connect(ask_database=True)
        build_engine_mock.assert_called()
        ask_database_mock.assert_not_called()
        AskCredentials.assert_not_called()

    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler."
           "ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_3(self, build_engine_mock, ask_database_mock,
                       AskCredentials):
        """Verify connect() depends on build_engine() to raise ValueError.
        """
        self.alchemist.connected = False
        build_engine_mock.side_effect = OperationalError("", "", "")

        with self.assertRaises(SQLCredentialsError):
            self.alchemist.connect()

        build_engine_mock.assert_called()
        ask_database_mock.assert_not_called()
        AskCredentials.assert_called()

    def build_engine_side_effect(self, mock_engine):
        """Helper function for side effect usage.
        """
        self.alchemist._engine = mock_engine

    @patch("pdm_utils.classes.alchemyhandler.MetaData")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_metadata_1(self, ask_database_mock, build_engine_mock,
                              metadata_mock):
        """Verify build_metadata() relies on AlchemyHandler properties.
        """
        self.alchemist.has_database = False
        self.alchemist.connected = False

        self.alchemist.build_metadata()

        ask_database_mock.assert_called()
        build_engine_mock.assert_called()
        metadata_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.MetaData")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_metadata_2(self, ask_database_mock, build_engine_mock,
                              metadata_mock):
        """Verify build_metadata() calls ask_database() and build_engine().
        """
        self.alchemist.has_database = True
        self.alchemist.connected = True
        self.alchemist.connected_database = True

        self.alchemist.build_metadata()

        ask_database_mock.assert_not_called()
        build_engine_mock.assert_not_called()
        metadata_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_1(self, build_metadata_mock, build_graph_mock):
        """Verify build_graph() calls querying.build_graph().
        """
        build_graph_mock.return_value = "Graph"

        self.alchemist._metadata = "Metadata"

        self.alchemist.build_graph()

        build_metadata_mock.assert_not_called()
        build_graph_mock.assert_called_with("Metadata")

        self.assertEqual(self.alchemist._graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.querying.build_graph")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_graph_2(self, build_metadata_mock, build_graph_mock):
        """Verify build_graph() calls build_metadata().
        """
        build_graph_mock.return_value = "Graph"

        self.alchemist._metadata = None

        self.alchemist.build_graph()

        build_metadata_mock.assert_called()
        build_graph_mock.assert_called_with(None)

        self.assertEqual(self.alchemist._graph, "Graph")

    @patch("pdm_utils.classes.alchemyhandler.sessionmaker")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_session_1(self, ask_database_mock, build_engine_mock,
                             sessionmaker_mock):
        """Verify build_session() relies on AlchemyHandler properties.
        """
        self.alchemist.has_database = False
        self.alchemist.connected = False

        self.alchemist.build_session()

        ask_database_mock.assert_called()
        build_engine_mock.assert_called()
        sessionmaker_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.sessionmaker")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    def test_build_session_2(self, ask_database_mock, build_engine_mock,
                             sessionmaker_mock):
        """Verify build_session() calls ask_database() and build_engine().
        """
        self.alchemist.has_database = True
        self.alchemist.connected = True

        self.alchemist.build_session()

        ask_database_mock.assert_not_called()
        build_engine_mock.assert_not_called()
        sessionmaker_mock.assert_called()

    @patch("pdm_utils.classes.alchemyhandler.automap_base")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_mapper_1(self, build_metadata_mock, automap_base_mock):
        """Verify build_mapper() calls automap_base().
        """
        base_mock = Mock()
        automap_base_mock.return_value = base_mock

        self.alchemist._metadata = "Metadata"

        self.alchemist.build_mapper()

        build_metadata_mock.assert_not_called()
        automap_base_mock.assert_called_with(metadata="Metadata")

        self.assertEqual(self.alchemist._mapper, base_mock)

    @patch("pdm_utils.classes.alchemyhandler.automap_base")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_metadata")
    def test_build_mapper_2(self, build_metadata_mock, automap_base_mock):
        """Verify build_mapper() calls build_metadata().
        """
        base_mock = Mock()
        automap_base_mock.return_value = base_mock

        self.alchemist._metadata = None

        self.alchemist.build_mapper()

        build_metadata_mock.assert_called()
        automap_base_mock.assert_called_with(metadata=None)

        self.assertEqual(self.alchemist._mapper, base_mock)


if __name__ == "__main__":
    unittest.main()
