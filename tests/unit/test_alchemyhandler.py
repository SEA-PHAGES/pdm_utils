from sqlalchemy import MetaData
from sqlalchemy.engine.base import Engine
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.schemagraph import SchemaGraph
from unittest.mock import patch, Mock, PropertyMock
import unittest

class TestAlchemyHandler(unittest.TestCase):
    def setUp(self):
        self.alchemist = AlchemyHandler()

    def test_constructor_1(self):
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
    
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_credentials")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.ask_database")
    @patch("pdm_utils.classes.alchemyhandler.AlchemyHandler.build_engine")
    def test_connect_1(self, BuildEngine, AskDatabase, AskCredentials):
        self.alchemist.username = "Test"
        self.alchemist.password = "Test"
        self.alchemist.connect()
        BuildEngine.assert_called()
        AskDatabase.assert_not_called()
        AskCredentials.assert_not_called()
         

if __name__ == "__main__":
    unittest.main()
