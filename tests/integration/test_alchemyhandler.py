import sys
import unittest
from pathlib import Path

from networkx import Graph
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.engine.base import Engine
from sqlalchemy.ext.declarative.api import DeclarativeMeta
from sqlalchemy.orm.session import Session

from pdm_utils.classes.alchemyhandler import (
                                AlchemyHandler, MySQLDatabaseError)

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

# pdm_anon, pdm_anon, and pdm_test_db
user = test_db_utils.USER
pwd = test_db_utils.PWD
db = test_db_utils.DB


class TestAlchemyHandler(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    def setUp(self):
        self.alchemist = AlchemyHandler()

    def test_validate_database_1(self):
        """Verify validate_database() detects good database access.
        """
        self.alchemist.username = user
        self.alchemist.password = pwd
        self.alchemist.build_engine()

        self.alchemist.database = db
        self.alchemist.validate_database()

    def test_validate_database_2(self):
        """Verify validate_database() detects bad database access.
        """
        self.alchemist.username = user
        self.alchemist.password = pwd
        self.alchemist.build_engine()

        self.alchemist.database = "not_database"
        with self.assertRaises(MySQLDatabaseError):
            self.alchemist.validate_database()

    def test_build_engine_1(self):
        """Verify build_engine() creates and stores Engine object.
        """
        self.alchemist.username = user
        self.alchemist.password = pwd
        self.alchemist.build_engine()

        self.assertTrue(isinstance(self.alchemist.engine, Engine))
        self.assertTrue(self.alchemist.connected)

    def test_build_engine_2(self):
        """Verify build_engine() connects to database if has_database.
        """
        self.alchemist.username = user
        self.alchemist.password = pwd
        self.alchemist.database = db

        self.alchemist.build_engine()

        self.assertTrue(self.alchemist.connected_database)

    def connect_to_pdm_test_db(self):
        """Sets alchemist credentials and database to connect to pdm_test_db.
        """
        self.alchemist.username = user
        self.alchemist.password = pwd
        self.alchemist.database = db

    def test_build_metadata_1(self):
        """Verify build_metadata() creates and stores MetaData and Engine.
        """
        self.connect_to_pdm_test_db()

        self.assertTrue(isinstance(self.alchemist.metadata, MetaData))
        self.assertTrue(isinstance(self.alchemist.engine, Engine))

        self.assertTrue(self.alchemist._graph is None)
        self.assertTrue(self.alchemist._session is None)
        self.assertTrue(self.alchemist._mapper is None)

    def test_build_graph_1(self):
        """Verify build_graph() creates and stores, Graph, MetaData, and Engine.
        """
        self.connect_to_pdm_test_db()

        self.assertTrue(isinstance(self.alchemist.graph, Graph))
        self.assertTrue(isinstance(self.alchemist.metadata, MetaData))
        self.assertTrue(isinstance(self.alchemist.engine, Engine))

        self.assertTrue(self.alchemist._session is None)
        self.assertTrue(self.alchemist._mapper is None)

    def test_build_session_1(self):
        """Verify build_session() creates and stores Engine and Session
        """
        self.connect_to_pdm_test_db()

        self.assertTrue(isinstance(self.alchemist.session, Session))
        self.assertTrue(isinstance(self.alchemist.engine, Engine))

        self.assertTrue(self.alchemist._graph is None)
        self.assertTrue(self.alchemist._metadata is None)
        self.assertTrue(self.alchemist._mapper is None)

    def test_build_mapper_1(self):
        """Verify build_mapper() creates and stores, Mapper, MetaData, and Engine.
        """
        self.connect_to_pdm_test_db()

        self.assertTrue(isinstance(self.alchemist.mapper, DeclarativeMeta))
        self.assertTrue(isinstance(self.alchemist.metadata, MetaData))
        self.assertTrue(isinstance(self.alchemist.engine, Engine))

        self.assertTrue(self.alchemist._session is None)
        self.assertTrue(self.alchemist._graph is None)

    def test_engine_1(self):
        """Verify AlchemyHandler extracts credentials from engine.
        """
        engine = create_engine(self.alchemist.construct_engine_string(
                                                    username=user,
                                                    password=pwd,
                                                    database=db))

        self.alchemist.engine = engine

        self.assertEqual(self.alchemist.username, user)
        self.assertEqual(self.alchemist.password, pwd)
        self.assertEqual(self.alchemist.database, db)

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()


if __name__ == "__main__":
    unittest.main()
