import sys

import sqlalchemy
import pymysql
from getpass import getpass
from networkx import Graph
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError
from sqlalchemy.ext.automap import automap_base

from pdm_utils.functions import cartography
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import parsing

#-----------------------------------------------------------------------------
#GLOBAL VARIABLES
CREDENTIALS_MSG = ("Credentials invalid and maximum login attempts reached. "
                   "Please check your MySQL credentials and try again.")

DATABASE_MSG = ("Unable to connect to database with valid credentials.\n"
                "Please check your MySQL database access, "
                "and/or your database availability.")

class AlchemyHandler:
    def __init__(self, database=None, username=None, password=None):
        self._database = database
        self._username = username
        self._password = password
        #Username, password, and database are included for two reasons:
        #The ability to set credentials 'pythonically'
        #An understanding of whether the credentials are valid for the engine.

        self._engine = None
        self._metadata = None
        self._graph = None
        self._mapper = None
        self._session = None
        self.echo = False

        self.connected = False
        self.has_credentials = False

        if username != None and password != None:
            self.has_credentials = True

        self.has_database = False
        self.connected_database = False
        self._databases = []

        if database != None:
            self.has_database = True



#-----------------------------------------------------------------------------
#ALCHEMYHANDLER PROPERTIES
#-----------------------------------------------------------------------------

    @property
    def database(self):
        """Returns the AlchemyHandler's set database.

        :returns: Returns a copy of the database attribute.
        :rtype: str
        """
        database = self._database
        return database

    @database.setter
    def database(self, database):
        if database == None:
            self.has_database = False
            self.connected_database = False
            return

        if not isinstance(database, str):
            raise TypeError("Entered database name is not a string.")

        self._database = database

        self.has_database = True
        self.connected_database = False

    @property
    def username(self):
        """Returns the AlchemynHandler's set username.

        :returns: Returns a copy of the username attribute.
        :rtype: str
        """
        username = self._username
        return username

    @username.setter
    def username(self, username):
        if username == None:
            self.has_credentials = False
            self.connected = False
            return

        if not isinstance(username, str):
            raise TypeError("Entered username is not a string")

        self._username = username

        if self._password != None:
            self.has_credentials = True
        self.connected = False

        self.clear()

    @property
    def password(self):
        """Returns the AlchemyHandler's set password.

        :returns: Returns a copy of the password attribute.
        :rtype: str
        """
        password = self._password
        return password

    @password.setter
    def password(self, password):
        if password == None:
            self.has_credentials = False
            self.connected = False
            return

        if not isinstance(password, str):
            raise TypeError("Entered password is not a string.")

        self._password = password

        if self._username != None:
            self.has_credentials = True
        self.connected = False

        self.clear()

    @property
    def login_attempts(self):
        """Returns the AlchemyHandler's number of login attempts for login.

        :returns: Returns the number of login attempts for login.
        :rtype: str
        """
        login_attempts = self._login_attempts
        return login_attempts

    @property
    def engine(self):
        """Returns the AlchemyHandler's stored engine object.

        :returns: Returns the AlchemyHandler's stored engine object.
        :rtype: Engine
        """
        if self._engine is None:
            self.build_engine()

        engine = self._engine
        return engine

    @engine.setter
    def engine(self, engine):
        if engine is None:
            self.connected = False
            self._engine = None
            self.clear()
            return

        if not isinstance(engine, Engine):
            raise TypeError

        self.clear()
        self._engine = engine
        self.get_mysql_dbs()
        self.connected = True

    @property
    def session(self):
        """Returns the AlchemyHandler's stored session object.
        """
        if self._session is None:
            self.build_session()

        session = self._session
        return session

    @property
    def metadata(self):
        """Returns the AlchemyHandler's stored metadata object.

        :returns: Returns the AlchemyHandler's stored engine object.
        :rtype: MetaData
        """
        if self._metadata is None:
            self.build_metadata()

        metadata = self._metadata
        return metadata

    @property
    def graph(self):
        """Returns the AlchemyHandler's stored graph object.

        :returns: Returns the AlchemyHandler's stored metadata graph object.
        :rtype: Graph
        """
        if self._graph is None:
            self.build_graph()

        graph = self._graph
        return graph

    @property
    def mapper(self):
        """Returns the AlchemyHandler's stored automapper object.

        :returns: Returns the AlchemyHandler's stored mapper object.
        """
        if self._mapper is None:
            self.build_mapper()

        mapper = self._mapper
        return mapper

    @property
    def databases(self):
        """Returns a copy of the databases available to the current credentials

        :returns: Returns the AlchemyHandler's available databases.
        :rtype: list[str]
        """
        databases = self._databases.copy()
        return databases

#-----------------------------------------------------------------------------
#CONNECTION METHODS
#-----------------------------------------------------------------------------

    def ask_database(self):
        """Ask for database input to store in AlchemyHandler.
        """
        self._database = input("MySQL database: ")

        self.has_database = True
        self.connected = False

    def ask_credentials(self):
        """Ask for username and password input to store in AlchemyHandler.
        """
        self._username = getpass(prompt="MySQL username: ")
        self._password = getpass(prompt="MySQL password: ")

        self.has_credentials = True
        self.connected = False

    def validate_database(self):
        """Validate access to database using stored MySQL credentials.
        """
        if not self.connected:
            raise ValueError("AlchemyHandler currently not connected to MySQL.")
        if not self.has_database:
            raise AttributeError("No database in AlchemyHandler to validate")

        if not self._databases:
            self.get_mysql_dbs()

        if self._database not in self._databases:
            raise ValueError("User does not have access to "
                            f"database {self._database}")

    def build_engine(self):
        """Create and store SQLAlchemy Engine object.
        """
        if not self.connected:
            if not self.has_credentials:
                raise AttributeError("AlchemyHandler missing credentials.\n"
                                     "Cannot connect to MySQL.")

            login_string = self.construct_engine_string(
                                            username=self._username,
                                            password=self._password)

            self._engine = sqlalchemy.create_engine(login_string,
                                                            echo=self.echo)
            self._engine.connect()

            self.connected = True

            self._metadata = None
            self._graph = None

            self.get_mysql_dbs()

        if self.has_database:
            database = self._database

            self.validate_database()

            login_string = self.construct_engine_string(
                                        username=self._username,
                                        password=self._password,
                                        database=self._database)

            self._engine = sqlalchemy.create_engine(login_string,
                                                            echo=self.echo)
            self._engine.connect()

            self.connected_database = True

    def connect(self, ask_database=False, login_attempts=5, pipeline=False):
        """Ask for input to connect to MySQL and MySQL databases.

        :param ask_database: Toggle whether to connect to a database.
        :type ask_database: Boolean
        :param login_attempts: Set number of total login attempts.
        :type login_attempts: int
        """
        attempts = 0
        if self.has_credentials:
            try:
                self.build_engine()
            except:
                pass

        while(not self.connected and attempts < login_attempts):
            attempts += 1
            self.ask_credentials()

            try:
                self.build_engine()
            except:
                pass

        if not self.connected:
            if not pipeline:
                raise ValueError(CREDENTIALS_MSG)
            else:
                print(CREDENTIALS_MSG)
                sys.exit(1)

        if ask_database:
            try:
                self.build_engine()
            except:
                pass

            while(not self.connected_database and attempts < login_attempts):
                attempts += 1
                self.ask_database()

                try:
                    self.build_engine()
                except:
                    pass

            if not self.connected_database:
                if not pipeline:
                    raise ValueError(DATABASE_MSG)
                else:
                    print(DATABASE_MSG)
                    sys.exit(1)

    def construct_engine_string(self, db_type="mysql", driver="pymysql",
                            username="", password="", database=""):
        """Construct a SQLAlchemy engine URL.

        :param db_type: Type of SQL database.
        :type db_type: str
        :param driver: Name of the Python DBAPI used to connect.
        :type driver: str
        :param username: Username to login to SQL database.
        :type username: str
        :param password: Password to login to SQL database.
        :type password: str
        :param database: Name of the database to connect to.
        :type database: str
        :returns: URL string to create SQLAlchemy engine.
        :rtype: str
        """
        engine_string = f"{db_type}+{driver}://{username}:{password}@localhost/{database}"
        return engine_string

    def get_mysql_dbs(self):
        """Retrieve database names from MySQL.

        :returns: List of database names.
        :rtype: list
        """
        if not self.connected:
            raise ValueError("AlchemyHandler currently not connected to MySQL.")
        databases = mysqldb_basic.get_mysql_dbs(self._engine)
        self._databases = list(databases)

    def clear(self):
        """Clear properties tied to MySQL credentials/database.
        """
        self._engine = None
        self._metadata = None
        self._graph = None
        self._mapper = None

        if not self._session is None:
            try:
                self._session.close()
            except:
                pass

        self._session = None

        self._databases = []

#-----------------------------------------------------------------------------
#SQLALCHEMY-RELATED OBJECT GENERATION METHODS

    def build_metadata(self):
        """Create and store SQLAlchemy MetaData object.
        """
        if not self.has_database:
            self.ask_database()

        if not self.connected:
            self.build_engine()

        self._metadata = MetaData(bind=self._engine)
        self._metadata.reflect()

    def build_session(self):
        """Create and store SQLAlchemy Session object.
        """
        if not self.has_database:
            self.ask_database()

        if not self.connected:
            self.build_engine()

        if not self._session is None:
            try:
                self._session.close()
            except:
                pass

        session_maker_obj = sessionmaker(bind=self._engine)
        self._session = session_maker_obj()

    def build_graph(self):
        """Create and store SQLAlchemy MetaData related NetworkX Graph object.
        """
        if self._metadata is None:
            self.build_metadata()

        self._graph = querying.build_graph(self._metadata)

    def build_mapper(self):
        """Create and store SQLAlchemy automapper Base object.
        """
        if self._metadata is None:
            self.build_metadata()

        self._mapper = automap_base(metadata=self._metadata)
        self._mapper.prepare()

    def build_all(self):
        """Create and store all relevant SQLAlchemy objects.
        """
        self.build_session()
        self.build_graph()
        self.build_mapper()

#-----------------------------------------------------------------------------
#SQLALCHEMY QUALITY-OF-LIFE FUNCTIONS

    def get_map(self, table):
        """Get SQLAlchemy ORM map object.
        """
        if self._mapper is None:
            self.build_mapper()

        table = parsing.translate_table(self._metadata, table)
        return self._mapper.classes[table]
