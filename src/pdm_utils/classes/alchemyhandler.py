import sqlalchemy
import pymysql
from getpass import getpass
from networkx import Graph
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError
from pdm_utils.functions import cartography
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb
from pdm_utils.functions import parsing

class AlchemyHandler:
    def __init__(self, database=None, username=None, password=None):
        self._database = database
        self._username = username
        self._password = password

        self._engine = None
        self.metadata = None
        self.graph = None
        self.session = None
            
        self.connected = False
        self.connected_database = False
        self.has_credentials = False
        self.has_database = False
  
        if database != None:
            self.has_database = True

        if username != None and password != None:
            self.has_credentials = True

    @property
    def database(self):
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
        """Returns the DatabaseHandler's set username.
        :return username:
            Returns a copy of the username attribute.
        :type username: str
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

    @property
    def password(self):
        """Returns the DatabaseHandler's set password.
        :return password:
            Returns a copy of the password attribute.
        :type password: str
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

    @property
    def login_attempts(self):
        login_attempts = self._login_attempts
        return login_attempts

    @property
    def engine(self):
        engine = self._engine
        return engine

    @engine.setter
    def engine(self, engine):
        if engine is None:
            self.connected = False
            self._engine = None
            return 

        if not isinstance(engine, Engine):
            raise TypeError

        self._engine = engine
        self.connected = True

    @property
    def tables(self):
        if not self.metadata:
            self.build_metadata()
           
        return self.metadata.tables

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

        proxy = self._engine.execute("SHOW DATABASES")   

        results = proxy.fetchall()

        databases = []
        for result in results:
            databases.append(result[0])

        if self._database not in databases:
            raise ValueError("User does not have access to "
                            f"database {self._database}")

    def build_engine(self):
        """Create and store SQLAlchemy Engine object.
        """
        if not self.connected:
            if not self.has_credentials:
                self.ask_credentials()

            login_string = mysqldb.construct_engine_string(
                                            username=self._username,
                                            password=self._password)

            self._engine = sqlalchemy.create_engine(login_string)
            self._engine.connect()
          
            self.connected = True

            self.metadata = None
            self.graph = None

        if self.has_database:
            database = self._database

            self.validate_database()

            login_string = mysqldb.construct_engine_string(
                                        username=self._username,
                                        password=self._password,
                                        database=self._database)

            self._engine = sqlalchemy.create_engine(login_string)
            self._engine.connect()
            
            self.connected_database = True 
        
    def connect(self, ask_database=False, login_attempts=5):
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
            raise ValueError("Maximum login attempts reached.  Please check "
                             "your MySQL credentials and try again.")
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
                raise ValueError("Maximum database login attempts reached. "
                                 "Please your database access and try again.")
        
    def execute(self, executable, return_dict=True):
        """Use SQLAlchemy Engine to execute a MySQL query.

        :param executable: Input a executable MySQL query.
        :type executable: Select
        :type executable: str
        :param return_dict: Toggle whether execute returns dict or tuple.
        :type return_dict: Boolean
        :returns: Results from execution of given MySQL query.
        :rtype: list[dict]
        :rtype: list[tuple]
        """
        if self.engine is None:
            self.build_engine()

        proxy = self.engine.execute(executable)

        results = proxy.fetchall()

        if return_dict:
            results_dicts = []
            for result in results:
                results_dicts.append(dict(result))

            results = results_dicts 

        return results

    def scalar(self, executable):
        """Use SQLAlchemy Engine to execute a MySQL query.

        :param executable: Input a executable MySQL query.
        :type executable: Select
        :type executable: str
        :param return_dict: Toggle whether execute returns dict or tuple.
        :type return_dict: Boolean
        :returns: Results from execution of given MySQL query.
        :rtype: list[dict]
        :rtype: list[tuple]
        """
        if self.engine is None:
            self.build_engine()

        proxy = self.engine.execute(executable)

        scalar = proxy.scalar()

        return scalar

    def build_metadata(self):
        """Create and store SQLAlchemy MetaData object.
        """
        if not self.has_database:
            self.ask_database()

        if not self.connected:
            self.build_engine()
        
        self.metadata = MetaData(bind=self.engine)
        self.metadata.reflect()
        
        return True

    #TODO Clean dependancies and remove
    def translate_table(self, raw_table): 
        if not self.metadata:
            self.build_metadata()

        return parsing.translate_table(self.metadata, raw_table)  

    #TODO Clean dependancies and remove
    def translate_column(self, raw_column):
        if not self.metadata:
            self.build_metadata()

        return parsing.translate_column(self.metadata, raw_column) 

    #TODO Clean dependancies and remove
    def get_table(self, table): 
        if not self.metadata:
            self.build_metadata() 

        return querying.get_table(self.metadata, table)

    #TODO Clean dependancies and remove
    def get_column(self, column):
        if not self.metadata:
            self.build_metadata() 

        return querying.get_column(self.metadata, column) 

    def build_graph(self):
        """Create and store SQLAlchemy MetaData related NetworkX Graph object.
        """
        if not self.metadata:
            self.build_metadata()
        
        self.graph = querying.build_graph(self.metadata)

    #For when necessary
    #def build_session(self):
    #    if not self.has_database:
    #        self.connect(ask_database=True)
    #    if not self.connected:
    #        self.build_engine()
            

    #    session_maker = sessionmaker()
    #    self.session = session_maker(bind=self.engine)
    #    return 
    
    #TODO Clean dependancies and remove
    def get_map(self, template):
        if not self.metadata:
            self.build_metadata()

        return cartography.get_map(self.metadata, template)


