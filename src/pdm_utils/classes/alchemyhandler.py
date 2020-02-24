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
        self.has_database = False
        self.has_credentials = False
  
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
            self.connected = False
            return

        if not isinstance(database, str):
            raise TypeError("Entered database name is not a string.")

        self._database = database

        self.has_database = True
        self.connected = False

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
        if engine == None:
            self.connected = False
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
        self._database = input("MySQL database: ")

        self.has_database = True
        self.connected = False

    def ask_credentials(self):
        self._username = getpass(prompt="MySQL username: ")
        self._password = getpass(prompt="MySQL password: ")

        self.has_credentials = True
        self.connected = False

    def build_engine(self):
        if self.connected:
            return

        if not self.has_credentials:
            self.ask_credentials()

        try: 
            database = ""
            if self.has_database:
                database = self._database

            login_string = mysqldb.construct_engine_string(
                                            username=self._username,
                                            password=self._password,
                                            database=database)

            self._engine = sqlalchemy.create_engine(login_string)
            self._engine.connect()
            
            self.connected = True

            self.metadata = None
            self.graph = None
            return

        except OperationalError: 
            raise

    def connect(self, ask_database=False, login_attempts=5):
        if ask_database:
            if not self.has_database:
                self.ask_database()

        if not self.has_credentials:
            self.ask_credentials()

        attempts = 1
        try:
            self.build_engine()
            self.connected=True
            return
        except:
            pass

        while(not self.connected and attempts < login_attempts):
            try:
                self.build_engine()
            except:
                pass

            attempts += 1
            self.ask_credentials()

        if not self.connected:
            print("Maximum logout attempts reached.\n"
                  "Please check your credentials and try again")
            exit(1)

    def execute(self, executable, return_dict=True):
        if not self.engine:
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
        if not self.engine:
            self.build_engine()

        proxy = self.engine.execute(executable)

        scalar = proxy.scalar()

        return scalar

    def build_metadata(self):
        if not self.has_database:
            self.ask_database()

        if not self.connected:
            self.build_engine()

        self.metadata = MetaData(bind=self.engine)
        self.metadata.reflect()
        return True 

    def translate_table(self, raw_table): 
        if not self.metadata:
            self.build_metadata()

        return parsing.translate_table(self.metadata, raw_table)  

    def translate_column(self, raw_column):
        if not self.metadata:
            self.build_metadata()

        return parsing.translate_column(self.metadata, raw_column) 

    def get_table(self, table): 
        if not self.metadata:
            self.build_metadata() 

        return querying.get_table(self.metadata, table)

    def get_column(self, column):
        if not self.metadata:
            self.build_metadata() 

        return querying.get_column(self.metadata, column) 

    def build_graph(self):
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
    
    def get_map(self, template):
        if not self.metadata:
            self.build_metadata()

        return cartography.get_map(self.metadata, template)


