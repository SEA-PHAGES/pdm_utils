import sqlalchemy
import pymysql
from getpass import getpass
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError
from pdm_utils.functions import cartography
from pdm_utils.classes.schemagraph import SchemaGraph
from pdm_utils.functions import querying

class AlchemyHandler:
    def __init__(self, database=None, username=None, password=None, 
                 attempts=5, verbose=False):

        self._database = database
        self._username = username
        self._password = password
        self._login_attempts = attempts

        self._engine = None
        self.metadata = None
        self.graph = None
        self.session = None

        self.connected = False
        self.has_database = False
        self.has_credentials = False
   
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

        if self.password:
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

        if self.username:
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

        self._engine = engine

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
            return self.engine

        if not self.has_credentials:
            self.ask_credentials()

        try:
            login_string = ("mysql+pymysql://"
                           f"{self.username}:{self.password}@localhost")
            if self.has_database:
                login_string = login_string + f"/{self.database}"

            self.engine = sqlalchemy.create_engine(login_string)
            self.engine.connect()
            
            self.connected = True

            self.metadata = None
            self.graph = None
            return self.engine

        except OperationalError: 
            return False

    def build_metadata(self):
        if not self.has_database:
            return False

        if not self.connected:
            if not self.build_engine():
                return False

        self.metadata = MetaData(bind=self.engine)
        self.metadata.reflect()
        return True 

    def get_map(self, template):
        if not self.metadata:
            if not self.build_metadata():
                return None

        return cartography.get_map(self.metadata, template)

    def build_schemagraph(self):
        if not self.metadata:
            if not self.build_metadata():
                return False

        graph = SchemaGraph()
        graph.setup(self.metadata)
        self.graph = graph
        
        return True

    def build_session(self):
        if not self.has_database:
            return False
        if not self.connected:
            if not self.build_engine():
                return False

        session_maker = sessionmaker()
        self.session = session_maker(bind=self.engine)
        return True

    def connect(self, database=True):
        if database:
            self.ask_database()
        else:
            self.database = None

        if not self.has_credentials:
            self.ask_credentials()

        attempts = 1
        connected = self.build_engine()
        while(not connected and attempts < self.login_attempts):
            connected = self.build_engine()
            attempts += 1
            self.ask_credentials()

        if not connected:
            print("Maximum logout attempts reached.\n"
                  "Please check your credentials and try again")
            exit(1)

    def select(self, columns, where_clause=None, order_by_clause=None):
        if not self.graph:
            if not self.build_schemagraph():
                return None

        return querying.build_select(self.graph, columns,
                                            where_clause=where_clause,
                                            order_by_clause=order_by_clause)

    def count(self, columns, where_clause=None, order_by_clause=None):
        if not self.graph:
            if not self.build_schemagraph():
                return None

        return querying.build_count(self.graph, columns,
                                            where_clause=where_clause,
                                            order_by_clause=order_by_clause)

    def distinct(self, columns, where_clause=None, order_by_clause=None):
        if not self.graph:
            if not self.build_schemagraph():
                return None

        return querying.build_distinct(self.graph, columns,
                                            where_clause=where_clause,
                                            order_by_clause=order_by_clause)

    def execute(self, excutable, return_dict=True):
        if not self.engine:
            if not self.build_engine():
                return None

        proxy = self.engine.execute(executable)

        results = proxy.fetchall()

        if return_dict:
            results_dicts = []
            for result in results:
                results_dicts.append(dict(result))

            results = results_dicts 

        return results

    def show_tables(self):
        if not self.graph:
            if not self.build_schemagraph():
                return None

        return self.graph.show_tables()

    def get_table_node(self, table):
        if not self.graph:
            if not self.build_schemagraph():
                return None
        
        return self.graph.get_table(table)
        
    def get_table(self, table):
        if not self.graph:
            if not self.build_schemagraph():
                return None

        table_node = self.graph.get_table(table)
        
        alchemy_table = None

        if table_node:
            alchemy_table = table_node.table

        return alchemy_table
        
