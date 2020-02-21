import sqlalchemy
import pymysql
from getpass import getpass
from networkx import Graph
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.engine.base import Engine
from sqlalchemy.exc import OperationalError
from pdm_utils.functions import querying
from pdm_utils.functions import cartography

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
            login_string = ("mysql+pymysql://"
                           f"{self.username}:{self.password}@localhost")
            if self.has_database:
                login_string = login_string + f"/{self.database}"

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

        for table in self.metadata.tables.keys():
            if table.lower() == raw_table.lower():
                return table

        raise ValueError(f"Table '{raw_table}' requested to be filtered "
                         f"is not in '{self.sql_handle.database}'")

    def translate_column(self, raw_column):
        parsed_column = querying.parse_column(raw_column)

        table = self.translate_table([0])
        table_obj = self[table]["table"]
        for column in table_obj.columns.keys():
            if column.lower() == raw_column.lower():
                return column

        raise ValueError(f"Field '{raw_column}' requested to be filtered"
                         f" is not in '{table_object.name}'")

    def get_table(self, table): 
        table = self.translate_table(table)
        return self.metadata.tables[table]

    def get_column(self, column):
        parsed_column = querying.parse_column(column)
        table = self.get_table(parsed_column[0])

        columns_dict = dict(table.columns)
        column = columns_dict[parsed_column[1]]

        return column

    def get_map(self, template):
        if not self.metadata:
            self.build_metadata()

        return cartography.get_map(self.metadata, template)

    def build_graph(self):
        if not self.metadata:
            self.build_metadata()
        
        graph = Graph()
        for table in self.metadata.tables.keys():
            table_object = self.metadata.tables[table]
            graph.add_node(table_object.name, table=table_object) 

        for target_table in self.metadata.tables.keys():
            target_table_obj  = self.metadata.tables[target_table]
        
            for foreign_key in target_table_obj.foreign_keys:   
                referent_column = foreign_key.column
                referent_table_obj = referent_column.table
                referent_table = referent_table_obj.name

                graph.add_edge(target_table, referent_table, 
                                             key=foreign_key)

        self.graph = graph
        
        return True

    def build_where(self, filter_expression):
        if not self.graph:
            self.build_graph()
         
        return querying.build_whereclause(self.graph, filter_expression)

    def build_select(self, columns, where=None, order_by=None):
        if not self.graph:
            self.build_graph()

        return querying.build_select(self.graph, columns,
                                                        where=where,
                                                        order_by=order_by)

    def build_count(self, columns, where=None, order_by=None):
        if not self.graph:
            self.build_graph()

        return querying.build_count(self.graph, columns,
                                                        where=where,
                                                        order_by=order_by)

    def build_distinct(self, columns, where=None, order_by=None):
        if not self.graph:
            self.build_graph()

        return querying.build_distinct(self.graph, columns,
                                                        where=where,
                                                        order_by=order_by)

    def build_session(self):
        if not self.has_database:
            raise
        if not self.connected:
            self.build_engine()
            

        session_maker = sessionmaker()
        self.session = session_maker(bind=self.engine)
        return 
  
