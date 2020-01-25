import sqlalchemy
import pymysql
from getpass import getpass
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError

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

    def build_session(self):
        if not self.has_database:
            return False
        if not self.connected:
            if not self.build_engine():
                return False

        session_maker = sessionmaker()
        self.session = session_maker(bind=self.engine)
        return True

    def connect(self, database=None):
        self.database = database

        attempts = 0
        connected = False
        while(not connected and attempts < self.login_attempts):
            self.ask_credentials()
            connected = self.build_engine()
            attempts += 1

        if not connected:
            print("Maximum logout attempts reached.\n"
                  "Please check your credentials and try again")
            exit(1)

if __name__ == "__main__":
    alchemist = AlchemyHandler()
