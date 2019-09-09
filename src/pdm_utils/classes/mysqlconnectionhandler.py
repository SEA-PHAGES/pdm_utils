import pymysql as pms
import getpass


class MySQLConnectionHandler:
    def __init__(self, username=None, password=None, database=None,
                 attempts=3):
        """
        This object is intended to handle creation of a connection to
        a MySQL database.  It gives a user 3 attempts to input valid
        MySQL credentials, and tests the database name to verify it
        exists.  If valid credentials are given and the database is
        real, a connection is established and can be used externally.
        :param database: the name of the MySQL database to connect to
        """
        # Error messages for user
        self.messages = {"bad user or pass": "Access denied for user '{}'@"
                                             "'localhost'. Please verify your "
                                             "account credentials and try "
                                             "again.",
                         "too many attempts": "Too many login attempts. Please"
                                              " verify your account "
                                              "credentials and try again.",
                         "bad database": "User '{}'@'localhost' does not have "
                                         "access to database '{}'. Please "
                                         "verify that the database name is "
                                         "spelled correctly and that your "
                                         "account has access to it.",
                         "already connected": "Already connected to MySQL. "
                                              "Multiple connections are not "
                                              "presently supported through a "
                                              "single instance of "
                                              "MySQLConnectionHandler.",
                         "no connection": "Cannot proceed without an open "
                                          "connection.",
                         "invalid syntax": "Invalid MySQL syntax for "
                                           "issued command '{}'.",
                         "validate login": "\nError: database validation can "
                                           "only be performed after "
                                           "validation of username and "
                                           "password.",
                         "database change": "\nError: database attribute "
                                            "cannot be changed with a "
                                            "connection open. Create a new "
                                            "MySQLConnectionHandler object "
                                            "to use database '{}'\n"}

        # Store pymysql connection object once it's created
        self.connection = None

        # Variables needed to establish a connection
        # Used by username, password, and database @properties
        self._username = username
        self._password = password
        self._database = database

        # How many login attempts can the user make?
        # Used by login_attempts @property
        self._login_attempts = attempts

        # Flag to track whether the user and password have been validated
        # Used by credential_status @property
        self._credential_status = False

        # Flag to track whether the database has been validated
        self._database_status = False
        
    @property
    def username(self):
        """
        Returns self.username
        :return: self.username
        """
        # Unit test passed
        return self._username

    @username.setter
    def username(self, value):
        """
        Sets self.username. Also resets self.credential_status to False.
        :param value: value to attach to self.username variable
        :return:
        """
        # Unit test passed
        self._username = value
        self._credential_status = False

    @property
    def password(self):
        """
        Returns self.password
        :return: self.password
        """
        # Unit test passed
        return self._password

    @password.setter
    def password(self, value):
        """
        Sets self.password. Also resets self.credential_status to False.
        :param value: value to attach to self.password variable
        :return:
        """
        # Unit test passed
        self._password = value
        self._credential_status = False

    @property
    def database(self):
        """
        Returns self.database
        :return: self.database
        """
        # Unit test passed
        return self._database

    @database.setter
    def database(self, value):
        """
        If self._database hasn't been set yet, sets it. Also resets
        self._database_status to False.
        :param value: value to attach to self.database variable
        :return:
        """
        # Unit test passed
        if self.connection is None:
            self._database = value
            self._database_status = False
        else:
            print(self.messages["database change"].format(value))
        
    def validate_database_access(self):
        """
        Tries to connect to the specified database using the verified
        username and password. If _username and _password haven't been
        verified yet it returns without doing anything.
        :return:
        """
        # Integration test passed
        if self.credential_status is True:
            try:
                con = pms.connect("localhost", self.username, self.password,
                                  self.database)
                con.close()
                self._database_status = True
            except pms.err.Error:
                self._database_status = False
        else:
            print(self.messages["validate login"])
        return

    @property
    def login_attempts(self):
        """
        Returns self.attempts_remaining
        :return: self.attempts_remaining
        """
        # Unit test passed
        return self._login_attempts

    @login_attempts.setter
    def login_attempts(self, value):
        """
        Sets self.attempts_remaining
        :param value: value to attach to self.attempts_remaining variable
        :return:
        """
        # Unit test passed
        self._login_attempts = value

    @property
    def credential_status(self):
        """
        Returns self.credential_status
        :return: self.credential_status
        """
        # Unit test passed
        return self._credential_status

    @credential_status.setter
    def credential_status(self, value):
        """
        Sets self.credential_status
        :param value: value to attach to self.credential_status
        :return:
        """
        # Unit tests passed
        # Make sure value is Boolean before setting the flag
        if value is True or value is False:
            self._credential_status = value
        else:
            print("{} is not an allowed type for this Boolean flag".format(
                type(value)))
            print("Using False")
            self._credential_status = False

    def get_credentials(self):
        """
        Gives user x attempts to input correct _username and _password,
        assuming valid credentials aren't already had.
        :return:
        """
        # Integration tests passed
        # If current credentials are unverified or
        while self.credential_status is False:
            if self.login_attempts > 0:
                self.ask_username_and_password()
                self.validate_credentials()
                # If self.credential_status is still false, bad user/pass
                if self.credential_status is False:
                    if self.login_attempts >= 1:
                        print(self.messages["bad user or pass"].format(
                            self.username))
                    else:
                        print(self.messages["too many attempts"])
                        return

    def ask_username_and_password(self):
        """
        Reverse increments (decrements? don't think that's a word...)
        the number of attempts remaining, then prompts the user for a
        MySQL _username and _password.
        :return:
        """
        # Integration test passed
        self.login_attempts = self.login_attempts - 1
        self.username = getpass.getpass(prompt="MySQL username: ")
        self.password = getpass.getpass(prompt="MySQL password: ")
        return

    def validate_credentials(self):
        """
        Tries to connect to MySQL localhost using the verified _username
        and _password. Successful connection triggers setting successful
        login flag to True. Otherwise, flag persists at False.
        :return:
        """
        # Integration tests passed
        try:
            con = pms.connect("localhost", self.username, self.password)
            con.close()
            self.credential_status = True
        except pms.err.Error:
            self.credential_status = False
        return

    def connection_status(self):
        """
        Returns True if the connection is open. False if the connection
        has been closed.
        :return:
        """
        # Integration tests (some tests joint with close_connection) passed
        if self.connection is not None:
            return self.connection.open
        else:
            return False

    def open_connection(self):
        """
        If connection status is False, open a new connection.
        :return:
        """
        # Integration tests passed
        # If a connection doesn't already exist (or if existing connection
        # was closed)
        if self.connection_status() is False:
            # If credentials have already been validated
            if self.credential_status is True:
                # Test database
                self.validate_database_access()
                # If database is valid
                if self._database_status is True:
                    # Create connection
                    self.connection = pms.connect("localhost",
                                                  self.username,
                                                  self.password,
                                                  self.database)
                # If database is invalid
                else:
                    # Print bad database message
                    print(self.messages["bad database"].format(self.username,
                                                               self.database))
            # If credentials have not been validated
            else:
                # Test them
                self.validate_credentials()
                # If tested credentials are not valid
                if self.credential_status is False:
                    # Prompt user for credentials with x attempts
                    self.get_credentials()
                # If credentials are now valid
                if self.credential_status is True:
                    # Test database
                    self.validate_database_access()
                    # If database is valid
                    if self._database_status is True:
                        # Create connection
                        self.connection = pms.connect("localhost",
                                                      self.username,
                                                      self.password,
                                                      self.database)
                    # If database is invalid
                    else:
                        # Print bad database message
                        print(self.messages["bad database"].format(
                            self.username, self.database))
                # If credentials are still invalid
                else:
                    # Don't create connection, leave flag False, return
                    return
        # If a connection does already exist
        else:
            # Print already connected message
            print(self.messages["already connected"])
        return

    def execute_query(self, query=""):
        """
        If connection exists and is open, attempts to attach a DictCursor
        to the connection and execute the input query. If connection doesn't
        exist or has been closed, tries to open a new one before creating
        DictCursor and executing query.
        :param query:
        :return: results list if available, None otherwise
        """
        # Integration tests passed
        if self.connection_status() is True:
            try:
                cursor = self.connection.cursor(pms.cursors.DictCursor)
                cursor.execute(query)
                results = cursor.fetchall()
                cursor.close()
                return results
            except pms.err.Error as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                return
        else:
            self.open_connection()
            if self.connection_status() is True:
                try:
                    cursor = self.connection.cursor(pms.cursors.DictCursor)
                    cursor.execute(query)
                    results = cursor.fetchall()
                    cursor.close()
                    return results
                except pms.err.Error as err:
                    print("Error {}: {}".format(err.args[0], err.args[1]))
                    return
            else:
                print(self.messages["no connection"])
                return

    def execute_transaction(self, statement_list=[]):
        """
        If connection exists and is open, attempts to attach a cursor to
        the connection and execute the commands in the input list. If
        connection doesn't exist or has been closed, tries to open a new one
        before creating DictCursor and executing query.
        :param statement_list: a list of any number of MySQL statements with
        no expectation that anything will return
        :return: 0 or 1 status code. 0 means no problems, 1 means problems
        """
        # TODO: integration tests
        if self.connection_status() is True:
            try:
                cursor = self.connection.cursor()
                cursor.execute("START TRANSACTION")
                for statement in statement_list:
                    try:
                        print(statement)
                        cursor.execute(statement)
                    except pms.err.ProgrammingError:
                        print(self.messages["invalid syntax"].format(
                            statement))
                        cursor.execute("ROLLBACK")
                        cursor.close()
                        return 1
                cursor.execute("COMMIT")
                cursor.close()
                return 0
            except pms.err.Error as err:
                print("Error {}: {}".format(err.args[0], err.args[1]))
                return 1
        else:
            self.open_connection()
            if self.connection_status() is True:
                try:
                    cursor = self.connection.cursor()
                    cursor.execute("START TRANSACTION")
                    for statement in statement_list:
                        try:
                            print(statement)
                            cursor.execute(statement)
                        except pms.err.ProgrammingError:
                            print(self.messages["invalid syntax"].format(
                                statement))
                            cursor.execute("ROLLBACK")
                            cursor.close()
                            return 1
                    cursor.execute("COMMIT")
                    cursor.close()
                    return 0
                except pms.err.Error as err:
                    print("Error {}: {}".format(err.args[0], err.args[1]))
                    return 1
            else:
                print(self.messages["no connection"])
                return 1

    def close_connection(self):
        """
        Calls close method on pymysql connection object
        :return:
        """
        # Integration tests passed
        if self.connection_status() is True:
            self.connection.close()
            self.connection = None
