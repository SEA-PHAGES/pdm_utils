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
                                              "MySQLConnectionHandler."}

        # Store pymysql connection object once it's created
        self.connection = None
        self.connection_open = False
        # Variables needed to establish a connection
        self.username = username
        self.password = password
        self.database = database

        # How many login attempts can the user make?
        self.attempts_remaining = attempts

        # valid_credentials flag lets us know whether user/pass work
        # valid_database flag lets us know whether user has access to db
        self.valid_credentials = False
        self.valid_database = False

    def set_username(self, value):
        """
        Sets username attribute to input value.
        :param value:
        :return:
        """
        self.username = value

    def set_password(self, value):
        """
        Sets password attribute to input value.
        :param value:
        :return:
        """
        self.password = value

    def set_database(self, value):
        """
        Sets database attribute to input value.
        :param value:
        :return:
        """
        self.database = value

    # TODO unit test.
    def get_credential_status(self):
        """
        Return True if username and password are valid, False if not.
        :return:
        """
        return self.valid_credentials

    # TODO unit test.
    def get_credentials(self):
        """
        Gives user x attempts to input correct username and password,
        assuming valid credentials aren't already had.
        :return:
        """
        while self.valid_credentials is False:
            if self.attempts_remaining > 0:
                self.ask_username_and_password()
                self.validate_credentials()
                # If self.valid_credentials is still false, bad user/pass
                if self.valid_credentials is False:
                    if self.attempts_remaining >= 1:
                        print(self.messages["bad user or pass"].format(
                            self.username))
                    else:
                        print(self.messages["too many attempts"])
                        return


    # TODO unit test.
    def create_connection(self):
        """
        If a connection doesn't already exist, attempts to create one.
        :return:
        """
        # If a connection doesn't already exist
        if self.connection is None or self.connection_open is False:
            # If credentials have already been validated
            if self.valid_credentials is True:
                # Test database
                self.validate_database_access()
                # If database is valid
                if self.valid_database is True:
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
                if self.valid_credentials is False:
                    # Prompt user for credentials with x attempts
                    self.get_credentials()
                # If credentials are now valid
                if self.valid_credentials is True:
                    # Test database
                    self.validate_database_access()
                    # If database is valid
                    if self.valid_database is True:
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

    # TODO unit test.
    def get_connection_status(self):
        """
        Returns True if the connection is open. False if the connection
        has been closed.
        :return:
        """
        if self.connection is not None:
            return self.connection_open

    # TODO unit test.
    def get_connection(self):
        """
        Returns the pymysql connection object
        :return:
        """
        return self.connection

    # TODO unit test.
    def close_connection(self):
        """
        Calls close method on pymysql connection object
        :return:
        """
        if self.connection:
            self.connection.close()
            self.connection = None

    # TODO unit test.
    def ask_username_and_password(self):
        """
        Reverse increments (decrements? don't think that's a word...)
        the number of attempts remaining, then prompts the user for a
        MySQL username and password.
        :return:
        """
        self.attempts_remaining -= 1
        self.username = getpass.getpass(prompt="MySQL username: ")
        self.password = getpass.getpass(prompt="MySQL password: ")
        return

    def validate_credentials(self):
        """
        Tries to connect to MySQL localhost using the verified username
        and password. Successful connection triggers setting successful
        login flag to True. Otherwise, flag persists at False.
        :return:
        """
        try:
            con = pms.connect("localhost", self.username, self.password)
            con.close()
            self.valid_credentials = True
        except pms.err.Error:
            self.valid_credentials = False
        return

    def validate_database_access(self):
        """
        Tries to connect to the specified database using the verified
        username and password. If username and password haven't been
        verified yet it returns without doing anything.
        :return:
        """
        if self.valid_credentials is True:
            try:
                con = pms.connect("localhost", self.username, self.password,
                                  self.database)
                con.close()
                self.valid_database = True
            except pms.err.Error:
                self.valid_database = False
        return
