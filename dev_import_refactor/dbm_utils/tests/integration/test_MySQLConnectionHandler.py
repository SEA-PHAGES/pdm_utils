"""Integration tests for the MySQLConnectionHandler class.
The tests assume the host is 'localhost' and there is a
'anonymous' user with 'anonymous' password set up."""

import unittest
from classes import MySQLConnectionHandler




class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler.MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pwd = "anonymous"
        self.invalid_user = "invalid"
        self.invalid_pwd = "invalid"
        self.valid_db = "Actino_Draft"




    def test_validate_credentials_1(self):
        """Verify that valid credentials are validated."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pwd
        self.handler.host = "localhost"
        self.handler.validate_credentials()
        self.assertTrue(self.handler.valid_credentials)

    def test_validate_credentials_2(self):
        """Verify that invalid username produces an error."""
        self.handler.username = self.invalid_user
        self.handler.password = self.valid_pwd
        self.handler.host = "localhost"
        self.handler.valid_credentials = True
        self.handler.validate_credentials()
        self.assertFalse(self.handler.valid_credentials)

    def test_validate_credentials_3(self):
        """Verify that invalid password produces an error."""
        self.handler.username = self.valid_user
        self.handler.password = self.invalid_pwd
        self.handler.host = "localhost"
        self.handler.valid_credentials = True
        self.handler.validate_credentials()
        self.assertFalse(self.handler.valid_credentials)




    def test_validate_database_access_checked_1(self):
        """Verify that valid credentials produce no error."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pwd
        self.handler.database = self.valid_db
        self.handler.valid_credentials = True
        self.handler.validate_database_access()
        self.assertTrue(self.handler.valid_database)

    def test_validate_database_access_checked_2(self):
        """Verify that valid credentials produce an error when the
        valid_credentials flag is False."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pwd
        self.handler.database = self.valid_db
        self.handler.valid_credentials = False
        self.handler.validate_database_access()
        self.assertFalse(self.handler.valid_database)


    def test_validate_database_access_checked_3(self):
        """Verify that invalid username produces an error when the
        valid_credentials flag is True."""
        self.handler.username = self.invalid_user
        self.handler.password = self.valid_pwd
        self.handler.database = self.valid_db
        self.handler.valid_credentials = True
        self.handler.validate_database_access()
        self.assertFalse(self.handler.valid_database)




    # TODO not working correctly since there is a connection.open attribute problem.
    def test_create_connection_1(self):
        """Verify that valid credentials produce no error."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pwd
        self.handler.database = self.valid_db
        self.handler.create_connection()
        self.assertTrue(self.handler.connection)




if __name__ == '__main__':
    unittest.main()
