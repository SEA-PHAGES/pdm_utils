"""Integration tests for the MySQLConnectionHandler class.
The tests assume the host is 'localhost' and there is a
'anonymous' user with 'anonymous' password set up."""

import unittest
from classes import MySQLConnectionHandler



class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler.MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pass = "anonymous"
        self.valid_db = "Actino_Draft"
        self.invalid_user = "invalid"
        self.invalid_pass = "invalid"
        self.invalid_db = "invalid"
        self.attempts = 5

    def test_validate_credentials_1(self):
        """Verify that valid credentials are validated properly."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.validate_credentials()
        self.assertTrue(self.handler.credential_status)

    def test_validate_credentials_2(self):
        """Verify that valid user and invalid password are invalid."""
        self.handler.username = self.valid_user
        self.handler.password = self.invalid_user
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_credentials_3(self):
        """Verify that invalid user and valid password are invalid."""
        self.handler.username = self.invalid_user
        self.handler.password = self.valid_pass
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_credentials_4(self):
        """Verify that invalid user and invalid password are invalid."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_user
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_database_access_1(self):
        """Verify that valid credentials and valid db produce no error."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler.validate_database_access()
        self.assertTrue(self.handler._database_status)

    def test_validate_database_access_2(self):
        """Verify that valid credentials produce "error" when credential
        status is False."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.validate_database_access()
        self.assertFalse(self.handler._database_status)

    # TODO disabled since it requires user response.
    # def test_get_credentials_1(self):
    #     """If credential_status is False, should ask for user/pass up to 3
    #     times."""
    #     self.handler.username = self.valid_user
    #     self.handler.password = self.valid_pass
    #     self.handler.get_credentials()
    #     self.assertTrue(self.handler.credential_status)
    #
    # def test_get_credentials_2(self):
    #     """If credential_status is True, should do nothing even is invalid
    #     credentials are used."""
    #     self.handler.username = self.invalid_user
    #     self.handler.password = self.invalid_pass
    #     self.handler.credential_status = True
    #     self.handler.get_credentials()
    #     self.assertTrue(self.handler.credential_status)

    # def test_ask_user_and_pass(self):
    #     """Credentials should be asked, login_attempts should decrease by 1"""
    #     self.handler.login_attempts = self.attempts
    #     self.handler.ask_username_and_password()
    #     with self.subTest():
    #         self.assertEqual(self.handler.username, self.valid_user)
    #     with self.subTest():
    #         self.assertEqual(self.handler.password, self.valid_pass)
    #     with self.subTest():
    #         self.assertEqual(self.handler.login_attempts, self.attempts - 1)

    def test_connection_status_1(self):
        """Connection status should be False if connection is None"""
        self.handler.connection = None
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())

    def test_close_connection_1(self):
        """Close connection should fail if connection is None"""
        self.handler.close_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())

    def test_open_connection_1(self):
        """Should create new connection readily if valid credentials and
        database are used and appropriate flags are set."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)

    def test_connection_status_2(self):
        """Connection status should be True with valid open connection"""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())

    def test_close_connection_2(self):
        """Close connection should close connection if valid one exists"""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        self.handler.close_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())

    def test_close_connection_3(self):
        """Close connection should fail to close connection if valid
        connection has already been closed."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        self.handler.close_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        self.handler.close_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())

    def test_open_connection_2(self):
        """Connection should not be made if valid credentials are used with
        an invalid database - unvalidated valid credentials."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.invalid_db
        self.handler._database_status = False
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)

    def test_open_connection_3(self):
        """Connection should not be made if valid credentials are used with
        an invalid database - validated valid credentials."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.invalid_db
        self.handler.credential_status = True
        self.handler._database_status = False
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)

    # TODO disabled since it requires user response.
    # def test_open_connection_4(self):
    #     """If invalid credentials are used with a valid database, should be
    #     prompted for correct credentials - unvalidated invalid credentials."""
    #     self.handler.username = self.invalid_user
    #     self.handler.password = self.invalid_pass
    #     self.handler.database = self.valid_db
    #     self.handler._database_status = True
    #     self.handler.open_connection()
    #     with self.subTest():
    #         self.assertIsNotNone(self.handler.connection)

    def test_open_connection_5(self):
        """Should fail to open connection if one already exists."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())

    def test_open_connection_6(self):
        """Should succeed in opening a new connection if valid connection
        was closed."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        self.handler.close_connection()
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())







































### TODO Old integration tests. Can probably be removed.

#
# class TestMySQLConnectionHandler(unittest.TestCase):
#     def setUp(self):
#         self.handler = MySQLConnectionHandler.MySQLConnectionHandler()
#         self.valid_user = "anonymous"
#         self.valid_pwd = "anonymous"
#         self.invalid_user = "invalid"
#         self.invalid_pwd = "invalid"
#         self.valid_db = "Actino_Draft"
#
#
#
#
#     def test_validate_credentials_1(self):
#         """Verify that valid credentials are validated."""
#         self.handler.username = self.valid_user
#         self.handler.password = self.valid_pwd
#         self.handler.host = "localhost"
#         self.handler.validate_credentials()
#         self.assertTrue(self.handler.valid_credentials)
#
#     def test_validate_credentials_2(self):
#         """Verify that invalid username produces an error."""
#         self.handler.username = self.invalid_user
#         self.handler.password = self.valid_pwd
#         self.handler.host = "localhost"
#         self.handler.valid_credentials = True
#         self.handler.validate_credentials()
#         self.assertFalse(self.handler.valid_credentials)
#
#     def test_validate_credentials_3(self):
#         """Verify that invalid password produces an error."""
#         self.handler.username = self.valid_user
#         self.handler.password = self.invalid_pwd
#         self.handler.host = "localhost"
#         self.handler.valid_credentials = True
#         self.handler.validate_credentials()
#         self.assertFalse(self.handler.valid_credentials)
#
#
#
#
#     def test_validate_database_access_checked_1(self):
#         """Verify that valid credentials produce no error."""
#         self.handler.username = self.valid_user
#         self.handler.password = self.valid_pwd
#         self.handler.database = self.valid_db
#         self.handler.valid_credentials = True
#         self.handler.validate_database_access()
#         self.assertTrue(self.handler.valid_database)
#
#     def test_validate_database_access_checked_2(self):
#         """Verify that valid credentials produce an error when the
#         valid_credentials flag is False."""
#         self.handler.username = self.valid_user
#         self.handler.password = self.valid_pwd
#         self.handler.database = self.valid_db
#         self.handler.valid_credentials = False
#         self.handler.validate_database_access()
#         self.assertFalse(self.handler.valid_database)
#
#
#     def test_validate_database_access_checked_3(self):
#         """Verify that invalid username produces an error when the
#         valid_credentials flag is True."""
#         self.handler.username = self.invalid_user
#         self.handler.password = self.valid_pwd
#         self.handler.database = self.valid_db
#         self.handler.valid_credentials = True
#         self.handler.validate_database_access()
#         self.assertFalse(self.handler.valid_database)
#
#
#
#
#     # TODO not working correctly since there is a connection.open attribute problem.
#     def test_create_connection_1(self):
#         """Verify that valid credentials produce no error."""
#         self.handler.username = self.valid_user
#         self.handler.password = self.valid_pwd
#         self.handler.database = self.valid_db
#         self.handler.create_connection()
#         self.assertTrue(self.handler.connection)




if __name__ == '__main__':
    unittest.main()
