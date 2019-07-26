"""Unit tests for the MySQLConnectionHandler class.
The tests assume the host is 'localhost' and there is a
'anonymous' user with 'anonymous' password set up."""

import unittest
from classes import MySQLConnectionHandler
import getpass


class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler.MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pwd = "anonymous"
        self.invalid_user = "invalid"
        self.invalid_pwd = "invalid"
        self.valid_db = "Actino_Draft"


    def test_set_database(self):
        """Verify the database is set correctly."""
        self.handler.set_database(self.valid_db)
        self.assertEqual(self.handler.database, self.valid_db)


    def test_set_username(self):
        """Verify the username is set correctly."""
        self.handler.set_username(self.valid_user)
        self.assertEqual(self.handler.username, self.valid_user)

    def test_set_password(self):
        """Verify the password is set correctly."""
        self.handler.set_password(self.valid_pwd)
        self.assertEqual(self.handler.password, self.valid_pwd)

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

    def test_create_connection_1(self):
        """Verify that valid credentials produce no error."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pwd
        self.handler.database = self.valid_db
        self.handler.connection = False
        self.handler.create_connection()
        self.assertTrue(self.handler.connection)




    # TODO develop more unit tests for create_connection() method and
    # for other methods.
    #
    # def test_create_connection_unchecked_bad_creds_good_db(self):
    #     """."""
    #     # Should be prompted to input valid credentials
    #     self.handler.set_username("anonymous")
    #     self.handler.set_password("anonymus")
    #     self.handler.set_database("Actino_Draft")
    #     self.handler.create_connection()
    #     with self.subTest():
    #         self.assertEqual(self.handler.database, "Actino_Draft")
    #         if self.handler.attempts_remaining > 0:
    #             self.assertTrue(self.handler.have_connection)
    #
    # def test_create_connection_unchecked_good_creds_bad_db(self):
    #     """."""
    #     # Should notify bad database for user
    #     self.handler.set_username("anonymous")
    #     self.handler.set_password("anonymous")
    #     self.handler.set_database("Actino_Daft")
    #     self.handler.create_connection()
    #     with self.subTest():
    #         self.assertEqual(self.handler.database, "Actino_Daft")
    #         self.assertEqual(self.handler.attempts_remaining, 3)
    #         self.assertEqual(self.handler.username, "anonymous")
    #         self.assertEqual(self.handler.password, "anonymous")
    #
    # def test_create_connection_checked_good_creds_good_db(self):
    #     """."""
    #     self.handler.set_username("anonymous")
    #     self.handler.set_password("anonymous")
    #     self.handler.set_database("Actino_Draft")
    #     self.handler.validate_credentials()
    #     self.handler.create_connection()
    #     with self.subTest():
    #         self.assertEqual(self.handler.database, "Actino_Draft")
    #         self.assertEqual(self.handler.attempts_remaining, 3)
    #         self.assertEqual(self.handler.username, "anonymous")
    #         self.assertEqual(self.handler.password, "anonymous")
    #         self.assertTrue(self.handler.have_connection)
    #         self.assertTrue(self.handler.valid_credentials)
    #         self.assertTrue(self.handler.valid_database)
    #
    # def test_close_connection(self):
    #     """."""
    #     self.handler.set_username("anonymous")
    #     self.handler.set_password("anonymous")
    #     self.handler.set_database("Actino_Draft")
    #     self.handler.validate_credentials()
    #     self.handler.create_connection()
    #     with self.subTest():
    #         self.assertEqual(self.handler.database, "Actino_Draft")
    #         self.assertEqual(self.handler.attempts_remaining, 3)
    #         self.assertEqual(self.handler.username, "anonymous")
    #         self.assertEqual(self.handler.password, "anonymous")
    #         self.assertTrue(self.handler.have_connection)
    #         self.assertTrue(self.handler.valid_credentials)
    #         self.assertTrue(self.handler.valid_database)
    #         self.handler.close_connection()
    #     with self.subTest():
    #         self.assertEqual(self.handler.database, "Actino_Draft")
    #         self.assertEqual(self.handler.attempts_remaining, 3)
    #         self.assertEqual(self.handler.username, "anonymous")
    #         self.assertEqual(self.handler.password, "anonymous")
    #         self.assertEqual(self.handler.database, "Actino_Draft")
    #         self.assertEqual(self.handler.attempts_remaining, 3)
    #         self.assertEqual(self.handler.username, "anonymous")
    #         self.assertEqual(self.handler.password, "anonymous")


if __name__ == '__main__':
    unittest.main()
