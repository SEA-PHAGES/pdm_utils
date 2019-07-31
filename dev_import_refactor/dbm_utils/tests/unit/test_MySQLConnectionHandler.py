"""Unit tests for the MySQLConnectionHandler class.
The tests assume the host is 'localhost' and there is a
'anonymous' user with 'anonymous' password set up."""

import unittest
from classes import MySQLConnectionHandler


class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler.MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pwd = "anonymous"
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





    # TODO develop more unit tests for create_connection() method and
    # for other methods.
    #
    #
    #
    # TODO below are Christian's tests.
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
