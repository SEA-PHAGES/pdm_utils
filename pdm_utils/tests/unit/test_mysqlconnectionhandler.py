"""
Unit tests for the MySQLConnectionHandler class.
Assumptions: 'anonymous'@'localhost' identified by 'anonymous' exists
with 'anonymous' password and has
select granted on all tables of all databases.
"""

import unittest
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler


class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pass = "anonymous"
        self.valid_db = "Actino_Draft"
        self.attempts = 5

    def test_database_property_1(self):
        """Verify database is set/get correctly."""
        self.handler.database = self.valid_db
        with self.subTest():
            self.assertEqual(self.handler.database, self.valid_db)
        with self.subTest():
            self.assertFalse(self.handler._database_status)

    def test_database_property_2(self):
        """Verify _database_status gets reset when database changes."""
        self.handler._database_status = True
        self.handler.database = self.valid_db
        with self.subTest():
            self.assertEqual(self.handler.database, self.valid_db)
        with self.subTest():
            self.assertFalse(self.handler._database_status)

    def test_username_property_1(self):
        """Verify username is set/get correctly."""
        self.handler.username = self.valid_user
        with self.subTest():
            self.assertEqual(self.handler.username, self.valid_user)
        with self.subTest():
            self.assertFalse(self.handler.credential_status)

    def test_username_property_2(self):
        """Verify credential_status is reset when username changes."""
        self.handler.credential_status = True
        self.handler.username = self.valid_user
        with self.subTest():
            self.assertEqual(self.handler.username, self.valid_user)
        with self.subTest():
            self.assertFalse(self.handler.credential_status)

    def test_password_property_1(self):
        """Verify password is set/get correctly."""
        self.handler.password = self.valid_pass
        with self.subTest():
            self.assertEqual(self.handler.password, self.valid_pass)
        with self.subTest():
            self.assertFalse(self.handler.credential_status)

    def test_password_property_2(self):
        """Verify credential_status is reset when password changes."""
        self.handler.credential_status = True
        self.handler.password = self.valid_pass
        with self.subTest():
            self.assertEqual(self.handler.password, self.valid_pass)
        with self.subTest():
            self.assertFalse(self.handler.credential_status)

    def test_login_attempts_property(self):
        """Verify login_attempts is set/get correctly."""
        self.handler.login_attempts = self.attempts
        self.assertEqual(self.handler.login_attempts, self.attempts)

    def test_credential_status_property_1(self):
        """Verify credential_status is set/get correctly."""
        self.handler.credential_status = True
        self.assertTrue(self.handler.credential_status)

    def test_credential_status_property_2(self):
        """Verify credential_status rejects non-boolean values"""
        self.handler.credential_status = "invalid"
        self.assertFalse(self.handler.credential_status)


if __name__ == "__main__":
    unittest.main()
