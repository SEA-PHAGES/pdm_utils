"""
Integration tests for the MySQLConnectionHandler class.
Assumptions: 'anonymous'@'localhost' identified by 'anonymous' exists
with 'anonymous' password and has
select granted on all tables of all databases.
"""

import unittest
from unittest.mock import patch
from classes.MySQLConnectionHandler import MySQLConnectionHandler


class TestMySQLConnectionHandler(unittest.TestCase):
    def setUp(self):
        self.handler = MySQLConnectionHandler()
        self.valid_user = "anonymous"
        self.valid_pass = "anonymous"
        self.valid_db = "Actino_Draft"
        self.invalid_user = "invalid"
        self.invalid_pass = "invalid"
        self.invalid_db = "invalid"
        self.attempts = 5
        self.valid_query = "SELECT database()"
        self.invalid_query = "INVALID"
        self.valid_transaction = ["SELECT COUNT(PhageID) FROM phage",
                                  "SELECT COUNT(GeneID) FROM gene",
                                  "SELECT COUNT(name) FROM pham"]
        self.invalid_transaction = ["SELECT COUNT(PhageID) FROM phage",
                                    "SELECT COUNT(GeneID) FROM gene",
                                    "SELECT COUNT(GeneID) FROM phage"]

    def test_validate_credentials_1(self):
        """Valid username and password should result in True
        credential_status."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.validate_credentials()
        self.assertTrue(self.handler.credential_status)

    def test_validate_credentials_2(self):
        """Valid username with invalid password should result in False
        credential_status."""
        self.handler.username = self.valid_user
        self.handler.password = self.invalid_user
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_credentials_3(self):
        """Invalid username with valid password should result in False
        credential_status."""
        self.handler.username = self.invalid_user
        self.handler.password = self.valid_pass
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_credentials_4(self):
        """Invalid user and invalid password should result in False
        credential_status."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_user
        self.handler.validate_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_validate_database_access_1(self):
        """Valid credentials and valid database should result in True
        database_status if credential_status is True."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        # bypass validate_credentials() because we know True
        self.handler.credential_status = True
        self.handler.validate_database_access()
        self.assertTrue(self.handler._database_status)

    def test_validate_database_access_2(self):
        """Valid credentials and valid database should result in False
        database_status if credential_status is False."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.validate_database_access()
        self.assertFalse(self.handler._database_status)

    def test_get_credentials_1(self):
        """If credential_status isn't True, should ask for
        user/pass...handled by patch and credential_status should be True."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        user_input = [self.valid_user, self.valid_pass]
        # patch getpass so we can feed good credentials without manual input
        with patch('getpass.getpass', side_effect=user_input):
            self.handler.get_credentials()
        self.assertTrue(self.handler.credential_status)

    def test_get_credentials_2(self):
        """Credential status False should ask for user/pass up to 3 times
        ...handled by patch. Give bad values 3x should say too many login
        attempts and credential status stays False."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        user_input = [self.invalid_user, self.invalid_pass,
                      self.invalid_user, self.invalid_pass,
                      self.invalid_user, self.invalid_pass]
        # patch getpass so we can feed good credentials without manual input
        with patch('getpass.getpass', side_effect=user_input):
            self.handler.get_credentials()
        self.assertFalse(self.handler.credential_status)

    def test_get_credentials_3(self):
        """Credential status False should ask for user/pass up to 3 times
        ...handled by patch. Give bad values 2x then correct values should
        be ok and credential_status changes to True."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        user_input = [self.invalid_user, self.invalid_pass,
                      self.invalid_user, self.invalid_pass,
                      self.valid_user, self.valid_pass]
        # patch getpass so we can feed good credentials without manual input
        with patch('getpass.getpass', side_effect=user_input):
            self.handler.get_credentials()
        self.assertTrue(self.handler.credential_status)

    def test_get_credentials_4(self):
        """If credential_status is True, should do nothing even is invalid
        credentials are used."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_pass
        self.handler.credential_status = True
        self.handler.get_credentials()
        self.assertTrue(self.handler.credential_status)

    def test_ask_user_and_pass(self):
        """Credentials should be asked, login_attempts should decrease by 1"""
        self.handler.login_attempts = self.attempts
        user_input = [self.valid_user, self.valid_pass]
        with patch('getpass.getpass', side_efect=user_input):
            self.handler.ask_username_and_password()
        self.assertEqual(self.handler.login_attempts, self.attempts - 1)

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

    def test_open_connection_4(self):
        """If invalid credentials are used with a valid database, should be
        prompted for correct credentials - unvalidated invalid credentials."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_pass
        self.handler.database = self.valid_db
        self.handler._database_status = True
        user_input = [self.valid_user, self.valid_pass]
        with patch('getpass.getpass', side_effect=user_input):
            self.handler.open_connection()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)

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

    def test_execute_query_1(self):
        """With a valid connection and valid query, cursor should be made and
        query should be executed"""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        return_list = self.handler.execute_query(self.valid_query)
        self.assertEqual(return_list[0]["database()"], self.valid_db.lower())

    def test_execute_query_2(self):
        """With a valid connection and invalid query, cursor should be made
        and query should fail to be executed. Error message should print"""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        return_list = self.handler.execute_query(self.invalid_query)
        self.assertIsNone(return_list)

    def test_execute_query_3(self):
        """With a valid query and credentials/db, but no connection, a
        connection should be made, cursor should be made and query should be
        executed."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_list = self.handler.execute_query(self.valid_query)
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        with self.subTest():
            self.assertEqual(return_list[0]["database()"],
                self.valid_db.lower())

    def test_execute_query_4(self):
        """With valid credentials/db but no connection and invalid query,
        a connection should be made, cursor should be made but query should
        fail to be executed. Error message should print."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_list = self.handler.execute_query(self.invalid_query)
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        with self.subTest():
            self.assertIsNone(return_list)

    def test_execute_query_5(self):
        """With invalid credentials, valid db, no connection and valid query, a
        connection should fail to be made - error message should print."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_pass
        self.handler.database = self.valid_db
        # Force connection attempt with bad credentials - otherwise
        # credentials would be asked
        self.handler.credential_status = True
        self.handler._database_status = False
        return_list = self.handler.execute_query(self.valid_query)
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        with self.subTest():
            self.assertIsNone(return_list)

    def test_execute_query_6(self):
        """With valid credentials, invalid db, no connection and a valid
        query, a connection should fail to be made - error message should
        print."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.invalid_db
        self.handler.credential_status = True
        # Force connection attempt with bad database
        self.handler._database_status = False
        return_list = self.handler.execute_query(self.valid_query)
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        with self.subTest():
            self.assertIsNone(return_list)

    def test_execute_query_7(self):
        """Valid everything but no query should be ok; error message will
        print and nothing will return."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_list = self.handler.execute_query()
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        with self.subTest():
            self.assertIsNone(return_list)

    def test_execute_transaction_1(self):
        """Valid everything should result in creation of cursor and execution
        of all statements in the transaction - return code 0."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        return_code = self.handler.execute_transaction(self.valid_transaction)
        self.assertEqual(return_code, 0)

    def test_execute_transaction_2(self):
        """Valid connection but invalid transaction should return code 1."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        self.handler.open_connection()
        return_code = self.handler.execute_transaction(self.invalid_transaction)
        self.assertEqual(return_code, 1)

    def test_execute_transaction_3(self):
        """No connection but valid credentials/db and valid transaction
        should return code 0."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_code = self.handler.execute_transaction(self.valid_transaction)
        self.assertEqual(return_code, 0)

    def test_execute_transaction_4(self):
        """No connection but valid credentials/db and invalid transaction
        should return code 1."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_code = self.handler.execute_transaction(self.invalid_transaction)
        with self.subTest():
            self.assertIsNotNone(self.handler.connection)
        with self.subTest():
            self.assertTrue(self.handler.connection_status())
        with self.subTest():
            self.assertEqual(return_code, 1)

    def test_execute_transaction_5(self):
        """No connection and invalid credentials with valid db and
        transaction should return code 1."""
        self.handler.username = self.invalid_user
        self.handler.password = self.invalid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_code = self.handler.execute_transaction(self.valid_transaction)
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        with self.subTest():
            self.assertEqual(return_code, 1)

    def test_execute_transaction_6(self):
        """No connection and invalid db with valid credentials and
        transaction should return code 1."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.invalid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_code = self.handler.execute_transaction(self.valid_transaction)
        with self.subTest():
            self.assertIsNone(self.handler.connection)
        with self.subTest():
            self.assertFalse(self.handler.connection_status())
        with self.subTest():
            self.assertEqual(return_code, 1)

    def test_execute_transaction_7(self):
        """Everything ok but no transaction should return 0."""
        self.handler.username = self.valid_user
        self.handler.password = self.valid_pass
        self.handler.database = self.valid_db
        self.handler.credential_status = True
        self.handler._database_status = True
        return_code = self.handler.execute_transaction()
        self.assertEqual(return_code, 0)


if __name__ == "__main__":
    unittest.main()
