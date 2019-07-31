"""
Integration tests for the MySQLConnectionHandler class.
Assumptions: 'anonymous'@'localhost' identified by 'anonymous' exists and
has select granted on all tables of all databases.
"""

import unittest
from main_scripts.python3.classes.MySQLConnectionHandler import MySQLConnectionHandler


class TestMySQLConnectionHandler(unittest.TestCase):
	def setUp(self):
		self.handler = MySQLConnectionHandler()
		self.valid_user = "anonymous"
		self.valid_pass = "anonymous"
		self.valid_db = "Actino_Draft"
		self.invalid_user = "invalid"
		self.invalid_pass = "invalid"
		self.invalid_db = "invalid"

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


if __name__ == "__main__":
	unittest.main()
