"""
Integration tests for the MySQLConnectionHandler class.
Assumptions: 'anonymous'@'localhost' identified by 'anonymous' exists and
has select granted on all tables of all databases.
"""

import unittest
from unittest.mock import patch
from main_scripts.python3.classes.MySQLConnectionHandler import \
	MySQLConnectionHandler


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

	@patch("classes.MySQLConnectionHandler"
		   ".get_credentials", return_value="anonymous")
	def test_get_credentials_1(self):
		"""If credential_status isn't True, even if valid username and password
		are initialized, should ask for user/pass up to 3 times."""
		self.handler.username = self.valid_user
		self.handler.password = self.valid_pass
		self.handler.get_credentials()
		self.assertTrue(self.handler.credential_status)


if __name__ == "__main__":
	unittest.main()
