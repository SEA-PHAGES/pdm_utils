"""
Unit tests for the MySQLConnectionHandler class.
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


if __name__ == "__main__":
	unittest.main()
