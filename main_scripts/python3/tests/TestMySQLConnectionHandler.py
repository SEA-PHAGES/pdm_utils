import unittest
from main_scripts.python3.classes.MySQLConnectionHandler import \
	MySQLConnectionHandler


class TestMySQLConnectionHandler(unittest.TestCase):
	def setUp(self):
		self.handler = MySQLConnectionHandler()

	def test_set_database(self):
		self.handler.set_database("Actino_Draft")
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertIsNone(self.handler.username)
			self.assertIsNone(self.handler.password)
			self.assertFalse(self.handler.have_connection)
			self.assertFalse(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_set_username(self):
		self.handler.set_username("anonymous")
		with self.subTest():
			self.assertIsNone(self.handler.database)
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertIsNone(self.handler.password)
			self.assertFalse(self.handler.have_connection)
			self.assertFalse(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_set_password(self):
		self.handler.set_password("anonymous")
		with self.subTest():
			self.assertIsNone(self.handler.database)
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertIsNone(self.handler.username)
			self.assertEqual(self.handler.password, "anonymous")
			self.assertFalse(self.handler.have_connection)
			self.assertFalse(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_test_username_password_good(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.test_username_and_password()
		with self.subTest():
			self.assertIsNone(self.handler.database)
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertFalse(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_test_username_password_bad(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymus")
		self.handler.test_username_and_password()
		with self.subTest():
			self.assertIsNone(self.handler.database)
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymus")
			self.assertFalse(self.handler.have_connection)
			self.assertFalse(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_test_database_checked_creds_good_db(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Draft")
		self.handler.test_username_and_password()
		self.handler.test_database()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertFalse(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertTrue(self.handler.valid_database)

	def test_test_database_unchecked_creds_good_db(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Draft")
		self.handler.test_database()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertFalse(self.handler.have_connection)
			self.assertFalse(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_test_database_checked_creds_bad_db(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Daft")
		self.handler.test_username_and_password()
		self.handler.test_database()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Daft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertFalse(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertFalse(self.handler.valid_database)

	def test_create_connection_unchecked_good_creds_good_db(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Draft")
		self.handler.create_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertTrue(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertTrue(self.handler.valid_database)

	def test_create_connection_unchecked_bad_creds_good_db(self):
		# Should be prompted to input valid credentials
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymus")
		self.handler.set_database("Actino_Draft")
		self.handler.create_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			if self.handler.attempts_remaining > 0:
				self.assertTrue(self.handler.have_connection)

	def test_create_connection_unchecked_good_creds_bad_db(self):
		# Should notify bad database for user
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Daft")
		self.handler.create_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Daft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")

	def test_create_connection_checked_good_creds_good_db(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Draft")
		self.handler.test_username_and_password()
		self.handler.create_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertTrue(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertTrue(self.handler.valid_database)

	def test_close_connection(self):
		self.handler.set_username("anonymous")
		self.handler.set_password("anonymous")
		self.handler.set_database("Actino_Draft")
		self.handler.test_username_and_password()
		self.handler.create_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertTrue(self.handler.have_connection)
			self.assertTrue(self.handler.valid_credentials)
			self.assertTrue(self.handler.valid_database)
		self.handler.close_connection()
		with self.subTest():
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")
			self.assertEqual(self.handler.database, "Actino_Draft")
			self.assertEqual(self.handler.attempts_remaining, 3)
			self.assertEqual(self.handler.username, "anonymous")
			self.assertEqual(self.handler.password, "anonymous")


if __name__ == '__main__':
	unittest.main()