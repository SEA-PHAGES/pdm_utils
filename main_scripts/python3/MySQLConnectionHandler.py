import sys
import pymysql as pms
import getpass


class MySQLConnectionHandler:
	def __init__(self, database):
		"""
		This object is intended to handle creation of a connection to
		a MySQL database.  It gives a user 3 attempts to input valid
		MySQL credentials, and tests the database name to verify it
		exists.  If valid credentials are given and the database is
		real, a connection is established and can be used externally.
		:param database: the name of the MySQL database to connect to
		"""
		self.connection = None
		self.attempts_remaining = 3
		self.username = None
		self.password = None
		self.database = database
		self.successful_login = False
		self.valid_database = False

	def create_connection(self):
		while not self.successful_login:
			if self.attempts_remaining > 0:
				self.get_username_and_password()
				self.test_username_and_password()
				if self.successful_login is False:
					if self.attempts_remaining >= 1:
						print("Incorrect username or password. You have {} "
							  "attempts remaining. Please try "
							  "again.".format(self.attempts_remaining))
					else:
						print("Too many attempts. Please verify that your "
							  "username and password are correct, and try "
							  "again.")
						sys.exit(1)
		self.test_database()
		if self.valid_database is True:
			connection = pms.connect("localhost", self.username,
									 self.password, self.database)
			self.connection = connection
		else:
			print("Invalid database selection. Please verify the name of "
				  "your database, and try again.")
		return

	def get_username_and_password(self):
		self.attempts_remaining -= 1
		self.username = getpass.getpass(prompt="MySQL username: ")
		self.password = getpass.getpass(prompt="MySQL password: ")
		return

	def test_username_and_password(self):
		try:
			con = pms.connect("localhost", self.username, self.password)
			con.close()
			self.successful_login = True
		except pms.err.Error:
			self.successful_login = False
		return

	def test_database(self):
		try:
			con = pms.connect("localhost", self.username, self.password,
							  self.database)
			con.close()
			self.valid_database = True
		except pms.err.Error:
			self.valid_database = False
		return
