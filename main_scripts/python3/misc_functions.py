import getpass
import pymysql as pms
import os


def expand_path(input_path):
	"""
	This function attempts to coerce input paths in any relative format
	into an absolute path. "~" should be expanded to the user's full
	home directory. "./" or "../" notation needs to be expanded.
	Directories should end in "/"
	:param input_path: the path to be expanded
	:return expanded_path: the expanded path to the file/dir indicated
	"""
	home_dir = os.path.expanduser("~")

	if input_path[0] == "~":
		expanded_path = home_dir + input_path[1:]
	else:
		expanded_path = input_path
	expanded_path = os.path.abspath(expanded_path)
	return expanded_path


def verify_path(expanded_path, kind=None):
	if kind == "file":
		if os.path.isfile(expanded_path) is True:
			return True
		else:
			return False
	elif kind == "dir":
		if os.path.isdir(expanded_path) is True:
			return True
		else:
			return False
	elif kind is None:
		if os.path.exists(expanded_path) is True:
			return True
		else:
			return False
	else:
		print("{} is not a valid kind for this function. Please try again "
			  "using one of (None, dir, file).")


def get_mysql_login():
	"""
	This function asks for a mysql username and password from the user.
	:return uname: MySQL username
	:return pword: MySQL password
	"""
	uname = getpass.getpass(prompt="MySQL username: ")
	pword = getpass.getpass(prompt="MySQL password: ")
	return uname, pword


def validate_mysql_login(username, password, database=None):
	"""
	This function tests a username and password (and database if
	supplied) to see whether a connection can be established using
	those parameters.
	:param username: MySQL username
	:param password: MySQL password
	:param database: MySQL database, default is None
	:return Boolean: True if connection is made, False if connection failed
	"""
	if database is not None:
		try:
			con = pms.connect("localhost", username, password, database)
			con.close()
			return True
		except pms.OperationalError:
			return False
	else:
		try:
			con = pms.connect("localhost", username, password)
			con.close()
			return True
		except pms.OperationalError:
			return False


def create_mysql_connection(username, password, database):
	"""
	This function creates a MySQL connection using the verified MySQL username, password, and database arguments. It
	returns that connection if established.
	:param username: the MySQL username
	:param password: the MySQL password
	:param database: the MySQL database
	:return con: pymysql connection object
	"""
	try:
		con = pms.connect("localhost", username, password, database)
		return con
	except pms.OperationalError:
		print("Could not establish connection to the database")
