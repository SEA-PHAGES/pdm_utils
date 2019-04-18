import pymysql as pms


class RandomFieldUpdateHandler:
	def __init__(self, connection):
		self.connection = connection
		self.table = ""
		self.field = ""
		self.value = ""
		self.key_name = ""
		self.key_value = ""
		self.valid_ticket = False

	def execute_ticket(self):
		if self.valid_ticket is False:
			return 0
		try:

			print("\nCommand to execute:")
			print("UPDATE {} SET {} = '{}' WHERE {} = '{}'".format(
				self.table, self.field, self.value, self.key_name,
				self.key_value))
			proceed = input("Do you wish to proceed? (y/n) ")
			if proceed.lower() in ["yes", "y"]:
				cur = self.connection.cursor()
				command = "UPDATE {} SET {} = '{}' WHERE {} = '{}'".format(
					self.table, self.field, self.value, self.key_name,
					self.key_value)
				cur.execute(command)
				cur.execute("COMMIT")
				cur.close()
			else:
				print("Skipping this ticket...")
				return 0
		except pms.err.Error as err:
			print("Error {}: {}".format(err[0], err[1]))
			return 0
		return 1

	def validate_ticket(self):
		# Validate the table
		valid_fields = []
		valid_keys = []
		try:
			cur = self.connection.cursor()
			command = "DESCRIBE {}".format(self.table)
			cur.execute(command)
			tuples = cur.fetchall()
			for t in tuples:
				valid_fields.append(t[0])
				if t[2] != "":
					valid_keys.append(t[0])
				else:
					pass
			command = "SELECT * FROM {} WHERE {} = '{}'".format(self.table,
																self.key_name,
																self.key_value)
			cur.execute(command)
			tuples = cur.fetchall()
			if len(tuples) > 0:
				if self.key_name in valid_keys:
					if self.field in valid_fields:
						self.valid_ticket = True
					else:
						print("\nInvalid field '{}'".format(self.field))
						self.valid_ticket = False
				else:
					print("\nInvalid key '{}'".format(self.key_name))
					self.valid_ticket = False
			else:
				print("\nInvalid key value '{}'".format(self.key_value))
				self.valid_ticket = False
		except pms.err.Error:
			print("\nInvalid table '{}'".format(self.table))
			self.valid_ticket = False
