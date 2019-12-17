import pymysql as pms
from pdm_utils.functions import mysqldb
from pdm_utils.functions import basic

# TODO unittest entire class.
class RandomFieldUpdateHandler:
    def __init__(self, connection):
        """
        This object is a ticket validation/execution machine for
        performing any number of simple, single-field database updates
        on a MySQL database.
        :param connection: the connection to the MySQL db
        """
        self.connection = connection        # MySQL connection
        self.table = ""                        # database table
        self.table_valid = False            # assume table invalid until proven
        self.field = ""                        # field whose value will be changed
        self.field_valid = False            # assume field invalid until proven
        self.value = ""                        # replacement value
        self.key_name = ""                    # e.g. PhageID, GeneID
        self.key_name_valid = False            # assume key invalid until proven
        self.key_value = ""                    # e.g. Phrann, SEA_PHRANN_29
        self.key_value_valid = False        # assume invalid until proven
        self.valid_ticket = False            # can't execute ticket unless True

    def validate_table(self):
        """
        This function attempts to validate the table by simply querying
        for the table's description.
        :return:
        """
        try:
            cur = self.connection.cursor()
            cur.execute("DESCRIBE {}".format(self.table))
            cur.close()
            self.table_valid = True
        except pms.err.Error:
            self.table_valid = False
            print("\nInvalid table '{}'".format(self.table))
        return

    def validate_field(self):
        """
        This function attempts to validate the replacement field by
        checking whether it's on the list of fields in the indicated
        table.
        :return:
        """
        # Can't run this function without a valid table.
        if self.table_valid is False:
            return
        try:
            valid_fields = []
            cur = self.connection.cursor()
            cur.execute("DESCRIBE {}".format(self.table))
            fields = cur.fetchall()
            for field in fields:
                valid_fields.append(field[0])
            if self.field in valid_fields:
                self.field_valid = True
            else:
                self.field_valid = False
                print("\nInvalid replacement field '{}' for table "
                      "'{}'".format(self.field, self.table))
        except pms.err.Error as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            self.table_valid = False
        return

    def validate_key_name(self):
        """
        This function attempts to validate the selection key by
        checking whether it's a field in the table marked as any kind
        of key.
        :return:
        """
        # Can't run this function without a valid table.
        if self.table_valid is False:
            return
        try:
            valid_keys = []
            cur = self.connection.cursor()
            cur.execute("DESCRIBE {}".format(self.table))
            fields = cur.fetchall()
            for field in fields:
                if field[2] != "":
                    valid_keys.append(field[0])
            if self.key_name in valid_keys:
                self.key_name_valid = True
            else:
                self.key_name_valid = False
                print("\nInvalid selection key '{}' for table '{}'".format(
                    self.key_name, self.table))
        except pms.err.Error as err:
            print("Error {}: {}".format(err.args[0], err.args[1]))
            self.key_name_valid = False
        return

    def validate_key_value(self):
        """
        This function attempts to validate the selection key's value by
        querying the database for the data associated with that key and
        value on the indicated table
        :return:
        """
        # Can't run this function without a valid table or selection key
        if self.table_valid is False:
            return
        elif self.key_name_valid is False:
            return
        try:
            cur = self.connection.cursor()
            cur.execute("SELECT * FROM {} WHERE {} = '{}'".format(
                self.table, self.key_name, self.key_value))
            tuples = cur.fetchall()
            if len(tuples) > 0:
                self.key_value_valid = True
            else:
                self.key_value_valid = False
                print("\nInvalid selection value '{}' for key '{}' in table "
                      "'{}'".format(self.key_value, self.key_name, self.table))
        except pms.err.Error as err:
            print ("Error {}: {}".format(err.args[0], err.args[1]))
            self.key_value_valid = False
        return

    def validate_ticket(self):
        """
        This function runs all 4 of the object's built-in ticket
        validation methods, and checks whether any of the ticket inputs
        were invalid.  If any are invalid, reject the ticket.  If none
        are invalid, accept the ticket.
        :return:
        """
        self.validate_table()
        self.validate_field()
        self.validate_key_name()
        self.validate_key_value()
        if False in [self.table_valid, self.field_valid, self.key_name_valid,
                self.key_value_valid]:
            self.valid_ticket = False
        else:
            self.valid_ticket = True
        return

    def execute_ticket(self):
        """
        This function checks whether the ticket is valid.  If it is not
        valid, the function returns with code 0, indicating failure to
        execute the ticket.  If the ticket is valid, request input from
        the user to verify that they actually want to proceed with the
        update they've proposed.  If response is in the affirmative,
        the ticket is executed.  Otherwise, indicate that this ticket
        will be skipped, and return 0 as the ticket was not executed.
        If an error is encountered during execution of the ticket,
        print error message and return 0.  If the ticket is executed
        without issue, return 1 indicating success.
        :return:
        """
        if self.valid_ticket is False:
            return 0
        try:
            command = mysqldb.create_update(self.table, self.field,
                        self.value, self.key_name, self.key_value)
            # print("\nCommand to execute:")
            # print(command)
            # prompt = "Do you wish to proceed? (y/n) "
            # result = basic.ask_yes_no(prompt=prompt, response_attempt=3)
            result = True
            if result == True:
                cur = self.connection.cursor()
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
