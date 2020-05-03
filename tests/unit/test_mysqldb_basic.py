"""Unit tests for schema-agnostic basic MySQL tasks."""

from pathlib import Path
import sys
import unittest

from pdm_utils.functions import mysqldb_basic

USER = "user"
PWD = "pwd"
DB = "db"

class TestMysqldbBasic1(unittest.TestCase):

    def test_mysqldump_command_1(self):
        """Verify mysqldump command list is correctly constructed."""
        result = mysqldb_basic.mysqldump_command(USER, PWD, DB)
        exp = ["mysqldump", "-u", USER, f"-p{PWD}", DB]
        with self.subTest():
            self.assertEqual(len(result), len(exp))

        for i in range(len(exp)):
            with self.subTest():
                self.assertEqual(result[i], exp[i])




    def test_mysql_login_command_1(self):
        """Verify mysql command list is correctly constructed."""
        result = mysqldb_basic.mysql_login_command(USER, PWD, DB)
        exp = ["mysql", "-u", USER, f"-p{PWD}", DB]
        with self.subTest():
            self.assertEqual(len(result), len(exp))

        for i in range(len(exp)):
            with self.subTest():
                self.assertEqual(result[i], exp[i])




    def test_convert_for_sql_1(self):
        """Verify non-empty value returned is encapsulated with "'"."""
        value = mysqldb_basic.convert_for_sql(
                                    "A", check_set={"Singleton"}, single=True)
        self.assertEqual(value, "'A'")

    def test_convert_for_sql_2(self):
        """Verify non-empty value returned is encapsulated with '"'."""
        value = mysqldb_basic.convert_for_sql(
                                    "A", check_set={"Singleton"}, single=False)
        self.assertEqual(value, '"A"')

    def test_convert_for_sql_3(self):
        """Verify empty value returned is NULL."""
        value = mysqldb_basic.convert_for_sql("", check_set={""}, single=True)
        self.assertEqual(value, "NULL")

    def test_convert_for_sql_4(self):
        """Verify 'Singleton' value returned is NULL."""
        value = mysqldb_basic.convert_for_sql(
                            "Singleton", check_set={"Singleton"}, single=True)
        self.assertEqual(value, "NULL")



if __name__ == '__main__':
    unittest.main()
