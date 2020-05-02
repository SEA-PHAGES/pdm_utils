"""Integration tests for basic MySQL tasks."""

from pathlib import Path
import sys
import shutil
import unittest
from unittest.mock import patch

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import mysqldb_basic
from pdm_utils.pipelines import convert_db



# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
DB2 = "pdm_test_2" # This db is sometimes created.
DB3 = "pdm_test_3" # This db is never created.

# Create the main test directory in which all files will be
# created and managed. Gets created once for all tests.
test_root_dir = Path("/tmp", "pdm_utils_tests_mysqldb_basic")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# Within main test folder, create new folder that will get created/removed
# for each test.
results_path = Path(test_root_dir, "output")

TABLES_QUERY = ("SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = '{}'")


class TestMysqldbBasic1(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_filled_test_db()
        self.alchemist = AlchemyHandler(username=USER, password=PWD)
        self.alchemist.build_engine()
        self.engine = self.alchemist.engine

    def tearDown(self):
        self.engine.dispose()
        # Remove 'pdm_test_db'
        exists = test_db_utils.check_if_exists()
        if exists:
            test_db_utils.remove_db()

        # Remove 'pdm_test_2'
        exists = test_db_utils.check_if_exists(db=DB2)
        if exists:
            test_db_utils.remove_db(db=DB2)



    def test_drop_db_1(self):
        """Verify existing database is dropped."""
        before = test_db_utils.check_if_exists()
        result = mysqldb_basic.drop_db(self.engine, DB)
        after = test_db_utils.check_if_exists()
        with self.subTest():
            self.assertTrue(before)
        with self.subTest():
            self.assertFalse(after)
        with self.subTest():
            self.assertEqual(result, 0)

    def test_drop_db_2(self):
        """Verify non-existing database is not dropped."""
        before = test_db_utils.check_if_exists(db=DB2)
        result = mysqldb_basic.drop_db(self.engine, DB2)
        after = test_db_utils.check_if_exists(db=DB2)
        with self.subTest():
            self.assertFalse(before)
        with self.subTest():
            self.assertFalse(after)
        with self.subTest():
            self.assertEqual(result, 1)




    def test_create_db_1(self):
        """Verify new database is created."""
        before = test_db_utils.check_if_exists(db=DB2)
        result = mysqldb_basic.create_db(self.engine, DB2)
        after = test_db_utils.check_if_exists(db=DB2)
        with self.subTest():
            self.assertFalse(before)
        with self.subTest():
            self.assertTrue(after)
        with self.subTest():
            self.assertEqual(result, 0)

    def test_create_db_2(self):
        """Verify already-existing database is not created."""
        before = test_db_utils.check_if_exists()
        result = mysqldb_basic.create_db(self.engine, DB)
        after = test_db_utils.check_if_exists()
        with self.subTest():
            self.assertTrue(before)
        with self.subTest():
            self.assertTrue(after)
        with self.subTest():
            self.assertEqual(result, 1)




    def test_drop_create_db_1(self):
        """Verify already-existing database is dropped and created."""
        before = test_db_utils.check_if_exists()
        before_tables = test_db_utils.execute(TABLES_QUERY.format(DB))
        result = mysqldb_basic.drop_create_db(self.engine, DB)
        after = test_db_utils.check_if_exists()
        after_tables = test_db_utils.execute(TABLES_QUERY.format(DB))
        with self.subTest():
            self.assertTrue(before)
        with self.subTest():
            self.assertTrue(after)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertTrue(len(before_tables) > 0)
        with self.subTest():
            self.assertTrue(len(after_tables) == 0)

    def test_drop_create_db_2(self):
        """Verify non-existing database is created."""
        before = test_db_utils.check_if_exists(DB2)
        result = mysqldb_basic.drop_create_db(self.engine, DB2)
        after = test_db_utils.check_if_exists(DB2)
        after_tables = test_db_utils.execute(TABLES_QUERY.format(DB2), db=DB2)
        with self.subTest():
            self.assertFalse(before)
        with self.subTest():
            self.assertTrue(after)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertTrue(len(after_tables) == 0)


    @patch("pdm_utils.functions.mysqldb_basic.drop_db")
    def test_drop_create_db_3(self, drop_mock):
        """Verify database is not created if there is an error during drop."""
        drop_mock.return_value = 1
        before = test_db_utils.check_if_exists()
        before_tables = test_db_utils.execute(TABLES_QUERY.format(DB))

        result = mysqldb_basic.drop_create_db(self.engine, DB)
        after = test_db_utils.check_if_exists()
        after_tables = test_db_utils.execute(TABLES_QUERY.format(DB))
        with self.subTest():
            self.assertTrue(before)
        with self.subTest():
            self.assertTrue(after)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertTrue(len(before_tables) > 0)
        with self.subTest():
            self.assertTrue(len(after_tables) > 0)
        with self.subTest():
            self.assertEqual(len(before_tables), len(after_tables))
        with self.subTest():
            drop_mock.assert_called()




class TestMysqldbBasic2(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_filled_test_db()
        test_db_utils.create_empty_test_db(db=DB2)

        stmt1 = "UPDATE version SET Version = 1"
        test_db_utils.execute(stmt1)
        stmt2 = "UPDATE version SET Version = 0"
        test_db_utils.execute(stmt2, db=DB2)

        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        self.engine = self.alchemist.engine

    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()
        test_db_utils.remove_db(db=DB2)




    def test_copy_db_1(self):
        """Verify data from database1 is copied to database2."""
        before_v1 = test_db_utils.get_data(test_db_utils.version_table_query)
        before_v2 = test_db_utils.get_data(test_db_utils.version_table_query, db=DB2)
        result = mysqldb_basic.copy_db(self.engine, DB2)
        after_v1 = test_db_utils.get_data(test_db_utils.version_table_query)
        after_v2 = test_db_utils.get_data(test_db_utils.version_table_query, db=DB2)
        with self.subTest():
            self.assertNotEqual(before_v1[0]["Version"], before_v2[0]["Version"])
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertEqual(after_v1[0]["Version"], after_v2[0]["Version"])

    def test_copy_db_2(self):
        """Verify no data is copied since databases are the same."""
        before_v1 = test_db_utils.get_data(test_db_utils.version_table_query)
        result = mysqldb_basic.copy_db(self.engine, DB)
        after_v1 = test_db_utils.get_data(test_db_utils.version_table_query)
        with self.subTest():
            self.assertEqual(before_v1[0]["Version"], after_v1[0]["Version"])
        with self.subTest():
            self.assertEqual(result, 0)

    def test_copy_db_3(self):
        """Verify no data is copied since new database does not exist."""
        result = mysqldb_basic.copy_db(self.engine, DB3)
        self.assertEqual(result, 1)

    @patch("pdm_utils.functions.mysqldb_basic.pipe_commands")
    def test_copy_db_4(self, pc_mock):
        """Verify no data is copied if en error is encountered during copying."""
        # Raise an error instead of calling pipe_commands() so that
        # the exception block is entered.
        pc_mock.side_effect = ValueError("Error raised")
        result = mysqldb_basic.copy_db(self.engine, DB2)
        self.assertEqual(result, 1)




    def test_install_db_1(self):
        """Verify new database is installed."""
        stmt1 = "UPDATE version SET Version = 0"
        test_db_utils.execute(stmt1)
        before = test_db_utils.get_data(test_db_utils.version_table_query)
        result = mysqldb_basic.install_db(self.engine, test_db_utils.TEST_DB_FILEPATH)
        after = test_db_utils.get_data(test_db_utils.version_table_query)
        with self.subTest():
            self.assertNotEqual(before[0]["Version"], after[0]["Version"])
        with self.subTest():
            self.assertTrue(after[0]["Version"] > 0)
        with self.subTest():
            self.assertEqual(result, 0)

    @patch("subprocess.check_call")
    def test_install_db_2(self, cc_mock):
        """Verify new database is not installed."""
        # Raise an error instead of calling check_call() so that
        # the exception block is entered.
        cc_mock.side_effect = ValueError("Error raised")
        stmt1 = "UPDATE version SET Version = 0"
        test_db_utils.execute(stmt1)
        before = test_db_utils.get_data(test_db_utils.version_table_query)
        result = mysqldb_basic.install_db(self.engine, test_db_utils.TEST_DB_FILEPATH)
        after = test_db_utils.get_data(test_db_utils.version_table_query)
        with self.subTest():
            self.assertEqual(before[0]["Version"], after[0]["Version"])
        with self.subTest():
            self.assertTrue(after[0]["Version"] == 0)
        with self.subTest():
            self.assertEqual(result, 1)




class TestMysqldbBasic3(unittest.TestCase):

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

if __name__ == '__main__':
    unittest.main()
