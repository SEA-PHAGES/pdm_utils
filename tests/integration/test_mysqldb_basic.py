"""Integration tests for schema-agnostic basic MySQL tasks."""

from pathlib import Path
import sys
import shutil
import unittest
from unittest.mock import patch

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.constants import constants
from pdm_utils.functions import mysqldb_basic



# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils
import test_data_utils

USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
DB2 = "pdm_test_2" # This db is sometimes created.
DB3 = "pdm_test_3" # This db is never created.

TABLE = "phage"
COLUMN = "PhageID"

COUNT_QUERY = "SELECT COUNT(*) FROM phage"
PHAGE_QUERY = "SELECT * FROM phage"
GENE_QUERY = "SELECT * FROM gene"

PHAGE = "phage"
GENE = "gene"

# TODO this should be moved to test_db_utils
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
        """Verify non-existing database is not dropped but is created."""
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



    def test_db_exists_1(self):
        """Verify db_exists() can detect existing databases.
        """
        self.assertTrue(mysqldb_basic.db_exists(self.engine, DB))

    def test_db_exists_2(self):
        """Verify db_exists() can detect non-existing databases.
        """
        self.assertFalse(mysqldb_basic.db_exists(self.engine, DB3))

    def test_db_exists_3(self):
        """Verify db_exists() can detect if a database was dropped.
        """
        self.assertTrue(mysqldb_basic.db_exists(self.engine, DB))

        mysqldb_basic.drop_db(self.engine, DB)
        self.assertFalse(mysqldb_basic.db_exists(self.engine, DB))

    def test_get_mysql_dbs_1(self):
        """Verify set of databases is retrieved when engine
        is not connected to a specific database."""
        databases = mysqldb_basic.get_mysql_dbs(self.engine)
        self.assertTrue(DB in databases)




    def test_get_tables_1(self):
        """Verify set of tables is retrieved when engine
        is not connected to a specific database."""
        tables = mysqldb_basic.get_tables(self.engine, DB)
        self.assertTrue(TABLE in tables)




    def test_get_columns_1(self):
        """Verify set of columns is retrieved when engine
        is not connected to a specific database."""
        columns = mysqldb_basic.get_columns(self.engine, DB, TABLE)
        self.assertTrue(COLUMN in columns)



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




    def test_get_mysql_dbs_2(self):
        """Verify set of databases is retrieved when engine
        is connected to a specific database."""
        databases = mysqldb_basic.get_mysql_dbs(self.engine)
        self.assertTrue(DB in databases)




    def test_get_tables_2(self):
        """Verify set of tables is retrieved when engine
        is connected to the same database."""
        tables = mysqldb_basic.get_tables(self.engine, DB)
        self.assertTrue(TABLE in tables)

    def test_get_tables_3(self):
        """Verify set of tables is retrieved when engine
        is connected to a different database."""
        tables = mysqldb_basic.get_tables(self.engine, DB2)
        self.assertTrue(TABLE in tables)




    def test_get_columns_2(self):
        """Verify set of columns is retrieved when engine
        is not connected to the same database."""
        columns = mysqldb_basic.get_columns(self.engine, DB, TABLE)
        self.assertTrue(COLUMN in columns)

    def test_get_columns_3(self):
        """Verify set of columns is retrieved when engine
        is not connected to a different database."""
        columns = mysqldb_basic.get_columns(self.engine, DB2, TABLE)
        self.assertTrue(COLUMN in columns)




class TestMysqldbBasic3(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db()
        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        self.engine = self.alchemist.engine

    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    def test_get_table_count_1(self):
        """Verify the correct number of phages is returned when
        the database is empty."""
        count = mysqldb_basic.get_table_count(self.engine, TABLE)
        self.assertEqual(count, 0)

    def test_get_table_count_2(self):
        """Verify the correct number of phages is returned when
        the database contains one genome."""
        phage_data = test_data_utils.get_trixie_phage_data()
        test_db_utils.insert_data(PHAGE, phage_data)
        count = mysqldb_basic.get_table_count(self.engine, TABLE)
        self.assertEqual(count, 1)




    def test_get_first_row_data_1(self):
        """Verify empty dictionary is returned when there is no data."""
        data = mysqldb_basic.get_first_row_data(self.engine, TABLE)
        self.assertEqual(len(data.keys()), 0)

    def test_get_first_row_data_2(self):
        """Verify dictionary is returned when there is one row of data."""
        phage_data = test_data_utils.get_trixie_phage_data()
        test_db_utils.insert_data(PHAGE, phage_data)
        data = mysqldb_basic.get_first_row_data(self.engine, TABLE)
        self.assertTrue(COLUMN in data.keys())

    def test_get_first_row_data_3(self):
        """Verify dictionary is returned when there are two rows of data."""
        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data1["PhageID"] = "Trixie"
        phage_data2["PhageID"] = "L5"
        test_db_utils.insert_data(PHAGE, phage_data1)
        test_db_utils.insert_data(PHAGE, phage_data2)
        # Get all data from table just to confirm there is more than one row.
        all_data = test_db_utils.get_data(test_db_utils.phage_table_query)
        data = mysqldb_basic.get_first_row_data(self.engine, TABLE)
        with self.subTest():
            self.assertEqual(len(all_data), 2)
        with self.subTest():
            self.assertTrue(COLUMN in data.keys())




    def test_first_1(self):
        """Verify dictionary is returned when there are two rows of data."""
        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data1["PhageID"] = "Trixie"
        phage_data2["PhageID"] = "L5"
        test_db_utils.insert_data(PHAGE, phage_data1)
        test_db_utils.insert_data(PHAGE, phage_data2)
        data = mysqldb_basic.first(self.engine, PHAGE_QUERY, return_dict=True)
        self.assertTrue(COLUMN in data.keys())

    def test_first_2(self):
        """Verify tuple is returned when there are two rows of data."""
        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data1["PhageID"] = "Trixie"
        phage_data2["PhageID"] = "L5"
        test_db_utils.insert_data(PHAGE, phage_data1)
        test_db_utils.insert_data(PHAGE, phage_data2)
        data = mysqldb_basic.first(self.engine, PHAGE_QUERY, return_dict=False)
        with self.subTest():
            self.assertIsInstance(data, tuple)
        with self.subTest():
            self.assertTrue(len(data) > 1)


    def test_scalar_1(self):
        """Verify dictionary is returned when there are two rows of data."""
        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data1["PhageID"] = "Trixie"
        phage_data2["PhageID"] = "L5"
        test_db_utils.insert_data(PHAGE, phage_data1)
        test_db_utils.insert_data(PHAGE, phage_data2)
        count = mysqldb_basic.scalar(self.engine, COUNT_QUERY)
        self.assertEqual(count, 2)


class TestMysqldbBasic4(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db(db=DB2)
        self.alchemist2 = AlchemyHandler(database=DB2, username=USER, password=PWD)
        self.alchemist2.build_engine()
        self.engine2 = self.alchemist2.engine

        test_db_utils.create_empty_test_db()
        self.alchemist1 = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist1.build_engine()
        self.engine1 = self.alchemist1.engine

        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data3 = test_data_utils.get_trixie_phage_data()

        phage_data1["PhageID"] = "L5"
        phage_data2["PhageID"] = "Trixie"
        phage_data3["PhageID"] = "D29"

        phage_data1["HostGenus"] = "Mycobacterium"
        phage_data2["HostGenus"] = "Mycobacterium"
        phage_data3["HostGenus"] = "Gordonia"

        phage_data1["Accession"] = "ABC123"
        phage_data2["Accession"] = "XYZ456"
        phage_data3["Accession"] = ""

        phage_data1["Cluster"] = "A"
        phage_data2["Cluster"] = "B"
        phage_data3["Cluster"] = "NULL"

        phage_data1["Subcluster"] = "A1"
        phage_data2["Subcluster"] = "NULL"
        phage_data3["Subcluster"] = "NULL"

        phage_data1["Sequence"] = "atcg"
        phage_data2["Sequence"] = "AATT"
        phage_data3["Sequence"] = "GGCC"

        phage_data1["Length"] = 6
        phage_data2["Length"] = 4
        phage_data3["Length"] = 5

        phage_data1["DateLastModified"] = constants.EMPTY_DATE
        phage_data2["DateLastModified"] = constants.EMPTY_DATE
        phage_data3["DateLastModified"] = constants.EMPTY_DATE

        phage_data_list = [phage_data1, phage_data2, phage_data3]
        for phage_data in phage_data_list:
            test_db_utils.insert_data(PHAGE, phage_data)

        gene_data1 = test_data_utils.get_trixie_gene_data()
        gene_data2 = test_data_utils.get_trixie_gene_data()
        gene_data3 = test_data_utils.get_trixie_gene_data()
        gene_data4 = test_data_utils.get_trixie_gene_data()

        gene_data1["PhageID"] = "Trixie"
        gene_data2["PhageID"] = "Trixie"
        gene_data3["PhageID"] = "Trixie"
        gene_data4["PhageID"] = "D29"

        gene_data1["GeneID"] = "Trixie_1"
        gene_data2["GeneID"] = "Trixie_2"
        gene_data3["GeneID"] = "Trixie_3"
        gene_data4["GeneID"] = "D29_1"

        gene_data_list = [gene_data1, gene_data2, gene_data3, gene_data4]
        for gene_data in gene_data_list:
            test_db_utils.insert_data(GENE, gene_data)


    def tearDown(self):
        self.engine1.dispose()
        self.engine2.dispose()
        test_db_utils.remove_db()
        test_db_utils.remove_db(db=DB2)




    def test_get_distinct_1(self):
        """Retrieve a set of all distinct values when data is not present."""
        result = mysqldb_basic.get_distinct(self.engine2, "phage", "PhageID")
        exp = set()
        self.assertEqual(result, exp)

    def test_get_distinct_2(self):
        """Retrieve a set of all distinct values when data is present."""
        result1 = mysqldb_basic.get_distinct(
                        self.engine1, "phage", "PhageID")
        result2 = mysqldb_basic.get_distinct(
                        self.engine1, "phage", "HostGenus", null="test")
        result3 = mysqldb_basic.get_distinct(
                        self.engine1, "phage", "Accession")
        result4 = mysqldb_basic.get_distinct(
                        self.engine1, "phage", "Cluster", null="Singleton")
        result5 = mysqldb_basic.get_distinct(
                        self.engine1, "phage", "Subcluster", null="none")

        exp_phage_id = {"L5", "Trixie", "D29"}
        exp_host_genus = {"Mycobacterium", "Gordonia"}
        exp_accession = {"ABC123", "XYZ456", ""}
        exp_cluster = {"A", "B", "Singleton"}
        exp_subcluster = {"A1", "none"}

        with self.subTest():
            self.assertEqual(result1, exp_phage_id)
        with self.subTest():
            self.assertEqual(result2, exp_host_genus)
        with self.subTest():
            self.assertEqual(result3, exp_accession)
        with self.subTest():
            self.assertEqual(result4, exp_cluster)
        with self.subTest():
            self.assertEqual(result5, exp_subcluster)




    def test_retrieve_data_1(self):
        """Verify that a dictionary of data is retrieved for a valid PhageID."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["L5"], query=PHAGE_QUERY)
        with self.subTest():
            self.assertEqual(len(result_list[0].keys()), 14)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "L5")

    def test_retrieve_data_2(self):
        """Verify that an empty list is retrieved for an invalid PhageID."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["EagleEye"], query=PHAGE_QUERY)
        self.assertEqual(len(result_list), 0)

    def test_retrieve_data_3(self):
        """Verify that dictionaries of data are retrieved for a list of two
        valid PhageIDs."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["L5","Trixie"], query=PHAGE_QUERY)
        self.assertEqual(len(result_list), 2)

    def test_retrieve_data_4(self):
        """Verify that dictionaries of data are retrieved for a list of three
        valid PhageIDs and one invalid PhageID."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["L5","Trixie","EagleEye","D29"],
                        query=PHAGE_QUERY)
        self.assertEqual(len(result_list), 3)

    def test_retrieve_data_5(self):
        """Verify that dictionaries of data are retrieved for multiple
        valid PhageIDs when no list is provided."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, query=PHAGE_QUERY)
        self.assertEqual(len(result_list), 3)

    def test_retrieve_data_6(self):
        """Verify that a list of CDS data is retrieved for a valid PhageID."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["Trixie"], query=GENE_QUERY)
        with self.subTest():
            self.assertEqual(len(result_list), 3)
        with self.subTest():
            self.assertEqual(len(result_list[0].keys()), 13)
        with self.subTest():
            self.assertEqual(result_list[0]["PhageID"], "Trixie")

    def test_retrieve_data_7(self):
        """Verify that an empty list of CDS data is retrieved
        for an invalid PhageID."""
        result_list = mysqldb_basic.retrieve_data(
                        self.engine1, column="PhageID",
                        id_list=["L5"], query=GENE_QUERY)
        self.assertEqual(len(result_list), 0)

    def test_retrieve_data_8(self):
        """Verify that a list of all CDS data is retrieved when no
        PhageID is provided."""
        result_list = mysqldb_basic.retrieve_data(
                                self.engine1, query=GENE_QUERY)
        self.assertEqual(len(result_list), 4)




if __name__ == '__main__':
    unittest.main()
