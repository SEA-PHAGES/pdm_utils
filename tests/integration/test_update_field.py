"""Integration tests for the update_field pipeline."""

import csv
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

from pdm_utils import run
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.pipelines import update_field

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

pipeline = "update"
DB = test_db_utils.DB
USER = test_db_utils.USER
PWD = test_db_utils.PWD


# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_update_field")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder and update table that will get created/removed for each test.
test_folder = Path(test_root_dir, "input")
update_table = Path(test_folder, "update_table.csv")


def create_update(table, field, value, phage_id=None):
    """Creates a MySQL UPDATE statement."""
    num_field_set = {"RetrieveRecord", "AnnotationAuthor"}
    if field not in num_field_set:
        value = f"'{value}'"
    statement = f"UPDATE {table} SET {field} = {value}"
    if phage_id is not None:
        statement = statement + f" WHERE PhageID = '{phage_id}'"
    return statement

def get_alice_ticket(field, value):
    """Returns a dictionary of update ticket data for Alice."""
    dict = {
        "table": "phage",
        "field": field,
        "value": value,
        "key_name": "PhageID",
        "key_value": "Alice"
        }
    return dict

def create_update_table(list_of_data_dicts, file_path):
    """Create an update ticket table."""

    headers = set()
    for dict in list_of_data_dicts:
        headers = headers | dict.keys()

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    with open(file_path, "w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        file_writer.writerow(headers_dict)
        for data_dict in list_of_data_dicts:
            file_writer.writerow(data_dict)

def get_unparsed_args(file=None, version=False):
    """Returns list of command line arguments to update fields."""
    unparsed_args = ["run.py", pipeline, DB]
    if file is not None:
        unparsed_args.extend(["-f", str(file)])
    if version:
        unparsed_args.append("-v")
    return unparsed_args

def phage_id_dict(list_of_dicts):
    """Convert list of MySQL data to dictionary where key = PhageID."""
    dict1 = {}
    for x in range(len(list_of_dicts)):
        dict2 = list_of_dicts[x]
        phage_id = dict2["PhageID"]
        dict1[phage_id] = dict2
    return dict1




class TestUpdate(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

    def setUp(self):

        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        test_folder.mkdir()

        # Standardize values in certain fields to define the data
        stmt1 = create_update("phage", "Status", "unknown")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "HostGenus", "unknown")
        test_db_utils.execute(stmt2)
        stmt3 = create_update("phage", "Accession", "")
        test_db_utils.execute(stmt3)
        stmt4 = create_update("phage", "Cluster", "Z")
        test_db_utils.execute(stmt4)
        stmt5 = create_update("phage", "Subcluster", "Z1")
        test_db_utils.execute(stmt5)
        stmt6 = "UPDATE version SET Version = 0"
        test_db_utils.execute(stmt6)

    def tearDown(self):
        shutil.rmtree(test_folder)




    @patch("pdm_utils.pipelines.update_field.establish_database_connection")
    def test_main_1(self, edc_mock):
        """Verify update runs with empty ticket table."""
        edc_mock.return_value = self.alchemist
        create_update_table([], update_table)
        unparsed_args = get_unparsed_args(file=update_table)
        run.main(unparsed_args)
        version_table = test_db_utils.get_data(test_db_utils.version_table_query)
        phage_table = test_db_utils.get_data(test_db_utils.phage_table_query)
        data_dict = phage_id_dict(phage_table)
        alice = data_dict["Alice"]
        trixie = data_dict["Trixie"]
        # Nothing should be different.
        with self.subTest():
            self.assertEqual(alice["HostGenus"], "unknown")
        with self.subTest():
            self.assertEqual(trixie["HostGenus"], "unknown")
        with self.subTest():
            self.assertEqual(version_table[0]["Version"], 0)

    @patch("pdm_utils.pipelines.update_field.establish_database_connection")
    def test_main_2(self, edc_mock):
        """Verify update runs with five tickets in ticket table."""
        edc_mock.return_value = self.alchemist
        host_genus = "Mycobacterium"
        cluster = "A"
        subcluster = "A2"
        status = "final"
        accession = "ABC123"
        tkt1 = get_alice_ticket("HostGenus", host_genus)
        tkt2 = get_alice_ticket("Cluster", cluster)
        tkt3 = get_alice_ticket("Subcluster", subcluster)
        tkt4 = get_alice_ticket("Status", status)
        tkt5 = get_alice_ticket("Accession", accession)
        tkts = [tkt1, tkt2, tkt3, tkt4, tkt5]
        create_update_table(tkts, update_table)
        unparsed_args = get_unparsed_args(file=update_table)
        run.main(unparsed_args)
        version_table = test_db_utils.get_data(test_db_utils.version_table_query)
        phage_table = test_db_utils.get_data(test_db_utils.phage_table_query)
        data_dict = phage_id_dict(phage_table)
        alice = data_dict["Alice"]
        trixie = data_dict["Trixie"]
        with self.subTest():
            self.assertEqual(alice["HostGenus"], host_genus)
        with self.subTest():
            self.assertEqual(alice["Cluster"], cluster)
        with self.subTest():
            self.assertEqual(alice["Subcluster"], subcluster)
        with self.subTest():
            self.assertEqual(alice["Accession"], accession)
        with self.subTest():
            self.assertEqual(alice["Status"], status)
        # Just confirm that only Alice data was changed.
        with self.subTest():
            self.assertEqual(trixie["HostGenus"], "unknown")
        with self.subTest():
            self.assertEqual(version_table[0]["Version"], 0)

    @patch("pdm_utils.pipelines.update_field.establish_database_connection")
    def test_main_3(self, edc_mock):
        """Verify version data is updated."""
        edc_mock.return_value = self.alchemist
        unparsed_args = get_unparsed_args(version=True)
        run.main(unparsed_args)
        version_table = test_db_utils.get_data(test_db_utils.version_table_query)
        phage_table = test_db_utils.get_data(test_db_utils.phage_table_query)
        data_dict = phage_id_dict(phage_table)
        alice = data_dict["Alice"]
        with self.subTest():
            self.assertEqual(version_table[0]["Version"], 1)
        # Just confirm that only version data was changed.
        with self.subTest():
            self.assertEqual(alice["HostGenus"], "unknown")

    @patch("pdm_utils.pipelines.update_field.establish_database_connection")
    def test_main_4(self, edc_mock):
        """Verify version data and phage table data are updated."""
        edc_mock.return_value = self.alchemist
        host_genus = "Mycobacterium"
        tkt = get_alice_ticket("HostGenus", host_genus)
        create_update_table([tkt], update_table)
        unparsed_args = get_unparsed_args(file=update_table, version=True)
        run.main(unparsed_args)
        version_table = test_db_utils.get_data(test_db_utils.version_table_query)
        phage_table = test_db_utils.get_data(test_db_utils.phage_table_query)
        data_dict = phage_id_dict(phage_table)
        alice = data_dict["Alice"]
        with self.subTest():
            self.assertEqual(alice["HostGenus"], host_genus)
        with self.subTest():
            self.assertEqual(version_table[0]["Version"], 1)


if __name__ == '__main__':
    unittest.main()
