"""Integration tests for the freeze pipeline."""

from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from pdm_utils import run
from pdm_utils.pipelines import freeze_db
from pdm_utils.classes.alchemyhandler import AlchemyHandler

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

pipeline = "freeze"
user = test_db_utils.USER
pwd = test_db_utils.PWD
db = test_db_utils.DB
db2 = "pdm_test_2"
COUNT_PHAGE = "SELECT COUNT(*) as count FROM phage"

def create_update(table, field, value, phage_id=None):
    """Creates a MySQL UPDATE statement."""
    statement = f"UPDATE {table} SET {field} = '{value}'"
    if phage_id is not None:
        statement = statement + f" WHERE PhageID = '{phage_id}'"
    return statement


def get_unparsed_freeze_args():
    """Returns list of command line arguments to freeze database."""
    unparsed_args = ["run.py", pipeline, db,
                      "-n", db2
                    ]
    return unparsed_args

class TestFreeze(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        # Remove 'pdm_test_db'
        test_db_utils.remove_db()

    def setUp(self):
        self.alchemist = AlchemyHandler(database=db, username=user, password=pwd)
        self.alchemist.build_engine()
        # Standardize values in certain fields to define the data
        stmt1 = create_update("phage", "Status", "draft")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "HostGenus", "Mycobacterium")
        test_db_utils.execute(stmt2)
        stmt3 = create_update("gene", "Notes", "repressor")
        test_db_utils.execute(stmt3)
        stmt4 = "UPDATE version SET Version = 1"
        test_db_utils.execute(stmt4)
        self.unparsed_args = get_unparsed_freeze_args()

    def tearDown(self):
        # Remove 'pdm_test_2'
        exists = test_db_utils.check_if_exists(db=db2)
        if exists:
            test_db_utils.remove_db(db=db2)




    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_1(self, edc_mock):
        """Verify frozen database is created from database with
        no change in genome count when no filters are provided."""
        edc_mock.return_value = self.alchemist
        run.main(self.unparsed_args)
        count1 = test_db_utils.get_data(COUNT_PHAGE, db=db)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        version = test_db_utils.get_data(test_db_utils.version_table_query,
                                         db=db2)
        with self.subTest():
            self.assertEqual(count1[0]["count"], count2[0]["count"])
        with self.subTest():
            self.assertEqual(version[0]["Version"], 1)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_2(self, edc_mock):
        """Verify frozen database is created from database with
        all genomes removed."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "phage.Status != draft"])
        run.main(self.unparsed_args)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 0)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_3(self, edc_mock):
        """Verify frozen database is created from database with
        only one 'final' genome."""
        edc_mock.return_value = self.alchemist
        stmt = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt)
        self.unparsed_args.extend(["-f", "phage.Status!=draft"])
        run.main(self.unparsed_args)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 1)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_4(self, edc_mock):
        """Verify frozen database is created from database with
        one genome based on two filters."""
        edc_mock.return_value = self.alchemist
        stmt = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt)
        filters = "phage.Status != draft AND phage.HostGenus = Mycobacterium"
        self.unparsed_args.extend(["-f", filters])
        run.main(self.unparsed_args)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 1)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_5(self, edc_mock):
        """Verify frozen database is created from database with
        no genomes based on two filters."""
        edc_mock.return_value = self.alchemist
        stmt = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt)
        filters = "phage.Status != draft AND phage.HostGenus = Gordonia"
        self.unparsed_args.extend(["-f", filters])
        run.main(self.unparsed_args)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 0)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_6(self, edc_mock):
        """Verify frozen database is created from database with
        one genome based on two filters from two tables."""
        edc_mock.return_value = self.alchemist
        stmt = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt)
        filters = "phage.Status != draft AND gene.Notes = repressor"
        self.unparsed_args.extend(["-f", filters])
        run.main(self.unparsed_args)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 1)

    @patch("sys.exit")
    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_7(self, edc_mock, exit_mock):
        """Verify pipeline exits with invalid filter column."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "phage.Invalid != draft"])
        run.main(self.unparsed_args)
        exit_mock.assert_called()

    @patch("sys.exit")
    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_8(self, edc_mock, exit_mock):
        """Verify pipeline exits with invalid filter table."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "Invalid.Invalid != draft"])
        run.main(self.unparsed_args)
        exit_mock.assert_called()

    @patch("sys.exit")
    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_9(self, edc_mock, exit_mock):
        """Verify pipeline exits with invalid filter string."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "Invalid != draft"])
        run.main(self.unparsed_args)
        exit_mock.assert_called()

    @patch("sys.exit")
    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_10(self, edc_mock, exit_mock):
        """Verify pipeline exits with invalid quoted filter value."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "phage.Status != 'draft'"])
        run.main(self.unparsed_args)
        exit_mock.assert_called()

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_11(self, edc_mock):
        """Verify data is changed when there is data in the database
        and reset = True."""
        edc_mock.return_value = self.alchemist
        stmt = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt)
        self.unparsed_args.extend(["-f", "phage.Status != draft", "-r"])
        run.main(self.unparsed_args)
        count = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        version = test_db_utils.get_data(test_db_utils.version_table_query,
                                         db=db2)
        with self.subTest():
            self.assertEqual(count[0]["count"], 1)
        with self.subTest():
            self.assertEqual(version[0]["Version"], 0)

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_12(self, edc_mock):
        """Verify data is changed when there is NO data in the database
        and reset = True."""
        edc_mock.return_value = self.alchemist
        self.unparsed_args.extend(["-f", "phage.Status != draft", "-r"])
        run.main(self.unparsed_args)
        count = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        version = test_db_utils.get_data(test_db_utils.version_table_query,
                                         db=db2)
        with self.subTest():
            self.assertEqual(count[0]["count"], 0)
        with self.subTest():
            self.assertEqual(version[0]["Version"], 0)

if __name__ == '__main__':
    unittest.main()
