"""Integration tests for the get_gb pipeline."""

from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

from pdm_utils import run
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.pipelines import get_gb_records

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils


# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_get_gb_records")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder that will get created/removed for each test.
test_folder = Path(test_root_dir, "output")

# How the output folder is named in the get_gb_records pipeline.
results_folder = Path(get_gb_records.RESULTS_FOLDER)
results_path = Path(test_folder, results_folder)

pipeline = "get_gb_records"
USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB

def create_update(table, field, value, phage_id=None):
    """Creates a MySQL UPDATE statement."""
    statement = f"UPDATE {table} SET {field} = '{value}'"
    if phage_id is not None:
        statement = statement + f" WHERE PhageID = '{phage_id}'"
    return statement

def get_unparsed_args():
    """Returns list of command line arguments to get gb records."""
    unparsed_args = ["run.py", pipeline, DB, "-o", str(test_folder)]
    return unparsed_args

D29_ACC = "AF022214"
L5_ACC = "Z18946"
TRIXIE_ACC = "JN408461"

def count_files(path_to_folder):
    """Count number of files in a folder."""
    count = 0
    for item in path_to_folder.iterdir():
        count += 1
    return count

class TestGetGBRecords(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

    def setUp(self):
        test_folder.mkdir()
        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        # Standardize values in certain fields to define the data
        stmt1 = create_update("phage", "Status", "draft")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "HostGenus", "Mycobacterium")
        test_db_utils.execute(stmt2)
        stmt3 = create_update("phage", "Accession", "")
        test_db_utils.execute(stmt3)
        stmt4 = create_update("gene", "Notes", "repressor")
        test_db_utils.execute(stmt4)
        stmt5 = "UPDATE version SET Version = 1"
        test_db_utils.execute(stmt5)
        self.unparsed_args = get_unparsed_args()

    def tearDown(self):
        shutil.rmtree(test_folder)


    @patch("pdm_utils.pipelines.get_gb_records.AlchemyHandler")
    def test_main_1(self, alchemy_mock):
        """Verify no GenBank record is retrieved."""
        alchemy_mock.return_value = self.alchemist
        run.main(self.unparsed_args)
        count = count_files(results_path)
        with self.subTest():
            self.assertTrue(results_path.exists())
        with self.subTest():
            self.assertEqual(count, 0)

    @patch("pdm_utils.pipelines.get_gb_records.AlchemyHandler")
    def test_main_2(self, alchemy_mock):
        """Verify one GenBank record is retrieved."""
        alchemy_mock.return_value = self.alchemist
        stmt = create_update("phage", "Accession", TRIXIE_ACC, "Trixie")
        test_db_utils.execute(stmt)
        run.main(self.unparsed_args)
        count = count_files(results_path)
        with self.subTest():
            self.assertTrue(results_path.exists())
        with self.subTest():
            self.assertEqual(count, 1)

    @patch("pdm_utils.pipelines.get_gb_records.AlchemyHandler")
    def test_main_3(self, alchemy_mock):
        """Verify no GenBank record is retrieved based on one filter."""
        alchemy_mock.return_value = self.alchemist
        stmt = create_update("phage", "Accession", TRIXIE_ACC, "Trixie")
        test_db_utils.execute(stmt)
        self.unparsed_args.extend(["-f", "phage.Status!=draft"])
        run.main(self.unparsed_args)
        count = count_files(results_path)
        with self.subTest():
            self.assertTrue(results_path.exists())
        with self.subTest():
            self.assertEqual(count, 0)

    @patch("pdm_utils.pipelines.get_gb_records.AlchemyHandler")
    def test_main_4(self, alchemy_mock):
        """Verify one GenBank record is retrieved based on one filter."""
        alchemy_mock.return_value = self.alchemist
        stmt1 = create_update("phage", "Accession", TRIXIE_ACC, "Trixie")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "Status", "final", "Trixie")
        test_db_utils.execute(stmt2)
        self.unparsed_args.extend(["-f", "phage.Status!=draft"])
        run.main(self.unparsed_args)
        count = count_files(results_path)
        with self.subTest():
            self.assertTrue(results_path.exists())
        with self.subTest():
            self.assertEqual(count, 1)

if __name__ == '__main__':
    unittest.main()
