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
db2 = "test_db_2"

SET_ALL_DRAFT = "UPDATE phage SET Status = 'draft'"
SET_ALL_FINAL = "UPDATE phage SET Status = 'final'"
COUNT_PHAGE = "SELECT COUNT(*) as count FROM phage"

def get_unparsed_freeze_args():
    """Returns list of command line arguments to freeze database."""
    unparsed_args = ["run.py", pipeline, db,
                      "-n", db2
                    ]
    return unparsed_args

class TestFreeze(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_filled_test_db()
        self.alchemist = AlchemyHandler(database=db, username=user, password=pwd)
        self.alchemist.build_engine()

    def tearDown(self):
        # Remove test_db and test_db_2
        test_db_utils.remove_db()
        exists = test_db_utils.check_if_exists(db=db2)
        if exists:
            test_db_utils.remove_db(db=db2)



    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_1(self, edc_mock):
        """Verify frozen database is created from database with
        no 'draft' genomes."""
        edc_mock.return_value = self.alchemist
        test_db_utils.execute(SET_ALL_FINAL)
        unparsed_args = get_unparsed_freeze_args()
        run.main(unparsed_args)
        count1 = test_db_utils.get_data(COUNT_PHAGE, db=db)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count1[0]["count"], count2[0]["count"])

    @patch("pdm_utils.pipelines.freeze_db.establish_database_connection")
    def test_main_2(self, edc_mock):
        """Verify frozen database is created from database with
        all 'draft' genomes."""
        edc_mock.return_value = self.alchemist
        test_db_utils.execute(SET_ALL_DRAFT)
        unparsed_args = get_unparsed_freeze_args()
        run.main(unparsed_args)
        count1 = test_db_utils.get_data(COUNT_PHAGE, db=db)
        count2 = test_db_utils.get_data(COUNT_PHAGE, db=db2)
        self.assertEqual(count2[0]["count"], 0)




if __name__ == '__main__':
    unittest.main()
