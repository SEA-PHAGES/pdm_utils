from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.pipelines import review

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils


USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
TEST_DIR = "/tmp/pdm_utils_tests_review"

class TestPhamReview(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

        self.test_dir = Path(TEST_DIR)
        if self.test_dir.is_dir():
            shutil.rmtree(TEST_DIR)

        self.test_dir.mkdir()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()
        shutil.rmtree(TEST_DIR)

    def setUp(self):
        self.alchemist = AlchemyHandler()
        self.alchemist.username=USER
        self.alchemist.password=PWD
        self.alchemist.database=DB
        self.alchemist.connect(ask_database=True, login_attempts=0)

        self.review_test_dir = self.test_dir.joinpath("review_test_dir")

    def tearDown(self):
        if self.review_test_dir.is_dir():
                shutil.rmtree(str(self.review_test_dir))

    def test_review_1(self):
        """Verify execute_review() creates new directory as expected.
        """
        review.execute_review(self.alchemist, self.test_dir, 
                              self.review_test_dir.name)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_review_2(self):
        """Verify execute_review() filter parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name, 
                              filters=("phage.Cluster='A' "
                                       "AND phage.Subcluster='A2'"))

        self.assertTrue(self.review_test_dir.is_dir())

    def test_review_3(self):
        """Verify execute_review() group parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              groups=["phage.Cluster"])

        self.assertTrue(self.review_test_dir.is_dir())

        clusterA_dir = self.review_test_dir.joinpath("A")
        self.assertTrue(clusterA_dir.is_dir())

    def test_review_4(self):
        """Verify execute_review() sort parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              sort=["gene.Name"])

        self.assertTrue(self.review_test_dir.is_dir())

    def test_review_5(self):
        """Verify execute_review() review parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              review=False)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_review_6(self):
        """Verify execute_review() pg_report parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              pg_report=True)

        self.assertTrue(self.review_test_dir.is_dir())

        gene_report_dir = self.review_test_dir.joinpath("GeneReports")
        self.assertTrue(gene_report_dir.is_dir())

if __name__ == "__main__":
    unittest.main()
