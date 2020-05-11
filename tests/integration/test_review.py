from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
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
        self.review_test_dir = self.test_dir.joinpath("review_test_dir")

        self.alchemist = AlchemyHandler()
        self.alchemist.username=USER
        self.alchemist.password=PWD
        self.alchemist.database=DB
        self.alchemist.connect(ask_database=True, login_attempts=0)

        self.db_filter = Filter(alchemist=self.alchemist)
        self.db_filter.add(review.BASE_CONDITIONALS)
        self.db_filter.key = "gene.PhamID"

    def tearDown(self):
        if self.review_test_dir.is_dir():
                shutil.rmtree(str(self.review_test_dir))

    def test_execute_review_1(self):
        """Verify execute_review() creates new directory as expected.
        """
        review.execute_review(self.alchemist, self.test_dir, 
                              self.review_test_dir.name)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_review_2(self):
        """Verify execute_review() filter parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name, 
                              filters=("phage.Cluster='A' "
                                       "AND phage.Subcluster='A2'"))

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_review_3(self):
        """Verify execute_review() group parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              groups=["phage.Cluster"])

        self.assertTrue(self.review_test_dir.is_dir())

        clusterA_dir = self.review_test_dir.joinpath("A")
        self.assertTrue(clusterA_dir.is_dir())

    def test_execute_review_4(self):
        """Verify execute_review() sort parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              sort=["gene.Name"])

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_review_5(self):
        """Verify execute_review() review parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              review=False)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_review_6(self):
        """Verify execute_review() pg_report parameter functions as expected.
        """
        review.execute_review(self.alchemist, self.test_dir,
                              self.review_test_dir.name,
                              pg_report=True)

        self.assertTrue(self.review_test_dir.is_dir())

        gene_report_dir = self.review_test_dir.joinpath("GeneReports")
        self.assertTrue(gene_report_dir.is_dir())

    def test_review_phams_1(self):
        """Verify review_phams() correctly identifies disrepencies.
        """
        self.db_filter.values = self.db_filter.build_values(
                          where=self.db_filter.build_where_clauses())

        review.review_phams(self.db_filter)

        self.assertFalse(39854 in self.db_filter.values)
        self.assertTrue(40481 in self.db_filter.values)

    def test_get_pf_data_1(self):
        """Verify get_pf_data() retrieves and returns data as expected.
        """
        self.db_filter.values = [40481]

        pf_data = review.get_pf_data(self.alchemist, self.db_filter)

        self.assertTrue(isinstance(pf_data, list))
       
        for header in review.PF_HEADER:
            with self.subTest(header=header):
                self.assertTrue(header in pf_data[0].keys())
                self.assertFalse(isinstance(pf_data[0][header], list))

    def test_get_pg_data_1(self):
        """Verify get_pg_data() retreives and retrusn data as expected.
        """
        self.db_filter.values = [40481]

        pg_data = review.get_pg_data(self.alchemist, self.db_filter, 40481)

        self.assertTrue(isinstance(pg_data, list))
       
        for header in review.PG_HEADER:
            with self.subTest(header=header):
                self.assertTrue(header in pg_data[0].keys())
                self.assertFalse(isinstance(pg_data[0][header], list))
    
    
if __name__ == "__main__":
    unittest.main()
