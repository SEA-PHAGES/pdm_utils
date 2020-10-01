

from pathlib import Path
import shutil
import sys
import unittest

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.pipelines import pham_review

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils


USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
TEST_DIR = "/tmp/pdm_utils_tests_pham_review"


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
        self.review_test_dir = self.test_dir.joinpath("pham_review_test_dir")

        self.alchemist = AlchemyHandler()
        self.alchemist.username = USER
        self.alchemist.password = PWD
        self.alchemist.database = DB
        self.alchemist.connect(ask_database=True, login_attempts=0)

        self.db_filter = Filter(alchemist=self.alchemist)
        self.db_filter.add(pham_review.BASE_CONDITIONALS)
        self.db_filter.key = "gene.PhamID"

    def tearDown(self):
        if self.review_test_dir.is_dir():
            shutil.rmtree(str(self.review_test_dir))

    def test_execute_pham_review_1(self):
        """Verify execute_pham_review() creates new directory as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_pham_review_2(self):
        """Verify execute_pham_review() filter parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   filters=("phage.Cluster='A' "
                                            "AND phage.Subcluster='A2'"))

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_pham_review_3(self):
        """Verify execute_pham_review() group parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   groups=["phage.Cluster"])

        self.assertTrue(self.review_test_dir.is_dir())

        clusterA_dir = self.review_test_dir.joinpath("A")
        self.assertTrue(clusterA_dir.is_dir())

    def test_execute_pham_review_4(self):
        """Verify execute_pham_review() sort parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   sort=["gene.Name"])

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_pham_review_5(self):
        """Verify execute_pham_review() review parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   no_review=True)

        self.assertTrue(self.review_test_dir.is_dir())

    def test_execute_pham_review_6(self):
        """Verify execute_pham_review() gr_reports parameter functions
        as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   gr_reports=True)

        self.assertTrue(self.review_test_dir.is_dir())

        pham_report_dir = self.review_test_dir.joinpath("PhamReports")
        self.assertTrue(pham_report_dir.is_dir())

    def test_execute_pham_review_7(self):
        """Verify execute_pham_review() s_report parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   s_report=True)

        self.assertTrue(self.review_test_dir.is_dir())

        summary_report_file = self.review_test_dir.joinpath(
                                                        "SummaryReport.txt")
        self.assertTrue(summary_report_file.is_file())

    def test_execute_pham_review_8(self):
        """Verify execute_pham_review() psr_reports parameter functions as expected.
        """
        pham_review.execute_pham_review(
                                   self.alchemist, folder_path=self.test_dir,
                                   folder_name=self.review_test_dir.name,
                                   psr_reports=True)

        self.assertTrue(self.review_test_dir.is_dir())

        pham_report_dir = self.review_test_dir.joinpath("PhamReports")
        self.assertTrue(pham_report_dir.is_dir())

    def test_review_phams_1(self):
        """Verify review_phams() correctly identifies disrepencies.
        """
        self.db_filter.values = self.db_filter.build_values(
                        where=self.db_filter.build_where_clauses())

        pham_review.review_phams(self.db_filter)

        self.assertFalse(39854 in self.db_filter.values)
        self.assertTrue(40481 in self.db_filter.values)

    def test_get_review_data_1(self):
        """Verify get_review_data() retrieves and returns data as expected.
        """
        self.db_filter.values = [40481]

        review_data = pham_review.get_review_data(
                                                self.alchemist, self.db_filter)

        self.assertTrue(isinstance(review_data, list))

        for header in pham_review.REVIEW_HEADER:
            with self.subTest(header=header):
                self.assertTrue(header in review_data[0].keys())
                self.assertFalse(isinstance(review_data[0][header], list))

    def test_get_summary_data_1(self):
        """Verify get_summary_data() retrieves and returns data as expected.
        """
        self.db_filter.values = [40481]

        summary_data = pham_review.get_summary_data(
                                                self.alchemist, self.db_filter)

        self.assertTrue(isinstance(summary_data, dict))

    def test_get_gr_data_1(self):
        """Verify get_g_data() retrieves and returns data as expected.
        """
        self.db_filter.values = [40481]

        gr_data = pham_review.get_gr_data(
                                        self.alchemist, self.db_filter, 40481)

        self.assertTrue(isinstance(gr_data, list))

        for header in pham_review.GR_HEADER:
            with self.subTest(header=header):
                self.assertTrue(header in gr_data[0].keys())
                self.assertFalse(isinstance(gr_data[0][header], list))


if __name__ == "__main__":
    unittest.main()
