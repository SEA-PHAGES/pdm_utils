import unittest
from pathlib import Path
from unittest import Mock
from unittest import patch
from unittest import PropertyMock

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.pipelines import resubmit
from pdm_utils.pipelines.review import PF_HEADER

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
TEST_DIR = "/tmp/pdm_utils_tests_resubmit"

TEST_DATA = [{"Pham" : ,"#Members" : ,"Clusters" : "#Functions" : , 
              "Functional Calls" : ,"Final Call" : },
             {"Pham" : ,"#Members" : ,"Clusters" : "#Functions" : , 
              "Functional Calls" : ,"Final Call" : },
             {"Pham" : ,"#Members" : ,"Clusters" : "#Functions" : , 
              "Functional Calls" : ,"Final Call" : },
             {"Pham" : ,"#Members" : ,"Clusters" : "#Functions" : , 
              "Functional Calls" : ,"Final Call" : }]


class TestGenbankResubmit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

        self.test_dir = Path(TEST_DIR)
        if self.test_dir.is_dir():
            shutil.rmtree(TEST_DIR)

        self.test_dir.mkdir()
        self.resubmit_form = self.test_dir.joinpath("resubmit_form.txt")

        basic.export_data_dicts(TEST_DATA, self.resubmit_form, PF_HEADERS, 
                                                        include_headers=True)

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()
        shutil.rmtree(TEST_DIR)

    def setUp(self):
        self.alchemist = AlchemyHandler()
        self.alchemist.username = USER
        self.alchemist.password = PWD
        self.alchemist.database = DB
        self.alchemist.connect(ask_database=True, login_attemptes=0)

        self.export_test_dir = self.test_dir.joinpath("export_test_dir")

    def tearDown(self):
        if self.export_test_dir.is_dir():
            shutil.rmtree(str(self.export_test_dir))

    def test_execute_resubmit_1(self):
        """Verify execute_resubmit creates new directory as expected.
        """
       
def write_test_resubmit_form(export_path, file_path):
    data_dicts = [{}] 

    

    

if __name__ == "__main__":
    unittest.main()
