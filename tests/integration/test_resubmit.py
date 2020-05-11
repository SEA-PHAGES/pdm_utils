from pathlib import Path
import shutil
import sys 
import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
from pdm_utils.functions import basic
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

TEST_DATA = [{"Pham" : 40481,"#Members" : 4,"Clusters" : "A", "#Functions" : 3,
              "Functional Calls" : "Hypothetical Protein;terminase;Terminase", 
              "Final Call" : "terminase"},
             {"Pham" : 25050,"#Members" : 4,"Clusters" : "A;N;None",
              "#Functions" : 2, 
              "Functional Calls" : "Hypothetical Protein;minor tail protein",
              "Final Call" : "Hypothetical Protein"},
             {"Pham" : 40880,"#Members" : 4,"Clusters" : "A;N","#Functions" : 3, 
              "Functional Calls" : "Hypothetical Protein;Holin;holin",
              "Final Call" : "holin"},
             {"Pham" : 39529,"#Members" : 5,"Clusters" : "A;None", 
              "#Functions" : 5, 
              "Functional Calls" : ("Hypothetical Protein;endonuclease;"
                                    "Endo VII;EndoVII;endonuclease VII"),
              "Final Call" : "endonuclease VII"}]


class TestGenbankResubmit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

        self.test_dir = Path(TEST_DIR)
        if self.test_dir.is_dir():
            shutil.rmtree(TEST_DIR)

        self.test_dir.mkdir()
        self.resubmit_form = self.test_dir.joinpath("resubmit_form.txt")

        basic.export_data_dict(TEST_DATA, self.resubmit_form, PF_HEADER, 
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
        self.alchemist.connect(ask_database=True, login_attempts=0)

        self.resubmit_test_dir = self.test_dir.joinpath("resubmit_test_dir")

    def tearDown(self):
        if self.resubmit_test_dir.is_dir():
            shutil.rmtree(str(self.resubmit_test_dir))

    def test_execute_resubmit_1(self):
        """Verify execute_resubmit creates new directory as expected.
        """
        resubmit.execute_resubmit(self.alchemist, TEST_DATA, self.test_dir,
                                    self.resubmit_test_dir.name)

        self.assertTrue(self.resubmit_test_dir.is_dir())
       
    
if __name__ == "__main__":
    unittest.main()
