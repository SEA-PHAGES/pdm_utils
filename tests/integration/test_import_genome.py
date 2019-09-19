"""Integration tests for the main import pipeline."""



import unittest
import os
import shutil
from pdm_utils.pipelines.db_import import import_genome
from pdm_utils.constants import constants
from pdm_utils.classes import bundle, genome, ticket
from pdm_utils.classes import mysqlconnectionhandler as mch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference



# The following integration tests user the 'tester' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "tester"
pwd = "tester"
db = "test_db"


class TestImportGenomeMain(unittest.TestCase):


    def setUp(self):


        self.base_dir = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/test_import")
        os.mkdir(self.base_dir)

        self.genome_folder = \
            os.path.join(self.base_dir, "genome_folder")
        os.mkdir(self.genome_folder)

        self.test_import_table = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_import_table_1.csv")


        self.sql_handle = mch.MySQLConnectionHandler()



    def test_prepare_tickets_1(self):
        """Verify dictionary is returned from a correct import table file."""
        eval_flags=constants.RUN_MODE_PHAGESDB
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table,
                        eval_flags=constants.RUN_MODE_PHAGESDB,
                        description_field="product")
        self.assertEqual(len(tkt_dict.keys()), 2)







    # TODO finish after building tests for other module functions.
    # def test_import_io_1(self):
    #     """."""
    #
    #     import_genome.import_io(sql_handle=self.sql_handle,
    #                             genome_folder=self.genome_folder,
    #                             import_table_file=self.test_import_table,
    #                             filename_flag=False, test_run=True,
    #                             description_field="product",
    #                             run_mode="phagesdb")



    def tearDown(self):
        shutil.rmtree(self.base_dir)






if __name__ == '__main__':
    unittest.main()
