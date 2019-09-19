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
            os.path.join(os.path.dirname(__file__),
            "test_files/test_import_table_1.csv")

        self.test_flat_file1 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_1.gb")

        self.test_flat_file2 = \
            os.path.join(os.path.dirname(__file__),
            "test_files/test_flat_file_2.gb")

        self.sql_handle = mch.MySQLConnectionHandler()


        self.eval_flags = {
            "check_locus_tag":True,
            "check_description_field":True,
            "check_replace":True,
            "check_trna":True,
            "import_locus_tag":True,
            "check_id_typo":True,
            "check_host_typo":True,
            "check_author":True,
            "check_description":True,
            "check_gene":True
            }
        self.tkt1 = ticket.GenomeTicket()
        self.tkt1.id = 1
        self.tkt1.type = "add"
        self.tkt1.phage_id = "L5"
        self.tkt1.run_mode = "phagesdb"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.host_genus = "Mycobacterium"
        self.tkt1.cluster = "A"
        self.tkt1.subcluster = "A2"
        self.tkt1.annotation_status = "draft"
        self.tkt1.annotation_author = 1
        self.tkt1.annotation_qc = 1
        self.tkt1.retrieve_record = 1
        self.tkt1.accession = "ABC123"

        self.tkt2 = ticket.GenomeTicket()


    def test_prepare_tickets_1(self):
        """Verify dictionary is returned from a correct import table file."""
        eval_flags=constants.RUN_MODE_PHAGESDB
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table,
                        eval_flags=constants.RUN_MODE_PHAGESDB,
                        description_field="product")
        self.assertEqual(len(tkt_dict.keys()), 2)










    def test_prepare_bundle_1(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no phagesdb data."""
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1)
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["add"]
        bndl_tkt = bndl.ticket
        ff_tkt_pair = bndl.genome_pair_dict["flat_file_add"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(bndl.id, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "L5")






    def test_prepare_bundle_2(self):
        """Verify bundle is returned from a flat file with:
        no record."""
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file2,
                    ticket_dict=tkt_dict, id=1)
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)




    def test_prepare_bundle_3(self):
        """Verify bundle is returned from a flat file with:
        one record, no ticket."""
        tkt_dict = {"L5x":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1)
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)




    def test_prepare_bundle_4(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, with phagesdb data."""
        self.tkt1.host_genus = "retrieve"
        self.tkt1.cluster = "B"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(filename=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1)
        ff_gnm = bndl.genome_dict["flat_file"]
        pdb_gnm = bndl.genome_dict["phagesdb"]
        ff_tkt_pair = bndl.genome_pair_dict["flat_file_add"]
        ff_pdb_pair = bndl.genome_pair_dict["flat_file_phagesdb"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")


    # TODO continue unit testing prepare_bundle().
    # Next test: use a 'replace' ticket.









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
