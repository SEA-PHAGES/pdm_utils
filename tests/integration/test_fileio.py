import csv
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions import fileio

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

# pdm_anon, pdm_anon, and pdm_test_db
USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB

unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
test_file_dir = Path(test_dir, "test_files")

TMPDIR_PREFIX = "pdm_utils_tests_fileio_"
# Can set TMPDIR_BASE to string such as "/tmp/" to track tmp directory location.
TMPDIR_BASE = "/tmp"

class TestFileIO(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        base_dir = Path(TMPDIR_BASE)
        self.test_dir = base_dir.joinpath(TMPDIR_PREFIX)
        test_db_utils.create_filled_test_db()

        if self.test_dir.is_dir():
            shutil.rmtree(self.test_dir)

        self.test_dir.mkdir()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()
        shutil.rmtree(self.test_dir)

    def setUp(self):
        self.alchemist = AlchemyHandler()
        self.alchemist.username=USER
        self.alchemist.password=PWD
        self.alchemist.database=DB
        self.alchemist.connect(ask_database=True, login_attempts=0)
        self.alchemist.build_graph()

        self.test_import_table_1 = Path(test_file_dir, "test_import_table_1.csv")
        self.tkt_dict1 = {"phage_id": "L5", "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"phage_id": "Trixie", "host_genus": "Mycobacterium"}

        self.fasta_dict_1 = {"Trixie_CDS_11" : ("MASIQGKLIALVLKYGISYLRKHPELLKEI"
                                                "SKHIPGKVDDLVLEVLAKLLGV")}
        self.fasta_dict_2 = {"TRIXIE_CDS_3" : ("MSGFDDKIVDQAQAIVPADDYDALPLAGPGR"
                                               "WAHVPGGLTLYTNDDTVLFAQGDMSTIESSY"
                                               "LFQAMEKLRLAGKTASQAFDILRLEADAISG"
                                               "DLSELAEE"),
                             "L5_CDS_3"     : ("MAQMQATHTIEGFLAVEVAPRAFVAENGHVL"
                                               "TRLSATKWGGGEGLEILNYEGPGTVEVSDEK"
                                               "LAEAQRASEVEAELRREVGKE")}
        
        self.fileio_test_dir = self.test_dir.joinpath("fileio_test_dir")
        self.fileio_test_dir.mkdir()
        self.data_dict_file =  self.fileio_test_dir.joinpath("table.csv")
        self.fasta_file = self.fileio_test_dir.joinpath("translations.fasta")

    def tearDown(self):
        shutil.rmtree(self.fileio_test_dir)

    def test_retrieve_data_dict_1(self):
        """Verify a correctly structured file can be opened."""
        list_of_data_dicts = \
            fileio.retrieve_data_dict(self.test_import_table_1)
        self.assertEqual(len(list_of_data_dicts), 2)

    def test_export_data_dict_1(self):
        """Verify data is exported correctly."""

        list_of_data = [self.tkt_dict1, self.tkt_dict2]
        headers = ["type", "phage_id", "host_genus", "cluster"]
        fileio.export_data_dict(list_of_data, self.data_dict_file,
                                    headers, include_headers=True)

        exp_success_tkts = []
        with open(self.data_dict_file,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_success_tkts.append(dict)

        with self.subTest():
            self.assertEqual(len(exp_success_tkts), 2)
        with self.subTest():
            self.assertEqual(set(exp_success_tkts[0].keys()), set(headers))

    def test_write_fasta_1(self): 
        """Verify write_fasta() creates readable fasta formatted file"""
        fileio.write_fasta(self.fasta_dict_1, self.fasta_file)

        record = SeqIO.read(self.fasta_file, "fasta")
        id = list(self.fasta_dict_1.keys())[0]
        seq = self.fasta_dict_1[id]

        self.assertEqual(record.id, id)
        self.assertEqual(str(record.seq), seq)
           
    def test_write_fasta_2(self):
        """Verify write_fasta() can properly concatenate fasta files"""
        fileio.write_fasta(self.fasta_dict_2, self.fasta_file)

        records = SeqIO.parse(self.fasta_file, "fasta")
        
        keys = list(self.fasta_dict_2.keys())
        
        for record in records:
            self.assertTrue(record.id in keys)

            seq = self.fasta_dict_2[record.id]
            self.assertEqual(str(record.seq), seq)

if __name__ == "__main__":
    unittest.main()
