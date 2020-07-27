import csv
import filecmp
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import Entrez, SeqIO

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.filter import Filter
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

def get_acc_id_dict(alchemist):
    """Test helper function to retrieve accessions of database entries.
    """
    db_filter = Filter(alchemist=alchemist)
    db_filter.key = "phage.PhageID"

    db_filter.values = db_filter.build_values() 
    groups = db_filter.group("phage.Accession")

    return groups

def copy_gb_ft_files(ncbi_handle, acc_id_dict, records_path):
    """Test helper functions to create duplicates of GenBank ft files.
    """
    file_lines = ncbi_handle.readlines()

    feature_format = re.compile(">Feature .+\|(\w+)(\..)\|\n")

    total_files_path = records_path.joinpath("TOTAL.tbl")
    total_files_handle = total_files_path.open(mode="w")
    
    file_handle = None
    for line in file_lines:
        total_files_handle.write(line) 

        if not re.match(feature_format, line) is None:
            if not file_handle is None:
                file_handle.close()

            accession_split = re.split(feature_format, line)
            accession = accession_split[1]
            phage_id = acc_id_dict[accession_split[1]][0]

            file_name = (f"{phage_id}.tbl") 
            file_path = records_path.joinpath(file_name)
            if file_path.is_file():
                raise Exception

            file_handle = file_path.open(mode="w")

            file_handle.write(line)
        else:
            if file_handle is None:
                continue
            file_handle.write(line) 

    if not file_handle is None:
        file_handle.close()

    total_files_handle.close()
    ncbi_handle.close()

class TestFeatureTableParser(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        base_dir = Path(TMPDIR_BASE)
        self.test_dir = base_dir.joinpath(TMPDIR_PREFIX)
        test_db_utils.create_filled_test_db()
       
        if self.test_dir.is_dir():
            shutil.rmtree(self.test_dir)

        self.test_dir.mkdir()

        self.alchemist = AlchemyHandler()
        self.alchemist.username=USER
        self.alchemist.password=PWD
        self.alchemist.database=DB
        self.alchemist.connect(ask_database=True, login_attempts=0)
        self.acc_id_dict = get_acc_id_dict(self.alchemist)

        accession_list = list(self.acc_id_dict.keys())
        ncbi_handle = Entrez.efetch(db="nucleotide", rettype="ft",
                                    id=",".join(accession_list), retmode="text")

        copy_gb_ft_files(ncbi_handle, self.acc_id_dict, self.test_dir)

    @classmethod 
    def tearDownClass(self):
        test_db_utils.remove_db()
        shutil.rmtree(self.test_dir)

    def setUp(self):
        self.indv_test_dir = self.test_dir.joinpath("RewrittenFeatureTables")
        self.indv_test_dir.mkdir()

    def tearDown(self):
        shutil.rmtree(self.indv_test_dir)

    def test_read_feature_table_1(self):
        """Verify read_feature_table() can read a number of differing tbl files.
        """
        for tbl_file in self.test_dir.iterdir():
            if tbl_file.is_file():
                file_name = tbl_file.name
                if file_name == "TOTAL.tbl":
                    continue

                with self.subTest(file_name=file_name):
                    with tbl_file.open(mode="r") as filehandle:
                        record = fileio.read_feature_table(filehandle)
                        
                        self.assertFalse(record is None)
                        self.assertTrue(record.id in self.acc_id_dict.keys())
                        
                        record_name = self.acc_id_dict[record.id][0]
                        self.assertEqual(".".join([record_name, "tbl"]), 
                                                   file_name)

    def test_read_feature_table_2(self):
        """Verify read_feature_table() returns None from an empty file.
        """
        empty_file = self.indv_test_dir.joinpath("empty.tbl")
        with empty_file.open(mode="w") as filehandle:
            filehandle.write("Nonsense\n")

        with empty_file.open(mode="r") as filehandle:
            record = fileio.read_feature_table(filehandle)

        self.assertTrue(record is None)

    def test_parse_feature_table_1(self):
        """Verify parse_feature_table() can read a concatenated feature table.
        """
        total_file_path = self.test_dir.joinpath("TOTAL.tbl")
        with total_file_path.open(mode="r") as filehandle:
            records = fileio.parse_feature_table(filehandle)

            for record in records:
                with self.subTest(accession=record.id):
                    self.assertFalse(record is None)
                    self.assertTrue(record.id in self.acc_id_dict.keys())
                        
                    record_name = self.acc_id_dict[record.id][0]
                    record_path = self.test_dir.joinpath(
                                         ".".join([record_name, "tbl"]))

                    self.assertTrue(record_path.is_file())

    def test_parse_write_feature_table_1(self):
        """Verify reading the feature table in and rewriting it exactly mimics 
        copied feature table files.
        """
        total_file_path = self.test_dir.joinpath("TOTAL.tbl")  

        records = []
        file_names = []
        with total_file_path.open(mode="r") as filehandle:
            parser = fileio.parse_feature_table(filehandle)
            for record in parser:
                record_name = self.acc_id_dict[record.id][0]
                record.name = record_name
               
                file_names.append(".".join([record_name, "tbl"]))
                
                records.append(record)

        fileio.write_feature_table(records, self.indv_test_dir)

        file_diffs = filecmp.cmpfiles(self.test_dir, self.indv_test_dir, 
                                      file_names)

        for discrepant_file in file_diffs[1]:
            print(f"TEST ERROR: data in {discrepant_file} is incorrect:")
            source_file = self.test_dir.joinpath(discrepant_file)
            rewritten_file = self.indv_test_dir.joinpath(discrepant_file)

            source_handle = source_file.open(mode="r")
            rewritten_handle = rewritten_file.open(mode="r")

            source_data = source_handle.readlines()
            rewritten_data = rewritten_handle.readlines()
           
            source_handle.close()
            rewritten_handle.close()

            for i in range(len(source_data)):
                source_line = source_data[i]
                print(source_line.rstrip("\n"))
                self.assertTrue(i <= len(rewritten_data) - 1)

                rewritten_line = rewritten_data[i]

                self.assertEqual(source_line, rewritten_line)

        self.assertTrue(len(file_diffs[1]) == 0)

if __name__ == "__main__":
    unittest.main()
