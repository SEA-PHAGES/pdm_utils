import csv
import tempfile
import unittest
from pathlib import Path

from pdm_utils.functions import fileio

unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
test_file_dir = Path(test_dir, "test_files")

TMPDIR_PREFIX = "pdm_utils_tests_fileio_"
# Can set TMPDIR_BASE to string such as "/tmp/" to track tmp directory location.
TMPDIR_BASE = "/tmp"

class TestFileIO(unittest.TestCase):


    def setUp(self):
        self.test_import_table_1 = Path(test_file_dir, "test_import_table_1.csv")
        self.tkt_dict1 = {"phage_id": "L5", "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"phage_id": "Trixie", "host_genus": "Mycobacterium"}
        self.tmpdir = tempfile.TemporaryDirectory(prefix=TMPDIR_PREFIX,
                                                  dir=TMPDIR_BASE)
        self.base_dir = Path(self.tmpdir.name)
        self.export_file = Path(self.base_dir, "table.csv")

    def tearDown(self):
        self.tmpdir.cleanup()




    def test_retrieve_data_dict_1(self):
        """Verify a correctly structured file can be opened."""
        list_of_data_dicts = \
            fileio.retrieve_data_dict(self.test_import_table_1)
        self.assertEqual(len(list_of_data_dicts), 2)




    def test_export_data_dict_1(self):
        """Verify data is exported correctly."""

        list_of_data = [self.tkt_dict1, self.tkt_dict2]
        headers = ["type", "phage_id", "host_genus", "cluster"]
        fileio.export_data_dict(list_of_data, self.export_file,
                                    headers, include_headers=True)

        exp_success_tkts = []
        with open(self.export_file,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_success_tkts.append(dict)

        with self.subTest():
            self.assertEqual(len(exp_success_tkts), 2)
        with self.subTest():
            self.assertEqual(set(exp_success_tkts[0].keys()), set(headers))
