"""Integration tests for misc. ticket functions"""

import unittest
import os
import csv
from pdm_utils.functions import tickets
from pathlib import Path
import shutil

class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):

        self.unittest_file = Path(__file__)
        self.unittest_dir = self.unittest_file.parent
        self.test_import_table_1 = Path(self.unittest_dir,
                                     "test_files/test_import_table_1.csv")
        self.base_dir = Path(self.unittest_dir, "test_files/test_tickets")
        self.base_dir.mkdir()

        self.tkt_dict1 = {"phage_id": "L5", "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"phage_id": "Trixie", "host_genus": "Mycobacterium"}

        self.export_file = Path(self.base_dir, "table.csv")


    def tearDown(self):
        shutil.rmtree(self.base_dir)



    def test_retrieve_ticket_data_1(self):
        """Verify a correctly structured file can be opened."""
        list_of_ticket_data = \
            tickets.retrieve_ticket_data(self.test_import_table_1)
        self.assertEqual(len(list_of_ticket_data), 2)




    def test_export_ticket_data_1(self):
        """Verify data is exported correctly."""

        list_of_data = [self.tkt_dict1, self.tkt_dict2]
        headers = ["type", "phage_id", "host_genus", "cluster"]
        tickets.export_ticket_data(list_of_data, self.export_file,
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




if __name__ == '__main__':
    unittest.main()
