"""Integration tests for misc. ticket functions"""

import unittest
import os
import csv
from pdm_utils.functions import tickets
class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.test_import_table_1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_import_table_1.csv")

        # self.test_import_table_2 = \
        #     os.path.join(os.path.dirname(__file__), \
        #     "test_files/test_import_table_2.csv")

    def test_retrieve_ticket_data_1(self):
        """Verify a correctly structured file can be opened."""
        list_of_ticket_data = \
            tickets.retrieve_ticket_data(self.test_import_table_1)
        self.assertEqual(len(list_of_ticket_data), 2)

    #
    # def test_retrieve_import_ticket_data_2(self):
    #     """Verify a correctly structured file can be opened."""
    #     list_of_ticket_data = \
    #         tickets.retrieve_import_ticket_data(self.test_import_table_2)
    #     self.assertEqual(len(list_of_ticket_data), 4)
    #     print(list_of_ticket_data[0])
    #     print(list_of_ticket_data[1])
    #     print(list_of_ticket_data[2])
    #     print(list_of_ticket_data[3])
    #     print(list_of_ticket_data[4])

if __name__ == '__main__':
    unittest.main()
