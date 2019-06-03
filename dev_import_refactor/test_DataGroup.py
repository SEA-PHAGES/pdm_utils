""" Unit tests for the DataGroup Class."""


import DataGroup
import unittest


class TestDataGroupClass(unittest.TestCase):


    def setUp(self):
        self.normal_ticket_list = ["add",
                                "Trixe",
                                "Mycobacterium",
                                "A",]








    # def test_parse_import_ticket_2(self):
    #     """Verify improperly structured data is not parsed."""
    #     ticket, eval = \
    #         FunctionsTicket.parse_import_ticket(self.short_ticket_list)
    #     with self.subTest():
    #         self.assertIsNotNone(eval)
    #     with self.subTest():
    #         self.assertIsNone(ticket)
    #     with self.subTest():
    #         self.assertEqual(eval.status, "error")



























if __name__ == '__main__':
    unittest.main()
