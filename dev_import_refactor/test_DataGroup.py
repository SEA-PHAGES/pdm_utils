""" Unit tests for the DataGroup Class."""


import DataGroup
import Genome
import Ticket
import Eval
import unittest


class TestDataGroupClass(unittest.TestCase):


    def setUp(self):

        self.data_group = DataGroup.DataGroup()
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.ticket = Ticket.ImportTicket()




    def test_set_evaluation_1(self):
        self.data_group.set_evaluation("none")
        self.assertEqual(len(self.data_group.evaluations), 1)

    def test_set_evaluation_2(self):
        self.data_group.set_evaluation("warning","message1")
        self.assertEqual(len(self.data_group.evaluations), 1)

    def test_set_evaluation_3(self):
        self.data_group.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.data_group.evaluations), 1)




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
