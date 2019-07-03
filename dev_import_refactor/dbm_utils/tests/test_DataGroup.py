""" Unit tests for the DataGroup Class."""


from classes import DataGroup
from classes import Genome
from classes import Ticket
from classes import Eval
import unittest


class TestDataGroupClass(unittest.TestCase):


    def setUp(self):

        self.data_group = DataGroup.DataGroup()
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.ticket = Ticket.ImportTicket()
        self.dict_name = "misc"
        self.data_group.evaluations_dict[self.dict_name] = []




    # TODO this are probably no longer needed.
    # def test_set_evaluation_1(self):
    #     """Check that empty EvalResult is appended to a evaluation list
    #     in the dictionary of evaluations."""
    #
    #     self.data_group.set_evaluation(self.dict_name, "none")
    #     self.assertEqual( \
    #         len(self.data_group.evaluations_dict[self.dict_name]), 1)
    #
    # def test_set_evaluation_2(self):
    #     """Check that 'warning' EvalResult is appended to a evaluation list
    #     in the dictionary of evaluations."""
    #
    #     self.data_group.set_evaluation(self.dict_name, "warning", "message1")
    #     self.assertEqual( \
    #         len(self.data_group.evaluations_dict[self.dict_name]), 1)
    #
    # def test_set_evaluation_3(self):
    #     """Check that 'error' EvalResult is appended to a evaluation list
    #     in the dictionary of evaluations."""
    #
    #     self.data_group.set_evaluation(self.dict_name, "error", "message1", "message2")
    #     self.assertEqual( \
    #         len(self.data_group.evaluations_dict[self.dict_name]), 1)






























if __name__ == '__main__':
    unittest.main()
