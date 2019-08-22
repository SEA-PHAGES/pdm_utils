""" Unit tests for the CdsPair Class."""


from classes import CdsPair
from classes import Eval
import unittest


class TestCdsPairClass(unittest.TestCase):


    def setUp(self):
        self.cds_pair = CdsPair.CdsPair()
        self.cds1 = Cds.Cds()
        self.cds2 = Cds.Cds()



    #
    #
    # def test_set_evaluation_1(self):
    #     """Set an empty evaluation object."""
    #     self.cds_pair.set_evaluation("none")
    #     self.assertEqual(len(self.cds_pair.evaluations), 1)
    #
    # def test_set_evaluation_2(self):
    #     """Set a warning evaluation object."""
    #     self.cds_pair.set_evaluation("warning","message1")
    #     self.assertEqual(len(self.cds_pair.evaluations), 1)
    #
    # def test_set_evaluation_3(self):
    #     """Set an error evaluation object."""
    #     self.cds_pair.set_evaluation("error","message1","message2")
    #     self.assertEqual(len(self.cds_pair.evaluations), 1)






























if __name__ == '__main__':
    unittest.main()
