""" Unit tests for the GenomePair Class."""


import GenomePair
import Genome
import Ticket
import Eval
import unittest


class TestGenomePairClass(unittest.TestCase):


    def setUp(self):
        self.genome_pair = GenomePair.GenomePair()
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.ticket = Ticket.ImportTicket()










    def test_set_evaluation_1(self):
        """Set an empty evaluation object."""
        self.genome_pair.set_evaluation("none")
        self.assertEqual(len(self.genome_pair.evaluations), 1)

    def test_set_evaluation_2(self):
        """Set a warning evaluation object."""
        self.genome_pair.set_evaluation("warning","message1")
        self.assertEqual(len(self.genome_pair.evaluations), 1)

    def test_set_evaluation_3(self):
        """Set an error evaluation object."""
        self.genome_pair.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.genome_pair.evaluations), 1)






























if __name__ == '__main__':
    unittest.main()
