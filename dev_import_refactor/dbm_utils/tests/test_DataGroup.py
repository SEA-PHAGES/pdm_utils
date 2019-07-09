""" Unit tests for the DataGroup Class."""


from classes import DataGroup
from classes import Genome
from classes import GenomePair
from classes import Ticket
from classes import Eval
import unittest


class TestDataGroupClass(unittest.TestCase):


    def setUp(self):

        self.data_group = DataGroup.DataGroup()
        self.genome1 = Genome.Genome()
        self.genome1.type = "flat_file"
        self.genome2 = Genome.Genome()
        self.genome2.type = "phamerator"
        self.ticket = Ticket.GenomeTicket()






    def test_set_genome_pair_1(self):
        """Check that a genome pair is set if both keys are present."""

        self.data_group.ticket = self.ticket
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.genome_dict[self.genome2.type] = self.genome2
        genome_pair = GenomePair.GenomePair()
        self.data_group.set_genome_pair(genome_pair, "phamerator", "flat_file")
        self.assertEqual(list(self.data_group.genome_pair_dict.keys())[0],
                            "phamerator_flat_file")

    def test_set_genome_pair_2(self):
        """Check that a genome pair is not set if one key is not present."""

        self.data_group.ticket = self.ticket
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.genome_dict[self.genome2.type] = self.genome2
        genome_pair = GenomePair.GenomePair()
        self.data_group.set_genome_pair(genome_pair, "invalid", "flat_file")
        self.assertEqual(len(self.data_group.genome_pair_dict.keys()), 0)







    def test_check_matched_genome_1(self):
        """Check that no error is produced when the genome type is present."""
        self.data_group.ticket = self.ticket
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.genome_dict[self.genome2.type] = self.genome2
        self.data_group.check_matched_genome("phamerator")
        self.assertEqual(self.data_group.evaluations[0].status, "correct")

    def test_check_matched_genome_2(self):
        """Check that an error is produced when the genome type is present."""
        self.data_group.ticket = self.ticket
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.genome_dict[self.genome2.type] = self.genome2
        self.data_group.check_matched_genome("invalid")
        self.assertEqual(self.data_group.evaluations[0].status, "error")


























if __name__ == '__main__':
    unittest.main()
