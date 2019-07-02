""" Unit tests for the GenomePair Class."""


from classes import GenomePair
from classes import Genome
from classes import Ticket
from classes import Eval
import unittest


class TestGenomePairClass(unittest.TestCase):


    def setUp(self):
        self.genome_pair = GenomePair.GenomePair()
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.ticket = Ticket.ImportTicket()









    # TODO this is probably no longer needed.
    # def test_set_evaluation_1(self):
    #     """Set an empty evaluation object."""
    #     self.genome_pair.set_evaluation("none")
    #     self.assertEqual(len(self.genome_pair.evaluations), 1)
    #
    # def test_set_evaluation_2(self):
    #     """Set a warning evaluation object."""
    #     self.genome_pair.set_evaluation("warning","message1")
    #     self.assertEqual(len(self.genome_pair.evaluations), 1)
    #
    # def test_set_evaluation_3(self):
    #     """Set an error evaluation object."""
    #     self.genome_pair.set_evaluation("error","message1","message2")
    #     self.assertEqual(len(self.genome_pair.evaluations), 1)




    def test_compare_genome_sequence_1(self):
        """Check that identical sequences produce no warning."""
        self.genome1.sequence = "ABCD"
        self.genome2.sequence = "ABCD"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_genome_sequence()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_genome_sequence_2(self):
        """Check that different sequences produce a warning."""
        self.genome1.sequence = "ABCD"
        self.genome2.sequence = "ABCDE"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_genome_sequence()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_genome_length_1(self):
        """Check that identical sequence lengths produce no warning."""
        self.genome1._length = 5
        self.genome2._length = 5
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_genome_length()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_genome_length_2(self):
        """Check that different sequence lengths produce a warning."""
        self.genome1._length = 5
        self.genome2._length = 6
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_genome_length()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_cluster_1(self):
        """Check that identical clusters produce no warning."""
        self.genome1.cluster = "A"
        self.genome2.cluster = "A"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_cluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_cluster_2(self):
        """Check that different clusters produce a warning."""
        self.genome1.cluster = "A"
        self.genome2.cluster = "B"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_cluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_subcluster_1(self):
        """Check that identical subclusters produce no warning."""
        self.genome1.subcluster = "A1"
        self.genome2.subcluster = "A1"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_subcluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_subcluster_2(self):
        """Check that different subclusters produce a warning."""
        self.genome1.subcluster = "A1"
        self.genome2.subcluster = "A2"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_subcluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_accession_1(self):
        """Check that identical accessions produce no warning."""
        self.genome1.accession = "ABC123"
        self.genome2.accession = "ABC123"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_accession()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_accession_2(self):
        """Check that different accessions produce a warning."""
        self.genome1.accession = "ABC1234"
        self.genome2.accession = "ABC123"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_accession()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_host_1(self):
        """Check that identical hosts produce no warning."""
        self.genome1.host = "Mycobacterium"
        self.genome2.host = "Mycobacterium"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_host()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_host_2(self):
        """Check that different hosts produce a warning."""
        self.genome1.host = "Mycobacterium"
        self.genome2.host = "Mycobacteriums"
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2
        self.genome_pair.compare_host()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




if __name__ == '__main__':
    unittest.main()
