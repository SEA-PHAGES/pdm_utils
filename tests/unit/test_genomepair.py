""" Unit tests for the GenomePair Class."""

from datetime import datetime
from pdm_utils.classes import genomepair
from pdm_utils.classes import genome
from pdm_utils.classes import ticket
from pdm_utils.classes import eval
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class TestGenomePairClass(unittest.TestCase):


    def setUp(self):
        self.genome1 = genome.Genome()
        self.genome2 = genome.Genome()
        self.tkt = ticket.ImportTicket()
        self.genome_pair = genomepair.GenomePair()
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2

        self.date_jan1 = datetime.strptime('1/1/2000', '%m/%d/%Y')
        self.date_feb1 = datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.date_feb1_b = datetime.strptime('2/1/2000', '%m/%d/%Y')




    def test_compare_attribute_1(self):
        """Verify no error is produced when both genomes have
        identical ids as expected."""
        self.genome1.id = "Trixie"
        self.genome2.id = "Trixie"
        self.genome_pair.compare_attribute("id", expect_same=True,
                                           eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")


    def test_compare_attribute_2(self):
        """Verify no error is produced when both genomes have
        different ids as expected."""
        self.genome1.id = "Trixie"
        self.genome2.id = "L5"
        self.genome_pair.compare_attribute("id", expect_same=False)
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.genome_pair.evaluations[0].id)


    def test_compare_attribute_3(self):
        """Verify an error is produced when both genomes have
        different ids unexpectedly."""
        self.genome1.id = "Trixie"
        self.genome2.id = "L5"
        self.genome_pair.compare_attribute("id", expect_same=True)
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_attribute_4(self):
        """Verify an error is produced when both genomes have
        identical ids unexpectedly."""
        self.genome1.id = "Trixie"
        self.genome2.id = "Trixie"
        self.genome_pair.compare_attribute("id", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_attribute_5(self):
        """Verify no error is produced when both genomes have
        different retrieve_records as expected."""
        self.genome1.retrieve_record = 1
        self.genome2.retrieve_record = 0
        self.genome_pair.compare_attribute("retrieve_record", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")


    def test_compare_attribute_6(self):
        """Verify no error is produced when both genomes have
        identical seqs as expected."""
        self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome2.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome_pair.compare_attribute("seq", expect_same=True)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_attribute_7(self):
        """Verify no error is produced when both genomes have
        different seqs as expected."""
        self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome2.seq = Seq("AAAT", IUPAC.ambiguous_dna)
        self.genome_pair.compare_attribute("seq", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_attribute_8(self):
        """Verify no test is performed when the attribute is invalid."""
        self.genome_pair.compare_attribute("invalid", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")




    def test_compare_date_1(self):
        """Verify no error is produced when
        genome1 is older than genome2 as expected."""
        self.genome1.date = self.date_jan1
        self.genome2.date = self.date_feb1
        self.genome_pair.compare_date("older", "eval_id")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")

    def test_compare_date_2(self):
        """Verify no error is produced when
        genome1 is newer than genome2 as expected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_jan1
        self.genome_pair.compare_date("newer")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.genome_pair.evaluations[0].id)

    def test_compare_date_3(self):
        """Verify no error is produced when
        genome1 is equal to genome2 as expected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_feb1_b
        self.genome_pair.compare_date("equal")
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_date_4(self):
        """Verify no error is produced when
        genome1 is older than genome2 unexpectedly."""
        self.genome1.date = self.date_jan1
        self.genome2.date = self.date_feb1
        self.genome_pair.compare_date("newer")
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_date_5(self):
        """Verify no test is performed when an
        invalid comparison is selected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_jan1
        self.genome_pair.compare_date("invalid")
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")






if __name__ == '__main__':
    unittest.main()
