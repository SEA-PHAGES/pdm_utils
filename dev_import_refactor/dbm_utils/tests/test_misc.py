""" Unit tests for misc. functions."""

from classes import DataGroup
from classes import Genome
from functions import misc
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestMiscFunctions(unittest.TestCase):


    def setUp(self):
        pass

        self.genome1 = Genome.Genome()
        self.genome1.phage_id = "Trixie"
        self.genome1.phage_name = "Genome1"
        self.genome1.type = "phamerator"

        self.genome2 = Genome.Genome()
        self.genome2.phage_id = "Trixie"
        self.genome2.phage_name = "Genome2"
        self.genome2.type = "flat_file"


        self.genome3 = Genome.Genome()
        self.genome3.phage_id = "L5"
        self.genome3.phage_name = "Genome3"
        self.genome3.type = "phamerator"

        self.genome4 = Genome.Genome()
        self.genome4.phage_id = "L5"
        self.genome4.phage_name = "Genome4"
        self.genome4.type = "flat_file"

        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup2 = DataGroup.DataGroup()








    def test_match_genome_by_phage_id_1(self):
        """Verify that one genome is matched correctly."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genome_by_phage_id(
            self.datagroup1, genomes_to_match, "phamerator")
        matched_genome = self.datagroup1.genome_dict["flat_file"]
        self.assertEqual(matched_genome.phage_name, "Genome2")

    def test_match_genome_by_phage_id_2(self):
        """Verify that no genome is matched since there is no
        reference genome in the DataGroup that matches the key."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genome_by_phage_id(
            self.datagroup1, genomes_to_match, "invalid")
        self.assertEqual(len(self.datagroup1.genome_dict.keys()), 1)

    def test_match_genome_by_phage_id_3(self):
        """Verify that no genome is matched since there is no
        genome in the genome dictionary for matching."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        genomes_to_match = {}
        misc.match_genome_by_phage_id(
            self.datagroup1, genomes_to_match, "phamerator")
        self.assertEqual(len(self.datagroup1.genome_dict.keys()), 1)

    def test_match_genome_by_phage_id_4(self):
        """Verify that one genome is matched correctly with the 'key2'
        parameter explicitly added."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genome_by_phage_id(
            self.datagroup1, genomes_to_match, "phamerator", "new_type")
        matched_genome = self.datagroup1.genome_dict["new_type"]
        self.assertEqual(matched_genome.phage_name, "Genome2")










    def test_match_genomes_1(self):
        """Verify that one genome is matched correctly."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        list1 = [self.datagroup1] # Trixie phamerator genome.
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genomes(list1, genomes_to_match, "phamerator", "new_genome")
        matched_genome = list1[0].genome_dict["new_genome"]
        self.assertEqual(matched_genome.phage_name, "Genome2")


    def test_match_genomes_2(self):
        """Verify that no genome is matched since there is no
        DataGroup object."""

        list1 = []
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genomes(list1, genomes_to_match, "phamerator", "new_genome")
        self.assertEqual(len(list1), 0)

    def test_match_genomes_3(self):
        """Verify that no genome is matched since there is no
        reference genome in the DataGroup."""

        list1 = [self.datagroup1]
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genomes(list1, genomes_to_match, "phamerator", "new_genome")
        self.assertEqual(len(list1[0].genome_dict.keys()), 0)

    def test_match_genomes_4(self):
        """Verify that no genome is matched since there is no
        genome in the genome dictionary for matching."""

        list1 = [self.datagroup1]
        genomes_to_match = {}
        misc.match_genomes(list1, genomes_to_match, "phamerator", "new_genome")
        self.assertEqual(len(list1[0].genome_dict.keys()), 0)

    def test_match_genomes_5(self):
        """Verify that two genomes are matched correctly."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        self.datagroup2.genome_dict[self.genome3.type] = self.genome3

        # Trixie phamerator genome, L5 phamerator genome.
        list1 = [self.datagroup1, self.datagroup2]

        genomes_to_match = {self.genome2.phage_id: self.genome2,
                            self.genome4.phage_id: self.genome4}
        misc.match_genomes(list1, genomes_to_match, "phamerator", "new_genome")
        matched_genome2 = list1[0].genome_dict["new_genome"]
        matched_genome4 = list1[1].genome_dict["new_genome"]
        with self.subTest():
            self.assertEqual(matched_genome2.phage_name, "Genome2")
        with self.subTest():
            self.assertEqual(matched_genome4.phage_name, "Genome4")

    def test_match_genomes_6(self):
        """Verify that one genome is matched correctly with the 'key2'
        parameter omitted."""

        self.datagroup1.genome_dict[self.genome1.type] = self.genome1
        list1 = [self.datagroup1] # Trixie phamerator genome.
        genomes_to_match = {self.genome2.phage_id: self.genome2}
        misc.match_genomes(list1, genomes_to_match, "phamerator")
        matched_genome = list1[0].genome_dict["flat_file"]
        self.assertEqual(matched_genome.phage_name, "Genome2")





    def test_create_fasta_seqrecord_1(self):
        """."""
        record = misc.create_fasta_seqrecord("Mycobacterium phage Trixie", "ATCGC")
        with self.subTest():
            self.assertTrue(isinstance(record, SeqRecord))
        with self.subTest():
            self.assertTrue(isinstance(record.seq, Seq))
        with self.subTest():
            self.assertEqual(record.description, "Mycobacterium phage Trixie")
        with self.subTest():
            self.assertEqual(record.seq, "ATCGC")



















if __name__ == '__main__':
    unittest.main()
