""" Unit tests for the DataGroup Class."""


from classes import DataGroup
from classes import Genome
from classes import GenomePair
from classes import Cds
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




    def test_check_genome_dictionary_1(self):
        """Check that no error is produced when a genome is present
        in the dictionary and is expected to be present."""
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.check_genome_dictionary("flat_file")
        self.assertEqual(self.data_group.evaluations[0].status, "correct")

    def test_check_genome_dictionary_2(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.check_genome_dictionary("flat_file", False)
        self.assertEqual(self.data_group.evaluations[0].status, "error")

    def test_check_genome_dictionary_3(self):
        """Check that no error is produced when a genome is not present
        in the dictionary and is not expected to be present."""
        self.data_group.check_genome_dictionary("flat_file", False)
        self.assertEqual(self.data_group.evaluations[0].status, "correct")

    def test_check_genome_dictionary_4(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.data_group.check_genome_dictionary("flat_file")
        self.assertEqual(self.data_group.evaluations[0].status, "error")




    def test_check_genome_pair_dictionary_1(self):
        """Check that no error is produced when a genome_pair is present
        in the dictionary and is expected to be present."""
        self.data_group.genome_pair_dict["flat_file_phamerator"] = ""
        self.data_group.check_genome_pair_dictionary("flat_file_phamerator")
        self.assertEqual(self.data_group.evaluations[0].status, "correct")

    def test_check_genome_pair_dictionary_2(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.data_group.genome_pair_dict["flat_file_phamerator"] = ""
        self.data_group.check_genome_pair_dictionary("flat_file_phamerator", False)
        self.assertEqual(self.data_group.evaluations[0].status, "error")

    def test_check_genome_pair_dictionary_3(self):
        """Check that no error is produced when a genome_pair is not present
        in the dictionary and is not expected to be present."""
        self.data_group.check_genome_pair_dictionary("flat_file", False)
        self.assertEqual(self.data_group.evaluations[0].status, "correct")

    def test_check_genome_pair_dictionary_4(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.data_group.check_genome_pair_dictionary("flat_file")
        self.assertEqual(self.data_group.evaluations[0].status, "error")





class TestDataGroupClass2(unittest.TestCase):


    def setUp(self):

        self.ticket1 = Ticket.GenomeTicket()
        self.cds1 = Cds.CdsFeature()
        self.cds2 = Cds.CdsFeature()
        self.genome1 = Genome.Genome()
        self.genome1.type = "flat_file"
        self.genome1.cds_features.append(self.cds1)
        self.genome1.cds_features.append(self.cds2)
        self.genome2 = Genome.Genome()
        self.genome2.type = "phamerator"
        self.genome_pair1 = GenomePair.GenomePair()
        self.genome_pair2 = GenomePair.GenomePair()
        self.data_group = DataGroup.DataGroup()
        self.data_group.ticket = self.ticket1
        self.data_group.genome_dict[self.genome1.type] = self.genome1
        self.data_group.genome_dict[self.genome2.type] = self.genome2
        self.data_group.genome_pair_dict["genome_pair1"] = self.genome_pair1
        self.data_group.genome_pair_dict["genome_pair2"] = self.genome_pair2
        self.eval1 = Eval.Eval()
        self.eval2 = Eval.Eval()




    def test_check_for_errors_1(self):
        """Check that no error is counted."""

        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_2(self):
        """Check that a DataGroup 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.data_group.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_3(self):
        """Check that a DataGroup 'error' eval is counted."""
        self.eval1.status = "error"
        self.data_group.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 1)




    def test_check_for_errors_4(self):
        """Check that a ticket 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.ticket1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_5(self):
        """Check that a ticket 'error' eval is counted."""
        self.eval1.status = "error"
        self.ticket1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 1)




    def test_check_for_errors_6(self):
        """Check that a genome 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.genome1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_7(self):
        """Check that a genome 'error' eval is counted."""
        self.eval1.status = "error"
        self.genome1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 1)

    def test_check_for_errors_8(self):
        """Check that two genome 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.genome1.evaluations.append(self.eval1)
        self.genome2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_9(self):
        """Check that two genome 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.genome1.evaluations.append(self.eval1)
        self.genome2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 2)

    def test_check_for_errors_10(self):
        """Check that a cds 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.cds1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_11(self):
        """Check that a cds 'error' eval is counted."""
        self.eval1.status = "error"
        self.cds1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 1)

    def test_check_for_errors_12(self):
        """Check that two cds 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.cds1.evaluations.append(self.eval1)
        self.cds2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_13(self):
        """Check that two cds 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.cds1.evaluations.append(self.eval1)
        self.cds2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 2)

    def test_check_for_errors_14(self):
        """Check that a genome_pair 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.genome_pair1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_15(self):
        """Check that a genome_pair 'error' eval is counted."""
        self.eval1.status = "error"
        self.genome_pair1.evaluations.append(self.eval1)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 1)

    def test_check_for_errors_16(self):
        """Check that two genome_pair 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.genome_pair1.evaluations.append(self.eval1)
        self.genome_pair2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 0)

    def test_check_for_errors_17(self):
        """Check that two genome_pair 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.genome_pair1.evaluations.append(self.eval1)
        self.genome_pair2.evaluations.append(self.eval2)
        self.data_group.check_for_errors()
        self.assertEqual(self.data_group._errors, 2)


if __name__ == '__main__':
    unittest.main()
