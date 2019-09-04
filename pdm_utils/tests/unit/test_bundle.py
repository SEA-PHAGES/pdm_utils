""" Unit tests for the Bundle Class."""


from classes import bundle
from classes import Genome
from classes import GenomePair
from classes import cds
from classes import Ticket
from classes import eval
import unittest


class TestBundleClass1(unittest.TestCase):


    def setUp(self):

        self.bndl = bundle.Bundle()
        self.genome1 = Genome.Genome()
        self.genome1.type = "flat_file"
        self.genome2 = Genome.Genome()
        self.genome2.type = "phamerator"
        self.ticket = Ticket.GenomeTicket()




    def test_set_genome_pair_1(self):
        """Check that a genome pair is set if both keys are present."""

        self.bndl.ticket = self.ticket
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        genome_pair = GenomePair.GenomePair()
        self.bndl.set_genome_pair(genome_pair, "phamerator", "flat_file")
        self.assertEqual(list(self.bndl.genome_pair_dict.keys())[0],
                            "phamerator_flat_file")

    def test_set_genome_pair_2(self):
        """Check that a genome pair is not set if one key is not present."""

        self.bndl.ticket = self.ticket
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        genome_pair = GenomePair.GenomePair()
        self.bndl.set_genome_pair(genome_pair, "invalid", "flat_file")
        self.assertEqual(len(self.bndl.genome_pair_dict.keys()), 0)




    def test_check_matched_genome_1(self):
        """Check that no error is produced when the genome type is present."""
        self.bndl.ticket = self.ticket
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        self.bndl.check_matched_genome("phamerator", "eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_matched_genome_2(self):
        """Check that an error is produced when the genome type is present."""
        self.bndl.ticket = self.ticket
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        self.bndl.check_matched_genome("invalid")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)




    def test_check_genome_dictionary_1(self):
        """Check that no error is produced when a genome is present
        in the dictionary and is expected to be present."""
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.check_genome_dictionary("flat_file", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_genome_dictionary_2(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.check_genome_dictionary("flat_file", False)
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)

    def test_check_genome_dictionary_3(self):
        """Check that no error is produced when a genome is not present
        in the dictionary and is not expected to be present."""
        self.bndl.check_genome_dictionary("flat_file", False)
        self.assertEqual(self.bndl.evaluations[0].status, "correct")

    def test_check_genome_dictionary_4(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.bndl.check_genome_dictionary("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "error")




    def test_check_genome_pair_dictionary_1(self):
        """Check that no error is produced when a genome_pair is present
        in the dictionary and is expected to be present."""
        self.bndl.genome_pair_dict["flat_file_phamerator"] = ""
        self.bndl.check_genome_pair_dictionary(
            "flat_file_phamerator", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_genome_pair_dictionary_2(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.bndl.genome_pair_dict["flat_file_phamerator"] = ""
        self.bndl.check_genome_pair_dictionary("flat_file_phamerator", False)
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)

    def test_check_genome_pair_dictionary_3(self):
        """Check that no error is produced when a genome_pair is not present
        in the dictionary and is not expected to be present."""
        self.bndl.check_genome_pair_dictionary("flat_file", False)
        self.assertEqual(self.bndl.evaluations[0].status, "correct")

    def test_check_genome_pair_dictionary_4(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.bndl.check_genome_pair_dictionary("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "error")





class TestBundleClass2(unittest.TestCase):


    def setUp(self):

        self.ticket1 = Ticket.GenomeTicket()
        self.cds1 = cds.Cds()
        self.cds1.id = "L5_1"
        self.cds2 = cds.Cds()
        self.cds2.id = "L5_2"
        self.cds3 = cds.Cds()
        self.cds3.id = "L5_3"
        self.genome1 = Genome.Genome()
        self.genome1.type = "flat_file"
        self.genome1.cds_features.append(self.cds1)
        self.genome1.cds_features.append(self.cds2)
        self.genome2 = Genome.Genome()
        self.genome2.type = "phamerator"
        self.genome_pair1 = GenomePair.GenomePair()
        self.genome_pair2 = GenomePair.GenomePair()
        self.bndl = bundle.Bundle()
        self.bndl.ticket = self.ticket1
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        self.bndl.genome_pair_dict["genome_pair1"] = self.genome_pair1
        self.bndl.genome_pair_dict["genome_pair2"] = self.genome_pair2
        self.eval1 = eval.Eval()
        self.eval2 = eval.Eval()




    def test_check_for_errors_1(self):
        """Check that no error is counted."""
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_2(self):
        """Check that a Bundle 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.bndl.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_3(self):
        """Check that a Bundle 'error' eval is counted."""
        self.eval1.status = "error"
        self.bndl.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_4(self):
        """Check that a ticket 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.ticket1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_5(self):
        """Check that a ticket 'error' eval is counted."""
        self.eval1.status = "error"
        self.ticket1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_6(self):
        """Check that a genome 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.genome1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_7(self):
        """Check that a genome 'error' eval is counted."""
        self.eval1.status = "error"
        self.genome1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_8(self):
        """Check that two genome 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.genome1.evaluations.append(self.eval1)
        self.genome2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_9(self):
        """Check that two genome 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.genome1.evaluations.append(self.eval1)
        self.genome2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)

    def test_check_for_errors_10(self):
        """Check that a cds 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.cds1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_11(self):
        """Check that a cds 'error' eval is counted."""
        self.eval1.status = "error"
        self.cds1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_12(self):
        """Check that two cds 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.cds1.evaluations.append(self.eval1)
        self.cds2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_13(self):
        """Check that two cds 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.cds1.evaluations.append(self.eval1)
        self.cds2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)

    def test_check_for_errors_14(self):
        """Check that a genome_pair 'correct' eval is not counted."""
        self.eval1.status = "correct"
        self.genome_pair1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_15(self):
        """Check that a genome_pair 'error' eval is counted."""
        self.eval1.status = "error"
        self.genome_pair1.evaluations.append(self.eval1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_16(self):
        """Check that two genome_pair 'correct' evals are not counted."""
        self.eval1.status = "correct"
        self.eval2.status = "correct"
        self.genome_pair1.evaluations.append(self.eval1)
        self.genome_pair2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_17(self):
        """Check that two genome_pair 'error' evals are counted."""
        self.eval1.status = "error"
        self.eval2.status = "error"
        self.genome_pair1.evaluations.append(self.eval1)
        self.genome_pair2.evaluations.append(self.eval2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)




    def test_get_evaluations_1(self):
        """Verify one eval is returned from Bundle eval list."""
        self.bndl.evaluations.append(self.eval1)
        eval_dict = self.bndl.get_evaluations()
        self.assertEqual(len(eval_dict["bundle"]), 1)

    def test_get_evaluations_2(self):
        """Verify one eval is returned from Ticket eval list."""
        self.ticket1.evaluations.append(self.eval1)
        eval_dict = self.bndl.get_evaluations()
        self.assertEqual(len(eval_dict["ticket"]), 1)

    def test_get_evaluations_3(self):
        """Verify one eval is returned from each genome eval list."""
        self.genome1.evaluations.append(self.eval1)
        self.genome1.evaluations.append(self.eval2)
        self.genome2.evaluations.append(self.eval1)
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_flat_file"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_phamerator"]), 1)

    def test_get_evaluations_4(self):
        """Verify one eval is returned from each Cds eval list in
        each genome."""
        self.cds1.evaluations.append(self.eval1)
        self.cds1.evaluations.append(self.eval2)
        self.cds2.evaluations.append(self.eval1)
        self.cds3.evaluations.append(self.eval1)
        self.genome1.cds_features = [self.cds1, self.cds2]
        self.genome2.cds_features = [self.cds3]
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_1"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_2"]), 1)
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_3"]), 1)

    def test_get_evaluations_5(self):
        """Verify one eval is returned from each genome_pair eval list."""
        self.genome_pair1.evaluations.append(self.eval1)
        self.genome_pair1.evaluations.append(self.eval2)
        self.genome_pair2.evaluations.append(self.eval1)
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_pair_genome_pair1"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_pair_genome_pair2"]), 1)




if __name__ == '__main__':
    unittest.main()
