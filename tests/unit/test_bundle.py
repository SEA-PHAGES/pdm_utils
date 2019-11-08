""" Unit tests for the Bundle Class."""


from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import genomepair
from pdm_utils.classes import source
from pdm_utils.classes import cds
from pdm_utils.classes import ticket
from pdm_utils.classes import eval
import unittest


class TestBundleClass1(unittest.TestCase):


    def setUp(self):

        self.bndl = bundle.Bundle()
        self.genome1 = genome.Genome()
        self.genome1.type = "flat_file"
        self.genome2 = genome.Genome()
        self.genome2.type = "phamerator"
        self.tkt = ticket.GenomeTicket()




    def test_set_genome_pair_1(self):
        """Check that a genome pair is set if both keys are present."""

        self.bndl.ticket = self.tkt
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        genome_pair = genomepair.GenomePair()
        self.bndl.set_genome_pair(genome_pair, "phamerator", "flat_file")
        self.assertEqual(list(self.bndl.genome_pair_dict.keys())[0],
                            "phamerator_flat_file")

    def test_set_genome_pair_2(self):
        """Check that a genome pair is not set if one key is not present."""

        self.bndl.ticket = self.tkt
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        genome_pair = genomepair.GenomePair()
        self.bndl.set_genome_pair(genome_pair, "invalid", "flat_file")
        self.assertEqual(len(self.bndl.genome_pair_dict.keys()), 0)




    def test_check_ticket_1(self):
        """Check that no error is produced when a ticket is present."""
        self.bndl.ticket = self.tkt
        self.bndl.check_ticket("eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_ticket_2(self):
        """Check that an error is produced when a ticket is not present."""
        self.bndl.ticket = None
        self.bndl.check_ticket("eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")








    def test_check_genome_dict_1(self):
        """Check that no error is produced when a genome is present
        in the dictionary and is expected to be present."""
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.check_genome_dict("flat_file", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_genome_dict_2(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.check_genome_dict("flat_file", False)
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)

    def test_check_genome_dict_3(self):
        """Check that no error is produced when a genome is not present
        in the dictionary and is not expected to be present."""
        self.bndl.check_genome_dict("flat_file", False)
        self.assertEqual(self.bndl.evaluations[0].status, "correct")

    def test_check_genome_dict_4(self):
        """Check that an error is produced when a genome is not present
        in the dictionary and is expected to be present."""
        self.bndl.check_genome_dict("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "error")




    def test_check_genome_pair_dict_1(self):
        """Check that no error is produced when a genome_pair is present
        in the dictionary and is expected to be present."""
        self.bndl.genome_pair_dict["flat_file_phamerator"] = ""
        self.bndl.check_genome_pair_dict(
            "flat_file_phamerator", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_genome_pair_dict_2(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.bndl.genome_pair_dict["flat_file_phamerator"] = ""
        self.bndl.check_genome_pair_dict("flat_file_phamerator", False)
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)

    def test_check_genome_pair_dict_3(self):
        """Check that no error is produced when a genome_pair is not present
        in the dictionary and is not expected to be present."""
        self.bndl.check_genome_pair_dict("flat_file", False)
        self.assertEqual(self.bndl.evaluations[0].status, "correct")

    def test_check_genome_pair_dict_4(self):
        """Check that an error is produced when a genome_pair is not present
        in the dictionary and is expected to be present."""
        self.bndl.check_genome_pair_dict("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "error")





class TestBundleClass2(unittest.TestCase):


    def setUp(self):

        self.ticket1 = ticket.GenomeTicket()
        self.src1 = source.Source()
        self.src1.id = "L5_SRC_1"
        self.src2 = source.Source()
        self.src2.id = "L5_SRC_2"
        self.src3 = source.Source()
        self.src3.id = "L5_SRC_3"
        self.cds1 = cds.Cds()
        self.cds1.id = "L5_CDS_1"
        self.cds2 = cds.Cds()
        self.cds2.id = "L5_CDS_2"
        self.cds3 = cds.Cds()
        self.cds3.id = "L5_CDS_3"
        self.genome1 = genome.Genome()
        self.genome1.type = "flat_file"
        self.genome1.cds_features.append(self.cds1)
        self.genome1.cds_features.append(self.cds2)
        self.genome1.source_features.append(self.src1)
        self.genome1.source_features.append(self.src2)
        self.genome2 = genome.Genome()
        self.genome2.type = "phamerator"
        self.genome_pair1 = genomepair.GenomePair()
        self.genome_pair2 = genomepair.GenomePair()
        self.bndl = bundle.Bundle()
        self.bndl.ticket = self.ticket1
        self.bndl.genome_dict[self.genome1.type] = self.genome1
        self.bndl.genome_dict[self.genome2.type] = self.genome2
        self.bndl.genome_pair_dict["genome_pair1"] = self.genome_pair1
        self.bndl.genome_pair_dict["genome_pair2"] = self.genome_pair2

        self.eval_correct1 = eval.Eval(status="correct")
        self.eval_correct2 = eval.Eval(status="correct")
        self.eval_error1 = eval.Eval(status="error")
        self.eval_error2 = eval.Eval(status="error")




    def test_check_for_errors_1(self):
        """Check that no error is counted."""
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_2(self):
        """Check that a Bundle 'correct' eval is not counted."""
        self.bndl.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_3(self):
        """Check that a Bundle 'error' eval is counted."""
        self.bndl.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_4(self):
        """Check that a ticket 'correct' eval is not counted."""
        self.ticket1.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_5(self):
        """Check that a ticket 'error' eval is counted."""
        self.ticket1.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_6(self):
        """Check that a bundle with no ticket is not counted."""
        self.bndl.ticket = None
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_7(self):
        """Check that a genome 'correct' eval is not counted."""
        self.genome1.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_8(self):
        """Check that a genome 'error' eval is counted."""
        self.genome1.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_9(self):
        """Check that two genome 'correct' evals are not counted."""
        self.genome1.evaluations.append(self.eval_correct1)
        self.genome2.evaluations.append(self.eval_correct2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_10(self):
        """Check that two genome 'error' evals are counted."""
        self.genome1.evaluations.append(self.eval_error1)
        self.genome2.evaluations.append(self.eval_error2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)

    def test_check_for_errors_11(self):
        """Check that a source feature 'correct' eval is not counted."""
        self.src1.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_12(self):
        """Check that a source feature 'error' eval is counted."""
        self.src1.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_13(self):
        """Check that two source feature 'correct' evals are not counted."""
        self.src1.evaluations.append(self.eval_correct1)
        self.src2.evaluations.append(self.eval_correct2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_14(self):
        """Check that two source feature 'error' evals are counted."""
        self.src1.evaluations.append(self.eval_error1)
        self.src2.evaluations.append(self.eval_error2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)

    def test_check_for_errors_15(self):
        """Check that a cds 'correct' eval is not counted."""
        self.cds1.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_16(self):
        """Check that a cds 'error' eval is counted."""
        self.cds1.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_17(self):
        """Check that two cds 'correct' evals are not counted."""
        self.cds1.evaluations.append(self.eval_correct1)
        self.cds2.evaluations.append(self.eval_correct2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_18(self):
        """Check that two cds 'error' evals are counted."""
        self.cds1.evaluations.append(self.eval_error1)
        self.cds2.evaluations.append(self.eval_error2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)

    def test_check_for_errors_19(self):
        """Check that a genome_pair 'correct' eval is not counted."""
        self.genome_pair1.evaluations.append(self.eval_correct1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_20(self):
        """Check that a genome_pair 'error' eval is counted."""
        self.genome_pair1.evaluations.append(self.eval_error1)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 1)

    def test_check_for_errors_21(self):
        """Check that two genome_pair 'correct' evals are not counted."""
        self.genome_pair1.evaluations.append(self.eval_correct1)
        self.genome_pair2.evaluations.append(self.eval_correct2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 0)

    def test_check_for_errors_22(self):
        """Check that two genome_pair 'error' evals are counted."""
        self.genome_pair1.evaluations.append(self.eval_error1)
        self.genome_pair2.evaluations.append(self.eval_error2)
        self.bndl.check_for_errors()
        self.assertEqual(self.bndl._errors, 2)




    def test_get_evaluations_1(self):
        """Verify one eval is returned from Bundle eval list."""
        self.bndl.evaluations.append(self.eval_correct1)
        eval_dict = self.bndl.get_evaluations()
        self.assertEqual(len(eval_dict["bundle"]), 1)

    def test_get_evaluations_2(self):
        """Verify one eval is returned from Ticket eval list."""
        self.ticket1.evaluations.append(self.eval_correct1)
        eval_dict = self.bndl.get_evaluations()
        self.assertEqual(len(eval_dict["ticket"]), 1)

    def test_get_evaluations_3(self):
        """Verify one eval is returned from each genome eval list."""
        self.genome1.evaluations.append(self.eval_correct1)
        self.genome1.evaluations.append(self.eval_correct2)
        self.genome2.evaluations.append(self.eval_error1)
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_flat_file"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_phamerator"]), 1)

    def test_get_evaluations_4(self):
        """Verify one eval is returned from each Source eval list in
        each genome."""
        self.src1.evaluations.append(self.eval_correct1)
        self.src1.evaluations.append(self.eval_correct2)
        self.src2.evaluations.append(self.eval_error1)
        self.src3.evaluations.append(self.eval_error2)
        self.genome1.source_features = [self.src1, self.src2]
        self.genome2.source_features = [self.src3]
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["src_L5_SRC_1"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["src_L5_SRC_2"]), 1)
        with self.subTest():
            self.assertEqual(len(eval_dict["src_L5_SRC_3"]), 1)

    def test_get_evaluations_5(self):
        """Verify one eval is returned from each Cds eval list in
        each genome."""
        self.cds1.evaluations.append(self.eval_correct1)
        self.cds1.evaluations.append(self.eval_correct2)
        self.cds2.evaluations.append(self.eval_error1)
        self.cds3.evaluations.append(self.eval_error2)
        self.genome1.cds_features = [self.cds1, self.cds2]
        self.genome2.cds_features = [self.cds3]
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_CDS_1"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_CDS_2"]), 1)
        with self.subTest():
            self.assertEqual(len(eval_dict["cds_L5_CDS_3"]), 1)

    def test_get_evaluations_6(self):
        """Verify one eval is returned from each genome_pair eval list."""
        self.genome_pair1.evaluations.append(self.eval_correct1)
        self.genome_pair1.evaluations.append(self.eval_correct2)
        self.genome_pair2.evaluations.append(self.eval_error1)
        eval_dict = self.bndl.get_evaluations()
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_pair_genome_pair1"]), 2)
        with self.subTest():
            self.assertEqual(len(eval_dict["genome_pair_genome_pair2"]), 1)




class TestBundleClass3(unittest.TestCase):


    def setUp(self):
        self.bndl = bundle.Bundle()
        self.tkt = ticket.GenomeTicket()
        self.gnm = genome.Genome()
        self.bndl.genome_dict["flat_file"] = self.gnm
        self.bndl.ticket = self.tkt




    def test_check_compatible_type_and_annotation_status_1(self):
        """Verify that no error is produced with "add" type and "draft"
        annotation_status."""
        self.tkt.type = "add"
        self.gnm.annotation_status = "draft"
        self.bndl.check_compatible_type_and_annotation_status("flat_file",
            eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].id, "eval_id")

    def test_check_compatible_type_and_annotation_status_2(self):
        """Verify that an error is produced with "add" type and "final"
        annotation_status."""
        self.tkt.type = "add"
        self.gnm.annotation_status = "final"
        self.bndl.check_compatible_type_and_annotation_status("flat_file")
        with self.subTest():
            self.assertEqual(self.bndl.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.bndl.evaluations[0].id)

    def test_check_compatible_type_and_annotation_status_3(self):
        """Verify that no error is produced with "replace" type and "final"
        annotation_status."""
        self.tkt.type = "replace"
        self.gnm.annotation_status = "final"
        self.bndl.check_compatible_type_and_annotation_status("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "correct")

    def test_check_compatible_type_and_annotation_status_4(self):
        """Verify that an error is produced with "replace" type and "draft"
        annotation_status."""
        self.tkt.type = "replace"
        self.gnm.annotation_status = "draft"
        self.bndl.check_compatible_type_and_annotation_status("flat_file")
        self.assertEqual(self.bndl.evaluations[0].status, "error")





if __name__ == '__main__':
    unittest.main()
