""" Unit tests for the ticket class"""

from pdm_utils.classes import ticket
import unittest

class TestGenomeTicketClass(unittest.TestCase):

    def setUp(self):

        # Empty ticket to test simple methods
        self.tkt = ticket.GenomeTicket()




    def test_set_type_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "ADD"
        self.tkt.set_type(value)
        self.assertEqual(self.tkt.type, "add")

    def test_set_type_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_type(value)
        self.assertEqual(self.tkt.type, "none")




    def test_set_description_field_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Product"
        self.tkt.set_description_field(value)
        self.assertEqual(self.tkt.description_field, "product")

    def test_set_description_field_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_description_field(value)
        self.assertEqual(self.tkt.description_field, "none")




    def test_set_run_mode_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "PhagesDB"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "phagesdb")

    def test_set_run_mode_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "none")




    def test_set_value_flag_1(self):
        """Verify that the 'retain' setting is set to True."""
        self.tkt._value_flag = False
        self.tkt.cluster = "retain"
        self.tkt.set_value_flag("retain")
        self.assertTrue(self.tkt._value_flag)

    def test_set_value_flag_2(self):
        """Verify that the 'retain' setting is set to False."""
        self.tkt._value_flag = True
        self.tkt.cluster = "A"
        self.tkt.set_value_flag("retain")
        self.assertFalse(self.tkt._value_flag)




    def test_set_field_trackers_1(self):
        """Check that keys are assigned to sets as expected."""
        self.tkt.data_dict = {"host_genus":"retrieve",
                              "cluster":"retain",
                              "subcluster":"A2"}
        self.tkt.set_field_trackers()
        with self.subTest():
            self.assertEqual(self.tkt.data_retrieve, {"host_genus"})
        with self.subTest():
            self.assertEqual(self.tkt.data_retain, {"cluster"})
        with self.subTest():
            self.assertEqual(self.tkt.data_ticket, {"subcluster"})

    def test_set_field_trackers_2(self):
        """Check that no keys are assigned to sets."""
        self.tkt.data_dict = {}
        self.tkt.set_field_trackers()
        with self.subTest():
            self.assertEqual(len(self.tkt.data_retrieve), 0)
        with self.subTest():
            self.assertEqual(len(self.tkt.data_retain), 0)
        with self.subTest():
            self.assertEqual(len(self.tkt.data_ticket), 0)


    def test_check_type_1(self):
        """Check that no error is produced if the
        type is present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.tkt.type = "add"
        self.tkt.check_type(set1, True, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_type_2(self):
        """Check that an error is produced if the
        type is not present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.tkt.type = "none"
        self.tkt.check_type(set1)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.phage_id = "Trixie"
        self.tkt.check_phage_id(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")



    def test_check_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the empty/null set."""
        set1 = set(["none"])
        self.tkt.phage_id = "none"
        self.tkt.check_phage_id(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_description_field_1(self):
        """Check that no error is produced if the
        description_field is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.description_field = "Mycobacterium"
        self.tkt.check_description_field(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_description_field_2(self):
        """Check that an error is produced if the
        description_field is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.description_field = "none"
        self.tkt.check_description_field(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_run_mode_1(self):
        """Check that no error is produced if the
        run_mode is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.run_mode = "Mycobacterium"
        self.tkt.check_run_mode(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_run_mode_2(self):
        """Check that an error is produced if the
        run_mode is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.run_mode = "none"
        self.tkt.check_run_mode(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_duplicate_id_1(self):
        """Check that no error is produced if the
        id is not present in the set of duplicated values."""
        dupe_set = set([0, 1])
        self.tkt.phage_id = 2
        self.tkt.check_duplicate_id(dupe_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_duplicate_id_2(self):
        """Check that an error is produced if the
        id is present in the set of duplicated values."""
        dupe_set = set([0, 1])
        self.tkt.id = 1
        self.tkt.check_duplicate_id(dupe_set)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_duplicate_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.tkt.phage_id = "D29"
        self.tkt.check_duplicate_phage_id(dupe_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_duplicate_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.tkt.phage_id = "Trixie"
        self.tkt.check_duplicate_phage_id(dupe_set)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)







    def test_check_compatible_type_and_annotation_status_1(self):
        """Verify that no error is produced with "add" type and "draft"
        annotation_status."""
        self.tkt.type = "add"
        self.tkt.annotation_status = "draft"
        self.tkt.check_compatible_type_and_annotation_status(
            eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_compatible_type_and_annotation_status_2(self):
        """Verify that an error is produced with "add" type and "final"
        annotation_status."""
        self.tkt.type = "add"
        self.tkt.annotation_status = "final"
        self.tkt.check_compatible_type_and_annotation_status()
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)

    def test_check_compatible_type_and_annotation_status_3(self):
        """Verify that no error is produced with "replace" type and "final"
        annotation_status."""
        self.tkt.type = "replace"
        self.tkt.annotation_status = "final"
        self.tkt.check_compatible_type_and_annotation_status()
        self.assertEqual(self.tkt.evaluations[0].status, "correct")

    def test_check_compatible_type_and_annotation_status_4(self):
        """Verify that no error is produced with "replace" type and "draft"
        annotation_status."""
        self.tkt.type = "replace"
        self.tkt.annotation_status = "draft"
        self.tkt.check_compatible_type_and_annotation_status()
        self.assertEqual(self.tkt.evaluations[0].status, "error")






if __name__ == '__main__':
    unittest.main()
