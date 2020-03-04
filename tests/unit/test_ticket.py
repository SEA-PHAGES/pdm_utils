""" Unit tests for the ticket class"""

from pdm_utils.classes import ticket
import unittest

class TestImportTicketClass(unittest.TestCase):

    def setUp(self):

        # Empty ticket to test simple methods
        self.tkt = ticket.ImportTicket()




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
        value = "Final"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "final")

    def test_set_run_mode_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "none")




    def test_check_attribute_1(self):
        """Verify no error is produced when the id
        is in the check_set and is expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.tkt.id = "Trixie"
        self.tkt.check_attribute("id", check_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_attribute_2(self):
        """Verify an error is produced when the id
        is not in the check_set and is expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.tkt.id = "D29"
        self.tkt.check_attribute("id", check_set, True)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)

    def test_check_attribute_3(self):
        """Verify no error is produced when the id
        is not in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.tkt.id = "D29"
        self.tkt.check_attribute("id", check_set, False)
        self.assertEqual(self.tkt.evaluations[0].status, "correct")

    def test_check_attribute_4(self):
        """Verify an error is produced when the id
        is in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.tkt.id = "Trixie"
        self.tkt.check_attribute("id", check_set, False)
        self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_attribute_5(self):
        """Verify no test is performed when the attribute is invalid."""
        check_set = set([1, 0])
        self.tkt.check_attribute("invalid", check_set, True)
        self.assertEqual(self.tkt.evaluations[0].status, "untested")




    def test_check_eval_flags_1(self):
        """Check that no error is produced if the
        eval_flags is empty and not expected to contain data."""
        self.tkt.eval_flag_dict = {}
        self.tkt.check_eval_flags(False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_eval_flags_2(self):
        """Check that no error is produced if the
        eval_flags is not empty and expected to contain data."""
        self.tkt.eval_flags = {"a":1, "b":2}
        self.tkt.check_eval_flags(True)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)

    def test_check_eval_flags_3(self):
        """Check that an error is produced if the
        eval_flags is not empty and not expected to contain data."""
        self.tkt.eval_flags = {"a":1, "b":2}
        self.tkt.check_eval_flags(False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_eval_flags_4(self):
        """Check that an error is produced if the
        eval_flags is empty and expected to contain data."""
        self.tkt.eval_flags = {}
        self.tkt.check_eval_flags(True)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")




    def test_check_compatible_type_and_data_retain_1(self):
        """Verify that no error is produced with "add" type and data_retain
        is empty."""
        self.tkt.type = "add"
        self.tkt.data_retain = set()
        self.tkt.check_compatible_type_and_data_retain(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_compatible_type_and_data_retain_2(self):
        """Verify that no error is produced with "replace" type and data_retain
        is not empty."""
        self.tkt.type = "replace"
        self.tkt.data_retain = set(["host_genus"])
        self.tkt.check_compatible_type_and_data_retain(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")

    def test_check_compatible_type_and_data_retain_3(self):
        """Verify that an error is produced with "add" type and data_retain
        is not empty."""
        self.tkt.type = "add"
        self.tkt.data_retain = set(["host_genus"])
        self.tkt.check_compatible_type_and_data_retain(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")




    def test_check_valid_data_source_1(self):
        """Check that no error is produced if
        data_add is empty and check_set is empty."""
        self.tkt.data_add = set()
        check_set = set()
        self.tkt.check_valid_data_source("data_add", check_set=check_set,
                                         eval_id="test")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "test")

    def test_check_valid_data_source_2(self):
        """Check that no error is produced if
        data_add is empty and check_set is not empty."""
        self.tkt.data_add = set()
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_add", check_set=check_set)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)

    def test_check_valid_data_source_3(self):
        """Check that no error is produced if
        data_add is not empty and valid and check_set is not empty."""
        self.tkt.data_add = set(["a"])
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_add", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "correct")

    def test_check_valid_data_source_4(self):
        """Check that error is produced if
        data_add is not empty and invalid and check_set is not empty."""
        self.tkt.data_add = set(["a", "c"])
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_add", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_valid_data_source_5(self):
        """Check that error is produced if
        data_retain is not empty and invalid and check_set is not empty."""
        self.tkt.data_retain = set(["a", "c"])
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_retain", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_valid_data_source_6(self):
        """Check that error is produced if
        data_retrieve is not empty and invalid and check_set is not empty."""
        self.tkt.data_retrieve = set(["a", "c"])
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_retrieve", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_valid_data_source_7(self):
        """Check that error is produced if
        data_parse is not empty and invalid and check_set is not empty."""
        self.tkt.data_parse = set(["a", "c"])
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("data_parse", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "error")

    def test_check_valid_data_source_8(self):
        """Check that error is produced if invalid ref_set_attr."""
        check_set = set(["a", "b"])
        self.tkt.check_valid_data_source("invalid", check_set=check_set)
        self.assertEqual(self.tkt.evaluations[0].status, "error")




if __name__ == '__main__':
    unittest.main()
