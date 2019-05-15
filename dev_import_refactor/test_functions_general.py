""" Unit tests for general functions."""


import functions_general
import unittest
import re
import Genome




class TestGeneralFunctionsClass(unittest.TestCase):


    def setUp(self):
        pass





    def test_remove_draft_suffix_1(self):
        """Verify '_Draft' suffix is removed."""
        new_value = functions_general.remove_draft_suffix("Phage_Draft")
        self.assertEqual(new_value, "Phage")


    def test_remove_draft_suffix_2(self):
        """Verify important data is not removed."""
        new_value = functions_general.remove_draft_suffix("PhageDraft")
        self.assertEqual(new_value, "PhageDraft")







    def test_reformat_strand_1(self):
        """Verify 'Forward' strand is converted to long format."""
        new_value = functions_general.reformat_strand("Forward", "long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_2(self):
        """Verify 'F' strand is converted to long format."""
        new_value = functions_general.reformat_strand("F", "long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_3(self):
        """Verify '1' strand is converted to long format."""
        new_value = functions_general.reformat_strand(1, "long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_4(self):
        """Verify 'Reverse' strand is converted to long format."""
        new_value = functions_general.reformat_strand("Reverse", "long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_5(self):
        """Verify 'R' strand is converted to long format."""
        new_value = functions_general.reformat_strand("R", "long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_6(self):
        """Verify '-1' strand is converted to long format."""
        new_value = functions_general.reformat_strand(-1, "long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_7(self):
        """Verify non-standard strand is not converted."""
        new_value = functions_general.reformat_strand("None", "long")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_8(self):
        """Verify 'Forward' strand is converted to short format."""
        new_value = functions_general.reformat_strand("Forward", "short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_9(self):
        """Verify 'F' strand is converted to short format."""
        new_value = functions_general.reformat_strand("F", "short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_10(self):
        """Verify '1' strand is converted to short format."""
        new_value = functions_general.reformat_strand(1, "short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_11(self):
        """Verify 'Reverse' strand is converted to short format."""
        new_value = functions_general.reformat_strand("Reverse", "short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_12(self):
        """Verify 'R' strand is converted to short format."""
        new_value = functions_general.reformat_strand("R", "short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_13(self):
        """Verify '-1' strand is converted to short format."""
        new_value = functions_general.reformat_strand(-1, "short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_14(self):
        """Verify non-standard strand is not converted."""
        new_value = functions_general.reformat_strand("None", "short")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_15(self):
        """Verify 'Forward' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand("Forward", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_16(self):
        """Verify 'F' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand("F", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_17(self):
        """Verify '1' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand(1, "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_18(self):
        """Verify 'Reverse' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand("Reverse", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_19(self):
        """Verify 'R' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand("R", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_20(self):
        """Verify '-1' strand is converted to numeric format."""
        new_value = functions_general.reformat_strand(-1, "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_21(self):
        """Verify non-standard strand is not converted."""
        new_value = functions_general.reformat_strand("None", "numeric")
        self.assertEqual(new_value, 0)

    def test_reformat_strand_22(self):
        """Verify non-standard format does not convert strand."""
        new_value = functions_general.reformat_strand("Forward", "other")
        self.assertEqual(new_value, "Forward")













    def test_reformat_description_1(self):
        """Verify 'None' description is converted."""
        description = None
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_2(self):
        """Verify description is converted properly."""
        description = " Hypothetical Protein "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "Hypothetical Protein")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_3(self):
        """Verify description is converted properly."""
        description = "\\N "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "\\N")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_4(self):
        """Verify description is converted properly."""
        description = " . "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, ".")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_5(self):
        """Verify description is converted properly."""
        description = " 1 "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "1")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_6(self):
        """Verify 'gp' description is converted properly."""
        description = " gp12345 "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "gp12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_7(self):
        """Verify 'gp' description is converted properly."""
        description = " GP12345a "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "GP12345a")
        with self.subTest():
            self.assertEqual(value2, "GP12345a")

    def test_reformat_description_8(self):
        """Verify 'orf' description is converted properly."""
        description = " orf12345 "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "orf12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_9(self):
        """Verify 'orf' description is converted properly."""
        description = " ORF12345a "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "ORF12345a")
        with self.subTest():
            self.assertEqual(value2, "ORF12345a")

    def test_reformat_description_10(self):
        """Verify 'gp' description is converted properly."""
        description = " gp 12345 "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "gp 12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_11(self):
        """Verify 'gp' description is converted properly."""
        description = " GP 12345a "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "GP 12345a")
        with self.subTest():
            self.assertEqual(value2, "GP 12345a")

    def test_reformat_description_12(self):
        """Verify 'putative protein' description is converted properly."""
        description = " Putative Protein12345 "
        value1, value2 = functions_general.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "Putative Protein12345")
        with self.subTest():
            self.assertEqual(value2, "")




    def test_find_expression_1(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["Trixie", "L5", "D29", "Trixie"]
        value = functions_general.find_expression(pattern, list_of_items)
        self.assertEqual(value, 2)

    def test_find_expression_2(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["L5", "D29"]
        value = functions_general.find_expression(pattern, list_of_items)
        self.assertEqual(value, 0)










if __name__ == '__main__':
    unittest.main()
