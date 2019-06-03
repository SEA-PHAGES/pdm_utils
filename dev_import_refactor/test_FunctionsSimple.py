""" Unit tests for general functions."""


import FunctionsSimple
import unittest
import re
#import Genome




class TestGeneralFunctions(unittest.TestCase):


    def setUp(self):
        pass





    def test_remove_draft_suffix_1(self):
        """Verify '_Draft' suffix is removed."""
        new_value = FunctionsSimple.remove_draft_suffix("Phage_Draft")
        self.assertEqual(new_value, "Phage")


    def test_remove_draft_suffix_2(self):
        """Verify important data is not removed."""
        new_value = FunctionsSimple.remove_draft_suffix("PhageDraft")
        self.assertEqual(new_value, "PhageDraft")







    def test_reformat_strand_1(self):
        """Verify 'Forward' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand("Forward", "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_2(self):
        """Verify 'F' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand("F", "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_3(self):
        """Verify '1' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand(1, "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_4(self):
        """Verify 'Reverse' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand("Reverse", "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_5(self):
        """Verify 'R' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand("R", "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_6(self):
        """Verify '-1' strand is converted to long format."""
        new_value = FunctionsSimple.reformat_strand(-1, "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_7(self):
        """Verify non-standard strand is not converted."""
        new_value = FunctionsSimple.reformat_strand("None", "fr_long")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_8(self):
        """Verify 'Forward' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand("Forward", "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_9(self):
        """Verify 'F' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand("F", "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_10(self):
        """Verify '1' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand(1, "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_11(self):
        """Verify 'Reverse' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand("Reverse", "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_12(self):
        """Verify 'R' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand("R", "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_13(self):
        """Verify '-1' strand is converted to short format."""
        new_value = FunctionsSimple.reformat_strand(-1, "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_14(self):
        """Verify non-standard strand is not converted."""
        new_value = FunctionsSimple.reformat_strand("None", "fr_short")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_15(self):
        """Verify 'Forward' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand("Forward", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_16(self):
        """Verify 'F' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand("F", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_17(self):
        """Verify '1' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand(1, "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_18(self):
        """Verify 'Reverse' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand("Reverse", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_19(self):
        """Verify 'R' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand("R", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_20(self):
        """Verify '-1' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand(-1, "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_21(self):
        """Verify non-standard strand is not converted."""
        new_value = FunctionsSimple.reformat_strand("None", "numeric")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_22(self):
        """Verify non-standard format does not convert strand."""
        new_value = FunctionsSimple.reformat_strand("Forward", "other")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_23(self):
        """Verify 'Watson' strand is converted to numeric format."""
        new_value = FunctionsSimple.reformat_strand("Watson", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_24(self):
        """Verify 'Crick' strand is converted to operator format."""
        new_value = FunctionsSimple.reformat_strand("Crick", "operator")
        self.assertEqual(new_value, "-")

    def test_reformat_strand_25(self):
        """Verify 'operator' strand is converted to wc_long format."""
        new_value = FunctionsSimple.reformat_strand("+", "wc_long")
        self.assertEqual(new_value, "watson")

    def test_reformat_strand_26(self):
        """Verify capitalization for strings."""
        new_value = FunctionsSimple.reformat_strand("+", "wc_long", True)
        self.assertEqual(new_value, "Watson")

    def test_reformat_strand_27(self):
        """Verify capitalization does not impact numeric format."""
        new_value = FunctionsSimple.reformat_strand("+", "numeric", True)
        self.assertEqual(new_value, 1)

    def test_reformat_strand_28(self):
        """Verify capitalization does not impact operator format."""
        new_value = FunctionsSimple.reformat_strand("w", "operator", True)
        self.assertEqual(new_value, "+")





    def test_reformat_description_1(self):
        """Verify 'None' description is converted."""
        description = None
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_2(self):
        """Verify description is converted properly."""
        description = " Hypothetical Protein "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "Hypothetical Protein")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_3(self):
        """Verify description is converted properly."""
        description = "\\N "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "\\N")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_4(self):
        """Verify description is converted properly."""
        description = " . "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, ".")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_5(self):
        """Verify description is converted properly."""
        description = " 1 "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "1")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_6(self):
        """Verify 'gp' description is converted properly."""
        description = " gp12345 "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "gp12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_7(self):
        """Verify 'gp' description is converted properly."""
        description = " GP12345a "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "GP12345a")
        with self.subTest():
            self.assertEqual(value2, "GP12345a")

    def test_reformat_description_8(self):
        """Verify 'orf' description is converted properly."""
        description = " orf12345 "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "orf12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_9(self):
        """Verify 'orf' description is converted properly."""
        description = " ORF12345a "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "ORF12345a")
        with self.subTest():
            self.assertEqual(value2, "ORF12345a")

    def test_reformat_description_10(self):
        """Verify 'gp' description is converted properly."""
        description = " gp 12345 "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "gp 12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_11(self):
        """Verify 'gp' description is converted properly."""
        description = " GP 12345a "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "GP 12345a")
        with self.subTest():
            self.assertEqual(value2, "GP 12345a")

    def test_reformat_description_12(self):
        """Verify 'putative protein' description is converted properly."""
        description = " Putative Protein12345 "
        value1, value2 = FunctionsSimple.reformat_description(description)
        with self.subTest():
            self.assertEqual(value1, "Putative Protein12345")
        with self.subTest():
            self.assertEqual(value2, "")




    def test_find_expression_1(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["Trixie", "L5", "D29", "Trixie"]
        value = FunctionsSimple.find_expression(pattern, list_of_items)
        self.assertEqual(value, 2)

    def test_find_expression_2(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["L5", "D29"]
        value = FunctionsSimple.find_expression(pattern, list_of_items)
        self.assertEqual(value, 0)







    def test_parse_fasta_file_1(self):
        """Verify it fasta file data is parsed correctly."""
        fasta_data = ">Trixie     \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        sequence = FunctionsSimple.parse_fasta_file(fasta_data)
        self.assertEqual(expected_sequence, sequence)
















if __name__ == '__main__':
    unittest.main()
