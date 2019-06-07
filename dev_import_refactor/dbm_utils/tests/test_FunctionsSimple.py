""" Unit tests for general functions."""


from functions import FunctionsSimple
import unittest
import re




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




    def test_identify_unique_items_1(self):
        """Verify the same list is returned."""
        input_list = ['a','b','c']
        unique_list, duplicate_list = \
            FunctionsSimple.identify_unique_items(input_list)
        unique_list.sort()
        with self.subTest():
            self.assertEqual(len(unique_list), 3)
        with self.subTest():
            self.assertEqual(unique_list, input_list)
        with self.subTest():
            self.assertEqual(len(duplicate_list), 0)

    def test_identify_unique_items_2(self):
        """Verify a unique list with no duplicates is returned, and a
        duplicate list is returned with one item."""
        input_list = ['a','b','c','c']
        expected_unique_list = ['a','b']
        expected_duplicate_list = ['c']
        unique_list, duplicate_list = \
            FunctionsSimple.identify_unique_items(input_list)
        unique_list.sort()
        with self.subTest():
            self.assertEqual(len(unique_list), 2)
        with self.subTest():
            self.assertEqual(unique_list, expected_unique_list)
        with self.subTest():
            self.assertEqual(duplicate_list, expected_duplicate_list)

    def test_identify_unique_items_3(self):
        """Verify a unique list with no items is returned, and a
        duplicate list is returned with three items."""
        input_list = ['a','b','c','c', 'a','b']
        expected_unique_list = []
        expected_duplicate_list = ['a','b','c']
        unique_list, duplicate_list = \
            FunctionsSimple.identify_unique_items(input_list)
        duplicate_list.sort()
        with self.subTest():
            self.assertEqual(len(unique_list), 0)
        with self.subTest():
            self.assertEqual(unique_list, expected_unique_list)
        with self.subTest():
            self.assertEqual(duplicate_list, expected_duplicate_list)

    def test_identify_unique_items_4(self):
        """Verify the function works with a different object type than
        strings."""
        input_list = [('a','b'), ('a','c'), ('c','a'), ('a','c')]
        expected_unique_list = [('a','b'), ('c','a')]
        expected_duplicate_list = [('a','c')]
        unique_list, duplicate_list = \
            FunctionsSimple.identify_unique_items(input_list)
        unique_list.sort()
        with self.subTest():
            self.assertEqual(len(unique_list), 2)
        with self.subTest():
            self.assertEqual(unique_list, expected_unique_list)
        with self.subTest():
            self.assertEqual(duplicate_list, expected_duplicate_list)





    def test_trim_characters_1(self):
        """Verify empty string is not changed."""
        input_string = ""
        output_string = FunctionsSimple.trim_characters(input_string)
        self.assertEqual(output_string, "")

    def test_trim_characters_2(self):
        """Verify string is not changed."""
        input_string = "abc"
        output_string = FunctionsSimple.trim_characters(input_string)
        self.assertEqual(output_string, input_string)

    def test_trim_characters_3(self):
        """Verify string with one leading character is trimmed."""
        input_string = ".abc"
        expected_output_string = "abc"
        output_string = FunctionsSimple.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)

    def test_trim_characters_4(self):
        """Verify string with two leading characters is trimmed."""
        input_string = ";.abc"
        expected_output_string = "abc"
        output_string = FunctionsSimple.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)

    def test_trim_characters_5(self):
        """Verify string with leading and trailing characters is trimmed."""
        input_string = ".,;-.,;-.,;-abc.,;-.,;-.,;-"
        expected_output_string = "abc"
        output_string = FunctionsSimple.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)








    def test_parse_names_from_record_field_1(self):
        """Test empty string."""
        string = ""
        expected_phage = ""
        expected_host = ""
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_2(self):
        """Verify host before 'phage' is identified."""
        string = "Mycobacterium phage"
        expected_phage = ""
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_3(self):
        """Verify host before 'virus' is identified."""
        string = "Mycobacterium virus"
        expected_phage = ""
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_4(self):
        """Verify phage after 'phage' is identified."""
        string = "phage Trixie"
        expected_phage = "Trixie"
        expected_host = ""
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_5(self):
        """Verify phage after 'virus' is identified."""
        string = "virus Trixie"
        expected_phage = "Trixie"
        expected_host = ""
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_6(self):
        """Verify both host and phage are identified."""
        string = "Mycobacterium phage Trixie"
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_7(self):
        """Verify both host and phage are identified from long description."""
        string = "skdj sjakl Mycobacterium phage Trixie skdjfl ksjd ksjdk"
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_8(self):
        """Verify host is identified from messy description."""
        string = ".;Mycobacterium.., phage"
        expected_phage = ""
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_9(self):
        """Verify phage is identified from messy description."""
        string = "phage .,;Trixie..;;,"
        expected_phage = "Trixie"
        expected_host = ""
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_10(self):
        """Verify host is identified from messy 'phage'."""
        string = ".;Mycobacterium.., .phage;"
        expected_phage = ""
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_11(self):
        """Verify host and phage are identified from string with
        multiple whitespace characters."""
        string = "  Mycobacterium       phage   Trixie   "
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_12(self):
        """Verify both host and phage are identified from
        long messy complex description."""
        string = " .kl;; .,Mycobacterium...  ;phage  Trixie,. complete genome."
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_13(self):
        """Verify host joined with 'phage' is identified from one word
        string."""
        string = "Mycobacteriophage"
        expected_phage = ""
        expected_host = "Mycobacterio"
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_14(self):
        """Verify phage is identified from short description."""
        string = "Trixie"
        expected_phage = "Trixie"
        expected_host = ""
        output_phage, output_host = \
            FunctionsSimple.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)









if __name__ == '__main__':
    unittest.main()
