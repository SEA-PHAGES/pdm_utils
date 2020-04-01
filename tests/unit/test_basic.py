"""Unit tests for general functions."""


from pdm_utils.functions import basic
from datetime import datetime
import unittest
import re




class TestBasicFunctions1(unittest.TestCase):

    def test_edit_suffix_1(self):
        """Verify '_Draft' suffix is removed."""
        new_value = basic.edit_suffix("Bonanza_Draft", "remove")
        self.assertEqual(new_value, "Bonanza")

    def test_edit_suffix_2(self):
        """Verify important data is not removed."""
        new_value = basic.edit_suffix("BonanzaDraft", "remove")
        self.assertEqual(new_value, "BonanzaDraft")

    def test_edit_suffix_3(self):
        """Verify '_Draft' suffix is added."""
        new_value = basic.edit_suffix("Bonanza", "add")
        self.assertEqual(new_value, "Bonanza_Draft")

    def test_edit_suffix_4(self):
        """Verify suffix is not added."""
        new_value = basic.edit_suffix("Bonanza_Draft", "add")
        self.assertEqual(new_value, "Bonanza_Draft")

    def test_edit_suffix_5(self):
        """Verify '_Draft' suffix is added."""
        new_value = basic.edit_suffix("BonanzaDraft", "add")
        self.assertEqual(new_value, "BonanzaDraft_Draft")

    def test_edit_suffix_6(self):
        """Verify value is not changed due to invalid option."""
        new_value = basic.edit_suffix("Bonanza_Draft", "invalid")
        self.assertEqual(new_value, "Bonanza_Draft")




    def test_create_indices_1(self):
        """Verify correct indices produced from batch_size of 1."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 1)
        with self.subTest():
            self.assertEqual(len(new_list),5)
        with self.subTest():
            self.assertEqual(new_list[4], (4,5))

    def test_create_indices_2(self):
        """Verify correct indices produced from batch_size of 2."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 2)
        with self.subTest():
            self.assertEqual(len(new_list),3)
        with self.subTest():
            self.assertEqual(new_list[0], (0,2))
        with self.subTest():
            self.assertEqual(new_list[2], (4,5))

    def test_create_indices_3(self):
        """Verify correct indices produced from batch_size of 3."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 3)
        with self.subTest():
            self.assertEqual(len(new_list),2)
        with self.subTest():
            self.assertEqual(new_list[0], (0,3))
        with self.subTest():
            self.assertEqual(new_list[1], (3,5))

    def test_create_indices_4(self):
        """Verify correct indices produced from batch_size of 4."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 4)
        with self.subTest():
            self.assertEqual(len(new_list),2)
        with self.subTest():
            self.assertEqual(new_list[0], (0,4))
        with self.subTest():
            self.assertEqual(new_list[1], (4,5))

    def test_create_indices_5(self):
        """Verify correct indices produced from batch_size of 5."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 5)
        with self.subTest():
            self.assertEqual(len(new_list),1)
        with self.subTest():
            self.assertEqual(new_list[0], (0,5))

    def test_create_indices_6(self):
        """Verify correct indices produced from batch_size of 6."""
        input_list = ["a", "b", "c", "d", "e"]
        new_list = basic.create_indices(input_list, 6)
        with self.subTest():
            self.assertEqual(len(new_list),1)
        with self.subTest():
            self.assertEqual(new_list[0], (0,5))




    def test_reformat_strand_1(self):
        """Verify 'Forward' strand is converted to long format."""
        new_value = basic.reformat_strand("Forward", "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_2(self):
        """Verify 'F' strand is converted to long format."""
        new_value = basic.reformat_strand("F", "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_3(self):
        """Verify '1' strand is converted to long format."""
        new_value = basic.reformat_strand(1, "fr_long")
        self.assertEqual(new_value, "forward")

    def test_reformat_strand_4(self):
        """Verify 'Reverse' strand is converted to long format."""
        new_value = basic.reformat_strand("Reverse", "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_5(self):
        """Verify 'R' strand is converted to long format."""
        new_value = basic.reformat_strand("R", "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_6(self):
        """Verify '-1' strand is converted to long format."""
        new_value = basic.reformat_strand(-1, "fr_long")
        self.assertEqual(new_value, "reverse")

    def test_reformat_strand_7(self):
        """Verify non-standard strand is not converted."""
        new_value = basic.reformat_strand("None", "fr_long")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_8(self):
        """Verify 'Forward' strand is converted to short format."""
        new_value = basic.reformat_strand("Forward", "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_9(self):
        """Verify 'F' strand is converted to short format."""
        new_value = basic.reformat_strand("F", "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_10(self):
        """Verify '1' strand is converted to short format."""
        new_value = basic.reformat_strand(1, "fr_short")
        self.assertEqual(new_value, "f")

    def test_reformat_strand_11(self):
        """Verify 'Reverse' strand is converted to short format."""
        new_value = basic.reformat_strand("Reverse", "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_12(self):
        """Verify 'R' strand is converted to short format."""
        new_value = basic.reformat_strand("R", "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_13(self):
        """Verify '-1' strand is converted to short format."""
        new_value = basic.reformat_strand(-1, "fr_short")
        self.assertEqual(new_value, "r")

    def test_reformat_strand_14(self):
        """Verify non-standard strand is not converted."""
        new_value = basic.reformat_strand("None", "fr_short")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_15(self):
        """Verify 'Forward' strand is converted to numeric format."""
        new_value = basic.reformat_strand("Forward", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_16(self):
        """Verify 'F' strand is converted to numeric format."""
        new_value = basic.reformat_strand("F", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_17(self):
        """Verify '1' strand is converted to numeric format."""
        new_value = basic.reformat_strand(1, "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_18(self):
        """Verify 'Reverse' strand is converted to numeric format."""
        new_value = basic.reformat_strand("Reverse", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_19(self):
        """Verify 'R' strand is converted to numeric format."""
        new_value = basic.reformat_strand("R", "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_20(self):
        """Verify '-1' strand is converted to numeric format."""
        new_value = basic.reformat_strand(-1, "numeric")
        self.assertEqual(new_value, -1)

    def test_reformat_strand_21(self):
        """Verify non-standard strand is not converted."""
        new_value = basic.reformat_strand("None", "numeric")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_22(self):
        """Verify non-standard format does not convert strand."""
        new_value = basic.reformat_strand("Forward", "other")
        self.assertEqual(new_value, "NA")

    def test_reformat_strand_23(self):
        """Verify 'Watson' strand is converted to numeric format."""
        new_value = basic.reformat_strand("Watson", "numeric")
        self.assertEqual(new_value, 1)

    def test_reformat_strand_24(self):
        """Verify 'Crick' strand is converted to operator format."""
        new_value = basic.reformat_strand("Crick", "operator")
        self.assertEqual(new_value, "-")

    def test_reformat_strand_25(self):
        """Verify 'operator' strand is converted to wc_long format."""
        new_value = basic.reformat_strand("+", "wc_long")
        self.assertEqual(new_value, "watson")

    def test_reformat_strand_26(self):
        """Verify capitalization for strings."""
        new_value = basic.reformat_strand("+", "wc_long", True)
        self.assertEqual(new_value, "Watson")

    def test_reformat_strand_27(self):
        """Verify capitalization does not impact numeric format."""
        new_value = basic.reformat_strand("+", "numeric", True)
        self.assertEqual(new_value, 1)

    def test_reformat_strand_28(self):
        """Verify capitalization does not impact operator format."""
        new_value = basic.reformat_strand("w", "operator", True)
        self.assertEqual(new_value, "+")

    def test_reformat_strand_29(self):
        """Verify 'F' strand is converted to first abbreviated format."""
        new_value = basic.reformat_strand("F", "fr_abbrev1")
        self.assertEqual(new_value, "for")

    def test_reformat_strand_30(self):
        """Verify 'F' strand is converted to second abbreviated format."""
        new_value = basic.reformat_strand("F", "fr_abbrev2")
        self.assertEqual(new_value, "fwd")








    def test_reformat_coordinates_1(self):
        """Verify 0-based half open interval is converted to
        1-based closed interval."""
        output_start, output_stop = \
            basic.reformat_coordinates(5, 10, "0_half_open", "1_closed")
        with self.subTest():
            self.assertEqual(output_start, 6)
        with self.subTest():
            self.assertEqual(output_stop, 10)

    def test_reformat_coordinates_2(self):
        """Verify 0-based half open interval is not converted."""
        output_start, output_stop = \
            basic.reformat_coordinates(5, 10, "0_half_open", "0_half_open")
        with self.subTest():
            self.assertEqual(output_start, 5)
        with self.subTest():
            self.assertEqual(output_stop, 10)

    def test_reformat_coordinates_3(self):
        """Verify 1-based closed interval is converted to
        0-based half open interval."""
        output_start, output_stop = \
            basic.reformat_coordinates(5, 10, "1_closed", "0_half_open")
        with self.subTest():
            self.assertEqual(output_start, 4)
        with self.subTest():
            self.assertEqual(output_stop, 10)

    def test_reformat_coordinates_4(self):
        """Verify 1-based closed interval is not converted."""
        output_start, output_stop = \
            basic.reformat_coordinates(5, 10, "1_closed", "1_closed")
        with self.subTest():
            self.assertEqual(output_start, 5)
        with self.subTest():
            self.assertEqual(output_stop, 10)

    def test_reformat_coordinates_5(self):
        """Verify invalid input format is not converted."""
        with self.assertRaises(ValueError):
            basic.reformat_coordinates(5, 10, "invalid", "1_closed")

    def test_reformat_coordinates_6(self):
        """Verify invalid output format is not converted."""
        with self.assertRaises(ValueError):
            basic.reformat_coordinates(5, 10, "1_closed", "invalid")




    def test_check_empty_1(self):
        """Verify 'a' is not in the empty set."""
        result = basic.check_empty("a")
        self.assertFalse(result)

    def test_check_empty_2(self):
        """Verify '' is in the empty set."""
        result = basic.check_empty("")
        self.assertTrue(result)

    def test_check_empty_3(self):
        """Verify 'none' is in the empty set."""
        result = basic.check_empty("none")
        self.assertTrue(result)

    def test_check_empty_4(self):
        """Verify 'null' is in the empty set."""
        result = basic.check_empty("null")
        self.assertTrue(result)

    def test_check_empty_5(self):
        """Verify None is in the empty set."""
        result = basic.check_empty(None)
        self.assertTrue(result)

    def test_check_empty_6(self):
        """Verify 'not applicable' is in the empty set."""
        result = basic.check_empty("not applicable")
        self.assertTrue(result)

    def test_check_empty_7(self):
        """Verify 'na' is in the empty set."""
        result = basic.check_empty("na")
        self.assertTrue(result)

    def test_check_empty_8(self):
        """Verify 'n/a' is in the empty set."""
        result = basic.check_empty("n/a")
        self.assertTrue(result)

    def test_check_empty_9(self):
        """Verify '0' is in the empty set."""
        result = basic.check_empty("0")
        self.assertTrue(result)

    def test_check_empty_10(self):
        """Verify 0 is in the empty set."""
        result = basic.check_empty(0)
        self.assertTrue(result)

    def test_check_empty_11(self):
        """Verify empty datetime object is in the empty set."""
        empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
        result = basic.check_empty(empty_date)
        self.assertTrue(result)

    def test_check_empty_12(self):
        """Verify 'NA' is in the empty set."""
        result = basic.check_empty("NA")
        self.assertTrue(result)

    def test_check_empty_13(self):
        """Verify 'NA' is now not in the empty set."""
        result = basic.check_empty("NA", lower = False)
        self.assertFalse(result)

    def test_check_empty_14(self):
        """Verify None is still in the empty set."""
        result = basic.check_empty(None, lower = False)
        self.assertTrue(result)

    def test_check_empty_15(self):
        """Verify 0 is still in the empty set."""
        result = basic.check_empty(0, lower = False)
        self.assertTrue(result)

    def test_check_empty_16(self):
        """Verify empty datetime object is still in the empty set."""
        empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
        result = basic.check_empty(empty_date, lower = False)
        self.assertTrue(result)




    def test_convert_empty_1(self):
        """Verify '' is converted to 'none'."""
        new_value = basic.convert_empty("", "none_string")
        self.assertEqual(new_value, "none")

    def test_convert_empty_2(self):
        """Verify '' is converted to 'null'."""
        new_value = basic.convert_empty("", "null_string")
        self.assertEqual(new_value, "null")

    def test_convert_empty_3(self):
        """Verify '' is converted to None."""
        new_value = basic.convert_empty("", "none_object")
        self.assertEqual(new_value, None)

    def test_convert_empty_4(self):
        """Verify '' is converted to 'not applicable'."""
        new_value = basic.convert_empty("", "na_long")
        self.assertEqual(new_value, "not applicable")

    def test_convert_empty_5(self):
        """Verify '' is converted to 'na'."""
        new_value = basic.convert_empty("", "na_short")
        self.assertEqual(new_value, "na")

    def test_convert_empty_6(self):
        """Verify '' is converted to '0'."""
        new_value = basic.convert_empty("", "zero_string")
        self.assertEqual(new_value, "0")

    def test_convert_empty_7(self):
        """Verify '' is converted to 0."""
        new_value = basic.convert_empty("", "zero_num")
        self.assertEqual(new_value, 0)

    def test_convert_empty_8(self):
        """Verify 'NA' is converted to 0."""
        new_value = basic.convert_empty("NA", "zero_num")
        self.assertEqual(new_value, 0)

    def test_convert_empty_9(self):
        """Verify 'NONE' is converted to 0."""
        new_value = basic.convert_empty("NONE", "zero_num")
        self.assertEqual(new_value, 0)

    def test_convert_empty_10(self):
        """Verify None is converted to 0."""
        new_value = basic.convert_empty(None, "zero_num")
        self.assertEqual(new_value, 0)

    def test_convert_empty_11(self):
        """Verify 0 is converted to '0'."""
        new_value = basic.convert_empty(0, "zero_string")
        self.assertEqual(new_value, "0")

    def test_convert_empty_12(self):
        """Verify 0 is converted to uppercased 'NA'."""
        new_value = basic.convert_empty(0, "na_short", upper = True)
        self.assertEqual(new_value, "NA")

    def test_convert_empty_13(self):
        """Verify 'NONE' is not converted due to invalid format."""
        new_value = basic.convert_empty("NONE", "asdf")
        self.assertEqual(new_value, "NONE")

    def test_convert_empty_14(self):
        """Verify 0 is not converted due to invalid format."""
        new_value = basic.convert_empty(0, "asdf")
        self.assertEqual(new_value, 0)

    def test_convert_empty_15(self):
        """Verify None is not converted due to invalid format."""
        new_value = basic.convert_empty(None, "asdf")
        self.assertEqual(new_value, None)

    def test_convert_empty_16(self):
        """Verify 'asdf' is not converted due to invalid input value."""
        new_value = basic.convert_empty("Asdf", "none_string")
        self.assertEqual(new_value, "Asdf")

    def test_convert_empty_17(self):
        """Verify 'NA' is converted to empty date."""
        empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
        new_value = basic.convert_empty("NA", "empty_datetime_obj")
        self.assertEqual(new_value, empty_date)

    def test_convert_empty_18(self):
        """Verify empty date is converted to "none"."""
        empty_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
        new_value = basic.convert_empty(empty_date, "none_string")
        self.assertEqual(new_value, "none")

    def test_convert_empty_19(self):
        """Verify 'NONE' is converted to 'N/A'."""
        new_value = basic.convert_empty("NONE", "n/a", upper = True)
        self.assertEqual(new_value, "N/A")

    def test_convert_empty_20(self):
        """Verify 'N/A' is converted to ""."""
        new_value = basic.convert_empty("N/A", "empty_string")
        self.assertEqual(new_value, "")




    def test_reformat_description_1(self):
        """Verify 'None' raw_description is converted."""
        raw_description = None
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_2(self):
        """Verify raw_description is converted properly."""
        raw_description = " Hypothetical Protein "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "Hypothetical Protein")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_3(self):
        """Verify raw_description is converted properly."""
        raw_description = "\\N "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "\\N")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_4(self):
        """Verify raw_description is converted properly."""
        raw_description = " . "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, ".")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_5(self):
        """Verify raw_description is converted properly."""
        raw_description = " 1 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "1")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_6(self):
        """Verify 'gp' raw_description is converted properly."""
        raw_description = " gp12345 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "gp12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_7(self):
        """Verify 'gp' raw_description is converted properly."""
        raw_description = " GP12345a "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "GP12345a")
        with self.subTest():
            self.assertEqual(value2, "GP12345a")

    def test_reformat_description_8(self):
        """Verify 'orf' raw_description is converted properly."""
        raw_description = " orf12345 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "orf12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_9(self):
        """Verify 'orf' raw_description is converted properly."""
        raw_description = " ORF12345a "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "ORF12345a")
        with self.subTest():
            self.assertEqual(value2, "ORF12345a")

    def test_reformat_description_10(self):
        """Verify 'gp' raw_description is converted properly."""
        raw_description = " gp 12345 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "gp 12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_11(self):
        """Verify 'gp' raw_description is converted properly."""
        raw_description = " GP 12345a "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "GP 12345a")
        with self.subTest():
            self.assertEqual(value2, "GP 12345a")

    def test_reformat_description_12(self):
        """Verify 'putative protein' raw_description is converted properly."""
        raw_description = " Putative Protein12345 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "Putative Protein12345")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_13(self):
        """Verify 'gp10.1' raw_description is converted properly."""
        raw_description = " gp10.1 "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "gp10.1")
        with self.subTest():
            self.assertEqual(value2, "")

    def test_reformat_description_14(self):
        """Verify 'gp10.1 protein' raw_description is converted properly."""
        raw_description = " gp10.1 protein "
        value1, value2 = basic.reformat_description(raw_description)
        with self.subTest():
            self.assertEqual(value1, "gp10.1 protein")
        with self.subTest():
            self.assertEqual(value2, "")




    def test_find_expression_1(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["Trixie", "L5", "D29", "Trixie"]
        value = basic.find_expression(pattern, list_of_items)
        self.assertEqual(value, 2)

    def test_find_expression_2(self):
        """Verify it returns correct number of items."""
        pattern = re.compile("^" + "Trixie" + "$")
        list_of_items = ["L5", "D29"]
        value = basic.find_expression(pattern, list_of_items)
        self.assertEqual(value, 0)




    def test_identify_unique_items_1(self):
        """Verify the same set is returned."""
        input_list = ['a','b','c']
        unique_set, duplicate_set = \
            basic.identify_unique_items(input_list)
        expected_unique_set = set(input_list)
        with self.subTest():
            self.assertEqual(len(unique_set), 3)
        with self.subTest():
            self.assertEqual(unique_set, expected_unique_set)
        with self.subTest():
            self.assertEqual(len(duplicate_set), 0)

    def test_identify_unique_items_2(self):
        """Verify a unique set with no duplicates is returned, and a
        duplicate set is returned with one item."""
        input_list = ['a','b','c','c']
        expected_unique_set = set(['a','b'])
        expected_duplicate_set = set(['c'])
        unique_set, duplicate_set = \
            basic.identify_unique_items(input_list)
        with self.subTest():
            self.assertEqual(len(unique_set), 2)
        with self.subTest():
            self.assertEqual(unique_set, expected_unique_set)
        with self.subTest():
            self.assertEqual(duplicate_set, expected_duplicate_set)

    def test_identify_unique_items_3(self):
        """Verify a unique set with no items is returned, and a
        duplicate set is returned with three items."""
        input_list = ['a','b','c','c', 'a','b']
        expected_unique_set = set([])
        expected_duplicate_set = set(['a','b','c'])
        unique_set, duplicate_set = \
            basic.identify_unique_items(input_list)
        with self.subTest():
            self.assertEqual(len(unique_set), 0)
        with self.subTest():
            self.assertEqual(unique_set, expected_unique_set)
        with self.subTest():
            self.assertEqual(duplicate_set, expected_duplicate_set)

    def test_identify_unique_items_4(self):
        """Verify the function works with a different object type than
        strings."""
        input_list = [('a','b'), ('a','c'), ('c','a'), ('a','c')]
        expected_unique_set = set([('a','b'), ('c','a')])
        expected_duplicate_set = set([('a','c')])
        unique_set, duplicate_set = \
            basic.identify_unique_items(input_list)
        with self.subTest():
            self.assertEqual(len(unique_set), 2)
        with self.subTest():
            self.assertEqual(unique_set, expected_unique_set)
        with self.subTest():
            self.assertEqual(duplicate_set, expected_duplicate_set)





    def test_trim_characters_1(self):
        """Verify empty string is not changed."""
        input_string = ""
        output_string = basic.trim_characters(input_string)
        self.assertEqual(output_string, "")

    def test_trim_characters_2(self):
        """Verify string is not changed."""
        input_string = "abc"
        output_string = basic.trim_characters(input_string)
        self.assertEqual(output_string, input_string)

    def test_trim_characters_3(self):
        """Verify string with one leading character is trimmed."""
        input_string = ".abc"
        expected_output_string = "abc"
        output_string = basic.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)

    def test_trim_characters_4(self):
        """Verify string with two leading characters is trimmed."""
        input_string = ";.abc"
        expected_output_string = "abc"
        output_string = basic.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)

    def test_trim_characters_5(self):
        """Verify string with leading and trailing characters is trimmed."""
        input_string = ".,;-.,;-.,;-abc.,;-.,;-.,;-"
        expected_output_string = "abc"
        output_string = basic.trim_characters(input_string)
        self.assertEqual(output_string, expected_output_string)




    def test_parse_names_from_record_field_1(self):
        """Test empty string."""
        string = ""
        expected_phage = ""
        expected_host = ""
        output_phage, output_host = \
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
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
            basic.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_15(self):
        """Verify 'Unclassified' is not identified."""
        string = "Unclassified."
        expected_phage = ""
        expected_host = ""
        output_phage, output_host = \
            basic.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)

    def test_parse_names_from_record_field_16(self):
        """Verify phage is identified."""
        string = "Trixie Unclassified."
        expected_phage = "Trixie"
        expected_host = ""
        output_phage, output_host = \
            basic.parse_names_from_record_field(string)
        with self.subTest():
            self.assertEqual(output_phage, expected_phage)
        with self.subTest():
            self.assertEqual(output_host, expected_host)




    def test_compare_sets_1(self):
        """Verify output when there is no intersection."""
        set1 = set(['a', 'b', 'c'])
        set2 = set(['d', 'e', 'f'])

        expected_set1_diff = set1
        expected_set2_diff = set2
        expected_set_int = set([])

        set_intersection, set1_diff, set2_diff = \
            basic.compare_sets(set1, set2)
        with self.subTest():
            self.assertEqual(set_intersection, expected_set_int)
        with self.subTest():
            self.assertEqual(set1_diff, expected_set1_diff)
        with self.subTest():
            self.assertEqual(set2_diff, expected_set2_diff)

    def test_compare_sets_2(self):
        """Verify output when there is one shared item."""
        set1 = set(['a', 'b', 'c'])
        set2 = set(['a', 'e', 'f'])

        expected_set1_diff = set(['b', 'c'])
        expected_set2_diff = set(['e', 'f'])
        expected_set_int = set(['a'])

        set_intersection, set1_diff, set2_diff = \
            basic.compare_sets(set1, set2)
        with self.subTest():
            self.assertEqual(set_intersection, expected_set_int)
        with self.subTest():
            self.assertEqual(set1_diff, expected_set1_diff)
        with self.subTest():
            self.assertEqual(set2_diff, expected_set2_diff)

    def test_compare_sets_3(self):
        """Verify output when there is all items in one set are shared."""
        set1 = set(['a', 'b', 'c'])
        set2 = set(['a', 'b', 'c'])

        expected_set1_diff = set([])
        expected_set2_diff = set([])
        expected_set_int = set1

        set_intersection, set1_diff, set2_diff = \
            basic.compare_sets(set1, set2)
        with self.subTest():
            self.assertEqual(set_intersection, expected_set_int)
        with self.subTest():
            self.assertEqual(set1_diff, expected_set1_diff)
        with self.subTest():
            self.assertEqual(set2_diff, expected_set2_diff)

    def test_compare_sets_4(self):
        """Verify output when there all items are shared."""
        set1 = set(['a', 'b', 'c'])
        set2 = set(['a', 'b', 'c', 'd'])

        expected_set1_diff = set([])
        expected_set2_diff = set(['d'])
        expected_set_int = set(['a', 'b', 'c'])

        set_intersection, set1_diff, set2_diff = \
            basic.compare_sets(set1, set2)
        with self.subTest():
            self.assertEqual(set_intersection, expected_set_int)
        with self.subTest():
            self.assertEqual(set1_diff, expected_set1_diff)
        with self.subTest():
            self.assertEqual(set2_diff, expected_set2_diff)




    def test_match_items_1(self):
        """Verify all unique items are matched."""
        input_list1 = ['a', 'b', 'c']
        input_list2 = ['a', 'b', 'c']

        exp_matched_unique = set(['a', 'b', 'c'])
        exp_set1_unmatched_unique = set([])
        exp_set2_unmatched_unique = set([])
        exp_set1_duplicate = set([])
        exp_set2_duplicate = set([])

        matched_unique, \
        set1_unmatched_unique, \
        set2_unmatched_unique, \
        set1_duplicate, \
        set2_duplicate = \
            basic.match_items(input_list1, input_list2)

        with self.subTest():
            self.assertEqual(matched_unique, exp_matched_unique)
        with self.subTest():
            self.assertEqual(set1_unmatched_unique, exp_set1_unmatched_unique)
        with self.subTest():
            self.assertEqual(set2_unmatched_unique, exp_set2_unmatched_unique)
        with self.subTest():
            self.assertEqual(set1_duplicate, exp_set1_duplicate)
        with self.subTest():
            self.assertEqual(set2_duplicate, exp_set2_duplicate)

    def test_match_items_2(self):
        """Verify no unique items are matched."""
        input_list1 = ['a', 'b', 'c']
        input_list2 = ['d', 'e', 'f']

        exp_matched_unique = set([])
        exp_set1_unmatched_unique = set(input_list1)
        exp_set2_unmatched_unique = set(input_list2)
        exp_set1_duplicate = set([])
        exp_set2_duplicate = set([])

        matched_unique, \
        set1_unmatched_unique, \
        set2_unmatched_unique, \
        set1_duplicate, \
        set2_duplicate = \
            basic.match_items(input_list1, input_list2)

        with self.subTest():
            self.assertEqual(matched_unique, exp_matched_unique)
        with self.subTest():
            self.assertEqual(set1_unmatched_unique, exp_set1_unmatched_unique)
        with self.subTest():
            self.assertEqual(set2_unmatched_unique, exp_set2_unmatched_unique)
        with self.subTest():
            self.assertEqual(set1_duplicate, exp_set1_duplicate)
        with self.subTest():
            self.assertEqual(set2_duplicate, exp_set2_duplicate)

    def test_match_items_3(self):
        """Verify no non-unique items are matched."""
        input_list1 = ['a', 'b', 'c'] + ['a', 'b', 'c']
        input_list2 = ['d', 'e', 'f'] + ['d', 'e', 'f']

        exp_matched_unique = set([])
        exp_set1_unmatched_unique = set([])
        exp_set2_unmatched_unique = set([])
        exp_set1_duplicate = set(input_list1)
        exp_set2_duplicate = set(input_list2)

        matched_unique, \
        set1_unmatched_unique, \
        set2_unmatched_unique, \
        set1_duplicate, \
        set2_duplicate = \
            basic.match_items(input_list1, input_list2)

        with self.subTest():
            self.assertEqual(matched_unique, exp_matched_unique)
        with self.subTest():
            self.assertEqual(set1_unmatched_unique, exp_set1_unmatched_unique)
        with self.subTest():
            self.assertEqual(set2_unmatched_unique, exp_set2_unmatched_unique)
        with self.subTest():
            self.assertEqual(set1_duplicate, exp_set1_duplicate)
        with self.subTest():
            self.assertEqual(set2_duplicate, exp_set2_duplicate)

    def test_match_items_4(self):
        """Verify unique and non-unique items are sorted."""
        input_list1 = ['a', 'b', 'c'] + ['a', 'b', 'c'] + ['z'] + ['x']
        input_list2 = ['d', 'e', 'f'] + ['d', 'e', 'f'] + ['z'] + ['y']

        exp_matched_unique = set(['z'])
        exp_set1_unmatched_unique = set(['x'])
        exp_set2_unmatched_unique = set(['y'])
        exp_set1_duplicate = set(['a', 'b', 'c'])
        exp_set2_duplicate = set(['d', 'e', 'f'])

        matched_unique, \
        set1_unmatched_unique, \
        set2_unmatched_unique, \
        set1_duplicate, \
        set2_duplicate = \
            basic.match_items(input_list1, input_list2)

        with self.subTest():
            self.assertEqual(matched_unique, exp_matched_unique)
        with self.subTest():
            self.assertEqual(set1_unmatched_unique, exp_set1_unmatched_unique)
        with self.subTest():
            self.assertEqual(set2_unmatched_unique, exp_set2_unmatched_unique)
        with self.subTest():
            self.assertEqual(set1_duplicate, exp_set1_duplicate)
        with self.subTest():
            self.assertEqual(set2_duplicate, exp_set2_duplicate)




    def test_split_string_1(self):
        """Verify non-numeric string is not split."""
        left, right = basic.split_string("ABCD")
        with self.subTest():
            self.assertEqual(left, "ABCD")
        with self.subTest():
            self.assertEqual(right, "")

    def test_split_string_2(self):
        """Verify numeric string is not split."""
        left, right = basic.split_string("1234")
        with self.subTest():
            self.assertEqual(left, "")
        with self.subTest():
            self.assertEqual(right, "1234")

    def test_split_string_3(self):
        """Verify alphanumeric string is split correctly."""
        left, right = basic.split_string("ABC0123")
        with self.subTest():
            self.assertEqual(left, "ABC")
        with self.subTest():
            self.assertEqual(right, "0123")

    def test_split_string_4(self):
        """Verify mixed alphanumeric string is not split."""
        left, right = basic.split_string("A1B2C3")
        with self.subTest():
            self.assertEqual(left, "")
        with self.subTest():
            self.assertEqual(right, "")

    def test_split_string_5(self):
        """Verify float string is not split."""
        left, right = basic.split_string("1234.5")
        with self.subTest():
            self.assertEqual(left, "")
        with self.subTest():
            self.assertEqual(right, "1234.5")

    def test_split_string_6(self):
        """Verify alphanumeric string with float is split correctly."""
        input = "ABC123.5"
        left, right = basic.split_string(input)
        with self.subTest():
            self.assertEqual(left, "ABC")
        with self.subTest():
            self.assertEqual(right, "123.5")




    def test_compare_cluster_subcluster_1(self):
        """Check that empty cluster and empty subcluster
        do not produce an error."""
        cluster = "none"
        subcluster = "none"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertTrue(result)

    def test_compare_cluster_subcluster_2(self):
        """Check that empty cluster and non-empty subcluster
        produce an error."""
        cluster = "none"
        subcluster = "AA1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_3(self):
        """Check that singleton cluster and empty subcluster
        do not produce an error."""
        cluster = "singleton"
        subcluster = "none"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertTrue(result)

    def test_compare_cluster_subcluster_4(self):
        """Check that singleton cluster and non-empty subcluster
        produce an error."""
        cluster = "singleton"
        subcluster = "AA1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_5(self):
        """Check that Unknown cluster and empty subcluster
        do not produce an error."""
        cluster = "UNK"
        subcluster = "none"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertTrue(result)

    def test_compare_cluster_subcluster_6(self):
        """Check that Unknown cluster and non-empty subcluster
        produce an error."""
        cluster = "UNK"
        subcluster = "AA1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_7(self):
        """Check that non-empty cluster and empty subcluster
        do not produce an error."""
        cluster = "AA"
        subcluster = "none"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertTrue(result)

    def test_compare_cluster_subcluster_8(self):
        """Check that non-empty cluster and matching non-empty subcluster
        do not produce an error."""
        cluster = "AA"
        subcluster = "AA1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertTrue(result)

    def test_compare_cluster_subcluster_9(self):
        """Check that non-empty cluster and non-matching non-empty subcluster
        produce an error."""
        cluster = "AA"
        subcluster = "A1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_10(self):
        """Check that non-empty cluster and non-matching non-empty subcluster
        produce an error."""
        cluster = "AA"
        subcluster = "1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_11(self):
        """Check that non-empty cluster and non-matching non-empty subcluster
        produce an error."""
        cluster = "AA"
        subcluster = "AAA1"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)

    def test_compare_cluster_subcluster_12(self):
        """Check that non-empty cluster and non-matching non-empty subcluster
        produce an error."""
        cluster = "AA"
        subcluster = "A123"
        result = basic.compare_cluster_subcluster(cluster, subcluster)
        self.assertFalse(result)




    def test_identify_one_list_duplicates_1(self):
        """Verify duplicate items generates an error."""
        duplicate_items = \
            basic.identify_one_list_duplicates(["Trixie", "Trixie"])
        self.assertEqual(len(duplicate_items), 1)

    def test_identify_one_list_duplicates_2(self):
        """Verify non-duplicate items do not generate an error."""
        duplicate_items = \
            basic.identify_one_list_duplicates(["Trixie", "L5"])
        self.assertEqual(len(duplicate_items), 0)

    def test_identify_one_list_duplicates_3(self):
        """Verify multiple duplicate items generate an error."""
        item_list = ["Trixie", "Trixie", "L5", "L5", "L5", "D29"]
        duplicate_items = \
            basic.identify_one_list_duplicates(item_list)
        self.assertEqual(len(duplicate_items), 2)




    def test_identify_two_list_duplicates_1(self):
        """Verify duplicate items in two lists generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["Trixie", "RedRock"]
        duplicate_items = \
            basic.identify_two_list_duplicates(list1, list2)
        self.assertEqual(len(duplicate_items), 1)

    def test_identify_two_list_duplicates_2(self):
        """Verify non-duplicate items in two lists do not generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["D29", "RedRock"]
        duplicate_items = \
            basic.identify_two_list_duplicates(list1, list2)
        self.assertEqual(len(duplicate_items), 0)

    def test_identify_two_list_duplicates_3(self):
        """Verify multiple duplicate items in two lists generate an error."""
        list1 = ["Trixie", "L5", "D29", "RedRock"]
        list2 = ["D29", "RedRock"]
        duplicate_items = \
            basic.identify_two_list_duplicates(list1, list2)
        self.assertEqual(len(duplicate_items), 2)




    def test_check_value_expected_in_set_1(self):
        """Check value that is in the set and expected to be in set."""
        set1 = set(["none", "empty"])
        value = "none"
        result = basic.check_value_expected_in_set(value, set1)
        self.assertTrue(result)

    def test_check_value_expected_in_set_2(self):
        """Check value that is in the set and not expected to be in set."""
        set1 = set(["none", "empty"])
        value = "none"
        result = basic.check_value_expected_in_set(value, set1, False)
        self.assertFalse(result)

    def test_check_value_expected_in_set_3(self):
        """Check value that is not in the set and expected to be in set."""
        set1 = set(["none", "empty"])
        value = "EagleEye"
        result = basic.check_value_expected_in_set(value, set1)
        self.assertFalse(result)

    def test_check_value_expected_in_set_4(self):
        """Check value that is not in the set and not expected to be in set."""
        set1 = set(["none", "empty"])
        value = "EagleEye"
        result = basic.check_value_expected_in_set(value, set1, False)
        self.assertTrue(result)




    def test_check_value_in_two_sets_1(self):
        """Check value that is in neither set."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["D29", "RedRock"])
        value = "EagleEye"
        result = basic.check_value_in_two_sets(value, set1, set2)
        self.assertEqual(result, "neither")

    def test_check_value_in_two_sets_2(self):
        """Check value that is in the first set."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["D29", "RedRock"])
        value = "Trixie"
        result = basic.check_value_in_two_sets(value, set1, set2)
        self.assertEqual(result, "first")

    def test_check_value_in_two_sets_3(self):
        """Check value that is in the second set."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["D29", "RedRock"])
        value = "D29"
        result = basic.check_value_in_two_sets(value, set1, set2)
        self.assertEqual(result, "second")

    def test_check_value_in_two_sets_4(self):
        """Check value that is in both sets."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["Trixie", "RedRock"])
        value = "Trixie"
        result = basic.check_value_in_two_sets(value, set1, set2)
        self.assertEqual(result, "both")




    def test_lower_case_1(self):
        """Check that a string value is not changed if it is not within
        the set of valid values."""
        input_value = "Trixie"
        output_value = basic.lower_case(input_value)
        self.assertEqual(output_value, input_value)

    def test_lower_case_2(self):
        """Check that a string value is changed if it is within
        the set of valid values."""
        input_value = "NONE"
        output_value = basic.lower_case(input_value)
        self.assertEqual(output_value, "none")

    def test_lower_case_2(self):
        """Check that an int value is changed if it is within
        the set of valid values."""
        input_value = 1
        output_value = basic.lower_case(input_value)
        self.assertEqual(output_value, 1)




    def test_truncate_value_1(self):
        """Verify short value is not changed."""
        input_value = "short"
        output_value = basic.truncate_value(input_value, 6, "...")
        self.assertEqual(output_value, input_value)

    def test_truncate_value_2(self):
        """Verify long value is changed with suffix shorter than max length."""
        input_value = "abcdefghijkl"
        output_value = basic.truncate_value(input_value, 6, "...")
        exp_output_value = "abc..."
        self.assertEqual(output_value, exp_output_value)

    def test_truncate_value_3(self):
        """Verify long value is changed with suffix longer than max length."""
        input_value = "abcdefghijkl"
        output_value = basic.truncate_value(input_value, 6, ".......")
        exp_output_value = "......"
        self.assertEqual(output_value, exp_output_value)



    def test_convert_list_to_dict_1(self):
        """Verify list with dictionaries with unique intended keys is
        converted correctly."""
        input_list = [{"name":"Trixie", "host":"Mycobacterium"},
                      {"name":"L5", "host":"Gordonia"}]
        output_dict = basic.convert_list_to_dict(input_list, "name")
        self.assertEqual(output_dict.keys(), set(["Trixie", "L5"]))

    def test_convert_list_to_dict_2(self):
        """Verify list with dictionaries with non-unique intended keys is
        not converted correctly."""
        input_list = [{"name":"Trixie", "host":"Mycobacterium"},
                      {"name":"Trixie", "host":"Gordonia"}]
        output_dict = basic.convert_list_to_dict(input_list, "name")
        self.assertEqual(len(output_dict.keys()), 0)




    def test_join_strings_1(self):
        """Verify two strings joined as expected."""
        input_list = ["string1", "string2"]
        output = basic.join_strings(input_list)
        expected = "string1 string2"
        self.assertEqual(output, expected)

    def test_join_strings_2(self):
        """Verify empty and non-empty strings are joined as expected."""
        input_list = ["", "string2"]
        output = basic.join_strings(input_list)
        expected = "string2"
        self.assertEqual(output, expected)

    def test_join_strings_3(self):
        """Verify non-empty and empty strings are joined as expected."""
        input_list = ["string1", ""]
        output = basic.join_strings(input_list)
        expected = "string1"
        self.assertEqual(output, expected)

    def test_join_strings_4(self):
        """Verify None and string are joined as expected."""
        input_list = [None, "string2"]
        output = basic.join_strings(input_list)
        expected = "string2"
        self.assertEqual(output, expected)

    def test_join_strings_5(self):
        """Verify multiple values are joined as expected."""
        input_list = ["", "string1", "", None, "string2"]
        output = basic.join_strings(input_list)
        expected = "string1 string2"
        self.assertEqual(output, expected)

    def test_join_strings_6(self):
        """Verify delimiter change."""
        input_list = ["string1", "string2"]
        output = basic.join_strings(input_list, delimiter="; ")
        expected = "string1; string2"
        self.assertEqual(output, expected)




    def test_is_float_1(self):
        """Verify integer string is True."""
        output = basic.is_float("1")
        self.assertTrue(output)

    def test_is_float_2(self):
        """Verify alpha string is False."""
        output = basic.is_float("ABC")
        self.assertFalse(output)

    def test_is_float_3(self):
        """Verify decimal string is False."""
        output = basic.is_float("1.1")
        self.assertTrue(output)



class TestBasicFunctions2(unittest.TestCase):

    def setUp(self):
        self.set1 = {1, 2}
        self.set2 = {3, 4}
        self.set3 = {1, 5}
        self.set4 = {6, 7}
        self.dict1 = {"a": self.set1, "b": self.set2}
        self.dict2 = {"a": self.set3, "c": self.set4}
        self.dict3 = {}

    def test_merge_set_dicts_1(self):
        """Verify new dict contains all keys and merged sets
        two non-empty dictionaries with non-identical keys."""
        dict4 = basic.merge_set_dicts(self.dict1, self.dict2)
        keys4 = {"a", "b", "c"}
        set5 = {1, 2, 5}
        set6 = {3, 4}
        set7 = {6, 7}
        with self.subTest():
            self.assertEqual(dict4.keys(), keys4)
        with self.subTest():
            self.assertEqual(dict4["a"], set5)
        with self.subTest():
            self.assertEqual(dict4["b"], set6)
        with self.subTest():
            self.assertEqual(dict4["c"], set7)
        with self.subTest():
            self.assertEqual(self.dict1.keys(), {"a", "b"})
        with self.subTest():
            self.assertEqual(self.dict1["a"], self.set1)
        with self.subTest():
            self.assertEqual(self.dict1["b"], self.set2)
        with self.subTest():
            self.assertEqual(self.dict2.keys(), {"a", "c"})
        with self.subTest():
            self.assertEqual(self.dict2["a"], self.set3)
        with self.subTest():
            self.assertEqual(self.dict2["c"], self.set4)

    def test_merge_set_dicts_2(self):
        """Verify new dict contains all keys and merged sets from
        one non-empty dictionary and one empty dictionary."""
        dict4 = basic.merge_set_dicts(self.dict1, self.dict3)
        keys3 = {"a", "b"}
        set5 = {1, 2}
        set6 = {3, 4}
        with self.subTest():
            self.assertEqual(dict4.keys(), keys3)
        with self.subTest():
            self.assertEqual(dict4["a"], set5)
        with self.subTest():
            self.assertEqual(dict4["b"], set6)
        with self.subTest():
            self.assertEqual(self.dict1.keys(), {"a", "b"})
        with self.subTest():
            self.assertEqual(self.dict1["a"], set5)
        with self.subTest():
            self.assertEqual(self.dict1["b"], set6)
        with self.subTest():
            self.assertEqual(self.dict3.keys(), set())




if __name__ == '__main__':
    unittest.main()
