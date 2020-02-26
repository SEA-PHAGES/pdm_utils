from pdm_utils.functions import parsing
from sqlalchemy import Column
from unittest.mock import Mock, patch, PropertyMock
import unittest
import re

class TestParsing(unittest.TestCase):
    def test_parse_out_ends_1(self):
        string = "Trixie"

        trimmed_string = parsing.parse_out_ends(string)

        self.assertEqual(string, trimmed_string)

    def test_parse_out_ends_2(self):
        string = " Trixie"

        trimmed_string = parsing.parse_out_ends(string)

        self.assertEqual("Trixie", trimmed_string)
    def test_parse_out_ends_3(self):
        string = " Trixie "

        trimmed_string = parsing.parse_out_ends(string)
        
        self.assertEqual("Trixie", trimmed_string)

    def test_parse_out_ends_4(self):
        with self.assertRaises(TypeError):
            parsing.parse_out_ends(["Fail"])

    def test_parse_cmd_string_1(self):
        filters_string = "phage.PhageID=Trixie AND gene.Notes=Antirepressor"
        
        parsed_string = parsing.parse_cmd_string(filters_string)

        self.assertTrue(len(parsed_string) == 2)
        self.assertEqual(parsed_string[0], ["phage.PhageID=Trixie"])
        self.assertEqual(parsed_string[1], ["gene.Notes=Antirepressor"])

    def test_parse_cmd_string_2(self):
        filters_string = ("phage.PhageID = Trixie OR "
                          "gene.Notes  = Antirepressor")

        parsed_string = parsing.parse_cmd_string(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID = Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes  = Antirepressor")

    def test_parse_cmd_string_3(self):
        filters_string = ("phage.PhageID LIKE Trixie OR "
                          "gene.Notes IS NOT Antirepressor")
 
        parsed_string = parsing.parse_cmd_string(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID LIKE Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes IS NOT Antirepressor")

    def test_parse_cmd_string_4(self):
        filters_string = ("phage.PhageID LIKE Trixie OR "
                          "gene.Notes IS NOT Antirepressor")
 
        parsed_string = parsing.parse_cmd_string(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID LIKE Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes IS NOT Antirepressor")

    def test_parse_cmd_string_5(self):
        filters_string = ("phage.PhageID LIKE Trixie OR")
        
        with self.assertRaises(ValueError):
            parsed_cmd_line = parsing.parse_cmd_string(filters_string)

    def test_parse_cmd_filter_5(self):
        filters_string = ("phage.PhageID LIKE Trixie")

        parsed_filters = parsing.parse_cmd_string(filters_string)
    

    def test_parse_column_1(self):
        column_string = "phage.PhageID"

        parsed_column = parsing.parse_column(column_string)
        
        self.assertEqual(parsed_column[0], "phage")
        self.assertEqual(parsed_column[1], "PhageID")

    def test_parse_column_2(self):
        column_string = " phage.PhageID "

        parsed_column = parsing.parse_column(column_string)
        
        self.assertEqual(parsed_column[0], "phage")
        self.assertEqual(parsed_column[1], "PhageID")

    def test_parse_column_3(self):
        column_string = "phage. PhageID"

        with self.assertRaises(ValueError):
            parsed_column = parsing.parse_column(column_string)
        
    def test_parse_filter_1(self):
        filter = "phage.PhageID=Trixie"

        parsed_filter = parsing.parse_filter(filter)

        self.assertTrue(len(parsed_filter) == 4)
        self.assertEqual(parsed_filter[0], "phage")
        self.assertEqual(parsed_filter[1], "PhageID")
        self.assertEqual(parsed_filter[2], "=")
        self.assertEqual(parsed_filter[3], "Trixie")

    def test_parse_filter_2(self):
        filter = "phage.PhageID = Trixie"

        parsed_filter = parsing.parse_filter(filter)

        self.assertTrue(len(parsed_filter) == 4)
        self.assertEqual(parsed_filter[0], "phage")
        self.assertEqual(parsed_filter[1], "PhageID")
        self.assertEqual(parsed_filter[2], "=")
        self.assertEqual(parsed_filter[3], "Trixie")

    def test_parse_filter_3(self):
        filter = "phage. PhageID = Trixie"

        with self.assertRaises(ValueError):
            parsed_filter = parsing.parse_filter(filter)

    def test_parse_filter3(self):
        filter = "phage.PhageID Trixie = "

        with self.assertRaises(ValueError):
            parsed_filter = parsing.parse_filter(filter)

    def test_check_operator_1(self):
        pass

    def test_check_operator_2(self):
        pass
    
    def test_check_operator_3(self):
        pass

    def test_check_operator_4(self):
        pass

    def test_check_operator_5(self):
        pass

    def test_check_operator_6(self):
        pass

if __name__ == "__main__":
    unittest.main()

