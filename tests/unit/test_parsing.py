from pdm_utils.functions import parsing
from unittest.mock import Mock, patch, PropertyMock
import unittest
import re

class TestParsing(unittest.TestCase):
    def test_parse_cmd_filters_1(self):
        filters_string = "phage.PhageID=Trixie AND gene.Notes=Antirepressor"
        
        parsed_string = parsing.parse_cmd_filters(filters_string)

        self.assertTrue(len(parsed_string) == 2)
        self.assertEqual(parsed_string[0], ["phage.PhageID=Trixie"])
        self.assertEqual(parsed_string[1], ["gene.Notes=Antirepressor"])

    def test_parse_cmd_filters_2(self):
        filters_string = ("phage.PhageID = Trixie OR "
                          "gene.Notes  = Antirepressor")

        parsed_string = parsing.parse_cmd_filters(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID = Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes  = Antirepressor")

    def test_parse_cmd_filters_3(self):
        filters_string = ("phage.PhageID LIKE Trixie OR "
                          "gene.Notes IS NOT Antirepressor")
 
        parsed_string = parsing.parse_cmd_filters(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID LIKE Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes IS NOT Antirepressor")

    def test_parse_cmd_filters_4(self):
        filters_string = ("phage.PhageID LIKE Trixie OR "
                          "gene.Notes IS NOT Antirepressor")
 
        parsed_string = parsing.parse_cmd_filters(filters_string)

        self.assertTrue(len(parsed_string) == 1)
        self.assertTrue(len(parsed_string[0]) == 2)
        self.assertEqual(parsed_string[0][0], "phage.PhageID LIKE Trixie")
        self.assertEqual(parsed_string[0][1], "gene.Notes IS NOT Antirepressor")


    def test_random(self):
        pass 





if __name__ == "__main__":
    unittest.main()

